/*
  wendy.c: One time-step of a one-dimensional N-body code
*/
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>
#include <wendy.h>
#include <bst.h>
bool _PRINT_DIAGNOSTICS= false;
bool _PRINT_ALL_DIAGNOSTICS= false;
double _solve_coll_quad(double c0, double c1, double c2){
  // Solves for collisions under quadratic motion: a t^2/2 + vt + x
  double mba;
  mba= -c1/c2;
  return 0.5 * ( mba + sqrt ( mba * mba - 4. * c0/c2) );
}
double _solve_coll_harm(double c0, double c1, double c2, double omega){
  // Solves for collisions under harmonic motion: A cos(omegaxt+phi)+a/omega^2
  double A, B, out;
  A= c0-c2;
  B= c1/omega;
  out= -asin(c2/sqrt(B*B+A*A))-atan2(A,B);
  if ( out < 0 ) {
    out= M_PI+asin(c2/sqrt(B*B+A*A))-atan2(A,B);
  }
  //printf("Solves equation22? %g,%g,%g,%g\n",A*cos(out),B*sin(out),C,
  // A * cos(out) + B * sin(out) - C);
  //fflush(stdout);
  return out/omega;
}
void _wendy_nbody_onestep(int N, double * x, double * v, double * a,
			  double * m, int * sindx,
			  int * cindx, double * next_tcoll, double * tcoll,
			  int M, double * xt, double * vt, int * stindx, 
			  double dt,
			  int maxcoll, int maxcoll_tp, int * err, 
			  int * ncoll, double * time_elapsed,
			  int * ncoll_tp){
  int cnt_coll,cnt_coll_tp,ii, tmpi, jj, kk, tmpN;
  int c_in_x_indx, c_in_x_next_indx, ctindx, newctindx, cindx_tp, c_in_xt_indx;
  double dv,tdt, dm, tmpd;
  double * next_tcoll_tm;
  struct node* minNode;
  struct node* bst_tcollt;
  int * ntp_left;
  double * tcollt_left, * tcollt_rite;
  struct node** bst_tcollt_left, ** bst_tcollt_rite;
  int tp_assign_low, tp_assign_high;
  double * tcollt_min;
  bool tp_left;
  bool tp_updated= false;
  double * t= (double *) malloc ( N * sizeof(double) );
#pragma omp parallel for schedule(static,chunk) private(ii)
  for (ii=0; ii < N; ii++) *(t+ii)= 0.;
  double * ttp;
  cnt_coll= 0;
  cnt_coll_tp= 0;
  // Build binary search tree for keeping track of collision times
  int * idx= (int *) malloc ( (N-1) * sizeof(int) );
  int * idx_tp;
#pragma omp parallel for schedule(static,chunk) private(ii)
  for (ii=0; ii < N-1; ii++) *(idx+ii)= ii;
  struct node* bst_tcoll= bst_build(N-1,idx,tcoll);
  free(idx);
  //bst_inorder(bst_tcoll);
  if ( M > 0 ) {
    // To keep track of time at which each test-particle is
    ttp= (double *) malloc ( M * sizeof(double) );
#pragma omp parallel for schedule(static,chunk) private(ii)
    for (ii=0; ii < M; ii++) *(ttp+ii)= 0.;
    // Build left/right binary search trees for test particles for every mass
    ntp_left= (int *) malloc ( N * sizeof(int) );
    tcollt_left= (double *) malloc ( M * sizeof(double) );
    tcollt_rite= (double *) malloc ( M * sizeof(double) );
    bst_tcollt_left= (struct node **) malloc ( N * sizeof(struct node*) );
    bst_tcollt_rite= (struct node **) malloc ( N * sizeof(struct node*) );
    idx_tp= (int *) malloc ( M * sizeof(int) );
#pragma omp parallel for schedule(static,chunk) private(ii)
    for (ii=0; ii < M; ii++) *(idx_tp+ii)= ii;
    tp_assign_high= 0;
    jj= 0;
    ii= 0;
    for (jj=0; jj < N+1; jj++) {
      tp_assign_low= tp_assign_high;
      while (  ( ii < M ) &&						\
	       ( ( jj < N && *(xt + *(stindx+ii)) < *(x + *(sindx+jj) ) ) \
		 || ( jj == N ) ) )
	ii++;
      tp_assign_high= ii;
      tmpN= tp_assign_high - tp_assign_low;
      if ( jj > 0 ) {
	for (kk=0; kk < tmpN; kk++)
	  *(tcollt_rite+tp_assign_low+kk) = _solve_coll_quad(		\
		*(xt + *(stindx+tp_assign_low+kk)) - *(x + *(sindx+jj-1)),
		*(vt + *(stindx+tp_assign_low+kk)) - *(v + *(sindx+jj-1)),
		-0.5 * *(m + *(sindx+jj-1))); //rel. a = indiv. mass
	*(bst_tcollt_rite+jj-1)= bst_build(tp_assign_high\
					   -tp_assign_low,
					   idx_tp+tp_assign_low,
					   tcollt_rite+tp_assign_low);
      }
      if ( jj < N ) {
	*(ntp_left+jj)= tp_assign_high;
	for (kk=0; kk < tmpN; kk++)
	  *(tcollt_left+tp_assign_low+kk) = _solve_coll_quad(		\
		*(xt + *(stindx+tp_assign_low+kk)) - *(x + *(sindx+jj)),
		*(vt + *(stindx+tp_assign_low+kk)) - *(v + *(sindx+jj)),
		0.5 * *(m + *(sindx+jj))); //rel. a = indiv. mass
	*(bst_tcollt_left+jj)=bst_build(tp_assign_high-tp_assign_low,
					idx_tp+tp_assign_low,
					tcollt_left+tp_assign_low);
      } 
    }
    //Now make a BST with the minimum times for each node
    idx= (int *) malloc ( 2*N * sizeof(int) );
#pragma omp parallel for schedule(static,chunk) private(ii)
    for (ii=0; ii < 2*N; ii++) *(idx+ii)= ii;
    tcollt_min= (double *) malloc ( 2*N * sizeof(double) );
    for (jj=0; jj < N; jj++) {
      if ( *(bst_tcollt_left+jj) )
	*(tcollt_min+jj)= *(bst_minValueNode(*(bst_tcollt_left+jj))->val);
      else
	*(tcollt_min+jj)= -1;
      if ( *(bst_tcollt_rite+jj) )
	*(tcollt_min+N+jj)= *(bst_minValueNode(*(bst_tcollt_rite+jj))->val);
      else
	*(tcollt_min+N+jj)= -1;
    }
    bst_tcollt= bst_build(2*N,idx,tcollt_min);
    next_tcoll_tm= (double *) malloc ( sizeof ( double ) );
    *next_tcoll_tm= *(bst_minValueNode(bst_tcollt)->val);
    free(idx);
    /*
    if ( _PRINT_DIAGNOSTICS ) {
	printf("Printing tp BST in order\n");
	fflush(stdout);
	bst_inorder(bst_tcollt);
    }
    */
  } else next_tcoll_tm= &dt; // just so that it is >= dt
  // Time how long the loop takes
  clock_t time_begin= clock();
  while ( ( *next_tcoll < dt || *next_tcoll_tm < dt ) \
	  && cnt_coll < maxcoll && cnt_coll_tp < maxcoll_tp ){
    /*
      if ( _PRINT_DIAGNOSTICS ) {
    printf("Next collisions? %g,%g\n",*next_tcoll,*next_tcoll_tm);
    fflush(stdout);
    */
    if ( M == 0 || *next_tcoll < *next_tcoll_tm ) {
      cnt_coll+= 1;
      // collide, update collided particles
      c_in_x_indx= *(sindx+ *cindx);
      c_in_x_next_indx= *(sindx+ *cindx+1);
      tdt= ( *next_tcoll - *(t + c_in_x_indx) );
      dv= *(a + c_in_x_indx) * tdt;
      *(x + c_in_x_indx)+= dv * tdt / 2. + *(v + c_in_x_indx) * tdt;
      *(v + c_in_x_indx)+= dv;
      *(t + c_in_x_indx)= *next_tcoll;
      tdt= ( *next_tcoll - *(t + c_in_x_next_indx) );
      dv= *(a + c_in_x_next_indx) * tdt;
      *(x + c_in_x_next_indx)+= dv * tdt / 2. + *(v + c_in_x_next_indx) * tdt;
      *(v + c_in_x_next_indx)+= dv;
      *(t + c_in_x_next_indx)= *next_tcoll;
      // swap
      tmpi= *(sindx+ *cindx);
      *(sindx+ *cindx)= *(sindx+ *cindx+1);
      *(sindx+ *cindx+1)= tmpi;
      // track mass and update accelerations
      dm= *(m + c_in_x_next_indx) - *(m + c_in_x_indx);
      *(a + c_in_x_indx)-= dm;
      *(a + c_in_x_next_indx)-= dm;
      tmpd= *(a + c_in_x_indx);
      *(a + c_in_x_indx)= *(a + c_in_x_next_indx);
      *(a + c_in_x_next_indx)= tmpd;
      // Update collision times, solution for delta x = 0
      //printf("Collide in %f\n",2.* (*(v + c_in_x_indx) - *(v + c_in_x_next_indx)) /
      //	   (*(a + c_in_x_next_indx) - *(a + c_in_x_indx)));
      //fflush(stdout);
      tmpd= 2.* (*(v + c_in_x_indx) - *(v + c_in_x_next_indx)) / 
	(*(a + c_in_x_next_indx) - *(a + c_in_x_indx));
      bst_tcoll= bst_deleteNode(bst_tcoll,tcoll+*cindx);
      *(tcoll + *cindx)= *next_tcoll+tmpd;
      bst_tcoll= bst_forceInsert(bst_tcoll,*cindx,tcoll+*cindx);
      // Abuse the c_in_x_indx and c_in_x_next_indx arrays
      if ( *cindx > 0 ){
	c_in_x_indx= *(sindx+ *cindx-1);
	c_in_x_next_indx= *(sindx+ *cindx);
	tdt= *(t + c_in_x_indx) - *next_tcoll;
	bst_tcoll= bst_deleteNode(bst_tcoll,tcoll+*cindx-1);
	*(tcoll + *cindx -1)= *next_tcoll +
	  _solve_coll_quad(*(x + c_in_x_indx) + *(a+c_in_x_indx) * tdt*tdt / 2. 
			   - *(v+c_in_x_indx) * tdt - *(x + c_in_x_next_indx),
			   *(v + c_in_x_indx) - *(a + c_in_x_indx) * tdt
			   - *(v + c_in_x_next_indx),
			   0.5*(*(a + c_in_x_indx) - *(a + c_in_x_next_indx)));
	bst_tcoll= bst_forceInsert(bst_tcoll,*cindx-1,tcoll+*cindx-1);
      }
      if ( *cindx < N-2 ){
	c_in_x_indx= *(sindx+ *cindx+2);
	c_in_x_next_indx= *(sindx+ *cindx+1);
	tdt= *(t + c_in_x_indx) - *next_tcoll;
	bst_tcoll= bst_deleteNode(bst_tcoll,tcoll+*cindx+1);
	*(tcoll + *cindx+1)= *next_tcoll +
	  _solve_coll_quad(*(x + c_in_x_indx) + *(a+c_in_x_indx) * tdt*tdt / 2. 
			   - *(v+c_in_x_indx) * tdt - *(x + c_in_x_next_indx),
			   *(v + c_in_x_indx) - *(a + c_in_x_indx) * tdt
			   - *(v + c_in_x_next_indx),
			   0.5*(*(a + c_in_x_indx) - *(a + c_in_x_next_indx)));
	bst_tcoll= bst_forceInsert(bst_tcoll,*cindx+1,tcoll+*cindx+1);
      }
      //Also update test particles associated with the swapped pair
      if ( M > 0 ) {
	//Assign left tp from left mp to left tp of rite mp
	if ( *cindx > 0 ) {
	  tmpN= *(ntp_left+*cindx)-*(ntp_left+*cindx-1);
	}
	else {
	  tmpN= *(ntp_left+*cindx);
	}
	if ( tmpN > 0 ) { // there is a left tree to move
	  if ( *cindx > 0 ) {
	    tmpi= *(ntp_left+*cindx-1);
	  }
	  else {
	    tmpi= 0;
	  }
	  c_in_x_indx= *(sindx+ *cindx+1); // bc already flipped
	  c_in_x_next_indx= *(sindx+ *cindx);
	  //First delete node in overall bst_tcollt, left of left, rite of rite
	  if ( *(bst_tcollt_left+*cindx) )
	    bst_tcollt= bst_deleteNode(bst_tcollt,tcollt_min+*cindx);
	  tmpd= *(a + c_in_x_next_indx) + *(m + c_in_x_next_indx); // tp accel.
	  for (kk=0; kk < tmpN; kk++) {
	    tdt= *next_tcoll - *(ttp + *(stindx+tmpi+kk));
	    *(tcollt_left+tmpi+kk) = *next_tcoll + _solve_coll_quad(	\
	    *(xt + *(stindx+tmpi+kk)) + *(vt + *(stindx+tmpi+kk)) * tdt \
	    + 0.5 * tmpd * tdt * tdt - *(x + c_in_x_next_indx),
	    *(vt + *(stindx+tmpi+kk)) + tmpd * tdt - *(v + c_in_x_next_indx),
	    0.5 * *(m + c_in_x_next_indx)); //rel. a = indiv. mas

	  if ( *(tcollt_left+tmpi+kk)-*next_tcoll < 0 ) {
	    //	    if ( _PRINT_DIAGNOSTICS ) {
	    printf("Printing tp bst in order at end\n");
	    bst_inorder(bst_tcollt);
	    printf("Next collisions? %g,%g in %g\n",
		   *next_tcoll,*next_tcoll_tm,dt);
	    fflush(stdout);
	    bst_inorder(*(bst_tcollt_left+*cindx));
	      printf("Left collision data: %i,%i,%i,%g,%g,%g,%g,%g,%g,%g,%g,%g\n",
		     *cindx,tmpi+kk,tmpN,
		     *(xt + *(stindx+tmpi+kk)),
		     *(vt + *(stindx+tmpi+kk)),
		     tdt,
		     tmpd,
		     *(x + c_in_x_next_indx),
		     *(v + c_in_x_next_indx),
		     *(xt + *(stindx+tmpi+kk))		\
		     + *(vt + *(stindx+tmpi+kk)) * tdt	\
		     + 0.5 * tmpd * tdt * tdt,
		     *(vt + *(stindx+tmpi+kk)) + tmpd * tdt,
		     *(ttp + *(stindx+tmpi+kk)));
	      printf("\033[1m\033[30m" " NEGATIVE LEFT COLLISION TIME %g \n" "\033[0m",*(tcollt_left+tmpi+kk)-*next_tcoll);
	      fflush(stdout);
	      //	    }
	      abort();
	  }

	  }
	  *(bst_tcollt_left+*cindx)= bst_build(tmpN,
					       idx_tp+tmpi,
					       tcollt_left+tmpi);
	  bst_destroy(*(bst_tcollt_rite+*cindx));
	  *(bst_tcollt_rite+*cindx)= NULL;
	  //insert new tcoll in overall bst_tcollt
	  *(tcollt_min+*cindx)=						\
	    *(bst_minValueNode(*(bst_tcollt_left+*cindx))->val);
	  bst_tcollt= bst_forceInsert(bst_tcollt,*cindx,tcollt_min+*cindx);
	  tp_updated= true;
	}
	//Assign rite tp from rite mp to rite tp of left mp
	if ( *cindx < N-2 ) {
	  tmpN= *(ntp_left+*cindx+2)-*(ntp_left+*cindx+1);
	}
	else {
	  tmpN= M-*(ntp_left+*cindx+1);
	}
	if ( tmpN > 0 ) { // there is a right tree to move
	  if ( *cindx < N-2 ) {
	    tmpi= *(ntp_left+*cindx+1);
	  }
	  else {
	    tmpi= *(ntp_left+*cindx+1);
	  }
	  c_in_x_indx= *(sindx+ *cindx+1); // bc already flipped
	  c_in_x_next_indx= *(sindx+ *cindx);
	  if ( *(bst_tcollt_rite+*cindx+1) )
	    bst_tcollt= bst_deleteNode(bst_tcollt,tcollt_min+N+*cindx+1);
	  tmpd= *(a + c_in_x_indx) - *(m + c_in_x_indx); // tp acceleration
	  for (kk=0; kk < tmpN; kk++) {
	    tdt= *next_tcoll - *(ttp + *(stindx+tmpi+kk));
	    *(tcollt_rite+tmpi+kk) = *next_tcoll + _solve_coll_quad(	\
	    *(xt + *(stindx+tmpi+kk)) + *(vt + *(stindx+tmpi+kk)) * tdt \
	    + 0.5 * tmpd * tdt * tdt- *(x + c_in_x_indx),
	    *(vt + *(stindx+tmpi+kk)) + tmpd * tdt - *(v + c_in_x_indx),
	    -0.5 * *(m + c_in_x_indx)); //rel. a = indiv. mass
	    if ( *(tcollt_rite+tmpi+kk)-*next_tcoll < 0 ) {
	      if ( _PRINT_DIAGNOSTICS ) {
		printf("Printing tp bst in order\n");
		bst_inorder(bst_tcollt);
		printf("Next collisions? %g,%g in %g\n",
		       *next_tcoll,*next_tcoll_tm,dt);
		fflush(stdout);
		bst_inorder(*(bst_tcollt_rite+*cindx+1));
		printf("Rite collision data: %i,%i,%i,%g,%g,%g,%g,%g,%g,%g,%g,%g\n",
		       *cindx,tmpi+kk,tmpN,
		       *(xt + *(stindx+tmpi+kk)),
		       *(vt + *(stindx+tmpi+kk)),
		       tdt,
		       tmpd,
		       *(x + c_in_x_next_indx),
		       *(v + c_in_x_next_indx),
		       *(xt + *(stindx+tmpi+kk))	\
		       + *(vt + *(stindx+tmpi+kk)) * tdt	\
		       + 0.5 * tmpd * tdt * tdt,
		       *(vt + *(stindx+tmpi+kk)) + tmpd * tdt,
		       *(ttp + *(stindx+tmpi+kk)));
		printf("\033[1m\033[30m" " NEGATIVE RITE COLLISION TIME %g\n" "\033[0m",*(tcollt_rite+tmpi+kk)-*next_tcoll);
		fflush(stdout);
	      }
	      abort();
	    }
	    
	  }
	  
	  *(bst_tcollt_rite+*cindx+1)= bst_build(tmpN,
						 idx_tp+tmpi,
						 tcollt_rite+tmpi);
	  bst_destroy(*(bst_tcollt_left+*cindx+1));
	  *(bst_tcollt_left+*cindx+1)= NULL;
	  //insert new tcoll in overall bst_tcollt
	  *(tcollt_min+N+*cindx+1)=					\
	    *(bst_minValueNode(*(bst_tcollt_rite+*cindx+1))->val);
	  bst_tcollt= bst_forceInsert(bst_tcollt,N+*cindx+1,
				      tcollt_min+N+*cindx+1);
	  tp_updated= true;
	}
	/*
	if ( _PRINT_ALL_DIAGNOSTICS ) {
	  for (jj=0; jj < N; jj++) {
	    printf("Printing new rite BST in order for mass %i\n",jj);
	    fflush(stdout);
	    bst_inorder(*(bst_tcollt_rite+jj));
	    printf("Printing new left BST in order for mass %i\n",jj);
	    fflush(stdout);
	    bst_inorder(*(bst_tcollt_left+jj));
	  }
	  printf("Printing tp BST in order\n");
	  fflush(stdout);
	  bst_inorder(bst_tcollt);
	}
	*/
      }
      // Find minimum
      minNode= bst_minValueNode(bst_tcoll);
      *cindx= minNode->idx;
      *next_tcoll= *minNode->val;
      //printf("Next one %f\n",*next_tcoll);
      //fflush(stdout);
    } else if ( *next_tcoll_tm < dt ){
      // Next handle test-particle -- massive-particle collisions
      // Locate colliding tp, delete from its own BSTs
      cnt_coll_tp+= 1;
      minNode= bst_minValueNode(bst_tcollt);
      cindx_tp= minNode->idx;
      bst_tcollt= bst_deleteNode(bst_tcollt,tcollt_min+cindx_tp);
      if ( cindx_tp < N ) { // collision from the left
	tp_left= true;
	if ( cindx_tp > 0 )
	  bst_tcollt= bst_deleteNode(bst_tcollt,tcollt_min+N+cindx_tp-1);
	minNode= bst_minValueNode(*(bst_tcollt_left+cindx_tp));
	ctindx= minNode->idx;
	*(bst_tcollt_left+cindx_tp)=bst_deleteNode(*(bst_tcollt_left+cindx_tp),
						   tcollt_left+ctindx);
	if ( cindx_tp > 0 )
	  *(bst_tcollt_rite+cindx_tp-1)=\
	    bst_deleteNode(*(bst_tcollt_rite+cindx_tp-1),tcollt_rite+ctindx);
      } else { // collision from the right
	tp_left= false;
	cindx_tp-= N;
	if ( cindx_tp < N-1 )
	  bst_tcollt= bst_deleteNode(bst_tcollt,tcollt_min+cindx_tp+1);
	minNode= bst_minValueNode(*(bst_tcollt_rite+cindx_tp));
	ctindx= minNode->idx;
	*(bst_tcollt_rite+cindx_tp)=bst_deleteNode(*(bst_tcollt_rite+cindx_tp),
						   tcollt_rite+ctindx);
	if ( cindx_tp < N-1 )
	  *(bst_tcollt_left+cindx_tp+1)=\
	    bst_deleteNode(*(bst_tcollt_left+cindx_tp+1),tcollt_left+ctindx);
      }
      //printf("Collision at %i from left: %i\n",cindx_tp,(int)tp_left);
      // Recompute min. left and rite and insert in overall tree
      if ( cindx_tp < N-1+(int)tp_left ) {
	if ( *(bst_tcollt_left+cindx_tp+1-(int)tp_left) ) {
	  *(tcollt_min+cindx_tp+1-(int)tp_left)=			\
	    *(bst_minValueNode(*(bst_tcollt_left+cindx_tp+1-(int)tp_left))->val);
	  bst_tcollt= bst_forceInsert(bst_tcollt,cindx_tp+1-(int)tp_left,
				      tcollt_min+cindx_tp+1-(int)tp_left);
	}
      }
      if ( cindx_tp > (int)tp_left-1 ) {
	if ( *(bst_tcollt_rite+cindx_tp-(int)tp_left) ) {
	  *(tcollt_min+N+cindx_tp-(int)tp_left)=			\
	    *(bst_minValueNode(*(bst_tcollt_rite+cindx_tp-(int)tp_left))->val);
	  bst_tcollt= bst_forceInsert(bst_tcollt,N+cindx_tp-(int)tp_left,
				      tcollt_min+N+cindx_tp-(int)tp_left);
	}
      }
      // Move colliding tp to impact point
      c_in_xt_indx= *(stindx+ctindx);
      tdt= ( *next_tcoll_tm - *(ttp + c_in_xt_indx) );
      dv= (*(a + *(sindx+cindx_tp)) \
	   + (2 * (int)tp_left - 1) * *(m + *(sindx+cindx_tp))  ) * tdt;
      /*
      if ( _PRINT_DIAGNOSTICS ) {
	printf("a: %g, dv: %g\n",*(a + *(sindx+cindx_tp)),dv);
	fflush(stdout);
      }
      */
      *(xt + c_in_xt_indx)+= dv * tdt / 2. + *(vt + c_in_xt_indx) * tdt;
      *(vt + c_in_xt_indx)+= dv;
      *(ttp + c_in_xt_indx)= *next_tcoll_tm;
      // move the collided particle to the appropriate edge, so the tp array
      // remains sorted in the between-massive-particle buckets; edit any
      // other use of the sorted indices
      fflush(stdout);
      *(ntp_left+cindx_tp)+= 1 - 2 * (int)tp_left;
      newctindx= *(ntp_left+cindx_tp)-1+(int)tp_left;
      if ( ctindx != newctindx ) { // If moving tps around
	tmpi= *(stindx+ newctindx);
	*(stindx+ newctindx)= *(stindx+ ctindx);
	*(stindx+ ctindx)= tmpi;
	if ( cindx_tp < N-1+(int)tp_left ) {
	  *(bst_tcollt_left+cindx_tp+1-(int)tp_left)=			\
	    bst_deleteNode(*(bst_tcollt_left+cindx_tp+1-(int)tp_left),
			   tcollt_left+newctindx);
	  *(tcollt_left+ctindx)= *(tcollt_left+newctindx);
	  *(bst_tcollt_left+cindx_tp+1-(int)tp_left)=			\
	    bst_forceInsert(*(bst_tcollt_left+cindx_tp+1-(int)tp_left),
			    ctindx,tcollt_left+ctindx);
	}
	if ( cindx_tp > -1+(int)tp_left ) {
	  *(bst_tcollt_rite+cindx_tp-(int)tp_left)=			\
	    bst_deleteNode(*(bst_tcollt_rite+cindx_tp-(int)tp_left),
			   tcollt_rite+newctindx);
	  *(tcollt_rite+ctindx)= *(tcollt_rite+newctindx);
	  *(bst_tcollt_rite+cindx_tp-(int)tp_left)=			\
	    bst_forceInsert(*(bst_tcollt_rite+cindx_tp-(int)tp_left),
			    ctindx,tcollt_rite+ctindx);
	}
      }
      // Compute left/rite collision times for new position and insert
      if ( cindx_tp > -(int)tp_left ) {
	if ( *(bst_tcollt_rite+cindx_tp-1+(int)tp_left) )
	  bst_tcollt= bst_deleteNode(bst_tcollt,
				     tcollt_min+N+cindx_tp-1+(int)tp_left);
	c_in_x_indx= *(sindx+cindx_tp-1+(int)tp_left);
	tdt= ( *next_tcoll_tm - *(t + c_in_x_indx) );
	if ( tp_left ) 
	  *(tcollt_rite+newctindx) = *next_tcoll_tm \
	    +2.* (*(vt + c_in_xt_indx) - *(v + c_in_x_indx)	\
		  - *(a + c_in_x_indx) * tdt) \
	    / *(m + c_in_x_indx);
	else
	  *(tcollt_rite+newctindx) = *next_tcoll_tm			\
	    +_solve_coll_quad(						\
          *(xt + c_in_xt_indx) - *(x + c_in_x_indx) - *(v + c_in_x_indx) * tdt\
	  - 0.5 * *(a + c_in_x_indx) * tdt * tdt,
	  *(vt + c_in_xt_indx) - *(v + c_in_x_indx) - *(a + c_in_x_indx) * tdt,
	  -0.5 * *(m + c_in_x_indx)); //rel. a = indiv. mass
	/*
	if ( _PRINT_DIAGNOSTICS ) {
	  printf("Next rite collision data %g, %g, %g, %g, %g, %g, %g, %g, %g, %i\n",
		 *(xt + c_in_xt_indx),*(vt + c_in_xt_indx),
		 *(x + c_in_x_indx),*(v + c_in_x_indx),*(a + c_in_x_indx),*(m + c_in_x_indx),tdt,*(t + c_in_x_indx),*next_tcoll_tm,cindx_tp+(int)tp_left);
	  printf("Next rite collision %g\n",*(tcollt_rite+newctindx)-*next_tcoll_tm);
	}
	*/
	if (*(tcollt_rite+newctindx)-*next_tcoll_tm < 0) {
	  printf("\033[1m\033[30m" " NEGATIVE RITE COLLISION TIME \n" "\033[0m");
	  fflush(stdout);
	  abort();
	}
	*(bst_tcollt_rite+cindx_tp-1+(int)tp_left)=			\
	  bst_forceInsert(*(bst_tcollt_rite+cindx_tp-1+(int)tp_left),
			  newctindx,tcollt_rite+newctindx);
	*(tcollt_min+N+cindx_tp-1+(int)tp_left)=			\
	  *(bst_minValueNode(*(bst_tcollt_rite+cindx_tp-1+(int)tp_left))->val);
	bst_tcollt= bst_forceInsert(bst_tcollt,N+cindx_tp-1+(int)tp_left,
				    tcollt_min+N+cindx_tp-1+(int)tp_left);
      }
      if ( cindx_tp < N-(int)tp_left ) {
	if ( *(bst_tcollt_left+cindx_tp+(int)tp_left) )
	  bst_tcollt= bst_deleteNode(bst_tcollt,
				     tcollt_min+cindx_tp+(int)tp_left);
	c_in_x_indx= *(sindx+cindx_tp+(int)tp_left);
	tdt= ( *next_tcoll_tm - *(t + c_in_x_indx) );
	if ( tp_left ) 
	  *(tcollt_left+newctindx) = *next_tcoll_tm			\
	    +_solve_coll_quad(						\
          *(xt + c_in_xt_indx) - *(x + c_in_x_indx) - *(v + c_in_x_indx) * tdt\
	  - 0.5 * *(a + c_in_x_indx) * tdt * tdt,
	  *(vt + c_in_xt_indx) - *(v + c_in_x_indx) - *(a + c_in_x_indx) * tdt,
	  0.5 * *(m + c_in_x_indx)); //rel. a = indiv. mass
	else
	  *(tcollt_left+newctindx) = *next_tcoll_tm			\
	    -2. * (*(vt + c_in_xt_indx) - *(v + c_in_x_indx) \
		   - *(a + c_in_x_indx) * tdt)		     \
	    / *(m + c_in_x_indx);
	/*
	if ( _PRINT_DIAGNOSTICS ) {
	  printf("Next left collision data %g, %g, %g, %g, %g, %g, %g, %g, %g, %i,%i,%i\n",
		 *(xt + c_in_xt_indx),*(vt + c_in_xt_indx),
		 *(x + c_in_x_indx),*(v + c_in_x_indx),*(a + c_in_x_indx),*(m + c_in_x_indx),tdt,*(t + c_in_x_indx),*next_tcoll_tm,cindx_tp,(int)tp_left,newctindx);
	  printf("Next left collision %g\n",
		 *(tcollt_left+newctindx)-*next_tcoll_tm);
	}
	*/
	if (*(tcollt_left+newctindx)-*next_tcoll_tm < 0) {
	  printf("\033[1m\033[30m" " NEGATIVE LEFT COLLISION TIME \n" "\033[0m");
	  fflush(stdout);
	  abort();
	}
	*(bst_tcollt_left+cindx_tp+(int)tp_left)=			\
	  bst_forceInsert(*(bst_tcollt_left+cindx_tp+(int)tp_left),
			  newctindx,tcollt_left+newctindx);
	*(tcollt_min+cindx_tp+(int)tp_left)=				\
	  *(bst_minValueNode(*(bst_tcollt_left+cindx_tp+(int)tp_left))->val);
	bst_tcollt= bst_forceInsert(bst_tcollt,cindx_tp+(int)tp_left,
				    tcollt_min+cindx_tp+(int)tp_left);
      }
      tp_updated= true;
    }
  
    if ( M > 0 && tp_updated ) {
      *next_tcoll_tm= *(bst_minValueNode(bst_tcollt)->val);
      tp_updated= false; // reset
    }
  }
  clock_t time_end= clock();
  *time_elapsed= (double) (time_end-time_begin) / CLOCKS_PER_SEC;
  //printf("Next %f\n",*next_tcoll-dt);
  //fflush(stdout);
  // Update all to next snapshot
#pragma omp parallel for schedule(static,chunk) private(ii,tdt,dv)
  for (ii=0; ii < N; ii++) {
    tdt= dt - *(t+ii);
    dv= *(a+ii) * tdt;
    *(x+ii)+= dv * tdt / 2. + *(v+ii) * tdt;
    *(v+ii)+= dv;
  }
#pragma omp parallel for schedule(static,chunk) private(ii)
  for (ii=0; ii < N-1; ii++) {
    *(tcoll+ii)-= dt;
  }
  *next_tcoll-= dt;
  if ( M > 0 ) {
    if ( _PRINT_ALL_DIAGNOSTICS ) {
      for (ii=0; ii < N; ii++)
	printf("ntp_left %i: %i\n",ii,*(ntp_left+ii));
      fflush(stdout);
    }

#pragma omp parallel for schedule(static,chunk) private(ii,jj,tmpd,tmpN,tdt,dv,c_in_xt_indx)
    for (ii=0; ii < N+1; ii++) {
      if ( ii == N ) {
	tmpN= M - *(ntp_left+N-1);
	tmpd= ( *(a + *(sindx+N-1) ) - *(m + *(sindx+N-1) ) );
	tmpi= *(ntp_left+N-1);
      } else {
	tmpd= ( *(a + *(sindx+ii) ) + *(m + *(sindx+ii) ) );
	if ( ii == 0 ) {
	  tmpN= *ntp_left;
	  tmpi= 0;
	} else {
	  tmpN= *(ntp_left+ii)-*(ntp_left+ii-1);
	  tmpi= *(ntp_left+ii-1);
	}
      }
      for (jj=0; jj < tmpN; jj++) {
	c_in_xt_indx= *(stindx + tmpi + jj );
	tdt= dt - *(ttp + c_in_xt_indx);
	if ( _PRINT_DIAGNOSTICS ) {
	  printf("da %i,%i: %g, dt: %g, x: %g, v: %g\n",ii,jj,tmpd,tdt,
		 *(xt + c_in_xt_indx),*(vt + c_in_xt_indx));
	  fflush(stdout);
	}
	dv= tmpd * tdt;
	*(xt + c_in_xt_indx)+= dv * tdt / 2. + *(vt + c_in_xt_indx) * tdt;
	*(vt + c_in_xt_indx)+= dv;
      }
    }
  }
  free(t);
  bst_destroy(bst_tcoll);
  if ( M > 0 ) {
    free(ttp);
    for (jj=0; jj < N; jj++) {
      bst_destroy(*(bst_tcollt_left+jj));
      bst_destroy(*(bst_tcollt_rite+jj));
    }
    free(bst_tcollt_left);
    free(bst_tcollt_rite);
    free(tcollt_left);
    free(tcollt_rite);
    free(ntp_left);
    free(idx_tp);
    bst_destroy(bst_tcollt);
    free(tcollt_min);
    free(next_tcoll_tm);
  }
  *ncoll= cnt_coll;
  *ncoll_tp= cnt_coll_tp;
  if ( cnt_coll == maxcoll )
    *err= -2;
  if ( cnt_coll_tp == maxcoll_tp )
    *err= -3;
}

void _wendy_nbody_harm_onestep(int N, double * x, double * v, double * a,
			       double * m, int * sindx,
			       int * cindx, double * next_tcoll,
			       double * tcoll,double dt,
			       int maxcoll, int * err, 
			       int * ncoll,double * time_elapsed,
			       double omega){
  int cnt_coll,ii, tmpi;
  int c_in_x_indx, c_in_x_next_indx;
  double tdt, dm, tmpd, cosot, sinot;
  struct node* minNode;
  double * t= (double *) malloc ( N * sizeof(double) );
#pragma omp parallel for schedule(static,chunk) private(ii)
  for (ii=0; ii < N; ii++) *(t+ii)= 0.;
  cnt_coll= 0;
  // Build binary search tree for keeping track of collision times
  int * idx= (int *) malloc ( (N-1) * sizeof(int) );
#pragma omp parallel for schedule(static,chunk) private(ii)
  for (ii=0; ii < N-1; ii++) *(idx+ii)= ii;
  struct node* bst_tcoll= bst_build(N-1,idx,tcoll);
  free(idx);
  //bst_inorder(bst_tcoll);
  // Time how long the loop takes
  clock_t time_begin= clock();
  while ( *next_tcoll < dt && cnt_coll < maxcoll ){
    //printf("Colliding in %f\n",*next_tcoll);
    //fflush(stdout);
    cnt_coll+= 1;
    // collide, update collided particles
    c_in_x_indx= *(sindx+ *cindx);
    c_in_x_next_indx= *(sindx+ *cindx+1);
    tdt= ( *next_tcoll - *(t + c_in_x_indx) );
    sinot = sin( omega * tdt );
    cosot= sqrt(1.-sinot*sinot);
    tmpd= *(x + c_in_x_indx) - *(a + c_in_x_indx);
    *(x + c_in_x_indx)= tmpd * cosot \
      + *(v + c_in_x_indx) / omega * sinot \
      + *(a + c_in_x_indx);
    *(v + c_in_x_indx)= -tmpd * omega * sinot \
      + *(v + c_in_x_indx) * cosot;
    *(t + c_in_x_indx)= *next_tcoll;
    tdt= ( *next_tcoll - *(t + c_in_x_next_indx) );
    sinot = sin( omega * tdt );
    cosot= sqrt(1.-sinot*sinot);
    tmpd= *(x + c_in_x_next_indx) - *(a + c_in_x_next_indx);
    *(x + c_in_x_next_indx)= tmpd * cosot \
      + *(v + c_in_x_next_indx) / omega * sinot \
      + *(a + c_in_x_next_indx) ;
    *(v + c_in_x_next_indx)= -tmpd * omega * sinot \
      + *(v + c_in_x_next_indx) * cosot;
    *(t + c_in_x_next_indx)= *next_tcoll;
    //printf("Collide? %g, %g, %g\n",*(x + c_in_x_indx),*(x + c_in_x_next_indx),
    //   *(x + c_in_x_indx)-*(x + c_in_x_next_indx));
    //fflush(stdout);
    // swap
    tmpi= *(sindx+ *cindx);
    *(sindx+ *cindx)= *(sindx+ *cindx+1);
    *(sindx+ *cindx+1)= tmpi;
    // track mass and update accelerations
    dm= *(m + c_in_x_next_indx) - *(m + c_in_x_indx);
    *(a + c_in_x_indx)-= dm;
    *(a + c_in_x_next_indx)-= dm;
    tmpd= *(a + c_in_x_indx);
    *(a + c_in_x_indx)= *(a + c_in_x_next_indx);
    *(a + c_in_x_next_indx)= tmpd;
    // Update collision times, solution for delta x = 0
    tmpd= (*(v + c_in_x_indx) - *(v + c_in_x_next_indx)) / \
      (*(a + c_in_x_next_indx) - *(a + c_in_x_indx)) / omega;
    tmpd*= tmpd;
    tmpd= (1.-tmpd)/(1.+tmpd);
    tmpd= acos(tmpd)/omega;
    bst_tcoll= bst_deleteNode(bst_tcoll,tcoll+*cindx);
    *(tcoll + *cindx)= *next_tcoll+tmpd;
    bst_tcoll= bst_forceInsert(bst_tcoll,*cindx,tcoll+*cindx);
    //printf("Collide in %f\n",*(tcoll + *cindx)-*next_tcoll);
    //fflush(stdout);
    // Abuse the c_in_x_indx and c_in_x_next_indx arrays
    if ( *cindx > 0 ){
      c_in_x_indx= *(sindx+ *cindx-1);
      c_in_x_next_indx= *(sindx+ *cindx);
      // Also forward the previous sheet to the collision time for convenience
      tdt= *next_tcoll - *(t + c_in_x_indx);
      sinot = sin( omega * tdt );
      cosot= sqrt(1.-sinot*sinot);
      tmpd= *(x + c_in_x_indx) - *(a + c_in_x_indx);
      *(x + c_in_x_indx)= tmpd * cosot	   \
	      + *(v + c_in_x_indx) / omega * sinot	\
	      + *(a + c_in_x_indx);
      *(v + c_in_x_indx)= -tmpd * omega * sinot		\
	+ *(v + c_in_x_indx) * cosot;
      *(t + c_in_x_indx)= *next_tcoll;
      bst_tcoll= bst_deleteNode(bst_tcoll,tcoll+*cindx-1);
      *(tcoll + *cindx -1)= *next_tcoll +
	_solve_coll_harm(*(x + c_in_x_indx) - *(x + c_in_x_next_indx),
			 *(v + c_in_x_indx) - *(v + c_in_x_next_indx),
			 *(a + c_in_x_indx) - *(a + c_in_x_next_indx),
			 omega);
      bst_tcoll= bst_forceInsert(bst_tcoll,*cindx-1,tcoll+*cindx-1);
    }
    if ( *cindx < N-2 ){
      c_in_x_indx= *(sindx+ *cindx+2);
      c_in_x_next_indx= *(sindx+ *cindx+1);
      // Also forward the next sheet to the collision time for convenience
      tdt= *next_tcoll - *(t + c_in_x_indx);
      sinot = sin( omega * tdt );
      cosot= sqrt(1.-sinot*sinot);
      tmpd= *(x + c_in_x_indx) - *(a + c_in_x_indx);
      *(x + c_in_x_indx)= tmpd * cosot	   \
	+ *(v + c_in_x_indx) / omega * sinot	\
	+ *(a + c_in_x_indx);
      *(v + c_in_x_indx)= -tmpd * omega * sinot	\
	+ *(v + c_in_x_indx) * cosot;
      *(t + c_in_x_indx)= *next_tcoll;
      bst_tcoll= bst_deleteNode(bst_tcoll,tcoll+*cindx+1);
      *(tcoll + *cindx+1)= *next_tcoll +
	_solve_coll_harm(*(x + c_in_x_indx) - *(x + c_in_x_next_indx),
			 *(v + c_in_x_indx) - *(v + c_in_x_next_indx),
			 *(a + c_in_x_indx) - *(a + c_in_x_next_indx),
			 omega);
      bst_tcoll= bst_forceInsert(bst_tcoll,*cindx+1,tcoll+*cindx+1);
    }
    // Find minimum
    minNode= bst_minValueNode(bst_tcoll);
    *cindx= minNode->idx;
    *next_tcoll= *minNode->val;
    //printf("Next one %f\n",*next_tcoll);
    //fflush(stdout);
  }
  clock_t time_end= clock();
  *time_elapsed= (double) (time_end-time_begin) / CLOCKS_PER_SEC;
  //printf("Next %f\n",*next_tcoll-dt);
  //fflush(stdout);
  // Update all to next snapshot
#pragma omp parallel for schedule(static,chunk) private(ii,tdt,dv,cosot,sinot,tmpd)
  for (ii=0; ii < N; ii++) {
    tdt= dt - *(t+ii);
    sinot = sin( omega * tdt );
    cosot= sqrt(1.-sinot*sinot);
    tmpd= *(x+ii) - *(a+ii);
    *(x+ii)= tmpd * cosot + *(v+ii) / omega * sinot + *(a+ii);
    *(v+ii)= -tmpd * omega * sinot + *(v+ii) * cosot;
  }
#pragma omp parallel for schedule(static,chunk) private(ii)
  for (ii=0; ii < N-1; ii++) {
    *(tcoll+ii)-= dt;
  }
  *next_tcoll-= dt;
  free(t);
  bst_destroy(bst_tcoll);
  *ncoll= cnt_coll;
  if ( cnt_coll == maxcoll )
    *err= -2;
}

// Approximate solution using leapfrog integration w/ exact forces
void leapfrog_leapq(int N, struct array_w_index *xi,double *v,double dt){
  int ii;
  for (ii=0; ii < N; ii++)
    (xi+ii)->val+= dt * *(v + (xi+ii)->idx);
}
void leapfrog_leappq(int N, struct array_w_index * xi,
		     double *v,double dt_kick,double dt_drift,
		     double *a){
  int ii;
  for (ii=0; ii< N; ii++) {
    *(v + (xi+ii)->idx)+= dt_kick * *(a + (xi+ii)->idx);
    (xi+ii)->val+= dt_drift * *(v + (xi+ii)->idx);
  }
}
int argsort_compare_function(const void *a,const void *b) {
  struct array_w_index *x = (struct array_w_index *) a;
  struct array_w_index *y = (struct array_w_index *) b;
  if (x->val < y->val) return -1;
  else if (x->val > y->val) return 1; 
  else return 0;
}
void _nbody_force(int N, struct array_w_index * xi, double * m, double * a,
		  double omega2, double * cumulmass,double * revcumulmass){
  int ii;
  // argsort
  qsort(xi,N,sizeof(struct array_w_index),argsort_compare_function);
  // Compute cumulative mass
  for (ii=0; ii< N-1; ii++)
    *(cumulmass+ii+1)= *(cumulmass+ii) + *(m+(xi+ii)->idx); 
  // Now compute acceleration from reverse-cumulative mass
  for (ii=0; ii< N-1; ii++)
    *(revcumulmass+N-ii-2)= *(revcumulmass+N-ii-1) + *(m+(xi+N-ii-1)->idx);
  if ( omega2 < 0 )
    for (ii=0; ii< N; ii++)
      *(a + (xi+ii)->idx)= *(revcumulmass+ii) - *(cumulmass+ii);
  else
    for (ii=0; ii< N; ii++)
      *(a + (xi+ii)->idx)= *(revcumulmass+ii) - *(cumulmass+ii) \
	- omega2 * (xi+ii)->val;
}
void _wendy_nbody_approx_onestep(int N, struct array_w_index * xi, 
				 double * x, double * v, 
				 double * m, double * a,
				 double dt, int nleap, double omega2,
				 int * err,double * time_elapsed,
				 double * cumulmass, double * revcumulmass){
  int ii;
  clock_t time_begin= clock();
  //drift half
  leapfrog_leapq(N,xi,v,dt/2.);
  //now drift full for a while
  for (ii=0; ii < (nleap-1); ii++){
    //kick+drift
    _nbody_force(N,xi,m,a,omega2,cumulmass,revcumulmass);
    leapfrog_leappq(N,xi,v,dt,dt,a);
  }
  //end with one last kick and drift
  _nbody_force(N,xi,m,a,omega2,cumulmass,revcumulmass);
  leapfrog_leappq(N,xi,v,dt,dt/2.,a);
  //de-sort
  for (ii=0; ii< N; ii++)
    *(x+(xi+ii)->idx)= (xi+ii)->val;
  clock_t time_end= clock();
  *time_elapsed= (double) (time_end-time_begin) / CLOCKS_PER_SEC;
}

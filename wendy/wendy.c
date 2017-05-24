/*
  wendy.c: One time-step of a one-dimensional N-body code
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <wendy.h>
#include <bst.h>
double _solve_quad_pos(double c0, double c1, double c2){
  // Solves for collisions under quadratic motion: a t^2/2 + vt + x
  double ca, mba, sqD, out;
  ca= c0/c2;
  mba= -c1/c2;
  if ( ca >= 0. && mba < 0. )
    return INFINITY;
  sqD= sqrt ( mba * mba - 4. * ca);
  if ( ca <= 0.){
    return 0.5 * ( mba + sqD );
  }
  out= 0.5 * (mba - sqD );
  if ( out < 1.E-16 )
    return 0.5 * ( mba + sqD );
  else
    return out;	     
}
double _solve_harm_pos(double c0, double c1, double c2, double omega){
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
			  double dt,
			  int maxcoll, int * err, int * ncoll){
  int cnt_coll,ii, tmpi;
  int c_in_x_indx, c_in_x_next_indx;
  double dv,tdt, dm, tmpd;
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
  while ( *next_tcoll < dt && cnt_coll < maxcoll ){
    //printf("Colliding in %f\n",*next_tcoll);
    //fflush(stdout);
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
    if ( tmpd <= 0. ) 
      tmpd= INFINITY;
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
	_solve_quad_pos(*(x + c_in_x_indx) + *(a+c_in_x_indx) * tdt * tdt / 2. 
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
	_solve_quad_pos(*(x + c_in_x_indx) + *(a+c_in_x_indx) * tdt * tdt / 2. 
			- *(v+c_in_x_indx) * tdt - *(x + c_in_x_next_indx),
			*(v + c_in_x_indx) - *(a + c_in_x_indx) * tdt
			- *(v + c_in_x_next_indx),
			0.5*(*(a + c_in_x_indx) - *(a + c_in_x_next_indx)));
      bst_tcoll= bst_forceInsert(bst_tcoll,*cindx+1,tcoll+*cindx+1);
    }
    // Find minimum
    minNode= bst_minValueNode(bst_tcoll);
    *cindx= minNode->idx;
    *next_tcoll= *minNode->val;
    //printf("Next one %f\n",*next_tcoll);
    //fflush(stdout);
  }
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
  free(t);
  bst_destroy(bst_tcoll);
  *ncoll= cnt_coll;
  if ( cnt_coll == maxcoll )
    *err= -2;
}

void _wendy_nbody_harm_onestep(int N, double * x, double * v, double * a,
			       double * m, int * sindx,
			       int * cindx, double * next_tcoll,
			       double * tcoll,double dt,
			       int maxcoll, int * err, int * ncoll,
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
	_solve_harm_pos(*(x + c_in_x_indx) - *(x + c_in_x_next_indx),
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
	_solve_harm_pos(*(x + c_in_x_indx) - *(x + c_in_x_next_indx),
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

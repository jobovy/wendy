/*
  wendy.c: One time-step of a one-dimensional N-body code
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <wendy.h>
#include <bst.h>
double _solve_quad_pos(double c0, double c1, double c2){
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
    // Update collision times
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
 

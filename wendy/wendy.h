/*
  C code for 1D N-body calculations
*/
#ifndef __WENDY_H__
#define __WENDY_H__
#ifndef CHUNK_PARALLEL_LEAPFROG
#define CHUNK_PARALLEL_LEAPFROG 1000
#endif
struct array_w_index
{
  int idx;
  double val;
};
double _solve_quad_pos(double,double,double);
double _solve_harm_pos(double,double,double,double);
void _wendy_nbody_onestep(int, double *, double *, double *,
			  double *,int *,
			  int *, double *, double *,
			  double,
			  int,int *,int *,double *);
void _wendy_nbody_harm_onestep(int, double *, double *, double *,
			       double *,int *,
			       int *, double *, double *,
			       double,
			       int,int *,int *,double *,double);
void _wendy_nbody_approx_onestep(int,struct array_w_index *, double *, 
				 double *,double *, double *,double,
				 double, int,double,int,int *,double *,
				 double *);
#endif /* wendy.h */

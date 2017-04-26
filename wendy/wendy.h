/*
  C code for 1D N-body calculations
*/
#ifndef __WENDY_H__
#define __WENDY_H__
double _solve_quad_pos(double,double,double);
void _wendy_nbody_onestep(int, double *, double *, double *,
			  double *,int *,
			  int *, double *, double *,
			  double,
			  int,int *,int *);
#endif /* wendy.h */

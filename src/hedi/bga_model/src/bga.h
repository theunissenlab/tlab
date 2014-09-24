#ifndef BGA_H
#define BGA_H

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_errno.h>
#include "bga.h"

typedef struct 
{
/* data */
	int nap; /* number of aalpha points */
	int nbp; /* number of beta point */
	int nfp; /* number of filter point */

	double * alpha; /* k values : interpolation */
	double * talpha; /* time points */ 
	gsl_spline * a_spline; /* spline interpolation of k*/
 	gsl_interp_accel *a_acc ; /* accelerator */

	double * beta;  /* b values */
	double * tbeta;
	gsl_spline * b_spline; /* spline interpolation of b*/
	gsl_interp_accel *b_acc ; /* accelerator */

	double * filter; /* filter valuer interpolation */
	double * freqs; /* frequency points */ 
	gsl_spline * f_spline; /* spline interpolation of k*/
 	gsl_interp_accel *f_acc ; /* accelerator */
 	
	double gamma;
	double x0;
	double y0; 
}  ODEparams;

int ODEfuncNS(double t, const double y[], double f[], void *params);
int ODEjacNS(double t, const double y[], double *dfdy, double dfdt[], void *params);
void ODErunNS(double *data, double tmax, double dt, double *pa);
int ODEfunc(double t, const double y[], double f[], void *params);
int ODEjac(double t, const double y[], double *dfdy, double dfdt[], void *params);


ODEparams * ODEinit(int nap, double *alpha, double *talpha, int nbp, double *beta, double *tbeta, int nfp, double *filter, double *freqs, double gamma, double x0, double y0);
void ODEupdate(ODEparams * pa, double *alpha, double *beta, double *filter);
void smBGAs(double * data, int size,double gamma, double alpha, double beta);
void smBGAsl(double * data, int size,double gamma, double * alpha, double * talpha, int nap,double * beta, double *tbeta, int nbp,double * filter, double * freq, int nfp);

#endif

#ifndef SM_H
#define SM_H

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <sndfile.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <gsl/gsl_multimin.h>


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

typedef struct 
{
	char fname[100];  /* wav filename */
	double samplerate;
	double cut_freq;
	double tmax;
	double dt;
	int nfft;
	int ntarget;
	int nslices;
	int target_ns;
	double * spectrum;
	double  * vdata;
	int nchunks;
	int n;


	ODEparams * ODEpa;

} OPTparams;


/* BGA ODE's without spline */

int ODEfuncNS(double t, const double y[], double f[],void *params);
int ODEjacNS(double t, const double y[], double *dfdy, double dfdt[], void *params);
void ODErunNS(double * data, double tmax, double dt,double * pa);
void ODEfilterNS(double * data, int nsample,int nfft,double * filter);

/* BGA ODE's with spline */

int ODEfunc (double t, const double y[], double f[],void *params);
int ODEjac (double t, const double y[], double *dfdy, double dfdt[], void *params);
void ODErun(double * data, double tmax, double dt,ODEparams * pa);

/* filter appication */
void ODEfilter(double * data, int nsample,int nfft,ODEparams * pa);

double getABSmax(double * data, int size);
void  scaleData(double scale,double * data, int size);
void save_wave_file(char * fname,double * data,int size);
double * load_wave_file(char * fname, int * ns);


void ODEtarget(double * x, double * data, int ns, int nfft, int ntarget, int nslices,int * chunks);
double * get_spectrum_data(double * data, int ns, int nfft, int ntarget, int nslices,int * chunks);

ODEparams * ODEinit(int nap,double * alpha,double * talpha, int nbp, double* beta, double *tbeta,int nfp,double * filter, double * freqs,double gamma, double x0,double y0);
void ODEupdate(ODEparams* pa,double * alpha,double* beta,double * filter);
OPTparams * OPTinit(char * fname,int nap,double * salpha,int nbp,double * sbeta, int nfp,double *sfilter);
double CHI2(double * a, double *b, int size);
void OPTresponse(double *p, double *x, int m, int n, void *data);
double OPTchi2(double *p, double *x, int m, int n, void *data);



#endif
# 
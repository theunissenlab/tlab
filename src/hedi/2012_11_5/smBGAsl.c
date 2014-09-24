#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>

#include "mex.h"


typedef struct 
{
/* data */
	int nap; /* number of aalpha points */
	int nbp; /* number of beta point */

	double * alpha; /* k values : interpolation */
	double * talpha; /* time points */ 
	gsl_spline * a_spline; /* spline interpolation of k*/
 	gsl_interp_accel *a_acc ; /* accelerator */

	double * beta;  /* b values */
	double * tbeta;
	gsl_spline * b_spline; /* spline interpolation of b*/
	gsl_interp_accel *b_acc ; /* accelerator */
 	
	double gamma;
	double x0;
	double y0; 
}  ODEparams;



int func (double t, const double y[], double f[],
      void *params)
{

  ODEparams * pa = (ODEparams *)params;
 
  double alpha=gsl_spline_eval(pa->a_spline,t, pa->a_acc);
  double beta=gsl_spline_eval(pa->b_spline,t, pa->b_acc);
  double gamma=pa->gamma,gamma2=gamma*gamma;
  double yy=y[1],x=y[0];
/*printf("%e %e %e\n",t,alpha,beta );*/
  f[0] = yy;
  f[1] = alpha*gamma2 - beta*gamma2*x-gamma2*x*x*x - gamma*x*x*yy+gamma2*x*x - gamma*x*yy;
  return GSL_SUCCESS;
}

int
jac (double t, const double y[], double *dfdy, 
     double dfdt[], void *params)
{

  ODEparams * pa = (ODEparams *)params;
  double alpha = gsl_spline_eval(pa->a_spline,t, pa->a_acc);
  double beta = gsl_spline_eval(pa->b_spline,t, pa->b_acc);
  double gamma=pa->gamma,gamma2=gamma*gamma;
  double yy=y[1],x=y[0];


  /* double * p = (double *)params;*/
  	gsl_matrix_view dfdy_mat = gsl_matrix_view_array (dfdy, 2, 2);
  	gsl_matrix * m = &dfdy_mat.matrix; 
	/*	0              1
		-k -2*c*y    -(p-b)-cx^2 
	*/
  	gsl_matrix_set (m, 0, 0, 0.0);
	gsl_matrix_set (m, 0, 1, 1.0);
  	gsl_matrix_set (m, 1, 0, -beta*gamma2 - 3*gamma2*x*x - 2*gamma*x*yy + 2*gamma2*x - gamma*yy);
  	gsl_matrix_set (m, 1, 1, -gamma*x*x - gamma*x);
	dfdt[0] = 0.0;
    dfdt[1] = 0.0;
  return GSL_SUCCESS;
}



ODEparams * init(int nap,double * alpha,double * talpha, int nbp, double* beta, double *tbeta,double gamma, double x0,double y0)
{
	ODEparams * pa= (ODEparams*)malloc(sizeof(ODEparams));

	pa->gamma=gamma;
	pa->x0=x0;
	pa->y0=y0;
	
	pa->nap=nap;
	pa->nbp=nbp;

	pa->alpha=calloc(nap,sizeof(double));
	pa->talpha=calloc(nap,sizeof(double));
	pa->beta=calloc(nbp,sizeof(double));
	pa->tbeta=calloc(nbp,sizeof(double));
	int i;
	memcpy(pa->alpha,alpha,nap*sizeof(double));
	memcpy(pa->talpha,talpha,nap*sizeof(double));
	memcpy(pa->beta,beta,nbp*sizeof(double));
	memcpy(pa->tbeta,tbeta,nbp*sizeof(double));
/*	int i;
/	for(i=0;i<nbp;i++) printf("%lg",pa->tb[i]);*/
    
	pa->a_acc=gsl_interp_accel_alloc ();
	pa->b_acc=gsl_interp_accel_alloc ();
	pa->a_spline = gsl_spline_alloc (gsl_interp_linear,nap);

	int err =	gsl_spline_init (pa->a_spline, pa->talpha, pa->alpha, pa->nap);
	if (err != GSL_SUCCESS)
		{
			mexErrMsgTxt("error: init vector interp - not good size / not increasing");
		}
	pa->b_spline = gsl_spline_alloc (gsl_interp_linear,nbp);
	err = gsl_spline_init (pa->b_spline, pa->tbeta, pa->beta, pa->nbp);
	if (err != GSL_SUCCESS)
		{
		mexErrMsgTxt("error: init vector interp - not good size / not increasing");
		}	
	pa->gamma=gamma;
	pa->x0=x0;
	pa->y0=y0;
	return pa;
}


void runOde(double * data, double tmax, double dt,ODEparams * pa)
{

	const gsl_odeiv_step_type * T = gsl_odeiv_step_rk8pd;
  	gsl_odeiv_step * s = gsl_odeiv_step_alloc (T, 2);
  	gsl_odeiv_control * c = gsl_odeiv_control_y_new (1e-10, 0.0);
  	gsl_odeiv_evolve * e  = gsl_odeiv_evolve_alloc (2);
  	gsl_odeiv_system sys = {func, jac, 2, pa};
  	int totalframe=ceil(tmax/dt);
  	int i=0;
	double t1=dt,t=0.0;  		
	double h = 1e-10;
	double y[2] = { pa->x0,pa->y0};
	

	while(t1<tmax)
    {
       while (t < t1)
  	    {  
	               int status = gsl_odeiv_evolve_apply (e, c, s,
					       &sys, 
					       &t, t1,
					       &h, y);
	  
        	  if (status != GSL_SUCCESS)
	               break;
      	}	  
	    i++;	
      	t1+=dt;
      	data[i]=y[0];
/*    	printf("%lf %lf\n",t1,y[0]); */

	}
	gsl_odeiv_evolve_free (e);
  	gsl_odeiv_control_free (c);
	gsl_odeiv_step_free (s);

}

void mexFunction(int nlhs,       mxArray *plhs[],
		 int nrhs, const mxArray *prhs[])
{

	int mrows,ncols;
	mrows = mxGetM(prhs[0]);
	ncols = mxGetN(prhs[0]);
  	if( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) ||
      !(mrows==1 && ncols==1) ) {
    mexErrMsgTxt("Input must be a noncomplex scalar double.");
  }
	double * ptmax=mxGetPr(prhs[0]);
	double * pgamma=mxGetPr(prhs[1]);



	int	marows = mxGetM(prhs[2]);
	int nacols = mxGetN(prhs[2]);
	if (nacols!=1)
		mexErrMsgTxt("a must be a one column vector");
	int nap=marows;

	int	mtarows = mxGetM(prhs[3]);
	int ntacols = mxGetN(prhs[3]);
	if (ntacols!=1)
		mexErrMsgTxt("ta must be a one column vector");
	if (mtarows != marows)
		mexErrMsgTxt("ta and k must be have the same length");			
	
	int	mbrows = mxGetM(prhs[4]);
	int nbcols = mxGetN(prhs[4]);
	if (nbcols!=1)
		mexErrMsgTxt("b must be a one column vector");
	int nbp=mbrows;
	

	int	mtbrows = mxGetM(prhs[5]);
	int ntbcols = mxGetN(prhs[5]);
	if (ntbcols!=1)
		mexErrMsgTxt("tb must be a one column vector");
	if (mtbrows != mbrows)
		mexErrMsgTxt("tb and b must be have the same length");			
	
	gsl_error_handler_t * err=gsl_set_error_handler_off ();

	double samplerate=44100.;

	double tmax=ptmax[0];/*0.200;*/
	double gamma=pgamma[0];

	double * alpha=mxGetPr(prhs[2]);	
	double * talpha=mxGetPr(prhs[3]);
	double * beta=mxGetPr(prhs[4]);
	double * tbeta=mxGetPr(prhs[5]);	

  	
	double dt=1/(samplerate);

	int totalframe=ceil(tmax/dt);

	ODEparams * pa=init(nap,alpha,talpha,nbp,beta,tbeta,gamma,0,1);


  	plhs[0] = mxCreateDoubleMatrix(1,totalframe, mxREAL);
  	double * data = mxGetPr(plhs[0]);
	runOde(data,tmax,dt,pa);
	
}
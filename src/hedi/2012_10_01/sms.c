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

	double k; /* k values : interpolation */
	double b;  /* b values */
	double p;
	double c;
	double x0;
	double y0; 
}  ODEparams;



int func (double t, const double y[], double f[],
      void *params)
{
  ODEparams * pa = (ODEparams *)params;
  double k=pa->k;
  double b=pa->b;
  double p=pa->p;
  double c=pa->c;

  f[0] = y[1];
  f[1] = -k*y[0] -(p-b)*y[1]- c*y[1]*y[0]*y[0];

  return GSL_SUCCESS;
}

int
jac (double t, const double y[], double *dfdy, 
     double dfdt[], void *params)
{

	ODEparams * pa = (ODEparams *)params;
	double k=pa->k;
	double b=pa->b;
  	double p=pa->p;
  	double c=pa->c;

  //double * p = (double *)params;
  	gsl_matrix_view dfdy_mat = gsl_matrix_view_array (dfdy, 2, 2);
  	gsl_matrix * m = &dfdy_mat.matrix; 
	/*	0              1
		-k -2*c*y    -(p-b)-cx^2 
	*/
  	gsl_matrix_set (m, 0, 0, 0.0);
	gsl_matrix_set (m, 0, 1, 1.0);
  	gsl_matrix_set (m, 1, 0, -2.0*c*y[0]*y[1] -k);
  	gsl_matrix_set (m, 1, 1, -(p-b)-c*y[0]*y[0]);
	dfdt[0] = 0.0;
    dfdt[1] = 0.0;
  return GSL_SUCCESS;
}



ODEparams * init(double k,double b, double c, double p,double x0, double y0)
{
	ODEparams * pa= (ODEparams*)malloc(sizeof(ODEparams));
	pa->c=c;
	pa->p=p;
	pa->k=k;
	pa->b=b;
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
	
	/* transients */ 

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
      	t1+=dt;     
//      	printf("%lf %lf\n",t1,y[0]);

	}
	t1=dt;
	t=0;
	i=0;
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
//      	printf("%lf %lf\n",t1,y[0]);

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
	double * pc=mxGetPr(prhs[1]);
	double * pp=mxGetPr(prhs[2]);
	double * pk=mxGetPr(prhs[3]);
	double * pb=mxGetPr(prhs[4]);
	double * px0=mxGetPr(prhs[5]);
	double * py0=mxGetPr(prhs[6]);

  	//	printf("%g %d %p %g\n",tmax,nrhs,ptmax,ptmax[0]);
	//tmax=ptmax[0];
	// c=*(mxGetPr(prhs[1]));
	// p=*(mxGetPr(prhs[2]));
	
 	double samplerate=44100.;
	double tmax=ptmax[0];//0.200;
	double c=pc[0]; //2000
	double p=pp[0];//13000;
	double k=pk[0];
	double b=pb[0];
	double x0=px0[0];
	double y0=py0[0];
	double dt=1/(samplerate);

	int totalframe=ceil(tmax/dt);
//	printf("%d %lg\n",totalframe,tmax);

	ODEparams * pa=init(k,b,c,p,x0,y0);


  	plhs[0] = mxCreateDoubleMatrix(1,totalframe, mxREAL);
  	double * data = mxGetPr(plhs[0]);
	runOde(data,tmax,dt,pa);
	
}
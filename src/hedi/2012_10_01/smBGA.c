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

	double alpha; /* k values : interpolation */
	double beta;  /* b values */
	double gamma;
	double x0;
	double y0; 
}  ODEparams;



int func (double t, const double y[], double f[],
      void *params)
{
  ODEparams * pa = (ODEparams *)params;
  double alpha=pa->alpha;
  double beta=pa->beta;
  double gamma=pa->gamma,gamma2=gamma*gamma;
  double yy=y[1],x=y[0];

  f[0] = yy;
  f[1] = alpha*gamma2 - beta*gamma2*x-gamma2*x*x*x - gamma*x*x*yy+gamma2*x*x - gamma*x*yy;
  return GSL_SUCCESS;
}

int
jac (double t, const double y[], double *dfdy, 
     double dfdt[], void *params)
{

  ODEparams * pa = (ODEparams *)params;
  double alpha=pa->alpha;
  double beta=pa->beta;
  double gamma=pa->gamma,gamma2=gamma*gamma;
  double yy=y[1],x=y[0];


  //double * p = (double *)params;
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



ODEparams * init(double alpha,double beta, double gamma, double x0,double y0)
{
	ODEparams * pa= (ODEparams*)malloc(sizeof(ODEparams));
	pa->alpha=alpha;
	pa->beta=beta;
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
	double * palpha=mxGetPr(prhs[1]);
	double * pbeta=mxGetPr(prhs[2]);
	double * pgamma=mxGetPr(prhs[3]);
	double * px0=mxGetPr(prhs[4]);
	double * py0=mxGetPr(prhs[5]);

  	//	printf("%g %d %p %g\n",tmax,nrhs,ptmax,ptmax[0]);
	//tmax=ptmax[0];
	// c=*(mxGetPr(prhs[1]));
	// p=*(mxGetPr(prhs[2]));
	
 	double samplerate=44100.;
	double tmax=ptmax[0];//0.200;
	double alpha=palpha[0]; //2000
	double beta=pbeta[0];//13000;
	double gamma=pgamma[0];
	double x0=px0[0];
	double y0=py0[0];
	double dt=1/(samplerate);

	int totalframe=ceil(tmax/dt);
//	printf("%d %lg\n",totalframe,tmax);

	ODEparams * pa=init(alpha,beta,gamma,x0,y0);


  	plhs[0] = mxCreateDoubleMatrix(1,totalframe, mxREAL);
  	double * data = mxGetPr(plhs[0]);
	runOde(data,tmax,dt,pa);
	
}
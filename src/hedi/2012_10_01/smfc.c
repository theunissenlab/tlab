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
	int nkp; /* number of k points */
	int nbp; /* number of p point */

	double * k; /* k values : interpolation */
	double * tk; /* time points */ 
	gsl_spline * k_spline; /* spline interpolation of k*/
 	gsl_interp_accel *k_acc ; /* accelerator */

	double * b;  /* b values */
	double * tb;
	gsl_spline * b_spline; /* spline interpolation of b*/
	gsl_interp_accel *b_acc ; /* accelerator */
 	
	double * p; /* k values : interpolation */
	double * tp; /* time points */ 
	gsl_spline * p_spline; /* spline interpolation of k*/
 	gsl_interp_accel *p_acc ; /* accelerator */

	double c; 
}  ODEparams;



int func (double t, const double y[], double f[],
      void *params)
{
  ODEparams * pa = (ODEparams *)params;
  double k=gsl_spline_eval(pa->k_spline,t, pa->k_acc);
  double b=gsl_spline_eval(pa->b_spline,t, pa->b_acc);
  double p=pa->p;
  double c=pa->c;
//  printf("[%g k=%g b=%g p=%g]\n",t,k,b,p);	
  f[0] = y[1];
  f[1] = -k*y[0] -(p-b)*y[1]- c*y[1]*y[0]*y[0];

  return GSL_SUCCESS;
}

int
jac (double t, const double y[], double *dfdy, 
     double dfdt[], void *params)
{

	ODEparams * pa = (ODEparams *)params;
	double k=gsl_spline_eval(pa->k_spline,t, pa->k_acc);
 	double b=gsl_spline_eval(pa->b_spline,t, pa->b_acc);
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



ODEparams * init(int nkp, double * k, double * kt,int nbp, double *b, double *bt, double c, double p)
{
	ODEparams * pa= (ODEparams*)malloc(sizeof(ODEparams));
	pa->c=c;
	pa->p=p;
	pa->nkp=nkp;
	pa->nbp=nbp;
	pa->k=calloc(nkp,sizeof(double));
	pa->tk=calloc(nkp,sizeof(double));
	pa->b=calloc(nbp,sizeof(double));
	pa->tb=calloc(nbp,sizeof(double));
	memcpy(pa->k,k,nkp*sizeof(double));
	memcpy(pa->tk,kt,nkp*sizeof(double));
	memcpy(pa->b,b,nbp*sizeof(double));
	memcpy(pa->tb,bt,nbp*sizeof(double));
//	int i;
//	for(i=0;i<nbp;i++) printf("%lg",pa->tb[i]);
	pa->k_acc=gsl_interp_accel_alloc ();
	pa->b_acc=gsl_interp_accel_alloc ();
	pa->k_spline = gsl_spline_alloc (gsl_interp_cspline,nkp);
	gsl_spline_init (pa->k_spline, pa->tk, pa->k, pa->nkp);
	pa->b_spline = gsl_spline_alloc (gsl_interp_cspline,nbp);
	gsl_spline_init (pa->b_spline, pa->tb, pa->b, pa->nbp);
	return pa;
}


void update(ODEparams * pa, double * k, double * kt, double *b, double *bt, double c, double p)
{
	pa->c=c;
	pa->p=p;
	memcpy(pa->k,k,pa->nkp*sizeof(double));
	memcpy(pa->tk,kt,pa->nkp*sizeof(double));
	memcpy(pa->b,b,pa->nbp*sizeof(double));
	memcpy(pa->tb,bt,pa->nbp*sizeof(double));	
	gsl_spline_init (pa->k_spline, pa->tk, pa->k, pa->nkp);
	gsl_spline_init (pa->b_spline, pa->tb, pa->b, pa->nbp);	
}

void cleanup(ODEparams * pa)
{
	gsl_spline_free (pa->k_spline);
	gsl_spline_free (pa->b_spline);
    gsl_interp_accel_free (pa->k_acc);
    gsl_interp_accel_free (pa->b_acc);
    free(pa->k);free(pa->tk);free(pa->b);free(pa->tb);
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
	double y[2] = { 1.00, 1.0 };

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
	double * pb=mxGetPr(prhs[2]);
	
	int	mkrows = mxGetM(prhs[3]);
	int nkcols = mxGetN(prhs[3]);
	if (nkcols!=1)
		mexErrMsgTxt("k must be a one column vector");
	int nkp=mkrows;

	int	mtkrows = mxGetM(prhs[4]);
	int ntkcols = mxGetN(prhs[4]);
	if (ntkcols!=1)
		mexErrMsgTxt("tk must be a one column vector");
	if (mtkrows != mkrows)
		mexErrMsgTxt("tk and k must be have the same length");			


	int	mbrows = mxGetM(prhs[5]);
	int nbcols = mxGetN(prhs[5]);
	if (nbcols!=1)
		mexErrMsgTxt("b must be a one column vector");
	int nbp=mbrows;

	int	mtbrows = mxGetM(prhs[6]);
	int ntbcols = mxGetN(prhs[6]);
	if (ntbcols!=1)
		mexErrMsgTxt("tb must be a one column vector");
	if (mtbrows != mbrows)
		mexErrMsgTxt("tb and b must be have the same length");			
	
	double * k=mxGetPr(prhs[3]);	
	double * kt=mxGetPr(prhs[4]);
	double * b=mxGetPr(prhs[5]);
	double * bt=mxGetPr(prhs[6]);	


  	//	printf("%g %d %p %g\n",tmax,nrhs,ptmax,ptmax[0]);
	//tmax=ptmax[0];
	// c=*(mxGetPr(prhs[1]));
	// p=*(mxGetPr(prhs[2]));
	
 	double samplerate=44100.;
	double tmax=ptmax[0];//0.200;
	double c=pc[0]; //2000
	double p=pb[0];//13000;
	double dt=1/(samplerate);

	/* final time check */
	if ((tmax>kt[mkrows-1])||(tmax>bt[mbrows-1]))
		mexErrMsgTxt("tmax must be < both tb and tk max times");		

	// double	kt[]={0, 0.05, 0.10 ,0.125,0.15,0.25};
	// double 	k[]= {4e7, 4e7  , 4e7 ,4e7  ,4e7    ,4e7 };
	// double	bt[]={0, 0.05, 0.10 ,0.125,0.15,0.25};
	// double 	b[]= {4e4, 4e4, 4e4 ,4e4  ,4e4    ,4e4 };
	
	


	int totalframe=ceil(tmax/dt);
	printf("%d %lg\n",totalframe,tmax);

	ODEparams * pa=init(nkp,k,kt,nbp,b,bt,c,p);


  	plhs[0] = mxCreateDoubleMatrix(1,totalframe, mxREAL);
  	double * data = mxGetPr(plhs[0]);
	runOde(data,tmax,dt,pa);
	cleanup(pa);
}
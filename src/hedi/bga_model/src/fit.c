#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <gsl/gsl_multimin.h>

#include "spectrum.h"
#include "optim.h"
#include "fit.h"


double minFun(const gsl_vector * v, void *params)
{
	FIT * fit=(FIT*)params;
	int nfilter=fit->nfilter;
	double fmin=fit->fmin;
	double fmax=fit->fmax;
	SPEC_DATA * spec=fit->target_spec;

	int frame=fit->current_frame;
	double alpha = gsl_vector_get(v,0);


	double beta = gsl_vector_get(v,1);


	double x[nfilter+2];
	double y[nfilter+2];

	int i;
//	for(i=0;i<nfilter;i++) {x[i]=gsl_vector_get(v,2+i); y[i]=gsl_vector_get(v,2+i+nfilter);}

	 memcpy(x, gsl_vector_const_ptr (v, 2),nfilter*sizeof(double));
//	 printf("%d\n",__LINE__);
	 memcpy(y, gsl_vector_const_ptr (v, 2+nfilter),nfilter*sizeof(double));
	x[nfilter]=-1000;
	x[nfilter+1]=10000;
	y[nfilter]=0;
	y[nfilter+1]=0;

	//for(i=0;i<nfilter;i++) printf("%lf %lf\n",x[i],y[i]);
	return spec_distance(alpha,beta, x,y,nfilter+2,frame, fmin,fmax, spec,0);	
}


double  fitOneFrame(FIT * fit,int frame, double * guess)
{


	int nfilter=fit->nfilter;

	double fmax=fit->fmax;
	double fmin=fit->fmin;
	fit->current_frame=frame;
	double alpha = guess[0];
	double beta = guess[1];
	double x[nfilter];
	double y[nfilter];
	memcpy(x,guess+2,nfilter*sizeof(double));
	memcpy(y,guess+2+nfilter,nfilter*sizeof(double));

 	const gsl_multimin_fminimizer_type *T =
	gsl_multimin_fminimizer_nmsimplex2;
	gsl_multimin_fminimizer *s = NULL;	

	gsl_multimin_function minex_func;
	size_t          iter = 0;
	int             status;
	double          size;

	int m = 2*nfilter+2;
	int i;

	gsl_vector  * start = gsl_vector_alloc(m);
	gsl_vector     *ss = gsl_vector_alloc(m);
	for (i = 0; i < m; i++)
 		gsl_vector_set(ss, i, 0.05);
 	for (i = 2; i < 2+nfilter; i++)
 		gsl_vector_set(ss, i, 1e3);	

 	for (i = 0; i < m; i++)
 		gsl_vector_set(start,i,guess[i]);

 	minex_func.n = m;
	minex_func.f = minFun;
	minex_func.params = fit;
	s = gsl_multimin_fminimizer_alloc(T, m);
	gsl_multimin_fminimizer_set(s, &minex_func, start, ss);

	do {
		iter++;
//		printf("%d\n",iter);
		status = gsl_multimin_fminimizer_iterate(s);

		if (status)
			break;

		size = gsl_multimin_fminimizer_size(s);
		status = gsl_multimin_test_size(size, 1e-2);

		// if (status == GSL_SUCCESS) {
		// 	printf("converged to minimum at\n");
		// }

		// for(i=0;i<m;i++) 
		// 	printf("%e ",gsl_vector_get(s->x,i));
		// printf("\n");
	}
	while (status == GSL_CONTINUE && iter < 500);
	
	double d=minFun(s->x,fit);
	memcpy(fit->best_current_frame, gsl_vector_const_ptr (s->x, 0),m*sizeof(double));	
	gsl_vector_free(start);
	gsl_vector_free(ss);
	gsl_multimin_fminimizer_free(s);
	return d;
}

double minFunv2(const gsl_vector * v, void *params)
{

	FIT * fit=(FIT*)params;

	int nfilter=fit->nfilter;
	double fmin=fit->fmin;
	double fmax=fit->fmax;

	SPEC_DATA * spec=fit->target_spec;

	int frame=fit->current_frame;
	double alpha = gsl_vector_get(v,0);
	double beta=fit->slope*alpha + fit->intercept;
	double x[nfilter+2];
	double y[nfilter+2];

	int i;
//	for(i=0;i<nfilter;i++) {x[i]=gsl_vector_get(v,2+i); y[i]=gsl_vector_get(v,2+i+nfilter);}

	 memcpy(x, gsl_vector_const_ptr (v, 1),nfilter*sizeof(double));
//	 printf("%d\n",__LINE__);
	 memcpy(y, gsl_vector_const_ptr (v, 1+nfilter),nfilter*sizeof(double));
	x[nfilter]=-1000;
	x[nfilter+1]=10000;
	y[nfilter]=0;
	y[nfilter+1]=0;

	//for(i=0;i<nfilter;i++) printf("%lf %lf\n",x[i],y[i]);
	return spec_distance(alpha,beta, x,y,nfilter+2,frame, fmin,fmax, spec,0);	
}


double  fitOneFramev2(FIT * fit,int frame, double slope, double intercept,double * guess)
{



	int nfilter=fit->nfilter;
	fit->intercept=intercept;
	fit->slope=slope;
	double fmax=fit->fmax;
	double fmin=fit->fmin;
	fit->current_frame=frame;
	double alpha = guess[0];
	double beta= alpha*slope+intercept;
	double x[nfilter];
	double y[nfilter];

	memcpy(x,guess+1,nfilter*sizeof(double));
	memcpy(y,guess+1+nfilter,nfilter*sizeof(double));

 	const gsl_multimin_fminimizer_type *T =
	gsl_multimin_fminimizer_nmsimplex2;
	gsl_multimin_fminimizer *s = NULL;	

	gsl_multimin_function minex_func;
	size_t          iter = 0;
	int             status;
	double          size;

	int m = 2*nfilter+1;
	int i;

	gsl_vector  * start = gsl_vector_alloc(m);
	gsl_vector     *ss = gsl_vector_alloc(m);
	for (i = 0; i < m; i++)
 		gsl_vector_set(ss, i, 0.05);
 	for (i = 1; i < 1+nfilter; i++)
 		gsl_vector_set(ss, i, 1e3);	

 	for (i = 0; i < m; i++)
 		gsl_vector_set(start,i,guess[i]);

 	minex_func.n = m;
	minex_func.f = minFunv2;
	minex_func.params = fit;
	s = gsl_multimin_fminimizer_alloc(T, m);
	gsl_multimin_fminimizer_set(s, &minex_func, start, ss);

	do {
		iter++;
		//printf("%d\n",iter);
		status = gsl_multimin_fminimizer_iterate(s);

		if (status)
			break;

		size = gsl_multimin_fminimizer_size(s);
		status = gsl_multimin_test_size(size, 1e-2);

		// if (status == GSL_SUCCESS) {
		// 	printf("converged to minimum at\n");
		// }
		double d=minFun(s->x,fit);
		printf("%lg ",d);
		for(i=0;i<m;i++) 
			printf("%e ",gsl_vector_get(s->x,i));
		printf("\n");
	}
	while (status == GSL_CONTINUE && iter < 500);
	
	double d=minFun(s->x,fit);
	memcpy(fit->best_current_frame, gsl_vector_const_ptr (s->x, 0),m*sizeof(double));	
	gsl_vector_free(start);
	gsl_vector_free(ss);
	gsl_multimin_fminimizer_free(s);
	return d;
}
// int
// main(int argc, char **argv)
// {

// 	/* multimin stuff */
// 	const gsl_multimin_fminimizer_type *T =
// 	gsl_multimin_fminimizer_nmsimplex2;
// 	gsl_multimin_fminimizer *s = NULL;
// 	gsl_vector     *ss, *x;
// 	gsl_multimin_function minex_func;
// 	size_t          iter = 0;
// 	int             status;
// 	double          size;

// 	/* loading and setting up */
// 	int             i, k;
// 	int             nap = 11;
// 	int             nbp = 4;
// 	int             nfp = 10;
// 	const int       m = nap + nbp + nfp;

// 	/* Starting point */
// 	x = gsl_vector_alloc(m);


	
// 	double          alpha[] = {-0.26, -0.1374, -0.1542, -0.1552, -0.1522, -0.1507, -0.1486, -0.1447, -0.1293, -0.1291, -0.1187};
// 	double          beta[] = {0.0105, 0.0105, 0.0105, 0.0105, 0.0105, 0.0105, 0.0105, 0.0105, 0.0105, 0.0105};
// 	double          filter[] = {1e-2, 1e-1, 1, 2, 4, 8, 1, 1.0, 0.0, 0.0};

// 	for (i = 0; i < nap; i++)
// 		gsl_vector_set(x, i, alpha[i]);
// 	for (i = 0; i < nbp; i++)
// 		gsl_vector_set(x, i + nap, beta[i]);
// 	for (i = 0; i < nfp; i++)
// 		gsl_vector_set(x, i + nap + nbp, filter[i]);

	
// 	 load file and params 

// 	OPTparams      *opt_pa = OPTinit(argv[1], nap, alpha, nbp, beta, nfp, filter);
// 	const int       n = opt_pa->n;

// 	/* Set initial step sizes */
// 	ss = gsl_vector_alloc(m);
// 	for (i = 0; i < nap; i++)
// 		gsl_vector_set(ss, i, 1e-1);
// 	for (i = 0; i < nbp; i++)
// 		gsl_vector_set(ss, i + nap, 1e-1);
// 	for (i = 0; i < nfp; i++)
// 		gsl_vector_set(ss, i + nap + nbp, 1);

// 	/* Initialize method and iterate */
// 	minex_func.n = m;
// 	minex_func.f = minFun;
// 	minex_func.params = opt_pa;

// 	s = gsl_multimin_fminimizer_alloc(T, m);
// 	gsl_multimin_fminimizer_set(s, &minex_func, x, ss);

// 	do {
// 		iter++;
// 		status = gsl_multimin_fminimizer_iterate(s);

// 		if (status)
// 			break;

// 		size = gsl_multimin_fminimizer_size(s);
// 		status = gsl_multimin_test_size(size, 1e-2);

// 		if (status == GSL_SUCCESS) {
// 			printf("converged to minimum at\n");
// 		}
// 		for(i=0;i<m;i++) 
// 			printf("%e ",gsl_vector_get(s->x,i));
// 		printf("\n");
// 	}
// 	while (status == GSL_CONTINUE && iter < 100);

// 	gsl_vector_free(x);
// 	gsl_vector_free(ss);
// 	gsl_multimin_fminimizer_free(s);


// }




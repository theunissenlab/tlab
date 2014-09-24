#include <stdio.h>


#include "sm.h"
#include "wave.h"


typedef struct 
{
	int NFREQ;
	int nfft;
	double samplerate;
	gsl_interp_accel *acc;
  	gsl_spline *spline;
  	double *z;
	/* data */
} param_set;


double minFun(const gsl_vector * v, void *params)
{
	param_set * pa=(param_set*)params;

	int k;
	int nfft=pa->nfft;
	int NFREQ=pa->NFREQ;
	double alpha=gsl_vector_get(v, 0);
	double beta=gsl_vector_get(v, 1);
	double *x=pa->z;
	double samplerate=pa->samplerate;
	
	double sfilter[NFREQ+1];
	double freqs[NFREQ+1];
	for (k = 0; k < NFREQ; k++)
	{
		sfilter[k] = gsl_vector_get(v, k+2);
		freqs[k]=8000./NFREQ*k;
	}		
	sfilter[NFREQ]=0;
	freqs[NFREQ]=samplerate/2;
	gsl_spline_init(pa->spline,freqs,sfilter,NFREQ+1);
	
	double filter[nfft];
	double z[nfft],y[2*nfft];
	for(k=0;k<nfft;k++) 
	{
		double ff=k/(1.0*nfft)*samplerate/2;
		filter[k]=gsl_spline_eval(pa->spline,ff,pa->acc);
	}

	
	double pab[]={alpha,beta,24000,1,0};	
	double tmax=2*nfft/samplerate;
	double dt=1/samplerate;

 	ODErunNS(y,tmax,dt,pab);	
	SLIfilter(y,nfft,filter);`
	
	
	SLItarget(z,y,nfft);

	double s=0;
	for(k=0;k<nfft;k++)
		s += (x[k]-z[k])*(x[k]-z[k]);
	return s;
}


int
main(int argc, char **argv)
{

	char fname[100];
	int nfft=512;
	
	int NFREQ=20;
	int FREQM=8000;
	int N=NFREQ+1;

	double samplerate=44100;
	
	double W=3;


	int i,k;
	double wav[nfft],awav[2*nfft],fawav[nfft],sfilter[nfft];

	const gsl_multimin_fminimizer_type *T =
	gsl_multimin_fminimizer_nmsimplex2;
	gsl_multimin_fminimizer *s = NULL;
	

	gsl_multimin_function minex_func;
	size_t          iter = 0;
	int             status;
	double          size;
	
	sprintf(fname,"%s",argv[1]);
	wave_array * wa = create_array(fname,nfft);
	param_set * params=(param_set*)malloc(sizeof(param_set));
	params->nfft=nfft;
	params->NFREQ=NFREQ;
	params->samplerate=samplerate;
	int nslices=wa->nslices;

	for(k=0;k<NFREQ;k++)
		sfilter[k]=exp(-(k-NFREQ/2.0)*(k-NFREQ/2.0)/W);
	sfilter[NFREQ]=0;

	params->acc = gsl_interp_accel_alloc ();
  	const gsl_interp_type *t = gsl_interp_linear; 
  	params->spline = gsl_spline_alloc (t, N);

	
	// double          sfilter[] = {1e-2, 1e-1, 1, 2, 4, 8, 1, 1.0, 0.0, 0.0};
	// double          freqs[] = {0, 1000, 2000,3000,4000,5000, 6000,7000,8000,samplerate/2};
	
	double filter[nfft];
	int m=NFREQ+2;
	gsl_vector * x=gsl_vector_alloc(NFREQ+2);
	gsl_vector * ss = gsl_vector_alloc(NFREQ+2);

	gsl_vector_set(ss, 0, 1e-1);
	gsl_vector_set(ss, 1, 1e-2);
	for (i = 0; i < NFREQ; i++)
		gsl_vector_set(ss, i + 2, 1);


	double best_alpha[nslices];
	double best_beta[nslices];
	double best_filter[nfft];
	for (k = 0; k < nfft; ++k)
	{
		best_filter[k]=0.;
	}


	for(i=0;i<wa->nslices;i++)
	{
		double alpha0=-0.15;
		double beta0 = 0.0105;

		
		
		gsl_vector_set(x,0,alpha0);
		gsl_vector_set(x,1,beta0);
		
		for (k = 0; k < NFREQ; k++)	
			gsl_vector_set(x, k+2, sfilter[k]);
		
		SLItarget(wav,wa->slices[i]->data,wa->nfft);
		params->z=wav;

		minex_func.n = NFREQ+2;
		minex_func.f = minFun;
		minex_func.params = params;

		s = gsl_multimin_fminimizer_alloc(T, m);
		gsl_multimin_fminimizer_set(s, &minex_func, x, ss);

		iter=0;
		status=GSL_CONTINUE;

		do {
			 iter++;
			status = gsl_multimin_fminimizer_iterate(s);			
			size = gsl_multimin_fminimizer_size(s);
			status = gsl_multimin_test_size(size, 1e-2);

		}
		while (status == GSL_CONTINUE && iter < 1000);
	
		//  printf("%d %e\n",i,minFun(s->x,params));	

		//  for(k=0;k<NFREQ+2;k++) 
		//  		printf("%e ",gsl_vector_get(s->x,k));
		//  printf("\n");
		 best_alpha[i]=gsl_vector_get(s->x,0);		
		best_beta[i]=gsl_vector_get(s->x,1);
		for(k=0;k<NFREQ;k++)
			best_filter[k]+= gsl_vector_get(s->x,k+2)/nslices;
	}	
	for(k=0;k<wa->nslices;k++)
		printf("%lf %lf\n",best_alpha[k],best_beta[k]);

}




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


double minFun(double alpha, double beta, double * filter,double * freqs, void *params)
{
	param_set * pa=(param_set*)params;

	int k;
	int nfft=pa->nfft;
	int NFREQ=pa->NFREQ;
	double *x=pa->z;
	double samplerate=pa->samplerate;

	double z[nfft],y[2*nfft];
	
	double pab[]={alpha,beta,24000,1,0};	
	double tmax=2*nfft/samplerate;
	double dt=1/samplerate;

 	ODErunNS(y,tmax,dt,pab);	
	SLIfilter(y,nfft,filter);
	y[0]=0;
	x[0]=0;
	double m =getABSmax(y,2*nfft);
	scaleData(1/m,y,2*nfft);

	SLItarget(z,y,nfft);
	m=getABSmax(z,nfft);
	double mm=getABSmax(x,nfft);
	scaleData(mm/m,z,nfft);

	z[0]=0;
	double s=0;
	for(k=0;k<nfft/2;k++)
	{
		s += (x[k]-z[k])*(x[k]-z[k]);
//		printf("%lf %lf %lf \n",x[k]*x[k],z[k]*z[k],y[k]);
	}
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

	
	sprintf(fname,"%s",argv[1]);
	wave_array * wa = create_array(fname,nfft);
	param_set * params=(param_set*)malloc(sizeof(param_set));
	params->nfft=nfft;
	params->NFREQ=NFREQ;
	params->samplerate=samplerate;

	int nslices=wa->nslices;
	params->acc = gsl_interp_accel_alloc ();
  	const gsl_interp_type *t = gsl_interp_linear; 
  	params->spline = gsl_spline_alloc (t, N);

	double best_alpha[nslices];
	double best_beta[nslices];

	for(i=0;i<wa->nslices;i++)
	{

		double sfilter[NFREQ+1];
		double d=1e6;
		double alpha,beta,balpha,bbeta;
		SLItarget(wav,wa->slices[i]->data,wa->nfft);

		params->z=wav;
	//	printf("%d\n",i);
		for(k=0;k<NFREQ;k++)
		{
			double fmax=8000;
			int index=floor((k*fmax/NFREQ)/(samplerate/2)*nfft);
		//	printf("%d %d %e\n",k,index,wav[index]*wav[index]);
			sfilter[k]=wav[index]*wav[index];
		}
		double freqs[NFREQ+1];
		for (k = 0; k < NFREQ; k++)
		{
			freqs[k]=8000./NFREQ*k;
		}		
		sfilter[NFREQ]=0;
		freqs[NFREQ]=samplerate/2;
		gsl_spline_init(params->spline,freqs,sfilter,NFREQ+1);
	
		double filter[nfft];
	
		for(k=0;k<nfft;k++) 
		{
			double ff=k/(1.0*nfft)*samplerate/2;
			filter[k]=gsl_spline_eval(params->spline,ff,params->acc);
		//	printf(" ** %e\n", filter[k]);
		}
		for(alpha=-0.2;alpha<-0.05;alpha+=0.001)
		 for(beta=0.00;beta<0.05;beta +=0.0005)
			{
			// alpha=-0.15;
			// beta=0.049;

				double D=minFun(alpha,beta,filter,freqs,params);
//				printf("%lf\n",D);
				 if (D<d)
				 {
				 	d=D;
					balpha=alpha;
					bbeta=beta;
				//	printf("%d %e %e %e\n",i,d,balpha,bbeta);
				}
			}			
		best_beta[i]=bbeta;
		best_alpha[i]=balpha;
		printf("%d %e %e\n",i,best_alpha[i],best_beta[i]);

	}
	// for(i=0;i<nslices;i++)
	// 	printf("%e %e\n",best_alpha[i],best_beta[i]);

}




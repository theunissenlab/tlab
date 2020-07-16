#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_spline.h>
#include "utils.h"
#include "optim.h"

gsl_spline * filter_spline; /* spline interpolation of k*/
gsl_interp_accel *filter_acc ; /* accelerator */
 	

int _compare(const void *a,const void *b)
{
	param_sort *A=(param_sort*)a;
	param_sort *B=(param_sort*)b;

	if (A->x>B->x)
		return 1;
	else return -1;
}

void interp_alloc(int nbfilter)
{
	filter_acc = gsl_interp_accel_alloc();
	filter_spline = gsl_spline_alloc(gsl_interp_linear, nbfilter);
}	

void _interp_init(int nparam,double *x, double *y)
{
	param_sort S[nparam];
	int i;
	for(i=0;i<nparam;i++)
	{
		S[i].y=y[i];
		S[i].x=x[i];
	}
	qsort(S,nparam,sizeof(param_sort),_compare);
	for(i=0;i<nparam;i++)
	{
		x[i]=S[i].x;
		y[i]=S[i].y;
//		printf("%lg %lg\n", x[i],y[i]);
	}	
	gsl_spline_init(filter_spline,x,y,nparam);
}

void get_vocal(double alpha,double beta,SPEC_DATA * fit_spec)
{
	double pa[]={alpha,beta,24000};
	double tmax=0.5; 
	double dt=1/fit_spec->sampleRate;
	int n=ceil(tmax/dt);
	double * data = (double *) calloc(n,sizeof(double));
	ODErunNS(data,tmax,dt,pa);
//	printf("%d %d %d %lf %lf\n",__LINE__,n,1,fit_spec->f,fit_spec->fband );
	TimeFreq(data, n,fit_spec->sampleRate,fit_spec->f,fit_spec->fband,fit_spec);
}


double * get_average_spectrum(SPEC_DATA*fit_spec,double fmax,double fmin)
{

	int frame;
	int tspecsize;
	frame=fit_spec->frameCount/2;
	double * output = GetSpec(fit_spec,frame, fmax,fmin,&tspecsize);
	for(frame=fit_spec->frameCount/2;frame<fit_spec->frameCount-1;frame++)
	{
		double * sp=GetSpec(fit_spec,frame, fmax,fmin,&tspecsize);
		int j;
		for(j=0;j<tspecsize;j++) 
			output[j]+=sp[j]/(fit_spec->frameCount/2)	;
		free(sp);
	}
	return output;

}
double spec_distance(double alpha, double beta, double *x, double *y, int nbfilter, int frame, double fmin, double fmax, SPEC_DATA * spec, int verbose)
{

	/* load target spectrum at the correct frame */
	int spsize;
	double * target_spectrum = GetSpec(spec,frame,fmax,fmin,&spsize);
	normalize(target_spectrum,spsize);
	double sampleRate = spec->sampleRate;
	/* compute vocal chord spectrum */
	SPEC_DATA  fit_spec;
	fit_spec.sampleRate=sampleRate;
	fit_spec.f=spec->f;
	fit_spec.fband=spec->fband;

//	printf("%d %lg\n",__LINE__,sampleRate);
	get_vocal(alpha,beta,&fit_spec);
//	printf("%d %lg\n",__LINE__,sampleRate);
	double * vocal_spectrum=get_average_spectrum(&fit_spec,fmax,fmin);
	free(fit_spec.data);
	/* compute interpolation parameter*/
	_interp_init(nbfilter,x,y);

	/* filter the spectrum */
	int i;
	for(i=0;i<spsize;i++)
	{
		double freq= i * sampleRate / spec->winLength;
		
		vocal_spectrum[i] *= fabs(gsl_spline_eval(filter_spline, freq, filter_acc));			
		
	}
	 normalize(vocal_spectrum,spsize);


	 /* compute the chi2 */
	double c= CHI2(vocal_spectrum,target_spectrum,spsize);
	for(i=0;i<spsize;i++)
	{
			double freq= i * sampleRate / spec->winLength;
	if (verbose)	printf("%lf %lf %lf %lf \n",freq,target_spectrum[i],vocal_spectrum[i],fabs(gsl_spline_eval(filter_spline, freq, filter_acc)));
	}
	//fprintf(stderr,"## %lg\n",c);
	free(target_spectrum);
	free(vocal_spectrum);
	
	return c;
}
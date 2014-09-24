#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <sndfile.h>
#include "src/bga.h"
#include "src/spectrum.h"
#include "src/utils.h"
#include "src/optim.h"
#include "src/fit.h"



int main(int argc,char ** argv)
{
		int i,k;
	char name[100];
	sprintf(name,"%s",argv[1]);
	int ns;
	double         * input=load_wave_file(name,&ns);
	SPEC_DATA spec;
	double f=1000; 
	double sampleRate=44100.;
	double fband=125;

	spec.f=f;
	spec.fband=fband;

	
	double alpha=atof(argv[2]);
	double beta=atof(argv[3]);
	double f1=atof(argv[4]);
	double f2=atof(argv[5]);
	double f3=atof(argv[6]);
	double a1=atof(argv[7]);
	double a2=atof(argv[8]);
	double a3=atof(argv[9]);

	/* spectrum of the target file */


	TimeFreq(input,ns,sampleRate,f, fband,&spec);


	/* starting frame */
	int frameCount=spec.frameCount;
	int frame=frameCount/2;

	int nfilter=3;
	int m=2+2*nfilter;
	double fmax=8000;
	double fmin=200;
	interp_alloc(nfilter+2);
	double best[m];

	FIT fit_data;	
	fit_data.target_spec=&spec;
	fit_data.nfilter=nfilter;
	
	fit_data.fmin=fmin;
	fit_data.fmax=fmax;
	fit_data.best_current_frame=best;

	double verybest[m];
	double guess[]={alpha,beta,f1,f2,f3,a1,a2,a3};
	memcpy(best,guess,m*sizeof(double));
	memcpy(verybest,guess,m*sizeof(double));
	
	double fit_values[(m+2)*frameCount];
	for(i=0;i<(m+2)*frameCount;i++) 
		fit_values[i]=0.0;

	for(frame=frameCount/2;frame<frameCount;frame+=10)
	{


		double d = fitOneFrame(&fit_data,frame,best);
		printf("%d %lg ",frame,d);
		for(i=0;i<m;i++)
			printf("%lg ",fit_data.best_current_frame[i]);
		printf("\n");		
		memcpy(best,fit_data.best_current_frame,m*sizeof(double));	
	}
	
	 
}


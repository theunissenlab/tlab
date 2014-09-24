#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "src/bga.h"
#include "src/spectrum.h"
#include "src/optim.h"
#include "src/utils.h"


void find_max_divisors(double * tab, int size, int nb,double * vmax, int * emax)
{
	
	int i,d;
	for(d=0;d<nb;d++)
	{
		int imax=-1;
		double m=-1;
		int offset=size/(1.0*nb)*d;
		for(i=0;i<size/(1.0*nb);i++)
		{
	
			if (fabs(tab[i+offset])>m)
			{
				m=tab[i+offset];
				imax=i+offset;
			}
		}
		emax[d]=imax;
		vmax[d]=m;
	}
	
}

int main(int argc,char ** argv)
{
		int i,k;
	char name[100];
	int ns;
	SPEC_DATA spec;
	double f=1000; 
	double sampleRate=44100.;
	double fband=25;
	double fmin=200;
	double fmax=8000;
	
	spec.f=f;
	spec.fband=fband;
	spec.sampleRate=sampleRate;


	sprintf(name,"%s",argv[1]);
	double         * input=load_wave_file(name,&ns);

	TimeFreq(input,ns,sampleRate,f,fband,&spec);
	double amin=-0.70, amax=0.0, da=0.01;
	double bmin=-0.70, bmax=0.0, db=0.001;

	int frame;
	int nb=10;
	int imax[nb];
	double vmax[nb];
	double vtmax;
	for(frame=0;frame<spec.frameCount;frame++)
	{
		int size;
		double * v = GetSpec(&spec,frame,fmax,fmin,&size);
		int fftLen=spec.winLength;
		int	imin=ceil(fmin/sampleRate*fftLen);
		int	ixmax=ceil(fmax/sampleRate*fftLen);
		int nfreq=ixmax-imin;

		find_max_divisors(v,nfreq,1,vmax,imax);
		printf("%d %lg %lg ",frame,fmin+(fmax-fmin)*imax[0]/nfreq,vmax[0]);
		vtmax=vmax[0];
		find_max_divisors(v,nfreq,nb,vmax,imax);
		for(i=0;i<nb;i++)
			printf("%lg %lg ",fmin+(fmax-fmin)*imax[i]/nfreq,vmax[i]/vtmax);
		printf("\n");

		free(v);
	}
	free(spec.data);
	free(spec.input);
	
}


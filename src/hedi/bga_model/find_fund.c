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
	SPEC_DATA spec;
	double f=1000; 
	double sampleRate=44100.;
	double fband=25;
	double fmin=200;
	double fmax=8000;
	
	spec.f=f;
	spec.fband=fband;
	spec.sampleRate=sampleRate;

	double amin=-0.70, amax=0.0, da=0.01;
	double bmin=-0.70, bmax=0.0, db=0.001;

	double alpha,beta;
	FILE * gg=fopen("fun.dat","w");
	for(alpha=amin;alpha<amax;alpha+=da)
	{
		for(beta=alpha;beta<alpha+0.3;beta+=db)
		{
			int nb=10;
			int imax[nb];
			double vmax[nb];
			double vtmax;
			fprintf(gg,"%lg %lg ",alpha,beta);
			get_vocal(alpha,beta,&spec);
			double *v=get_average_spectrum(&spec,fmax,fmin);
			
			int fftLen=spec.winLength;
			int	imin=ceil(fmin/sampleRate*fftLen);
			int	ixmax=ceil(fmax/sampleRate*fftLen);
			int nfreq=ixmax-imin;


			find_max_divisors(v,nfreq,1,vmax,imax);
			fprintf(gg,"%lg %lg ",fmin+(fmax-fmin)*imax[0]/nfreq,vmax[0]);
			vtmax=vmax[0];
			find_max_divisors(v,nfreq,nb,vmax,imax);
			for(i=0;i<nb;i++)
				fprintf(gg,"%lg %lg ",fmin+(fmax-fmin)*imax[i]/nfreq,vmax[i]/vtmax);
			fprintf(gg,"\n");


			fprintf(stdout,"%lg %lg ",alpha,beta);
			fprintf(stdout,"%lg %lg\n",fmin+(fmax-fmin)*imax[0]/nfreq,vmax[0]);


			free(v);
			free(spec.data);
			free(spec.input);
		}
	}
	fclose(gg);
}


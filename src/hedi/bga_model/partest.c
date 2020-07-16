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
	
	


	int nfilter=3;
	int m=2+2*nfilter;
	double fmax=8000;
	double fmin=200;
	double alpha=atof(argv[2]);
	double beta=atof(argv[3]);
	double f1=atof(argv[4]);
	double f2=atof(argv[5]);
	double f3=atof(argv[6]);
	double a1=atof(argv[7]);
	double a2=atof(argv[8]);
	double a3=atof(argv[9]);
	spec.f=f;
	spec.fband=fband;

	double slope=-0.6025;
	double intercept=-0.5979;
	beta=slope*alpha+intercept;


	TimeFreq(input,ns,sampleRate,f,fband,&spec);
	
	int frame=spec.frameCount/2;
	interp_alloc(nfilter+2);
	
	double guess[]={alpha,beta,-1000,f1,f2,f3,10000.,0.0,a1,a2,a3,0.0};
	spec_distance(guess[0], guess[1], guess+2, guess+7, 5,frame, fmin, fmax, &spec,1);

}


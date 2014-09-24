#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_complex.h>
#include "spectrum.h"
     
int _iGaussianSpectrum(double *  data,double * input, int inputsize,int increment, int winLength,int frameCount)
{
	int i,k;
	int fftLen=winLength;
	double ppinput[inputsize+winLength];
	for(i=0;i<inputsize+winLength;i++) ppinput[i]=0.;
	/* complex conversion for fft*/
	for(i=0;i<inputsize;i++) ppinput[i+winLength/2]=input[i];	

 	gsl_fft_complex_wavetable * wavetable;
    gsl_fft_complex_workspace * workspace;

	int nstd = 6;     
	double ws[fftLen];
	double wvar=(winLength/nstd)*(winLength/nstd);	

	for(i=0;i<winLength;i++)
		ws[i]=exp(-0.5*((i+1)-(winLength+1.)/2)*((i+1)-(winLength+1.)/2)/wvar);
	
	wavetable = gsl_fft_complex_wavetable_alloc (fftLen);
    workspace = gsl_fft_complex_workspace_alloc (fftLen);
    double buffer[2*fftLen];

	for(i=0;i<frameCount;i++)
	{
		int start=i*increment;
		int last=start+winLength;
		for(k=0;k<fftLen;k++)
		{
			buffer[2*k]=ws[k]*ppinput[start+k];
			buffer[2*k+1]=0.;
		}
		gsl_fft_complex_forward (buffer, 1, fftLen, wavetable, workspace);		
 		int offset=i*2*fftLen;
		memcpy(data+offset,buffer,sizeof(double)*2*fftLen);
	}   
    gsl_fft_complex_wavetable_free (wavetable);
    gsl_fft_complex_workspace_free (workspace);
}


int GaussianSpectrum(double * input, int inputsize,int increment, int winLength,SPEC_DATA* spec)
{
	/* compute the gaussian spectrum according to the tlab matlab file*/

	if (winLength%2==1) winLength +=1;
	int i,k; 
	int newLength=inputsize+winLength;
	int fftLen=winLength;
	int frameCount=floor((newLength-winLength)/increment)+1;
	double * data=(double*)calloc(2*(fftLen)*frameCount,sizeof(double));
	_iGaussianSpectrum(data,input,inputsize,increment,winLength,frameCount);
    spec->frameCount=frameCount;
    spec->winLength=winLength;
    spec->size=inputsize+winLength;
    spec->increment=increment;
    spec->data=data;
    return 0;
}


double*  GaussianSpectrumS(double * input, int inputsize,int increment, int winLength,int  * fc, int *fft)
{
	/* compute the gaussian spectrum according to the tlab matlab file*/

	if (winLength%2==1) winLength +=1;

	int newLength=inputsize+winLength;
	*fft=winLength;
	int fftLen=winLength;
	int frameCount=floor((newLength-winLength)/increment)+1;
	*fc=frameCount;
	double * data=(double*)calloc(2*(fftLen)*frameCount,sizeof(double));
	_iGaussianSpectrum(data,input,inputsize,increment,winLength,frameCount);	

    return data;
}

double *  _iInvertAndAdd(double * spec,int frameCount,int winLength,int increment,  double * data)
{
	int fftLen=winLength;
	int i,k;
	/* gaussian window preparation */	
	int nstd = 6;     
	double ws[fftLen];
	double wvar=(winLength/nstd)*(winLength/nstd);
	for(i=0;i<winLength;i++)
		ws[i]=exp(-0.5*((i+1)-(winLength+1)/2)*((i+1)-(winLength+1)/2)/wvar);
	/* fft stuff */
	gsl_fft_complex_wavetable * wavetable;
    gsl_fft_complex_workspace * workspace;	
	wavetable = gsl_fft_complex_wavetable_alloc (fftLen);
    workspace = gsl_fft_complex_workspace_alloc (fftLen);

    double buffer[2*fftLen];
    double win[fftLen];
    double y[frameCount*increment+winLength];
    double c[frameCount*increment+winLength];


	for(i=0;i<frameCount;i++)
	{
		int offset=i*2*fftLen;
		int first = i*increment;

		memcpy(buffer,spec+offset,sizeof(double)*2*fftLen);

		gsl_fft_complex_inverse(buffer, 1, fftLen, wavetable, workspace);		
		

		for(k=0;k<fftLen;k++)
			win[k]=buffer[2*k]*ws[k];

		for(k=0;k<fftLen;k++)
			{
				y[first+k] += win[k];
				c[first+k] += ws[k]*ws[k];
			}

	}

	for(i=0;i<frameCount*increment+winLength;i++)
		{			
			y[i] /=c[i];

		}
	memcpy(data,y+winLength/2,frameCount*increment*sizeof(double));

	return data;
}


double *  InvertAndAdd(SPEC_DATA * spec_data)
{
	int frameCount=spec_data->frameCount;
	int winLength=spec_data->winLength;
	int increment=spec_data->increment;
	double * spec=spec_data->data;
	double * data=(double*)calloc(frameCount*increment,sizeof(double));
	return _iInvertAndAdd(spec,frameCount,winLength,increment,data);
}

void InvertAndAddS(double * spec,int frameCount,int winLength,int increment,double * data)
{
	 _iInvertAndAdd(spec,frameCount,winLength,increment,data);
	 
}




int TimeFreq(double * input, int inputsize,double sampleRate, double f, double fband, SPEC_DATA* spec)
{
	int nstd=6;
//	double fband=125;
	double twindow = nstd/(fband*2.0*M_PI);
    int winLength = round(twindow*sampleRate);
    winLength = round(winLength/2)*2; 
	int    increment = round(sampleRate/f); 
    spec->sampleRate=sampleRate;
    spec->input=input;
//    printf("tw %lg wl %d inc %d\n",twindow,winLength,increment);
 //   printf("%d %d %d\n",__LINE__,inputsize,increment );
	GaussianSpectrum(input, inputsize, increment, winLength, spec); 

}

double * GetSpec(SPEC_DATA*spec_data, int frame, double fmax, double fmin,int * spsize)
{
	int frameCount=spec_data->frameCount;
	int winLength=spec_data->winLength;
	int fftLen=winLength;
	int increment=spec_data->increment;
	double sampleRate=spec_data->sampleRate;
	double * data=spec_data->data;
	int i,k;

	if ((frame<0)||(frame>=frameCount)) return 0x0;
	int offset=frame*fftLen;
	int imax,imin;
	imin=ceil(fmin/sampleRate*fftLen);
	imax=ceil(fmax/sampleRate*fftLen);
	double * sp=(double*)calloc((imax-imin),sizeof(double));
	*spsize=imax-imin;
	
	for(i=0;i<imax-imin;i++)
	{
		double dx=data[2*(i+imin+offset)],dy=data[2*(i+imin+offset)+1];
		double v=sqrt(dx*dx+dy*dy);
		sp[i]=v;
	}	
	return sp;
}
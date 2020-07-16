#ifndef SPEC_H
#define SPEC_H

typedef struct 
{
	int frameCount;
    int winLength;
    int size;
    int increment;
    double sampleRate;
    double * data;
    double * input;
    double fband;
	double f;  
} SPEC_DATA;

int GaussianSpectrum(double * input, int inputsize,int increment, int winLength,SPEC_DATA* spec);
double * GaussianSpectrumS(double * input, int inputsize,int increment, int winLength,int  * fc, int *fft);
int TimeFreq(double * input, int inputsize,double sampleRate, double f, double fband,SPEC_DATA* spec);
double * GetSpec(SPEC_DATA*spec, int frame, double fmax, double fmin, int * size);
double *  _iInvertAndAdd(double * spec,int frameCount,int winLength,int increment,  double * data);
double *  InvertAndAdd(SPEC_DATA * spec_data);
void InvertAndAddS(double * spec,int frameCount,int winLength,int increment,double * data);

#endif
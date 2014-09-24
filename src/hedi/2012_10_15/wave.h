#ifndef WAVE_H
#define WAVE_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sndfile.h>



typedef struct 
{
	double * data; // whole data (including surroundings)

} wave_slice;

typedef struct 
{
	wave_slice ** slices;
	double *  wave_data; // formatted data (zero padded before and after)
	int ns;  // actual samples # then real data is nfft/2:nfft/2+ns 
	int nslices; 
	int nfft;
} wave_array;

void save_wave_file(char * fname,double * data,int size);
double * load_wave_file(char * fname, int * ns);



wave_slice * create_slice(double * wave_data,int start,int nfft);
wave_array * create_array(char * fname,int nfft);


#endif
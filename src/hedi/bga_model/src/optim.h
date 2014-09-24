#ifndef OPTIM_H
#define OPTIM_H

#include "spectrum.h"
typedef struct 
{
	double y;
	double x;
	/* data */
} param_sort;



void interp_alloc(int nbfilter);
double spec_distance(double alpha, double beta, double *x, double *y, int nbfilter, int frame, double fmin, double fmax, SPEC_DATA * spec, int verbose);
double * get_average_spectrum(SPEC_DATA*fit_spec,double fmax,double fmin);
void get_vocal(double alpha,double beta,SPEC_DATA * fit_spec);
#endif
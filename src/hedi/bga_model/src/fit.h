#ifndef FIT_H
#define FIT_H

typedef struct 
{
	SPEC_DATA * target_spec;
	int nfilter;
	int current_frame;
	double * best_current_frame;
	double fmin;
	double fmax;
	double slope;
	double intercept;
}FIT;


double fitOneFrame(FIT * fit,int frame, double * guess);
double fitOneFramev2(FIT * fit,int frame, double slope, double intercept,double * guess);

#endif 
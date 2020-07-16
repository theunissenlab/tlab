#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <sndfile.h>
#include "utils.h"


void
save_wave_file(char *fname, double *data, int size)
{
	SNDFILE        *file;
	SF_INFO         sfinfo;

	sfinfo.samplerate = 44100;
	sfinfo.frames = size;
	sfinfo.channels = 1;
	sfinfo.format = (SF_FORMAT_WAV | SF_FORMAT_PCM_16);
//	double          Mx = getABSmax(data, size);
	//scaleData(1 / Mx, data, size);
	file = sf_open(fname, SFM_WRITE, &sfinfo);
	sf_write_double(file, data, sfinfo.channels * size);
	sf_close(file);
}

double         *
load_wave_file(char *fname, int *ns)
{
	SNDFILE        *file;
	SF_INFO         sfinfo;
	file = sf_open(fname, SFM_READ, &sfinfo);
	(*ns) = sfinfo.frames;

	double         *data = (double *) calloc((*ns), sizeof(double));
	int k = sf_readf_double(file, data, (*ns));
	int i;


	sf_close(file);
	return data;
}

double
getABSmax(double *data, int size)
{
	int             i;
	double          Max = -1;
	for (i = 0; i < size; i++)
		if (fabs(data[i]) > Max)
			Max = fabs(data[i]);
	return Max;
}
void
scaleData(double scale, double *data, int size)
{
	int             i;
	double          Max = -1;
	for (i = 0; i < size; i++)
		data[i] *= scale;
}

void normalize(double *data,int size)
{
	double m=getABSmax(data, size);
	scaleData(1/m,data,size);
}

double
CHI2(double *a, double *b, int size)
{
	double          d = 0;
	int             i;
	for (i = 0; i < size; i++) {
		//printf("%d %e %e\n", i, a[i], b[i]);
		d += (a[i] - b[i]) * (a[i] - b[i]);
	}
	return d;
}
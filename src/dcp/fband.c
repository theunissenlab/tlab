#include <stdio.h>
#include <malloc.h>
#include <string.h>
#include <math.h>
#include <nrutil.h>
#include <nr2.h>
#include <stdlib.h>


#include "dcp_ft.h"

#define GAUSSIAN /* Overlapping Gaussian window */
/* GAUSSIAN or HANNING or BOXCAR must be defined */

#define log2(x) log(x)/log(2.0)

/* ----------------------------- stim_fband() ------------------------------- */
stim_fband(pstim, nband, len, twindow, pow_level, f_low, f_high, f_width, f_step, amp_samp)
DCP_STIM *pstim;
int *nband;
int *len;
int twindow;             /* window of silence to add to stimulus on both sides in ms */
double ***pow_level;
float *f_low;
float *f_high;
float *f_width;
float *f_step;
float amp_samp;
{
int i, j, nb, c_n;
int nframes, nbands;
int istart, nwindow;
int i_low, i_high;
double n_step;
double *pr;
float *c_song, *c_temp;
double f2, fres, df;
double fresamp2, fresampmean;
double expval;
double x, y, tstep;
double fval, fmean, i_val;
double a_val, a_low, a_high;
char fname[128];
int wlen;
FILE *f_save, *f_debug;
/* #define DEBUG */

   /* Allocate memory for amplitude enveloppes */

	nframes = pstim->length;            /* Length of stim in points */
	*nband = nbands = ((*f_high)-(*f_low))/(*f_step);  /* Number of frequency bands */
	tstep = amp_samp/pstim->samprate;
	*len = (int)ceil(nframes*tstep);   /* Length of stim in sampling rate */
	wlen = (int)floor(twindow*amp_samp/1000.0);

	if ( tstep > 1.0 ) 
	{
		 printf("Error: Sampling rate is too low\n");
		 exit(1);
	}
	*pow_level = (double **)calloc(nbands+1, sizeof(double *) );
	for ( nb = -1; ++nb <= nbands; ) (*pow_level)[nb] = (double *)calloc(*len+2*wlen, sizeof(double));

	/* Allocate memory for complex array */
	nwindow = (int)floor(twindow*pstim->samprate/1000.0);
	c_n = (int)(log2((double)(nframes+2*nwindow))+0.5);
	c_n = (int)(pow(2.0,(double)c_n)+0.1);
	if ( c_n < (nframes+2*nwindow) ) c_n *=2;
	fres = ((float)pstim->samprate)/c_n;
   istart = (c_n-nframes)/2;
	printf("c_n = %d istart = %d nframes = %d nwindow=%d\n", c_n, istart, nframes, nwindow);

   c_song = (float *)calloc(c_n*2, sizeof(float));
	c_temp = (float *)calloc(c_n*2, sizeof(float));

	/* Stuff complex array */
	for ( i = -1; ++i < nframes; ) c_song[2*(i+istart)] = pstim->stimulus[i];

	/* Go to the fourier domain and zero out negative freqs */
	four1(c_song-1, (unsigned long)c_n, (int)1);

	for ( i = c_n; ++i < 2*c_n; ) c_song[i] = 0.0;

	/* Calculate analytical signal for each band */
	fres = ((float)pstim->samprate)/c_n;
	*f_step = (*f_high-*f_low)/nbands;
	f2=(*f_width)*(*f_width);
	printf ("New frequency step: low = %g high = %g step=%g\n", *f_low, *f_high, *f_step);

	for ( nb = -1; ++nb < nbands; )
	{
		bcopy(c_song, c_temp, sizeof(float)*c_n*2);
		fmean = *f_low + (nb+0.5)*(*f_step);
		
		for ( i = -1; ++i <= c_n/2; )
		{
			fval = i*fres;
			df = fval - fmean;
#ifdef GAUSSIAN
			expval = -0.5*df*df/f2;
			if ( expval > -10 )
			{
				expval = exp(expval);
				c_temp[2*i] *= expval;
				c_temp[2*i+1] *= expval;
			}
			else
			{
				c_temp[2*i] = 0.0;
				c_temp[2*i+1] = 0.0;
			}
#endif
#ifdef BOXCAR
        if ( fabs(df) >= *f_width/2 )
			{
				c_temp[2*i] = 0.0;
				c_temp[2*i+1] = 0.0;
			}
#endif
#ifdef HANNING
        if ( fabs(df) >= *f_width/2 )
			{
				c_temp[2*i] = 0.0;
				c_temp[2*i+1] = 0.0;
			}
			else
			{
				expval = 0.5*(1+cos((2.0*M_PI)*(df/(*f_width))));
				c_temp[2*i] *= expval;
				c_temp[2*i+1] *= expval;
			}
#endif
		}
		four1(c_temp-1, (unsigned long)c_n, (int)-1);

		/* Calculate amplitude enveloppe */
		for ( i = -1; ++i < c_n; )
		{
			x = c_temp[2*i];
			y = c_temp[2*i+1];
			c_temp[2*i] = x*x+y*y;
			c_temp[2*i+1] = 0.0;
		}

#ifdef DEBUG
		f_debug = fopen("before.txt","w");
		for ( i = -1; ++i < c_n; ) fprintf(f_debug,"%g\n",c_temp[2*i]);
		fclose(f_debug);
#endif

		/* Low pass filter Amplitude before resampling if needed */

		if ( amp_samp  < 10.0*(*f_width) )
		{
 /* I need to implement a different type of filter - some type of MA filter. Also need to figure out whether I need to filter in the frequency domain.. */
			four1(c_temp-1, (unsigned long)c_n, (int)1);
			fresampmean = amp_samp/2;
			fresamp2=25*25;
			for ( i = -1; ++i <= c_n/2; )
			{
				fval = i*fres;
				if ( fval < fresampmean ) continue;
				df = fval - fresampmean;
				expval = -0.5*df*df/fresamp2;
				if ( expval > -100 )
				{
					expval = exp(expval);
					c_temp[2*i] *= expval;
					c_temp[2*i+1] *= expval;
					if ( i != 0 && i != c_n/2 )
					{
						c_temp[2*c_n-2*i] *= expval;
						c_temp[2*c_n-2*i+1] *= expval;
					}
				}
				else
				{
					c_temp[2*i] = 0.0;
					c_temp[2*i+1] = 0.0;
					if ( i != 0 && i != c_n/2 )
					{
						c_temp[2*c_n-2*i] = 0.0;
						c_temp[2*c_n-2*i+1] = 0.0;
					}
				}
			}
			four1(c_temp-1, (unsigned long)c_n, (int)-1);
			for ( i = -1; ++i < c_n; ) c_temp[2*i] /= (double)c_n;
#ifdef DEBUG
		f_debug = fopen("after.txt","w");
		for ( i = -1; ++i < c_n; ) fprintf(f_debug,"%g\n",c_temp[2*i]);
		fclose(f_debug);
		printf("Debug files printed\n");
		printf("Enter character to continue\n");
		getchar();
#endif
		}

		/* Stuff Amplitude */
		for ( j = -1; ++j < *len + 2*wlen; )
		{
			i_val = (j-wlen)/tstep;
			i_low = (int)floor(i_val) + istart;
			i_high = (int)ceil(i_val) + istart;
			i_val += istart;

			x = c_temp[2*i_low];
			if ( x < 0.0 )
			{
				printf("i_low = %d i_val = %f x = %f\n",i_low,i_val,x);
				x = 0.0;
			}
			a_low = x;
			a_val = (1.0+i_low-i_val)*a_low;

			if (i_low != i_high)
			{
				x = c_temp[2*i_high];
				if ( x < 0.0 ) x = 0.0;
				a_high = x;
				a_val += (1.0+i_val-i_high)*a_high;
			}
			(*pow_level)[nb][j] = sqrt(a_val);
			(*pow_level)[nbands][j] += sqrt(a_val); 
		}
	}
	free(c_temp);
	free(c_song);
}

#include <stdio.h>
#include <malloc.h>
#include <string.h>
#include <math.h>
#include <nrutil.h>
#include <nr2.h>
#include <stdlib.h>


#include "dcp_ft.h"

#define log2(x) log(x)/log(2.0)

/* ----------------------------- stim_fband() ------------------------------- */
stim_fband(pstim, nband, len, twindow, pow_level, f_low, f_high, f_width)
DCP_STIM *pstim;
int *nband;
int *len;
int twindow;             /* window of silence to add to stimulus on both sides in ms */
double ***pow_level;
float *f_low;
float *f_high;
float *f_width;
{
#define GAUSSIAN

int i, j, nb, c_n;
int nframes, nwindow, nbands;
int istart, n_step;
double *pr;
float *c_song, *c_temp;
double fres;
int nb_pos, n_group, ind_f_low, ind_f_start;
double ind_f_mean, ind_f_std, expval;
double x, y, frac;
double fstep, fval;
char fname[128];
FILE *f_save;



   /* Allocate memory for amplitude enveloppes */
	nframes = pstim->length;
	*nband = nbands = ((*f_high)-(*f_low))/(*f_width);
	n_step = pstim->samprate/1000.0;
	*len = nframes/n_step;
	nwindow = twindow*n_step;
	if ( nframes%n_step)  *len += 1;
	*pow_level = (double **)calloc(nbands+1, sizeof(double *) );
	for ( nb = -1; ++nb <= nbands; ) (*pow_level)[nb] = (double *)calloc(*len+2*twindow, sizeof(double));

   /* Check to see if computation is on disk */
	if ( strcmp(pstim->type,"song") == 0 || (strcmp(pstim->type,"seg") == 0 && strcmp(pstim->segcmd,"1-end") == 0) )
	{
		sprintf(fname,"%s/%s.fband",pstim->songpath,pstim->songfile);
		printf("%s:\n",fname); 
	}

	/* Allocate memory for complex array */
	c_n = (int)(log2((double)(nframes+2*nwindow))+0.5);
	c_n = (int)(pow(2.0,(double)c_n)+0.1);
	if ( c_n < (nframes+2*nwindow) ) c_n *=2;
	fres = ((float)pstim->samprate)/c_n;
	nb_pos = ((*f_high)-(*f_low))/fres;
	n_group = (int)(nb_pos/nbands + 0.5);
	while ( n_group == 0 )
	{
		c_n *=2;
		fres = ((float)pstim->samprate)/c_n;
		nb_pos = ((*f_high)-(*f_low))/fres;
		n_group = (int)(nb_pos/nbands + 0.5);
	}

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
	nb_pos = ((*f_high) - (*f_low))/fres;
	n_group = (int)(nb_pos/nbands + 0.5);
	ind_f_low = (*f_low)/fres;
	printf ("Old frequency bounds : low = %g high = %g step %g\n", (*f_low), (*f_high), *f_width);
	*f_low = ind_f_low*fres;
	*f_high = (ind_f_low+nbands*n_group)*fres;
	fstep = (*f_high-*f_low)/nbands;
	*f_width = fstep;
	printf ("New frequency bounds : low = %g high = %g step=%g\n", *f_low, *f_high, *f_width);

	for ( nb = -1; ++nb < nbands; )
	{
		bcopy(c_song, c_temp, sizeof(float)*c_n*2);
		
		ind_f_start = ind_f_low + nb*n_group;
		ind_f_mean = ind_f_start + (n_group-1)/2.0;
		ind_f_std = n_group;  /* The std with the best approx to square */
		for ( i = -1; ++i <= c_n/2; )
		{
		   if ( i < ind_f_start )
			{
#ifdef GAUSSIAN
				expval = -0.5*(i-ind_f_mean)*(i-ind_f_mean)/(ind_f_std*ind_f_std);
				if ( expval > -10 )
				{
					c_temp[2*i] *= exp(expval);
					c_temp[2*i+1] *= exp(expval);
				}
				else
				{
					c_temp[2*i] = 0.0;
					c_temp[2*i+1] = 0.0;
				}
#else
				c_temp[2*i] = 0.0;
				c_temp[2*i+1] = 0.0;
#endif
			}
		   else if  ( i < ind_f_start+n_group )  
			{
#ifdef HANNING
				/* Hanning window in frequency */
				expval =  0.5*(1-cos(2.0*M_PI*(i-ind_f_start)/(double)n_group));
				c_temp[2*i] *= expval;
				c_temp[2*i+1] *= expval;
#endif
#ifdef WELCH
            /* Welch window in frequency */
				frac = i-ind_f_start-n_group/2;
				frac = frac*frac;
				frac /= n_group*n_group/4;
				c_temp[2*i] *= 1.0-frac;
				c_temp[2*i+1] *= 1.0-frac;
#endif
#ifdef GAUSSIAN
				expval = -0.5*(i-ind_f_mean)*(i-ind_f_mean)/(ind_f_std*ind_f_std);
				if ( expval > -10 )
				{
					c_temp[2*i] *= exp(expval);
					c_temp[2*i+1] *= exp(expval);
				}
				else
				{
					c_temp[2*i] = 0.0;
					c_temp[2*i+1] = 0.0;
				}
#endif
			}
			else 
			{
#ifdef GAUSSIAN
				expval = -0.5*(i-ind_f_mean)*(i-ind_f_mean)/(ind_f_std*ind_f_std);
				if ( expval > -10 )
				{
					c_temp[2*i] *= exp(expval);
					c_temp[2*i+1] *= exp(expval);
				}
				else
				{
					c_temp[2*i] = 0.0;
					c_temp[2*i+1] = 0.0;
            }
#else
				c_temp[2*i] = 0.0;
				c_temp[2*i+1] = 0.0;
#endif
			}
		}
		four1(c_temp-1, (unsigned long)c_n, (int)-1);

		/* Stuff Amplitude */
		fval = fstep*(nb+0.5) + *f_low;
		for ( i = -1; ++i < nframes + 2*nwindow; )
		{
			 
			if (i%n_step == 0 )
			{
				j = i/n_step;
				x = c_temp[2*(i+istart-nwindow)];
				y = c_temp[2*(i+istart-nwindow)+1];
				(*pow_level)[nb][j] = sqrt(x*x + y*y);
				(*pow_level)[nbands][j] += x*x + y*y; 
			}
		}
	}
	free(c_temp);
	free(c_song);
}

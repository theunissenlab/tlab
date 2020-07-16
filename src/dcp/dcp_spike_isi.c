#include <stdio.h>
#include <malloc.h>
#include <string.h>
#include <math.h>
#include <nrutil.h>
#include <nr2.h>
#include <time.h>
#include <stdlib.h>


#include "dcp_ft.h"
DCP_DATA *read_raw();

main(argc, argv)
int argc;
char *argv[];
{
int it, ib, st;
int i, j, k;
int bin0;
#define TWINDOW 700
DCP_DATA *pin;
DCP_TRIAL *ptin;
DCP_DATA *read_data();
float f_low, f_high;
double time_spike, time_spike2;
double **stim_env;        /* Set of stimulus envelopes in ms units */
double **stim_avg;        /* Average Set of stimulus envelopes in ms units */
double **avg_env;
double **rand_env;
double *cross_spike;
double **filter;
double *average_rand;
int nband, nlen, trial_spikes;
int nspikes=0;
int *spike_ns, *stim_ns, *rand_ns;
int tzero2, tzero;
int init_flg = 0;
int rawfile = 0;
int count_stim;
int moredat;
FILE *f_spike, *f_spike2;

	if ( argc < 2 )
	{
		printf("Error: dcp_spike_avg takes at least 1 argument\n");
		printf("dcp_spike_avg {-raw} <file1> <file2> ... \n");
		printf("    -raw reads raw file. otherwise reads dcp file.\n");
		exit(1);
	}

   argc --;
	argv ++;

	while ( argv[0][0] == '-' )
	{
		if ( strncmp(argv[0],"-raw") == 0 ) rawfile = 1;
		else 
		{
			printf("Error: Unknown flag %s\n", argv[0]);
			exit(1);
		}
		argc --;
		argv ++;
	}

	init_flg = 1;
	tzero2 = tzero = 0;
	cross_spike = (double *)calloc(2*TWINDOW, sizeof(double));
	stim_ns = (int *)calloc(2*TWINDOW, sizeof(int));
	rand_ns = (int *)calloc(2*TWINDOW, sizeof(int));
	spike_ns = (int *)calloc(2*TWINDOW, sizeof(int));
	srand48((long)time(NULL));

   /* Calculate spike triggered averages and auto-correlation of spikes */
	f_spike = fopen("stim.avg","w");
	f_spike2 = fopen("avg_stim.avg","w");
   while (argc)
	{
		moredat = 1;
		while ( moredat )
		{
			if ( rawfile == 0 )
			{
				pin = read_data(argv[0]);
				if ( pin == NULL )
				{
					printf("Error reading dcp file %s\n", argv[0]);
					exit(1);
				}
				if ( pin->ssflag )
				{
					printf("Not yet dealing with spike sorted data\n");
					exit(1);
				}
				strcpy(stim_spec, pin->stim_spec);
				stim = pin->pstim = setup_stimulus(stim_spec);
				if ( stim == NULL )
				{
					printf("Error setting up stimulus for %s:%s\n", argv[1], pin->stim_spec);
					exit(1);
				}
				moredat = 0;
			}
			else
			{
				pin = read_raw(argv[0], &moredat);
				if (moredat==0) 
				{
               /* The following code is repeated down below */
					for (it = -1; ++it < nlen; )
					{
						fprintf(f_spike2, "%d ", tzero2+it);
						for ( ib = -1; ++ib <= nband; )
							fprintf(f_spike2, "%g ", stim_avg[ib][it]/count_stim);
						fprintf(f_spike2, "\n");
					}
					tzero2 += nlen;
					for ( ib = -1; ++ib <= nband; ) free(stim_avg[ib]);
					free(stim_avg);
					break;
				}
            stim = pin->pstim;
			}

			/* Pre-process the stimulus */
			stim_fband(stim,&nband,&nlen,&stim_env,&f_low,&f_high);

			if ( init_flg )
			{
				avg_env = (double **)calloc(nband+1, sizeof(double *));
				rand_env = (double **)calloc(nband+1, sizeof(double *));
				for ( ib = -1; ++ib <= nband; ) 
				{
					avg_env[ib] = (double *)calloc(2*TWINDOW, sizeof(double));
					rand_env[ib] = (double *)calloc(2*TWINDOW, sizeof(double));
				}
			}

			if ( moredat != 2 )
			{
				stim_avg = (double **)calloc(nband+1, sizeof(double *));
				for ( ib = -1; ++ib <= nband; ) 
					stim_avg[ib] = (double *)calloc(nlen, sizeof(double));
				count_stim = 0;
			}
         for ( ib = -1; ++ib <= nband; )
         {
				for ( it = -1; ++it < nlen; )
					if ( stim_env[ib][it] != 0.0 )
						stim_avg[ib][it] += log(stim_env[ib][it]);
         }
			count_stim ++;

			for ( i = -1; ++i < pin->n; )
			{
				ptin = pin->ptrial + i;
				trial_spikes = 0;
				for ( j = -1; ++j < ptin->nspikes; )
				{
					time_spike = ptin->time[j]*1.0e-3;
					if ( time_spike < pin->pre - TWINDOW ) continue;
					if ( time_spike >= pin->pre + stim->dur + TWINDOW ) break;
					bin0 = nint(time_spike) - pin->pre - TWINDOW;
					nspikes ++;
					trial_spikes ++;
					for ( it = -1; ++it < 2*TWINDOW; )
					{
						st = bin0 + it;
						if ( st < -TWINDOW || (st + pin->pre) < 0 ) continue;
						if ( st >= stim->dur+TWINDOW || st >= (pin->cd-pin->pre) ) break;
						spike_ns[it] ++;
               }
					for ( k = -1; ++k < ptin->nspikes; )
					{
						time_spike2 = ptin->time[k]*1.0e-3;
						if ( time_spike2 < pin->pre - TWINDOW ) continue;
						if ( time_spike2 >= pin->pre + stim->dur + TWINDOW ) break;
						st = nint(time_spike2 - time_spike) + TWINDOW;
						if ( st < 0 ) continue;
						if ( st >= 2*TWINDOW ) break;
						cross_spike[st] += 1.0;
					}
					for ( it = -1 ; ++it < 2*TWINDOW; )
					{
						st = bin0 + it;
						if ( st < 0 ) continue;
						if ( st >= nlen ) break;
						stim_ns[it] ++;
						for ( ib = -1; ++ib <= nband; )
							if ( stim_env[ib][st] != 0.0 )
								avg_env[ib][it] += log(stim_env[ib][st]);
					}
				}
				for ( j = -1; ++j < trial_spikes; )
				{
					time_spike = pin->pre - TWINDOW + drand48()*(stim->dur + 2*TWINDOW);
					bin0 = nint(time_spike) - pin->pre - TWINDOW;
					for ( it = -1 ; ++it < 2*TWINDOW; )
					{
						st = bin0 + it;
						if ( st < 0 ) continue;
						if ( st >= nlen ) break;
						rand_ns[it] ++;
						for ( ib = -1; ++ib <= nband; )
								if ( stim_env[ib][st] != 0.0 )
									rand_env[ib][it] += log(stim_env[ib][st]);
					}
				}
			}
			free_dcp_memory(&pin);
			if ( moredat == 0 || init_flg) 
			{
				for (it = -1; ++it < nlen; )
				{
					fprintf(f_spike, "%d ", tzero+it);
					for ( ib = -1; ++ib <= nband; )
						fprintf(f_spike, "%g ", stim_env[ib][it]);
					fprintf(f_spike, "\n");
				}
				tzero += nlen;
				init_flg = 0;
			}
			free_env(nlen,nband,stim_env);
			if (moredat==0) 
			{
				for (it = -1; ++it < nlen; )
				{
					fprintf(f_spike2, "%d ", tzero+it);
					for ( ib = -1; ++ib <= nband; )
						fprintf(f_spike2, "%g ", stim_avg[ib][it]/count_stim);
					fprintf(f_spike2, "\n");
				}
				tzero2 += nlen;
				for ( ib = -1; ++ib <= nband; ) free(stim_avg[ib]);
				free(stim_avg);
				break;
			}
		}
		argc --;
		argv ++;
	}
	fclose(f_spike);
	fclose(f_spike2);
				 
	f_spike = fopen("spike.avg","w");
	for (it = -1; ++it < 2*TWINDOW; )
	{
		fprintf(f_spike, "%d ", it-TWINDOW);
		for ( ib = -1; ++ib <= nband; )
		{
			avg_env[ib][it] /= (double)stim_ns[it];
			fprintf(f_spike, "%g ", avg_env[ib][it]);
		}
		fprintf(f_spike, "\n");
	}
	fclose(f_spike);

	f_spike = fopen("random.avg","w");
	for (it = -1; ++it < 2*TWINDOW; )
	{
		fprintf(f_spike, "%d ", it-TWINDOW);
		for ( ib = -1; ++ib <= nband; )
		{
			rand_env[ib][it] /= (double)rand_ns[it];
			fprintf(f_spike, "%g ", rand_env[ib][it]);
		}
		fprintf(f_spike, "\n");
	}
	fclose(f_spike);

	f_spike = fopen("cross.avg", "w");
	for ( it = -1; ++ it < 2*TWINDOW; )
	{
		cross_spike[it] /= (double)spike_ns[it];
		fprintf(f_spike, "%g\n", cross_spike[it]);
	}
	fclose(f_spike);

	/* Calculate the optimal rev reconstruction filter */
	filter = (double **)calloc(nband+1, sizeof(double *));
	average_rand = (double *)calloc(nband+1, sizeof(double) );
	for ( ib = -1; ++ib <= nband; )
	{
		filter[ib] = (double *)calloc(2*TWINDOW, sizeof(double));
		average_rand[ib] = 0.0;
		for (it = -1; ++it < 2*TWINDOW; ) average_rand[ib] += rand_env[ib][it];
		average_rand[ib] /= 2.0*TWINDOW;
		for (it = -1; ++it < 2*TWINDOW; ) filter[ib][it] = avg_env[ib][it] - average_rand[ib];
	}
	rev_filter(filter, cross_spike, 2*TWINDOW, nband+1);

	f_spike = fopen("filter.avg","w");
	for (it = -1; ++it < 2*TWINDOW; )
	{
		fprintf(f_spike, "%d ", it-TWINDOW);
		for ( ib = -1; ++ib <= nband; )
			fprintf(f_spike, "%g ", filter[ib][it]);
		fprintf(f_spike, "\n");
	}
	fclose(f_spike);

   /* free other stuff */
	/* free_filter() */
}
/* ----------------------------- rev_filter() ------------------------------- */
rev_filter( filter, cross_spike, nt, nb)
double **filter;      /* Returns reverse recon filter */
double *cross_spike;  /* Auto-cross correlation of stim and spike */
int nt;               /* number of time points */
int nb;               /* number of bands */
{
int it, ib, istart;
int c_n;
float *c_stim, *c_spike;
float div;
/* #define DEBUG */
#ifdef DEBUG
FILE *f_debug;
#endif


	c_n = (int)(log2((double)nt)+0.5);
	c_n = (int)(pow(2.0,(double)c_n)+0.1);
	if ( c_n < nt ) c_n *= 2;

#ifdef DEBUG
	f_debug = fopen("debug.avg","w");
#endif

   istart = (c_n-nt)/2;

   c_spike = (float *)calloc(c_n*2, sizeof(float));
	c_stim = (float *)calloc(c_n*2, sizeof(float));

	/* Stuff complex array */
	for ( it = -1; ++it < nt; ) c_spike[2*(it+istart)] = cross_spike[it]*cos((it-nt/2)*M_PI/nt);
	/* for ( it = -1; ++it < nt; ) c_spike[2*(it+istart)] = cross_spike[it]; */
#ifdef DEBUG
	for ( it = -1; ++it < nt; ) fprintf(f_debug,"%g ", cross_spike[it]);
	fprintf(f_debug,"\n");
#endif
	four1(c_spike-1, (unsigned long)c_n, (int)1);

   for ( ib = -1; ++ib < nb; )
	{
		for ( it = -1; ++it < c_n*2; ) c_stim[it] = 0.0;
		for ( it = -1; ++it < nt; ) c_stim[2*(it+istart)] = filter[ib][it]*cos((it-nt/2)*M_PI/nt); 
		/* for ( it = -1; ++it < nt; ) c_stim[2*(it+istart)] = filter[ib][it]; */
#ifdef DEBUG
      if ( ib == nb-1 )
			for ( it = -1; ++it < nt; ) fprintf(f_debug,"%g ", filter[ib][it]);
		fprintf(f_debug,"\n");
#endif
		four1(c_stim-1, (unsigned long)c_n, (int)1);
      for ( it = -1; ++it < c_n; ) 
		{ 
			div = sqrt( c_spike[2*it]*c_spike[2*it] + c_spike[2*it+1]*c_spike[2*it+1]);
			if ( div != 0.0 )
			{
				c_stim[2*it] /= div;
				c_stim[2*it+1] /= div;
			}
			else
			{ 
				c_stim[2*it] = 0.0;
				c_stim[2*it+1] = 0.0;
			}
		}
		four1(c_stim-1, (unsigned long)c_n, (int)-1);
		for ( it = -1; ++it < nt; ) filter[ib][it] = c_stim[2*(it+istart)]/c_n;
#ifdef DEBUG
      if ( ib == nb-1 )
			for ( it = -1; ++it < nt; ) fprintf(f_debug,"%g ", filter[ib][it]);
		fprintf(f_debug,"\n");
#endif
	}

#ifdef DEBUG
	fclose(f_debug);
#endif

   free((void *)c_spike);
   free((void *)c_stim);

}
/* ----------------------------- stim_fband() ------------------------------- */
stim_fband(pstim, nband, len, pow_level, f_low, f_high)
DCP_STIM *pstim;
int *nband;
int *len;
double ***pow_level;
float *f_low;
float *f_high;
{
#define FLOW 250.0
#define FHIGH 8000.0
#define FWIDTH 250.0
#define SAMPRATE 32000.0
#define GAUSSIAN

int i, j, nb, c_n;
int nframes, nbands;
int istart, n_step;
double *pr;
float *c_song, *c_temp;
double fres;
int nb_pos, n_group, ind_f_low, ind_f_start;
double ind_f_mean, ind_f_std, expval;
double x, y, frac;
double fstep, fval;


   /* Allocate memory for amplitude enveloppes */
	nframes = pstim->length;
	*nband = nbands = (FHIGH-FLOW)/FWIDTH;
	n_step = SAMPRATE/1000.0;
	*len = nframes/n_step;
	if ( nframes%n_step)  *len += 1;
	*pow_level = (double **)calloc(nbands+1, sizeof(double *) );
	for ( nb = -1; ++nb <= nbands; ) (*pow_level)[nb] = (double *)calloc(*len, sizeof(double));

	/* Allocate memory for complex array */
	c_n = (int)(log2((double)nframes)+0.5);
	c_n = (int)(pow(2.0,(double)c_n)+0.1);
	if ( c_n < nframes ) c_n *=2;
	fres = SAMPRATE/c_n;
	nb_pos = (FHIGH-FLOW)/fres;
	n_group = (int)(nb_pos/nbands + 0.5);
	while ( n_group == 0 )
	{
		c_n *=2;
		fres = SAMPRATE/c_n;
		nb_pos = (FHIGH-FLOW)/fres;
		n_group = (int)(nb_pos/nbands + 0.5);
	}

   istart = (c_n-nframes)/2;
	printf("c_n = %d istart = %d nframes = %d\n", c_n, istart, nframes);

   c_song = (float *)calloc(c_n*2, sizeof(float));
	c_temp = (float *)calloc(c_n*2, sizeof(float));

	/* Stuff complex array */
	for ( i = -1; ++i < nframes; ) c_song[2*(i+istart)] = pstim->stimulus[i];

	/* Go to the fourier domain and zero out negative freqs */
	four1(c_song-1, (unsigned long)c_n, (int)1);

	for ( i = c_n; ++i < 2*c_n; ) c_song[i] = 0.0;

	/* Calculate analytical signal for each band */
	fres = SAMPRATE/c_n;
	nb_pos = (FHIGH - FLOW)/fres;
	n_group = (int)(nb_pos/nbands + 0.5);
	ind_f_low = FLOW/fres;
	printf ("Old frequency bounds : low = %g high = %g\n", FLOW, FHIGH);
	*f_low = ind_f_low*fres;
	*f_high = (ind_f_low+nbands*n_group)*fres;
	fstep = (f_high-f_low)/nbands;
	printf ("New frequency bounds : low = %g high = %g\n", *f_low, *f_high);

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
		for ( i = -1; ++i < nframes; )
		{
			 
			if (i%n_step == 0 )
			{
				j = i/n_step;
				x = c_temp[2*(i+istart)];
				y = c_temp[2*(i+istart)+1];
				(*pow_level)[nb][j] = sqrt(x*x + y*y);
				(*pow_level)[nbands][j] += x*x + y*y; 
			}
		}
	}
	free(c_temp);
	free(c_song);
}

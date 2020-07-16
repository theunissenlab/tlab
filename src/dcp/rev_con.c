#include <stdio.h>
#include <malloc.h>
#include <string.h>
#include <math.h>
#include <nrutil.h>
#include <nr2.h>
#include <time.h>
#include <stdlib.h>


#include "dcp_ft.h"
#define SAMPRATE 32000.0
DCP_DATA *read_raw();

main(argc, argv)
int argc;
char *argv[];
{
int it, ib, st;
int i, j, k;
int bin0;
char stim_spec[MAX_SPEC_LEN];
DCP_DATA *pin;
DCP_TRIAL *ptin;
DCP_STIM *stim;
DCP_DATA *read_data();
double time_spike;
double **filter, **random;
double **rand_env;
double **stim_recon;
double *average_rand;
double *spike_trials;
double *recon_trials;
double **stim_env;
float f_low, f_high;
int nband, nlen, trial_spikes;
int tzero = 0, n_step;
int init_flg = 0;
int rawfile = 0;
int count_trials;
int moredat;
int nt, TWINDOW;
FILE *f_recon;
FILE *f_stim, *f_spike, *f_recontrials;

   if ( argc < 2 )
   {
      printf("Error: rev_con takes at least 1 argument\n");
      printf("rev_con {-raw} <file1> <file2> ... \n");
      printf("    -raw reads raw file. otherwise reads dcp file.\n");
      exit(1);
   }

   argc --;
   argv ++;

   while ( argv[0][0] == '-' )
   {
      if ( strncmp(argv[0],"-raw",4) == 0 ) rawfile = 1;
      else
      {
         printf("Error: Unknown flag %s\n", argv[0]);
         exit(1);
      }
      argc --;
      argv ++;
   }

	read_filter("random.avg", &random, &nband, &nt);
	average_rand = (double *)calloc(nband+1, sizeof(double));
	for ( ib = -1; ++ib <=nband; )
	{
		for ( it = -1; ++it < nt; ) average_rand[ib] += random[ib][it];
		average_rand[ib] /= (double)nt;
	}
	
	read_filter("filter.avg", &filter, &nband, &nt);
	TWINDOW = nt/2;

   /*  Reverse reconstruction */
	f_recon = fopen("recon.avg", "w");
	f_spike = fopen("spike.trials","w");
	f_recontrials = fopen("recon.trials","w");
	f_stim = fopen("stim.trial","w");
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
				init_flg = 1;
			}
			else
			{
				if ( moredat == 1 ) init_flg = 1;
				pin = read_raw(argv[0], &moredat);
				if ( pin != NULL ) stim = pin->pstim;
			}

			if ( pin != NULL )
			{
			/* Pre-process the stimulus */
			stim_fband(stim,&nband,&nlen,&stim_env,&f_low,&f_high);
/* This was for the version with no stim processing
			n_step = SAMPRATE/1000.0;
			nlen = stim->length/n_step;
			if ( stim->length%n_step ) nlen += 1;
*/
			if ( init_flg )
			{
				stim_recon = (double **)calloc(nband+1, sizeof(double));
				for ( ib = -1; ++ib <= nband; )
					stim_recon[ib] = (double *)calloc(nlen, sizeof(double));
				recon_trials = (double *)calloc(nlen, sizeof(double));
				spike_trials = (double *)calloc(nlen, sizeof(double));
				init_flg = 0;
				count_trials = 0;
			}
			for ( it = -1; ++it < nlen; ) 
				if ( stim_env[nband][it] != 0.0 ) 
					fprintf(f_stim,"%g ", log(stim_env[nband][it])); 
               
				else
					fprintf(f_stim,"0.0 ");
			fprintf(f_stim,"\n");
			fflush(f_stim);
			for ( i = -1; ++i < pin->n; )
			{
				ptin = pin->ptrial + i;
				for ( it = -1; ++it < nlen; ) 
				{
                spike_trials[it] = 0.0;
					 recon_trials[it] = 0.0;
				}
				for ( j = -1; ++j < ptin->nspikes; )
				{
					time_spike = ptin->time[j]*1.0e-3;
					if ( time_spike < pin->pre - TWINDOW ) continue;
					if ( time_spike > pin->pre + stim->dur + TWINDOW ) break;
					bin0 = nint(time_spike) - pin->pre - TWINDOW;
					spike_trials[bin0 + TWINDOW] = 1.0;
					for ( it = -1 ; ++it < 2*TWINDOW; )
					{
						st = bin0 + it;
						if ( st >=0  && st < nlen-1 )
						{
							for ( ib = -1; ++ib <= nband; )
								stim_recon[ib][st] += filter[ib][it];
							recon_trials[st] += filter[nband][it];
						}
					}
				}
				for ( it = -1; ++it < nlen; ) 
				{
					fprintf(f_recontrials,"%g ", recon_trials[it]);
					fprintf(f_spike,"%g ", spike_trials[it]);
				}
				fprintf(f_recontrials,"\n");
				fprintf(f_spike,"\n");
				fflush(f_recontrials);
				fflush(f_spike);
				count_trials ++;
			}
			free_dcp_memory(&pin);
			}
			if ( moredat == 0) 
			{
				for (it = -1; ++it < nlen; )
				{
					fprintf(f_recon, "%d ", it + tzero);
					for (ib = -1; ++ib <= nband; )
					{
						stim_recon[ib][it] =  average_rand[ib] + stim_recon[ib][it]/count_trials;
						fprintf(f_recon, "%g ", stim_recon[ib][it]);
					}
					fprintf(f_recon,"\n");
					fflush(f_recon);
				}
				tzero += nlen;
				free(stim_recon);
			}
		}
		argc --;
		argv ++;
	}
	fclose(f_recon);

   /* free other stuff */
	/* free_filter() */
}
/* ----------------------------- read_filter() ----------------------------- */
int read_filter(fname, filter, nband, nt)
/* Reads the filte file create by dcp_spike_avg */
char *fname;
double ***filter;
int *nband;
int *nt;
{
int ib;
int it;
int l;
FILE *f_filter;
#define BUFFSIZE 10000
char buff[BUFFSIZE];
char *ptr;

   f_filter = fopen(fname,"r");
	if ( f_filter == NULL ) 
	{
		printf("Error: could not open file %s\n", fname);
		exit(1);
	}

	/* Read first line to get number of bands */
	if ( fgets(buff, BUFFSIZE, f_filter ) == NULL ) 
	{
		printf("Error reading %s\n",fname);
		exit(1);
	}
	if ( (l=strlen(buff)) == BUFFSIZE ) 
	{
		printf("Line too long in %s\n",fname);
		exit(1);
	}
	buff[l-1]=0;
	it = 1;
	if ( strtok(buff," ") == NULL )
	{
		printf("Error reading %s\n",fname);
		exit(1);
	}
	ib = 1;
	while ( strtok(NULL," ") ) ib ++;
	*nband = ib - 2;
	while ( fgets(buff, BUFFSIZE, f_filter) ) it += 1;
	*nt = it;

   printf("Number of bands %d. Number of time points %d\n", *nband, *nt);
	/* allocate space for filter */
	*filter = (double **)calloc(*nband+1, sizeof(double *));
   for ( ib = -1; ++ib <= *nband; )
      (*filter)[ib] = (double *)calloc((*nt), sizeof(double));

   rewind(f_filter);


	it = 0;
	while ( fgets(buff, BUFFSIZE, f_filter) ) 
	{
		ib = 0;
		l = strlen(buff);
		buff[l-1]=0;
      strtok(buff," "); /* first element it time */
		while ( (ptr=strtok(NULL," ")) != NULL )
		{
			(*filter)[ib][it] = atof(ptr);
			ib ++;
		}
		it ++;
	}

	fclose(f_filter);
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

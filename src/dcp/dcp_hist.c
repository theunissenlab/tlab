#include <stdio.h>
#include <malloc.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "dcp_ft.h"

main(argc, argv)
int argc;
char *argv[];
{
int i,j;
char buff[256];
DCP_DATA *pin;
DCP_DATA *read_data();
DCP_STIM *pstim;
DCP_TRIAL *ptin;
FILE *f_out;   /* Output file that has number of bins for pre, stim and post */
FILE *f_pre;
FILE *f_stim;
FILE *f_post;
FILE *f_stim_power;
FILE *f_fix;
FILE *f_dp;
float tbin;
int ibin;
int nbins_pre, nbins_stim, nbins_post;
int pre, dur, first_dur, stim_start;
double *hist_pre, *hist_stim, *hist_post, *hist_power;
double spike_val;
double spike_sum, spike_sum2, spike_trial;
double mean_pow_dprime;
float time_base;
float time_spike;
char **song_name;
int *syll_beg, *syll_end;
int nlines;
int songind;
int stim_flg=0;   /* if stim flag then make a histogram from stimulus in power units */
int fix_flg=0;    /* if fix_flg is on then the begining of stim is adjusted for silences */
int dp_flg=0;     /* if dp_flg is on then sum of rates, sum square of rates, N, are 
						 * calculated for dprime calculations. dp_flg is then 1 for first, 
						 * 2 for rest and 3 for all */
int lat=10;       /* Latency in ms */
   
   argc --;
   argv ++;

   /* Deal with arguments */
   while ( argv[0][0] == '-' )
   {
      if ( strncmp(argv[0],"-stim",5) == 0 )
      {
         stim_flg = 1;
      }
      else if ( strncmp(argv[0], "-dprime", 7) == 0 )
		{
			argc --;
			argv ++;
			if ( strncmp(argv[0], "first", 5) == 0 )
			{
				dp_flg = 1;
				fix_flg = 1;
			}
			else if ( strncmp(argv[0], "rest", 4) == 0 )
			{
				dp_flg = 2;
				fix_flg = 1;
			}
			else if ( strncmp(argv[0], "all", 3) == 0 )
			{
				dp_flg = 3;
			}
			else
			{
				printf("Error: Unknown modifier for dprime flag: %s\n", argv[0]);
				printf("Possible values are first, rest or all\n");
				exit(1);
			}
		}
		else if ( strncmp(argv[0], "-fixstart", 8) == 0 )
		{
			fix_flg = 1;
		}
      else if ( strncmp(argv[0], "-lat", 7) == 0 )
		{
			argc --;
			argv ++;
			lat = atoi(argv[0]);
		}
		else
      {
         printf("Error: Unknown flag %s\n", argv[0]);
         exit(1);
      }
      argc --;
      argv ++;
   }

	if ( argc != 2 )
	{
		printf("Error: dcp_hist takes at least two arguments\n");
		printf("dcp_hist <-stim> <-fixstart> <-dprime (first, rest or last)> <time bin in ms>  <dcp_file>\n");
		exit(1);
	}


	/* Read bin width in ms */
	tbin = atof(argv[0]);
	if ( tbin <= 0 )
	{
		printf("Error: binwidth of %f is not a real number\n", tbin);
		exit(1);
	}

	/* Read in fix table if fixstart is set or if dprime is set to first or rest */
	if ( fix_flg )
	{
		f_fix = fopen("/auto/fhome/fet/progs/dcp_util/fix.table","r");
		if ( f_fix == NULL ) 
		{
			printf("Error opening fix.table to find beginining of first syllable\n");
			exit(1);
		}
		nlines = 0;
		while (fgets(buff, 256, f_fix)) nlines ++;
		rewind(f_fix);
		song_name = (char **) calloc(nlines, sizeof(char *));
		syll_beg = (int *) calloc(nlines, sizeof(int));
		syll_end = (int *) calloc(nlines, sizeof(int));
		for (i = -1; ++i < nlines; ) 
		{
			song_name[i] = (char *)calloc( 80, sizeof(char));
			fgets(buff, 256, f_fix);
			sscanf(buff, "%s %d %d", song_name[i], syll_beg+i, syll_end+i);
		}
		fclose(f_fix);

	}


   /* Read file */
	pin = read_data(argv[1]);
	if ( pin == NULL ) 
	{
		printf("Error reading dcp file %s\n", argv[2]);
		exit(1);
	}

	if ( pin->ssflag )
	{
		printf("Not yet dealing with spike sorted data\n");
		exit(1);
	}

   pstim=pin->pstim = setup_stimulus(pin->stim_spec);
	if ( pstim == NULL )
	{
		printf("Error setting up stimulus for %s:%s\n", argv[1], pin->stim_spec);
		exit(1);
	}

	/* If fixflag is set or dprime first then check if stimulus is song or seg */
	if ( fix_flg )
	{
		 if ( strcmp(pstim->type, "song") && strcmp(pstim->type, "seg") ) 
		 {
			 printf("Warning: The stimulus is not a song or a segment but fix_flag or dprime flag is set to first or rest\n");
			 printf("Fix flags will be ignored. dprime will be calculated for entire stim");
			 fix_flg = 0;
			 if (dp_flg ) dp_flg =3;
		}
		else
		{
			songind = -1;
			for (i=-1; ++i < nlines; )
			{
				if ( !strcmp(pstim->songfile, song_name[i]) || !strcmp(pstim->segfile, song_name[i]) )
				{
					songind = i;
					break;
				}
			}
			if (songind == -1)
			{
				printf("Song %s or Segment %s was not found in the list of songs\n", pstim->songfile, pstim->segfile);
			   printf("Fix flags will be ignored. dprime will be calculated for entire stim");
			   fix_flg = 0;
			   if (dp_flg ) dp_flg =3;
			}
			else
			{
				 printf("Found song information for d-prime calculation or for fix start\n");
				 printf("Song: %s\n", song_name[i]);
			}
		}
	}

   if ( fix_flg )
	{
	   pre= pin->pre + syll_beg[songind];
	   dur = pstim->dur - syll_beg[songind];
		first_dur = syll_end[songind]-syll_beg[songind];
		stim_start = syll_beg[songind];
	}
	else
	{
		pre= pin->pre;
		dur= pstim->dur;
		stim_start = 0;
	}

	nbins_pre = (int)floor(pre/tbin);
	nbins_stim = (int)floor(dur/tbin);
	nbins_post = (int)floor((pin->cd-pre-dur)/tbin);
	
	hist_pre = (double *)calloc(nbins_pre, sizeof(double));
	hist_stim = (double *)calloc(nbins_stim, sizeof(double));
	hist_post = (double *)calloc(nbins_post, sizeof(double));

   /* Calculate stimulus power if requested */
   if ( stim_flg )
   {
      hist_power = (double *)calloc(nbins_stim, sizeof(double));
      power_stim(pstim, nbins_stim, hist_power, tbin, stim_start, dp_flg, first_dur, &mean_pow_dprime);
		f_stim_power = fopen("power.hist","w");
		if ( f_stim_power == NULL )
		{
			printf("Error: could not open output data file power.hist\n");
			exit(1);
      }
		for ( i = -1; ++i < nbins_stim; ) fprintf(f_stim_power,"%g\n", hist_power[i]);
		fclose(f_stim_power);
      free((void *)hist_power);
   }


	f_out = fopen("out.hist","w");
	if ( f_out == NULL )
	{
		printf("Error: could not open output data file out.hist\n");
		exit(1);
	}
	fprintf(f_out, "%d\t%d\t%d\n", nbins_pre, nbins_stim, nbins_post);
	fclose(f_out);

   /* Make histograms */
	spike_val = 1000/(tbin*pin->n);   /* Make histograms in units of spikes/s */
	spike_sum = 0;
	spike_sum2 = 0;
	for ( i = -1; ++i < pin->n; )
	{
		ptin = pin->ptrial + i;
		spike_trial = 0;
		for ( j = -1; ++ j < ptin->nspikes; )
		{
			time_spike = ptin->time[j]*1.0e-3;
			if ( time_spike < pre )
			{
				ibin = (int)floor(time_spike/tbin);
				if (ibin < nbins_pre) hist_pre[ibin] += spike_val;
			}
			else if (time_spike < pre + dur)
			{
				ibin = (int)floor((time_spike-pre)/tbin);
				if (ibin < nbins_stim) hist_stim[ibin] += spike_val;
				if (dp_flg)
				{
					if (dp_flg == 1 && time_spike > (pre+lat) && time_spike <= (pre+first_dur+lat) )
						spike_trial ++;
					else if (dp_flg == 2 && time_spike > (pre+first_dur+lat) )
						spike_trial ++;
					else if (dp_flg == 3 && time_spike > (pre+lat) )
						spike_trial ++;
				}
			}
			else
			{
				if (dp_flg)
				{
					if (dp_flg == 1 && time_spike > (pre+lat) && time_spike <= (pre+first_dur+lat) )
						spike_trial ++;
					else if (dp_flg == 2 && time_spike > (pre+first_dur+lat) && time_spike <= (pre+dur+lat))
						spike_trial ++;
					else if (dp_flg == 3 && time_spike <= (pre+dur+lat) )
						spike_trial ++;
				}
				ibin = (int)floor((time_spike-pre-dur)/tbin);
				if (ibin < nbins_post) hist_post[ibin] += spike_val;
			}

		}
		if (dp_flg )
		{
			if (dp_flg == 1)
				spike_trial /= (double)first_dur/1000.0;
			else if (dp_flg == 2)
				spike_trial /= (double)(dur-first_dur)/1000.0;
			else
				spike_trial /= (double)dur/1000.0;
			spike_sum += spike_trial;
			spike_sum2 += spike_trial*spike_trial;
		}
	}

   /* Write out histograms */
	f_pre = fopen("pre.hist","w");
	if ( f_pre == NULL )
	{
		printf("Error could not open output file pre.hist\n");
		exit(1);
	}
	for ( i = -1; ++i < nbins_pre; ) fprintf(f_pre,"%g\n", hist_pre[i]);
	fclose(f_pre);

	f_stim = fopen("stim.hist","w");
	if ( f_stim == NULL )
	{
		printf("Error could not open output file stim.hist\n");
		exit(1);
	}
	for ( i = -1; ++i < nbins_stim; ) fprintf(f_stim,"%g\n", hist_stim[i]);
	fclose(f_stim);

	f_post = fopen("post.hist","w");
	if ( f_post == NULL )
	{
		printf("Error could not open output file post.hist\n");
		exit(1);
	}
	for ( i = -1; ++i < nbins_post; ) fprintf(f_post,"%g\n", hist_post[i]);
	fclose(f_post);

	/* Write out sum and sum square values for dprime calculation */
	if (dp_flg)
	{
		f_dp = fopen("dprime.hist", "w");
		if ( f_dp == NULL )
		{
			printf("Error could not open output file dprime.hist\n");
			exit(1);
		}
		fprintf(f_dp,"%d %g %g %g\n", pin->n, spike_sum, spike_sum2, mean_pow_dprime);
		fclose(f_dp);
	}


}
/* ------------------------------- power_stim() ------------------------------*/
power_stim(pstim, nbins, power, tbin, start_time, dp_flg, first_dur, mean_pow)
DCP_STIM *pstim;
int nbins;
double *power;
float tbin;
int start_time;  /* Actual begining of stimulus in ms */
int dp_flg;    /* dp_flag: 0 none, 1 first , 2 rest, 3 all */
int first_dur; /* duration of first syllable in ms */
double *mean_pow;  /* returns the mean power for dprime calculation */
{
int it, ibin;
double xbin, msbin, val;
double time_ms;
int start_it;
int count_mean;

    xbin = 1000.0/(pstim->samprate*tbin);
	 msbin = 1000.0/pstim->samprate;
	 start_it = (int)(start_time*pstim->samprate/1000.0);
	 *mean_pow = 0;
	 count_mean = 0;
    for (it = start_it-1; ++it < pstim->length; )
    {
        ibin = (int)floor((it-start_it)*xbin);
		  time_ms = (it-start_it)*msbin;
        if (ibin >= nbins ) break;
        val =  pstim->stimulus[it]*pstim->stimulus[it];
        power[ibin] += val;
		  if (dp_flg == 1 && time_ms < first_dur )
		  {
			  count_mean ++;
			  *mean_pow += val;
		  }
		  else if (dp_flg == 2 && time_ms > first_dur)
		  {
			  count_mean ++;
			  *mean_pow += val;
			}  
			else if (dp_flg == 3)
			{
			  count_mean ++;
			  *mean_pow += val;
			}
    }
	 *mean_pow /= (double)count_mean;
    for (it = -1; ++it < nbins; ) power[it] *= xbin/32768;
}

#include <stdio.h>
#include <malloc.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>


#include "dcp_ft.h"

main(argc, argv)
int argc;
char *argv[];
{
int i, j;
int it, ib, st, iJN, irec;
int bin0;
char stim_spec[MAX_SPEC_LEN];
DCP_DATA *pin;
DCP_TRIAL *ptin;
DCP_DATA *read_data();
double time_spike;
double **stim_env;       /* Set of stimulus envelopes in ms units */
double *stim_avg;        /* Average Set of stimulus envelopes in ms units */
double **avg_env, ***avg_env_JN;
double stim_val;
double **psth, avg_psth;
int *argvlen, nlen, tot_len;
char **argvrec;
int file_trials=0, *file_trials_JN;
int nJN, nband, nrec;
int nspikes=0, nfiles=0;
int *stim_ns, **stim_ns_JN;
int lin_flg, sil_flg, end_flg = 0;
int sil_window, end_window=0;
FILE *f_count, *f_spike, *f_check, *f_JN, *f_avg, *f_spect;
char fname[100];
char cbuff[1000];

int TWINDOW; 
float f_low;
float f_high;
float f_width;
float f_step;

   argc --;
	argv ++;

	while ( argc && argv[0][0] == '-' )
	{
		if ( strncmp(argv[0],"-end",4) == 0 )
		{
			end_flg = 1;
         argc --;
         argv ++;
         end_window = atoi(argv[0]);
      }
		else 
		{
			printf("Error: Unknown flag %s\n", argv[0]);
			exit(1);
		}
		argc --;
		argv ++;
	}

	if ( argc != 0 )
	{
		printf("Error: dcp_stim_spike takes no argument\n");
		exit(1);
	}

   /* Read Parameters */
	f_count = fopen("stim_init_count.avg","r");
   if ( f_count == NULL )
   {
      printf("Error: cannot open parameter file stim_init_count.avg\n");
      exit(1);
   }
   fscanf(f_count,"%d %d %d %f %f %f %d %d %f", &nrec, &nband, &TWINDOW, &f_low, &f_high, &f_step, &lin_flg, &sil_window, &f_width);
	printf("Starting dcp_stim_spike. Parameters:\n");
   printf("nrec=%d nband=%d TWINDOW=%d f_low=%g f_high=%g f_step=%g f_width=%g lin_flg=%d sil_window=%d\n",
          nrec, nband, TWINDOW, f_low, f_high, f_step, f_width, lin_flg, sil_window);
   fclose(f_count);

   /* Read the files from dcp_stim_init */
	f_check = fopen("stim_init.rec","r");
	if ( f_check == NULL )
	{
		 printf("Could not read rec file from stim_init: stim_int.rec\n");
		 exit(1);
	}

   argvrec = (char **)calloc(nrec, sizeof(char *));
   for ( irec=-1; ++irec < nrec; ) argvrec[irec]=(char *)calloc(256, sizeof(char));
   argvlen = (int *)calloc(nrec, sizeof(int));

   irec = 0;
	printf("Records to be processed:\n");
   while ( fgets( cbuff, 1000, f_check) ) 
   {
     sscanf(cbuff,"%s %d",argvrec[irec], argvlen+irec);
	  printf("Record %d: %s %d\n", irec, argvrec[irec], argvlen[irec]);
     irec ++;
   }
   if ( irec != nrec )
   { 
      printf("Error: number of lines in stim_init.rec does not match nrec\n");
      exit(1);
   }
   printf("done\n");
   fclose(f_check);

   /* Read and allocate stim_avg */
   f_avg = fopen("stim_init_avg.avg","r");
   if ( f_avg == NULL )
   {
      printf("Error: cannot open average stim file stim_init_avg.avg\n");
      exit(1);
   }
   stim_avg = (double *)calloc(nband, sizeof(double));
   ib=0;
   while ( fgets(cbuff, 1000, f_avg) )
   {
     stim_avg[ib] = atof(cbuff);
     ib++;
   }
   if ( ib != nband )
   {
      printf ("Error: missmatch between number of bands and averages\n");
      exit(1);
   }
   fclose(f_avg);


   /* Initialize memory  */
	stim_ns = (int *)calloc(2*TWINDOW+1, sizeof(int));
	nJN = nrec;
	stim_ns_JN = (int **)calloc(nJN, sizeof(int *));
	avg_env_JN = (double ***)calloc(nJN, sizeof(double **));
	file_trials_JN = (int *)calloc(nJN, sizeof(int));
	for ( i = -1; ++i < nJN; )
		stim_ns_JN[i] = (int *)calloc(2*TWINDOW+1, sizeof(int));
	avg_env = (double **)calloc(nband, sizeof(double *));
	for ( ib = -1; ++ib < nband; ) 
		avg_env[ib] = (double *)calloc(2*TWINDOW+1, sizeof(double));
	for ( iJN = -1; ++iJN < nJN; )
	{
		avg_env_JN[iJN] =  (double **)calloc(nband, sizeof(double *));
		for ( ib = -1; ++ib < nband; ) 
			avg_env_JN[iJN][ib] =  (double *)calloc(2*TWINDOW+1, sizeof(double));
	}


   /* Calculate the psth for all stim  */
	psth = (double **)calloc(nrec, sizeof(double *));
	for ( irec = -1; ++irec < nrec; )
	{
		pin = read_data(argvrec[irec]);
		if ( pin == NULL )
		{
         printf("Error reading dcp file %s\n", argvrec[irec]);
         exit(1);
      }
      nlen = argvlen[irec];
		psth[irec]=(double *)calloc(nlen+end_window, sizeof(double) );
		for ( i = -1; ++i < pin->n; )
		{
			ptin = pin->ptrial + i;
			for ( j = -1; ++j < ptin->nspikes; )
			{
				time_spike = ptin->time[j]*1.0e-3;
            bin0 = nint(time_spike) - pin->pre + sil_window;
            if ( bin0 < 0 ) continue;
            if (bin0 >= nlen+end_window ) break;
            psth[irec][bin0] += 1.0;
				nspikes ++;
         }
		}
		for ( j = -1; ++j < nlen+end_window; ) psth[irec][j] /= (double)pin->n;
		free_dcp_memory(&pin);
	}
	avg_psth = 0.0;
	tot_len = 0;
	for ( irec = -1; ++irec < nrec; ) 
	{
		 nlen = argvlen[irec];
		 for ( j = -1; ++j < nlen+end_window; ) avg_psth += psth[irec][j];
		 tot_len += nlen+end_window;
	}
	avg_psth /= (double)tot_len;

	printf("Stim-Spike cross-correlation is based on %d spikes\n", nspikes);
	printf("\tMean Firing rate = %f (spikes/s) Total stim durarion = %f (s)\n", avg_psth*1000.0, tot_len/1000.0);



   /* Calculate cross-correlation between psth and stim */
	for ( irec = -1; ++irec < nrec; )
	{

      /* Allocate memory for spectrogram and read it */
      sprintf(fname,"stim_spect%d.dat",irec);
      f_spect = fopen(fname, "r");
      stim_env = (double **)calloc(nband, sizeof(double *));
      nlen = argvlen[irec];
      for ( ib = -1; ++ib < nband; )
      {
         stim_env[ib] = (double *)calloc(nlen, sizeof(double));
         fread(stim_env[ib], sizeof(double), nlen, f_spect);
      }
      fclose(f_spect);

      iJN = irec;

		for ( j = -1; ++j < nlen + end_window; )
		{
			for ( it = -1 ; ++it <= 2*TWINDOW; )
			{
				st = j + it - TWINDOW;
				if ( st < 0) continue;
				else if ( st >= nlen ) break; 
				stim_ns_JN[iJN][it] ++;
				stim_ns[it] ++;
				for ( ib = -1; ++ib < nband; )
				{
					stim_val = (stim_env[ib][st]-stim_avg[ib])*(psth[irec][j]-avg_psth);
					avg_env[ib][it] += stim_val;
					avg_env_JN[iJN][ib][it] += stim_val;
				}
			}
		}
		free_env(nlen,nband-1,stim_env);
		printf("Done with record %d\n",irec);
	}

				 
	/* Calculate JN estimates */
	for ( iJN = -1; ++iJN < nJN; )
	{
		sprintf(fname,"stim_spike_JN%d.avg", iJN+1);
		f_JN = fopen(fname,"w");
		if ( f_JN == NULL )
		{
			printf("Error opening file %\n",fname);
			exit(1);
		}
		for (it = -1; ++it <= 2*TWINDOW; ) stim_ns_JN[iJN][it] = stim_ns[it] - stim_ns_JN[iJN][it];
		for ( ib = -1; ++ib < nband; )
		{
			for (it = -1; ++it <= 2*TWINDOW; )
			{
				avg_env_JN[iJN][ib][it] = avg_env[ib][it] - avg_env_JN[iJN][ib][it];
				if ( stim_ns_JN[iJN][it] > 0 )
					avg_env_JN[iJN][ib][it] /= (double)stim_ns_JN[iJN][it];
			}
			fwrite(avg_env_JN[iJN][ib], sizeof(double), 2*TWINDOW+1, f_JN);
		}
		fclose(f_JN);
	}


	f_spike = fopen("stim_spike.avg","w");
	for ( ib = -1; ++ib < nband; )
	{
		for (it = -1; ++it <= 2*TWINDOW; )
		{
			if ( stim_ns[it] > 0 )
				avg_env[ib][it] /= (double)stim_ns[it];
		}
		fwrite(avg_env[ib], sizeof(double), 2*TWINDOW+1, f_spike);
	}
	fclose(f_spike);
	printf("dcp_stim_spike is finished\n");

}

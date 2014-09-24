#include <stdio.h>
#include <malloc.h>
#include <string.h>
#include <math.h>
#include <nrutil.h>
#include <nr2.h>
#include <stdlib.h>


#include "dcp_ft.h"

main(argc, argv)
int argc;
char *argv[];
{
int it, ib, st, irec;
int i, j, k;
int bin0;
char stim_spec[MAX_SPEC_LEN];
int TWINDOW;  
DCP_DATA *pin;
DCP_TRIAL *ptin;
DCP_DATA *read_data();
double time_spike, time_spike2;
double *cross_spike;
int nband, nlen, trial_spikes, nrec;
int nspikes=0;
int sil_flg, lin_flg;
int sil_window;
float f_low, f_high, f_step, f_width;
int *spike_ns;
FILE *f_spike, *f_check, *f_count;
char **argvrec;
int *argvlen;
int file_trials=0;
int spike_count_avg=0;
char cbuff[1000];

   argc --;
	argv ++;

	if ( argc != 0 )
	{
		printf("Error: dcp_spike does not take any arguments\n");
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
   while ( fgets( cbuff, 1000, f_check) ) 
   {
     sscanf(cbuff,"%s %d",argvrec[irec], argvlen+irec);
	  printf("Line %d: %s %d\n", irec, argvrec[irec], argvlen[irec]);
     irec ++;
   }
   if ( irec != nrec )
   { 
      printf("Error: number of lines in stim_init.rec does not match nrec\n");
      exit(1);
   }
   printf("done\n");
   fclose(f_check);


	cross_spike = (double *)calloc(2*TWINDOW+1, sizeof(double));
	spike_ns = (int *)calloc(2*TWINDOW+1, sizeof(int));

   /* Calculate auto-correlation of spikes */
	for ( irec=-1; ++irec < nrec; )
	{
		pin = read_data(argvrec[irec]);
		nlen = argvlen[irec];
		if ( pin == NULL )
		{
			printf("Error reading dcp file %s\n", argvrec[irec]);
			exit(1);
		}
		if ( pin->ssflag )
		{
			printf("Not yet dealing with spike sorted data\n");
			exit(1);
		}

		for ( i = -1; ++i < pin->n; )
		{
			ptin = pin->ptrial + i;
			trial_spikes = 0;
			for ( j = -1; ++j < ptin->nspikes; )
			{
				time_spike = ptin->time[j]*1.0e-3;
				if ( time_spike < pin->pre - sil_window ) continue;
				if ( time_spike > pin->pre + nlen - sil_window ) break;
				bin0 = nint(time_spike) - pin->pre - sil_window;
				nspikes ++;
				trial_spikes ++;
#ifdef BIASED
				for ( it = -1; ++it <= 2*TWINDOW; )
				{
					st = bin0 + it;
					if ( st < -sil_window || (st + pin->pre) < 0 ) continue;
					if ( st > nlen - sil_window || st >= (pin->cd-pin->pre) ) break;
					spike_ns[it] ++;
				}
#endif
				for ( k = -1; ++k < ptin->nspikes; )
				{
					time_spike2 = ptin->time[k]*1.0e-3;
					if ( time_spike2 < pin->pre - sil_window ) continue;
					if ( time_spike2 > pin->pre + nlen - sil_window ) break;
					st = nint(time_spike2 - time_spike) + TWINDOW;
					if ( st < 0 ) continue;
					if ( st > 2*TWINDOW ) break;
					cross_spike[st] += 1.0;
				}
			}
		}
		file_trials += pin->n;
      spike_count_avg += nlen*pin->n;
		free_dcp_memory(&pin);
	}
				 
	f_spike = fopen("spike.avg", "w");
	for ( it = -1; ++ it <= 2*TWINDOW; )
	{
		if ( it < TWINDOW )
			spike_ns[it] = spike_count_avg - (file_trials)*(TWINDOW-it);
		else
			spike_ns[it] = spike_count_avg - (file_trials)*(it-TWINDOW);

		cross_spike[it] /= (double)spike_ns[it];
		fprintf(f_spike, "%g\n", cross_spike[it]);
	}
	fclose(f_spike);

}

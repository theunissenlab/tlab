#include <stdio.h>
#include <malloc.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>


#include "dcp_ft.h"
DCP_DATA *read_raw();

main(argc, argv)
int argc;
char *argv[];
{
int it1, it2, ib1, ib2, st, xb, irec;
char stim_spec[MAX_SPEC_LEN];
DCP_DATA *pin;
DCP_STIM *stim;
DCP_DATA *read_data();
double **stim_env;        /* Set of stimulus envelopes in ms units */
double **env_corr;
double ***JN_env_corr;
double *stim_avg;
double stim1, stim2;
float f_low, f_high, f_step, f_width;
float amp_samp;
char **argvrec;
int *argvlen;
int lin_flg;
int dc_flg;
int **ns_env_corr;
int ***JN_ns_env_corr;
int ns;
int nrec;
double *tot_corr;
int nband, nlen, ncorr;
int sil_window;
FILE *f_spect;
FILE *f_stim, *f_check, *f_count, *f_avg;
int TWINDOW;
int twindow_len;
char cbuff[1000];
char fname[64];


   argc --;
	argv ++;

   dc_flg = 0;
   while ( argc != 0 && argv[0][0] == '-' )
   {
      if ( strncmp(argv[0],"-dc",4) == 0 )
      {
         dc_flg = 1;
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
		printf("Error: dcp_stim_stat does not take any arguments\n");
		exit(1);
	}


   /* Read Parameters */
	f_count = fopen("stim_init_count.avg","r");
   if ( f_count == NULL )
   {
      printf("Error: cannot open parameter file stim_init_count.avg\n");
      exit(1);
   }
   fscanf(f_count,"%d %d %d %f %f %f %d %d %f %f", &nrec, &nband, &TWINDOW, &f_low, &f_high, &f_step, &lin_flg, &sil_window, &f_width, &amp_samp);
	twindow_len = (int)ceil(TWINDOW*amp_samp/1000);
	printf("Starting dcp_stim. Parameters:\n");
   printf("nrec=%d nband=%d TWINDOW=%d (%d pts) f_low=%g f_high=%g f_step=%g f_width=%g lin_flg=%d sil_window=%d amp_samp=%f\n",
          nrec, nband, TWINDOW, twindow_len, f_low, f_high, f_step, f_width, lin_flg, sil_window, amp_samp);
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
	printf("Records to be processed\n");
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
   fclose(f_check);

   /* Read and allocate stim_avg */
   f_avg = fopen("stim_init_avg.avg","r");
   if ( f_avg == NULL )
   {
      printf("Error: cannot open average stim file stim_init_avg.avg\n");
      exit(1);
   }
   stim_avg = (double *)calloc(nband, sizeof(double));
   ib1=0;
   while ( fgets(cbuff, 1000, f_avg) )
   {
     stim_avg[ib1] = atof(cbuff);
     printf("Band %d stim_avg = %g\n", ib1, stim_avg[ib1]);
     ib1++;
   }
   if ( ib1 != nband )
   {
      printf ("Error: missmatch between number of bands and averages\n");
      exit(1);
   }
   fclose(f_avg);


   /* Allocate memory for the cross-correlation */
	tot_corr = (double *)calloc(2*twindow_len+1, sizeof(double));
	ncorr = nband;
	JN_env_corr = (double ***)calloc(nrec, sizeof(double **));
	JN_ns_env_corr = (int ***)calloc(nrec, sizeof(int **));
	env_corr = (double **)calloc(ncorr, sizeof(double *));
	ns_env_corr = (int **)calloc(ncorr, sizeof(int *));
	for ( irec = -1; ++irec < nrec; ) 
	{
		JN_env_corr[irec] = (double **)calloc(ncorr, sizeof(double *));
		JN_ns_env_corr[irec] = (int **)calloc(ncorr, sizeof(int *));
	}

	for ( ib1 = -1; ++ib1 < ncorr; ) 
	{
		env_corr[ib1] = (double *)calloc(2*twindow_len+1, sizeof(double));
		ns_env_corr[ib1] = (int *)calloc(2*twindow_len+1, sizeof(int));
		for ( irec = -1; ++irec < nrec; ) 
		{
			JN_env_corr[irec][ib1] = (double *)calloc(2*twindow_len+1, sizeof(double));
			JN_ns_env_corr[irec][ib1] = (int *)calloc(2*twindow_len+1, sizeof(int));
		}
	}

   printf("Starting cross-correlation calculation ...\n");
	for ( irec = -1; ++irec < nrec; )
	{
		pin = read_data(argvrec[irec]);
      if ( pin == NULL )
      {
         printf("Error reading dcp file %s\n", argvrec[irec]);
         exit(1);
      }

      /* Allocate memory for spectrogram and read it */
		sprintf(fname,"stim_spect%d.dat",irec);
		f_spect = fopen(fname, "r");
		stim_env = (double **)calloc(nband, sizeof(double *));
      nlen = argvlen[irec];
		for ( ib1 = -1; ++ib1 < nband; ) 
      {
			stim_env[ib1] = (double *)calloc(nlen, sizeof(double));
         fread(stim_env[ib1], sizeof(double), nlen, f_spect);
      }
      fclose(f_spect);

		for ( ib1 = -1; ++ib1 < nband; )
		{
			for ( ib2 = ib1-1; ++ib2 < nband; )
			{
				xb = ib2-ib1;
				for ( it1 = -1; ++ it1 < nlen; )
				{
					stim1 = stim_env[ib1][it1];
					if ( dc_flg == 0 ) stim1 -= stim_avg[ib1];

					for ( it2 = it1-twindow_len-1; ++it2 <= it1+twindow_len; )
					{  
						if ( it2 < 0 ) continue;
						else if ( it2 < nlen) 
                  {
							stim2 = stim_env[ib2][it2];
							if ( dc_flg == 0 ) stim2 -= stim_avg[ib2];
						}
						else break;

						st = it2 - it1 + twindow_len;
						env_corr[xb][st] += (stim2*stim1)*pin->n;
						ns_env_corr[xb][st] += pin->n;
						JN_env_corr[irec][xb][st] += (stim2*stim1)*pin->n;
						JN_ns_env_corr[irec][xb][st] += pin->n;

						if ( stim1*stim2 < -1.0e20 || stim1*stim2 > 1.0e20 )
						{
							printf("Warning: nb1 = %d nb2 = %d it1 = %d it2 = %d\n",
							ib1,ib2, it1, it2);
							printf("stim1*stim2 is out of ordinary: stim1 = %g stim2 = %g stim_avg1 = %g stim_avg2 = %g \n", stim1, stim2, stim_avg[ib1], stim_avg[ib2]);
						}

					}
				}
			}
		}

		free_dcp_memory(&pin);
		free_env(nlen,nband-1,stim_env);
		printf("Done with record %d\n", irec);
	}

	printf("Writing output files to disk\n", irec);
   /* Write out results for JN files - these are the delete one values */
	for ( irec = -1; ++ irec < nrec; )
	{
		sprintf(fname, "stim_stat%d.avg", irec);
		f_stim = fopen(fname,"w");
		if ( f_stim == NULL )
		{
			printf("Error: Could not write output file %s\n", fname);
			exit(1);
		}

		for (xb = -1; ++xb < ncorr; )
		{
			for ( st = -1; ++st <= 2*twindow_len; )
			{
				ns = ns_env_corr[xb][st] - JN_ns_env_corr[irec][xb][st];
				if ( ns != 0 )
					JN_env_corr[irec][xb][st] = (env_corr[xb][st] - JN_env_corr[irec][xb][st])/(double)ns;
				else
					JN_env_corr[irec][xb][st] = 0.0;

			}
			fwrite(JN_env_corr[irec][xb] , sizeof(double), 2*twindow_len+1, f_stim);
		}
		fclose(f_stim);
	}

   /* Now write the overall mean */
   f_stim = fopen("stim_stat.avg","w");
	if ( f_stim == NULL )
	{
		printf("Error: Could not write output file stim_stat.avg\n");
		exit(1);
	}

	for (xb = -1; ++xb < ncorr; )
	{
		for ( st = -1; ++st <= 2*twindow_len; )
		{
			if ( ns_env_corr[xb][st] != 0 )
				env_corr[xb][st] /= (double)ns_env_corr[xb][st];
		}
		fwrite(env_corr[xb] , sizeof(double), 2*twindow_len+1, f_stim);
	}
	fclose(f_stim);


	printf("dcp_stim_stat is done\n");

				 
}

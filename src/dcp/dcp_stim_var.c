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
int it1, ib1, irec;
char stim_spec[MAX_SPEC_LEN];
DCP_DATA *pin;
DCP_STIM *stim;
DCP_DATA *read_data();
double **stim_env;        /* Set of stimulus envelopes in ms units */
double *stim_avg;
double *stim_var;
double stim1;
float f_low, f_high, f_width, f_step;
char **argvrec;
int *argvlen;
int lin_flg;
long ns_var;
int nrec;
double *tot_corr;
int nband, nlen;
int sil_window;
FILE *f_spect;
FILE *f_stim, *f_check, *f_count, *f_avg;
int TWINDOW;
char cbuff[1000];
char fname[64];


   argc --;
	argv ++;

	if ( argc != 0 )
	{
		printf("Error: dcp_stim does not take any arguments\n");
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
	printf("Starting dcp_stim_var. Parameters:\n");
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


   /* Allocate memory for the variance */
	stim_var = (double *)calloc(nband, sizeof(double));
	ns_var = 0;

   printf("Starting variance calculation...\n");
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

		ns_var += pin->n*nlen;
		for ( ib1 = -1; ++ib1 < nband; )
		{
			for ( it1 = -1; ++ it1 < nlen; )
			{
				stim1 = stim_env[ib1][it1]-stim_avg[ib1];
				stim_var[ib1] += (stim1*stim1)*pin->n;
			}
		}

		free_dcp_memory(&pin);
		free_env(nlen,nband-1,stim_env);
		printf("Done with record %d\n", irec);
	}

   f_stim = fopen("stim_var.avg","w");
	if ( f_stim == NULL )
	{
		printf("Error: Could not write output file stim_var.avg\n");
		exit(1);
	}

	for ( ib1 = -1; ++ib1 < nband; )
	{
		if ( ns_var != 0 )
				stim_var[ib1] /= (double)ns_var;
		fprintf(f_stim, "%g\n", stim_var[ib1]);
	}
	fclose(f_stim);
	printf("dcp_stim_var is done\n");

				 
}

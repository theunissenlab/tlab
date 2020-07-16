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
int it, ib;
int i, j, k;
int bin0;
char stim_spec[MAX_SPEC_LEN];
int nt, nf;
DCP_DATA *pin;
DCP_TRIAL *ptin;
DCP_STIM *stim;
DCP_DATA *read_data();
double time_spike;
double **stim_env;       /* Set of stimulus envelopes in ms units */
int *spike_trials;
int  nband, nlen ;
char fname[100];
FILE *f_stim, *f_spike, *f_count;

/* These are the default parameters */
float f_low=250.0, f_high=8000.0, f_width=250.0;
int TWINDOW=0;

	if ( argc < 2 )
	{
		printf("Error: dcp_show takes at least 1 argument\n");
		printf("dcp_show  <file1> <file2> ... \n");
		exit(1);
	}

   argc --;
	argv ++;
   
	while ( argv[0][0] == '-' )
	{
	   if ( strncmp(argv[0], "-window", 7) == 0 )
      {
         argc --;
         argv ++;
         TWINDOW = atoi(argv[0]);
      }
      else if ( strncmp(argv[0], "-f_high", 6) == 0 )
      {
         argc --;
         argv ++;
         f_high = atof(argv[0]);
      }
      else if ( strncmp(argv[0], "-f_low", 5) == 0 )
      {
         argc --;
         argv ++;
         f_low = atof(argv[0]);
      }
      else if ( strncmp(argv[0], "-f_width", 7) == 0 )
      {
         argc --;
         argv ++;
         f_width = atof(argv[0]);
      }
		else 
		{
			printf("Error: Unknown flag %s\n", argv[0]);
			exit(1);
		}
		argc --;
		argv ++;
	}

	/* Open output files */
	sprintf(fname,"spike_all.dat",nf);
	f_spike = fopen(fname,"w");
	sprintf(fname,"stim_all.dat",nf);
	f_stim = fopen(fname,"w");
	sprintf(fname,"count_all.dat",nf);
	f_count = fopen(fname,"w");
	

   nf = 0;
   while (argc)
	{
		pin = read_data(argv[0]);
		nf ++;
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


		/* Pre-process the stimulus */
		stim_fband(stim,&nband,&nlen,TWINDOW,&stim_env,&f_low,&f_high,&f_width);


		spike_trials = (int *)calloc(nlen+2*TWINDOW, sizeof(int));
		for ( it = -1; ++it < nlen+2*TWINDOW; )
		{  
			fprintf(f_stim,"%g ", stim_env[nband][it]/(double)nband);
		}

		for ( i = -1; ++i < pin->n; )
		{
			ptin = pin->ptrial + i;

			for ( j = -1; ++j < ptin->nspikes; )
			{
				time_spike = ptin->time[j]*1.0e-3;
				if ( time_spike < pin->pre - TWINDOW ) continue;
				if ( time_spike > pin->pre + stim->dur + TWINDOW ) break;
				bin0 = nint(time_spike) - pin->pre + TWINDOW ;
				spike_trials[bin0] += 1;
			}
		}
		for ( it = -1; ++it < nlen+2*TWINDOW; ) 
			fprintf(f_spike,"%d ", spike_trials[it]);
		free(spike_trials);

		free_dcp_memory(&pin);
		free_env(nlen+2*TWINDOW,nband,stim_env);
		argc --;
		argv ++;
	}

   /* Write output to files */
	fprintf(f_count,"%d %d %d %f %f %f\n", nband, TWINDOW, nlen, f_low, f_high, f_width);
	fclose(f_count);
	fclose(f_spike);
	fclose(f_stim);

}

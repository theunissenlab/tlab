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
float f_step;
int TWINDOW=0;
int lin_flg=0;

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
      if ( strncmp(argv[0],"-lin",4) == 0 )
      {
         lin_flg = 1;
      }
	   else if ( strncmp(argv[0], "-window", 7) == 0 )
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

		/* Open output files */
		sprintf(fname,"spike%d.dat",nf);
		f_spike = fopen(fname,"w");
		sprintf(fname,"stim%d.dat",nf);
		f_stim = fopen(fname,"w");
		sprintf(fname,"count%d.dat",nf);
		f_count = fopen(fname,"w");
      

		/* Pre-process the stimulus */
		f_step = f_width;
		stim_fband(stim,&nband,&nlen,TWINDOW,&stim_env,&f_low,&f_high,&f_width,&f_step);

      fprintf(f_count,"%d %d %d %f %f %f\n", nband, TWINDOW, nlen, f_low, f_high, f_width);
		fclose(f_count);

      /* Take log of enveloppes */
      if ( lin_flg == 0 )
      {
			for ( ib = -1; ++ib < nband; )
				for ( it = -1; ++it < nlen+2*TWINDOW; )
					stim_env[ib][it] = log(stim_env[ib][it]+1.0);
		}

		spike_trials = (int *)calloc(nlen+2*TWINDOW, sizeof(int));
		for ( ib = -1; ++ib < nband; )
		{
			for ( it = -1; ++it < nlen+2*TWINDOW; )
			{  
				fprintf(f_stim,"%g ", stim_env[ib][it]);
			}
			fprintf(f_stim,"\n");
		}
		fclose(f_stim);
		printf("Stim dur = %d\n", stim->dur);

		if ( pin->n == 0 )
		{
		   printf("Warning: no trials in %s\n", argv[0]);
		}

		for ( i = -1; ++i < pin->n; )
		{
			ptin = pin->ptrial + i;
			for ( it = -1; ++it < nlen+2*TWINDOW; ) spike_trials[it] = 0;

			for ( j = -1; ++j < ptin->nspikes; )
			{
				time_spike = ptin->time[j]*1.0e-3;
				if ( time_spike < pin->pre - TWINDOW ) continue;
				if ( time_spike >= pin->pre + stim->dur + TWINDOW ) break;
				bin0 = nint(time_spike) - pin->pre + TWINDOW ;
				if (bin0 < 0 || bin0 >= nlen+2*TWINDOW ) continue;
				spike_trials[bin0]++;
			}
			for ( it = -1; ++it < nlen+2*TWINDOW; ) 
				fprintf(f_spike,"%d ", spike_trials[it]);
			fprintf(f_spike,"\n");
		}
		free(spike_trials);
		fclose(f_spike);

		free_dcp_memory(&pin);
		free_env(nlen+2*TWINDOW,nband,stim_env);
		argc --;
		argv ++;
	}


}

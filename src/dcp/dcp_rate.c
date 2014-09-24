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
DCP_DATA *pin;
DCP_DATA *read_data();
DCP_STIM *pstim;
DCP_TRIAL *ptin;
FILE *f_rate;
int pre, dur;
double *rate_trial;
double time_spike;
int spike_trial;
   
   argc --;
   argv ++;

   /* Deal with arguments */
	if ( argc != 1 )
	{
		printf("Error: dcp_rate takes one argument\n");
		printf("dcp_rate <dcp_file>\n");
		exit(1);
	}


   /* Read file */
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

   pstim=pin->pstim = setup_stimulus(pin->stim_spec);
	if ( pstim == NULL )
	{
		printf("Error setting up stimulus for %s:%s\n", argv[1], pin->stim_spec);
		exit(1);
	}

	pre= pin->pre;
	dur= pstim->dur;

	rate_trial = (double *)calloc(pin->n, sizeof(double));

	f_rate = fopen("rate.trial","w");
	if ( f_rate == NULL )
	{
		printf("Error: could not open output data file rate.trial\n");
		exit(1);
	}

   /* Find rate for each trial */
	for ( i = -1; ++i < pin->n; )
	{
		ptin = pin->ptrial + i;
		spike_trial = 0;
		for ( j = -1; ++ j < ptin->nspikes; )
		{
			time_spike = ptin->time[j]*1.0e-3;
			if (time_spike > pre && time_spike <= pre + dur)
			{
						spike_trial ++;
			}
		}
		rate_trial[i] = (1000.0*spike_trial)/(double)dur;
	}

   /* Write out rate */
	for ( i = -1; ++i < pin->n; )
	{
	   fprintf(f_rate,"%g\n", rate_trial[i]);
	}

	fclose(f_rate);

}

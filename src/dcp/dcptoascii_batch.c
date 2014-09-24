/* dcp_stim_init.c - main file for dcp_stim_init. dcp_stim_init should be called
before dcp_stim, dcp_stim_spike, dcp_forward, etc.  It calculates the 
spectrogram of all the stimuli. It outputs 4 types of files:
stim_init.rec		- has file name + length of stim
stim_spect(file number).dat  - the spectrograms in binary format
stim_init_count.avg - the parameters 
stim_init_avg.avg - the average stim in each band
*/


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
int i, j;
char stim_spec[MAX_SPEC_LEN];
DCP_DATA *pin;
DCP_STIM *stim;
DCP_DATA *read_data();
DCP_TRIAL *ptin;
double stim_len;       
int lin_flg;
int argc_start;
char **argv_start;
char fname[80];
int nf;
FILE *f_check, *f_spike;
float time_spike;
char delims[] = "/";
char *last_name = NULL, *temp_name=NULL;

/* These are the default parameters that can be changed on the command line */

   argc --;
	argv ++;

	/* Deal with arguments */
	lin_flg = 0;
	while ( argv[0][0] == '-' )
	{
		if ( strncmp(argv[0],"-lin",4) == 0 )
		{
			printf("Linear flag ignored\n");
			lin_flg = 1;
		}
		else 
		{
			printf("Error: Unknown flag %s\n", argv[0]);
			exit(1);
		}
		argc --;
		argv ++;
	}
	if ( argc == 0 )
	{
		printf("Error: dcp_spike_avg takes at least 1 argument\n");
		printf("dcptoascii_batch <file1> <file2> ... \n");
		exit(1);
	}


   /* Dump data to ascii files */
	f_check = fopen("stim_init.rec","w");
	nf = 0;
	argc_start = argc;
	argv_start = argv;

   while (argc)
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
		if ( pin->n == 0 )
		{
		   printf("No trials in dcp file %s\n", pin->n);
			exit(1);
		}
		strcpy(stim_spec, pin->stim_spec);
		stim = pin->pstim = setup_stimulus(stim_spec);
		if ( stim == NULL )
		{
			printf("Error setting up stimulus for %s:%s\n", argv[0], pin->stim_spec);
			exit(1);
		}
		stim_len = (1000.0*stim->length)/stim->samprate; /* Stim length in ms */
		fprintf(f_check,"%s %g\n", argv[0], stim_len);
		  
		/* Write the spike data as spike arrival times */
		last_name = argv[0];
		temp_name = strtok(argv[0], delims);
		while( temp_name != NULL ) 
		{
		   last_name = temp_name;
		 	temp_name  = strtok( NULL, delims );
	   } 

		sprintf(fname,"%s.txt", last_name);
		f_spike = fopen(fname,"w");
		if ( f_spike == NULL )
		{
			printf("Error opening spike file %s\n", fname);
			exit(1);
		}

		/* Dump spikes */
		for ( i = -1; ++i < pin->n; )
		{
			ptin = pin->ptrial + i;
			for ( j = -1; ++ j < ptin->nspikes; )
			{
				time_spike = ptin->time[j]*1.0e-3 - pin->pre;  /* time_spike is in milli seconds */
				/* if ( time_spike <=0 ) continue; */
				fprintf(f_spike,"%.3f ", time_spike);
			}
			fprintf(f_spike,"\n");
		}
		fclose(f_spike);
		free_dcp_memory(&pin);
		argc --;
		argv ++;
		nf ++;
	}
	fclose(f_check);


}


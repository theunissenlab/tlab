#include <stdlib.h>
#include <stdio.h>
#include <malloc.h>
#include <string.h>
#include <sndfile.h>

#include "dcp_ft.h"

main(argc, argv)
int argc;
char *argv[];
{
int i,j;
char f_name[256];
char stim_spec[MAX_SPEC_LEN];
DCP_DATA *pin;
DCP_DATA *read_data();
DCP_TRIAL *ptin;
DCP_STIM *stim;
FILE *f_spike;
SNDFILE *f_stim ;
SF_INFO sfinfo ;
		  
float time_spike;

	if ( argc != 3 )
	{
		printf("Error: dcptostrfpak takes two arguments\n");
		printf("dcptostrfpak <dcp_file> <id for output files>\n");
		exit(1);
	}

   /* Read file */
	pin = read_data(argv[1]);
	if ( pin == NULL ) 
	{
		printf("Error reading dcp file %s\n", argv[1]);
		exit(1);
	}
	if ( pin->ssflag )
	{
		printf("Not yet dealing with spike sorted data\n");
		exit(1);
	}

	/* Write the stimulus as a wave file */
	strcpy(stim_spec, pin->stim_spec);
	stim = pin->pstim = setup_stimulus(stim_spec);
   if ( stim == NULL )
   {
       printf("Error setting up stimulus for %s:%s\n", argv[1], pin->stim_spec);
       exit(1);
   }

	/* write the wav file ... */
	memset (&sfinfo, 0, sizeof (sfinfo)) ;

	sfinfo.samplerate  = stim->samprate;
	sfinfo.frames     = stim->length;
	sfinfo.channels    = 1 ;
	sfinfo.format      = (SF_FORMAT_WAV | SF_FORMAT_PCM_16) ;

	sprintf(f_name, "stim%d.wav", atoi(argv[2]));
   if (! (f_stim = sf_open (f_name, SFM_WRITE, &sfinfo)))
   {
		printf ("Error : Not able to open output file %s.\n", f_name) ;
      return 1 ;
   } 

	if (sf_write_short (f_stim, stim->stimulus, stim->length) != stim->length)
	{
		puts (sf_strerror (f_stim)) ;
		printf("Error writing the data to the sound file\n");
	}
	sf_close (f_stim) ;

	/* Write the spike data as spike arrival times */
	sprintf(f_name,"spike%d", atoi(argv[2]));
	f_spike = fopen(f_name,"w");
	if ( f_spike == NULL )
	{
		printf("Error opening spike file %s\n", f_name);
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
}


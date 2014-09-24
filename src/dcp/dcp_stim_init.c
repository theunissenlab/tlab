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
int it, ib;
char stim_spec[MAX_SPEC_LEN];
DCP_DATA *pin;
DCP_STIM *stim;
DCP_DATA *read_data();
double **stim_env;        /* Set of stimulus envelopes in ms units */
double *stim_avg;
int nband, nlen, ncorr;
int init_flg, count_avg, lin_flg;
int argc_start;
char **argv_start;
char fname[80];
int nf;
FILE *f_stim, *f_check, *f_count;

/* These are the default parameters that can be changed on the command line */
float f_low=250.0;
float f_high=8000.0;
float f_width=250.0;
float f_step = -1;
float amp_samp = 1000;
int sil_window=0;
int twindow=200;

   argc --;
	argv ++;

	/* Deal with arguments */
	lin_flg = 0;
	while ( argv[0][0] == '-' )
	{
		if ( strncmp(argv[0],"-lin",4) == 0 )
		{
			lin_flg = 1;
		}
		else if ( strncmp(argv[0],"-sil",4) == 0 )
		{
			argc --;
			argv ++;
			sil_window = atoi(argv[0]);
		}
		else if ( strncmp(argv[0],"-win",4) == 0 )
		{
			argc --;
			argv ++;
			twindow = atoi(argv[0]);
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
		else if ( strncmp(argv[0], "-f_step", 6) == 0 )
		{
			argc --;
			argv ++;
			f_step = atof(argv[0]);
		}
		else if ( strncmp(argv[0], "-amp_samp", 5) == 0 )
		{
			argc --;
			argv ++;
			amp_samp = atof(argv[0]);
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
		printf("dcp_stim <file1> <file2> ... \n");
		exit(1);
	}

	/* Set the default f_step to be equal to f_width. Probably want to change this later... */
	if ( f_step < 0.0 ) f_step = f_width;


   /* Calculate stimulus auto correlations */
	f_check = fopen("stim_init.rec","w");
	init_flg = 1;
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
			printf("Error setting up stimulus for %s:%s\n", argv[1], pin->stim_spec);
			exit(1);
		}

		/* Pre-process the stimulus */
		stim_fband(stim,&nband,&nlen,sil_window,&stim_env,&f_low,&f_high,&f_width,&f_step, amp_samp);

		if ( init_flg )
		{
			stim_avg = (double *)calloc(nband, sizeof(double));
			count_avg = 0;
			init_flg = 0;
		}

		for ( ib = -1; ++ib < nband; )
		{
			for ( it = -1; ++ it < nlen+2*sil_window; )
			{
				if ( lin_flg == 0 )
				{
					if ( stim_env[ib][it] < 0.0 ) 
					{
						printf("Error: negative amplitude enveloppe\n");
						printf("ib=%d it=%d env=%g\n", ib, it, stim_env[ib][it]);
						exit(1);
					}
					stim_env[ib][it] = log(stim_env[ib][it]+1.0);
				}
				stim_avg[ib] += stim_env[ib][it]*pin->n;
			}
		}
		count_avg += (nlen+2*sil_window)*pin->n; 

      /* Write out spectrogram */
		fprintf(f_check,"%s %d\n", argv[0], nlen+2*sil_window);
		sprintf(fname,"stim_spect%d.dat",nf);
		f_stim = fopen(fname,"w");
		for ( ib = -1; ++ib < nband; )
			fwrite(stim_env[ib], sizeof(double), nlen+2*sil_window, f_stim);
		fclose(f_stim);

		free_dcp_memory(&pin);
		free_env(nlen+2*sil_window,nband,stim_env);
		argc --;
		argv ++;
		nf ++;
	}
	fclose(f_check);

   /* Print out parameters */
	f_count = fopen("stim_init_count.avg","w");
   fprintf(f_count,"%d %d %d %f %f %f %d %d %f %f", nf, nband, twindow, f_low, f_high, f_step, lin_flg, sil_window, f_width, amp_samp);
   fclose(f_count);

   /* Print out averages */
	f_count = fopen("stim_init_avg.avg","w");
	for ( ib = -1; ++ib < nband; ) 
	{
		stim_avg[ib] /= count_avg;
		fprintf(f_count,"%g\n", stim_avg[ib]);
	}
	fclose(f_count);

}

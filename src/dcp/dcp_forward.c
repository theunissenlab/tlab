#include <stdio.h>
#include <malloc.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <sys/types.h>
#include <time.h>



#include "dcp_ft.h"


#define MAXTRIALCOUNT 20

main(argc, argv)
int argc;
char *argv[];
{
int it, ib, st, dt, nt, nb;
int i, j, k;
int bin0;
int nrec, irec;
char **argvrec;
int *argvlen;
char stim_spec[MAX_SPEC_LEN];
DCP_DATA *pin;
DCP_TRIAL *ptin;
DCP_DATA *pout;
DCP_TRIAL *ptout;
DCP_STIM *stim;
DCP_DATA *read_data();
double time_spike;
double **stim_env;       /* Set of stimulus envelopes in ms units */
double *stim_avg;        /* Average Set of stimulus envelopes in ms units */
double *stim_avg_test;   /* Average Set of stimulus envelopes in ms units */
double *est_spike, *avg_spike2, *avg_spike1;  /* The psth is divided in two */
double avg_est=0, avg_est2=0;
int totlen=0;
int *spike_trials;
int count_avg, trial_count, max_trial_count, current_max_value;
int  nband, nlen, nlast;
char fname[100];
char cbuff[1000];
float ***forward;
int matchflg;
char **argv_start;
int argc_start, nmax;
int *ntrials_hist;
int *ntrials_value;
FILE *f_in, *f_out_psth1, *f_out_psth2, *f_out_pre, *f_count;
FILE *f_ntrials, *f_check, *f_spect;
FILE *f_avg;
int lin_flg, gen_flg, sil_flg, end_flg;
int sil_window, end_window;
double min_est, max_est, mean_est, rate_factor;
double spike_thres;
int *fake_spike;
double stim_val;
int ns;
float f_low, f_high, f_step, f_width;
int TWINDOW;

#define MAX_JN 10


   argc --;
	argv ++;
   
   matchflg = 1;
	gen_flg = 0;
	end_window = 0;
	srand48((long)time(NULL));

	while ( argc && argv[0][0] == '-' )
	{
		if ( strncmp(argv[0],"-end",4) == 0 ) 
		{
			argc --;
			argv ++;
			end_window = atoi(argv[0]);
		}
		else if ( strncmp(argv[0],"-nomatch",8) == 0 ) 
		{
			matchflg = 0;
		}
		else if ( strncmp(argv[0],"-gen",4) == 0 ) 
		{
			gen_flg = 1;
		}
		else
		{
			printf("Error: Unknown flag %s\n", argv[0]);
			exit(1);
		}
		argc --;
		argv ++;
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
	nt = 2*TWINDOW+1;
   fclose(f_count);

   /* Read the files from dcp_stim_init */
	if ( matchflg )
	{
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
	}
	else
	{
		if ( argc == 0 )
		{
			printf("Error: missing file names. nomatch flag requires arguments.\n");
			exit(1);
		}
		nrec = argc;
		argvrec = (char **)calloc(nrec, sizeof(char *));
		for ( irec=-1; ++irec < nrec; ) 
		{
			argvrec[irec]=(char *)calloc(256, sizeof(char));
		   strcpy( argvrec[irec], argv[irec] );
		   printf("Argument %d: %s\n", irec, argvrec[irec]);
		}
	}


   /* Read and allocate stim_avg */
   f_avg = fopen("stim_init_avg.avg","r");
   if ( f_avg == NULL )
   {
      printf("Error: cannot open average stim file stim_init_avg.avg\n");
      exit(1);
   }
   stim_avg = (double *)calloc(nband, sizeof(double));
   stim_avg_test = (double *)calloc(nband, sizeof(double));
   ib=0;
   while ( fgets(cbuff, 1000, f_avg) )
   {
     stim_avg[ib] = atof(cbuff);
     printf("Band %d stim_avg = %g\n", ib, stim_avg[ib]);
     ib++;
   }
   if ( ib != nband )
   {
      printf ("Error: missmatch between number of bands and averages\n");
      exit(1);
   }
   fclose(f_avg);

	/* Read in forward filter file */
   if ( matchflg )
	{
		forward = (float ***)calloc(nrec,sizeof(float ***));
		for (i=-1; ++i< nrec; )
		{
			sprintf(fname,"forwardJN%d.filt",i);
			f_in = fopen(fname,"r");
			if ( f_in == NULL )
			{
				printf("Error: cannot open filter file %s\n",fname);
				exit(1);
			}
			forward[i] = (float **)calloc(nband, sizeof(float **));
			for (ib = -1; ++ib < nband; )
			{
				forward[i][ib]= (float *)calloc(nt, sizeof(float *));
				fread(forward[i][ib],sizeof(float),nt,f_in);
			}
			fclose(f_in);
		}
	}
	else
	{
		f_in = fopen("forward.filt","r");
		if ( f_in == NULL )
		{
			printf("Error: cannot open filter file forward.filt\n");
			exit(1);
		}
		forward = (float ***)calloc(1,sizeof(float ***));
		forward[0] = (float **)calloc(nband, sizeof(float **));
		for (ib = -1; ++ib < nband; )
		{
			forward[0][ib]= (float *)calloc(nt, sizeof(float *));
			fread(forward[0][ib],sizeof(float),nt,f_in);
		}
		fclose(f_in);
	}


   /* Check files and count max number of samples */
	f_ntrials = fopen("spike_psth.count", "w");
	if ( f_ntrials == NULL )
	{
		printf("Error: cannot open data file spike_psth.count\n");
		exit(1);
	}

	
	/*allocate space for the vector that keeps track of number of trials and
    its corresponding histogram.  Also outputs number of trials to spike_psth.count*/
    ntrials_value = (int *)calloc(MAXTRIALCOUNT,sizeof(int));
    ntrials_hist = (int *)calloc(MAXTRIALCOUNT,sizeof(int));
    trial_count=0;
	for ( irec = -1; ++irec < nrec; )
	{
		pin = read_data(argvrec[irec]);
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
      nlast = (pin->n/2)*2;
      nt=0;
      while(nt<trial_count)
        {
        if (nlast==ntrials_value[nt])
	  {
            ntrials_hist[nt]++;         
            break;	 
	  }
	nt++;
	if (nt >= MAXTRIALCOUNT ) 
	{
		printf("Too many different number of trials in each file\n"); 
		exit(1);
	}
         }   
      if (nt==trial_count)
      /*this means that a new number of trials occurred*/
          {
        ntrials_value[nt]=nlast;
        trial_count++;
            }
		fprintf(f_ntrials,"%d\n", nlast);
	}

		
	/*now trial_count is equal to the number of different trials that occurred
	ntrials_value and ntrials_hist are of length(trial_count).  This loop finds
	the largest, most common number of trials.*/
	max_trial_count=0;
	current_max_value=0;
	for (nt=-1; ++nt<trial_count;)
	     if (ntrials_hist[nt]>=max_trial_count)
	         if (ntrials_value[nt]>current_max_value)
	            {
	             max_trial_count=ntrials_hist[nt];
	             current_max_value=ntrials_value[nt];
	             }
	
	free(ntrials_hist);
	free(ntrials_value);
    
	/*the last number in spike_psth.count will be the numbers of trials used */
	fprintf(f_ntrials, "%d\n", current_max_value);
	fclose(f_ntrials);

	for ( irec = -1; ++irec < nrec; )
	{
		pin = read_data(argvrec[irec]);
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
        nlast = (pin->n/2)*2;
		if ( matchflg )   /* read spectrogram from file */
		{
			sprintf(fname,"stim_spect%d.dat",irec);
			f_spect = fopen(fname, "r");
			if ( f_spect == NULL )
			{
				printf("Error: cannot open spectrogram data file %s\n",fname);
				exit(1);
			}
			/* Allocated for nband+1 to match the allocation in fband */
			stim_env = (double **)calloc(nband+1, sizeof(double *));
			nlen = argvlen[irec];
			for ( ib = -1; ++ib < nband; )
			{
				stim_env[ib] = (double *)calloc(nlen, sizeof(double));
				fread(stim_env[ib], sizeof(double), nlen, f_spect);
			}
			fclose(f_spect);
		}
		else
      {
			strcpy(stim_spec, pin->stim_spec);
			stim = pin->pstim = setup_stimulus(stim_spec);
			if ( stim == NULL )
			{
				printf("Error setting up stimulus for %s:%s\n", argv[1], pin->stim_spec);
				exit(1);
			}
			/* Pre-process the stimulus */
			stim_fband(stim,&nb,&nlen,sil_window,&stim_env,&f_low,&f_high,&f_width,&f_step);
         nlen += 2*sil_window;

	 printf("nlen is %d\n", nlen);
			if ( lin_flg == 0)
				for ( ib = -1; ++ib < nb; )
					for ( it = -1; ++it < nlen; )
						stim_env[ib][it] = log(stim_env[ib][it]+1.0);

			if ( nb != nband )
			{
				printf("Missmatch between filter and stimulus #bands:%d vs %d\n",
						  nb, nband);
				exit(1);
			}
      }
      
      
     
		/* Open output files if number of trials is correct*/
	if (nlast>=current_max_value)
	{
		sprintf(fname,"spike_psth1_%d.dat",irec);
		f_out_psth1 = fopen(fname,"w");
		sprintf(fname,"spike_psth2_%d.dat",irec);
		f_out_psth2 = fopen(fname,"w");
		sprintf(fname,"spike_pre_%d.dat",irec);
		f_out_pre = fopen(fname,"w");
		if (f_out_psth1 == NULL || f_out_psth2 == NULL || f_out_pre == NULL )
		{
			printf("Error: Could not open output files\n");
			exit(1);
		}
        
		est_spike = (double *)calloc(nlen+end_window, sizeof(double));
		avg_spike1 = (double *)calloc(nlen+end_window, sizeof(double));
		avg_spike2 = (double *)calloc(nlen+end_window, sizeof(double));

      /* For debbuging purposes calculate average for that file */
		for ( ib = -1; ++ib < nband; )
		{
			stim_avg_test[ib] = 0.0;
			for ( it = -1; ++it < nlen; )
				stim_avg_test[ib] += stim_env[ib][it];
			stim_avg_test[ib] /= (double)nlen;
			printf("file %d Band %d average filter %g average file %g\n", irec, ib, stim_avg[ib], stim_avg_test[ib]);
		}

		/* Convolve stim and filter */
		for (ib = -1; ++ib < nband; )
		{
			for ( it = -1; ++it < nlen+end_window; )
			{
			   for (st=-1; ++st < nlen; )
				{
				   dt = it-st+TWINDOW;
					if ( dt < 0 ) break;
					if ( dt > 2*TWINDOW ) continue;
				       
					stim_val = stim_env[ib][st]-stim_avg[ib];
					
					if ( matchflg )
						est_spike[it] += stim_val*forward[irec][ib][dt];
					else
						est_spike[it] += stim_val*forward[0][ib][dt];
             }
			}
		}

		/* Calculate average rates */
		for (it=-1; ++it < nlen+end_window; )
		{
			avg_est += est_spike[it];
			avg_est2 += est_spike[it]*est_spike[it];
		}
		totlen += nlen+end_window;

      /* Ignore last trials if the trial number is greater than current_max_value*/
		for ( i = -1; ++i < current_max_value; )
		{
			ptin = pin->ptrial + i;

			for ( j = -1; ++j < ptin->nspikes; )
			{
				time_spike = ptin->time[j]*1.0e-3;
				bin0 = nint(time_spike) - pin->pre + sil_window;
				if ( bin0 < 0 ) continue;
				if (bin0 >= nlen+end_window ) break;
				if (i%2)
					avg_spike1[bin0] += 1.0;
				else
					avg_spike2[bin0] += 1.0;
			}
		}

		for ( it= -1; ++it < nlen + end_window; )
		{
			avg_spike1[it] /= (current_max_value/2);	
			avg_spike2[it] /= (current_max_value/2);	
		}

		fwrite(avg_spike1, sizeof(double), nlen+end_window, f_out_psth1);
		fwrite(avg_spike2, sizeof(double), nlen+end_window, f_out_psth2);
		fwrite(est_spike, sizeof(double), nlen+end_window, f_out_pre);
		
		fclose(f_out_psth1);
		fclose(f_out_psth2);
		fclose(f_out_pre);
		
		free_dcp_memory(&pin);
		free_env(nlen,nband,stim_env);
		free(est_spike);
		free(avg_spike1);
		free(avg_spike2);
		
	}
	}
	


	/* Generates a fake dcp file assuming a given rate and variance */
	if ( gen_flg )
	{
	  /* dcp_dummy_gen();*/
	  avg_est /= totlen;
	  avg_est2 = 1/((double)(totlen -1))*(avg_est2 - ((double)(totlen))*avg_est*avg_est);

	  for (irec=-1; ++irec < nrec; )
	  {
		pin = read_data(argvrec[irec]);
		if ( matchflg )   /* read spectrogram from file */
		{
			sprintf(fname,"stim_spect%d.dat",irec);
			f_spect = fopen(fname, "r");
			if ( f_spect == NULL )
			{
				printf("Error: cannot open spectrogram data file %s\n",fname);
				exit(1);
			}
			/* Allocated for nband+1 to match the allocation in fband */
			stim_env = (double **)calloc(nband+1, sizeof(double *));
			nlen = argvlen[irec];
			for ( ib = -1; ++ib < nband; )
			{
				stim_env[ib] = (double *)calloc(nlen, sizeof(double));
				fread(stim_env[ib], sizeof(double), nlen, f_spect);
			}
			fclose(f_spect);
		}
		else
      {
			strcpy(stim_spec, pin->stim_spec);
			stim = pin->pstim = setup_stimulus(stim_spec);
			if ( stim == NULL )
			{
				printf("Error setting up stimulus for %s:%s\n", argv[1], pin->stim_spec);
				exit(1);
			}
			/* Pre-process the stimulus */
			stim_fband(stim,&nb,&nlen,sil_window,&stim_env,&f_low,&f_high,&f_width,&f_step);
         nlen += 2*sil_window;

	 printf("nlen is %d\n", nlen);
			if ( lin_flg == 0)
				for ( ib = -1; ++ib < nb; )
					for ( it = -1; ++it < nlen; )
						stim_env[ib][it] = log(stim_env[ib][it]+1.0);

			if ( nb != nband )
			{
				printf("Missmatch between filter and stimulus #bands:%d vs %d\n",
						  nb, nband);
				exit(1);
			}
      }

		est_spike = (double *)calloc(nlen+end_window, sizeof(double));
		/* Convolve stim and filter */
		for (ib = -1; ++ib < nband; )
		{
			for ( it = -1; ++it < nlen+end_window; )
			{
			   for (st=-1; ++st < nlen; )
				{
				   dt = it-st+TWINDOW;
					if ( dt < 0 ) break;
					if ( dt > 2*TWINDOW ) continue;
				       
					stim_val = stim_env[ib][st]-stim_avg[ib];
					
					if ( matchflg )
						est_spike[it] += stim_val*forward[irec][ib][dt];
					else
						est_spike[it] += stim_val*forward[0][ib][dt];
             }
			}
		}

	  dcp_gen(est_spike, nlen+end_window, argvrec[irec], pin, irec, sil_window, avg_est, avg_est2);
		free_dcp_memory(&pin);
		free_env(nlen,nband,stim_env);
		free(est_spike);
	}
	}

		
	free(stim_avg);

}
int dcp_dummy_gen()
{
	return 1;
}

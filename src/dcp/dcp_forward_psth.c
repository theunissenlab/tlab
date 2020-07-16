#include <stdio.h>
#include <malloc.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <sys/types.h>
#include <time.h>



#include "dcp_ft.h"

main(argc, argv)
int argc;
char *argv[];
{
int it, ib, st, dt, nt, nb, cumtimepre;
double bg_mean, stim_onset_mean, stim_mean, onset_mean, bg_var, stim_onset_var, stim_var, onset_var;
int N_trials_stim;
double  bg_mean_temp, stim_onset_mean_temp, stim_mean_temp, onset_mean_temp;
double cov_stim_bg, cov_stim_onset_bg;
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
int matchflg;
char **argv_start;
int argc_start, nmax;
int *ntrials_hist;
int *ntrials_value;
FILE *f_out_psth1, *f_out_psth2, *f_out_pre, *f_count;
FILE *f_ntrials, *f_trial_number, *f_check, *f_spect, *f_means;
FILE *f_avg, *f_bg, *f_stim_means, *f_stim_onset_means, *f_onset_means;
int lin_flg, gen_flg, sil_flg, end_flg;
int sil_window, end_window;
double min_est, max_est, mean_est, rate_factor;
double spike_thres;
int *fake_spike;
int ns;
float f_low, f_high, f_step, f_width;
int TWINDOW;

#define MAX_JN 10
#define ONSET_TIME 100

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



   /* Check files and count max number of samples stim_list.count is the number of the data files that gets calculated (for use in update_filt.m)*/
        f_trial_number=fopen("stim_list.count", "w");
	f_ntrials = fopen("spike_psth.count", "w");
	if ( f_ntrials == NULL )
	{
		printf("Error: cannot open data file spike_psth.count\n");
		exit(1);
	}
	
	/*allocate space for the vector that keeps track of number of trials and
    its corresponding histogram.  Also outputs number of trials to spike_psth.count*/
    ntrials_value = (int *)calloc(20,sizeof(int));
    ntrials_hist = (int *)calloc(20,sizeof(int));
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
         }   
      if (nt==trial_count)
      /*this means that a new number of trials occurred*/
          {
        ntrials_value[nt]=nlast;
        trial_count++;
            }
		fprintf(f_ntrials,"%d\n", nlast);
		fprintf(f_trial_number, "%d\n", irec);
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
	fclose(f_trial_number);
	printf("current_max_value is %d\n", current_max_value);
    bg_mean = 0;
    stim_onset_mean = 0;
    stim_mean= 0;
    onset_mean = 0;
    bg_var = 0;
    stim_onset_var = 0;
    stim_var= 0;
    onset_var = 0;
    N_trials_stim=0;
    cov_stim_bg=0;
    cov_stim_onset_bg=0;	      
	f_bg=fopen("background_means", "w");	
	f_stim_onset_means= fopen("stim_onset_means", "w");
	f_stim_means=fopen("stim_means", "w");

	for ( irec = -1; ++irec < nrec; )
	  {

		pin = read_data(argvrec[irec]);
		cumtimepre=0;
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
     
		/* Open output files if number of trials is correct*/
	  
	  if (matchflg)
	    {
	    nlen = argvlen[irec];	
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

	    }
	        nlast = (pin->n/2)*2;
		/*only do this if there are at least current_max_value number of trials*/
    	if (nlast>=current_max_value)      
	{
		sprintf(fname,"spike_psth1_%d.dat",irec);
		f_out_psth1 = fopen(fname,"w");
		sprintf(fname,"spike_psth2_%d.dat",irec);
		f_out_psth2 = fopen(fname,"w");

		if (f_out_psth1 == NULL || f_out_psth2 == NULL )
		{
			printf("Error: Could not open output files\n");
			exit(1);
		}
        

		avg_spike1 = (double *)calloc(nlen+end_window, sizeof(double));
		avg_spike2 = (double *)calloc(nlen+end_window, sizeof(double));

	}	
		totlen += nlen+end_window;
		N_trials_stim+=pin->n;
	
		for ( i = -1; ++i < pin->n; )
		  {
		     bg_mean_temp = 0;
		    stim_onset_mean_temp = 0;
		    stim_mean_temp= 0;
		    onset_mean_temp = 0;
			ptin = pin->ptrial + i;		
			for ( j = -1; ++j < ptin->nspikes; )
			{
				time_spike = ptin->time[j]*1.0e-3;		     
				bin0 = nint(time_spike) - pin->pre + sil_window;                          
				if ( bin0 < 0 ) 
				  {
			
				  bg_mean_temp +=1.0;
				  continue;
				  }
				if (bin0 >= nlen+end_window ) 
				  {
				    break;
				  }
				if (bin0 < ONSET_TIME)
				  {
				    onset_mean_temp +=1.0;
				    stim_onset_mean_temp +=1.0;		       
				  }
				else
				  {
				    stim_mean_temp +=1.0;
				    stim_onset_mean_temp +=1.0;
				  }
			      /* Ignore last trials if the trial number is greater than current_max_value*/
				if(nlast>=current_max_value && i<current_max_value)
				  {

				if (i%2)
					avg_spike1[bin0] += 1.0;
				else
					avg_spike2[bin0] += 1.0;
				  }
		       
			}
			onset_mean_temp/=(double)ONSET_TIME;
			stim_mean_temp /=(double)(nlen+end_window-ONSET_TIME);
			stim_onset_mean_temp/=(double)nlen+end_window;
			bg_mean_temp/= (double)(pin->pre + sil_window);
			onset_mean += onset_mean_temp;
			stim_mean += stim_mean_temp;
			stim_onset_mean +=stim_onset_mean_temp;
			bg_mean+=bg_mean_temp;
			onset_var += pow(onset_mean_temp,2);
			stim_var += pow(stim_mean_temp,2);
			stim_onset_var +=pow(stim_onset_mean_temp,2);
			bg_var+= pow(bg_mean_temp,2);
			cov_stim_onset_bg+=bg_mean_temp*stim_onset_mean_temp;
			cov_stim_bg+=bg_mean_temp*stim_mean_temp;
			fprintf(f_bg, "%f\n", bg_mean_temp);
			fprintf(f_stim_onset_means,"%f\n", stim_onset_mean_temp);
			fprintf(f_stim_means,"%f\n", stim_mean_temp);
		  }
	
		if (nlast >=current_max_value)
		  {
		for ( it= -1; ++it < nlen + end_window; )
		{
			avg_spike1[it] /= (current_max_value/2);	
			avg_spike2[it] /= (current_max_value/2);	
		}

		fwrite(avg_spike1, sizeof(double), nlen+end_window, f_out_psth1);
		fwrite(avg_spike2, sizeof(double), nlen+end_window, f_out_psth2);
		fclose(f_out_psth1);
		fclose(f_out_psth2);
		free(avg_spike1);
		free(avg_spike2);
		  }
		
		free_dcp_memory(&pin);
		
	
	cumtimepre=+nlen;
	  }
	fclose(f_bg);
	fclose(f_stim_means);
	fclose(f_stim_onset_means);

			onset_mean /= (double)N_trials_stim;
			stim_mean /= (double)N_trials_stim;
			stim_onset_mean /= (double)N_trials_stim;
			bg_mean/= (double)N_trials_stim;

		
			onset_var=1/((double)N_trials_stim-1.0)*(onset_var-(double)N_trials_stim*pow(onset_mean, 2));
			stim_var=1/((double)N_trials_stim-1.0)*(stim_var-(double)N_trials_stim*pow(stim_mean, 2));
			stim_onset_var=1/((double)N_trials_stim-1.0)*(stim_onset_var-(double)N_trials_stim*pow(stim_onset_mean, 2));
			bg_var=1/((double)N_trials_stim-1.0)*(bg_var-(double)N_trials_stim*pow(bg_mean, 2));
			onset_var=sqrt(onset_var);
			stim_var=sqrt(stim_var);
			stim_onset_var=sqrt(stim_onset_var);
			bg_var=sqrt(bg_var);
			cov_stim_onset_bg=1/((double)N_trials_stim-1.0)*(cov_stim_onset_bg-(double)N_trials_stim*bg_mean*stim_onset_mean);
			cov_stim_bg=1/((double)N_trials_stim-1.0)*(cov_stim_bg-(double)N_trials_stim*bg_mean*stim_mean);
			printf("cov_stim_bg is %f\n", cov_stim_bg);
	f_means=fopen("mean_sd_values.dat", "w");
	  fprintf(f_means,"%f %f %f %f %f %f\n", bg_mean*1000, onset_mean*1000, stim_mean*1000, stim_onset_mean*1000, 0, 0);
	fprintf(f_means,"%f %f %f %f %f %f\n", bg_var*1000, onset_var*1000, stim_var*1000, stim_onset_var*1000, cov_stim_bg*1000*1000, cov_stim_onset_bg*1000*1000);
	
	free(stim_avg);

}


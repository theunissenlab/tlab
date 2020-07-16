#include <stdio.h>
#include <malloc.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

/* Calculates 2 order moments between spike trains from a series of pair of
	data files. Ouput:
*/

#include "dcp_ft.h"

main(argc, argv)
int argc;
char *argv[];
{
int i, j, k, id, nt, si;
int id_in, out_n;
int argc_sav;
char **argv_sav;
char stim_spec[MAX_SPEC_LEN];
char col_spec[MAX_SPEC_LEN];
DCP_DATA *pin1, *pin2;
DCP_TRIAL *ptin1, *ptin2, *ptin1_1, *ptin2_1;
DCP_STIM *stim1, *stim2;
char ret_char;
DCP_DATA *read_data();
#define BIN_SIZE 10.0
#define TMAX 500.0
int nbin, ibin0, ibin;
int last_ibin, last_ibin_shift;
double rand_spike;

double *cross_back, *cond_cross_back, *cond_cross_stim, *cross_stim;
double *cross_back_rand, *cross_stim_rand;
double *cond_cross_back_rand, *cond_cross_stim_rand;

double *auto_back, *auto_stim;
double *auto_back_rand, *auto_stim_rand;

int back_count, stim_count;
double back_rate, stim_rate;

double *auto2_back, *auto2_stim;
double *auto2_back_rand, *auto2_stim_rand;
 
int back2_count, stim2_count;
double back2_rate, stim2_rate;

double *xctime;
double **JNcross_back;
double **JNcross_stim;
double **JNcross_back_rand;
double **JNcross_stim_rand;
double **JNauto_back, **JNauto_stim;
double **JNauto2_back, **JNauto2_stim;
double **JNauto_back_rand, **JNauto_stim_rand;
double **JNauto2_back_rand, **JNauto2_stim_rand;

int *JNstim_count, *JNstim2_count;
int *JNback_count, *JNback2_count;

float time_spike1, time_spike2, dt;
float time_shift, time_spike2_shift;
int stim_dur, ntrials;
double stim_len, back_len;
double *JNstim_len, *JNback_len;
double *stim_norm, *back_norm;
double **JNstim_norm, **JNback_norm;
double tres,len;
char *option;
FILE *fout;

   argv ++;
   argc --;

   /* Read arguments */
   tres = BIN_SIZE;
   while ( argc && argv[0][0] == '-' )
   {
      option = argv[0] + 1;
      if ( strcmp(option,"tres") == 0 )
      {
         argv ++;
         argc --;
         if ( argc != 0 ) tres = atof(argv[0]);
      }
      else
      {
			printf("Error: Unknown option %s\n", option);
			printf("dcp_corr {-tres <val>} <file_A.1> <file_B.1> <file_A.2> <file_B.2> ...\n");
			exit(1);
		}
		argv ++;
		argc --;
	}
      

	if ( argc%2 != 0 || argc == 0 )
	{
		printf("Error: dcp_corr takes at an even number arguments\n");
		printf("dcp_corr {-tres <val>} <file_A.1> <file_B.1> <file_A.2> <file_B.2> ...\n");
		exit(1);
	}

   argc_sav = argc;
	argv_sav = argv;
	ntrials = 0;

	/* Check files */
   while ( argc > 1 )
	{
		pin1 = read_data(argv[0]);
		if ( pin1 == NULL )
		{
			printf("Error reading dcp file %s\n", argv[0]);
			exit(1);
		}
		if ( pin1->ssflag )
		{
			printf("Not yet dealing with spike sorted data\n");
			exit(1);
		}
		ntrials += pin1->n;
		free_dcp_memory(&pin1);
		argc -= 2;
		argv++;
		argv++;
	}

   /* Allocate space for all arrays */

   nbin = 2*(int)(TMAX/tres+1)+1;
	ibin0 = (nbin-1)/2 + 1;
	back_count = 0;
	stim_count = 0;
	stim2_count = 0;
	back2_count = 0;
	srand48(time(NULL));
   xctime = (double *)calloc(nbin, sizeof(double));
   cross_back = (double *)calloc(nbin, sizeof(double));
   cond_cross_back = (double *)calloc(nbin, sizeof(double));
   auto_back = (double *)calloc(nbin, sizeof(double));
   auto2_back = (double *)calloc(nbin, sizeof(double));
	cross_stim = (double *)calloc(nbin, sizeof(double));
	cond_cross_stim = (double *)calloc(nbin, sizeof(double));
	auto_stim = (double *)calloc(nbin, sizeof(double));
	auto2_stim = (double *)calloc(nbin, sizeof(double));
   back_norm = (double *)calloc(nbin, sizeof(double));
	stim_norm = (double *)calloc(nbin, sizeof(double));
	cross_back_rand  = (double *)calloc(nbin, sizeof(double));
	cond_cross_back_rand  = (double *)calloc(nbin, sizeof(double));
	auto_back_rand  = (double *)calloc(nbin, sizeof(double));
	auto2_back_rand  = (double *)calloc(nbin, sizeof(double));
	cross_stim_rand = (double *)calloc(nbin, sizeof(double));
	cond_cross_stim_rand = (double *)calloc(nbin, sizeof(double));
	auto_stim_rand = (double *)calloc(nbin, sizeof(double));
	auto2_stim_rand = (double *)calloc(nbin, sizeof(double));
	JNback_norm = (double **)calloc(ntrials, sizeof(double *));
	JNstim_norm = (double **)calloc(ntrials, sizeof(double *));
	JNcross_back = (double **)calloc(ntrials, sizeof(double *));
	JNcross_stim = (double **)calloc(ntrials, sizeof(double *));
	JNauto_back = (double **)calloc(ntrials, sizeof(double *));
	JNauto_stim = (double **)calloc(ntrials, sizeof(double *));
	JNauto2_back = (double **)calloc(ntrials, sizeof(double *));
	JNauto2_stim = (double **)calloc(ntrials, sizeof(double *));
	JNauto_back_rand = (double **)calloc(ntrials, sizeof(double *));
	JNauto_stim_rand = (double **)calloc(ntrials, sizeof(double *));
	JNauto2_back_rand = (double **)calloc(ntrials, sizeof(double *));
	JNauto2_stim_rand = (double **)calloc(ntrials, sizeof(double *));
	JNcross_back_rand  = (double **)calloc(ntrials, sizeof(double *));
	JNcross_stim_rand = (double **)calloc(ntrials, sizeof(double *));
	JNstim_count  = (int *)calloc(ntrials, sizeof(int));
	JNback_count = (int *)calloc(ntrials, sizeof(int));
	JNstim2_count  = (int *)calloc(ntrials, sizeof(int));
	JNback2_count = (int *)calloc(ntrials, sizeof(int));
	JNstim_len  = (double *)calloc(ntrials, sizeof(double));
	JNback_len = (double *)calloc(ntrials, sizeof(double));
	for ( nt = -1; ++nt < ntrials; )
	{
		JNback_norm[nt] = (double *)calloc(nbin, sizeof(double));
		JNstim_norm[nt] = (double *)calloc(nbin, sizeof(double));
		JNcross_back[nt] = (double *)calloc(nbin, sizeof(double));
		JNcross_stim[nt] = (double *)calloc(nbin, sizeof(double));
		JNcross_back_rand[nt]  = (double *)calloc(nbin, sizeof(double));
		JNcross_stim_rand[nt] = (double *)calloc(nbin, sizeof(double));
		JNauto_stim[nt] = (double *)calloc(nbin, sizeof(double));
		JNauto_back[nt] = (double *)calloc(nbin, sizeof(double));
		JNauto2_stim[nt] = (double *)calloc(nbin, sizeof(double));
		JNauto2_back[nt] = (double *)calloc(nbin, sizeof(double));
		JNauto_stim_rand[nt] = (double *)calloc(nbin, sizeof(double));
		JNauto_back_rand[nt] = (double *)calloc(nbin, sizeof(double));
		JNauto2_stim_rand[nt] = (double *)calloc(nbin, sizeof(double));
		JNauto2_back_rand[nt] = (double *)calloc(nbin, sizeof(double));
	}

   argc = argc_sav;
	argv = argv_sav;
	nt = 0;
	stim_len = 0.0;
	back_len = 0.0;
   while ( argc > 1 )
	{
		pin1 = read_data(argv[0]);
		if ( pin1 == NULL )
		{
			printf("Error reading dcp file %s\n", argv[0]);
			exit(1);
		}
		if ( pin1->ssflag )
		{
			printf("Not yet dealing with spike sorted data\n");
			exit(1);
		}
		strcpy(stim_spec, pin1->stim_spec);
		stim1=pin1->pstim = setup_stimulus(stim_spec);
		if ( stim1 == NULL )
		{
			printf("Error setting up stimulus for %s:%s\n", argv[0], pin1->stim_spec);
			exit(1);
		}

		pin2 = read_data(argv[1]);
		if ( pin2 == NULL )
		{
			printf("Error reading dcp file %s\n", argv[1]);
			exit(1);
		}
		if ( pin2->ssflag )
		{
			printf("Not yet dealing with spike sorted data\n");
			exit(1);
		}
		strcpy(stim_spec, pin2->stim_spec);
		stim2=pin2->pstim = setup_stimulus(stim_spec);
		if ( stim2 == NULL )
		{
			printf("Error setting up stimulus for %s:%s\n", argv[1], pin2->stim_spec);
			exit(1);
		}
		if ( pin1->n != pin2->n )
		{
			printf("The two files do not have the same number of trials\n");
			exit(1);
		}
		if ( strcmp(pin1->stim_spec,pin2->stim_spec) )
		{
			printf("The two files have different stimuli:\n");
			printf("\tFile %s: %s\n", argv[0], pin1->stim_spec);
			printf("\tFile %s: %s\n", argv[1], pin2->stim_spec);
			exit(1);
		}

		stim_dur = stim1->dur;
		rand_spike = 1.0/(double)(pin1->n - 1);

		/* Loop through trials */
		for ( i = -1; ++i < pin1->n; )
		{
			stim_len += stim_dur;
			back_len += pin1->pre;
			JNstim_len[nt] += stim_dur;
			JNback_len[nt] += pin1->pre;
			ptin1 = pin1->ptrial + i;
			ptin2 = pin2->ptrial + i;
			if ( i+1 < pin1->n )
			{
			  ptin1_1 = pin1->ptrial + (i+1);
			  ptin2_1 = pin2->ptrial + (i+1);
			}
			else
			{
			  ptin1_1 = pin1->ptrial;
			  ptin2_1 = pin2->ptrial;
			}

			/* Count how many hits for cell1 */
			for ( j = -1; ++j < ptin1->nspikes; )
			{
				time_spike1 = ptin1->time[j]*1.0e-3;
				if ( time_spike1 < pin1->pre )
				{
					back_count++; 
					JNback_count[nt]++;
				}
				else if ( time_spike1 < pin1->pre+stim_dur )
				{
					stim_count++; 
					JNstim_count[nt]++;
				}
         }
         /* ... and cell 2 */
		   for ( j = -1; ++j < ptin2->nspikes; )
		   {
			   time_spike2 = ptin2->time[j]*1.0e-3;
				if ( time_spike2 < pin1->pre )
				{
					back2_count++; 
					JNback2_count[nt]++;
				}
				else if ( time_spike2 < pin1->pre+stim_dur )
				{
					stim2_count++; 
					JNstim2_count[nt]++; 
				}
			}

			if ( pin1->pre < stim_dur ) time_shift = drand48()*pin1->pre;
			else time_shift = drand48()*stim_dur;

			for ( j = -1; ++j < ptin1->nspikes; )
			{
				time_spike1 = ptin1->time[j]*1.0e-3;
				last_ibin = -1;
				last_ibin_shift = -1;
				/* Calculate the cross correlation */
				for ( k = -1; ++k < ptin2->nspikes; )
				{
					time_spike2 = ptin2->time[k]*1.0e-3;
					time_spike2_shift = time_spike2 + time_shift;
					if ( time_spike2 <= pin1->pre+stim_dur && time_spike2 > pin1->pre)
					{
						if (time_spike2_shift > pin1->pre+stim_dur ) time_spike2_shift -= stim_dur;
					}
					else if ( time_spike2 <= pin1->pre)
					{
						if (time_spike2_shift > pin1->pre) time_spike2_shift -= pin1->pre;
					}
					dt = time_spike2 - time_spike1;
					if ( fabs(dt) < TMAX )
					{
						ibin = ibin0 + (int)(nint(dt/tres));
						if ( time_spike1 > pin1->pre && time_spike1 <= pin1->pre+stim_dur
						&& time_spike2 > pin1->pre && time_spike2 <= pin1->pre+stim_dur)
						{
							cross_stim[ibin]++;
							if (ibin != last_ibin ) cond_cross_stim[ibin] ++;
							JNcross_stim[nt][ibin]++;
						}
						else if ( time_spike1 <= pin1->pre && time_spike2 <= pin2->pre )
						{
							cross_back[ibin]++;
							if (ibin != last_ibin ) cond_cross_back[ibin] ++;
							JNcross_back[nt][ibin]++;
						}
						last_ibin = ibin;
				   }
					dt = time_spike2_shift - time_spike1;
					if ( fabs(dt) < TMAX )
					{
						ibin = ibin0 + (int)(nint(dt/tres));
						if ( time_spike1 > pin1->pre && time_spike1 <= pin1->pre+stim_dur
						&& time_spike2 > pin1->pre && time_spike2 <= pin1->pre+stim_dur)
						{
							if (ibin != last_ibin_shift ) cond_cross_stim_rand[ibin] ++;
						}
						else if ( time_spike1 <= pin1->pre && time_spike2 <= pin2->pre )
						{
							if (ibin != last_ibin_shift ) cond_cross_back_rand[ibin] ++;
						}
						last_ibin_shift = ibin;
				   }
					else
					{
						last_ibin_shift = -1;
					}
				}

				/* And the auto_correlation of 1 */
				for ( k = -1; ++k < ptin1->nspikes; )
				{
					time_spike2 = ptin1->time[k]*1.0e-3;
					dt = time_spike2 - time_spike1;
					if ( fabs(dt) < TMAX )
					{
						ibin = ibin0 + (int)(nint(dt/tres));
						if ( time_spike1 > pin1->pre && time_spike1 <= pin1->pre+stim_dur
						&& time_spike2 > pin1->pre && time_spike2 <= pin1->pre+stim_dur)
						{
							auto_stim[ibin]++;
							JNauto_stim[nt][ibin]++; 
						}
						else if ( time_spike1 <= pin1->pre && time_spike2 <= pin2->pre )
						{
							auto_back[ibin]++;
							JNauto_back[nt][ibin]++;
						}
					}
				}
				/* and the shuffled quantities */
				for ( si = -1; ++si < pin1->n; )
				{
					if ( si == i ) continue;
				   ptin1_1 = pin1->ptrial + si;
					
					for ( k = -1; ++k < ptin1_1->nspikes; )
					{
						time_spike2 = ptin1_1->time[k]*1.0e-3;
						dt = time_spike2 - time_spike1;
						if ( fabs(dt) < TMAX )
						{
							ibin = ibin0 + (int)(nint(dt/tres));
							if ( time_spike1 > pin1->pre && time_spike1 <= pin1->pre+stim_dur
							&& time_spike2 > pin1->pre && time_spike2 <= pin1->pre+stim_dur)
							{
								auto_stim_rand[ibin]+=rand_spike;
								JNauto_stim_rand[nt][ibin]+=rand_spike;
							}
							else if ( time_spike1 <= pin1->pre && time_spike2 <= pin2->pre )
							{
								auto_back_rand[ibin]+=rand_spike;
								JNauto_back_rand[nt][ibin]+=rand_spike;
							}
						}
					}
				}

				for ( si = -1; ++si < pin1->n; )
				{
					if ( si == i ) continue;
				   ptin2_1 = pin2->ptrial + si;
					for ( k = -1; ++k < ptin2_1->nspikes; )
					{
						time_spike2 = ptin2_1->time[k]*1.0e-3;
						dt = time_spike2 - time_spike1;
						if ( fabs(dt) < TMAX )
						{
							ibin = ibin0 + (int)(nint(dt/tres));
							if ( time_spike1 > pin1->pre && time_spike1 <= pin1->pre+stim_dur
							&& time_spike2 > pin1->pre && time_spike2 <= pin1->pre+stim_dur)
							{
								cross_stim_rand[ibin]+=rand_spike;
								JNcross_stim_rand[nt][ibin]+=rand_spike;
							}
							else if ( time_spike1 <= pin1->pre && time_spike2 <= pin2->pre )
							{
								cross_back_rand[ibin]+=rand_spike;
								JNcross_back_rand[nt][ibin]+=rand_spike;
							}
						}
					}
				}
			}

			/* The auto correlation of 2 */
			for ( j = -1; ++j < ptin2->nspikes; )
			{
				time_spike1 = ptin2->time[j]*1.0e-3;
				for ( k = -1; ++k < ptin2->nspikes; )
				{
					time_spike2 = ptin2->time[k]*1.0e-3;
					dt = time_spike2 - time_spike1;
					if ( fabs(dt) < TMAX )
					{
						ibin = ibin0 + (int)(nint(dt/tres));
						if ( time_spike1 > pin1->pre && time_spike1 <= pin1->pre+stim_dur
						&& time_spike2 > pin1->pre && time_spike2 <= pin1->pre+stim_dur)
						{
							auto2_stim[ibin]++;
							JNauto2_stim[nt][ibin]++; 
						}
						else if ( time_spike1 <= pin1->pre && time_spike2 <= pin2->pre )
						{
							auto2_back[ibin]++;
							JNauto2_back[nt][ibin]++;
						}
					}
				}
				/* and the shuffled quantities */
				for ( si = -1; ++si < pin1->n; )
				{
					if ( si == i ) continue;
				   ptin2_1 = pin2->ptrial + si;
					for ( k = -1; ++k < ptin2_1->nspikes; )
					{
						time_spike2 = ptin2_1->time[k]*1.0e-3;
						dt = time_spike2 - time_spike1;
						if ( fabs(dt) < TMAX )
						{
							ibin = ibin0 + (int)(nint(dt/tres));
							if ( time_spike1 > pin1->pre && time_spike1 <= pin1->pre+stim_dur
							&& time_spike2 > pin1->pre && time_spike2 <= pin1->pre+stim_dur)
							{
								auto2_stim_rand[ibin]+=rand_spike;
								JNauto2_stim_rand[nt][ibin]+=rand_spike;
							}
							else if ( time_spike1 <= pin1->pre && time_spike2 <= pin2->pre )
							{
								auto2_back_rand[ibin]+=rand_spike;
								JNauto2_back_rand[nt][ibin]+=rand_spike;
							}
						}
					}
				}
			}
			nt ++;
		}
		free_dcp_memory(&pin1);
		free_dcp_memory(&pin2);
		argc -= 2;
		argv++;
		argv++;
	}
	 
   /* First calculate the conditional correlation */
   if ( stim_count < stim2_count )
	{
		for (i=-1;++i<nbin; )
		{
			cond_cross_stim[i] /= (double)stim_count;
			cond_cross_stim_rand[i] /= (double)stim_count;
		}
	}
	else
	{
		for (i=-1;++i<nbin; )
		{
			cond_cross_stim[i] /= (double)stim2_count;
			cond_cross_stim_rand[i] /= (double)stim2_count;
		}
	}

   if ( back_count < back2_count )
	{
		for (i=-1;++i<nbin; )
		{
			cond_cross_back[i] /= (double)back_count;
			cond_cross_back_rand[i] /= (double)back_count;
		}
	}
	else
	{
		for (i=-1;++i<nbin; )
		{
			cond_cross_back[i] /= (double)back2_count;
			cond_cross_back_rand[i] /= (double)back2_count;
		}
	}
	

   /* Calculate the JN Variables, the cross covariances, and the right units */
	/* First the JN variables */
	for ( nt = -1; ++nt < ntrials; )
	{
		JNstim_len[nt] = stim_len - JNstim_len[nt];
		JNback_len[nt] = back_len - JNback_len[nt];

		JNstim_count[nt]  = stim_count - JNstim_count[nt];
		JNback_count[nt] = back_count - JNback_count[nt];
		JNstim2_count[nt]  = stim2_count - JNstim2_count[nt];
		JNback2_count[nt] = back2_count - JNback2_count[nt];
		for ( i = -1; ++i < nbin; )
		{
			JNcross_back[nt][i] = cross_back[i] - JNcross_back[nt][i];
			JNcross_stim[nt][i] = cross_stim[i] - JNcross_stim[nt][i];
			JNauto_back[nt][i] = auto_back[i] - JNauto_back[nt][i];
			JNauto_stim[nt][i] = auto_stim[i] - JNauto_stim[nt][i];
			JNauto2_back[nt][i] = auto2_back[i] - JNauto2_back[nt][i];
			JNauto2_stim[nt][i] = auto2_stim[i] - JNauto2_stim[nt][i];
			JNcross_back_rand[nt][i]  = (cross_back_rand[i] - JNcross_back_rand[nt][i]);
			JNcross_stim_rand[nt][i] = (cross_stim_rand[i] - JNcross_stim_rand[nt][i]);
			JNauto_stim_rand[nt][i] = (auto_stim_rand[i] - JNauto_stim_rand[nt][i]);
			JNauto2_stim_rand[nt][i] = (auto2_stim_rand[i] - JNauto2_stim_rand[nt][i]);
			JNauto_back_rand[nt][i] = (auto_back_rand[i] - JNauto_back_rand[nt][i]);
			JNauto2_back_rand[nt][i] = (auto2_back_rand[i] - JNauto2_back_rand[nt][i]);
		}
	}
	 
   /* Calculate mean rates */
	stim_rate = (double)stim_count/(stim_len/1000);
	stim2_rate = (double)stim2_count/(stim_len/1000);
	back_rate = (double)back_count/(back_len/1000);
	back2_rate = (double)back2_count/(back_len/1000);

   /* Now calculate the normalization values to get the correct units */
	for ( i = -1; ++i < nbin; ) 
	{
	   len = stim_len - abs(i-ibin0)*tres*ntrials;
		if ( len > 0 )
			stim_norm[i] = (1000.0/tres)*(1000.0/len);
		else
			stim_norm[i] = 0.0;
		len = back_len - abs(i-ibin0)*tres*ntrials;
		if ( len > 0 )
			back_norm[i] = (1000.0/tres)*(1000.0/len);
		else
			back_norm[i] = 0.0;
		for ( nt = -1; ++nt < ntrials; )
		{
			len = JNstim_len[nt] - abs(i-ibin0)*tres*ntrials;
			if ( len > 0 )
				JNstim_norm[nt][i] = (1000.0/tres)*(1000.0/len);
			else
				JNstim_norm[nt][i] = 0.0;
			len = JNback_len[nt] - abs(i-ibin0)*tres*ntrials;
			if ( len > 0 )
				JNback_norm[nt][i] = (1000.0/tres)*(1000.0/len);
			else
				JNback_norm[nt][i] = 0.0;
		}

	}

	for ( nt = -1; ++nt < ntrials; )
	{
		for ( i = -1; ++i < nbin; )
		{
			JNcross_back[nt][i] = (JNcross_back[nt][i]-JNcross_back_rand[nt][i])*JNback_norm[nt][i];
			JNauto_back[nt][i] = (JNauto_back[nt][i]-JNauto_back_rand[nt][i])*JNback_norm[nt][i];
			JNcross_back_rand[nt][i] *= JNback_norm[nt][i];
			JNauto_back_rand[nt][i] *= JNback_norm[nt][i];
			JNcross_stim[nt][i] = (JNcross_stim[nt][i]-JNcross_stim_rand[nt][i])*JNstim_norm[nt][i];
			JNauto_stim[nt][i] = (JNauto_stim[nt][i]-JNauto_stim_rand[nt][i])*JNstim_norm[nt][i];
			JNcross_stim_rand[nt][i] *= JNstim_norm[nt][i];
			JNauto_stim_rand[nt][i] *= JNstim_norm[nt][i];
			JNauto2_back[nt][i] = (JNauto2_back[nt][i]-JNauto2_back_rand[nt][i])*JNback_norm[nt][i];
			JNauto2_back_rand[nt][i] *= JNback_norm[nt][i];
			JNauto2_stim[nt][i] = (JNauto2_stim[nt][i]-JNauto2_stim_rand[nt][i])*JNstim_norm[nt][i];
			JNauto2_stim_rand[nt][i] *= JNstim_norm[nt][i];
		}
	}

   /* Calculate cross-covariances */
   for ( i = -1; ++i < nbin; )
	{
		 cross_stim[i] = (cross_stim[i] - cross_stim_rand[i])*stim_norm[i];
		 auto_stim[i] = (auto_stim[i] - auto_stim_rand[i])*stim_norm[i];
		 auto_stim_rand[i] *= stim_norm[i];
		 cross_stim_rand[i] *= stim_norm[i];

		 cross_back[i] = (cross_back[i] - cross_back_rand[i])*back_norm[i];
		 auto_back[i] = (auto_back[i] - auto_back_rand[i])*back_norm[i];
		 auto_back_rand[i] *= back_norm[i];
		 cross_back_rand[i] *= back_norm[i];

		 auto2_stim[i] = (auto2_stim[i] - auto2_stim_rand[i])*stim_norm[i];
		 auto2_stim_rand[i] *= stim_norm[i];
		 auto2_back[i] = (auto2_back[i] - auto2_back_rand[i])*back_norm[i];
		 auto2_back_rand[i] *= back_norm[i];
	}

   fout = fopen("cross.dat","w");
	if ( fout == NULL )
	{
		printf("Error: Could not open output file cross.dat\n");
		exit(1);
	}
	
	fwrite(&nbin,sizeof(int),1, fout);
	fwrite(&ntrials,sizeof(int), 1, fout);
	fwrite(&stim_rate,sizeof(double), 1, fout);
	fwrite(&stim2_rate, sizeof(double), 1, fout);
	fwrite(&back_rate, sizeof(double), 1, fout);
	fwrite(&back2_rate, sizeof(double), 1, fout);

	for ( i = -1; ++i < nbin; ) xctime[i] = (i-ibin0)*tres;

	fwrite(xctime,sizeof(double), nbin, fout);

   /* Write the 4 conditional corellation variables */
	fwrite(cond_cross_stim,sizeof(double), nbin, fout);
	fwrite(cond_cross_back,sizeof(double), nbin, fout);
	fwrite(cond_cross_stim_rand,sizeof(double), nbin, fout);
	fwrite(cond_cross_back_rand,sizeof(double), nbin, fout);


   /* Write the 12 variables. Used in the coherence */
	fwrite(auto_stim,sizeof(double), nbin, fout);
	fwrite(auto_back,sizeof(double), nbin, fout);
	fwrite(auto2_stim,sizeof(double), nbin, fout);
	fwrite(auto2_back,sizeof(double), nbin, fout);
	fwrite(cross_stim,sizeof(double), nbin, fout);
	fwrite(cross_back,sizeof(double), nbin, fout);
	fwrite(auto_stim_rand,sizeof(double), nbin, fout);
	fwrite(auto_back_rand,sizeof(double), nbin, fout);
	fwrite(auto2_stim_rand,sizeof(double), nbin, fout);
	fwrite(auto2_back_rand,sizeof(double), nbin, fout);
	fwrite(cross_stim_rand,sizeof(double), nbin, fout);
	fwrite(cross_back_rand,sizeof(double), nbin, fout);

   /* And write the 12 corresping JN variables */
	for (nt = -1; ++nt < ntrials; )
		fwrite(JNauto_stim[nt],sizeof(double), nbin, fout);
	for (nt = -1; ++nt < ntrials; )
		fwrite(JNauto_back[nt],sizeof(double), nbin, fout);
	for (nt = -1; ++nt < ntrials; )
		fwrite(JNauto2_stim[nt],sizeof(double), nbin, fout);
	for (nt = -1; ++nt < ntrials; )
		fwrite(JNauto2_back[nt],sizeof(double), nbin, fout);
	for (nt = -1; ++nt < ntrials; )
		fwrite(JNcross_stim[nt],sizeof(double), nbin, fout);
	for (nt = -1; ++nt < ntrials; )
		fwrite(JNcross_back[nt],sizeof(double), nbin, fout);
	for (nt = -1; ++nt < ntrials; )
		fwrite(JNauto_stim_rand[nt],sizeof(double), nbin, fout);
	for (nt = -1; ++nt < ntrials; )
		fwrite(JNauto_back_rand[nt],sizeof(double), nbin, fout);
	for (nt = -1; ++nt < ntrials; )
		fwrite(JNauto2_stim_rand[nt],sizeof(double), nbin, fout);
	for (nt = -1; ++nt < ntrials; )
		fwrite(JNauto2_back_rand[nt],sizeof(double), nbin, fout);
	for (nt = -1; ++nt < ntrials; )
		fwrite(JNcross_stim_rand[nt],sizeof(double), nbin, fout);
	for (nt = -1; ++nt < ntrials; )
		fwrite(JNcross_back_rand[nt],sizeof(double), nbin, fout);

	fclose(fout);

}

#ifdef SMOOTHING
#define NSMOOTH 11
double triang[NSMOOTH] = { 0.1667, 0.3333, 0.5000, 0.6667, 0.8333, 1.0000,
									0.8333, 0.6667, 0.5000, 0.3333, 0.1667 };
double *cross_stim_rate, *cross_back_rate;
double *auto_stim_rate, *auto_back_rate;
double totmult_stim, totmult_back;
	cross_stim_rate  = (double *)calloc(nbin, sizeof(double));
	cross_back_rate = (double *)calloc(nbin, sizeof(double));
	auto_stim_rate = (double *)calloc(nbin, sizeof(double));
	auto_back_rate  = (double *)calloc(nbin, sizeof(double));
	auto2_stim_rate = (double *)calloc(nbin, sizeof(double));
	auto2_back_rate  = (double *)calloc(nbin, sizeof(double));
/* I don't know if smoothing is important */
   /* Estimate cross <x1>2, <x2>2 and <x1x2> - by smoothing shift predictors */
   for ( i = -1; ++i < nbin; )
	{
		totmult_stim=0.0;
		totmult_back=0.0;
		totmult2_stim=0.0;
		totmult2_back=0.0;
		for ( j = -1; ++ j < NSMOOTH;)
		{
			id = i - NSMOOTH/2 + j;
         if ( id < 0 ) continue;
			if ( id == nbin ) break;
			if ( stim_count[id] != 0 )
			{
				cross_stim_rate[i] += triang[j]*(double)cross_stim_rand[id]/(double)stim_count[id];
				auto_stim_rate[i] += triang[j]*(double)auto_stim_rand[id]/(double)stim_count[id];
				totmult_stim += triang[j];
			}
			if ( stim2_count[id] != 0 )
			{
				auto2_stim_rate[i] += triang[j]*(double)auto2_stim_rand[id]/(double)stim_count[id];
				totmult2_stim += triang[j];
			}
			if ( back_count[id] != 0 )
			{
				cross_back_rate[i] += triang[j]*(double)cross_back_rand[id]/(double)back_count[id];
				auto_back_rate[i] += triang[j]*(double)auto_back_rand[id]/(double)back_count[id];
				totmult_back += triang[j];
			}
			if ( back2_count[id] != 0 )
			{
				auto2_back_rate[i] += triang[j]*(double)auto2_back_rand[id]/(double)back_count[id];
				totmult2_back += triang[j];
			}
		}
		cross_stim_rate[i] /= totmult_stim;
		auto_stim_rate[i] /= totmult_stim;
		auto2_stim_rate[i] /= totmult2_stim;
		cross_back_rate[i] /= totmult_back;
		auto_back_rate[i] /= totmult_back;
		auto2_back_rate[i] /= totmult2_back;
	}
	
	/* Do the JN statistics */
	JNMeancross_back[i] += JNcross_back[nt][i];
	JNMeancross_stim[i] += JNcross_stim[nt][i];
	JNMeancross_back_rand[i] += JNcross_back_rand[nt][i];
	JNMeancross_stim_rand[i] += JNcross_stim_rand[nt][i];
	for ( i = -1; ++ i < nbin; )
	{
		JNMeancross_back[i] /= ntrials;
		JNMeancross_stim[i] /= ntrials;
		JNMeancross_back_rand[i] /= ntrials;
		JNMeancross_stim_rand[i] /= ntrials;
		for ( nt = -1; ++nt < ntrials; )
		{
			JNStdcross_back[i] += (JNcross_back[nt][i] - JNMeancross_back[i])*
										 (JNcross_back[nt][i] - JNMeancross_back[i]);
			JNStdcross_stim[i] += (JNcross_stim[nt][i] - JNMeancross_stim[i])*
										 (JNcross_stim[nt][i] - JNMeancross_stim[i]);
			JNStdcross_back_rand[i] += (JNcross_back_rand[nt][i] - JNMeancross_back_rand[i])*
										 (JNcross_back_rand[nt][i] - JNMeancross_back_rand[i]);
			JNStdcross_stim_rand[i] += (JNcross_stim_rand[nt][i] - JNMeancross_stim_rand[i])*
										 (JNcross_stim_rand[nt][i] - JNMeancross_stim_rand[i]);
		}
		JNStdcross_back[i] = sqrt(ntrials/(ntrials-1)*JNStdcross_back[i]);
		JNStdcross_stim[i] = sqrt(ntrials/(ntrials-1)*JNStdcross_stim[i]);
		JNStdcross_back_rand[i] = sqrt(ntrials/(ntrials-1)*JNStdcross_back_rand[i]);
		JNStdcross_stim_rand[i] = sqrt(ntrials/(ntrials-1)*JNStdcross_stim_rand[i]);
	}

#endif

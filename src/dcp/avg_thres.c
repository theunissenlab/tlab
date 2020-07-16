#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "nr2.h"
#include "dcp_ft.h"

float AD_to_mV();

/* values for input line 1 */
float cal_m=13374.804688;  /* AD tics / V input */
float cal_b=-9.578730;

main(argc, argv)
int argc;
char *argv[];
{
DCP_DATA *pin;
DCP_TRIAL *ptin;
char ret_char;
DCP_DATA *read_data();
int i, j;
int max_t, min_t, dt;
int stim_dur, lat_points, stim_end;
int npoints, pre_points, stim_points;
int lat = 5;						/* Latency in ms */
int post = 20;                /* Post period in ms */
float *wave_avg;
float guv_vscale=10.0;
float back_avg, back_var;
float stim_avg, stim_var;
float stim_min, stim_max;
double pmax, pmin, pmaxmin;
float dmax, dmin;
char *dur_str;
FILE *fout;
#define NBINS 100
float back_max, back_prob[NBINS];
float xval, pxgauss;
int ind;


	if ( argc != 2 )
	{
		printf ("Error: wrong number of arguments. Usage:\n");
		printf ("wave_thresh  <dcp_filename>\n");
	}

   if ( (pin = read_data(argv[1])) == NULL ) exit(1);

	printf("%d trials found\n", pin->n);

	if ( pin-> waveform != 1 )
	{
	  printf("No waveform data for dcp file %s\n");
	  exit(1);
	}

   npoints = (pin->cd*pin->adfreq)/1000;
	wave_avg = (float *)calloc(npoints, sizeof(float));
	ptin = pin->ptrial;
	if ( ptin->npoints != npoints )
	{
		printf("Internal error\n");
		exit(1);
	}
	for ( j = -1; ++j < npoints; )
	  wave_avg[j] = (float)ptin->waveform[j];

   fout = fopen("wave.avg", "w");
	for ( j = -1; ++j < npoints; )
	{
	  wave_avg[j]=AD_to_mV(wave_avg[j], guv_vscale, (float)pin->adgain); 
     fprintf(fout,"%f %f\n", (j*1000.0)/pin->adfreq, wave_avg[j]);
	}
	fclose(fout);

	/* Calculate mean and variance for times before stim */
   pre_points = (pin->pre*pin->adfreq)/1000;
	back_avg = 0.0;
	back_var = 0.0;
	back_max = 0.0;
	for ( j = -1; ++j < pre_points; )
	{
		back_avg += wave_avg[j];
		back_var += wave_avg[j]*wave_avg[j];
		if ( wave_avg[j] > back_max ) back_max = wave_avg[j];
		if ( -wave_avg[j] > back_max ) back_max = -wave_avg[j];
	}
	for ( j = -1; ++j < NBINS; ) back_prob[j] = 0.0;
	for ( j = -1; ++j < pre_points; )
	{
		ind = (int)((wave_avg[j]/back_max)*NBINS/2 + NBINS/2 - 0.5);
		if ( ind < 0 ) ind = 0;
		if (ind > NBINS ) ind = NBINS;
		back_prob[ind] +=1.0/pre_points;
	}
	back_avg /= (float)pre_points;
	back_var = ((float)pre_points/(pre_points-1))*(back_var/(float)(pre_points) -back_avg*back_avg);
   printf("Background = %f +- %f\n",back_avg, sqrt(back_var));
	fout=fopen("prop.bak","w");
	for ( j = -1; ++j < NBINS; )
	{
		xval = (j-NBINS/2+0.5)*back_max*2.0/NBINS;
		pxgauss = (1.0/sqrt(2.0*M_PI*back_var))*exp(-0.5*(xval-back_avg)*(xval-back_avg)/back_var)*(2.0*back_max/NBINS);
		fprintf(fout,"%g %g %g\n",xval, pxgauss, back_prob[j]);
	}
	fclose(fout);

	dur_str = strtok(pin->stim_spec," ");
	while ( dur_str[0] != 'd' || dur_str[1] != '=' ) 
	{
		dur_str = strtok(NULL," ");
		if ( dur_str == NULL )
		{
			printf("Error parsing stimulus spec\n");
			exit(1);
		}
	}

	sscanf(dur_str,"d=%d",&stim_dur);
   stim_points = ((stim_dur+post)*pin->adfreq)/1000;
	lat_points = (lat*pin->adfreq)/1000;
	stim_avg = 0.0;
	stim_var = 0.0;
	stim_max = wave_avg[pre_points+lat_points];
	stim_min = wave_avg[pre_points+lat_points];
	max_t = min_t = pre_points+lat_points;
	stim_end = pre_points+lat_points+stim_points;
	if ( stim_end > npoints ) stim_end = npoints;
	for ( j = pre_points+lat_points; ++j < stim_end; )
	{
		stim_avg += wave_avg[j];
		stim_var += wave_avg[j]*wave_avg[j];
		if ( wave_avg[j] > stim_max ) 
		{
			stim_max = wave_avg[j];
			max_t = j;
		}
		if ( wave_avg[j] < stim_min )
		{
			stim_min = wave_avg[j];
			min_t = j;
		}
	}
	stim_avg /= (float)stim_points;
	stim_var = ((float)stim_points/(stim_points-1))*(stim_var/(float)(stim_points) -stim_avg*stim_avg);
	dt = min_t - max_t;
	if ( dt < 0 ) dt = -dt;
	dmax = stim_max-back_avg;
	if ( dmax < 0.0 ) dmax = -dmax;
	dmin = stim_min-back_avg;
	if ( dmin < 0.0 ) dmin = -dmin;
	pmax = 0.5*erff(dmax/sqrt(back_var)) + 0.5;
	pmax = 1.- pow(pmax,(double)dt);
	pmin =  0.5*erff(dmin/sqrt(back_var)) + 0.5;
	pmin = 1.- pow(pmin,(double)dt);
	printf("Stim = %f +- %f\n", stim_avg, sqrt(stim_var));
	printf("\tMax = %f T_max = %f P_max = %g\n", stim_max, (max_t-pre_points)*1000.0/(pin->adfreq), pmax) ;
	printf("\tMin = %f T_min = %f P_min = %g\n", stim_min, (min_t-pre_points)*1000.0/(pin->adfreq), pmin) ;
	printf("\tPmax&Pmin=%g\n", pmax*pmin);
	if ( pmax < 0.001 || pmin < 0.001 ) printf("ABOVE THRESHOLD\n");
	else printf("BELOW THRESHOLD\n");
	  
}
/* ------------------------ AD_to_mV() ----------------------------------- */
float AD_to_mV(xAD,vscale,ADgain)
     float xAD;         /* value (A/D units) */
     float vscale;      /* voltage gain (output V) / ( input V) */
     float ADgain;      /* (output V) / (input V) */
{
  float mV,AD,V;

  V = ADgain*(xAD - cal_b)/cal_m;
  mV = (V*1000.0)/vscale;
  return (mV);
}

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
DCP_DATA *psig, *ptest;
DCP_TRIAL *ptsig, *pttest;
DCP_DATA *read_data();
int  j, k;
int npoints;
float lat=5, post=20;                     /* Latency and post stim in ms */
int stim_dur;
int  stim_points, lat_points, stim_end, pre_points;
char *dur_str;
float *wave_sig, *wave_test, *wave_std;
double avg_sig, avg_test;
double std_test, avg_std_test;
double cross_val, cross_err, cross_norm;


	if ( argc != 3 )
	{
		printf ("Error: wrong number of arguments. Usage:\n");
		printf ("noise_corr_avg  <sig_file> <test_file>\n");
	}

   if ( (psig = read_data(argv[1])) == NULL ) exit(1);
   if ( (ptest = read_data(argv[2])) == NULL ) exit(1);

	if ( psig->n != 2 || ptest->n != 2 )
	{
		printf("Unexpected number of trials\n");
		exit(1);
	}

	if ( psig->waveform != 1 || ptest->waveform != 1 )
	{
	  printf("No waveform data\n");
	  exit(1);
	}

   npoints = (psig->cd*psig->adfreq)/1000;
   pre_points = (psig->pre*psig->adfreq)/1000;
	wave_sig = (float *)calloc(npoints, sizeof(float));
	wave_test = (float *)calloc(npoints, sizeof(float));
	wave_std = (float *)calloc(npoints, sizeof(float));

	pttest = ptest->ptrial;
	if ( pttest->npoints != npoints )
	{
		printf("Internal error\n");
		exit(1);
	}
	avg_test = 0.0;
	avg_std_test = 0.0;
	for ( j = -1; ++j < npoints; )
	{
	  wave_test[j] = (float)pttest->waveform[j];
	  wave_std[j] = (float)pttest[1].waveform[j];
	}

	for ( j = -1; ++j < pre_points; )
	{
	  avg_test += wave_test[j];
	  avg_std_test += wave_std[j];
	}
	avg_test /= (double)pre_points;
	avg_std_test /= (double)pre_points;

	std_test = 0.0;
	for ( j = -1; ++j < pre_points; )
	{
	  std_test += (wave_test[j]-avg_test)*(wave_test[j]-avg_test);
	}
	std_test /= (double)(pre_points-1);
	std_test = sqrt(std_test);

	ptsig = psig->ptrial;
	if ( ptsig->npoints != npoints )
	{
		printf("Internal error\n");
		exit(1);
	}
	avg_sig = 0.0;
	for ( j = -1; ++j < npoints; )
	  wave_sig[j] = (float)ptsig->waveform[j];
	for ( j = -1; ++j < pre_points; )
	  avg_sig += wave_sig[j];
	avg_sig /= (double)pre_points;

	/* Calculate cross correlation and expected variance */
	cross_val = 0.0;
	cross_err = 0.0;
	cross_norm = 0.0;

   dur_str = strtok(psig->stim_spec," ");
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
   stim_points = ((stim_dur+post)*psig->adfreq)/1000;
   lat_points = (lat*psig->adfreq)/1000;
   stim_end = pre_points+lat_points+stim_points;
   if ( stim_end > npoints ) stim_end = npoints;
   for ( j = pre_points+lat_points, k=-1; ++k < pre_points+lat_points; j++)
	{
		cross_val += (wave_sig[j]-avg_sig)*(wave_test[k]-avg_test);
      cross_norm += (wave_sig[j]-avg_sig)*(wave_sig[j]-avg_sig);
	}
	cross_val /= cross_norm;
	printf("%g\n", cross_val);

}

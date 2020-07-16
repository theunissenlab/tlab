#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "dcp_ft.h"

main(argc, argv)
int argc;
char *argv[];
{
#define NDIV 10
DCP_DATA *pin;
DCP_TRIAL *ptin;
DCP_DATA *read_data();
int npoints;
int i, j;
float *wave_avg;
float *wave_std;
float pow_trial;
float pow_avg, pow_std;
float pow_min, pow_max;
char name_out[256];
char str_spec[256];
int nbins;
int *pow_hist, pow_bin, pow_bin_th;
int *bin_no, out_n, cum_count;


	if ( argc != 2 )
	{
		printf ("Error: wrong number of arguments. Usage:\n");
		printf ("wave_avg  <dcp_filename>\n");
	}

   if ( (pin = read_data(argv[1])) == NULL ) exit(1);

	printf("%d trials found\n", pin->n);

	if ( pin->waveform != 1 )
	{
	  printf("No waveform data for dcp file %s\n");
	  exit(1);
	}

   npoints = (pin->cd*pin->adfreq)/1000;
	wave_avg = (float *)calloc(npoints, sizeof(float));
	wave_std = (float *)calloc(npoints, sizeof(float));
	pow_avg = 0.0;
	nbins = pin->n/NDIV;
	for ( i = -1; ++i < pin->n; )
	{
		ptin = pin->ptrial+i;
		if ( ptin->npoints != npoints )
		{
			printf("Internal error\n");
			exit(1);
		}
		pow_trial = 0.0;
		for ( j = -1; ++j < npoints; )
		  pow_trial += (float)ptin->waveform[j]*ptin->waveform[j];
		pow_avg += pow_trial;
		if ( i )
		{
			if ( pow_trial < pow_min ) pow_min = pow_trial;
			if ( pow_trial > pow_max ) pow_max = pow_trial;
		}
		else
		{
			pow_min = pow_max = pow_trial;
		}
	}
	pow_avg /= (float)pin->n;

   pow_std = 0.0;
	pow_hist = (int *)calloc(nbins, sizeof(int));
	bin_no = (int *)calloc(pin->n, sizeof(int));

	for ( i = -1; ++i < pin->n; )
	{
		ptin = pin->ptrial+i;
		pow_trial = 0.0;
		for ( j = -1; ++j < npoints; )
		  pow_trial += (float)ptin->waveform[j]*ptin->waveform[j];
		pow_std += (pow_avg-pow_trial)*(pow_avg-pow_trial);
		pow_bin = (int)((pow_trial-pow_min)/(pow_max-pow_min)*nbins);
		if ( pow_bin >= nbins ) pow_bin = nbins-1;
		bin_no[i] = pow_bin;
		pow_hist[pow_bin] ++;
	}
	pow_std /= (float)(pin->n-1);
   pow_std = sqrt(pow_std);

   printf("\n\n");
	for ( i = -1; ++i < nbins; ) 
		printf("%d ",pow_hist[i]);
   printf("\n\n");

   pow_bin_th = nbins;
	cum_count = 0;
	for ( i = -1; ++i < nbins-2; ) 
	{
		cum_count += pow_hist[i];
		if ( cum_count < pin->n/2) continue;
		if (pow_hist[i] == 0 && pow_hist[i+1] == 0 && pow_hist[i+2] == 0 ) { pow_bin_th = i; break;}
	}

	out_n = pin->n;
	for ( i = -1; ++i < pin->n; )
	{
		ptin = pin->ptrial+i;
      if ( bin_no[i] > pow_bin_th ) 
		{
			printf("Trial %d has extranous power spectra\n", i+1);
			out_n --;
		}
		else
		{
			for ( j = -1; ++j < npoints; ) 
			  wave_avg[j] += (float)ptin->waveform[j];
		}
	}
	for ( j = -1; ++j < npoints; ) wave_avg[j] /= (float)out_n;
	for ( i = -1; ++i < pin->n; )
	{
		ptin = pin->ptrial+i;
      if ( bin_no[i] <= pow_bin_th ) 
		{
			for ( j = -1; ++j < npoints; ) 
			  wave_std[j] += ((float)ptin->waveform[j] - wave_avg[j])*((float)ptin->waveform[j] - wave_avg[j]);
		}
	}

	for ( j = -1; ++j < npoints; ) 
	{
		wave_std[j] /= (float)((out_n)*(out_n-1));
		wave_std[j] = sqrt(wave_std[j]);
		pin->ptrial[0].waveform[j] = (short)(wave_avg[j]);
		pin->ptrial[1].waveform[j] = (short)(wave_std[j]);
	}

   sprintf( str_spec," Avg and Std n=%d", out_n);
   strcat( pin->col_spec, str_spec);
	pin->n = 2;
   sprintf(name_out,"%s.avg", argv[1]);
	write_dcp_data( pin, name_out);
}

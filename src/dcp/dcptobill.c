#include <stdio.h>
#include <malloc.h>
#include <string.h>

#include "dcp_ft.h"

main(argc, argv)
int argc;
char *argv[];
{
int i,j;
char file_name[256];
DCP_DATA *pin;
DCP_DATA *read_data();
DCP_TRIAL *ptin;
FILE *f_index;
FILE *f_bill;
float time_base;
float time_spike;

	if ( argc != 2 )
	{
		printf("Error: dcptobill takes one argument\n");
		printf("dcptobill <dcp_file>\n");
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

	f_index = fopen("index.bill","a");
	fprintf(f_index, "%s\t%s\t%d\t%d\n", argv[1], pin->stim_spec, pin->pre, pin->cd);
	fclose(f_index);

   sprintf(file_name,"%s.bill", argv[1]);
	f_bill = fopen(file_name,"w");
	if ( f_bill == NULL )
	{
		printf("Error could not open output file %s\n", file_name);
		exit(1);
	}

   /* Dump spikes */
	time_base = 0.0;
	for ( i = -1; ++i < pin->n; )
	{
		ptin = pin->ptrial + i;
		for ( j = -1; ++ j < ptin->nspikes; )
		{
			time_spike = ptin->time[j]*1.0e-6;
			fprintf(f_bill,"%.6f\n", time_spike+time_base);
		}
		time_base += 10.0;
	}
	fclose(f_bill);
}


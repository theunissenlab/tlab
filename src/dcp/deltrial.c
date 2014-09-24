#include <stdio.h>
#include <malloc.h>

#include "dcp_ft.h"

main(argc, argv)
int argc;
char *argv[];
{
int i, j;
DCP_DATA *pin;
DCP_DATA *pout;
DCP_TRIAL *ptin, *ptout;
char ret_char;
int deltrial;
DCP_DATA *read_data();

	if ( argc != 4 )
	{
		printf("Error: deltrial takes 3 arguments\n");
		printf("deltrial <trial number> <file_in> <file_out> \n");
		exit(1);
	}

   /* Argument check */
	deltrial = atoi(argv[1]);
	if ( file_exists(argv[3]) ) 
	{
		printf("Output file %s exists. Overwrite? (y/n)\n", argv[3]);
		ret_char = (char)getchar();
		if ( ret_char != 'y' ) exit(1);
	}

	pin = read_data(argv[2]);
	if ( pin == NULL ) 
	{
		printf("Error reading data\n");
		exit(1);
	}

   if ( deltrial < 1 || deltrial > pin->n )
	{
		printf("Trial number %d does not exists in input file\n", deltrial);
	   free_dcp_trial_memory(pin);
	   free(pin);
		exit(1);
	}

   /* Copy data except for the given trial. Note that time and waveform
		are not copied */ 
   pout = (DCP_DATA *)calloc(1,sizeof(DCP_DATA));
   memcpy(pout, pin, sizeof(DCP_DATA));
	pout->ptrial = (DCP_TRIAL *)calloc(pin->n-1, sizeof(DCP_TRIAL));
	pout->n = pin->n - 1;
	for ( i = -1,j=0; ++i < pin->n; )
	{
		if ( i == (deltrial-1) ) continue;
		ptin = pin->ptrial + i;
		ptout = pout->ptrial + j;
		memcpy(ptout, ptin, sizeof(DCP_TRIAL));
		j ++;
	}

   if ( write_dcp_data(pout, argv[3]) == RET_ERR )
		printf("Error writing out data\n");
	else
		printf("Wrote new data file %s\n", argv[3]);


}


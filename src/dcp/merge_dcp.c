#include <stdio.h>
#include <malloc.h>
#include <string.h>

#include "dcp_ft.h"

main(argc, argv)
int argc;
char *argv[];
{
int i, j;
int id_in, out_n;
char stim_spec[MAX_SPEC_LEN];
char col_spec[MAX_SPEC_LEN];
DCP_DATA *pin;
DCP_DATA *pout;
DCP_TRIAL *ptin, *ptout;
char ret_char;
DCP_DATA *read_data();

	if ( argc < 3 )
	{
		printf("Error: merge_dcp takes at least 3 arguments\n");
		printf("merge_dcp <file_in1> <file_in2>.. <file_out> \n");
		exit(1);
	}

   /* Argument check */
	if ( file_exists(argv[argc-1]) ) 
	{
		printf("Output file %s exists. Overwrite? (y/n)\n", argv[argc-1]);
		ret_char = (char)getchar();
		if ( ret_char != 'y' ) exit(1);
	}

   /* Do a first loop to get number of trials */
	out_n = 0;
   for ( id_in = 0; ++id_in < argc-1; )
	{
		pin = read_data(argv[id_in]);
		if ( pin == NULL ) 
		{
			printf("Error reading data for file %s\n", argv[id_in]);
			exit(1);
		}
		out_n += pin->n;
		if ( id_in == 1 )
		{
			strcpy(stim_spec, pin->stim_spec);
			strcpy(col_spec,"merged files:");
			strcat(col_spec, argv[id_in]);
		}
		else
		{
			if ( strcmp(stim_spec, pin->stim_spec) )
			{
				printf("Warning: Missmatch of stimulus specification in your files\n");
				printf("First file %s had:\n%s\n", argv[1], stim_spec);
				printf("Current file %s has:\n%s\n", argv[id_in], pin->stim_spec);
				printf("Do you want to proceed?(y/n)\n");
			   if ( ret_char != 'y' ) exit(1);
				printf("Merge files will use first stim as spec.\n");
			}
			if ( strlen(col_spec)+strlen(argv[id_in]) < MAX_SPEC_LEN-2 )
				strcat(col_spec, argv[id_in]);
			else if ( strlen(col_spec) < MAX_SPEC_LEN-5 )
				strcat(col_spec,"...");
			if ( strlen(col_spec) < MAX_SPEC_LEN-2 )
				strcat(col_spec," ");
		}
	   free_dcp_trial_memory(pin);
	   free(pin);
	}

   /* Merge all data into one data file */
   pout = (DCP_DATA *)calloc(1,sizeof(DCP_DATA));
   j = 0;
   for ( id_in = 0; ++id_in < argc-1; )
	{
		pin = read_data(argv[id_in]);
		if ( id_in == 1 )
		{
			memcpy(pout, pin, sizeof(DCP_DATA));
			pout->ptrial = (DCP_TRIAL *)calloc(out_n, sizeof(DCP_TRIAL));
			pout->n = out_n;
			strcpy(pout->col_spec, col_spec);
		}
		for ( i = -1; ++i < pin->n; )
		{
			ptin = pin->ptrial + i;
			ptout = pout->ptrial + j;
			memcpy(ptout, ptin, sizeof(DCP_TRIAL));
			j ++;
		}
		/* pin is not freed here to keep data stuff around */
	}

   if ( write_dcp_data(pout, argv[argc-1]) == RET_ERR )
		printf("Error writing out data\n");
	else
		printf("Wrote new data file %s\n", argv[argc-1]);


}


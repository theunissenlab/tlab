#include <stdio.h>
#include <malloc.h>
#include <string.h>

#include "dcp_ft.h"

main(argc, argv)
int argc;
char *argv[];
{
char new_spec[MAX_SPEC_LEN];
char temp_spec[MAX_SPEC_LEN];
char *new_path;
char file_name[256];
DCP_DATA *pin;
char ret_char;
char *str_tok;
DCP_DATA *read_data();
int froggy_exception;

	if ( argc != 2 )
	{
		printf("Error: dcp_fixpath takes one argument\n");
		printf("dcp_fixpath_adaptive <dcp_file> \n");
		exit(1);
	}

   /* Read file */
	pin = read_data(argv[1]);
	if ( pin == NULL ) 
	{
		printf("Error reading dcp file %s\n", argv[1]);
		exit(1);
	}

	strcpy(temp_spec, pin->stim_spec);
	str_tok = strtok( temp_spec, " ");
	if ( strcmp(str_tok, "stim") != 0 ) 
	{
		printf("Internal Error\n");
		exit(1);
	}
	strcpy(new_spec, str_tok);
	str_tok = strtok( NULL, " ");
	if ( strcmp(str_tok,"song") != 0 && strcmp(str_tok,"seg") != 0 )
	{
		printf("Stimulus spec is not a song or a seg but a %s\n", str_tok);
		exit(1);
	}
	strcat(new_spec, " ");
	strcat(new_spec, str_tok);

	str_tok = strtok( NULL, " ");
	if  (strncmp( str_tok, "p=/data/sarah/a", 15) == 0 )
	{
		printf("Old song path is: %s\n", str_tok+2);
		froggy_exception = 0;
	}
	else if (strncmp( str_tok, "p=/data/sarah/bf_zf/stims/bfsongs", 31) == 0 )
	{
		printf("Old song path is: %s\n", str_tok+2);
		froggy_exception = 1;
	}
	else if (strncmp( str_tok, "p=/data/sarah/bf_zf/stims/zfsongs", 31) == 0 )
	{
		printf("Old song path is: %s\n", str_tok+2);
		froggy_exception = 2;
	}
	else if (strncmp( str_tok, "p=/data/sarah/bf_zf/stims/mlnoise", 31) == 0 )
	{
		printf("Old song path is: %s\n", str_tok+2);
		froggy_exception = 3;
	}
	else
	{
		printf("Error: Old song path is not accounted for");
		printf("Old song path is: %s\n", str_tok+2);
		exit(1);
	}
	if ( froggy_exception == 0 )
	{
		new_path = "/auto/fsong/nandini/adaptive_stim/fullname_but_really_no_syl";
	}
	else if ( froggy_exception == 1 )
	{
		new_path = "/auto/fsong/nandini/bf_zf/bfsongs";
	}
	else if ( froggy_exception == 2 )
	{
		new_path = "/auto/fsong/nandini/bf_zf/zfsongs";
	}
	else if ( froggy_exception == 3 )
	{
		new_path = "/auto/fsong/nandini/bf_zf/mlnoise";
	}

	str_tok = strtok(NULL," ");
	/* Check to see if file exists */
	strcpy(file_name, new_path);
	strcat(file_name, "/");
	strcat(file_name, str_tok);

	if  ( fopen(file_name, "r") == NULL )
	{
		printf("Warning: Could not open song file %s\n", file_name);
		/* printf("Write modified dcp file anyway?\n");
		ret_char = (char)getchar();
		if ( ret_char != 'y' ) exit(1);
		*/
	}
	strcat(new_spec, " ");
	strcat(new_spec, "p=");
	strcat(new_spec, new_path);
	strcat(new_spec, " ");
	strcat(new_spec, str_tok);
	str_tok = strtok(NULL," ");
	while ( str_tok != NULL )
	{
		strcat(new_spec, " ");
		strcat(new_spec, str_tok);
		str_tok = strtok(NULL," ");
	}
	strcpy(pin->stim_spec, new_spec);
	printf("The new stimulus spec is: %s\n", pin->stim_spec);
				
   if ( write_dcp_data(pin, argv[1]) == RET_ERR )
		printf("Error writing out updated dcp file\n");

}


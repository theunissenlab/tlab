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

	if ( argc != 4 )
	{
		printf("Error: dcp_newpath takes three argument\n");
		printf("dcp_newsong <dcp_file> <new_path> <new_song>\n");
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
	if  (strncmp( str_tok, "p=", 2) == 0)
	{
		printf("Old song path is: %s\n", str_tok+2);
	}
	new_path = argv[2];

	str_tok = strtok(NULL," ");  /* THis points to old file */
	/* Check to see if file exists */
	strcpy(file_name, new_path);
	strcat(file_name, "/");
	strcat(file_name, argv[3]);

	if  ( fopen(file_name, "r") == NULL )
	{
		printf("Warning: Could not open song file %s\n", file_name);
		printf("Write modified dcp file anyway?\n");
		ret_char = (char)getchar();
		if ( ret_char != 'y' ) exit(1);
	}
	strcat(new_spec, " ");
	strcat(new_spec, "p=");
	strcat(new_spec, new_path);
	strcat(new_spec, " ");
	strcat(new_spec, argv[3]);
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


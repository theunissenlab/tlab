#include <stdio.h>
#include <malloc.h>
#include <string.h>

typedef struct {
	char stim_name[1024];
	char stim_id[256];
	char band[256];
	char phase[256];
	char syll[256];
	} STIMSUB;

main(argc, argv)
int argc;
char *argv[];
{
FILE *f_in;
FILE *f_sub;
char buff[1024];
char cp_buff[1024];
char *pstim;
STIMSUB *subtable;
int ntable, nt, l;
int nl, nc;
	

	if ( (f_in=fopen("tt", "r")) == NULL )
	{
		printf("Error: Could not open stat input file tt\n");
		exit(1);
	}

	if ( argc != 7 ) 
	{
		printf("Error: Wrong number of arguments\n");
		exit(1);
	}

	/* Read the stimulus substitution table */
	if ( (f_sub=fopen("subtable","r")) == NULL )
	{
		printf("Error: Could not open substitution table subtable\n");
		exit(1);
	}

	ntable = 0;
	while ( fgets(buff,1024,f_sub) ) ntable ++;
	subtable = (STIMSUB *)calloc(ntable, sizeof(STIMSUB));
	rewind(f_sub);
	nt = 0;
	/* No error checking here ... */
	while ( fgets(buff,1024,f_sub) )
	{
		 strcpy(subtable[nt].stim_name,strtok(buff,"	"));
		 strcpy(subtable[nt].stim_id,strtok(NULL,"	"));
		 strcpy(subtable[nt].band,strtok(NULL,"	"));
		 strcpy(subtable[nt].phase,strtok(NULL,"	"));
		 strcpy(subtable[nt].syll,strtok(NULL,"	"));
		 l=strlen(subtable[nt].syll);
		 subtable[nt].syll[l-1]=0;
		 nt++;
	}

	fgets(buff, 1024, f_in);
	printf("Bird\tUnit id\tUnit Class\tD Con\tD Rev\tD Noise\tStim id\tBand\tPhase\tSyll\t%s", buff);
	while ( fgets(buff, 1024, f_in) )
	{
		printf("%s\t%s\t%s\t", argv[1], argv[2], argv[3]);
		if ( argv[4][0] != '?') printf("%s\t", argv[4]);
		else printf("\t");
		if ( argv[5][0] != '?') printf("%s\t", argv[5]);
		else printf("\t");
		if ( argv[6][0] != '?') printf("%s\t", argv[6]);
		else printf("\t");

		/* Find name of stimulus file for substitutions */
		strcpy(cp_buff,buff);
		strtok(cp_buff,"	");
		pstim = strtok(NULL,"	");
		/* just look for first stim file of composite for the moment */
		nl = strlen(pstim);
		for ( nc = -1; ++nc < nl; )
		if ( pstim[nc] == '+' )
		{
			pstim[nc] = 0;
			break;
		}
		for ( nt = -1; ++ nt < ntable; )
			if ( strcmp(pstim,subtable[nt].stim_name) == 0 ) break;
		if ( nt == ntable )
			printf("Other\t\t\t\t");
		else
		{
			printf("%s\t", subtable[nt].stim_id);
			if ( subtable[nt].band[0] != '?' )
				printf("%s\t", subtable[nt].band);
			else
				printf("\t");
			if ( subtable[nt].phase[0] != '?' )
				printf("%s\t", subtable[nt].phase);
			else
				printf("\t");
			if ( subtable[nt].syll[0] != '?' )
				printf("%s\t", subtable[nt].syll);
			else
				printf("\t");
		}
		printf("%s",buff);
	}
   fclose(f_in);
}

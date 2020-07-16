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
FILE *f_sg;
FILE *f_sub;
char buff[1024];
char cp_buff[1024];
char *pstim;
STIMSUB *subtable;
int ntable, nt1, nt2, l;
int sfsigma;
int nl, nc;
	

	if ( (f_in=fopen("tt", "r")) == NULL )
	{
		printf("Error: Could not open stat input file tt\n");
		exit(1);
	}

	if ( (f_sg=fopen("sfsigma", "r")) == NULL )
	{
		printf("Error: Could not open sfsigma file\n");
		exit(1);
	}
	fscanf(f_sg,"%d", &sfsigma);
	fclose(f_sg);

	if ( argc != 4 ) 
	{
		printf("Error: Missing arguments\n");
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
	nt1 = 0;
	/* No error checking here ... */
	while ( fgets(buff,1024,f_sub) )
	{
		 strcpy(subtable[nt1].stim_name,strtok(buff,"	"));
		 strcpy(subtable[nt1].stim_id,strtok(NULL,"	"));
		 strcpy(subtable[nt1].band,strtok(NULL,"	"));
		 strcpy(subtable[nt1].phase,strtok(NULL,"	"));
		 strcpy(subtable[nt1].syll,strtok(NULL,"	"));
		 l=strlen(subtable[nt1].syll);
		 subtable[nt1].syll[l-1]=0;
		 nt1++;
	}

	fgets(buff, 1024, f_in);
	/* printf("Bird\tUnit id\tUnit Class\tsfsigma\tSType 1\tStype 2\t%s", buff); */

   nt1 = ntable;
	nt2 = ntable;
	while ( fgets(buff, 1024, f_in) )
	{
		printf("%s\t%s\t%s\t%d\t", argv[1], argv[2], argv[3], sfsigma);

		/* Find name of stimulus file for substitutions */
		if ( buff[0] != '\t' )
		{
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
			for ( nt1 = -1; ++ nt1 < ntable; )
				if ( strcmp(pstim,subtable[nt1].stim_name) == 0 ) break;
			pstim = strtok(NULL,"	");
			pstim = strtok(NULL,"	");
			nl = strlen(pstim);
			for ( nc = -1; ++nc < nl; ) 
				if ( pstim[nc] == '+' ) 
				{
					pstim[nc] = 0;
				   break;
				}
			for ( nt2 = -1; ++ nt2 < ntable; )
				if ( strcmp(pstim,subtable[nt2].stim_name) == 0 ) break;
		}
		if ( nt1 == ntable )
			printf("Other\t");
		else
			printf("%s\t", subtable[nt1].stim_id);
		if ( nt2 == ntable )
			printf("Other\t");
		else
			printf("%s\t", subtable[nt2].stim_id);
		printf("%s",buff);
	}
   fclose(f_in);
}

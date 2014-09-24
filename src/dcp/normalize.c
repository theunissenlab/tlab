#include <stdio.h>
#include <malloc.h>
#include <string.h>

typedef struct {
	char bird[16];
	char unit_id[16];
	float zbos;
	} ZBOS;

double atof();

main(argc, argv)
int argc;
char *argv[];
{
#define NVAR 30
FILE *f_in;
char buff[2048];
char cbuff[2048];
char *cvar[NVAR];
ZBOS *bosdata;
int nbos, nb, nv;
int found, l, off;
	

	if ( argc != 2 ) 
	{
		printf("Error: Wrong number of arguments\n");
		exit(1);
	}

	if ( (f_in=fopen(argv[1], "r")) == NULL )
	{
		printf("Error: Could not open input file %s\n", argv[1]);
		exit(1);
	}

	/* Read the file to count occurrences of BOS */
	fgets(buff,2048,f_in);
	nbos = 0;
	while ( fgets(buff,2048,f_in) ) 
	{
		nv = 0;
		cvar[nv] = strtok(buff, "	");
		l = strlen(cvar[nv]);
		off = 1;
		while ( cvar[nv][l+off] == '	' ) 
		{
			cvar[nv+off]=NULL;
			off++;
		}
		nv += off-1;
		for ( ; ++nv < NVAR; )
		{
			if ( (cvar[nv] = strtok(NULL,"	")) == NULL ) 
			{
				printf("Missing fields in input file. nv=%d\n", nv);
				exit(1);
			}
			l = strlen(cvar[nv]);
			off = 1;
			while ( cvar[nv][l+off] == '	')  
			{
				if ( nv+off < NVAR )
				{
					cvar[nv+off]=NULL;
					off++;
				}
				else
					break;
			}
			nv += off-1;
		}
		if ( !strcmp(cvar[6],"Bos") ) nbos ++;
	}

	bosdata = (ZBOS *)calloc(nbos, sizeof(ZBOS));
	rewind(f_in);
	nb = 0;
	/* No error checking here ... */
	fgets(buff,2048,f_in);
	while ( fgets(buff,2048,f_in) )
	{
		nv = 0;
		cvar[nv] = strtok(buff, "	");
		l = strlen(cvar[nv]);
		off = 1;
		while ( cvar[nv][l+off] == '	')  
		{
			cvar[nv+off]=NULL;
			off ++;
		}
		nv += off-1;
		for ( ; ++nv < NVAR; ) 
		{
			cvar[nv] = strtok(NULL,"	");
			l = strlen(cvar[nv]);
			off = 1;
			while ( cvar[nv][l+off] == '	') 
			{
				if ( nv+off < NVAR )
				{
					cvar[nv+off]=NULL;
					off ++;
				}
				else break;
			}
			nv += off-1;
		}
		if ( !strcmp(cvar[6],"Bos") ) 
		{
			strcpy(bosdata[nb].bird, cvar[0]);
			strcpy(bosdata[nb].unit_id, cvar[1]);
			bosdata[nb].zbos = atof(cvar[24]);
			nb++;
		}
	}

	rewind(f_in);
	fgets(buff, 2048, f_in);
	l = strlen(buff);
	buff[l-1] = 0;
	printf("%s\tPercent Z\n",buff);
	while ( fgets(buff, 1024, f_in) )
	{
		strcpy(cbuff,buff);
		l = strlen(cbuff);
		cbuff[l-1] = 0;

		nv = 0;
		cvar[nv] = strtok(buff, "	");
		l = strlen(cvar[nv]);
		off = 1;
		while ( cvar[nv][l+off] == '	') 
		{
			cvar[nv+off]=NULL;
			off ++;
		}
		nv += off -1;
		for ( ;++nv < NVAR; ) 
		{
			cvar[nv] = strtok(NULL,"	");
			l = strlen(cvar[nv]);
			off = 1;
			while ( cvar[nv][l+off] == '	') 
			{
				if ( nv+off < NVAR )
				{
					cvar[nv+off]=NULL;
					off++;
				}
				else
					break;
			}
			nv += off-1;
		}
		found = 0;
		for ( nb = -1; ++nb < nbos; )
		{
			if ( !strcmp(bosdata[nb].bird, cvar[0]) && !strcmp(bosdata[nb].unit_id,cvar[1]) )
			{
				if ( found ) 
				{
					printf("Error: Multiple BOS entries for bird %s unit %s\n", bosdata[nb].bird, bosdata[nb].unit_id);
					exit(1);
				}
				found ++;
				printf("%s\t%f\n", cbuff, atof(cvar[24])/bosdata[nb].zbos );
			}
		}
	}
   fclose(f_in);
}

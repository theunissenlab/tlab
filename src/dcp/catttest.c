#include <stdio.h>
#include <string.h>


main(argc,argv)
int argc;
char *argv[];
{
char sys_cmd[120];
char f_name[120];
char buffer[1024];
int nf, nl, l;
FILE *f_in, *f_list;

	if ( argc < 2 ) exit(1);

   system("rm -f catttest.list");
	while ( --argc )
	{
		sprintf(sys_cmd,"find . -name \"%s\" -print >> catttest.list",argv[argc]);
		system(sys_cmd);
	}

	f_list = fopen("catttest.list","r");

   nf = 0;
	while( fgets(f_name, 120, f_list) )
	{
		l = strlen(f_name);
		f_name[l-1]=0;
		if ( (f_in = fopen(f_name,"r")) == NULL )
		{
			printf("Error could not open file %s\n", f_name);
			exit(1);
		}
		nl = 0;
		while ( fgets(buffer, 1024, f_in) )
		{
			if ( nl == 0 )
			{
				if ( nf == 0 ) printf("%s", buffer);
			}
			else
				printf("%s",buffer);
			nl++;
		}
		nf++;
		fclose(f_in);
	}

}


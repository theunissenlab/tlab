#include <stdio.h>


main()
{
int ib;
int nt=401, nband=31;
float **forward;
FILE *f_out;


   f_out=fopen("forward.test","w");
   if ( f_out == NULL )
	{
		 printf("Could not open forward.test\n");
		 exit(1);
	}
	forward = (float **)calloc(nband, sizeof(float **));
	for (ib = -1; ++ib < nband; )
	{
		forward[ib]= (float *)calloc(nt, sizeof(float *));
		if (ib == 9 ) forward[ib][200]=1;
		fwrite(forward[ib],sizeof(float),nt,f_out);
	}
	fclose(f_out);
}


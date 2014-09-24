#include <stdio.h>

main()
{
int nt=401, nb=31;
float **forward;
FILE *f_in, *f_out;
int it, ib;


   f_in = fopen("forward.filt","r");
	forward = (float **)calloc(nb, sizeof(float **));
	for (ib = -1; ++ib < nb; )
	{
		forward[ib]= (float *)calloc(nt, sizeof(float *));
		fread(forward[ib],sizeof(float),nt,f_in);
   }

   f_out = fopen("test.dat","w");
	for (it=-1;++it<nt;)
	{
		fprintf(f_out,"%g\n",forward[2][it]);
	}
	fclose(f_out);
}

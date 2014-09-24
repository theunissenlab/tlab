#include <stdio.h>



#include "sm.h"


double
minFun(const gsl_vector * v, void *params)
{

	int             i;
	OPTparams      *OPTpa = (OPTparams *) params;
	ODEparams      *ODEpa = OPTpa->ODEpa;
	double         *vdata = OPTpa->vdata;
	int             nap = ODEpa->nap, nbp = ODEpa->nbp, nfp = ODEpa->nfp;
	int n=OPTpa->n;
	double          alpha[nap], beta[nbp], filter[nfp];
	for (i = 0; i < nap; i++)
		alpha[i] = gsl_vector_get(v, i);
	for (i = 0; i < nbp; i++)
		beta[i] = gsl_vector_get(v, i + nap);
	for (i = 0; i < nfp; i++)
		filter[i] = gsl_vector_get(v, i + nap + nbp);

	double          x[n];

	ODEupdate(ODEpa, alpha, beta, filter);

	ODErun(vdata, OPTpa->tmax, OPTpa->dt, ODEpa);
	ODEfilter(vdata, OPTpa->target_ns, OPTpa->nfft, ODEpa);
	scaleData(1 / getABSmax(vdata, OPTpa->target_ns), vdata, OPTpa->target_ns);
	
	static int wn=0;
	char fname[100];
	sprintf(fname,"waves/wave%d.wav",wn);
	save_wave_file(fname,vdata, OPTpa->target_ns);
	wn++;

	ODEtarget(x, vdata, OPTpa->target_ns, OPTpa->nfft, OPTpa->ntarget, OPTpa->nslices, &OPTpa->nchunks);
	return CHI2(OPTpa->spectrum, x, OPTpa->n);
}

int
main(int argc, char **argv)
{

	/* multimin stuff */
	const gsl_multimin_fminimizer_type *T =
	gsl_multimin_fminimizer_nmsimplex2;
	gsl_multimin_fminimizer *s = NULL;
	gsl_vector     *ss, *x;
	gsl_multimin_function minex_func;
	size_t          iter = 0;
	int             status;
	double          size;

	/* loading and setting up */
	int             i, k;
	int             nap = 11;
	int             nbp = 4;
	int             nfp = 10;
	const int       m = nap + nbp + nfp;

	/* Starting point */
	x = gsl_vector_alloc(m);


	
	double          alpha[] = {-0.26, -0.1374, -0.1542, -0.1552, -0.1522, -0.1507, -0.1486, -0.1447, -0.1293, -0.1291, -0.1187};
	double          beta[] = {0.0105, 0.0105, 0.0105, 0.0105, 0.0105, 0.0105, 0.0105, 0.0105, 0.0105, 0.0105};
	double          filter[] = {1e-2, 1e-1, 1, 2, 4, 8, 1, 1.0, 0.0, 0.0};

	for (i = 0; i < nap; i++)
		gsl_vector_set(x, i, alpha[i]);
	for (i = 0; i < nbp; i++)
		gsl_vector_set(x, i + nap, beta[i]);
	for (i = 0; i < nfp; i++)
		gsl_vector_set(x, i + nap + nbp, filter[i]);

	
	/* load file and params */

	OPTparams      *opt_pa = OPTinit(argv[1], nap, alpha, nbp, beta, nfp, filter);
	const int       n = opt_pa->n;

	/* Set initial step sizes */
	ss = gsl_vector_alloc(m);
	for (i = 0; i < nap; i++)
		gsl_vector_set(ss, i, 1e-1);
	for (i = 0; i < nbp; i++)
		gsl_vector_set(ss, i + nap, 1e-1);
	for (i = 0; i < nfp; i++)
		gsl_vector_set(ss, i + nap + nbp, 1);

	/* Initialize method and iterate */
	minex_func.n = m;
	minex_func.f = minFun;
	minex_func.params = opt_pa;

	s = gsl_multimin_fminimizer_alloc(T, m);
	gsl_multimin_fminimizer_set(s, &minex_func, x, ss);

	do {
		iter++;
		status = gsl_multimin_fminimizer_iterate(s);

		if (status)
			break;

		size = gsl_multimin_fminimizer_size(s);
		status = gsl_multimin_test_size(size, 1e-2);

		if (status == GSL_SUCCESS) {
			printf("converged to minimum at\n");
		}
		for(i=0;i<m;i++) 
			printf("%e ",gsl_vector_get(s->x,i));
		printf("\n");
	}
	while (status == GSL_CONTINUE && iter < 100);

	gsl_vector_free(x);
	gsl_vector_free(ss);
	gsl_multimin_fminimizer_free(s);


}




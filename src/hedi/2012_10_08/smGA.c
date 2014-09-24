#include <stdio.h>



#include "sm.h"
#include "../ga/utils.h"
#include "../ga/ga.h"

double
fitness(double * p, void *params)
{

	int             i;
	OPTparams      *OPTpa = (OPTparams *) params;
	ODEparams      *ODEpa = OPTpa->ODEpa;
	double         *vdata = OPTpa->vdata;
	int             nap = ODEpa->nap, nbp = ODEpa->nbp, nfp = ODEpa->nfp;
	int n=OPTpa->n;
	double          alpha[nap], beta[nbp], filter[nfp];
	for (i = 0; i < nap; i++)
		alpha[i] = p[i];
	for (i = 0; i < nbp; i++)
		beta[i] = p[i +nap];
	for (i = 0; i < nfp; i++)
		filter[i] =p[i + nap + nbp];

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

	/* loading and setting up */
	int             i, k;
	int             nap = 30;
	int             nbp = 30;
	int             nfp = 10;
	const int       m = nap + nbp + nfp;

	/* Starting point */
	double p[m];
	
	double          alpha[] = {-0.26, -0.1374, -0.1542, -0.1552, -0.1522, -0.1507, -0.1486, -0.1447, -0.1293, -0.1291, -0.1187};
	double          beta[] = {0.0105, 0.0105, 0.0105, 0.0105, 0.0105, 0.0105, 0.0105, 0.0105, 0.0105, 0.0105};
	double          filter[] = {1e-2, 1e-1, 1, 2, 4, 8, 1, 1.0, 0.0, 0.0};

	for (i = 0; i < nap; i++)
		p[i]=-0.15;//alpha[i];
	for (i = 0; i < nbp; i++)
		p[i + nap]=0.0105;//beta[i];
	for (i = 0; i < nfp; i++)
		p[i + nap + nbp]= filter[i];

	
	/* load file and params */

	OPTparams      *opt_pa = OPTinit(argv[1], nap, alpha, nbp, beta, nfp, filter);

	/* population */
	int npop=100; 

	gaPop * pop=gaCreatePopulation(npop,m,5/(1.*m),0.5);

	initialize_rng(1);

	double mean[m],std[m],mustd[m];
	for(i=0;i<m;i++) 
		{ mean[i]=p[i]; std[i]=1.0*fabs(p[i]); mustd[i]=1.0*fabs(p[i]);}
	gaInitPopulation(pop,mean,std);
	int generation=0;

	for(generation=0;generation<100;generation++)
	{
		/* compute fitness*/
		double mf=0,bf;
		for(i=0;i<pop->n;i++)
		{
			gaGen * g=pop->current[i];
			pop->fitness[i].index=i;
			pop->fitness[i].fitness=fitness(g->gen,opt_pa);
			mf += pop->fitness[i].fitness;
		}
		mf /= pop->n;
		gaSelection(pop);
		bf=pop->fitness[0].fitness;
		int k,id=pop->fitness[0].index;
		printf("%e %e \n",mf,bf);
		gaMutation(pop,mustd);
		gaCrossOver(pop);
	}

}




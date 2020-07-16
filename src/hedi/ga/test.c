#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "utils.h"
#include "ga.h"


double fitness(double * p, int size)
{
	int i;
	double f=0;
	for(i=0;i<size;i++)
	{
		double x=p[i];
		f += x*x*x*x -x*x*x-20*x*x+x+1;
	}
	return f;

}

int main()
{

	int n=100; 
	int m=20;
	gaPop * pop=gaCreatePopulation(n,m,1/(1.*m),0.5);

	initialize_rng(1);
	int i;
	double mean[m],std[m],mustd[m];
	for(i=0;i<m;i++) 
		{ mean[i]=5.0; std[i]=5.0; mustd[i]=0.5;}
	gaInitPopulation(pop,mean,std);
	int generation=0;
	for(generation=0;generation<100;generation++)
	{
		/* compute fitness*/
		double mf=0,bf;
		for(i=0;i<n;i++)
		{
			gaGen * g=pop->current[i];
			pop->fitness[i].index=i;
			pop->fitness[i].fitness=fitness(g->gen,m);
//			printf("%e ",pop->fitness[i].fitness);
			mf += pop->fitness[i].fitness;
		}
		mf /= n;
		gaSelection(pop);
		bf=pop->fitness[0].fitness;
		int k,id=pop->fitness[0].index;
		printf("%e %e \n",mf,bf);
//		for(k=0;k<m;k++) printf("%e ",pop->current[id]->gen[k]);
//		printf("%e \n",fitness(pop->current[id]->gen,m));
		gaMutation(pop,mustd);
//		gaCrossOver(pop);
	}

}
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "utils.h"
#include "ga.h"


// typedef struct
// {
// 	int index;
// 	double fitness;

// } gaIdx;

// typedef struct 
// {
// 	double * gen;
// } gaGen;

// typedef struct 
// {
// 	int n; // population size */
// 	int m; // gene size */
// 	gaGen ** current; // current population*/
// 	gaGen ** newest; // new population */
// 	gaIdx * fitness; // current fitness list */
// 	double mutation;
// 	double crossover;
// } gaPop;




int compFitness(const void * a,const void *b)
{
	gaIdx * ag=(gaIdx*)a;
	gaIdx * bg=(gaIdx*)b;

	if (ag->fitness<bg->fitness)
		return -1;
	else return 1;
}


gaPop *  gaCreatePopulation(int n, int m, double mut,double cro)
{
	int i;
	gaPop * pop=(gaPop*)malloc(sizeof(gaPop));
	pop->n=n;
	pop->m=m;
	pop->mutation=mut;
	pop->crossover=cro;

	pop->current=(gaGen**)malloc(sizeof(gaGen*)*n);
	pop->newest=(gaGen**)malloc(sizeof(gaGen*)*n);
	pop->fitness=(gaIdx*)malloc(sizeof(gaIdx)*n);
	return pop;	
}

void gaInitPopulation(gaPop * pop,double * mean, double * std)
{
	int i,k;
	for(i=0;i<pop->n;i++)
	{
		pop->current[i]=(gaGen*)malloc(sizeof(gaGen));
		pop->current[i]->gen=(double*)calloc(pop->m,sizeof(double));
		pop->newest[i]=(gaGen*)malloc(sizeof(gaGen)); /* dummy population */
		pop->newest[i]->gen=(double*)calloc(pop->m,sizeof(double));

		for(k=0;k<pop->m;k++)
			pop->current[i]->gen[k]=mean[k]+std[k]*white_noise();
	}

}

void gaMutation(gaPop* pop, double * mustd)
{
	int i,k;
	for(i=0;i<pop->n;i++)
		for(k=0;k<pop->m;k++)
		 if (box_noise()<pop->mutation)
				pop->current[i]->gen[k]+= mustd[k]*white_noise();
}

void gaCrossOver(gaPop * pop)
{
	/* newest and current have been inverted so newesrt contained the unselected previous population */
	int i,k;
	for(i=0;i<pop->n;i++)
		if (box_noise()<pop->crossover)
		{
			k=floor(box_noise()*pop->n);
			gaGen * h=pop->current[i];
			gaGen * g=pop->newest[k];
			int cut=floor(box_noise()*pop->m);
			if ((cut>=1)&&(cut<pop->m-1))
				if (box_noise()<0.5)
				{
					/* we copy the first segment */
					memcpy(h->gen,g->gen,cut*sizeof(double));
				}
				else
				{
				 /* we copy the second segment */
					memcpy(h->gen+cut,g->gen+cut,(pop->m-cut)*sizeof(double));
				}
		}
}

void gaSelection(gaPop* pop)
{
	/* fitness has been computed let's sort it */
	qsort(pop->fitness,pop->n,sizeof(gaIdx),compFitness);

	int child;
	for(child=0;child<pop->n;child++)
	{
		int n=pop->n;
		// select one guy at random 
		int i=0,s=n,p=n*(n-1)/2*box_noise();
		while(s<p)
		{
			i++;
			s+=n-i;
		}
		if (i>=n) i=n-1;

		int parent=pop->fitness[i].index;
	//	printf("%d %d %e\n",i,parent,pop->fitness[i].fitness);
		
		free(pop->newest[child]->gen);
		free(pop->newest[child]);

		pop->newest[child]=(gaGen*)malloc(sizeof(gaGen));
		pop->newest[child]->gen=(double*)calloc(pop->m,sizeof(double));
		memcpy(pop->newest[child]->gen,pop->current[parent]->gen,sizeof(double)*pop->m);
	}
	gaGen ** tmp=pop->current;
	pop->current=pop->newest;
	pop->newest=tmp;
}








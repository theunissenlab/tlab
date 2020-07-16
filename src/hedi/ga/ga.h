#ifndef GA_H
#define GA_H


typedef struct
{
	int index;
	double fitness;

} gaIdx;

typedef struct 
{
	double * gen;
} gaGen;

typedef struct 
{
	int n; /* population size */
	int m; /* gene size */
	gaGen ** current; /* current population*/
	gaGen ** newest; /* new population */
	gaIdx * fitness; /* current fitness list */
	double mutation;
	double crossover;
} gaPop;

gaPop *  gaCreatePopulation(int n, int m, double mut,double cro);
void gaInitPopulation(gaPop * pop,double * mean, double * std);
void gaMutation(gaPop* pop, double * mustd);
void gaCrossOver(gaPop * pop);
void gaSelection(gaPop* pop);






#endif
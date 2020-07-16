#include <malloc.h>
#include <string.h>
#include <math.h>
#include <nrutil.h>
#include <nr2.h>
#include <time.h>
#include <stdlib.h>


#include "dcp_ft.h"
DCP_DATA *read_raw();

/* ----------------------------- free_env() --------------------------------- */
free_env(nt, nbands, pow_level)
int nt;
int nbands;
double **pow_level;
{
int nb;

	for ( nb = -1; ++nb <= nbands; ) free(pow_level[nb]);
	free(pow_level);
}
/* --------------------- read_raw() -------------------------------------- */
DCP_DATA *read_raw(datafile, moredat)
char *datafile;
int *moredat;
{
int i;
static FILE *fp=NULL;
DCP_DATA *pp;
DCP_TRIAL *ppt;
DCP_STIM *stim;
int nspike;
int cd;

  if ( *moredat == 1 ) 
  {
	  if (!(fp = fopen(datafile,"r")))
	  {
		  printf("Error could not open raw datafile %s\n", datafile);
		  return NULL;
	  }
  }

  (*moredat) += 1;
  pp = (DCP_DATA *)calloc(1,sizeof(DCP_DATA));
  pp->n=1;
  ppt = pp->ptrial = (DCP_TRIAL *)calloc(pp->n, sizeof(DCP_TRIAL));
  pp->pre = 0;
  pp->waveform = NULL;

  if ( fread(&nspike,sizeof(int),1,fp) != 1 ) 
  {
     *moredat = 0;
     fclose(fp);
     fp = NULL;
     return NULL;
  }

  ppt->nspikes = nspike;
  printf("nspikes is %d\n", nspike);
  ppt->time = (unsigned int *)calloc(nspike,sizeof(unsigned int));
  if ( fread((char *)ppt->time, sizeof(unsigned int), nspike, fp) != nspike)
  {
		fprintf(stderr, "Error reading spike events in file %s.\n",datafile);
      exit(1);
  }
  for ( i = -1; ++i < nspike; ) ppt->time[i] = ppt->time[i]*1000.0/32;
  
  /* Now read the stimulus */
  stim = (DCP_STIM *)calloc(1, sizeof(DCP_STIM));
  pp->pstim = stim;
  if ( fread(&cd,sizeof(int),1,fp) != 1 ) 
  {
		fprintf(stderr, "Premature end in file %s.\n",datafile);
      exit(1);
  }
  strcpy(stim->type,"song");
  stim->length = cd;
  printf("npoints is %d\n", cd);
  stim->dur = stim->length*1000/da_freq();
  pp->cd = stim->dur;
  stim->stimulus = (short *)calloc(cd,sizeof(short));
  if ( fread((char *)stim->stimulus, sizeof(short), cd, fp) != cd )
  {
     fprintf(stderr, "Error reading stim in file %s.\n",datafile);
     exit(1);
  }
	
  return pp;

}

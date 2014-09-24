#include <stdio.h>
#include <malloc.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <sys/types.h>
 
#include "dcp_ft.h"

unsigned int factorial(unsigned int n);

void dcp_gen(double *est_spike,  int len, char *recname, DCP_DATA *pin, int irec, int sil_window, double avg_est, double avg_est2)
{
      /* Generates a fake dcp file with no noise with given rate in subdirectory rt# and creates do_stim, do_spike, do_stim_spike, and do_for script files in subdirectory rt#/analysis */
 
     int i, it, ns, nspk, ntrials;
     double rate;
     FILE *f_test; 
     FILE *f_rate;
     char ratename[1000];
     char fname[1000];
     char  dname[1000];
     char  dostimspikecall[1000];
     char  dostimcall[1000];
     char  dospikecall[1000];
     char  doforcall[1000];
     char  dostiminitcall[1000];
     char  systemcall[1000];
     int startflag;
     FILE *fdo_for, *fdo_stim, *fdo_spike, *fdostim_spike, *fdostim_init;
     int *fake_spike, *mean_fake_spike;
     DCP_DATA *pout;
     double mean_est, var_est, rate_factor, mean_est2;
     double spike_thres;
     DCP_TRIAL *ptout;
     double fireprob[20];
     double currentrate, currentexp;
     int fact, fireflag;
     int a, amax;
     int facti, bugflag;
     double cumfireprob, randomnum;

   for(rate = 0.04; rate<=.04;rate=rate+.01)
       for(ntrials = 100; ntrials <=100;ntrials = ntrials +500)
	 {
                         sprintf(fname,"%s%g.trial%d.gen", recname, rate, ntrials);
                       
			 sprintf(dname, "rate%gntrials%d", rate, ntrials);
			 printf("%s\n", dname);
		     
	     if (irec == 0)
                {
		     strcpy(dostimcall, dname);
		     strcat(dostimcall, "/do_stim");
		     strcpy(dospikecall, dname);
		     strcat(dospikecall, "/do_spike");
		     strcpy(doforcall, dname);
		     strcat(doforcall, "/do_for");
		     strcpy(dostimspikecall, dname);
		     strcat(dostimspikecall, "/do_stim_spike");
		     strcpy(dostiminitcall, dname);
		     strcat(dostiminitcall, "/do_stim_init");

			 /*creates directory for given rate and opens up do_stim, do_spike, do_stim_spike script files in those directories*/
		
		
		    
		     strcpy(systemcall, "mkdir ");
		     strcat(systemcall, dname);
		     printf("%s\n", systemcall);
                     system(systemcall);
		  
		    
		     fdo_for=fopen(doforcall, "w");
		     if (fdo_for == NULL)
		       {printf("error: opening do_for\n");
		       exit(1);
		       }
		     fputs("/home/fet/progs/dcp_util/dcp_forward\n", fdo_for);
		     fclose(fdo_for);
		     fdo_stim=fopen(dostimcall, "w"); 
		     if (fdo_stim == NULL)
		       {printf("error: opening do_stim\n");
		       exit(1);
		       }
		     fputs("/home/fet/progs/dcp_util/dcp_stim\n", fdo_stim);
		     fclose(fdo_stim);
		     fdo_spike = fopen(dospikecall, "w");
		     if (fdo_spike == NULL)
		       {printf("error: opening do_spike\n");
		       exit(1);
		       }
		     fputs("/home/fet/progs/dcp_util/dcp_spike\n", fdo_spike);
		     fclose(fdo_spike);
		     fdostim_spike=fopen(dostimspikecall, "w");
		     if (fdostim_spike == NULL)
		       {printf("error: opening do_stim_spike\n");
		       exit(1);
		       }
		     fputs("/home/fet/progs/dcp_util/dcp_stim_spike\n", fdostim_spike);
		     fclose(fdostim_spike);
		       fdostim_init = fopen(dostiminitcall, "w");
		     if (fdostim_init == NULL)
		       {printf("error: opening do_stim_init\n");
		       exit(1);
		       }
		     fputs("/home/fet/progs/dcp_util/dcp_stim_init", fdostim_init);
		     fclose(fdostim_init);
		     
		     system("pwd");
		   }




		 
		    fake_spike = (int *)calloc(len, sizeof(int));
		    mean_fake_spike = (int *)calloc(len, sizeof(int));
			 pout = (DCP_DATA *)calloc(1,sizeof(DCP_DATA));
			 memcpy(pout, pin, sizeof(DCP_DATA));
			 pout->n = ntrials;
			 pout->ptrial = (DCP_TRIAL *)calloc(ntrials, sizeof(DCP_TRIAL));
			 
			
			 /*min and max are not used.  the following creates a spike rate of variance 10 (spikes/s)^2*/
			 mean_est = 0.0;
			 mean_est2=0.0;
			 var_est = 0.0;
			 
			 rate_factor = 0.01/sqrt(avg_est2);  /* 10 Spikes/s */
			 printf("mean_est and rate_factor are %g %g\n", mean_est, rate_factor);
			 /*the following outputs the r(t) for each rate*/	
			 sprintf(ratename, "rate%g",rate *1000);
			 if (irec == 0)
			   f_rate = fopen(ratename,"w");
			 else
			   f_rate = fopen(ratename, "a");
			
			  for ( it = -1; ++it <len; )
			     fprintf(f_rate, "%g ", rate_factor*est_spike[it] + rate );
			
			  fclose(f_rate);
			   printf("rate made\n");
			 for ( i = -1; ++i < ntrials; )
			 {
				 ptout = pout->ptrial + i;
				 ptout->nspikes = 0;
				 for ( it = -1; ++it < len; )
				 {
					 fake_spike[it] = 0;
					 for( a = 0; a < 10; a++)
					   fireprob[a] = 0.0;
					 fireflag = 1;
					 currentrate  = 0;
					 if (rate_factor*est_spike[it] + rate > 0)
					   currentrate =  rate_factor*est_spike[it] + rate;
					
					 for( a = 0; fireflag; a++)
					   { facti = factorial(a);
					   currentexp = exp(-currentrate);
					   if (a ==0)
					     fireprob[a] =currentexp/facti;
					   else
					   fireprob[a] = pow(currentrate,a)*currentexp/facti;
					   if (fireprob[a] < .001)
					     {amax = a;
					     /* fireprob[a] is so small that no more probs need be calculated*/                                             fireflag = 0;
				
					      }
				   }
					 
					
					 cumfireprob = 0.0;
					 for(a=amax; a>=0 ; a--)
					   { cumfireprob += fireprob[a];
					     randomnum = drand48();
					 if ( randomnum <= cumfireprob )
	       				       {     
					   
    		            			 fake_spike[it] +=a;
						 (ptout->nspikes) +=a;
					        }	
					   }
			}
				 bugflag  = 0;
				 if (i == ntrials)
				 printf("Trial %d has %d spikes\n", i, ptout->nspikes);
				 ptout->time = (unsigned int *)calloc(ptout->nspikes,sizeof(unsigned int));
				 ns=-1;
				 for ( it = -1; ++it < len; )
					 if ( fake_spike[it] >0) 
					   {
					   for(nspk = 1; nspk <=fake_spike[it]; nspk ++)
					     {  
						 ptout->time[++ns] = (unsigned int)(pout->pre+it-sil_window)*1000 + 10*(nspk -1);
						 /* if (fake_spike[it] >1)
						    printf("nsp is %d time is %u\n",nspk,  ptout->time[ns]);*/
						      mean_fake_spike[it] += fake_spike[it];
					 }
					   }
			 }

 f_test = fopen("test.file","w");
			 for ( it = -1; ++it < len; )
					 fprintf(f_test,"%g ", (est_spike[it]-spike_thres)*rate_factor);
			 fprintf(f_test,"\n");
			 for ( it = -1; ++it < len; )
					 fprintf(f_test,"%g ", mean_fake_spike[it]/(double)pout->n);
			 fprintf(f_test,"\n");
			 fclose(f_test);

			 write_dcp_data(pout,fname);

			 /* writes out filename for do_whatever files*/
			
		       
			 fdostim_init = fopen(dostiminitcall, "a");
			 if (fdostim_init==NULL)
			   {printf("Error:can't open do_stim_init\n");
			   exit(1);
			   }
                         fprintf(fdostim_init, "%s ", fname);
                         fclose(fdostim_init);
                        


			 pout->pstim = NULL;
			 free_dcp_memory(&pout);
			 free(fake_spike);
      }
}


unsigned int factorial(unsigned int n)
{
  if (n == 0)
    return 1;
  else return n*factorial(n-1);
}

/*
 * Routines for specifying and generating stimuli for dcp like programs.
 */
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <malloc.h>

#include "dcp_ft.h"
#define DELIMITERS " \t"

double da_freq();

DCP_STIM *setup_stimulus(stim_spec)
char *stim_spec;
{
  int c=0,i,seed;
  char tone_string[160],fullname[160], pathname[160], segname[160];
  DCP_STIM *stim=NULL;
  char *av;

  av = strtok(stim_spec, DELIMITERS);
  if ( av == NULL || strcmp(av,"stim") )  return stim;
  av = strtok(NULL, DELIMITERS);

  /* Allocate space for stim */
  stim = (DCP_STIM *)calloc(1, sizeof(DCP_STIM));
  stim->samprate = 32000;    /* by default for non song stim */

  if (!strcmp(av,"tone")) {
    strcpy(stim->type,"tone");
	 av = strtok(NULL, DELIMITERS);
    while (av) {
      if (!strncmp(av,"f=",2)) stim->tone_freq = atoi(&av[2]);
      else if (!strncmp(av,"d=",2)) stim->dur=atoi(&av[2]);
      else if (!strncmp(av,"rt=",3)) stim->rise_time = atoi(&av[3]);
      else if (!strncmp(av,"ft=",3)) stim->fall_time = atoi(&av[3]);
      else 
      {
			printf("Warning: unknown tone parameter %s\n", av);
			free(stim);
			return NULL;
		}
		av = strtok(NULL, DELIMITERS);
    }
    allocate_stimulus(stim);
    make_tone(stim->stimulus, stim->length, (float)stim->tone_freq);
    ShapeWaveform(stim->stimulus, stim->length, stim->rise_time,stim->fall_time);
    
    /* This is a hack until we normalize to the power */
    for (i=0; i<stim->length; i++)
      stim->stimulus[i] *= 0.6;

  }


  else if (!strcmp(av,"sweep")) {
    strcpy(stim->type,"sweep");
	 av = strtok(NULL, DELIMITERS);
    while (av) {
      if (!strncmp(av,"f1=",3)) stim->sweep_f1 = atoi(&av[3]);
      else if (!strncmp(av,"f2=",3)) stim->sweep_f2=atoi(&av[3]);
      else if (!strncmp(av,"d=",2)) stim->dur=atoi(&av[2]);
      else if (!strncmp(av,"rt=",3)) stim->rise_time = atoi(&av[3]);
      else if (!strncmp(av,"ft=",3)) stim->fall_time = atoi(&av[3]);
      else 
      {
			printf("Warning: unknown sweep parameter %s\n", av);
			free(stim);
			return NULL;
		}
		av = strtok(NULL, DELIMITERS);
    }
    allocate_stimulus(stim);
    make_sweep(stim->stimulus, stim->length, (float)stim->sweep_f1,
	       (float)stim->sweep_f2);
    ShapeWaveform(stim->stimulus, stim->length,
		  stim->rise_time,stim->fall_time);
    
    /* This is a hack until we normalize to the power */
    for (i=0; i<stim->length; i++)
      stim->stimulus[i] *= 0.6;
  }
  else if (!strcmp(av,"tones")) {
    strcpy(stim->type,"tones");
	 av = strtok(NULL, DELIMITERS);
    while (av) {
      if (!strncmp(av,"n=",2)) {
			stim->nfreqs = atoi(&av[2]);
			if (stim->nfreqs > MAX_TONE_COMBS) {
			  printf("The maximum number of tones is %d.", MAX_TONE_COMBS);
			  free(stim);
			  return NULL;
			}
			tone_string[0]='\0';
			for (i=0; i<stim->nfreqs; i++) {
			  av = strtok(NULL, DELIMITERS);
			  stim->freq_array[i] = atoi(av);
			  strcat(tone_string,av);
			  if (i<stim->nfreqs-1) strcat(tone_string," ");
			}
      }
      else if (!strncmp(av,"d=",2)) stim->dur=atoi(&av[2]);
      else if (!strncmp(av,"rt=",3)) stim->rise_time = atoi(&av[3]);
      else if (!strncmp(av,"ft=",3)) stim->fall_time = atoi(&av[3]);
      else 
      {
			printf("Warning: unknown tones parameter %s\n", av);
			free(stim);
			return NULL;
		}
		av = strtok(NULL, DELIMITERS);
    }
    allocate_stimulus(stim);
    make_ToneCombo(stim->stimulus, stim->length, 
		   stim->freq_array, stim->nfreqs);
    ShapeWaveform(stim->stimulus, stim->length, 
		  stim->rise_time, stim->fall_time);

    /* This is a hack until we normalize to the power */
    for (i=0; i<stim->length; i++)
      stim->stimulus[i] *= 0.5;
  }
  else if (!strcmp(av,"noise")) {
    strcpy(stim->type,"noise");
	 av = strtok(NULL, DELIMITERS);
    while (av) {
      if (!strncmp(av,"d=",2)) stim->dur = atoi(&av[2]);
      else if (!strncmp(av,"rt=",3)) stim->rise_time = atoi(&av[3]);
      else if (!strncmp(av,"ft=",3)) stim->fall_time = atoi(&av[3]);
      else 
      {
			printf("Warning: unknown tones parameter %s\n", av);
			free(stim);
			return NULL;
		}
		av = strtok(NULL, DELIMITERS);
    }

    /* Initialize fstnoise if it hasn't been yet */
    /* gaussian noise with width of 2 std deviations */
    init_gaussian_table(2.0);
    seed = -1;
    init_random_sequence(&seed);

    allocate_stimulus(stim);
    fast_noise(stim->stimulus, stim->length, 1);
    ShapeWaveform(stim->stimulus, stim->length,
						stim->rise_time, stim->fall_time);

    /* This is a hack until we normalize to the power */
    for (i=0; i<stim->length; i++)
     stim->stimulus[i] *= 0.1;
  }

  else if (!strcmp(av,"song")) {
    strcpy(stim->type,"song");
    av = strtok(NULL, DELIMITERS);
    while (av) {
      if (!strncmp(av,"p=",2)) strcpy(stim->songpath, &av[2]);
      else {
			strcpy(stim->songfile, av);
			break;
		}
		av = strtok(NULL, DELIMITERS);
    }

    /* Try special path before looking in the global path list */
    combine_path_and_file(stim->songpath,stim->songfile,fullname);
    if (!file_exists(fullname)) {
		 printf("Could not open song file %s\n", fullname);
		 free(stim);
		 return NULL;
    } else {
      strcpy(pathname, stim->songpath);
    }

    /* Check to see if a segment file exist with the same name
       and a .seg extension */
    strcpy(segname,fullname);
    strcat(segname,".seg");
    if ( !file_exists(segname) )
    { 
		 if (ReadSongFile(&stim->stimulus, fullname, &stim->length, &stim->samprate ) == RET_ERR) {
			free(stim);
			return NULL;
		 } 
		 stim->dur = (stim->length*1000)/stim->samprate;
    }
    else
    {
		 /* Treat it as 1-end segment file */
		 strcpy(stim->type,"seg");
		 strcpy(segname,stim->songfile);
		 strcat(segname,".seg");
       strcpy(stim->segfile, segname);
		 strcpy(stim->segcmd,"1-end");
		 strcpy(segname,"1-end");
		 if (parse_segment_command(stim,segname) == RET_ERR) {
			free(stim);
			return NULL;
		 }
	  }	

  } else if (!strcmp(av,"seg")) {
    strcpy(stim->type,"seg");
	 av = strtok(NULL, DELIMITERS);

    if (!strncmp(av,"p=",2)) 
	 {
		 strcpy(stim->songpath, &av[2]);
		 av = strtok(NULL, DELIMITERS);
	 }
    strcpy(stim->segfile, av);
	 av = strtok(NULL, DELIMITERS);
    while (av) {
      if (!strncmp(av,"rt=",3)) stim->rise_time = atoi(&av[3]);
      else if (!strncmp(av,"ft=",3)) stim->fall_time = atoi(&av[3]);
      else break;	/* anything else is the segment sequence spec */
		av = strtok(NULL, DELIMITERS);
    }
    strcpy(stim->segcmd,av);
    if (parse_segment_command(stim,av) == RET_ERR) {
		free(stim);
		return NULL;
    }


  } else if (!strcmp(av,"vicmd")) {
	  printf("vicmd not dealt with\n");
	  free(stim);
	  return NULL;

  } else if (!strcmp(av,"burst_test")) {
	  printf("burst_test not dealt with\n");
	  free(stim);
	  return NULL;
  } else {
    printf("Unknown stimulus '%s'\n",av);
	 free(stim);
	 return NULL;
  }
  return stim;
}

ShapeWaveform(buffer, len, RiseTime, FallTime)
     short *buffer;	/* data to be shaped */
     int len;		/* length of buffer */
     int RiseTime;	/* rise time in ms */
     int FallTime;	/* fall time in ms */
{
  int time, i, last_member;
  double slope, x, rt, ft, sf=da_freq();

  if (!buffer) return RET_ERR;

  rt = (double)RiseTime;
  ft = (double)FallTime;
  last_member=len-1;

  if (rt != 0.0) {
    time = (int) (slope = rt/1000.0 * sf);
    slope = 1.0 / slope;
    for (i = 0; i < time && i < len; i++) {
      x = (double) (buffer[i]);
      x *= slope * i;
      buffer[i] = (short) (x);
    }
  }

  if (ft != 0.0) {
    time = (int) (slope = ft / 1000.0 * sf);
    slope = 1.0 / slope;
    for (i = 0; i < time && i < len; i++) {
      x = (double) (buffer[last_member - i]);
      x *= slope * i;
      buffer[last_member - i] = (short) (x);
    }
  }
  return RET_OK;
}

allocate_stimulus(stim)
DCP_STIM *stim;
{
  /* compute the length of the array */
  stim->length = (stim->dur*stim->samprate)/1000;
  stim->stimulus = (short *)calloc(stim->length,sizeof(short));
  return RET_OK;
}

double da_freq()
{
   return 200000.0;
}


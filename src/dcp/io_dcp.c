/* program to read a dcp file and do do stuff with it */

#include <stdio.h>
#include <malloc.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "digital.h"
#include "dcp_ft.h"
#include "swap.h"

/* Times before which files did not have the indicated specification */
#define NO_TIME_STAMP	689325417	/* Nov  4 23:16 */
#define NO_PROT_SPEC	692487767	/* Dec 11 13:42 */
#define USING_PROPORT	708371555	/* Jun 12 10:52 (16 bits) */
#define CAL_PROPORT	709271215	/* Jun 22 20:46 (calibrated range) */


/* A lean version of the dcp read_data routine from fileio.c of dcp */
DCP_DATA *read_data(datafile)
     char *datafile;
{
  int i;
  FILE *fp;
  struct stat file_info;
  time_t last_mod_time;
  char filetype[80];
  DCP_DATA *pp;

  if (!(fp = fopen(datafile,"r"))) 
  {
	  printf("Error could not open dcp datafile %s\n", datafile);
	  return NULL;
  }

  if (stat(datafile, &file_info)) return NULL;
  last_mod_time = file_info.st_mtime;

  pp = (DCP_DATA *)calloc(1,sizeof(DCP_DATA));
  fgets(filetype,80,fp);
  if (strstr(filetype,"dcp version") != NULL) {
    sscanf(filetype,"dcp version %s\n",pp->uv_version);
  } else {
    printf("Error old version of dcp file or not a dcp file\n");
	 fclose(fp);
	 free(pp);
	 return NULL;
  }
  if ( read_dcp_data(pp,fp,datafile) == RET_ERR)
  {
	  printf("Error reading datafile %s\n",datafile);
	  fclose(fp);
	  free_dcp_trial_memory(pp);
	  free(pp);
	  return NULL;
  }
  fclose(fp);
  return pp;
}

write_dcp_header(pp,fp)
     DCP_DATA *pp;
     FILE *fp;
{
  fprintf(fp,"dcp version %s\n",pp->uv_version);
  /* Commands must be terminated by newlines since the default is
   * for the strings to end with '\0'
   */  
  fprintf(fp,"%s\n",pp->stim_spec);
  fprintf(fp,"%s\n",pp->col_spec);
  fprintf(fp,"%s\n",pp->prot_spec);
}

skip_dcp_header(fp)
     FILE *fp;
{
  char tempstr[1024];
  if (fgets(tempstr,sizeof(tempstr),fp) == NULL)	/* stim spec */ 
    return RET_ERR;
  if (fgets(tempstr,sizeof(tempstr),fp) == NULL)	/*  col spec */ 
    return RET_ERR;
  if (fgets(tempstr,sizeof(tempstr),fp) == NULL)	/* prot spec */ 
    return RET_ERR;
  return RET_OK;
}

/* This reads a header after the filetype line has been read. */
read_dcp_header(pp,fp)
     DCP_DATA *pp;
     FILE *fp;
{
  /* In dcp, command strings are assumed to have no terminating newline */
  fgets(pp->stim_spec,MAX_SPEC_LEN,fp);
  pp->stim_spec[strlen(pp->stim_spec)-1] = '\0';
  fgets(pp->col_spec,MAX_SPEC_LEN,fp);
  pp->col_spec[strlen(pp->col_spec)-1] = '\0';
  fgets(pp->prot_spec,MAX_SPEC_LEN,fp);
  pp->prot_spec[strlen(pp->prot_spec)-1] = '\0';
}

/* This parses the prots_spec command. */
prot_parse(pp, sbuff)
DCP_DATA *pp;
char *sbuff;
{
int retval;

   retval = sscanf(sbuff,"n=%d pre=%d cd=%d ici=%d rand_ici=%d adfreq=%d adgain=%d waveform=%d\
 prot_iscale=%f prot_vscale=%f ssflag=%d prot_ra_num=%d prot_ra_cmd=%f prot_ra_start=%f\
 prot_ra_dur=%f prot_ra_delay=%f",
					 &(pp->n), &(pp->pre), &(pp->cd), &(pp->ici), &(pp->rand_ici), &(pp->adfreq), &(pp->adgain),
					 &(pp->waveform), &(pp->prot_iscale), &(pp->prot_vscale), &(pp->ssflag),
					 &(pp->prot_ra_num), &(pp->prot_ra_cmd), &(pp->prot_ra_start), 
					 &(pp->prot_ra_dur), &(pp->prot_ra_delay) );
  
   if ( retval != 16 )
   {
	   printf("Error in parsing protocol description in the header\n");
		printf("sbuff=\n%s\n",sbuff);
		printf("Number of values read = %d\n", retval);
	   return RET_ERR;
	}
	else
	  return RET_OK;
}
generate_prot(pp)
DCP_DATA *pp;
{
   sprintf(pp->prot_spec,"n=%d pre=%d cd=%d ici=%d rand_ici=%d adfreq=%d adgain=%d waveform=%d\
 prot_iscale=%f prot_vscale=%f ssflag=%d prot_ra_num=%d prot_ra_cmd=%f prot_ra_start=%f\
 prot_ra_dur=%f prot_ra_delay=%f",
					 (pp->n), (pp->pre), (pp->cd), (pp->ici), (pp->rand_ici), (pp->adfreq), (pp->adgain),
					 (pp->waveform), (pp->prot_iscale), (pp->prot_vscale), (pp->ssflag),
					 (pp->prot_ra_num), (pp->prot_ra_cmd), (pp->prot_ra_start), 
					 (pp->prot_ra_dur), (pp->prot_ra_delay) );
}

free_dcp_trial_memory(pp)
DCP_DATA *pp;
{
int i;
DCP_TRIAL *pptrial;

	if ( pp->ptrial != NULL )
	{
		for ( i = -1; ++i < pp->n; )
		{
			pptrial = pp->ptrial+i;
		   if ( pptrial->time != NULL  ) 
			{
				free(pptrial->time);
			   pptrial->time = NULL;
			}
			if ( pptrial->waveform != NULL )
			{
				free(pptrial->waveform);
				pptrial->waveform = NULL;
			}
		}
		free(pp->ptrial);
		pp->ptrial = NULL;
	}
}

allocate_dcp_trial_memory(pp)
DCP_DATA *pp;
{
int i;
int size;
DCP_TRIAL *pptrial;

	free_dcp_trial_memory(pp);
	pp->ptrial = (DCP_TRIAL *)calloc(pp->n, sizeof(DCP_TRIAL));
	if ( pp->waveform )
	{
		size = (pp->adfreq/1000)*pp->cd;
		for ( i = -1; ++i < pp->n; )
		{
			pptrial = pp->ptrial +i;
			pptrial->npoints = size;
			pptrial->waveform = (short *)calloc(size, sizeof(short));
		}
	}
}

free_dcp_stimulus(pp)
DCP_DATA *pp;
{
DCP_STIM *pstim;

   if ( pp != NULL )
	{
		pstim = pp->pstim;
		if ( pstim != NULL)
		{
		  if ( pstim->stimulus != NULL ) free(pstim->stimulus);
		  if ( pstim->song_seg != NULL ) free(pstim->song_seg);
		  free(pstim);
		  pp->pstim = NULL;
		}
	}

}

free_dcp_memory(pp)
DCP_DATA **pp;
{

	free_dcp_trial_memory(*pp);
	free_dcp_stimulus(*pp);
	if ( *pp != NULL ) free(*pp);
	*pp=NULL;
}

read_dcp_data(pp,fp,filename)
     FILE *fp;
     DCP_DATA *pp;
     char *filename;
{
  short *swp, *ewp;
  int i, j;
  int init_flg;
  unsigned int size;
  DCP_TRIAL *pptrial;

  read_dcp_header(pp,fp);

  /* Setup the appropriate stimulus. */
  if ( prot_parse(pp, pp->prot_spec) == RET_ERR ) return RET_ERR;

  allocate_dcp_trial_memory(pp);

  /* Read in spikes and waveforms */
  for (i=0; i<pp->n; i++) {
	 pptrial = pp->ptrial + i;
    fgets(pptrial->t_stamp,MAX_SPEC_LEN,fp);

    if (pp->ssflag) 
    {
      if (!feof(fp)) {
      if ( i == 0 ) init_flg = 1;
      else init_flg = 0;
		if (read_ss(pp,pptrial,fp,init_flg) == RET_ERR) {
        fprintf(stderr,"Resetting nreps=%d\n",i);
        pp->n=i;
        pptrial->ss=NULL;
        break;
     }
     }
	 }
	 else
	 {
		 fread((char *)&size, sizeof(unsigned int), 1, fp);
		 swap_u32bit(&size); 
		 pptrial->nspikes=size;
		 pptrial->time = NULL;
		 if (size > 0) 
		 {
			 pptrial->time = (unsigned int *)calloc(size,sizeof(unsigned int));
			 if ( pptrial->time == NULL  )
			 {
				 fprintf(stderr, "Error allocating memory for spikes in file %s.\n",filename);
				 fprintf(stderr, "	Collection %d asking for %d spikes\n",i,size);
				 return RET_ERR;
			 }
			 if (!feof(fp)) 
			 {
				 if (fread((char *)pptrial->time, sizeof(unsigned int), size, fp) != size)
				 {
					fprintf(stderr, "Error reading events in file %s.\n",filename);
					return RET_ERR;
				  }
				  else
				  {
					  for ( j = -1; ++j < size; ) swap_u32bit(&(pptrial->time[j])); 
				  }
			  } 
			  else 
			  {
				  fprintf(stderr,"\nCollection terminated prematurely in file %s.  nreps=%d\n",filename,i);
				  pp->n=i;
				  break;
			  }
		 }
	 }

    if (pp->waveform) {
      /* read waveform buffer */
		pptrial->npoints = (pp->adfreq/1000)*pp->cd;
      if (fread((char *)pptrial->waveform, sizeof(short),
		    pptrial->npoints, fp) != pptrial->npoints) 
		{
			fprintf(stderr, "Error reading waveform buffer in file %s\n",
				filename);
			return RET_ERR;
      }
		for ( j = -1; ++j < pptrial->npoints; ) swap_s16bit(&(pptrial->waveform[j])); 
      /* Data is stored as it comes out of the amp (with positive current
       * going into the electrode).  If we want a different convention,
       * we change it once.  To change it back, it must be reread.
       */
#ifdef IFNEDDED
      if (guv_current_sign == -1 && pp->uv_prot_iscale) {
	swp = &pp->waveform[i][0];
	ewp = &pp->waveform[i][pp->waveform_length];
	while (swp != ewp) *swp++ *= -1;
      }
#endif
    }

    /* read extcmds */
#ifdef FORLATER
    read_extcmds(pp,i,fp);
#endif
	 fread((char *)&size, sizeof(unsigned int), 1, fp);
	 swap_u32bit(&size); 
	 if ( size != 0 )
	 {
		 printf("Unexpected extcmds : %d\n", size);
		 return RET_ERR;
	 }
    
  }

#ifdef IFNEDDED
  if (pp->uv_ssflag) {
    set_default_sspikes(pp);
    copy_selected_spikes(pp);
  }
#endif

  return RET_OK;
}

/* Returns 1 if file exists, 0 if not found */
char *
file_last_modify_time(filename)
     char *filename;
{
  struct stat file_info;
  static char time[26];
  time[0]='\0';

  if (stat(filename, &file_info)) {
    /* Wasn't able to get file_info */
    return &time[0];
  }
  return ctime(&file_info.st_mtime);;
}

file_exists(filename)
     char *filename;
{
  struct stat file_info;

  if (stat(filename, &file_info)) {
    /* Wasn't able to get file_info */
    return 0;
  }

  if ((file_info.st_mode & S_IFMT) == S_IFDIR) {
    fprintf(stderr, "Requested file is a directory!\n");
    return 1;
  }
  return 1;
}

remove_old_file(filename)
     char *filename;
{
  char command[256];
  if (!file_exists(filename))
    return RET_ERR;
  
  sprintf(command,"/bin/rm %s",filename);
  system(command);
  return RET_OK;
}


write_time(pptrial,fp)
     DCP_TRIAL *pptrial;
     FILE *fp;
{
  long clock;

  clock = time((long *)0);
  strcpy(pptrial->t_stamp,ctime(&clock));
  
  /* ctime places a linefeed at the end of time string */
  fputs(pptrial->t_stamp,fp);
}

write_spike(pptrial,fp)
     DCP_TRIAL *pptrial;
     FILE *fp;
{
  int i;
  unsigned int ns, stime;

  ns = pptrial->nspikes;
  swap_u32bit(&ns); 
  if (fwrite((char *)&ns,sizeof(unsigned int),1,fp)
      < 1) {
    fprintf(stderr, "Error writing spike count.\n");
    return RET_ERR;
  }

  if (pptrial->nspikes) {
    /* write out the event buffer */
	 for ( i = -1; ++i < pptrial->nspikes; ) 
	 {
		 stime = pptrial->time[i];
		 swap_u32bit(&stime);
       if (fwrite((char *)&stime, sizeof(unsigned int), 1, fp) != 1) 
		 {
			fprintf(stderr, "Error writing spike buffer.\n");
			return RET_ERR;
		 }
	 }
  }
  return RET_OK;
}

read_ss(pp,pptrial,fp, initflg)
DCP_DATA *pp;
DCP_TRIAL *pptrial;
FILE *fp;
int initflg;
{
  char ssflag[7];
  if (initflg) {
    ss_SetSamplingRate(pp->adfreq);
    ss_Initialize();

    fgets(&ssflag[0],7,fp);
    if (strcmp(&ssflag[0],"SSFLAG") != 0) return RET_ERR;

    if (pp->spike_set) {
      ss_FreeSpikeSet(pp->spike_set);
      pp->spike_set = NULL;
    }

    if (ss_ReadSpikeSet(&pp->spike_set,fp) == 0) return RET_ERR;
  }
  pptrial->ss=NULL;
  if (ss_ReadSpikeList(&(pptrial->ss),fp) == 0) return RET_ERR;

  return RET_OK;
}

write_ss(pp,pptrial,fp,initflg)
DCP_DATA *pp;
DCP_TRIAL *pptrial;
FILE *fp;
int initflg;
{
  if (initflg) {
    fputs("SSFLAG",fp);
    if (ss_WriteSpikeSet(pp->spike_set,fp) == 0) return RET_ERR;
  }
  if (ss_WriteSpikeList(pptrial->ss,fp) == 0) return RET_ERR;

  return RET_OK;
}

write_waveform(pptrial,fp)
     DCP_TRIAL *pptrial;
     FILE *fp;
{
  int i;
  short *swp, *ewp;
  short waveform;

  /* write out waveform buffer */
  for ( i = -1; ++i < pptrial->npoints; )
  {
	  waveform = pptrial->waveform[i];
	  swap_s16bit(&waveform);
	  if (fwrite((char *)&waveform, sizeof(short), 1, fp) != 1) 
	  {
		 fprintf(stderr, "Error writing waveform buffer.\n");
		 return RET_ERR;
	  }
  }
#ifdef FORLATER
  /* Data is stored as it comes out of the amp (with positive current
   * going into the electrode).  If we want a different convention,
   * we change it once.  To change it back, it must be reread.
   */
  if (guv_current_sign == -1 && pp->uv_prot_iscale) {
    swp = &pp->waveform[i][0];
    ewp = &pp->waveform[i][pp->waveform_length];
    while (swp != ewp) *swp++ *= -1;
  }
#endif
  return RET_OK;
}

#ifdef FORLATER
read_extcmds(pp,rnum,fp)
     protocol_struct *pp;
     int rnum;	/* raster number */
     FILE *fp;
{
  int i,j;

  if (strcmp(pp->uv_version,"1.1.0") == 0) {
    fread(&pp->num_ext_cmds,sizeof(int),1,fp);

    if (rnum == 0) {		/* Only malloc the first time */
      pp->ext_cmd = (ext_cmd_struct *)
	malloc(pp->num_ext_cmds*sizeof(ext_cmd_struct));
    }
    /* Version 1.1.0 stored all extcmds with every waveform! */
    for (i=0; i<pp->num_ext_cmds; i++) {
      fread(&pp->ext_cmd[i].n,sizeof(int),1,fp);
      
      pp->ext_cmd[i].t = (float *)malloc(pp->ext_cmd[i].n*sizeof(float));
      pp->ext_cmd[i].y = (float *)malloc(pp->ext_cmd[i].n*sizeof(float));

      fread(pp->ext_cmd[i].t,sizeof(float),pp->ext_cmd[i].n,fp);
      fread(pp->ext_cmd[i].y,sizeof(float),pp->ext_cmd[i].n,fp);
    }
  } else {
    /* Version 1.1.1 and greater */
    /* Read the original number of extcmds, but discard it for now */
    fread(&pp->num_ext_cmds,sizeof(int),1,fp);

    if (pp->num_ext_cmds == 0) 
      return;

    /* Store what the external command was for each waveform, 
     * even if there is duplication. */
    pp->num_ext_cmds=pp->uv_nreps;

    if (rnum == 0) {
      pp->ext_cmd = (ext_cmd_struct *)
	malloc(pp->uv_nreps*sizeof(ext_cmd_struct));
    }
    i=rnum;
    fread(&pp->ext_cmd[i].n,sizeof(int),1,fp);

    pp->ext_cmd[i].t = (float *)malloc(pp->ext_cmd[i].n*sizeof(float));
    pp->ext_cmd[i].y = (float *)malloc(pp->ext_cmd[i].n*sizeof(float));

    fread(pp->ext_cmd[i].t,sizeof(float),pp->ext_cmd[i].n,fp);
    fread(pp->ext_cmd[i].y,sizeof(float),pp->ext_cmd[i].n,fp);
  }
}

write_extcmds(pp,rnum,fp)
     protocol_struct *pp;
     int rnum;
     FILE *fp;
{
  int i,j,cmdnum;
  fwrite(&pp->num_ext_cmds,sizeof(int),1,fp);

  if (pp->num_ext_cmds > 0)
    j=rnum%pp->num_ext_cmds;
  else 
    return;

  fwrite(&pp->ext_cmd[j].n,sizeof(int),1,fp);
  fwrite(pp->ext_cmd[j].t,sizeof(float),pp->ext_cmd[j].n,fp);
  fwrite(pp->ext_cmd[j].y,sizeof(float),pp->ext_cmd[j].n,fp);
}
#endif

/* writes data files in the current format */
write_dcp_data(pp,filename)
     DCP_DATA *pp;
     char *filename;
{
  int i;
  FILE *datafilep;
  DCP_TRIAL *pptrial;
  int ext_cmd=0;
  int init_flg;

  fprintf(stderr,"Writing dcp file %s.\n",filename);
  if (!(datafilep = fopen(filename,"w")))
  {
	  printf("Error opening output file %s\n", filename);
	  return RET_ERR;
  }

  generate_prot(pp);
  write_dcp_header(pp,datafilep);

  for (i=0; i<pp->n; i++) {

    pptrial = pp->ptrial + i;
    if (pptrial->t_stamp[0] == 0)
      write_time(pptrial,datafilep);
    else
      fputs(pptrial->t_stamp,datafilep);

	 if (pp->ssflag) {
		 if ( i ==  0 ) init_flg =1;
		 else init_flg = 0;
		 if (write_ss(pp,pptrial,datafilep, init_flg) == RET_ERR) return RET_ERR;
	  } else {
		 if (write_spike(pptrial,datafilep) == RET_ERR) return RET_ERR;
	 }
    if (pp->waveform) {
      if (write_waveform(pptrial,datafilep) == RET_ERR) return RET_ERR;
    }
#ifdef FORLATER
    write_extcmds(pp,i,datafilep);	/* save the votage output commands */
#endif
	 fwrite(&ext_cmd,sizeof(int),1,datafilep);
  }
  fclose(datafilep);
  return RET_OK;
}

combine_path_and_file(path, file, fullname)
     char *path, *file, *fullname;
{
  fullname[0]='\0';
  if (path[0]=='\0')
    getwd(path);     /* getwd requires the ucb universe */
  strcat(fullname,path);
  if (path[strlen(path)-1] != '/') strcat(fullname,"/");
  strcat(fullname,file);
  return RET_OK;
}

ReadSongFile(buffer,file,len,samprate)
     short **buffer;
     char *file;
     int *len;		/* duration in msecs (determined from file) */
	  int *samprate;
{
  FILE *fp;
  char *file_temp;
  int file_len;	
  
  file_len = strlen(file);
  file_temp = (char *)calloc(file_len+1, sizeof(char));
  strcpy(file_temp, file);

  if ( strcmp(file+file_len-4,".trf") == 0)
  {
      file_temp[file_len-4]=0;
	   if (!(fp = fopen(file_temp, "r"))) 
		{
		  if (!(fp = fopen(file, "r"))) 
		  {
			  free(file_temp);
			  return RET_ERR;
		  }
		  else
		  {
			  printf("Warning stimulus file was trf %s but non-trf file %s could not be found.\n", file, file_temp);	
		  }
		}
		else
		{
		  printf("Stimulus file was trf %s. Non-trf file %s will be used for analysis.\n", file, file_temp);	
		}
		free(file_temp);
  }
  if (!(fp = fopen(file, "r"))) return RET_ERR;
  return readin_songdata(fp,buffer,len,samprate);
}

readin_songdata(fp,buffer,len,samprate)
     FILE *fp;
     short **buffer;
     int *len;
	  int *samprate;
{
  char adstring[80];
  short *song,a,b;
  int i,ja,jb,newlen,song_ad_freq=0;
  int rlen;
  float scale;

  fgets(adstring,80,fp);
  if (strstr(adstring,"AD_FREQ:") == NULL) {
    rewind(fp);
    song_ad_freq=20000;		/* Default for backward compatability
				 * with masscomp song files */
  } else {
    sscanf(adstring,"AD_FREQ: %d Hz\n",&song_ad_freq);
  }
  *samprate=song_ad_freq;


  if ((fread((char *)len, sizeof(int), 1, fp)) == 0) return RET_ERR;
  swap_s32bit(len); 
  if (*buffer) {
    free((char *)*buffer); *buffer=NULL;
  }
 *buffer = (short *)malloc((*len+100)*sizeof(short));

  if ((rlen=fread((char *)*buffer, sizeof(short), *len, fp)) != *len) {
    fprintf(stderr,"Error reading songfile. Adf = %d len = %d read = %d\n", song_ad_freq, *len, rlen);
    return RET_ERR;
  }
  for ( i = -1; ++i < *len; ) swap_s16bit(*buffer+i); 
  fclose(fp);

  return RET_OK;
}

copy_selected_spikes(pp)
DCP_DATA *pp;
{
  DCP_TRIAL *pptrial;
  short *Mp;
  int i,j,c,nspikes,out;
  unsigned int *spikep;
  float *ssp;

  /* copy the selected spikes into the original single spike array */
  
  if ( pp->ssflag == 0 || pp->spike_set == NULL) 
    return RET_ERR;


  for (i=0; i<pp->n; i++) {
	 pptrial = pp->ptrial + i;
	 if ( pptrial->ss == NULL ) break;
    ssp = pptrial->ss->t;
    Mp = pptrial->ss->M;
    nspikes = pptrial->ss->totalspikes;
    pptrial->time
      = (unsigned int *)malloc((nspikes+2)*sizeof(unsigned int));
    spikep = &(pptrial->time[1]);
    c=0;
    for (j=0; j<nspikes; j++,Mp++,ssp++) {
      if (!sspike_selected(pp,*Mp)) continue;
      c++;
      *spikep = (int)(1000.0*(*ssp));
      spikep++;
    }
    pptrial->time[0] = c;
  }
}
sspike_selected(pp,m)
     DCP_DATA *pp;
     short m;			/* model code from spikelist->M[i] */
{
  int i,*sl, outlier_flag=pp->spike_list[pp->spike_set->nmodels+1];

  /*
   * The spike list model code is as follows:
   *    Let K equal the number of spike spike models.  The number of
   *    selected models will be <= K.
   *
   *    m = spikelist->M[i]:	
   *       0:K-1	- AP classified as spike model m
   *       K:2*K-1	- AP is an outlier but most probably spike model m
   *       -(K+1):-1	- AP is classified as spike model -(m+1)
   *                    which is not selected
   */
  
  if (m >= 2*pp->spike_set->nmodels) {
    return 0;
  } else if (m<0) {
    /* m<0 is the code for a model that was not selected during the sorting.
     * The data from these unselected models, however, is still stored.
     * Since the sorter needs to classify those models anyway.
     */
    m = -m - 1;
    return(pp->spike_list[m] && !outlier_flag);
  } else if (m >= pp->spike_set->nmodels) {
    return(pp->spike_list[m - pp->spike_set->nmodels] && outlier_flag);
  } else {
    return(pp->spike_list[m] && !outlier_flag);
  }
}


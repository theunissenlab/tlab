#include <stdio.h>
#include <string.h>
#include <math.h>

#include "swap.h"
#include "dcp_ft.h"

double da_freq();

parse_segment_command(stim,seg_cmd)
     DCP_STIM *stim;
     char *seg_cmd;
{
  char *range,*range_list,*end,*cp,*ep;
  short *stim_buffer,*temp;
  int i,j,k,c,first,last,reverse_flag,silence_flag,amplitude_ramp,
  	stim_idx,stim_len,dur_idx,
	endk,num_segs,seg_len,range_start,range_end,range_len;
  float dur;

  if (read_prot_segfile(stim) == RET_ERR) return RET_ERR;
  if (read_prot_segment_songfile(stim) == RET_ERR) return RET_ERR;

  if (!seg_cmd) return;

  range_list = (char *)malloc(2*(strlen(seg_cmd)+1)*sizeof(char));
  range = (char *)malloc(2*(strlen(seg_cmd)+1)*sizeof(char));

  strcpy(range_list,seg_cmd);

  if (end=strchr(range_list,',')) {
    *end = '\0';
    strcpy(range,range_list);
    for (cp=range_list,ep=end+1;; cp++,ep++) {
      *cp = *ep;
      if (*ep == '\0')
	break;
    }
  } else {
    strcpy(range,range_list);
    range_list[0] = '\0';
  }

  stim_len=stim_idx=0;
  while (range[0]) {
    reverse_flag=silence_flag=amplitude_ramp=0;
    dur=0.0;

    /* check for a special character (we can't use 'e' for anything) */
    switch(range[0]) {
      case 'r' :
	reverse_flag=1;
	for (cp=range; *cp!='\0'; cp++)
	  *cp = *(cp+1);
	break;
      case 'a' :
	amplitude_ramp=1;
	for (cp=range; *cp!='\0'; cp++)
	  *cp = *(cp+1);
	break;
      case 'd' :
	silence_flag=1;
	dur = atof(&range[1]);
	for (cp=range; *cp!='\0'; cp++)
	  *cp = *(cp+1);
	break;
      case 's' :
	silence_flag=1;
	for (cp=range; *cp!='\0'; cp++)
	  *cp = *(cp+1);
	break;
      case 'N' :
	fprintf(stderr,"Noise substitution not implemented yet.\n");

	if (stim_buffer) free((char *)stim_buffer);
	if (range) free((char *)range);
	if (range_list) free((char *)range_list);
	return RET_ERR;
      case 'n' :
	fprintf(stderr,"Noise envelope not implemented yet.\n");

	if (stim_buffer) free((char *)stim_buffer);
	if (range) free((char *)range);
	if (range_list) free((char *)range_list);
	return RET_ERR;
      case 'L' :
	fprintf(stderr,"Linear AM not implemented yet.\n");

	if (stim_buffer) free((char *)stim_buffer);
	if (range) free((char *)range);
	if (range_list) free((char *)range_list);
	return RET_ERR;
    }

    if (dur==0.0) {
      if (end=strrchr(range,'-')) {
	if (!strcmp(end+1,"end")) last = stim->nsegs;
	else last = atoi(end+1);
	*end = '\0';
	if (!strcmp(range,"end")) first = stim->nsegs;
	else first = atoi(range);
      } else {
	if (!strcmp(range,"end")) first = stim->nsegs;
	else first = atoi(range);
	last = first;
      }

      if (first < 0 || first > stim->nsegs || last < 0 || last > stim->nsegs) {
	fprintf(stderr,"Invalid segment range.\n");

	if (stim_buffer) free((char *)stim_buffer);
	if (range) free((char *)range);
	if (range_list) free((char *)range_list);
	return RET_ERR;
      }

      num_segs = fabs((float)(last-first))+1;
      i=first-1;
    } else {
      num_segs=1;
    }
    c=0;
    range_start=stim_idx;
    while (c<num_segs) {
      if (dur==0.0) {
	seg_len=stim->song_seg[i].len+1;	/* Larry's length is one off! */
	stim_len+=seg_len;
	if (stim_idx==0) {
	  stim_buffer=(short *)calloc((stim_len+1), sizeof(short));
	} else {
	  stim_buffer=
	    (short *)realloc(stim_buffer,(stim_len+1)*sizeof(short));
	  /* zero the new segment */
	  for (k = stim_len - seg_len; k <= stim_len; k++) {
	    stim_buffer[k] = 0;
	  }
	}

	/* copy syllable into stimulus buffer */
	if (silence_flag) {
	  for (k=stim->song_seg[i].begin; k<=stim->song_seg[i].end; k++)
	    stim_buffer[stim_idx++]=0;
	} else {
	  for (k=stim->song_seg[i].begin; k<=stim->song_seg[i].end; k++)
	    stim_buffer[stim_idx++]=stim->stimulus[k];

	  if (amplitude_ramp) {
	    ShapeWaveform(&stim_buffer[range_start],seg_len,
			  stim->rise_time,stim->fall_time);
	  }
	}
      } else {
	/* copy a fixed duration into the stimulus buffer */
	dur_idx = dur*da_freq()/1000.0;
	endk = stim_idx + dur_idx;
	stim_len += dur_idx;

	if (stim_idx==0) {
	  stim_buffer=(short *)malloc(stim_len*sizeof(short));
	} else {
	  stim_buffer=(short *)realloc(stim_buffer,stim_len*sizeof(short));
	}

	while (stim_idx <= endk) {
	  stim_buffer[stim_idx++]=0;
	}
      }

      if (first<=last) i++;
      else i--;
      c++;
    }
    range_end=stim_idx-1;

    if (reverse_flag && !silence_flag) {
      range_len = range_end-range_start+1;
      temp = (short *)malloc(range_len*sizeof(short));
      /* reverse the range just copied */
      for (j=0,k=range_start; k<=range_end; j++,k++)
	temp[j] = stim_buffer[k];
      for (j=0,k=range_end; k>=range_start; j++,k--)
	stim_buffer[k]=temp[j];
      free((char *)temp);
    }

    /* get next range from range list */
    if (range_list == NULL) break;

    if (end=strchr(range_list,',')) {
      *end = '\0';
      strcpy(range,range_list);
      for (cp=range_list,ep=end+1;; cp++,ep++) {
	*cp = *ep;
	if (*ep == '\0')
	  break;
      }
    } else {
      strcpy(range,range_list);
      range_list[0] = '\0';
    }
  }

  /* copy the stim buffer into the protocol stim buffer */

  stim->length=stim_len;
  stim->dur = (stim->length*1000)/stim->samprate;
  free(stim->stimulus);
  stim->stimulus = (short *)calloc(stim->length,sizeof(short));
  
  for (k=0; k<stim_len; k++)
    stim->stimulus[k] = stim_buffer[k];

  if (stim_buffer) free((char *)stim_buffer);
  if (range) free((char *)range);
  if (range_list) free((char *)range_list);

  return RET_OK;
}

read_prot_segfile(stim)
     DCP_STIM *stim;
{
  char str[80], fullname[160], altpath[160], *label, firstname[160];
  char *filename, segsong[160], pathname[160];
  int i;
  FILE *segfp, *songfile;
  seg_struct *seg;

  /* Try special path before looking in the global path list */
  combine_path_and_file(stim->songpath,stim->segfile,fullname);
  if (!file_exists(fullname)) {
      return RET_ERR;
  } else {
    strcpy(pathname, stim->songpath);
  }

  if ((segfp = fopen(fullname, "r")) == NULL) {
    fprintf(stderr, "Unable to open file \"%s\".", fullname);
    return RET_ERR;
  }

  /* read songfile name that this segment file was created from */
  fgets(segsong, 80, segfp);
  segsong[strlen(segsong)-1] = '\0';

  /* remove preceeding pathnames from the filename */
  if (strrchr(segsong,'/')) {
    filename=strrchr(segsong,'/');
    strcpy(stim->songfile,filename+1);
  } else {
    strcpy(stim->songfile,segsong);
  }
  

  fread((char *)&(stim->nsegs), sizeof(int), 1, segfp);
  swap_s32bit(&(stim->nsegs)); 
  stim->song_seg = (seg_struct *)malloc(stim->nsegs*sizeof(seg_struct));

  for (i=0; i<stim->nsegs; i++)
  {
    fread((char *)&stim->song_seg[i], sizeof(seg_struct), 1, segfp);
	 seg = stim->song_seg + i;
	 
	 swap_s32bit(&(seg->begin));
	 swap_s32bit(&(seg->end));
	 swap_s32bit(&(seg->len));
	 swap_s32bit(&(seg->subbegin));
	 swap_s32bit(&(seg->subend));
	 swap_s32bit(&(seg->sublen));
  }

  fclose(segfp);

  return RET_OK;
}

read_prot_segment_songfile(stim)
     DCP_STIM *stim;
{
  static char fullname[160], pathname[160];

  combine_path_and_file(stim->songpath,stim->songfile,fullname);
  if (!file_exists(fullname)) {
	  return RET_ERR;
  } else {
      strcpy(pathname, stim->songpath);
  }

  if (ReadSongFile(&stim->stimulus, fullname, &stim->length, &stim->samprate) == RET_ERR) {
    stim->length=0;
    return RET_ERR;
  } 
  return RET_OK;
}

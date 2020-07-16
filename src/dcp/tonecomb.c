#include <stdio.h>
#include <math.h>
#include <malloc.h>

#include "digital.h"

make_ToneCombo(stim_array,len,freq,nfreqs)
     short *stim_array;
     int *freq;
{
  short *temp=NULL;
  float *array=NULL;
  float ev, val;
  double scale;
  register int i, j;
  
  temp = (short *)malloc((unsigned)len*sizeof(short));
  array = (float *)calloc((unsigned)len,sizeof(float));

  for (i=0; i<nfreqs; i++) {
    make_tone(temp, len, (float)freq[i]);
    for (j=0; j<len; j++)
      array[j] += (float)temp[j];
  }
  /* find the extreme value of the tone combination to normalize */
  ev = 0.0;
  for (i = 0; i<len; i++) {
    val = array[i];
    if (ev < fabs(val)) ev = fabs(val);
  }
  /* We should normalize to achieve some average power value */
  /* normalize the tone_combo to max */
  scale = 0.5*(float)AD_MAX_VALUE/(double)ev;
  for (i=0; i<len; i++) 
    stim_array[i] = (short)(scale*(double)array[i]);
  if (temp) {
    free((char *)temp); temp=NULL;
  }
  if (array) {
    cfree((char *)array); array=NULL;
  }
  return nfreqs;
}

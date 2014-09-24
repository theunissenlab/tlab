#include <stdio.h>
#include <math.h>

#include "digital.h"

#define TWOPI 6.283185307

double da_freq();

make_tone(tone_array, len, w)
short *tone_array; /* pointer to array to store the sine wave */
int len;	/* length of the tone_array */
float w;	/* requested frequency in Hz */
{
  int i;
  float t,samp_period,A;

  A=0.5*(float)AD_MAX_VALUE;
  samp_period = 1.0/(float)da_freq();
  for (i=0,t=0.0; i<len; i++,t+=samp_period)
    tone_array[i] = (short)(A*sin(TWOPI*w*t));
}

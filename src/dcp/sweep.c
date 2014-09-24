#include <stdio.h>
#include <math.h>

#include "digital.h"

#define TWOPI 6.283185307

double da_freq();

make_sweep(stim_array, len, w1,w2)
short *stim_array; /* pointer to array to store the sine wave */
int len;	/* length of the tone_array */
float w1,w2;	/* freq start and freq end in Hz */
{
  int i;
  float t,T,A,m,phi;

  A=0.5*(float)AD_MAX_VALUE;
  m = (w2-w1)*da_freq()/(float)len;	/* sweep rate */
  T = 1.0/(float)da_freq();
  phi=0.0;
  for (i=0,t=0.0; i<len; i++,t+=T) {
    stim_array[i] = (short)(A*sin(TWOPI*(w1 + m*t)*t + phi));
    phi -= TWOPI*m*t*T;
  }
}

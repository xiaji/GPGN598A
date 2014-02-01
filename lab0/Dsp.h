#ifndef _DSP_H
#define _DSP_H

#include <stdio.h>

class Dsp {
public:
  void smooth(float a, float* x, float* y, int n) {
    float b = 1.0f-a;
    float xnm1 = x[n-1];
    int i;
    // forward
    y[0] = x[0];
    for (i=1; i<n; ++i) 
      y[i] = a*y[i-1]+b*x[i];
    // reverse
    y[n-1] = (y[n-1]+a*xnm1)/(1.0f+a);
    for (i=n-2; i>=0; --i) 
      y[i] = a*y[i+1]+b*y[i];
  }
  float mean(float* x, int n) {
    float sum = 0.0f;
    for (int i=0; i<n; ++i) 
      sum += x[i];
    return sum/n;
  }
};

#endif

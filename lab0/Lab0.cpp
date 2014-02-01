#include "Stopwatch.h"
#include "Dsp.h"

int main() {
  printf("c++\n");
  int n = 1001;
  double maxtime = 2.0;
  float a = 0.99f;
  float* x = new float[n];
  float* y = new float[n];
  for (int i=0; i<n; ++i) x[i]=0;
  x[0] = x[n/2] = x[n-1] = 1;
  int nsmooth;
  Stopwatch sw;
  Dsp dsp;
  sw.start();
  for (nsmooth=0; sw.time()<maxtime; ++nsmooth)
    dsp.smooth(a,x,y,n);
  sw.stop();
  printf("nsmooth = %d\n", nsmooth);
  printf("   mean = %12.8f\n", dsp.mean(y,n));
  printf("   time = %12.8f\n", sw.time());
  printf(" mflops = %d\n", (int)(6.0e-6*n*nsmooth/sw.time()));
  return 0;
}

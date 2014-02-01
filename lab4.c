/* 
* recursive 2D serial & parallel smoothing filter in C
*
*/

#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <math.h>
#include <sys/time.h>
#include <assert.h>
#define M_PI 3.14159265358979323846

typedef struct watch *Stopwatch;

struct watch{

    double time, start;
    int running;

};

double GetTime()
{
    struct timeval tv;
    gettimeofday(&tv,NULL);
    return (double)(tv.tv_sec + tv.tv_usec/1000000.0);
}

Stopwatch newStopwatch(){

    Stopwatch sw;
    sw = (Stopwatch) malloc(sizeof(* sw));
    return sw;
}

void Stopwatch_start(Stopwatch sw){

    if(!sw->running){
      sw->running = 1;
      sw->start = GetTime();
    }
}

void Stopwatch_stop(Stopwatch sw){

    if(sw->running){
      sw->time = GetTime()-sw->start;
      sw->running = 0;
    }
}

double Stopwatch_time(Stopwatch sw){

    if(sw->running)
        return sw->time + GetTime()-sw->start;
    else
        return GetTime();
}

void Stopwatch_reset(Stopwatch sw){

    Stopwatch_stop(sw);
    sw->time = 0.0;
}

void Stopwatch_restart(Stopwatch sw){

    Stopwatch_reset(sw);
    Stopwatch_start(sw);
}

void Stopwatch_delete(Stopwatch sw){

    free(sw);
}

float* zeroFloat1(int n){

    float *x;
    int i;

    x = (float *) malloc(n*sizeof(float));
    for(i=0; i<n; ++i)
      x[i] = 0.0;
    return x;
}

float** zeroFloat2(int n1, int n2){

    float **x;
    int i2;

    x = (float **) malloc(n2*sizeof(float *));
    for(i2=0; i2<n2; ++i2){
      x[i2] = zeroFloat1(n1);
    }
    return x;
}

float drand()   /* uniform distribution, (0..1] */
{
  return (rand()+1.0)/(RAND_MAX+1.0);
}

float random_normal() /* normal distribution, centered on 0, std dev 1 */
{
  return sqrt(-2*log(drand())) * cos(2*M_PI*drand());
}

float* randomFloat1(int n1, float mean, float std){

    float *x;
    int i1;

    x = (float *) malloc(n1*sizeof(float));
    for(i1=0; i1<n1; ++i1)
      x[i1] = std + mean*random_normal();
    return x;
}
float** randomFloat2(int n1, int n2, float mean, float std){

  float **x;
  int i2;

  x = (float **) malloc(n2*sizeof(float *));
  for(i2=0; i2<n2; ++i2)
    x[i2] = randomFloat1(n1,mean,std);
  return x;
}

void floatFree2(float **x, int n){

    int i;
    for(i=0; i<n; ++i)
      free(x[i]);
    free(x);
}

void smooth(float a, float *x, float *y, int n){

    float b = 1.0 - a;
    float yi;
    int i;

    yi = y[0] = x[0];
    /* forward */
    for(i=1; i<n-1; ++i){
      y[i] = yi = a*yi+b*x[i];
    }
    y[n-1] = (a*y[n-2]+x[n-1])/(1.0+a);

    /* backward */
    for(i=n-2; i>=0; --i)
      y[i] = a*y[i+1]+b*y[i];
}

void smooth1(float a, float **x, float **y, int n1, int n2){

    int i;
    for (i=0; i<n2; ++i)
      smooth(a,x[i],y[i],n1);
}

void smooth1P(float a, float **x, float **y, int n1, int n2, int nthreads){

    int i;

    omp_set_num_threads(nthreads);
#pragma omp parallel for schedule(dynamic)
    for (i=0; i<n2; ++i){
      smooth(a,x[i],y[i],n1);
    }
}

void smooth2(float a, float **x, float **y, int n1, int n2){

    float b = 1.0-a;
    int i1, i2;

    for(i1=0; i1<n1; ++i1)
      y[0][i1] = x[0][i1];

    /* forward */
    for(i2=1; i2<n2-2; ++i2){
      for(i1=0; i1<n1; ++i1)
        y[i2][i1] = a*y[i2-1][i1] + b*x[i2][i1];
    }
    for(i1=0; i1<n1; ++i1)
      y[n2-1][i1] = (a*y[n2-2][i1]+x[n2-1][i1])/(1.0+a);

    /* backward */
    for(i2=n2-2; i2>=0; --i2){
      for(i1=0; i1<n1; ++i1)
        y[i2][i1] = a*y[i2+1][i1] + b*y[i2][i1];
    }
}

void smooth2P(float a, float **x, float **y, int n1, int n2, int nthreads){

    float b = 1.0-a;
    int i1, i2, ithread, i1first, i1last, mi;

    omp_set_num_threads(nthreads);

    mi = 1+(n1-1)/nthreads;
    for(i1=0; i1<n1; ++i1)
      y[0][i1] = x[0][i1];
#pragma omp parallel shared(a,x,y,n1,n2,mi) private(i1first,i1last,ithread,i1,i2)
  {
    ithread = omp_get_thread_num();
    i1first = ithread*mi;
    i1last  = (i1first+mi)<n1?(i1first+mi):n1;

    /* forward */
    for(i2=1; i2<n2-2; ++i2){
      for(i1=i1first; i1<i1last; ++i1)
          y[i2][i1] = a*y[i2-1][i1] + b*x[i2][i1];
    }
    for(i1=i1first; i1<i1last; ++i1)
      y[n2-1][i1] = (a*y[n2-2][i1]+x[n2-1][i1])/(1.0+a);

    /* backward */
    for(i2=n2-2; i2>=0; --i2){
      for(i1=i1first; i1<i1last; ++i1)
          y[i2][i1] = a*y[i2+1][i1] + b*y[i2][i1];
    }
  }
}

void smooth2P2(float a, float **x, float **y, int n1, int n2, int nthreads){

    float b = 1.0-a;
    int i1, i2;

    omp_set_num_threads(nthreads);

#pragma omp parallel for schedule(dynamic)
    for(i1=0; i1<n1; ++i1)
      y[0][i1] = x[0][i1];

    /* forward */
    for(i2=1; i2<n2-2; ++i2){
#pragma omp parallel for schedule(dynamic)
      for(i1=0; i1<n1; ++i1)
          y[i2][i1] = a*y[i2-1][i1] + b*x[i2][i1];
    }
#pragma omp parallel for schedule(dynamic)
    for(i1=0; i1<n1; ++i1)
      y[n2-1][i1] = (a*y[n2-2][i1]+x[n2-1][i1])/(1.0+a);

    /* backward */
    for(i2=n2-2; i2>=0; --i2){
#pragma omp parallel for schedule(dynamic)
      for(i1=0; i1<n1; ++i1)
          y[i2][i1] = a*y[i2+1][i1] + b*y[i2][i1];
    }
}

float mean1(float *x, int n){

    float sum = 0.0;
    int i;
    for(i=0; i<n; ++i)
      sum += x[i];
    return sum/n;
}

float mean2(float **x, int n1, int n2){

    float sum = 0.0;
    int i;
    for(i=0; i<n2; ++i){
      sum += mean1(x[i],n1);
    }
    return sum/n2;
}

void assertEqual1(float *x, float *y, int n1){
 
    int i1;
    for(i1=0; i1<n1; ++i1)
      assert(x[i1]==y[i1]);
}

void assertEqual2(float **x, float **y, int n1, int n2){

    int i1, i2;
    for(i2=0; i2<n2; ++i2)
      assertEqual1(x[i2],y[i2],n1);
}

void bench1(float **x, float **ys, float **yp, int n1, int n2){

    float a = 0.99;
    float maxtime = 2.0;
    float mflopsS, mflopsP;
    int nsmooth, nthreads, nth;
    Stopwatch sw;
    sw = newStopwatch();

#pragma omp parallel
    {
      nthreads = omp_get_num_threads();
    }

    Stopwatch_start(sw);
    for(nsmooth = 0; Stopwatch_time(sw)<maxtime; ++nsmooth)
      smooth1(a,x,ys,n1,n2);
    Stopwatch_stop(sw);
    mflopsS = 0.000006*n1*n2*nsmooth/sw->time;
    fprintf(stderr,"smooth1  mean=%10f, mflops=%12f\n",mean2(ys,n1,n2),mflopsS);

    for(nth=1; nth<=nthreads; ++nth){
      Stopwatch_restart(sw);
      for(nsmooth = 0; Stopwatch_time(sw)<maxtime; ++nsmooth)
        smooth1P(a,x,yp,n1,n2,nth);
      Stopwatch_stop(sw);
      mflopsP = 0.000006*n1*n2*nsmooth/sw->time;
      fprintf(stderr,"smooth1P mean=%10f, mflops=%12f, nthreads=%2d, speedup=%10f \n",mean2(yp,n1,n2),mflopsP,nth,mflopsP/mflopsS);
      assertEqual2(yp,ys,n1,n2);
    }

    Stopwatch_delete(sw);
    fprintf(stderr,"**********************************************************************\n");
}

void bench2(float **x, float **ys, float **yp, int n1, int n2){

    float a = 0.99;
    float maxtime = 2.0;
    float mflopsS, mflopsP;
    int nsmooth, nthreads, nth;
    Stopwatch sw;
    sw = newStopwatch();

#pragma omp parallel
    {
      nthreads = omp_get_num_threads();
    }

    Stopwatch_restart(sw);
    for(nsmooth = 0; Stopwatch_time(sw)<maxtime; ++nsmooth)
      smooth2(a,x,ys,n1,n2);
    Stopwatch_stop(sw);
    mflopsS = 0.000006*n1*n2*nsmooth/sw->time;
    fprintf(stderr,"smooth2  mean=%10f, mflops=%12f\n",mean2(ys,n1,n2),mflopsS);

    for(nth=1; nth<=nthreads; ++nth){
      Stopwatch_restart(sw);
      for(nsmooth = 0; Stopwatch_time(sw)<maxtime; ++nsmooth)
        smooth2P(a,x,yp,n1,n2,nth);
      Stopwatch_stop(sw);
      mflopsP = 0.000006*n1*n2*nsmooth/sw->time;
      fprintf(stderr,"smooth2P mean=%10f, mflops=%12f, nthreads=%2d, speedup=%10f \n",mean2(yp,n1,n2),mflopsP,nth,mflopsP/mflopsS);
      assertEqual2(yp,ys,n1,n2);
    }
    Stopwatch_delete(sw);
}

int main(){

    float **x, **yp, **ys;
    int n1 = 1001, n2 = 1003;

    x  = randomFloat2(n1,n2,0.5,0.5);
    yp = zeroFloat2(n1,n2);
    ys = zeroFloat2(n1,n2);

    fprintf(stderr,"Benchmark serial & parallel smoothing of a 2D %d X %d image in C\n\n",n1,n2);
    bench1(x,ys,yp,n1,n2);
    bench2(x,ys,yp,n1,n2);

    floatFree2(x,n2);
    floatFree2(ys,n2);
    floatFree2(yp,n2);

    return 0;
}



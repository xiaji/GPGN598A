import sys
from java.awt import *
from java.lang import *
from java.util import *
from java.nio import *
from javax.swing import *

from edu.mines.jtk.awt import *
from edu.mines.jtk.io import *
from edu.mines.jtk.mosaic import *
from edu.mines.jtk.util import Stopwatch
from edu.mines.jtk.util.ArrayMath import *

from lab3 import *

#############################################################################

def main(args):
  #smooth2()
  #goBenchmark1()
  goBenchmark2()
  #goBenchmark2S()
  #goBenchmark2T()

def goBenchmark1():
  goBenchmark(1)

def goBenchmark2():
  goBenchmark(2)

def goBenchmark2S():
  goBenchmark(3)

def goBenchmark2T():
  goBenchmark(4)


def smooth2():
  a = 0.9
  x = readImage()
  y  = like(x)
  yp = like(x)
  Dsp.smooth2(a,x,y)
  Dsp.smooth2P(a,x,yp)
  plot(y,"smooth2")
  plot(yp,"smooth2P")
  plot(sub(y,yp),"dif")

def goBenchmark(t):
  nthread = Runtime.getRuntime().availableProcessors()
  print "The machine has",nthread,"threads"
  def bench(method,name,nth):
    a = 0.9
    x = readImage()
    y = like(x)
    maxtime = 2
    sw = Stopwatch()
    sw.start()
    nsmooth = 0
    while sw.time()<maxtime:
      method(a,x,y,nth)
      nsmooth += 1
    sw.stop()
    print "----------------------------------------------------------------"
    print "thread num = ",nth
    mflops = 6.0e-6*len(x[0])*len(x)*nsmooth/sw.time()
    print name+"mflops = ", (int)(mflops+0.5)
    print "mean of y = ", Dsp.mean(y)
    return mflops
  if (t==1):
    print "Now benchmark smooth1"
    mflopsS = bench(Dsp.smooth1,"smooth1 ",1)
    for nth in range(nthread):
      mflopsP = bench(Dsp.smooth1P,"smooth1P",nth+1)
      print "speed up = ", mflopsP/mflopsS
  if (t==2):
    print "Now benchmark smooth2"
    mflopsS = bench(Dsp.smooth2,  "smooth2  ",1)
    for nth in range(nthread):
      mflopsP = bench(Dsp.smooth2P, "smooth2P ",nth+1)
      print "speed up = ", mflopsP/mflopsS
  if (t==3):
    print "Now benchmark smooth2S"
    mflopsS = bench(Dsp.smooth2S, "smooth2S ",1)
    for nth in range(nthread):
      mflopsP = bench(Dsp.smooth2SP,"smooth2SP",nth+1)
      print "speed up = ", mflopsP/mflopsS
  if (t==4):
    print "Now benchmark smooth2T"
    mflopsS = bench(Dsp.smooth2T, "smooth2T ",1)
    for nth in range(nthread):
      mflopsP = bench(Dsp.smooth2TP,"smooth2TP",nth+1)
      print "speed up = ", mflopsP/mflopsS

#############################################################################

def like(x):
  return zerofloat(len(x[0]),len(x))

def testImage1():
  n1,n2 = 51,41
  x = zerofloat(n1,n2)
  x[n2/2][n1/2] = 1.0
  return x

def testImage2():
  n1,n2 = 51,41
  x = zerofloat(n1,n2)
  x[0   ][0   ] = 1.0
  x[0   ][n1-1] = 1.0
  x[n2-1][0   ] = 1.0
  x[n2-1][n1-1] = 1.0
  return x

def readImage():
  #n1,n2 = 750,600
  n1,n2 = 500,576
  x = zerofloat(n1,n2)
  ais = ArrayInputStream("luming576_500.dat",ByteOrder.BIG_ENDIAN)
  #ais = ArrayInputStream("dave.dat",ByteOrder.BIG_ENDIAN)
  ais.readFloats(x)
  ais.close()
  return x

def plot(x,title):
  sp = SimplePlot.asPixels(x)
  sp.setTitle(title)
  sp.addColorBar()
  sp.setSize(800,850)
  sp.plotPanel.setColorBarWidthMinimum(100)

#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())


#ifndef _STOPWATCH_H
#define _STOPWATCH_H

#include <stdio.h>
#include "time.h"

class Stopwatch {
public:
  Stopwatch() {
    _running = false;
    _time = 0;
  }
  void start() {
    if (!_running) {
      _running = true;
      _start = clock();
    }
  }

  void stop() {
    if (_running) {
      _time += clock()-_start;
      _running = false;
    }
  }

  void reset() {
    stop();
    _time = 0;
  }

  void restart() {
    reset();
    start();
  }

  double time() {
    if (_running) {
      return (double)(_time+(clock()-_start))/CLOCKS_PER_SEC;
    } else {
      return (double)(_time)/CLOCKS_PER_SEC;
    }
  }

  private: bool _running;
  long _start;
  long _time;
};

#endif

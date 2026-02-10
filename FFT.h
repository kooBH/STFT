#ifndef _H_FFT_BASE_
#define _H_FFT_BASE_

#include <string.h>

class FFTbase {
private : 

public:
  virtual ~FFTbase(){};
  virtual void FFT(double*) = 0;
  virtual void FFT(double**) = 0;
  virtual void FFT(double**,int) = 0;

  virtual void iFFT(double*) = 0;
  virtual void iFFT(double**) = 0;

};

#endif
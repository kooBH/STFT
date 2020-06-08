#ifndef _H_FFT_BASE_
#define _H_FFT_BASE_

/* If defined Use FFT from Numerical recipes 
 * else Use Ooura FFT */
#define _Numerical_FFT_

class FFT_Base {
protected:
  int frame_size;
  int channels;

public:
  FFT_Base(int fr, int ch) {
    frame_size = fr;
    channels = ch;
  }
  virtual ~FFT_Base(){}

  // 2D
  virtual void FFT(double **) {}
  virtual void iFFT(double **) {}
  // 1D
  virtual void FFT(double *) {}
  virtual void iFFT(double *) {}

  virtual void SingleFFT(double *) {}
  virtual void SingleiFFT(double *) {}
};

#endif

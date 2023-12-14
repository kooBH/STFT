#ifndef _H_DFT_
#define _H_DFT_

#include <string.h>
#include <math.h>

#include "FFT.h"
// It is okay to not inherit from FFTbase (FFT.h)
class DFTbrute :public FFTbase{
private :
  int n_fft;
  int n_hfft;
  int n_channels;

  double* buf;

  // Common factor in Talyor Formula
  double** factor_re;
  double** factor_im;

  double** ifactor_re;
  double** ifactor_im;

  const double DFT_pi = 3.14159265358979323846;

public :
  inline DFTbrute(int n_fft, int n_channels);
  inline ~DFTbrute();

  inline void DFT(double*);
  inline void DFT(double**);
  inline void DFT(double **, int target_channels);

  inline void iDFT(double*);
  inline void iDFT(double**);

  // Not FFT, auxilary function to interface with STFT
  inline void FFT(double*);
  inline void FFT(double**);
  inline void FFT(double**, int target_channels);
  inline void iFFT(double*);
  inline void iFFT(double**);
};

DFTbrute::DFTbrute(int n_fft_, int n_channels_) {
  this->n_fft = n_fft_;
  this->n_hfft = n_fft_ / 2 + 1;
  this->n_channels = n_channels_;
  buf = new double[n_fft_ + 2];

  factor_re = new double* [n_fft];
  for (int i = 0; i < n_fft; i++) {
    factor_re[i] = new double[n_fft];
    for (int j = 0; j < n_fft; j++)
      factor_re[i][j] = cos(-2.0 * DFT_pi * i * j / n_fft);
  }

  factor_im = new double* [n_fft];
  for (int i = 0; i < n_fft; i++) {
    factor_im[i] = new double[n_fft];
    for (int j = 0; j < n_fft; j++)
      factor_im[i][j] = sin(-2.0 * DFT_pi * i * j / n_fft);
  }

  ifactor_re = new double* [n_fft];
  for (int i = 0; i < n_fft; i++) {
    ifactor_re[i] = new double[n_fft];
    for (int j = 0; j < n_fft; j++)
      ifactor_re[i][j] = cos(2.0 * DFT_pi * i * j / n_fft);
  }

  ifactor_im = new double* [n_fft];
  for (int i = 0; i < n_fft; i++) {
    ifactor_im[i] = new double[n_fft];
    for (int j = 0; j < n_fft; j++)
      ifactor_im[i][j] = sin(2.0 * DFT_pi * i * j / n_fft);
  }
}

DFTbrute::~DFTbrute() {
  delete[] buf;
  for (int i = 0; i < n_fft; i++) {
    delete[] factor_re[i];
    delete[] factor_im[i];
    delete[] ifactor_re[i];
    delete[] ifactor_im[i];
  }
  delete[] factor_re;
  delete[] factor_im;
  delete[] ifactor_re;
  delete[] ifactor_im;
}

void DFTbrute::DFT(double** x) {
  for (int i = 0; i < n_channels; i++) {
    DFT(x[i]);
  }
}

void DFTbrute::iDFT(double** x) {
  for (int i = 0; i < n_channels; i++) {
    iDFT(x[i]);
  }
}

void DFTbrute::DFT(double** x, int target_channels) {
  for (int i = 0; i < target_channels; i++) {
    DFT(x[i]);
  }
}


/*
  x : input signal shape of [n_hfft*2] -> complex DFT output
  e.g.
  x = [s_1 s_2 s_3 .... s_n_fft, 0.0, 0.0]
  -> x = [real_1, imag_1, real_2, imag_2, ... , real_n_hfft, imag_n_hfft]
*/
void DFTbrute::DFT(double* x) {
  memcpy(buf, x, sizeof(double) * n_fft);
  int re, im;
  for (int k = 0; k < n_hfft; k++) {
    re = 2*k;
    im = 2*k + 1;
    x[re] = 0;
    x[im] = 0;
    for (int n = 0; n < n_fft; n++) {
      x[re] += buf[n] * factor_re[k][n];
      x[im] += buf[n] * factor_im[k][n];
    }
  }
}

/*
  x : complex DFT input[n_fft+2] -> output signal shape of [n_fft] 
  e.g.
  x = [real_1, imag_1, real_2, imag_2, ... , real_n_hfft, imag_n_hfft]
  ->x = [s_1 s_2 s_3 .... s_n_fft, .0.0, 0.0]

  n_fft = 8 example

0 : 3.18345404225869 + 0.00000000000000i
1 : 0.312254329339716 + 0.664395308387280i
2 : 0.407888857519389 - 0.398019401837019i
3 : 0.501028376018184 + 0.0657693633001526i
4 : -1.07520788103573 + 0.00000000000000i
5 : 0.501028376018184 - 0.0657693633001526i
6 : 0.407888857519389 + 0.398019401837019i
7 : 0.312254329339716 - 0.664395308387280i
*/
void DFTbrute::iDFT(double* x){
  memcpy(buf, x, sizeof(double) * (n_fft+2));
  double factor = 2* DFT_pi / n_fft;
  double e_re, e_im;
  int re, im;
  int i_l, i_r;

  int k;
  int k2;

  for (int n = 0; n < n_fft; n++) {
    x[n] = 0.0;

    // k = 0
    k = 0;
    re = 2 * k;
    im = 2 * k + 1;
    x[n] += buf[re]*ifactor_re[k][n];
    x[n] += buf[re]* ifactor_im[k][n];
    
    // 1 <= k < n_hfft -1
    for (int k = 1; k < n_hfft - 1; k++) {
      re = 2 * k;
      im = 2 * k + 1;
      x[n] += buf[re] * ifactor_re[k][n] - buf[im] * ifactor_im[k][n];
      x[n] += buf[im] * ifactor_re[k][n] + buf[re] * ifactor_im[k][n];
    }
    // k = n_hfft -1
    k = n_hfft - 1;
    re = 2 * k;
    im = 2 * k + 1;
    x[n] += buf[re]*ifactor_re[k][n];
    x[n] += buf[re]*ifactor_im[k][n];

    // n_hfft <= k < n_fft
    for (int k = n_hfft; k < n_fft; k++) {
      k2 = k - n_hfft;
      re = 2 * (n_hfft - k2 - 2);
      im = 2 * (n_hfft - k2 - 2)+1;
    
      x[n] += buf[re]*ifactor_re[k][n] + buf[im] * ifactor_im[k][n];
      x[n] += -buf[im] * ifactor_re[k][n] + buf[re] * ifactor_im[k][n];
    }
    x[n] /= n_fft;
  }
}

void DFTbrute::FFT(double* x) {
    DFT(x);
}
void DFTbrute::FFT(double** x) {
    DFT(x);
}
void DFTbrute::FFT(double** x,int target_channel) {
    DFT(x,target_channel);
}

void DFTbrute::iFFT(double* x) {
    iDFT(x);
}

void DFTbrute::iFFT(double** x) {
    iDFT(x);
}



#endif
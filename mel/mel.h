#ifndef _H_MEL_
#define _H_MEL_

/*
 by https://github.com/kooBH
 
 librosa - Slaney-style
 https://librosa.org/doc/main/_modules/librosa/filters.html#mel
 confirmed exact output with librosa default mel-filterbank

 + Slaney-style
 Slaney, M. Auditory Toolbox: A MATLAB Toolbox for Auditory
        Modeling Work. Technical Report, version 2, Interval Research Corporation, 1998.
 + HTK-style
 Young, S., Evermann, G., Gales, M., Hain, T., Kershaw, D., Liu, X.,
        Moore, G., Odell, J., Ollason, D., Povey, D., Valtchev, V., & Woodland, P.
        The HTK book, version 3.4. Cambridge University, March 2009.
*/

#include <cmath>
#include <stdio.h>

class mel{

private:
  int n_mels;
  int sr;
  int half_sr;
  int nfft;
  int nhfft;

  double* freq_fft; // [nhfft]
  double* freq_mel; // [n_mels+2]
  double* diff_mel; // [n_mels+1]

  double** ramps; //[n_mels+2][nhfft]
  double** filter_bank; ////[n_mels][nhfft]
  
  const double f_sp = 200.0 / 3;
  const double logstep = 0.068751777420949; // log(6.4)/27

  inline void hz_to_mel(double* f,int size, bool htk=false);
  inline double hz_to_mel(double f, bool htk=false);
  inline void mel_to_hz(double* m,int size, bool htk=false);
  inline double mel_to_hz(double m, bool htk=false);
public : 
  inline mel(int sr_ = 16000, int nfft_ = 512, int n_mels_ = 40, int f_min_ = -1, int f_max_ = -1, bool htk = false);
  inline ~mel();

  // TODO : more variety usage

  inline void filter(double*stft,double*out);

  inline void Export(const char* path);


};

mel::mel(int sr_,int nfft_, int n_mels_, int f_min_, int f_max_, bool htk) {
  int f_min, f_max;

  sr = sr_;
  nfft = nfft_;
  n_mels = n_mels_;

  nhfft = int(nfft / 2) + 1;

  freq_fft = new double[nhfft];
  freq_mel = new double[n_mels+2];

  if(f_min_ == -1) f_min = 0;
  else f_min = f_min_;

  if(f_max_ == -1) f_max = sr / 2;
  else f_max = f_max_;

  //[M] fftfreqs  = linspace(0, fs / 2, int32(1 + nfft / 2));
  for (int i = 0; i < nhfft; i++) {
    freq_fft[i] = i * ((sr/2.0) / (double)(nhfft-1));
  }

  // # Center freqs of each FFT bin
  // fftfreqs = fft_frequencies(sr = sr, n_fft = n_fft)
  {
    double min_mel = hz_to_mel(f_min, htk);
    double max_mel = hz_to_mel(f_max, htk);
    //[M] mels  = linspace(min_mel,max_mel,n_mels+2);
    for (int i = 0; i < n_mels + 2; i++) {
      freq_mel[i] = min_mel + i*(max_mel - min_mel) / (n_mels + 2 - 1) ;
    }
    mel_to_hz(freq_mel, n_mels + 2, htk);
  }

  //[M] fdiff = diff(mel_f);
  diff_mel = new double[n_mels + 1];
  for (int i = 0; i < n_mels + 1; i++) {
    diff_mel[i] = freq_mel[i + 1] - freq_mel[i];
  }
 
  //[P] ramps = np.subtract.outer(mel_f, fftfreqs)
  //[M] ramps  = mel_f' - fftfreqs;
  ramps = new double* [n_mels + 2];
  for (int i = 0; i < n_mels + 2; i++) {
    ramps[i] = new double[nhfft];
  }

  for (int i = 0; i < nhfft; i++) {
    for (int j = 0; j < n_mels + 2; j++) {
      ramps[j][i] = freq_mel[j] - freq_fft[i];
    }
  }

  filter_bank = new double* [n_mels];
  for (int i = 0; i < n_mels; i++)
    filter_bank[i] = new double[nhfft];

  double lower = 0;
  double upper = 0;
  for (int i = 0; i < n_mels; i++) {
    for (int j = 0; j < nhfft; j++) {
      // lower and upper slopes for all bins
      lower = -ramps[i][j] / diff_mel[i];
      upper = ramps[i + 2][j] / diff_mel[i + 1];

      // .. then intersect them with each other and zero
      filter_bank[i][j] = std::fmax(0.0, std::fmin(lower, upper));
    }
  }

  // librosa : If ¡®slaney¡¯, divide the triangular mel weights by the width of the mel band (area normalization).
  // torch : Slaney - style mel is scaled to be approx constant energy per channel
  //if (!htk) {
  if (!htk) {
    double enorm = 0.0;
    for (int i = 0; i < n_mels; i++) {
      enorm = 2.0 / (freq_mel[i + 2] - freq_mel[i]);
      for (int j = 0; j < nhfft; j++)
        filter_bank[i][j] *= enorm;
    }
  }

  /*
  for (int i = 0; i < sr; i += 100)
    printf("hz_to_mel(%d) : %lf, %lf\n",i,hz_to_mel(i,true),hz_to_mel(i,false));

  for (int i = 0; i < n_mels; i +=10)
    printf("mel_to_hz(%d) : %lf, %lf\n",i,mel_to_hz(i,true),mel_to_hz(i,false));
  */

  //printf("%lf %lf\n", hz_to_mel(1500, true), mel_to_hz(30, true));
  //printf("%lf %lf\n", hz_to_mel(1500, false), mel_to_hz(30, false));
}

mel::~mel() {
  delete[] freq_fft;
  delete[] freq_mel;
  delete[] diff_mel;

  for (int i = 0; i < n_mels + 2;i++)
    delete[] ramps[i];
  delete[] ramps;
  for (int i = 0; i < n_mels;i++)
    delete[] filter_bank[i];
  delete[] filter_bank;



}


void mel::filter(double* mag, double* out) {
  // TODO : sparse matrix operation
  for (int i = 0; i < n_mels; i++) {
    out[i] = 0;
    for (int j = 0; j < nhfft; j++) {
      out[i] += mag[j] * filter_bank[i][j];
    }
  }

}

void mel::hz_to_mel(double* f, int size, bool htk ) {
  if (htk) {
    for (int i = 0; i < size; i++) 
      f[i] = 2595.0 * std::log10(1.0 + f[i] / 700.0);
    return;
  }

  // fill in the linear part
  double f_min = 0.0;

  //  Fill in the log-scale part
  double min_log_hz = 1000.0; // beginning of log region (Hz)
  double min_log_mel = (min_log_hz - f_min) / f_sp; // same (Melz)

  for (int i = 0; i < size; i++) {
    if (f[i] >= min_log_hz) {
      f[i] = min_log_mel + std::log(f[i] / min_log_hz) / logstep;
    }
    else {
      f[i] = (f[i] - f_min) / f_sp;
    
    }
  }
}

double mel::hz_to_mel(double f, bool htk ) {

  if (htk)
      return 2595.0 * std::log10(1.0 + f/700.0);

  // fill in the linear part
  double f_min = 0.0;

  //  Fill in the log-scale part
  double min_log_hz = 1000.0; // beginning of log region (Hz)
  double min_log_mel = (min_log_hz - f_min) / f_sp; // same (Melz)

  if (f >= min_log_hz) {
    return min_log_mel + std::log(f / min_log_hz) / logstep;
  }
  else {
    return (f - f_min) / f_sp;
  }
}

void mel::mel_to_hz(double* m, int size, bool htk ) {
  if (htk) {
    for (int i = 0; i < size; i++)
      m[i] = 700.0 * (std::pow(10.0, (m[i] / 2595.0)) - 1.0);
    return;
  }

  // Fill in the linear scale
  double f_min = 0.0;

  // And now the nonlinear scale
  double min_log_hz = 1000.0;  // beginning of log region (Hz)
  double min_log_mel = (min_log_hz - f_min) / f_sp;  // same (Mels)

  for (int i = 0; i < size; i++) {
    if (m[i] >= min_log_mel) {
      m[i] = min_log_hz * std::exp(logstep * (m[i] - min_log_mel));
    }
    else {
      m[i] = f_min + f_sp * m[i];
    }
  }
}

double mel::mel_to_hz(double m, bool htk ) {
  if (htk) 
    return 700.0 * (std::pow(10.0, (m / 2595.0)) - 1.0);

    // Fill in the linear scale
    double f_min = 0.0;

    // And now the nonlinear scale
    double min_log_hz = 1000.0;  // beginning of log region (Hz)
    double min_log_mel = (min_log_hz - f_min) / f_sp;  // same (Mels)

    if (m >= min_log_mel) {
      return min_log_hz * std::exp(logstep * (m - min_log_mel));
    }
    else {
      return f_min + f_sp * m;
    }
}

void mel::Export(const char* path) {
  FILE*fp = fopen(path, "w");
  fprintf(fp,"mel sr %d, n_mels %d, n_hfft %d\n",sr,n_mels,nhfft);
  for (int i = 0; i < n_mels; i++) {
    for (int j = 0; j < nhfft; j++) {
      fprintf(fp,"%lf ",filter_bank[i][j]);
    }
  }
  fclose(fp);
}


#endif

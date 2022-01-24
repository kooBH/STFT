/*
 auther : https://github.com/kooBH
 
 librosa - Slaney-style
 https://librosa.org/doc/main/_modules/librosa/filters.html#mel
 confirmed exact output with librosa default mel-filterbank

 Note that there are different ways to compute MFCC.

 + reference : 
  @article{article,
  author = {Ganchev, Todor and Fakotakis, Nikos and George, Kokkinakis},
  year = {2005},
  month = {01},
  pages = {},
  title = {Comparative evaluation of various MFCC implementations on the speaker verification task},
  volume = {1},
  journal = {Proceedings of the SPECOM}
  }

*/

#include <cmath>

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

  inline void hz_to_mel(double* f,int size);
  inline double hz_to_mel(double f);
  inline void mel_to_hz(double* m,int size);
  inline double mel_to_hz(double m);
public : 
  inline mel(int,int,int);
  inline ~mel();

  // TODO : more variety usage
  inline int filter(double*stft,double*out);

};

mel::mel(int sr_=16000,int nfft_=512, int n_mels_=40) {
  sr = sr_;
  nfft = nfft_;
  n_mels = n_mels_;

  nhfft = int(nfft / 2) + 1;

  freq_fft = new double[nhfft];
  freq_mel = new double[n_mels+2];

  //[M] fftfreqs  = linspace(0, fs / 2, int32(1 + nfft / 2));
  for (int i = 0; i < nhfft; i++) {
    freq_fft[i] = i * ((sr/2) / (double)(nhfft-1));
  }

  double min_mel = hz_to_mel(0);
  double max_mel = hz_to_mel(sr/2);
  //[M] mels  = linspace(min_mel,max_mel,n_mels+2);
  for (int i = 0; i < n_mels + 2; i++) {
    freq_mel[i] = min_mel + i * max_mel / (n_mels + 2 -1);
  }
  mel_to_hz(freq_mel, n_mels + 2);

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

  // Slaney - style mel is scaled to be approx constant energy per channel
  double enorm = 0.0;
  for (int i = 0; i < n_mels; i++) {
    enorm = 2.0 / (freq_mel[i + 2] - freq_mel[i]);
    for (int j = 0; j < nhfft; j++)
      filter_bank[i][j] *= enorm;
  }
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

int mel::filter(double* stft, double* out) {
  // TODO : sparse matrix operation
  for (int i = 0; i < n_mels; i++) {
    out[i] = 0;
    for (int j = 0; j < nhfft; j++) {
      out[i] += stft[j] * filter_bank[i][j];
    }
  }

  // TODO : BLAS 

}

void mel::hz_to_mel(double* f, int size) {

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

double mel::hz_to_mel(double f) {
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

void mel::mel_to_hz(double* m, int size) {

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

double mel::mel_to_hz(double m) {

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


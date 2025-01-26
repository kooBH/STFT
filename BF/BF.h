#ifndef _H_BF_H_
#define _H_BF_H_

#include "BiquadFilter.h"

#include <vector>
/*
%% 1 - HPF
fc = 120;
BW = 2;
Q = (sqrt(2^BW))/(2^BW-1);
gain = 0;
type = 6;
[b1,a1] = BiquadFilterCoeff(type, fc, fs, Q, gain);
[H1,f] = freqz(b1,a1,fftSize,fs);


%% 2 - Bandpass (peaking)
fc = 250;
BW = 0.7;
Q = (sqrt(2^BW))/(2^BW-1);
gain = -6;
type = 9;
[b2,a2] = BiquadFilterCoeff(type, fc, fs, Q, gain);
[H2,f] = freqz(b2,a2,fftSize,fs);


%% 3 - Bandpass (peaking)
fc = 5500;
BW = 0.9;
Q = (sqrt(2^BW))/(2^BW-1);
gain = 6;
type = 9;
[b3,a3] = BiquadFilterCoeff(type, fc, fs, Q, gain);
[H3,f] = freqz(b3,a3,fftSize,fs);
*/
class BF {
private : 
  std::vector<BiquadFilter*>f1, f2, f3;
  int n_channels;
public : 
  BF(int n_channels);
  ~BF();
  void Process(double**x, int n_sample);
};

BF::BF(int n_channels_){
  double Fs = 16000, Fc;
  double BW, Q, gain;
  n_channels = n_channels_;

  Fc = 120;
  BW = 2;
  Q = (std::sqrt(std::pow(2,BW)))/(std::pow(2,BW)-1);
  gain = 0;
  for (int i = 0; i < n_channels; i++) {
    BiquadFilter* tmp = new BiquadFilter(BiquadFilter::TYPE::highpass, Fc, Fs, Q, gain);
    f1.push_back(tmp);
  }

  Fc = 250;
  BW = 0.7;
  Q = (std::sqrt(std::pow(2,BW)))/(std::pow(2,BW)-1);
  gain = -6;
  for (int i = 0; i < n_channels; i++) {
    BiquadFilter* tmp = new BiquadFilter(BiquadFilter::TYPE::peak, Fc, Fs, Q, gain);
    f2.push_back(tmp);
  }

  Fc = 5500;
  BW = 0.9;
  Q = (std::sqrt(std::pow(2,BW)))/(std::pow(2,BW)-1);
  gain = 6;
  for (int i = 0; i < n_channels; i++) {
    BiquadFilter* tmp = new BiquadFilter(BiquadFilter::TYPE::peak, Fc, Fs, Q, gain);
    f3.push_back(tmp);
  }
}

BF::~BF() {
  for (auto item: f1)
    delete item;
  for (auto item: f2)
    delete item;
  for (auto item: f3)
    delete item;
  f1.clear();
  f2.clear();
  f3.clear();
}

void BF::Process(double **x, int n_sample) {
  n_sample = 10;
  for (int i = 0; i < 10; i++)
    x[0][i] = i + 1;

  for (int i = 0; i < n_channels; i++) {
    f1[i]->Filter(x[i], n_sample);
    f2[i]->Filter(x[i], n_sample);
    f3[i]->Filter(x[i], n_sample);
  }
  for (int i = 0; i < 10; i++)
    printf("%d : %lf\n",i+1,x[0][i]);
  printf("--\n");
}




#endif
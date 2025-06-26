#ifndef _H_DRP_H_
#define _H_DRP_H_

#include <stdio.h>
#include <string.h>
#include <cmath>

class DRP{
private :
  double eps = 2.2204e-16;

  int n_hop;
  int fs;
  int nWin;
  int nChannel;

  double* buf1,*buf2, *gain;
  bool* index;

  double threshold = -40;
  double ratio  = 10;
  double attackTime = 0.01;
  double releaseTime = 0.02;
  double holdTime = 0;
  double kneeWidth = 0;

  double makeUpGain = 0;

  double levelDetectionState = 0;
  double holdTimeState = 0;
  double holdTimeSamples = 0;
  // factors for attack
  double alphaA;
  double betaA;
  // factors for release
  double betaR;
  double alphaR;
  double T;

  void gain_DRC(double* x);
  void gain_DRL(double* x);
  void gain_NG(double* x);
  void gain_smoothing_DRC_DRL();
  void gain_smoothing_DRE_NG();


  //
  void init_buffer(int n_hop);
  void init_param();


public :
  enum class MODE {NONE,DRC, DRL, NG};

  DRP(int n_hop, MODE);
  DRP(int n_hop, int fs, int nWin, int nChannel, double attackTime, double releaseTime, double kneeWidth);
  ~DRP();

  void Process(double* x);

  void DRC(double* x);
  void DRL(double* x);
  void NG(double* x);

  void Set_threshold(double val);
  void Set_ratio(double val);
  void Set_kneeWidth(double val);
  void Set_makeUpGain(double val);


private: MODE mode;
};

#endif
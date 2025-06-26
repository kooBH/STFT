#include "DRP.h"

DRP::DRP(int n_hop_, int fs_, int nWin_, int nChannel_, double attackTime_, double releaseTime_, double kneeWidth_) {
  init_buffer(n_hop_);
  this->n_hop = n_hop_; 
  this->fs = fs_; 
  this->nWin = nWin_;
  this->nChannel = nChannel_;
  this->attackTime = attackTime_;
  this->releaseTime = releaseTime_;
  this->kneeWidth = kneeWidth_;

  mode = MODE::NONE;

  init_param();
}

DRP::DRP(int n_hop_, MODE mode_) {
  init_buffer(n_hop_);
  this->n_hop = n_hop_; 
  this->fs = 16000;
  this->nWin = n_hop;
  this->nChannel = 1;
  this->mode = mode_;

  if (mode == MODE::DRC) {
    //attackTime = 1e-3;
    attackTime = 1e-4;
    releaseTime = 0.2;
    kneeWidth = 5;
    threshold =  20 * std::log10(0.15);
    //ratio = 7;
    ratio = 14;
    makeUpGain = 15;
  }
  else if (mode == MODE::DRL) {
    attackTime = 1e-4;
    releaseTime = 0.1;
    kneeWidth = 0;
    threshold =  20 * std::log10(0.9);
  }
  else if (mode == MODE::NG) {
    attackTime = 0.01;
    releaseTime = 0.02;
    holdTime = 0.01;
    kneeWidth = 0;
    threshold =  20 * std::log10(0.01);
    holdTimeSamples = holdTime * fs;
  }
  else {
    printf("ERROR::Mode(%d) is not supported\n",mode_);
  }
  
  init_param();
}

void DRP::init_buffer(int n_hop_) {
  n_hop = n_hop_;

  buf1 = new double[n_hop+1];
  memset(buf1, 0, sizeof(double) * (n_hop+1));

  buf2 = new double[n_hop];
  memset(buf2, 0, sizeof(double) * n_hop);

  gain = new double[n_hop];
  memset(gain, 0, sizeof(double) * n_hop);

  index = new bool[n_hop];
  memset(index, 0, sizeof(bool) * n_hop);

}

void DRP::init_param() {
  // [M] 10.^(obj.threshold/20);            
  T = std::pow(10, threshold / 20);

  alphaA = std::exp(-std::log(9) / (attackTime * fs));
  alphaR = std::exp(-std::log(9)/ (releaseTime * fs));
  betaA = 1 - alphaA;
  betaR = 1 - alphaR;
}

DRP::~DRP() {
  delete[] buf1;
  delete[] buf2;
  delete[] gain;
  delete[] index;
}

void DRP::Process(double* x) {
  if (mode == MODE::DRC)      DRC(x);
  else if (mode == MODE::DRL) DRL(x);
  else if (mode == MODE::NG)  NG(x);
  else {
    printf("ERROR::Mode(%d) is not specified.\n",mode);
  }
}

void DRP::gain_DRC(double* x) {
  // [M] buf1 = 20 * log10(max(abs(xFrame),eps));                    
  for (int i = 0; i < n_hop; i++)
    buf1[i] = 20 * std::log10(std::fmax(std::abs(x[i]), eps));

  // [M] buf2 = buf1;
  memcpy(buf2,buf1,sizeof(double) * n_hop);

  // hard knee
  // [M] index = 2*(buf1-obj.threshold)>obj.kneeWidth;
  for (int i = 0; i < n_hop; i++)
    index[i] = 2 * (buf1[i] - threshold) > kneeWidth;
  //[M] buf2(index)  =  obj.threshold + (buf1(index) - obj.threshold)/obj.ratio;
  for (int i = 0; i < n_hop; i++)
    if (index[i])
      buf2[i] = threshold + (buf1[i] - threshold) / ratio;

  // soft knee
  // [M]  if (obj.kneeWidth ~=0)  
  //      index = 2 * abs((buf1 - obj.threshold)) <= obj.kneeWidth;
  //      buf2(index) = buf1(index) + ((1 / obj.ratio - 1) * (buf1(index) - obj.threshold + obj.kneeWidth / 2). ^ 2 . / (2 * obj.kneeWidth));
  //      end
  if (kneeWidth != 0) {
    for (int i = 0; i < n_hop; i++) {
      index[i] = 2 * std::abs((buf1[i] - threshold)) <= kneeWidth;
      if (index[i])
        buf2[i] = buf1[i] + ((1 / ratio - 1) * std::pow((buf1[i] - threshold + kneeWidth / 2),2) / (2 * kneeWidth));
    }
  }
  // gain  = buf2 - buf1;
    for (int i = 0; i < n_hop; i++)
      gain[i] = buf2[i] - buf1[i];
}

void DRP::gain_DRL(double* x) {
  // [M] buf1 = 20 * log10(max(abs(xFrame),eps));                    
  for (int i = 0; i < n_hop; i++)
    buf1[i] = 20 * std::log10(std::fmax(std::abs(x[i]), eps));

  // [M] buf2 = buf1;
  memcpy(buf2, buf1, sizeof(double) * n_hop);

  // hard knee
  // [M] index = 2*(buf1-obj.threshold)>obj.kneeWidth;
  for (int i = 0; i < n_hop; i++)
    index[i] = 2 * (buf1[i] - threshold) > kneeWidth;
  //[M] buf2(index)  = obj.threshold;
  for (int i = 0; i < n_hop; i++)
    if (index[i])
      buf2[i] = threshold;

  // soft knee
  // [M]  if (obj.kneeWidth ~=0)  
  // index = 2*abs((buf1-obj.threshold))<=obj.kneeWidth;
  // buf2(index) = buf1(index) + ((buf1(index) - obj.threshold + obj.kneeWidth / 2). ^ 2 . / (2 *   obj.kneeWidth));
  //      end
  if (kneeWidth != 0) {
    for (int i = 0; i < n_hop; i++) {
      index[i] = 2 * std::abs((buf1[i] - threshold)) <= kneeWidth;
      if (index[i])
        buf2[i] = buf1[i] + (std::pow((buf1[i] - threshold + kneeWidth / 2), 2) / (2 * kneeWidth));
    }
  }
  // gain  = buf2 - buf1;
  for (int i = 0; i < n_hop; i++)
    gain[i] = buf2[i] - buf1[i];

}

void DRP::gain_NG(double* x) {

  //[M] xAbs = abs(xFrame);
  for (int i = 0; i < n_hop; i++)
    buf1[i] = std::abs(x[i]);
   
  //[M] T = 10. ^ (obj.threshold / 20);
  // pre-computed in the constructor

  //[M] gain = zeros(obj.nWin, obj.nChannel);
  memset(gain, 0, sizeof(double) * n_hop);

  //[M] index = (xAbs - T) > 0;
  for (int i = 0; i < n_hop; i++)
    index[i] = (buf1[i] - T) > 0.0;
  //[M] gain(index) = 1;
  for (int i = 0; i < n_hop; i++)
    if (index[i]) gain[i] = 1;
}

void DRP::gain_smoothing_DRC_DRL() {
  //[M] gainSm = zeros(obj.nWin + 1, obj.nChannel);
  memset(buf1, 0, sizeof(double) * (n_hop + 1));

  // init state
  //[M]  gainSm(1, :) = obj.levelDetectionState;
  buf1[0] = levelDetectionState;

  // [M]for sample = 1 : obj.nWin
  //      for channel = 1 : obj.nChannel
  //        if gain(sample, channel) <= gainSm(sample, channel)
  //           gainSm(sample + 1, channel) = obj.alphaA * gainSm(sample, channel) + obj.betaA * gain(sample, channel);
  //        else
  //           gainSm(sample + 1, channel) = obj.alphaR * gainSm(sample, channel) + obj.betaR * gain(sample, channel);
  //        end
  //      end
  //    end
  for (int i = 0; i < n_hop; i++) {
     if(gain[i] <= buf1[i])
        buf1[i+1] = alphaA * buf1[i] + betaA * gain[i];
      else
        buf1[i+1] = alphaR * buf1[i] + betaR * gain[i];
  }

  //Update state
  //[M]gainSm = gainSm(2:end, : );
  memcpy(gain, buf1+1, sizeof(double) * n_hop);
  
  //[M]obj.levelDetectionState = gainSm(end, :);
  levelDetectionState = buf1[n_hop];
}

void DRP::gain_smoothing_DRE_NG() {
  //[M] gainSm = zeros(obj.nWin + 1, obj.nChannel);
  memset(buf1, 0, sizeof(double) * (n_hop + 1));
  //[M] gainSm(1, :) = obj.levelDetectionState;
  buf1[0] = levelDetectionState;

  //[M] attackCount = obj.holdTimeState(1, :);
  double attackCount = holdTimeState;
  //[M] lim = obj.holdTimeSamples;
  double lim = holdTimeSamples;


  //[M]  for sample = 1 : obj.nWin
  //       for channel = 1 : obj.nChannel
  //         if gain(sample, channel) == gainSm(sample, channel)
  //           gainSm(sample + 1, channel) = gainSm(sample, channel);
  //           continue;
  //         end
  //         if gain(sample, channel) < gainSm(sample, channel) % attack
  //           if attackCount(channel) < lim
  //             % Still in hold mode
  //             attackCount(channel) = attackCount(channel) + 1;
  //             gainSm(sample + 1, channel) = gainSm(sample, channel);
  //           else
  //             gainSm(sample + 1, channel) = obj.alphaA * gainSm(sample, channel) + obj.betaA * gain(sample, channel);
  //           end
  //         else
  //         % release
  //           attackCount(channel) = 0;
  //           gainSm(sample + 1, channel) = obj.alphaR * gainSm(sample, channel) + obj.betaR * gain(sample, channel);
  //         end
  //       end
  //     end
  // omit channel iteration
  for (int i = 0; i < n_hop; i++) {
    if (gain[i] == buf1[i]) {
      buf1[i + 1] = buf1[i];
      continue;
    }
    // attack
    if (gain[i] < buf1[i]) {
      if (attackCount < lim) {
        // still in hold mode
        attackCount++;
        buf1[i + 1] = buf1[i];
      }
      else {
        buf1[i + 1] = alphaA * buf1[i] + betaA * gain[i];
      }
    }
    //release
    else {
      attackCount = 0;
      buf1[i + 1] = alphaR * buf1[i] + betaR * gain[i];
    }
  }
  
  //[M] gainSm = gainSm(2:end, : );
  memcpy(gain, buf1+1, sizeof(double) * n_hop);

  // % Update state
  //[M] obj.levelDetectionState = gainSm(end, :);
  levelDetectionState = buf1[n_hop];
  //[M] obj.holdTimeState = attackCount;
  holdTimeState = attackCount;
}


// Dynamic Range Compressor
void DRP::DRC(double* x) {

  //[M] [G, obj] = mpDRP_Gain_DRC(obj, xFrame);
  gain_DRC(x);
  //[M] [G, obj] = mpDRP_Gain_Smoothing_DRC_DRL(obj, G);
  gain_smoothing_DRC_DRL();

  for (int i = 0; i < n_hop; i++) {
    //[M] G = G + obj.makeUpGain;
    gain[i] += makeUpGain;

    // %% Convert gain to linear
    //[M] linearG = 10. ^ (G / 20);
    buf2[i] = std::pow(10, gain[i] / 20);

  // %% Apply gain
  //[M] yFrame = xFrame.*linearG;
    x[i] *= buf2[i];
  }
}

// Dynamic_Range_Limiter
void DRP::DRL(double* x) {
  //[M] [G, obj] = mpDRP_Gain_DRL(obj, xFrame);
  gain_DRL(x);

  //[M] [G, obj] = mpDRP_Gain_Smoothing_DRC_DRL(obj, G);
  gain_smoothing_DRC_DRL();

  for (int i = 0; i < n_hop; i++) {
    //[M] G = G + obj.makeUpGain;
    gain[i] += makeUpGain;

    // %% Convert gain to linear
    //[M] linearG = 10. ^ (G / 20);
    buf2[i] = std::pow(10, gain[i] / 20);

    // %% Apply gain
    //[M]  yFrame = xFrame.*linearG;
    x[i] *= buf2[i];
  }

}

// Noise_Gate
void DRP::NG(double* x) {
  //[M] [G, obj] = mpDRP_Gain_NG(obj, xFrame);
  gain_NG(x);
  //[M] [G, obj] = mpDRP_Gain_Smoothing_DRE_NG(obj, G);
  gain_smoothing_DRE_NG();

  //%% Apply gain
  //[M]  yFrame = xFrame.*G;
  for (int i = 0; i < n_hop; i++) {
    x[i] *= gain[i];
  }

}


void DRP::Set_threshold(double val) {
  this->threshold = val;
}

void DRP::Set_ratio(double val) {
  this->ratio = val;
}

void DRP::Set_kneeWidth(double val) {
  this->kneeWidth = val;
}

void DRP::Set_makeUpGain(double val) {
  this->makeUpGain = val;
}


#ifndef _H_STFT_
#define _H_STFT_

#include "Ooura_FFT.h"
#include "HannWindow.h"
#include "PostProcessor.h"

class STFT{
  private : 
    HannWindow *hw;
    Ooura_FFT *fft;
    PostProcessor *ap;

    int channels;
    int frame_size;
    int shift_size;
    int ol;

    double**buf;

  public :
    inline STFT(int channels,int frame,int shift);
    inline ~STFT();
    /* in from input device or file
    
      in : raw buffer from wav or mic
      length : shift_size * channels   (for not fully occupied input)
      out : STFTed buffer [channels][frame_size + 2] (half FFT with complex)
      */
    inline void stft(short*in,int length,double**out);
    inline void istft(double**in,short*out);

    /* Single-Channel    
      in : 1 x shift
      out : 1 x frame_size + 2 (half FFT with complex)
    */
    inline void stft(short* in, double* out);
    inline void stft(double* in, double* out);
};


STFT::STFT(int channels_,int frame_,int shift_){
  int i;
  channels = channels_;
  frame_size = frame_;
  shift_size = shift_;
  ol = frame_size - shift_size;

  hw = new HannWindow(frame_size, shift_size);
  fft= new Ooura_FFT(frame_size, channels);
  ap = new PostProcessor(frame_size, shift_size, channels);

  buf =  new double*[channels];
  for(i=0;i<channels;i++){
    buf[i] = new double[frame_size];
    memset(buf[i],0,sizeof(double)*frame_size);
  }
}

STFT::~STFT(){
  int i;
  delete hw;
  delete fft;
  delete ap;
  for(i=0;i<channels;i++)
    delete[] buf[i];
  delete[] buf;

}

void STFT::stft(short*in,int length,double**out){
  int i,j;
  /*** Shfit & Copy***/
  for (j = 0; j < channels; j++) {
    for (i = 0; i < ol; i++) {
      buf[j][i] = buf[j][i + shift_size];
    }
  }
  //// EOF
  if(length!=shift_size*channels){

    length = length/channels;
    for (i = 0; i < length; i++) {
      for (j = 0; j < channels; j++)
        buf[j][i + ol]
          =  (double)(in[i * channels+ j]);
    }
    for (i = length; i < shift_size; i++) {
      for (j = 0; j < channels; j++)
        buf[j][i + ol] = 0;
    }
    //// continue
  }else{
    for (i = 0; i < shift_size; i++) {
      for (j = 0; j < channels; j++){
        buf[j][i + ol] 
          = (double)(in[i * channels+ j]);
      }
    }
  }
  /*** Copy input -> hann_input buffer ***/
  for (i = 0; i < channels; i++)
    memcpy(out[i], buf[i], sizeof(double) * frame_size);

  // scaling for precision
  for (i = 0; i < channels; i++)
    for (j = 0; j < frame_size; j++)
      out[i][j] /= 32767.0;

  /*** Window ***/
  hw->Process(out, channels);

  /*** FFT ***/
  fft->FFT(out);
}


void STFT::istft(double**in,short*out){
  /*** iFFT ***/
  fft->iFFT(in);

  /*** Window ***/
  hw->Process(in, channels);

  /*** Output ***/
  memcpy(out,ap->Overlap(in),sizeof(short)*shift_size*channels);
}


// Single-Channel double
void STFT::stft(short* in, double* out){
	int i;
    /*** Shfit & Copy***/
    for (i = 0; i < ol; i++) {
        buf[0][i] = buf[0][i + shift_size];
    }
    for (i = 0; i < shift_size; i++)
        buf[0][ol + i] = static_cast<double>(in[i]);

    memcpy(out, buf, sizeof(double) * frame_size);
    
    /*** Window ***/
    hw->Process(out);

    /*** FFT ***/
    fft->FFT(out);
}
void STFT::stft(double* in, double* out) {
    int i;
	/*** Shfit & Copy***/
    for (i = 0; i < ol; i++) {
        buf[0][i] = buf[0][i + shift_size];
    }
    for (i = 0; i < shift_size; i++)
        buf[0][ol + i] = in[i];

    memcpy(out, buf, sizeof(double) * frame_size);

    /*** Window ***/
    hw->Process(out);

    /*** FFT ***/
    fft->FFT(out);
}

#endif

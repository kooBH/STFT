#ifndef _H_STFT_
#define _H_STFT_

#include "Ooura_FFT.h"
#include "HannWindow.h"
#include "PostProcessor.h"

class STFT {
private:
    const double MATLAB_scale = 32768;
    bool opt_scale = true;

    HannWindow* hw;
    Ooura_FFT* fft;
    PostProcessor* ap;

    int channels;
    int frame_size;
    int shift_size;
    int ol;

    double** buf;

public:
    STFT(int channels, int frame, int shift);
    ~STFT();

    void set_opt_scale(bool);
    /* in from input device or file

      in : raw buffer from wav or mic
      length : shift_size * channels   (for not fully occupied input)
      out : STFTed buffer [channels][frame_size + 2] (half FFT in complex)
      */
    void stft(short* in, int length, double** out);
    void istft(double** in, short* out);

    void stft(short* in, int length, double** out, int target_channels);

    /* 2-D raw input STFT
       in  : [channels][shift_size]   raw data in double
       out : [channels][frame_size+2]
    */
    void stft(double** in, double** out);
    void istft(double** in, double** out);
    void stft(double** in, double** out, int target_channels);

    /* Single-Channel STFT
      in : 1 x shift
      out : 1 x frame_size + 2 (half FFT in complex)
    */
    void stft(short* in, double* out);
    void stft(double* in, double* out);

    /* Single-Channel ISTFT
      in : 1 x frame_size + 2 (half FFT in complex)
      out : 1 x shift_size     */
    void istft(double* in, short* out);
    void istft(double* in, double* out);

    //for separated 3-channels wav
    void stft(short* in_1, short* in_2, short* in_3, int length, double** out);
};

#endif
#include "WAV.h"
#include "STFT.h"
#include "LPF.h"

int main() {
  int n_channels = 1;
  int n_hop = 256;
  int n_fft = 512;
  int n_overlap = n_fft - n_hop;

  WAV wav_in,wav_out(n_channels,16000);
  double **x, **X, **y;
  short* buf_io;

  LPF *lpf = new LPF(n_channels);

  Ooura_FFT* fft = new Ooura_FFT(n_fft, n_channels);
  SineWindow* hw = new SineWindow(n_fft, n_hop,-1,false);
  OA*ap = new OA(n_fft, n_hop, n_channels);

  x = new double* [n_channels];
  X = new double* [n_channels];
  for (int i = 0; i < n_channels; i++) {
    x[i] = new double[n_fft];
    X[i] = new double[n_fft + 2];

    memset(x[i], 0, sizeof(double) * n_fft);
    memset(X[i], 0, sizeof(double) * (n_fft+2));
  }
  buf_io = new short[n_channels * n_hop];

  wav_in.OpenFile("input.wav");
  wav_out.NewFile("output.wav");

  while (!wav_in.IsEOF()) {
    wav_in.ReadUnit(buf_io, n_hop * n_channels);

    // STFT

    /*** Shfit & Copy***/
    for (int j = 0; j < n_channels; j++) {
      for (int i = 0; i < n_overlap; i++) {
        x[j][i] = x[j][i + n_hop];
      }
      for (int i = 0; i < n_hop; i++) {
        x[j][n_overlap + i] = static_cast<double>(buf_io[i * n_channels + j] / 32768.0);
      }
      memcpy(X[j], x[j], sizeof(double) * n_fft);
    }

    hw->Process(X, n_channels);

    /*** FFT ***/
    fft->FFT(X);


    // iSTFT
    /*** iFFT ***/
    fft->iFFT(X);

    // IIR Fiter in iSTFT
    lpf->Process(X, n_fft);

    /*** Window ***/
    hw->Process(X, n_channels);

    // iSTFT - overlap-add
    ap->Overlap(X);
    for (int i = 0; i < n_channels; i++)
      for (int j = 0; j < n_hop; j++)
        buf_io[j * n_channels + i] = static_cast<short>(ap->Get_buf()[i][j] * 32768.0);
    wav_out.Append(buf_io, n_channels * n_hop);
  }
  wav_in.Finish();
  wav_out.Finish();


  return 0;
}
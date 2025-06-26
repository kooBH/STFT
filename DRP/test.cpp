#include "WAV.h"
#include "DRP.h"

int main() {
  int n_channels = 1;
  int n_hop = 256;

  WAV wav_in,wav_out(n_channels,16000);
  double* x;
  short* buf_io;

  x = new double[n_hop];
  buf_io = new short[n_channels * n_hop];

  wav_in.OpenFile("input.wav");
  wav_out.NewFile("output.wav");

  DRP DRC(n_hop, DRP::MODE::DRC);
  DRP DRL(n_hop, DRP::MODE::DRL);
  DRP NG(n_hop, DRP::MODE::NG);


  while (!wav_in.IsEOF()) {
    wav_in.ReadUnit(buf_io, n_hop * n_channels);

    for (int i = 0; i < n_hop; i++) {
      x[i] = static_cast<double>(buf_io[i])/32768.0;
    }

    DRC.Process(x);
    DRL.Process(x);
    NG.Process(x);


    for (int i = 0; i < n_hop; i++) {
      buf_io[i] = static_cast<short>(x[i] * 32768);
    }

    wav_out.Append(buf_io, n_channels * n_hop);
  }
  wav_in.Finish();
  wav_out.Finish();


  return 0;
}
#include "cpp/STFT.h"
#include "WAV/WAV.h"

int main(){
  
  const int ch = 6;
  const int rate = 16000;
  const int frame = 512;
  const int shift = 128;  
  int length;

  WAV input;
  WAV output(ch,rate);
  STFT process(ch,frame,shift);

  input.OpenFile("input.wav");
  output.NewFile("output.wav");

  short buf_in[ch*shift];
  double **data;
  short buf_out[ch*shift];

  data = new double*[ch];
  for(int i=0;i<ch;i++){
    data[i] = new double[frame+2];
    memset(data[i],0,sizeof(double)*(frame+2));
  }

  while(!input.IsEOF()){
    length = input.ReadUnit(buf_in,shift*ch);
    process.stft(buf_in,length,data);
    process.istft(data,buf_out);
    output.Append(buf_out,shift*ch);
  }

  for(int i=0;i<ch;i++)
    delete[] data[i];
  delete[] data;

  return 0;
}

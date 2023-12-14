#include "cpp/STFT.h"
#include "WAV/WAV.h"

int main(){
  
  const int ch = 1;
  const int rate = 16000;
  const int frame = 320;
  const int shift = 80;  
  int length;




  WAV input;
  WAV output(ch,rate);

  FILE* f_stft;

  STFT process(ch,frame,shift);

  input.OpenFile("../data/male_1.wav");
  output.NewFile("../data/male_1_output.wav");
  f_stft = fopen("../data/male_1_stft.bin", "wb");


  //short buf_in[ch*shift];
  short buf_in[512];
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
    fwrite(data[0],sizeof(double),frame+2,f_stft);

    process.istft(data,buf_out);
    output.Append(buf_out,shift*ch);
  }

  fclose(f_stft);

  for(int i=0;i<ch;i++)
    delete[] data[i];
  delete[] data;

  return 0;
}

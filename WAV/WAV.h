#ifndef _H_WAV_
#define _H_WAV_

#include <stdlib.h>
#include <stdio.h>
//for types
#include <stdint.h>

#include <cstdlib>
#include <string.h>
#include <vector>

// Data Type is fixed as short

class WAV {
  /*
   * http://www-mmsp.ece.mcgill.ca/Documents/AudioFormats/WAVE/WAVE.html
   * */
private:
  bool non_pcm;
  // 4bytes fixed size header infomation -> uint32_t
  char riff_id[4];    // riff string
  uint32_t riff_size; // overall size of file in bytes
  char wave_id[4];    // wave string
  char fmt_id[4];     // fmt string with trailing null char

  //          | pcm  | non-pcm
  // fmt_size |  16  |     18
  //
  uint32_t fmt_size;  // format chunk size 16,18,40
  short
      fmt_type; // format type 1-PCM 3-IEEE float 6- 8bit A law, 7- 8bit ,u law
  unsigned short channels;       // no of channel
  uint32_t sample_rate; // SampleRate(blocks per second)
  uint32_t byte_rate;   // ByteRate = SampleRate * NumChannels * BitsPerSample/8
  short block_align;    // NumChannels * BitsPerSample/8
  short bit_per_sample; // bits per sample, 8 - 8bits, 16-16bits etc

  /*   (if non -pcm )*/
  uint32_t cb_size;     //size of the extension
  char fact_id[4];
  uint32_t fact_size;
  uint32_t dwSampleLength; 

  char data_id[4];      // DATA string or FLLR string
  uint32_t data_size;   // NumSamples * NumChannels * BitsPerSample/8 - size of
                        // the nex chunk that will be read
  FILE *file;
  bool IsOpen;
  const char *file_name;
  // For Input usage only
  bool use_buf;
  int frame_size;
  int shift_size;

  int size_unit;

  void* buf; 
//  short* buf; 

public:
  WAV();
  WAV(short _ch, uint32_t _rate);
  WAV(short _ch, uint32_t _rate, int frame_size, int shift_size);
  ~WAV();
  int NewFile(const char *_file_name);
  int OpenFile(const char *_file_name);
  int Append(short *app_data, unsigned int app_size);
  int Append(float*app_data, unsigned int app_size);
  void WriteHeader();
  void Finish();

  void ReadHeader();

  /* There might be compile error for ReadUnit() in Visual Studio.
   * in this case, try to update your VS to most recent version. */
  size_t ReadUnit(short*dest,int unit);
  size_t ReadUnit(float*dest,int unit);
  int IsEOF() const;

  void Print() const;
  void Rewind();

  int Convert2ShiftedArray(double **raw);
  int Convert2ShiftedArray(double *raw);

  // Split 2 channel Wav into two 1 channel wav files.
  void SplitBy2(const char* f1,const char* f2);
  void SetSizes(int frame,int shift);

  int GetChannels();
  bool GetIsOpen();
  uint32_t GetSize(); 
  uint32_t GetSizeUnit(); 
  uint32_t GetSampleRate();
  short GetFmtType();
  void UseBuf(int frame_size,int shift_size);
  bool checkValidHeader();
  void* GetBuf();

  /*Split Wav file into each channel */
  void Split(char* );
};



#endif

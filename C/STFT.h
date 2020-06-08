#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "ooura_fft.h"

struct STFT{

  int channels;
  int frame;
  int shift;

  double**raw,**data;

  // hann window
  double*hann;

  //fft
  double *a,*w;
  int *ip;

}

void STFT_init(struct STFT*stft,int ch_,int frame_,int shift_){

  int channels = ch_;
  int frame = frame_;
  int shift = shift_;
  int ol= frame - shift;

  double**raw,**data;
  double*hann;

  int i,j;

  stft->channels = channels;
  stft->frame= frame_;
  stft->shift= shift_;

  hann = stft->hann;
  raw = stft->raw;
  data = stft->data;

  raw = (double**)malloc(sizeof(double*)*ch);
  for(i=0;i<ch;i++){
     raw[i] = (double*)malloc(sizeof(double)*frame);
     memset(raw[i],0,sizeof(double)*frame);
  }
  data = (double**)malloc(sizeof(double*)*ch);
  for(i=0;i<ch;i++){
     data[i] = (double*)malloc(sizeof(double)*(frame+2));
     memset(data[i],0,sizeof(double)*(frame+2));
  }

  hann = (double*)malloc(sizeof(double)*frame);
  switch((int)(frame/shift)){
    case 4:
        hann[0] = 0.0;
        for (i = 1; i < frame; ++i)
            hann[i] = 0.5 * (1.0 - cos(2.0 * M_PI * (double)i / (double)frame));

        tmp = sqrt((double)2 / 3);
        for (i = 1; i < frame; i++)
            hann[i] *= tmp;

        break;

    case 2:
        for (i = 0; i < frame; i++) {
            hann[i] = sin(M_PI * (i + 0.5) / frame);
        }

        break;
    default:
        printf("ERROR::frame/shift(%d) is Not Supported\n",
               frame / shift);
        exit(-1);
  }

  stft->a = (double*)malloc(sizeof(double)*frame);
  stft->w = (double*)malloc(sizeof(double)*frame);
  stft->ip = (int*)malloc(sizeof(int)*(int)(sqrt(frame/2)+1));

}

void STFT_free(struct STFT *stft){
  int i;

  free(stft->a);
  free(stft->w);
  free(stft->ip);

  for(i=0;i<stft->channels){
    free(stft->raw[i]);
  }
  free(stft->raw);
  free(stft->hann);
}

/* buf_in : input of STFT, wav file or record input stream format
 *          ex) for 4 channels input
 *              1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 ...
 * length : length of input
 * data   : FFT output, [channels][frame]
 * */
void STFT(struct STFT *stft,short * buf_in,int length,double**data){
  int i,j;
  int channels,frame,shift,ol;
  double **raw,*hann;

  /*** INIT ***/
  channels = stft->channels;
  frame = stft->frame;
  shift = stft->shift;
  ol = frame - shift;

  raw = stft->raw;
  hann = stft->hann;


  /*** SHIFT ***/
  for (j = 0; j < channels; j++) {
    for (i = 0; i < ol; i++) {
      raw[j][i] = raw[j][i + shift];
    }
  }

  /*** COPY ***/
  //// EOF
  if(length!=shift*channels){

    length = length/channels;
    for (i = 0; i < length; i++) {
      for (j = 0; j < channels; j++)
        raw[j][i + ol]
          =  (double)(buf_in[i * channels+ j]);
    }
    for (i = length; i < shift; i++) {
      for (j = 0; j < channels; j++)
        raw[j][i + ol] = 0;
    }
    //// continue
  }else{
    for (i = 0; i < shift; i++) {
      for (j = 0; j < channels; j++){
        raw[j][i + ol] 
          = (double)(buf_in[i * channels+ j]);
      }
    }
  }

  //shifted raw -> data
  for(i=0;i<channels;i++)
    memcpy(data[i],raw[i],sizeof(double)*frame);

  /*** Windowing ***/
  for(i=0;i<channels;i++){
    for(j=0;j<frame;j++)
      data[i][j] *= hann[j];
  }

  /*** FFT ***/
  for (j = 0; j < channels; j++) {
    ip[0] = 0;
    for (i = 0; i < frame; i++)
      a[i] = data[j][i];

    ooura_rdft(frame, 1, a, ip, w);

    for (i = 0; i < frame; i += 2) {
      data[j][i] = a[i];
      data[j][i + 1] = -a[i + 1];
    }
    data[j][1] = 0;
    data[j][frame] = a[1];
    data[j][frame+ 1] = 0;
  }

}

void ISTFT(struct STFT *stft,double **data, short * buf_out){

  /*** IFFT ***/
    for (j = 0; j < channels; j++) {
        double *t;
        t = data[j];
        ip[j][0] = 0;
        for (int i = 0; i < frame_size; i += 2) {
            a[j][i] = t[i];
            a[j][i + 1] = -t[i + 1];
        }
        a[j][1] = t[frame_size];

        rdft(frame_size, -1, a[j], ip[j], w[j]);
        for (int i = 0; i < frame_size; i++) {
            a[j][i] *= 2.0;
            a[j][i] /= frame_size;
        }

        for (int i = 0; i < frame_size; i++) {
            t[i] = a[j][i];
        }
    }

  /*** iFFT ***/
  fft_out->iFFT(pool);

  /*** Window ***/
  hw->Process(pool, channels);

  /*** Output ***/
  if (do_mvdr)
    out->Append(ap->Overlap(pool), shift_size * 1);
  else if(do_gcp)
    out->Append(ap->OverlapSingle(pool), shift_size);
  else
    out->Append(ap->Overlap(pool), shift_size * channels);
  }

}

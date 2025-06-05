# STFT

STFT(Short Time Fourier Transform), ISTFT(Inverse - Short Time Fourier Transform) for wav,mic input.  
And other pre-processings. 

+ Provides 25%,50% overlap STFT process.   
+ DFT, iDFT for the case that n_fft(frame_size) is not power of 2.  
+ Mel-filterbank for both HTK-style(pytorch default) and slaney-style(librosa-default)  
+ LPF for perceptual quality. 

## TODO
- ERB filterbank & normalization  
- Vorbis window  & spectrogram normalization  
- C++ wrapper for C  

## NOTE

```git clone --recursive https://github.com/kooBH/STFT.git```

To build test code you need to clone ```--recursive``` to use submodule  ```koobh/WAV.git```
If you already cloned, then,
```
git submoudle init
git submodule update
```
to use submodule.

## About FFT  
I'm currently using FFT of [Ooura](http://www.kurims.kyoto-u.ac.jp/~ooura/fft.html).  
Since, it is fastest FFT in a single header file.   
But, sometimes (usually  not), there are errors between MATLAB FFT output and Ooura FFT output.     
If you need to perfectly same output as MATLAB, you have to use other FFT library.  

## CMake

```CMAKE
...
set(SRC)
set(INCL)

# it will append sources and paths to SRC and INCL
include(${PATH_TO_STFT}/STFT/cpp/STFT_cpp.cmake )  

...

add_executable(${PROJECT_NAME} ${SRC})
target_include_directories(${PROJECT_NAME}    PRIVATE    ${INCL}  )

```

## EXAMPLE

+ Directly from buffer  

```cpp

#include "cpp/STFT.h"

...

// frame/shift must be 4 or 2
STFT process(channels,frame_size,shift_size);
WAV input;
WAV output(channels,sample_rate);
input.OpenFile('input.wav');

short* buf_in; //[channels * shift_size]
double **data; //[channels][frame_size]
short* buf_out; //[channels * shift_size]

while(!input.IsEOF()){
    length = input.ReadUnit(buf_in,shift*channels);
    process.stft(buf_in,length,data);
    process.istft(data,buf_out);
    output.Append(buf_out,shift*channels);
  }

```

+ From 2D array

```cpp

#include "cpp/STFT.h"

...

// frame/shift must be 4 or 2
STFT process(channels,frame_size,shift_size);
WAV input;
WAV output(channels,sample_rate);
input.OpenFile('input.wav');

double **raw; //[channels][shift_size]
double **data; //[channels][frame_size]
short* buf_out; //[channels * shift_size]

while(!input.IsEOF()){
    input.Convert2Array(raw);
    process.stft(raw,data);
    process.istft(data,buf_out);
    output.Append(buf_out,shift*channels);
  }

```

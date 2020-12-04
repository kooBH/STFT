# STFT

STFT(Short Time Fourier Transform), ISTFT(Inverse - Short Time Fourier Transform) for wav,mic input  

provide 25%,50% overlap STFT process.  

## NOTE

```git clone --recursive https://github.com/kooBH/STFT.git```

Need to clone ```--recursive``` to use submodule for test ```koobh/WAV.git```
If you alread cloned then,
```
git submoudle init
git submodule update
```
to use submodule

---

## EXAMPLE

```cpp

#include "cpp/STFT.h"

...

// frame/shift must be 4 or 2
STFT process(channels,frame,shift);

while(!input.IsEOF()){
    length = input.ReadUnit(buf_in,shift*channels);
    process.stft(buf_in,length,data);
    process.istft(data,buf_out);
    output.Append(buf_out,shift*channels);
  }

```

## STATUS

C[WIP] : only STFT is implemented, iSTFT is not implemented yet 
cpp    : verified (same output as MATLAB routine)

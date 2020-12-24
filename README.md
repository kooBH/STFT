# STFT

STFT(Short Time Fourier Transform), ISTFT(Inverse - Short Time Fourier Transform) for wav,mic input  

provides 25%,50% overlap STFT process.  

## NOTE

```git clone --recursive https://github.com/kooBH/STFT.git```

To build test code you need to clone ```--recursive``` to use submodule  ```koobh/WAV.git```
If you already cloned, then,
```
git submoudle init
git submodule update
```
to use submodule.

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

C       : only STFT is implemented, iSTFT is not implemented yet   
C++    : verified (same output as MATLAB routine)  

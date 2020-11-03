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
STFT process(ch,frame,shift);

while(!input.IsEOF()){
    length = input.ReadUnit(buf_in,shift*ch);
    process.stft(buf_in,length,data);
    process.istft(data,buf_out);
    output.Append(buf_out,shift*ch);
  }

```

## STATUS

C[WIP] : only STFT is implemented, iSTFT is not implemented  
cpp    : verified

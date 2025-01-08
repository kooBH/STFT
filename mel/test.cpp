#include <crtdbg.h>

#include "mel.h"
#include "STFT.h"

int main() {
  _CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF);
  //_CrtSetBreakAlloc(18);

  int f_min = 20;
  int f_max = 7600;

  {
    mel m(16000, 512, 80, f_min, f_max, true);
    m.Export("cpp_mel_htk.txt");
  }
  {
    mel m(16000, 512, 80, f_min, f_max, false);
    m.Export("cpp_mel_slaney.txt");
  }

  return 0;
}
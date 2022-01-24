#include <crtdbg.h>

#include "mel.h"
#include "STFT.h"

int main() {
  _CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF);
  //_CrtSetBreakAlloc(18);

  mel m(16000, 512, 40);

  return 0;
}
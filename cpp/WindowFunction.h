#ifndef _H_WINDOW_FUNCTION_
#define _H_WINDOW_FUNCTION_

#include <cmath>
#include <cstdio>

/* Implemented
- hamming
- hann(inherits hamming)
- sine
*/

class WindowFunction {
protected:
  const double pi = 3.141592653589793;

  int frame_size, shift_size;
  bool norm;
  bool sym;

  // if sym==false, window is analysis window
  double* window;
  // if sym==true, not used
  double* synthesis_window;

  void normalize();

public:
  inline WindowFunction(int frame_size, int shift_size, bool norm = true, bool sym = true);
  inline ~WindowFunction();
  // 2D
  inline void Process(double** buf, int channels);
  // 1D - multi channel
  inline void Process(double* buf, int channels);
  // 1D - single channel 
  inline void Process(double* buf);
  // 2D
  inline void WindowWithScaling(double** buf, int channels);
  // 1D - multi channel
  inline void WindowWithScaling(double* buf, int channels);
  // 1D - single channel 
  inline void WindowWithScaling(double* buf);
};

class HammingWindow : public WindowFunction {
public:
  HammingWindow(int frame_size_, int shift_size_ = -1, double alpha_ = 0.54, double beta_ = 0.46, bool norm_ = true) :WindowFunction(frame_size_,shift_size_, norm_,true) {
    int i;
    double tmp = 0;
    // default 25% overlap
    if (shift_size == -1) {
      shift_size = frame_size / 4;
    }
    for (i = 0; i < frame_size; i++)
      window[i] = alpha_ - beta_ * cos(2.0 * pi * (i / (double)frame_size));
    if (norm)normalize();
  }
};

class HannWindow : public HammingWindow {

public:
  HannWindow(int frame_size, int shift_size = -1, bool norm = true) :HammingWindow(frame_size, shift_size, 0.5, 0.5, norm) {
  }
};

class SineWindow : public WindowFunction {
private:

public:
  SineWindow(int frame_size, int shift_size = -1, bool norm = true) : WindowFunction(frame_size, norm) {
    int i;
    // default 50% overlap
    if (shift_size == -1) {
      shift_size = frame_size / 2;
    }

    for (i = 0; i < frame_size; i++)
      window[i] = sin(pi * ((double)i + 0.5) / (double)frame_size);
    if(norm)normalize();

  }
  ~SineWindow() {
  }
};


inline WindowFunction::WindowFunction(int frame_size_, int shift_size_, bool norm_, bool sym_) {
  int i;
  double tmp = 0;
  frame_size = frame_size_;
  shift_size = shift_size_;
  norm = norm_;
  sym = sym_;

  window = new double[frame_size];
  memset(window, 0, sizeof(double) * frame_size);

  if (!sym) {
    synthesis_window = new double[frame_size];
    memset(synthesis_window, 0, sizeof(double) * frame_size);
  }

}

inline WindowFunction::~WindowFunction() {
  delete[] window;
  if (!sym)
    delete[] synthesis_window;
}

inline void WindowFunction::normalize() {
  double tmp = 0.0;
  for (int i = 0; i < frame_size; i++)
    tmp += window[i] * window[i];
  tmp /= shift_size;
  tmp = std::sqrt(tmp);
  for (int i = 0; i < frame_size; i++)
    window[i] /= tmp;

}


inline void WindowFunction::Process(double** buffer,
  int n_channels) {
  int i, j;
  for (i = 0; i < n_channels; i++) {
#pragma omp parallel for
    for (j = 0; j <frame_size; j++) {
      buffer[i][j] *= window[j];
    }
  }
}
inline void WindowFunction::Process(double* buffer,
  int cwindowels) {
  int i, j;
  for (i = 0; i < cwindowels; i++) {
#pragma omp parallel for
    for (j = 0; j < frame_size; j++) {
      buffer[i * (frame_size + 2) + j] *= window[j];
    }
  }
}

inline void WindowFunction::Process(double* buffer) {
  int j;
  for (j = 0; j < frame_size; j++) {
    buffer[j] *= window[j];
  }
}

inline void WindowFunction::WindowWithScaling(double** buffer,
  int n_channels) {
  int i, j;
  for (i = 0; i < n_channels; i++) {
#pragma omp parallel for
    for (j = 0; j < frame_size; j++) {
      buffer[i][j] *= window[j];
      buffer[i][j] /= 32767.0;
    }
  }
}
inline void WindowFunction::WindowWithScaling(double* buffer) {
  int j;
  for (j = 0; j < frame_size; j++) {
    buffer[j] *= window[j];
    buffer[j] /= 32767.0;
  }
}
inline void WindowFunction::WindowWithScaling(double* buffer,
  int cwindowels) {
  int i, j;
  for (i = 0; i < cwindowels; i++) {
#pragma omp parallel for
    for (j = 0; j < frame_size; j++) {
      buffer[i * (frame_size + 2) + j] *= window[j];
      buffer[i * (frame_size + 2) + j] /= 32767.0;
    }
  }
}

#endif

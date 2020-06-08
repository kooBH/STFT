#ifndef _H_HANN_WINDOW_
#define _H_HANN_WINDOW_

#include <cmath>
#include <cstdio>

#define M_PI 3.14159265358979323846

class HannWindow {
private:
    double *hann;
    int shift_size;
    int frame_size;


public:
    inline HannWindow(int _frame_size, int _shift_size);
    inline ~HannWindow();
    // 2D
    inline void Process(double ** buf, int channels);
    // 1D - multi channel
    inline void Process(double * buf, int channels);
    // 1D - single channel 
    inline void Process(double * buf);
    // 2D
    inline void WindowWithScaling(double ** buf, int channels);
    // 1D - multi channel
    inline void WindowWithScaling(double * buf, int channels);
    // 1D - single channel 
    inline void WindowWithScaling(double * buf);
};

inline HannWindow::HannWindow(int _frame_size, int _shift_size) {
    int i;
    double temp;

    shift_size = _shift_size;
    frame_size = _frame_size;

    // From WPE_Online
    /*
     *
     Nfft = BufferSize * 4;
     Nfreq = Nfft / 2 + 1;
     Nwin = BufferSize * 4;
     *
     */
    hann = new double[frame_size];

    switch (frame_size / shift_size) {
    case 4:
        hann[0] = 0.0;
        for (i = 1; i < frame_size; ++i)
            hann[i] = 0.5 * (1.0 - cos(2.0 * M_PI * (double)i / (double)frame_size));

        temp = sqrt((double)2 / 3);
        for (i = 1; i < frame_size; i++)
            hann[i] *= temp;

        break;

    case 2:
        for (i = 0; i < frame_size; i++) {
            hann[i] = sin(M_PI * (i + 0.5) / frame_size);
        }

        break;
    default:
        printf("ERROR::frame_size/shift_size(%d) is Not Supported\n",
               frame_size / shift_size);
        exit(-1);
    }
#ifndef NDEBUG
/*
  printf("INFO::HannWindow\n");
  for(i=0;i<frame_size;i++)
  printf("%lf ",hann[i]);
  printf("\n");
*/
#endif
}

inline HannWindow::~HannWindow() { delete[] hann; }

inline void HannWindow::Process(double **buffer,
                                   int channels) {
    int i, j;
    for (i = 0; i < channels; i++) {
#pragma omp parallel for
        for (j = 0; j < frame_size; j++) {
            buffer[i][j] *= hann[j];
        }
    }
}
inline void HannWindow::Process(double *buffer,
                                   int channels) {
    int i, j;
    for (i = 0; i < channels; i++) {
#pragma omp parallel for
        for (j = 0; j < frame_size; j++) {
            buffer[i*(frame_size+2) + j] *= hann[j];
        }
    }
}

inline void HannWindow::Process(double *buffer){
    int j;
        for (j = 0; j < frame_size; j++) {
           buffer[j] *= hann[j];
        }
}

inline void HannWindow::WindowWithScaling(double **buffer,
                                   int channels) {
    int i, j;
    for (i = 0; i < channels; i++) {
#pragma omp parallel for
        for (j = 0; j < frame_size; j++) {
            buffer[i][j] *= hann[j];
            buffer[i][j] /= 32767.0;
        }
    }
}
inline void HannWindow::WindowWithScaling(double *buffer){
    int j;
        for (j = 0; j < frame_size; j++) {
            buffer[j] *= hann[j];
            buffer[j] /= 32767.0;
        }
}
inline void HannWindow::WindowWithScaling(double *buffer,
                                   int channels) {
    int i, j;
    for (i = 0; i < channels; i++) {
#pragma omp parallel for
        for (j = 0; j < frame_size; j++) {
            buffer[ i*(frame_size+2) + j] *= hann[j];
            buffer[ i*(frame_size+2) + j] /= 32767.0;
        }
    }
}

#endif

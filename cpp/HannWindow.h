#ifndef _H_HANN_WINDOW_
#define _H_HANN_WINDOW_

#include <cmath>
#include <cstdio>

class HannWindow {
    int cnt = 0;
private:
    // c standard defind 
    //#define M_PI (3.14159265358979323846)::

    // MATLAB 'pi'
    const double MATLAB_pi = 3.141592653589793;
    double* hann;
    int shift_size;
    int frame_size;


public:
    HannWindow(int _frame_size, int _shift_size);
    ~HannWindow();
    // 2D
    void Process(double** buf, int channels);
    // 1D - multi channel
    void Process(double* buf, int channels);
    // 1D - single channel 
    void Process(double* buf);
    // 2D
    void WindowWithScaling(double** buf, int channels);
    // 1D - multi channel
    void WindowWithScaling(double* buf, int channels);
    // 1D - single channel 
    void WindowWithScaling(double* buf);
};

#endif
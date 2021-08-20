#ifndef _H_AFTER_PROCESSOR_
#define _H_AFTER_PROCESSOR_

#include <cstdint>
#include <cstdlib>
#include <cstring>

class PostProcessor {
private:
    uint32_t frame_size;
    uint32_t shift_size;
    uint32_t channels;
    // frame_size / shift_size;
    uint32_t num_block;
    short *output;
    double **buf;
    uint32_t buf_offset;

public:
    inline PostProcessor(uint32_t _frame_size,
                          uint32_t _shift_size,
                          uint32_t _channels);
    inline ~PostProcessor();

    inline short *Overlap(double **in);
    inline short *Overlap(double *in);
    inline short *Array2WavForm(double **in);

    inline short *Frame2Wav(double *in);
};

inline PostProcessor::PostProcessor(uint32_t _frame_size,
                                      uint32_t _shift_size,
                                      uint32_t _channels) {
    int i;
    channels = _channels;
    frame_size = _frame_size;
    shift_size = _shift_size;

    buf_offset = 0;

    num_block = frame_size / shift_size;

    output = new short[shift_size * channels];
    buf = new double *[channels];
    for (i = 0; i < static_cast<int>(channels); i++)
        // buf[i] = new double[ num_block * frame_size ];
        buf[i] = new double[frame_size];
    for (i = 0; i < static_cast<int>(channels); i++)
        // memset(buf[i],0,  num_block * frame_size*sizeof(double)  );
        memset(buf[i], 0, frame_size * sizeof(double));
}

inline PostProcessor::~PostProcessor() {
    for (int i = 0; i < static_cast<int>(channels); i++)
        delete[] buf[i];
    delete[] buf;
    delete[] output;
}

/*
 * buffer         output
 * 0  0  0  a1    0
 * 0  0  a2 b1    0
 * 0  a3 b2 c1    0
 * a4 b3 c2 d1    a4 + a3 + a2 + a1
 * b4 c3 d2 e1    b4 + b3 + b2 + b1
 *
 * */
inline short *PostProcessor::Overlap(double **in) {
    int i, j;
    for (j = 0; j < static_cast<int>(channels); j++) {

        // Shift
        for (i = 0; i < static_cast<int>(frame_size - shift_size); i++)
            buf[j][i] = buf[j][i + shift_size];

        // Emptying Last Block
        memset(buf[j] + shift_size * (num_block - 1), 0,
               sizeof(double) * shift_size);

        // Sum
        for (i = 0; i < static_cast<int>(frame_size); i++) {
            buf[j][i] += in[j][i];
        }
    }
    // Distribution for Wav format
    for (i = 0; i < static_cast<int>(shift_size); i++) {
#pragma ivdep
        for (j = 0; j < static_cast<int>(channels); j++) {
            //output[i * channels + j] = (short)(buf[j][i] * 32767);
            output[i * channels + j] = static_cast<short>(buf[j][i]);
#ifndef NDEBUG
// printf("Overlap::buf[%d][%d]= %lf\n",j,i,buf[j][i]);
#endif
        }
    }
    return output;
}


inline short *PostProcessor::Overlap(double *in) {
  int i;
    // Shift
    for (i = 0; i < static_cast<int>(frame_size - shift_size); i++)
      buf[0][i] = buf[0][i + shift_size];

    // Emptying Last Block
    memset(buf[0] + shift_size * (num_block - 1), 0,
        sizeof(double) * shift_size);

    // Sum
    for (i = 0; i < static_cast<int>(frame_size); i++) {
      buf[0][i] += in[i];
    }

  // Distribution for Wav format
  for (i = 0; i < static_cast<int>(shift_size); i++) {
      output[i] = static_cast<short>(buf[0][i]);
  }

  return output;
}

inline short *PostProcessor::Array2WavForm(double **in) {
    int i, j;
    for (i = 0; i < static_cast<int>(shift_size); i++) {
        for (j = 0; j < static_cast<int>(channels); j++) {
            output[i * channels + j] = (short)in[j][i];
        }
    }
    return output;
}

inline short *PostProcessor::Frame2Wav(double *in) {
    int i, j;
    for (i = 0; i < static_cast<int>(frame_size); i++) {
        for (j = 0; j < static_cast<int>(channels); j++) {
            output[i * channels + j] = (short)in[j*frame_size + i];
        }
    }
    return output;
}

#endif

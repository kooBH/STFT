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
    PostProcessor(uint32_t _frame_size,
                          uint32_t _shift_size,
                          uint32_t _channels);
    ~PostProcessor();

    short *Overlap(double **in);
    short *Overlap(double *in);
    short *Array2WavForm(double **in);

    short *Frame2Wav(double *in);
    short *Get_output(); 
    double** Get_buf(); 
};

#endif
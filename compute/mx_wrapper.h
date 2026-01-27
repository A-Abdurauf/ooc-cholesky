#ifndef MX_WRAPPER_H
#define MX_WRAPPER_H

extern "C" {

// Launch quantization kernel for FP32 → MX
void mx_quantize_fp32_to_mx(
    const float* d_input,   // device ptr
    float* d_output,        // device ptr
    const float* d_maxvals, // device ptr
    long total_size,
    int axis_size,
    int post_axis_size,
    int scale_bits,
    int elem_ebits,
    int elem_mbits,
    float elem_max_norm,
    bool flush_fp32_subnorms,
    int rounding_mode
);

}

#endif

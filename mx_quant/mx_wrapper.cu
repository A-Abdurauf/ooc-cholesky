#include "mx_wrapper.h"
#include "mx.cuh"
#include <cuda_runtime.h>

extern "C" {

void mx_quantize_fp32_to_mx(
    const float* d_input,
    float* d_output,
    const float* d_maxvals,
    long total_size,
    int axis_size,
    int post_axis_size,
    int scale_bits,
    int elem_ebits,
    int elem_mbits,
    float elem_max_norm,
    bool flush_fp32_subnorms,
    int rounding_mode
) {
    int threads = 256;
    int blocks  = (total_size + threads - 1) / threads;

    quantize_mx_cuda_kernel<float><<<blocks, threads>>>(
        d_input,
        scale_bits,
        elem_ebits,
        elem_mbits,
        elem_max_norm,
        d_maxvals,
        total_size,
        axis_size,
        post_axis_size,
        flush_fp32_subnorms,
        (RoundingMode)rounding_mode,
        d_output
    );
}

} // extern "C"

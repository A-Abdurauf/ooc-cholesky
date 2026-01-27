#ifndef MX_SHARED_EXP_CUH
#define MX_SHARED_EXP_CUH

#include <cuda.h>
#include <cuda_runtime.h>
#include "common.cuh"

// max norm = (2^(elem_ebits-1) - 1)
// Example: E4M3 → 2^(4-1) - 1 = 7
__device__ __forceinline__
float mx_get_shared_scale(int shared_exp,
                          int scale_bits,
                          float elem_max_norm)
{
    // Let bias = 127 for FP32
    int fp32_bias = 127;
    int unbiased_exp = shared_exp - fp32_bias + 1;

    float max_val = ldexpf(1.0f, unbiased_exp);  // 2^exp
    float scale = max_val / elem_max_norm;
    return scale;
}

#endif

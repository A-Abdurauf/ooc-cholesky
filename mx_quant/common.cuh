#ifndef MX_COMMON_CUH
#define MX_COMMON_CUH

#include <cuda.h>
#include <cuda_runtime.h>
#include <stdint.h>
#include <math.h>

//--------------------------------------------
// Rounding mode enum
//--------------------------------------------
enum RoundingMode {
    ROUND_TO_NEAREST_EVEN = 0,
    ROUND_TOWARD_ZERO     = 1,
    ROUND_AWAY_FROM_ZERO  = 2
};

//--------------------------------------------
// Fast float-as-int reinterpret
//--------------------------------------------
__device__ __forceinline__ uint32_t as_uint(float x) {
    return reinterpret_cast<uint32_t&>(x);
}

__device__ __forceinline__ float as_float(uint32_t x) {
    return reinterpret_cast<float&>(x);
}

//--------------------------------------------
// Extract biased exponent of FP32
//--------------------------------------------
__device__ __forceinline__
int get_biased_exponent(float x)
{
    uint32_t bits = as_uint(x);
    int exp = (bits >> 23) & 0xFF;
    return exp;   // already biased
}

#endif

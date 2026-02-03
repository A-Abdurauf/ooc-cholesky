#ifndef MX_QUANTIZE_CUH
#define MX_QUANTIZE_CUH

#include <cuda.h>
#include <cuda_fp16.h>
#include "common.cuh"

//------------------------------------------------------
// Quantize a FP32 value into pseudo-MX (FP32 emulated)
//------------------------------------------------------
__device__ __forceinline__
float quantize_mx_elem(float x,
                       float scale,
                       bool flush_tile,
                       int ebits,
                       int mbits,
                       float elem_max_norm,
                       RoundingMode rm)
{
    // Special-case MX_FP16 to match Eigen path:
    // scale -> FP16 round-to-nearest-even -> rescale (no explicit clamp).
    if (ebits == 5 && mbits == 10) {
        float scaled = x / scale;
        __half h = __float2half_rn(scaled);
        float q = __half2float(h);
        return q * scale;
    }

    float scaled = x / scale;

    // saturate
    if (scaled > elem_max_norm) scaled = elem_max_norm;
    if (scaled < -elem_max_norm) scaled = -elem_max_norm;

    // flush subnormals
    // Note: The denominator in ldexpf is 1 - ebits, not ebits.
    if (flush_tile && fabsf(scaled) < ldexpf(1.f, 1 - ebits))
        scaled = 0.f;

    // Round
    float q;
    switch (rm) {
        case ROUND_TO_NEAREST_EVEN:
            q = nearbyintf(scaled);
            break;
        case ROUND_TOWARD_ZERO:
            q = truncf(scaled);
            break;
        case ROUND_AWAY_FROM_ZERO:
            q = (scaled > 0 ? ceilf(scaled) : floorf(scaled));
            break;
    }

    // --- CRITICAL FIX: BIT MASKING TO EMULATE MX PRECISION ---
    // Multiply the rounded integer 'q' by the scale factor.
    float result = q * scale;
    
    // The scale factor is 2^X / 15.0f. Multiplying by this non-power-of-2 factor
    // results in a full-precision FP32 result with a 'messy' mantissa.
    // To emulate the M=3 precision loss, we must manually truncate the FP32 mantissa.
    
    // Total FP32 mantissa bits is 23.
    const int fp32_mantissa_bits = 23;
    
    // We want to keep 'mbits' (e.g., 3) bits, so we zero out the rest.
    const int bits_to_truncate = fp32_mantissa_bits - mbits; // 23 - 3 = 20 bits
    
    // Only perform truncation if mbits is less than FP32 precision (i.e., mbits < 23)
    if (bits_to_truncate > 0) {
        // 1. Get the bit representation of the result
        uint32_t bits = as_uint(result);
        
        // 2. Create a mask: ~((1 << bits_to_truncate) - 1)
        // This is a mask of 1s in the upper bits and 0s in the lower bits_to_truncate.
        // E.g., for 20 bits: ...1110000000000000000000000000000
        uint32_t mask = ~((1 << bits_to_truncate) - 1);
        
        // 3. Apply the mask to the mantissa part (lower 23 bits).
        // The mantissa starts at bit 0. We must only mask the mantissa.
        // Bits: S (31) | EXP (30:23) | MANTISSA (22:0)
        
        // Ensure the mask only applies to the mantissa area.
        // We only need to shift the mask so that the lower bits are zeroed out.
        // Since we only care about the mantissa, we apply the mask directly to the whole integer.
        // This implicitly assumes the exponent of the result is not zero (not subnormal).
        bits &= mask;
        
        // 4. Convert back to float
        result = as_float(bits);
    }

    return result;
}

#endif
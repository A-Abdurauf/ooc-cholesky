#ifndef PYT_MX_MX_CUH     // legacy guard
#define PYT_MX_MX_CUH

#include "common.cuh"
#include "shared_exp.cuh"
#include "quantize.cuh"

//-----------------------------------------------------------------------
// quantize_mx_cuda_kernel (generic)
//-----------------------------------------------------------------------
template<typename T>
__global__ void quantize_mx_cuda_kernel(
    const T* __restrict__ input,
    const int scale_bits,
    const int elem_ebits,
    const int elem_mbits,
    const float elem_max_norm,
    const T* __restrict__ max_values,
    const long total_size,
    const int axis_size,
    const int post_axis_size,
    const bool flush_fp32_subnorms,
    const RoundingMode rounding_mode,
    T* __restrict__ output
) {
    const long offset = blockDim.x * blockIdx.x + threadIdx.x;
    if (offset >= total_size) return;

    const long post_axis_i = offset % post_axis_size;
    const long pre_axis_i = offset / (post_axis_size * axis_size);

    const long m_i = pre_axis_i * post_axis_size + post_axis_i;
    // CORRECT: Just read the float value directly
    const float max_elem = max_values[m_i];
    int shared_exp = get_biased_exponent(max_elem);
    bool flush_tile = (shared_exp == 0 && flush_fp32_subnorms);

    const float scale =
        mx_get_shared_scale(shared_exp, scale_bits, elem_max_norm);

    output[offset] = quantize_mx_elem(
        input[offset], scale, flush_tile,
        elem_ebits, elem_mbits, elem_max_norm, rounding_mode);
}

//-----------------------------------------------------------------------
// quantize_mx_innermost_cuda_kernel (warp tile)
//-----------------------------------------------------------------------
template<typename T>
__global__ void quantize_mx_innermost_cuda_kernel (
    const T* __restrict__ in,
    const int scale_bits,
    const int elem_ebits,
    const int elem_mbits,
    const float elem_max_norm,
    const long total_size,
    const int tile_size,
    const bool flush_fp32_subnorms,
    const RoundingMode rounding_mode,
    T* __restrict__ out
) {
    const long offset = blockDim.x * blockIdx.x + threadIdx.x;
    if (offset >= total_size) return;

    float elem = in[offset];
    float abs_elem = fabsf(elem);

    // warp-level reduction
    for (int mask = tile_size/2; mask > 0; mask /= 2) {
        float v = __shfl_xor_sync(0xFFFFFFFF, abs_elem, mask);
        abs_elem = fmaxf(abs_elem, v);
    }

    int shared_exp = get_biased_exponent(abs_elem);
    bool flush_tile = (shared_exp == 0 && flush_fp32_subnorms);

    float scale =
        mx_get_shared_scale(shared_exp, scale_bits, elem_max_norm);

    out[offset] = quantize_mx_elem(
        elem, scale, flush_tile,
        elem_ebits, elem_mbits, elem_max_norm, rounding_mode
    );
}

//-----------------------------------------------------------------------
// quantize_mx_by_tile_cuda_kernel
//-----------------------------------------------------------------------
template<typename T>
__global__ void quantize_mx_by_tile_cuda_kernel (
    const T* __restrict__ in,
    const int scale_bits,
    const int elem_ebits,
    const int elem_mbits,
    const float elem_max_norm,
    const int total_tiles,
    const int tile_size,
    const int num_tiles,
    const int axis_size,
    const int post_axis_size,
    const bool flush_fp32_subnorms,
    const RoundingMode rounding_mode,
    T* __restrict__ out
) {
    const long offset = blockDim.x * blockIdx.x + threadIdx.x;
    if (offset >= total_tiles) return;

    const long post_axis_i = offset % post_axis_size;
    const long num_tiles_i = (offset / post_axis_size) % num_tiles;
    const long pre_axis_i = offset / (num_tiles * post_axis_size);

    int adjusted_tile_size =
        ((num_tiles_i + 1) * tile_size > axis_size)
        ? (axis_size % tile_size)
        : tile_size;

    float abs_max_value = 0.0f;
    for (int i = 0; i < adjusted_tile_size; i++) {
        long idx =
            pre_axis_i * axis_size * post_axis_size +
            (num_tiles_i * tile_size + i) * post_axis_size +
            post_axis_i;

        float v = fabsf(in[idx]);
        abs_max_value = fmaxf(abs_max_value, v);
    }

    int shared_exp = get_biased_exponent(abs_max_value);
    bool flush_tile = (shared_exp == 0 && flush_fp32_subnorms);

    float scale =
        mx_get_shared_scale(shared_exp, scale_bits, elem_max_norm);

    for (int i = 0; i < adjusted_tile_size; i++) {
        long idx =
            pre_axis_i * axis_size * post_axis_size +
            (num_tiles_i * tile_size + i) * post_axis_size +
            post_axis_i;

        out[idx] = quantize_mx_elem(
            in[idx], scale, flush_tile,
            elem_ebits, elem_mbits, elem_max_norm, rounding_mode
        );
    }
}

#endif

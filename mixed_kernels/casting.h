#pragma once

#include <cublas_v2.h>

namespace mixed_kernels {
const void *casting1D(const cudaStream_t &stream, const void *from,
                      cudaDataType fromType, void *to, cudaDataType toType,
                      size_t size);

const void *casting2DRowMajor(const cudaStream_t &stream, int m, int n,
                              const void *from, cudaDataType fromType,
                              size_t ldFrom, void *to, cudaDataType toType,
                              size_t ldTo);

// to = from + beta * to
void plus2DRowMajor(const cudaStream_t &stream, int m, int n, double beta,
                    const void *from, cudaDataType fromType, size_t ldFrom,
                    void *to, cudaDataType toType, size_t ldTo);

bool operator<(cudaDataType_t lhs, cudaDataType_t rhs);

size_t cudaDataSize(cudaDataType type);
}  // namespace mixed_kernels

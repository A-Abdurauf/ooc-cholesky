#pragma once

#include <cublas_v2.h>

CUBLASAPI cublasStatus_t CUBLASWINAPI mixedKernelsGemmExRowMajorNTN(
    cublasHandle_t handle, void* cublasLtMeta, int m, int n, int k,
    double alpha, const void* A, cudaDataType Atype, int lda, const void* B,
    cudaDataType Btype, int ldb, double beta, void* C, cudaDataType Ctype,
    int ldc, cublasComputeType_t computeType, cublasGemmAlgo_t algo,
    void* buffer, size_t bufferSize);

CUBLASAPI cublasStatus_t CUBLASWINAPI mixedKernelsSyrkExRowMajorNTN(
    cublasHandle_t handle, void* cublasLtMeta, cublasFillMode_t uplo, int n,
    int k, double alpha, const void* A, cudaDataType Atype, int lda,
    double beta, void* C, cudaDataType Ctype, int ldc,
    cublasComputeType_t computeType, cublasGemmAlgo_t algo, void* buffer,
    size_t bufferSize);

namespace mixed_kernels::cublasLtFp8RowMajorNTNMeta {
CUBLASAPI cublasStatus_t CUBLASWINAPI
create(void** instance, uint64_t m, uint64_t n, uint64_t k, int64_t lda,
       int64_t ldb, int64_t ldc, std::size_t workspaceSize, void* workspace);

CUBLASAPI cublasStatus_t CUBLASWINAPI destroy(void* instance);

static CUBLASAPI cublasStatus_t CUBLASWINAPI
matmul(void* instance, cudaStream_t stream, double alpha, const void* A,
       const void* B, double beta, void* C, cudaDataType Ctype,
       void* castingBuffer, std::size_t castingBufferSize);
}  // namespace mixed_kernels::cublasLtFp8RowMajorNTNMeta
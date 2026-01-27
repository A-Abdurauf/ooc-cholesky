#include <cublas_v2.h>

#include <cstdio>

#include "casting.h"
#include "mixed_kernels.h"

CUBLASAPI cublasStatus_t CUBLASWINAPI mixedKernelsSyrkExRowMajorNTN(
    cublasHandle_t handle, void* cublasLtMeta, cublasFillMode_t uplo, int n,
    int k, double alpha, const void* A, cudaDataType Atype, int lda,
    double beta, void* C, cudaDataType Ctype, int ldc,
    cublasComputeType_t computeType, cublasGemmAlgo_t algo, void* buffer,
    size_t bufferSize) {
  const double alphaFP64 = alpha;
  const double betaFP64 = beta;
  const float alphaFP32 = alpha;
  const float betaFP32 = beta;
  switch (uplo) {
    case CUBLAS_FILL_MODE_UPPER:
      uplo = CUBLAS_FILL_MODE_LOWER;
      break;
    case CUBLAS_FILL_MODE_LOWER:
      uplo = CUBLAS_FILL_MODE_UPPER;
      break;
    default:
      break;
  }
  if (Ctype != CUDA_R_64F) {
    return CUBLAS_STATUS_INVALID_VALUE;
  }
  if (bufferSize < sizeof(double) * n * k) {
    printf("ERROR: bufferSize too small\n");
    return CUBLAS_STATUS_INVALID_VALUE;
  }
  cudaStream_t stream;
  cublasGetStream_v2(handle, &stream);
  cublasStatus_t status;
  switch (Atype) {
    case CUDA_R_8F_E4M3:
    case CUDA_R_16F:
    case CUDA_R_16BF:
    case CUDA_R_32F:
      status = mixedKernelsGemmExRowMajorNTN(
          handle, cublasLtMeta, n, n, k, alpha, A, Atype, lda, A, Atype, lda, 0,
          buffer, CUDA_R_32F, ldc, CUBLAS_COMPUTE_32F, algo, nullptr, 0);
      mixed_kernels::plus2DRowMajor(stream, n, n, beta, buffer, CUDA_R_32F,
                                    sizeof(float) * ldc, C, Ctype,
                                    mixed_kernels::cudaDataSize(Ctype) * ldc);
      return status;
    case CUDA_R_64F:
      return cublasDsyrk_v2(handle, uplo, CUBLAS_OP_T, n, k, &alphaFP64,
                            static_cast<const double*>(A), lda, &betaFP64,
                            static_cast<double*>(C), ldc);
    default:
      return CUBLAS_STATUS_INVALID_VALUE;
  }
}

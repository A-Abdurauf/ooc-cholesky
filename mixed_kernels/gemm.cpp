#include <cublasLt.h>

#include <cstdint>
#include <cstdio>

#include "casting.h"

namespace mixed_kernels {
#define CHECK_CUBLAS(statement)                        \
  do {                                                 \
    auto status = statement;                           \
    if (status != CUBLAS_STATUS_SUCCESS) {             \
      printf("Failed at %s:%d\n", __FILE__, __LINE__); \
    }                                                  \
                                                       \
  } while (0)

namespace cublasLtFp8RowMajorNTNMeta {
namespace {
struct Impl {
  const uint64_t m, n, k;
  const int64_t lda, ldb, ldc;
  const std::size_t workspaceSize;
  void* workspace;
  cublasLtHandle_t handle{nullptr};
  const cublasOperation_t transa{CUBLAS_OP_T}, transb{CUBLAS_OP_N};
  cublasLtMatmulPreference_t preference{nullptr};
  cublasLtMatrixLayout_t Adesc{nullptr}, Bdesc{nullptr};
  cublasLtMatmulDesc_t operationDesc{nullptr};
  struct {
    cublasLtMatrixLayout_t Cdesc{nullptr}, Ddesc{nullptr};
    cublasLtMatmulHeuristicResult_t heuristicResult{};
  } fp8fp8fp8{}, fp8fp8fp16{}, fp8fp8fp32{};
};
}  // namespace

CUBLASAPI cublasStatus_t CUBLASWINAPI
create(void** instance, uint64_t m, uint64_t n, uint64_t k, int64_t lda,
       int64_t ldb, int64_t ldc, std::size_t workspaceSize, void* workspace);

CUBLASAPI cublasStatus_t CUBLASWINAPI destroy(void* instance);

static CUBLASAPI cublasStatus_t CUBLASWINAPI
matmul(void* instance, cudaStream_t stream, double alpha, const void* A,
       const void* B, double beta, void* C, cudaDataType Ctype,
       void* castingBuffer, std::size_t castingBufferSize);
}  // namespace cublasLtFp8RowMajorNTNMeta

cublasStatus_t cublasLtFp8RowMajorNTNMeta::create(
    void** instance, uint64_t m, uint64_t n, uint64_t k, int64_t lda,
    int64_t ldb, int64_t ldc, std::size_t workspaceSize, void* workspace) {
  cudaDeviceProp deviceProp;
  cudaGetDeviceProperties_v2(&deviceProp, 0);
  if (deviceProp.major < 9) {
    return CUBLAS_STATUS_SUCCESS;
  }
  auto impl = new Impl{m, n, k, lda, ldb, ldc, workspaceSize, workspace};
  CHECK_CUBLAS(cublasLtCreate(&impl->handle));
  CHECK_CUBLAS(cublasLtMatmulPreferenceCreate(&impl->preference));
  CHECK_CUBLAS(cublasLtMatmulPreferenceSetAttribute(
      impl->preference, CUBLASLT_MATMUL_PREF_MAX_WORKSPACE_BYTES,
      &workspaceSize, sizeof(workspaceSize)));
  CHECK_CUBLAS(
      cublasLtMatrixLayoutCreate(&impl->Adesc, CUDA_R_8F_E4M3, k, m, lda));
  CHECK_CUBLAS(
      cublasLtMatrixLayoutCreate(&impl->Bdesc, CUDA_R_8F_E4M3, k, n, ldb));
  CHECK_CUBLAS(cublasLtMatmulDescCreate(&impl->operationDesc,
                                        CUBLAS_COMPUTE_32F, CUDA_R_32F));
  CHECK_CUBLAS(cublasLtMatmulDescSetAttribute(
      impl->operationDesc, CUBLASLT_MATMUL_DESC_TRANSA, &impl->transa,
      sizeof(impl->transa)));
  CHECK_CUBLAS(cublasLtMatmulDescSetAttribute(
      impl->operationDesc, CUBLASLT_MATMUL_DESC_TRANSB, &impl->transb,
      sizeof(impl->transb)));

  int returnedResults;
  CHECK_CUBLAS(cublasLtMatrixLayoutCreate(&impl->fp8fp8fp8.Cdesc, CUDA_R_16F, n,
                                          m, ldc));
  CHECK_CUBLAS(cublasLtMatrixLayoutCreate(&impl->fp8fp8fp8.Ddesc,
                                          CUDA_R_8F_E4M3, n, m, ldc));
  CHECK_CUBLAS(cublasLtMatmulAlgoGetHeuristic(
      impl->handle, impl->operationDesc, impl->Bdesc, impl->Adesc,
      impl->fp8fp8fp8.Cdesc, impl->fp8fp8fp8.Ddesc, impl->preference, 1,
      &impl->fp8fp8fp8.heuristicResult, &returnedResults));
  if (returnedResults == 0) {
    CHECK_CUBLAS(CUBLAS_STATUS_NOT_SUPPORTED);
  }

  CHECK_CUBLAS(cublasLtMatrixLayoutCreate(&impl->fp8fp8fp16.Cdesc, CUDA_R_16F,
                                          n, m, ldc));
  CHECK_CUBLAS(cublasLtMatrixLayoutCreate(&impl->fp8fp8fp16.Ddesc, CUDA_R_16F,
                                          n, m, ldc));
  CHECK_CUBLAS(cublasLtMatmulAlgoGetHeuristic(
      impl->handle, impl->operationDesc, impl->Bdesc, impl->Adesc,
      impl->fp8fp8fp16.Cdesc, impl->fp8fp8fp16.Ddesc, impl->preference, 1,
      &impl->fp8fp8fp16.heuristicResult, &returnedResults));
  if (returnedResults == 0) {
    CHECK_CUBLAS(CUBLAS_STATUS_NOT_SUPPORTED);
  }

  CHECK_CUBLAS(cublasLtMatrixLayoutCreate(&impl->fp8fp8fp32.Cdesc, CUDA_R_32F,
                                          n, m, ldc));
  CHECK_CUBLAS(cublasLtMatrixLayoutCreate(&impl->fp8fp8fp32.Ddesc, CUDA_R_32F,
                                          n, m, ldc));
  CHECK_CUBLAS(cublasLtMatmulAlgoGetHeuristic(
      impl->handle, impl->operationDesc, impl->Bdesc, impl->Adesc,
      impl->fp8fp8fp32.Cdesc, impl->fp8fp8fp32.Ddesc, impl->preference, 1,
      &impl->fp8fp8fp32.heuristicResult, &returnedResults));
  if (returnedResults == 0) {
    CHECK_CUBLAS(CUBLAS_STATUS_NOT_SUPPORTED);
  }
  *instance = impl;
  return CUBLAS_STATUS_SUCCESS;
}

cublasStatus_t cublasLtFp8RowMajorNTNMeta::matmul(
    void* instance, cudaStream_t stream, double alpha, const void* A,
    const void* B, double beta, void* C, cudaDataType Ctype,
    void* castingBuffer, std::size_t castingBufferSize) {
  auto impl = static_cast<Impl*>(instance);
  float alphaFP32 = alpha;
  float betaFP32 = beta;
  switch (Ctype) {
    case CUDA_R_32F:
      return cublasLtMatmul(impl->handle, impl->operationDesc, &alphaFP32, B,
                            impl->Bdesc, A, impl->Adesc, &betaFP32, C,
                            impl->fp8fp8fp32.Cdesc, C, impl->fp8fp8fp32.Ddesc,
                            &impl->fp8fp8fp32.heuristicResult.algo,
                            impl->workspace, impl->workspaceSize, stream);
    case CUDA_R_16F:
      return cublasLtMatmul(impl->handle, impl->operationDesc, &alphaFP32, B,
                            impl->Bdesc, A, impl->Adesc, &betaFP32, C,
                            impl->fp8fp8fp16.Cdesc, C, impl->fp8fp8fp16.Ddesc,
                            &impl->fp8fp8fp16.heuristicResult.algo,
                            impl->workspace, impl->workspaceSize, stream);
    case CUDA_R_8F_E4M3: {
      if (castingBufferSize < impl->m * impl->n * sizeof(nv_half)) {
        return CUBLAS_STATUS_INVALID_VALUE;
      }
      auto castingBuffer_ =
          casting2DRowMajor(stream, impl->m, impl->n, C, CUDA_R_8F_E4M3,
                            sizeof(uint8_t) * impl->ldc, castingBuffer,
                            CUDA_R_16F, sizeof(nv_half) * impl->ldc);
      return cublasLtMatmul(
          impl->handle, impl->operationDesc, &alphaFP32, B, impl->Bdesc, A,
          impl->Adesc, &betaFP32, castingBuffer_, impl->fp8fp8fp8.Cdesc, C,
          impl->fp8fp8fp8.Ddesc, &impl->fp8fp8fp8.heuristicResult.algo,
          impl->workspace, impl->workspaceSize, stream);
    }
    default:
      return CUBLAS_STATUS_NOT_SUPPORTED;
  }
}

cublasStatus_t cublasLtFp8RowMajorNTNMeta::destroy(void* instance) {
  cudaDeviceProp deviceProp;
  cudaGetDeviceProperties_v2(&deviceProp, 0);
  if (deviceProp.major < 9) {
    return CUBLAS_STATUS_SUCCESS;
  }
  auto impl = static_cast<Impl*>(instance);
  CHECK_CUBLAS(cublasLtMatrixLayoutDestroy(impl->fp8fp8fp32.Ddesc));
  CHECK_CUBLAS(cublasLtMatrixLayoutDestroy(impl->fp8fp8fp32.Cdesc));
  CHECK_CUBLAS(cublasLtMatrixLayoutDestroy(impl->fp8fp8fp16.Ddesc));
  CHECK_CUBLAS(cublasLtMatrixLayoutDestroy(impl->fp8fp8fp16.Cdesc));
  CHECK_CUBLAS(cublasLtMatrixLayoutDestroy(impl->fp8fp8fp8.Ddesc));
  CHECK_CUBLAS(cublasLtMatrixLayoutDestroy(impl->fp8fp8fp8.Cdesc));
  CHECK_CUBLAS(cublasLtMatmulDescDestroy(impl->operationDesc));
  CHECK_CUBLAS(cublasLtMatrixLayoutDestroy(impl->Adesc));
  CHECK_CUBLAS(cublasLtMatrixLayoutDestroy(impl->Bdesc));
  CHECK_CUBLAS(cublasLtMatmulPreferenceDestroy(impl->preference));
  CHECK_CUBLAS(cublasLtDestroy(impl->handle));
  return CUBLAS_STATUS_SUCCESS;
}
}  // namespace mixed_kernels

template <int size>
void* alignUp(void* ptr) {
  auto p = reinterpret_cast<std::uintptr_t>(ptr);
  std::uintptr_t aligned = (p + size - 1) & ~(size - 1);
  return reinterpret_cast<void*>(aligned);
}

CUBLASAPI cublasStatus_t CUBLASWINAPI mixedKernelsGemmExRowMajorNTN(
    cublasHandle_t handle, void* cublasLtMeta, int m, int n, int k,
    double alpha, const void* A, cudaDataType Atype, int lda, const void* B,
    cudaDataType Btype, int ldb, double beta, void* C, cudaDataType Ctype,
    int ldc, cublasComputeType_t computeType, cublasGemmAlgo_t algo,
    void* buffer, size_t bufferSize) {
  const double alphaFP64 = alpha;
  const double betaFP64 = beta;
  const float alphaFP32 = alpha;
  const float betaFP32 = beta;
  const half alphaFP16 = alpha;
  const half betaFP16 = beta;
  // easy part that we do not have to do any casting
  // If A, B and C have the same type and are not FP8, we can call cublas simply
  if (Atype == Btype && Atype == Ctype && Atype != CUDA_R_8F_E4M3 &&
      Atype != CUDA_R_8F_E5M2) {
    switch (computeType) {
      case CUBLAS_COMPUTE_16F:
      case CUBLAS_COMPUTE_16F_PEDANTIC:
        return cublasGemmEx(handle, CUBLAS_OP_T, CUBLAS_OP_N, n, m, k,
                            &alphaFP16, B, Btype, ldb, A, Atype, lda, &betaFP16,
                            C, Ctype, ldc, computeType, algo);

      case CUBLAS_COMPUTE_32F_FAST_16F:
      case CUBLAS_COMPUTE_32F_FAST_16BF:
      case CUBLAS_COMPUTE_32F_FAST_TF32:
      case CUBLAS_COMPUTE_32F:
      case CUBLAS_COMPUTE_32F_PEDANTIC:
        return cublasGemmEx(handle, CUBLAS_OP_T, CUBLAS_OP_N, n, m, k,
                            &alphaFP32, B, Btype, ldb, A, Atype, lda, &betaFP32,
                            C, Ctype, ldc, computeType, algo);

      case CUBLAS_COMPUTE_64F:
      case CUBLAS_COMPUTE_64F_PEDANTIC:
        return cublasGemmEx(handle, CUBLAS_OP_T, CUBLAS_OP_N, n, m, k,
                            &alphaFP64, B, Btype, ldb, A, Atype, lda, &betaFP64,
                            C, Ctype, ldc, computeType, algo);
      default:
        return CUBLAS_STATUS_INVALID_VALUE;
    }
  }
  // For D(FP32) = alpha A(FP16/BF16)B(FP16/BF16) + beta C(FP32)
  // Notice here we do not have to deal D is a 16-bit FP, because in this case
  // A, B and C have the same types, and will be dealt above
  if (Atype == Btype && (Atype == CUDA_R_16F || Atype == CUDA_R_16BF) &&
      (Ctype == CUDA_R_32F)) {
    switch (computeType) {
      case CUBLAS_COMPUTE_32F:
      case CUBLAS_COMPUTE_32F_PEDANTIC:
        return cublasGemmEx(handle, CUBLAS_OP_T, CUBLAS_OP_N, n, m, k,
                            &alphaFP32, B, Btype, ldb, A, Atype, lda, &betaFP32,
                            C, Ctype, ldc, computeType, algo);

      default:
        return CUBLAS_STATUS_INVALID_VALUE;
    }
  }

  // casting needed code
  cudaStream_t stream;
  cublasGetStream_v2(handle, &stream);

  if (Atype == Btype && Atype == CUDA_R_8F_E4M3 && (Ctype == CUDA_R_16F || Ctype == CUDA_R_32F)) {
    return mixed_kernels::cublasLtFp8RowMajorNTNMeta::matmul(
            cublasLtMeta, stream, alpha, A, B, beta, C, Ctype, nullptr, 0);
  }

  if (bufferSize < sizeof(double) * (m + n) * k) {
    printf("ERROR: bufferSize too small\n");
    return CUBLAS_STATUS_INVALID_VALUE;
  }

  void* bufferA = buffer;
  void* bufferB = static_cast<void*>(reinterpret_cast<std::byte*>(buffer) +
                                     (sizeof(double) * k * m));
  using mixed_kernels::operator<;
  switch (Ctype) {
    case CUDA_R_64F: {
      if (Atype == CUDA_R_8F_E4M3 && Btype == CUDA_R_8F_E4M3) {
        if (bufferSize < sizeof(float) * m * n) {
          printf("ERROR: bufferSize too small\n");
          return CUBLAS_STATUS_INVALID_VALUE;
        }
        auto status = mixed_kernels::cublasLtFp8RowMajorNTNMeta::matmul(
            cublasLtMeta, stream, alpha, A, B, 0 /* set to zero here */, buffer,
            CUDA_R_32F, nullptr, 0);
        mixed_kernels::plus2DRowMajor(stream, m, n, beta, buffer, CUDA_R_32F,
                                      sizeof(float) * ldc, C, CUDA_R_64F,
                                      sizeof(double) * ldc);
        return status;
      }
      if (Atype < CUDA_R_32F && Btype < CUDA_R_32F) {
        if (bufferSize <
            sizeof(float) * m * n + sizeof(nv_half) * k * (m + n)) {
          printf("ERROR: bufferSize too small\n");
          return CUBLAS_STATUS_INVALID_VALUE;
        }
        bufferB = alignUp<16>(static_cast<nv_half*>(bufferA) +
                              static_cast<std::size_t>(k) * m);
        auto bufferCasting = alignUp<16>(static_cast<nv_half*>(bufferB) +
                                         static_cast<std::size_t>(k) * n);
        auto bufferA_ = mixed_kernels::casting2DRowMajor(
            stream, m, k, A, Atype, mixed_kernels::cudaDataSize(Atype) * lda,
            bufferA, CUDA_R_16F, sizeof(nv_half) * lda);
        auto bufferB_ = mixed_kernels::casting2DRowMajor(
            stream, n, k, B, Btype, mixed_kernels::cudaDataSize(Btype) * ldb,
            bufferB, CUDA_R_16F, sizeof(nv_half) * lda);
        float zeroFP32 = 0;
        auto status = cublasGemmEx(
            handle, CUBLAS_OP_T, CUBLAS_OP_N, n, m, k, &alphaFP32, bufferB_,
            CUDA_R_16F, ldb, bufferA_, CUDA_R_16F, lda, &zeroFP32,
            bufferCasting, CUDA_R_32F, ldc, CUBLAS_COMPUTE_32F, algo);
        mixed_kernels::plus2DRowMajor(stream, m, n, beta, bufferCasting,
                                      CUDA_R_32F, sizeof(float) * ldc, C,
                                      CUDA_R_64F, sizeof(double) * ldc);
        return status;
      }
      auto bufferA_ = mixed_kernels::casting2DRowMajor(
          stream, m, k, A, Atype, mixed_kernels::cudaDataSize(Atype) * lda,
          bufferA, CUDA_R_64F, sizeof(double) * lda);
      auto bufferB_ = mixed_kernels::casting2DRowMajor(
          stream, n, k, B, Btype, mixed_kernels::cudaDataSize(Btype) * ldb,
          bufferB, CUDA_R_64F, sizeof(double) * ldb);
      return cublasGemmEx(handle, CUBLAS_OP_T, CUBLAS_OP_N, n, m, k, &alphaFP64,
                          bufferB_, CUDA_R_64F, ldb, bufferA_, CUDA_R_64F, lda,
                          &betaFP64, C, CUDA_R_64F, ldc, computeType, algo);
    }
    case CUDA_R_32F: {
      if (Atype < CUDA_R_32F && Btype < CUDA_R_32F) {
        auto bufferA_ = mixed_kernels::casting2DRowMajor(
            stream, m, k, A, Atype, mixed_kernels::cudaDataSize(Atype) * lda,
            bufferA, CUDA_R_16F, sizeof(nv_half) * lda);
        auto bufferB_ = mixed_kernels::casting2DRowMajor(
            stream, n, k, B, Btype, mixed_kernels::cudaDataSize(Btype) * ldb,
            bufferB, CUDA_R_16F, sizeof(nv_half) * ldb);
        return cublasGemmEx(handle, CUBLAS_OP_T, CUBLAS_OP_N, n, m, k,
                            &alphaFP32, bufferB_, CUDA_R_16F, ldb, bufferA_,
                            CUDA_R_16F, lda, &betaFP32, C, CUDA_R_32F, ldc,
                            CUBLAS_COMPUTE_32F, algo);
      }
      auto bufferA_ = mixed_kernels::casting2DRowMajor(
          stream, m, k, A, Atype, mixed_kernels::cudaDataSize(Atype) * lda,
          bufferA, CUDA_R_32F, sizeof(float) * lda);
      auto bufferB_ = mixed_kernels::casting2DRowMajor(
          stream, n, k, B, Btype, mixed_kernels::cudaDataSize(Btype) * ldb,
          bufferB, CUDA_R_32F, sizeof(float) * ldb);
      return cublasGemmEx(handle, CUBLAS_OP_T, CUBLAS_OP_N, n, m, k, &alphaFP32,
                          bufferB_, CUDA_R_32F, ldb, bufferA_, CUDA_R_32F, lda,
                          &betaFP32, C, CUDA_R_32F, ldc, computeType, algo);
    }
    case CUDA_R_16F: {
      auto bufferA_ = mixed_kernels::casting2DRowMajor(
          stream, m, k, A, Atype, mixed_kernels::cudaDataSize(Atype) * lda,
          bufferA, CUDA_R_16F, sizeof(nv_half) * lda);
      auto bufferB_ = mixed_kernels::casting2DRowMajor(
          stream, n, k, B, Btype, mixed_kernels::cudaDataSize(Btype) * ldb,
          bufferB, CUDA_R_16F, sizeof(nv_half) * ldb);
      switch (computeType) {
        case CUBLAS_COMPUTE_32F:
        case CUBLAS_COMPUTE_32F_PEDANTIC:
          return cublasGemmEx(handle, CUBLAS_OP_T, CUBLAS_OP_N, n, m, k,
                              &alphaFP32, bufferB_, CUDA_R_16F, ldb, bufferA_,
                              CUDA_R_16F, lda, &betaFP32, C, CUDA_R_16F, ldc,
                              computeType, algo);
        case CUBLAS_COMPUTE_16F:
        case CUBLAS_COMPUTE_16F_PEDANTIC:
          return cublasGemmEx(handle, CUBLAS_OP_T, CUBLAS_OP_N, n, m, k,
                              &alphaFP16, bufferB_, CUDA_R_16F, ldb, bufferA_,
                              CUDA_R_16F, lda, &betaFP16, C, CUDA_R_16F, ldc,
                              computeType, algo);
        default:
          return CUBLAS_STATUS_INVALID_VALUE;
      }
    }
    case CUDA_R_8F_E4M3: {
      auto bufferA_ = mixed_kernels::casting2DRowMajor(
          stream, m, k, A, Atype, mixed_kernels::cudaDataSize(Atype) * lda,
          bufferA, CUDA_R_8F_E4M3, sizeof(uint8_t) * lda);
      auto bufferB_ = mixed_kernels::casting2DRowMajor(
          stream, n, k, B, Btype, mixed_kernels::cudaDataSize(Btype) * ldb,
          bufferB, CUDA_R_8F_E4M3, sizeof(uint8_t) * ldb);
      auto castingBuffer = alignUp<16>(static_cast<std::uint8_t*>(bufferA) +
                                       static_cast<std::size_t>(k) * m);
      std::size_t castingBufferSize = static_cast<std::uint8_t*>(bufferB) -
                                      static_cast<std::uint8_t*>(castingBuffer);
      switch (computeType) {
        case CUBLAS_COMPUTE_32F:
        case CUBLAS_COMPUTE_32F_PEDANTIC:
          return mixed_kernels::cublasLtFp8RowMajorNTNMeta::matmul(
              cublasLtMeta, stream, alpha, bufferA_, bufferB_, beta, C,
              CUDA_R_8F_E4M3, castingBuffer, castingBufferSize);
        default:
          return CUBLAS_STATUS_INVALID_VALUE;
      }
    }
    default:
      return CUBLAS_STATUS_INVALID_VALUE;
  }
}

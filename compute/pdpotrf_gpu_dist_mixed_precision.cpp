/**
 *
 * @file pdpotrf.c
 *
 *  PLASMA auxiliary routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Jakub Kurzak
 * @author Hatem Ltaief
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @generated d Fri Apr  1 11:02:57 2016
 *
 **/
#include <cublasLt.h>

#include <atomic>
#include <memory>
#include <mutex>
#include <thread>
#include <unordered_map>
#include <vector>

#include "common.h"
#include "mixed_kernels.h"
#include "mixed_precision.h"
#include "plasma_d_mixed.h"

namespace {

inline cusolverStatus_t CUSOLVERAPI cusolverDnPotrfRowMajor(
    cusolverDnHandle_t handle, cusolverDnParams_t params, cublasFillMode_t uplo,
    int64_t n, cudaDataType dataTypeA, void *A, int64_t lda,
    cudaDataType computeType, void *pBuffer, size_t workspaceInBytes,
    int *info) {
  if (uplo == CUBLAS_FILL_MODE_LOWER) {
    uplo = CUBLAS_FILL_MODE_UPPER;
  } else {
    return CUSOLVER_STATUS_INVALID_VALUE;
  }
  return cusolverDnPotrf(handle, params, uplo, n, dataTypeA, A, lda,
                         computeType, pBuffer, workspaceInBytes, info);
}

CUBLASAPI cublasStatus_t CUBLASWINAPI mixedKernelsTrsmExRowMajor(
    cublasHandle_t handle, cublasSideMode_t side, cublasFillMode_t uplo,
    cublasOperation_t trans, cublasDiagType_t diag, int m, int n,
    const double alpha, const void *A, cudaDataType Atype, int lda, void *B,
    cudaDataType Btype, int ldb, void *castBuffer) {
  cudaStream_t stream;
  cublasGetStream_v2(handle, &stream);
  if (CUBLAS_OP_T == trans) {
    castBuffer =
        castToBuffer(m, n, B, Btype, ldb, castBuffer, Atype, ldb, stream);
  } else {
    return CUBLAS_STATUS_INVALID_VALUE;
  }
  if (CUBLAS_SIDE_RIGHT == side) {
    side = CUBLAS_SIDE_LEFT;
  } else {
    return CUBLAS_STATUS_INVALID_VALUE;
  }

  if (CUBLAS_FILL_MODE_LOWER == uplo) {
    uplo = CUBLAS_FILL_MODE_UPPER;
  } else {
    return CUBLAS_STATUS_INVALID_VALUE;
  }

  cublasStatus_t status;
  switch (Atype) {
    case CUDA_R_64F:
      status = cublasDtrsm_v2(handle, side, uplo, trans, diag, n, m,
                              static_cast<const double *>(&alpha),
                              static_cast<const double *>(A), lda,
                              static_cast<double *>(castBuffer), ldb);
      break;
    case CUDA_R_32F: {
      const float alphaReal = alpha;
      status = cublasStrsm_v2(handle, side, uplo, trans, diag, n, m,
                              static_cast<const float *>(&alphaReal),
                              static_cast<const float *>(A), lda,
                              static_cast<float *>(castBuffer), ldb);
      break;
    }
    default:
      printf("Not supported combination at %s:%d\n", __FILE__, __LINE__);
      status = CUBLAS_STATUS_INVALID_VALUE;
      break;
  }
  castToBuffer(m, n, castBuffer, Atype, ldb, B, Btype, ldb, stream);
  return status;
}

}  // namespace

inline MixedPrecisionTile &getTile(MixedPrecisionTiledArray &array, int m,
                                   int n) {
  return array.tiles[n + m * array.nt];
}

inline const MixedPrecisionTile &getTile(const MixedPrecisionTiledArray &array,
                                         int m, int n) {
  return array.tiles[n + m * array.nt];
}

inline void *getAddr(const MixedPrecisionTiledArray &array, int m, int n) {
  return getTile(array, m, n).data;
}
#define A(m, n) (double *)getAddr(A, m, n)
#undef BLKLDD
inline size_t BLKLDD(const MixedPrecisionTiledArray &array, int m) {
  return array.tiles[m + (array.nt - 1) * array.nt].ld;
}
/***************************************************************************/
/**
 *  Parallel tile Cholesky factorization - static scheduling
 **/

namespace {
void cudaDeleter(double *p) { CHECK_CUDA(cudaFree(p)); }
void *cudaAllocator(size_t size) {
  void *p;
  CHECK_CUDA(cudaMalloc(&p, size));
  return p;
}
}  // namespace

namespace {
struct DataTable {
  union Coord {
    std::pair<int, int> pair;
    long dummy;
  };
  static_assert(sizeof(Coord) == 2 * sizeof(int),
                "We assume long is doubled int, please report this bug");
  std::mutex mutex{};
  using MapValue = std::tuple<std::shared_ptr<double>, bool>;
  using Map = std::unordered_map<long, MapValue>;
  Map map{};
  std::shared_ptr<double> diagonal;
  int diagIdx{-1};
  int mapBound{-1};
  enum { SHARED_PTR = 0, FLAG = 1 };
  bool canAllocateNew{true};
};

std::tuple<std::shared_ptr<double>, bool, bool *> lookupTable(
    DataTable::Map &map, int mapBound, std::pair<int, int> pair, size_t size,
    bool canAllocateNew = true) {
  DataTable::Coord coord{.pair = pair};
  std::shared_ptr<double> buffer{nullptr};
  auto find = map.find(coord.dummy);
  bool exist = false;
  if (find != map.end()) {
    buffer = std::get<DataTable::SHARED_PTR>(find->second);
    exist = true;
  } else {
    if (!canAllocateNew || map.size() >= mapBound) {
      do {
        for (auto &[key, value] : map) {
          if (std::get<DataTable::SHARED_PTR>(value).use_count() == 1) {
            buffer = std::get<DataTable::SHARED_PTR>(value);
            map.insert(
                {coord.dummy,
                 {std::move(std::get<DataTable::SHARED_PTR>(value)), false}});
            find = map.find(coord.dummy);
            map.erase(key);
            break;
          }
        }
        // plasma_yield();
      } while (buffer == nullptr);
    } else {
      buffer.reset((double *)cudaAllocator(size), cudaDeleter);
      map.insert({coord.dummy, {buffer, false}});
      find = map.find(coord.dummy);
    }
  }
  return {std::move(buffer), exist, &std::get<DataTable::FLAG>(find->second)};
}

void setLandedTable(DataTable::Map &map, std::pair<int, int> pair, bool flag) {
  DataTable::Coord coord{.pair = pair};
  auto find = map.find(coord.dummy);
  if (find != map.end()) {
    std::get<DataTable::FLAG>(find->second) = flag;
  } else {
    printf("Unexpected situation happened at %s:%d, please report the bug\n",
           __FILE__, __LINE__);
  }
}

int getMapBound(int nb) {
  size_t freeSize, totalSize;
  cudaMemGetInfo(&freeSize, &totalSize);
  // NOTE: Jie: Use 90% free memory
  return freeSize / sizeof(double) * 9 / nb / nb / 10;
}
}  // namespace

void plasma_pdpotrf_gpu_reuse_data_mixed_precision_dist(
    plasma_context_t *plasma) {
  const int rank = PLASMA_RANK;
  double start, end;
  const int device = plasma->rankMapping[rank];
  cudaSetDevice(device);
  if (rank == 0) {
    plasma->volumeCPU2GPU = 0;
    plasma->volumeGPU2CPU = 0;
  }
  static_assert(sizeof(std::atomic<decltype(plasma->volumeCPU2GPU)>) ==
                sizeof(decltype(plasma->volumeCPU2GPU)));
  static_assert(sizeof(std::atomic<decltype(plasma->volumeGPU2CPU)>) ==
                sizeof(decltype(plasma->volumeGPU2CPU)));
  PLASMA_enum uplo;
  MixedPrecisionTiledArray A;
  PLASMA_sequence *sequence;
  PLASMA_request *request;

  int k, m, n;
  int next_k;
  int next_m;
  int next_n;
  int ldak, ldam, ldan;
  int *infoDevice;
  int info;
  int tempkn, tempmn;

  double zone = (double)1.0;
  double mzone = (double)-1.0;

  plasma_unpack_args_4(uplo, A, sequence, request);
  if (uplo == PlasmaUpper) {
    printf("Upper is not supported yet\n");
    return;
  }
  if (sequence->status != PLASMA_SUCCESS) return;
  ss_init(A.nt, A.nt, 0);
  // synced in ss_init
#define volumeCPU2GPU                                                \
  (reinterpret_cast<std::atomic<decltype(plasma->volumeCPU2GPU)> &>( \
      plasma->volumeCPU2GPU))
#define volumeGPU2CPU                                                \
  (reinterpret_cast<std::atomic<decltype(plasma->volumeGPU2CPU)> &>( \
      plasma->volumeGPU2CPU))

  k = 0;
  m = PLASMA_RANK;
  while (m >= A.nt) {
    k++;
    m = m - A.nt + k;
  }
  n = 0;
  int kDiag = -1;

  size_t workspaceInBytes;
  CHECK_CUSOLVER(cusolverDnPotrf_bufferSize(
      plasma->cusolver_handle[rank], plasma->cusolver_params[rank],
      CUBLAS_FILL_MODE_LOWER, A.nb, CUDA_R_64F, NULL, A.nb, CUDA_R_64F,
      &workspaceInBytes));
  double *buffer, *bufferA, *bufferB, *diagonal, *localDestination;
  const size_t bufferSizeInByte = (unsigned long)(A.nb) * A.nb * sizeof(double);
  cudaMalloc((void **)&buffer, 4 * bufferSizeInByte);
  bufferA = buffer;
  bufferB = bufferA + (unsigned long)(A.nb) * A.nb;
  diagonal = bufferB + (unsigned long)(A.nb) * A.nb;
  localDestination = diagonal + (unsigned long)(A.nb) * A.nb;
  void *pBuffer;
  CHECK_CUDA(
      cudaMalloc(&pBuffer, workspaceInBytes +
                               (sizeof(int) - workspaceInBytes % sizeof(int)) +
                               sizeof(int)));
  void *cublasLtMeta;
  mixed_kernels::cublasLtFp8RowMajorNTNMeta::create(&cublasLtMeta, A.nb, A.nb,
                                                    A.nb, A.nb, A.nb, A.nb,
                                                    workspaceInBytes, pBuffer);
  infoDevice = (int *)((long)pBuffer + workspaceInBytes +
                       (sizeof(int) - workspaceInBytes % sizeof(int)));

  double *castBufferA;
  CHECK_CUDA(cudaMalloc(reinterpret_cast<void **>(&castBufferA),
                        (unsigned long)(A.nb) * A.nb * sizeof(double) * 2));

  double one = 1., none = -1.;

  while (k < A.nt && m < A.nt && !ss_aborted()) {
    next_n = n;
    next_m = m;
    next_k = k;

    next_n++;
    if (next_n > next_k) {
      next_m += PLASMA_SIZE;
      while (next_m >= A.nt && next_k < A.nt) {
        next_k++;
        next_m = next_m - A.nt + next_k;
      }
      next_n = 0;
    }

    tempkn = k == A.nt - 1 ? A.n - k * A.nb : A.nb;
    tempmn = m == A.nt - 1 ? A.n - m * A.nb : A.nb;

    ldak = BLKLDD(A, k);
    ldan = BLKLDD(A, n);
    ldam = BLKLDD(A, m);

    if (m == k) {
      if (n == k) {
        const auto &tile = getTile(A, k, k);
        const auto sizeofTileElement = getSizeofTileElement(tile.dtype);
        if (k == 0) {
          start = get_current_time();
          CHECK_CUDA(cudaMemcpy2DAsync(
              localDestination, ldak * sizeofTileElement, A(k, k),
              ldak * sizeofTileElement, tempkn * sizeofTileElement, tempkn,
              cudaMemcpyHostToDevice, plasma->cuda_stream[rank]));
          volumeCPU2GPU += ldak * sizeofTileElement * tempkn;
          CHECK_CUDA(cudaStreamSynchronize(plasma->cuda_stream[rank]));
          end = get_current_time();
          log_event(rank, device, EVENT_C2G, start, end);
        }
        start = get_current_time();
        /*
         *  PlasmaLower
         */
        if (uplo == PlasmaLower) {
          CHECK_CUSOLVER(cusolverDnPotrfRowMajor(
              plasma->cusolver_handle[rank], plasma->cusolver_params[rank],
              CUBLAS_FILL_MODE_LOWER, tempkn, (cudaDataType)tile.dtype,
              localDestination, ldak, (cudaDataType)tile.dtype, pBuffer,
              workspaceInBytes, infoDevice));
          //                    CORE_dpotrf(
          //                        PlasmaLower,
          //                        tempkn,
          //                        A(k, k), ldak,
          //                        &info);
        }
        /*
         *  PlasmaUpper
         */
        else {
          CORE_dpotrf(PlasmaUpper, tempkn, A(k, k), ldak, &info);
        }
        CHECK_CUDA(cudaStreamSynchronize(plasma->cuda_stream[rank]));
        end = get_current_time();
        log_event(rank, device, EVENT_COMPUTE, start, end);
        start = get_current_time();
        CHECK_CUDA(cudaMemcpy2DAsync(
            A(k, k), ldak * sizeofTileElement, localDestination,
            ldak * sizeofTileElement, tempkn * sizeofTileElement, tempkn,
            cudaMemcpyDeviceToHost, plasma->cuda_stream[rank]));
        volumeGPU2CPU += ldak * sizeofTileElement * tempkn;
        CHECK_CUDA(cudaMemcpyAsync(&info, infoDevice, sizeof(int),
                                   cudaMemcpyDeviceToHost,
                                   plasma->cuda_stream[rank]));
        volumeGPU2CPU += sizeof(int);
        CHECK_CUDA(cudaStreamSynchronize(plasma->cuda_stream[rank]));
        if (info != 0) {
          plasma_request_fail(sequence, request, info + A.nb * k);
          ss_abort();
        }
        end = get_current_time();
        log_event(rank, device, EVENT_G2C, start, end);
        ss_cond_set(k, k, 1);
      } else {
        const auto &tileLocalDestination = getTile(A, k, k);
        const auto sizeofTileElementLocalDestination =
            getSizeofTileElement(tileLocalDestination.dtype);
        if (n == 0) {
          start = get_current_time();
          CHECK_CUDA(cudaMemcpy2DAsync(
              localDestination, ldak * sizeofTileElementLocalDestination,
              A(k, k), ldak * sizeofTileElementLocalDestination,
              tempkn * sizeofTileElementLocalDestination, tempkn,
              cudaMemcpyHostToDevice, plasma->cuda_stream[rank]));
          volumeCPU2GPU += ldak * sizeofTileElementLocalDestination * tempkn;
          end = get_current_time();
          log_event(rank, device, EVENT_C2G, start, end);
        }

        const auto &tileBufferA = getTile(A, k, n);
        const auto sizeofTileElementBufferA =
            getSizeofTileElement(tileBufferA.dtype);
        ss_cond_wait(k, n, 1);
        start = get_current_time();
        CHECK_CUDA(cudaMemcpy2DAsync(
            bufferA, ldan * sizeofTileElementBufferA, A(k, n),
            ldan * sizeofTileElementBufferA, A.nb * sizeofTileElementBufferA,
            tempkn, cudaMemcpyHostToDevice, plasma->cuda_stream[rank]));
        volumeCPU2GPU += ldan * sizeofTileElementBufferA * A.nb;
        CHECK_CUDA(cudaStreamSynchronize(plasma->cuda_stream[rank]));
        end = get_current_time();
        log_event(rank, device, EVENT_C2G, start, end);

        start = get_current_time();
        /*
         *  PlasmaLower
         */
        if (uplo == PlasmaLower) {
          CHECK_CUBLAS(mixedKernelsSyrkExRowMajorNTN(
              plasma->cublas_handle[rank], cublasLtMeta, CUBLAS_FILL_MODE_LOWER,
              tempkn, A.nb, none, bufferA, (cudaDataType)tileBufferA.dtype,
              ldan, one, localDestination,
              (cudaDataType)tileLocalDestination.dtype, ldak,
              (cublasComputeType_t)getComputeType(tileLocalDestination.dtype),
              CUBLAS_GEMM_DEFAULT, castBufferA,
              (unsigned long)(A.nb) * A.nb * sizeof(double) * 2));
          //                    CORE_dsyrk(
          //                         PlasmaLower, PlasmaNoTrans,
          //                         tempkn, A.nb,
          //                         -1.0, A(k, n), ldak,
          //                          1.0, A(k, k), ldak);
        }
        /*
         *  PlasmaUpper
         */
        else {
          CORE_dsyrk(PlasmaUpper, PlasmaTrans, tempkn, A.nb, -1.0, A(n, k),
                     ldan, 1.0, A(k, k), ldak);
        }
        CHECK_CUDA(cudaStreamSynchronize(plasma->cuda_stream[rank]));
        end = get_current_time();
        log_event(rank, device, EVENT_COMPUTE, start, end);
      }
    } else {
      if (n == k) {
        const auto &tileLocalDestination = getTile(A, m, k);
        const auto sizeofTileElementLocalDestination =
            getSizeofTileElement(tileLocalDestination.dtype);
        if (k == 0) {
          start = get_current_time();
          CHECK_CUDA(cudaMemcpy2DAsync(
              localDestination, ldak * sizeofTileElementLocalDestination,
              A(m, k), ldak * sizeofTileElementLocalDestination,
              A.nb * sizeofTileElementLocalDestination, tempmn,
              cudaMemcpyHostToDevice, plasma->cuda_stream[rank]));
          volumeCPU2GPU += ldak * sizeofTileElementLocalDestination * A.nb;
          end = get_current_time();
          log_event(rank, device, EVENT_C2G, start, end);
        }

        const auto &tileDiagonal = getTile(A, k, k);
        const auto sizeofTileElementDiagonal =
            getSizeofTileElement(tileDiagonal.dtype);
        ss_cond_wait(k, k, 1);
        if (kDiag != k) {
          start = get_current_time();
          CHECK_CUDA(cudaMemcpy2DAsync(
              diagonal, ldak * sizeofTileElementDiagonal, A(k, k),
              ldak * sizeofTileElementDiagonal,
              A.nb * sizeofTileElementDiagonal, A.nb, cudaMemcpyHostToDevice,
              plasma->cuda_stream[rank]));
          volumeCPU2GPU += ldak * sizeofTileElementDiagonal * A.nb;
          kDiag = k;
          CHECK_CUDA(cudaStreamSynchronize(plasma->cuda_stream[rank]));
          end = get_current_time();
          log_event(rank, device, EVENT_C2G, start, end);
        }

        start = get_current_time();
        /*
         *  PlasmaLower
         */
        if (uplo == PlasmaLower) {
          CHECK_CUBLAS(mixedKernelsTrsmExRowMajor(
              plasma->cublas_handle[rank], CUBLAS_SIDE_RIGHT,
              CUBLAS_FILL_MODE_LOWER, CUBLAS_OP_T, CUBLAS_DIAG_NON_UNIT, tempmn,
              A.nb, zone, diagonal, (cudaDataType)tileDiagonal.dtype, ldak,
              localDestination, (cudaDataType)tileLocalDestination.dtype, ldak,
              castBufferA));
          //                    CORE_dtrsm(
          //                        PlasmaRight, PlasmaLower, PlasmaTrans,
          //                        PlasmaNonUnit, tempmn, A.nb, zone, A(k, k),
          //                        ldak,
          //                              A(m, k), ldam);
        }
        /*
         *  PlasmaUpper
         */
        else {
          CORE_dtrsm(PlasmaLeft, PlasmaUpper, PlasmaTrans, PlasmaNonUnit, A.nb,
                     tempmn, zone, A(k, k), ldak, A(k, m), ldak);
        }
        CHECK_CUDA(cudaStreamSynchronize(plasma->cuda_stream[rank]));
        end = get_current_time();
        log_event(rank, device, EVENT_COMPUTE, start, end);
        start = get_current_time();
        CHECK_CUDA(cudaMemcpy2DAsync(
            A(m, k), ldak * sizeofTileElementLocalDestination, localDestination,
            ldak * sizeofTileElementLocalDestination,
            A.nb * sizeofTileElementLocalDestination, tempmn,
            cudaMemcpyDeviceToHost, plasma->cuda_stream[rank]));
        volumeGPU2CPU += ldak * sizeofTileElementLocalDestination * A.nb;
        CHECK_CUDA(cudaStreamSynchronize(plasma->cuda_stream[rank]));
        end = get_current_time();
        log_event(rank, device, EVENT_G2C, start, end);
        ss_cond_set(m, k, 1);
      } else {
        const auto &tileLocalDestination = getTile(A, m, k);
        const auto sizeofTileElementLocalDestination =
            getSizeofTileElement(tileLocalDestination.dtype);
        if (n == 0) {
          start = get_current_time();
          CHECK_CUDA(cudaMemcpy2DAsync(
              localDestination, ldak * sizeofTileElementLocalDestination,
              A(m, k), ldak * sizeofTileElementLocalDestination,
              A.nb * sizeofTileElementLocalDestination, tempmn,
              cudaMemcpyHostToDevice, plasma->cuda_stream[rank]));
          volumeCPU2GPU += ldak * sizeofTileElementLocalDestination * A.nb;
          CHECK_CUDA(cudaStreamSynchronize(plasma->cuda_stream[rank]));
          end = get_current_time();
          log_event(rank, device, EVENT_C2G, start, end);
        }

        const auto &tileBufferB = getTile(A, k, n);
        const auto sizeofTileElementBufferB =
            getSizeofTileElement(tileBufferB.dtype);
        ss_cond_wait(k, n, 1);
        start = get_current_time();
        CHECK_CUDA(cudaMemcpy2DAsync(
            bufferB, ldan * sizeofTileElementBufferB, A(k, n),
            ldan * sizeofTileElementBufferB, A.nb * sizeofTileElementBufferB,
            A.nb, cudaMemcpyHostToDevice, plasma->cuda_stream[rank]));
        volumeCPU2GPU += ldan * sizeofTileElementBufferB * A.nb;
        CHECK_CUDA(cudaStreamSynchronize(plasma->cuda_stream[rank]));
        end = get_current_time();
        log_event(rank, device, EVENT_C2G, start, end);

        const auto &tileBufferA = getTile(A, m, n);
        const auto sizeofTileElementBufferA =
            getSizeofTileElement(tileBufferA.dtype);
        ss_cond_wait(m, n, 1);
        start = get_current_time();
        CHECK_CUDA(cudaMemcpy2DAsync(
            bufferA, ldan * sizeofTileElementBufferA, A(m, n),
            ldan * sizeofTileElementBufferA, A.nb * sizeofTileElementBufferA,
            tempmn, cudaMemcpyHostToDevice, plasma->cuda_stream[rank]));
        volumeCPU2GPU += ldan * sizeofTileElementBufferA * A.nb;
        CHECK_CUDA(cudaStreamSynchronize(plasma->cuda_stream[rank]));
        end = get_current_time();
        log_event(rank, device, EVENT_C2G, start, end);
        start = get_current_time();
        /*
         *  PlasmaLower
         */
        if (uplo == PlasmaLower) {
          CHECK_CUBLAS(mixedKernelsGemmExRowMajorNTN(
              plasma->cublas_handle[rank], cublasLtMeta, tempmn, A.nb, A.nb,
              mzone, bufferA, (cudaDataType)tileBufferA.dtype, ldan, bufferB,
              (cudaDataType)tileBufferB.dtype, ldan, zone, localDestination,
              (cudaDataType)tileLocalDestination.dtype, ldak,
              (cublasComputeType_t)getComputeType(tileLocalDestination.dtype),
              CUBLAS_GEMM_DEFAULT, castBufferA,
              (unsigned long)(A.nb) * A.nb * sizeof(double) * 2));
          //                    CORE_dgemm(
          //                        PlasmaNoTrans, PlasmaTrans,
          //                        tempmn, A.nb, A.nb,
          //                        mzone, A(m, n), ldam,
          //                               A(k, n), ldak,
          //                         zone, A(m, k), ldam);
        }
        /*
         *  PlasmaUpper
         */
        else {
          CORE_dgemm(PlasmaTrans, PlasmaNoTrans, A.nb, tempmn, A.nb, mzone,
                     A(n, k), ldan, A(n, m), ldan, zone, A(k, m), ldak);
        }
        CHECK_CUDA(cudaStreamSynchronize(plasma->cuda_stream[rank]));
        end = get_current_time();
        log_event(rank, device, EVENT_COMPUTE, start, end);
      }
    }
    n = next_n;
    m = next_m;
    k = next_k;
  }
  CHECK_CUDA(cudaFree(pBuffer));
  CHECK_CUDA(cudaFree(buffer));
  CHECK_CUDA(cudaFree(castBufferA));
  CHECK_CUBLAS(
      mixed_kernels::cublasLtFp8RowMajorNTNMeta::destroy(cublasLtMeta));
  ss_finalize();
#undef volumeCPU2GPU
#undef volumeGPU2CPU
}

void plasma_pdpotrf_gpu_reuse_data_table_mixed_precision_dist(
    plasma_context_t *plasma) {
  const int rank = PLASMA_RANK;
  const int device = plasma->rankMapping[rank];
  cudaSetDevice(device);
  double start, end;
  if (rank == 0) {
    plasma->volumeCPU2GPU = 0;
    plasma->volumeGPU2CPU = 0;
  }

  if (rank == plasma->rootRankOfDevice[device]) {
    ((DataTable **)plasma->dataTable)[device] = new DataTable{};
  }
  static_assert(sizeof(std::atomic<decltype(plasma->volumeCPU2GPU)>) ==
                sizeof(decltype(plasma->volumeCPU2GPU)));
  static_assert(sizeof(std::atomic<decltype(plasma->volumeGPU2CPU)>) ==
                sizeof(decltype(plasma->volumeGPU2CPU)));
  PLASMA_enum uplo;
  MixedPrecisionTiledArray A;
  PLASMA_sequence *sequence;
  PLASMA_request *request;

  int k, m, n;
  int next_k;
  int next_m;
  int next_n;
  int ldak, ldam, ldan;
  int *infoDevice;
  int info;
  int tempkn, tempmn;

  double zone = (double)1.0;
  double mzone = (double)-1.0;

  plasma_unpack_args_4(uplo, A, sequence, request);
  if (uplo == PlasmaUpper) {
    printf("Upper is not supported yet\n");
    return;
  }
  if (sequence->status != PLASMA_SUCCESS) return;
#define volumeCPU2GPU                                                \
  (reinterpret_cast<std::atomic<decltype(plasma->volumeCPU2GPU)> &>( \
      plasma->volumeCPU2GPU))
#define volumeGPU2CPU                                                \
  (reinterpret_cast<std::atomic<decltype(plasma->volumeGPU2CPU)> &>( \
      plasma->volumeGPU2CPU))

  k = 0;
  m = PLASMA_RANK;
  while (m >= A.nt) {
    k++;
    m = m - A.nt + k;
  }
  n = 0;

  size_t workspaceInBytes;
  CHECK_CUSOLVER(cusolverDnPotrf_bufferSize(
      plasma->cusolver_handle[rank], plasma->cusolver_params[rank],
      CUBLAS_FILL_MODE_LOWER, A.nb, CUDA_R_64F, NULL, A.nb, CUDA_R_64F,
      &workspaceInBytes));
  //  double *buffer, *bufferA, *bufferB, *bufferC;
  const size_t bufferSizeInByte = (unsigned long)(A.nb) * A.nb * sizeof(double);
  std::shared_ptr<double> localDestination{
      (double *)cudaAllocator(bufferSizeInByte), cudaDeleter};
  //  cudaMalloc((void **)&buffer, 3 * bufferSizeInByte);
  //  bufferA = buffer;
  //  bufferB = bufferA + (unsigned long)(A.nb) * A.nb;
  //  bufferC = bufferB + (unsigned long)(A.nb) * A.nb;
  void *pBuffer;
  CHECK_CUDA(
      cudaMalloc(&pBuffer, workspaceInBytes +
                               (sizeof(int) - workspaceInBytes % sizeof(int)) +
                               sizeof(int)));
  void *cublasLtMeta;
  mixed_kernels::cublasLtFp8RowMajorNTNMeta::create(&cublasLtMeta, A.nb, A.nb,
                                                    A.nb, A.nb, A.nb, A.nb,
                                                    workspaceInBytes, pBuffer);
  infoDevice = (int *)((long)pBuffer + workspaceInBytes +
                       (sizeof(int) - workspaceInBytes % sizeof(int)));

  double *castBufferA;
  CHECK_CUDA(cudaMalloc(reinterpret_cast<void **>(&castBufferA),
                        (unsigned long)(A.nb) * A.nb * sizeof(double) * 2));

  double one = 1., none = -1.;
  ss_init(A.nt, A.nt, 0);
  // synced in ss_init
  auto dataTable = ((DataTable **)plasma->dataTable)[device];
  if (rank == plasma->rootRankOfDevice[device]) {
    // Do not have to sync because the first potrf is done by rank-0
    dataTable->diagonal.reset((double *)cudaAllocator(bufferSizeInByte),
                              cudaDeleter);
    const auto freeBlocks = getMapBound(A.nb);
    // Do not have to sync because the first potrf is done by rank-0
    dataTable->mapBound = max(freeBlocks / plasma->uniqueDevice,
                              min(plasma->uniqueDevice * 3, freeBlocks));
    if (dataTable->mapBound < plasma->uniqueDevice * 3) {
      printf("GPU out of memory, mostly because NB is too large\n");
      ss_abort();
    }
  }

  while (k < A.nt && m < A.nt && !ss_aborted()) {
    next_n = n;
    next_m = m;
    next_k = k;

    next_n++;
    if (next_n > next_k) {
      next_m += PLASMA_SIZE;
      while (next_m >= A.nt && next_k < A.nt) {
        next_k++;
        next_m = next_m - A.nt + next_k;
      }
      next_n = 0;
    }

    tempkn = k == A.nt - 1 ? A.n - k * A.nb : A.nb;
    tempmn = m == A.nt - 1 ? A.n - m * A.nb : A.nb;

    ldak = BLKLDD(A, k);
    ldan = BLKLDD(A, n);
    ldam = BLKLDD(A, m);

    if (m == k) {
      if (n == k) {
        const auto &tile = getTile(A, k, k);
        const auto sizeofTileElement = getSizeofTileElement(tile.dtype);
        if (k == 0) {
          start = get_current_time();
          CHECK_CUDA(cudaMemcpy2DAsync(
              localDestination.get(), ldak * sizeofTileElement, A(k, k),
              ldak * sizeofTileElement, tempkn * sizeofTileElement, tempkn,
              cudaMemcpyHostToDevice, plasma->cuda_stream[rank]));
          volumeCPU2GPU += ldak * sizeofTileElement * tempkn;
          CHECK_CUDA(cudaStreamSynchronize(plasma->cuda_stream[rank]));
          end = get_current_time();
          log_event(rank, device, EVENT_C2G, start, end);
        }
        start = get_current_time();
        /*
         *  PlasmaLower
         */
        if (uplo == PlasmaLower) {
          CHECK_CUSOLVER(cusolverDnPotrfRowMajor(
              plasma->cusolver_handle[rank], plasma->cusolver_params[rank],
              CUBLAS_FILL_MODE_LOWER, tempkn, (cudaDataType)tile.dtype,
              localDestination.get(), ldak, (cudaDataType)tile.dtype, pBuffer,
              workspaceInBytes, infoDevice));
          //                    CORE_dpotrf(
          //                        PlasmaLower,
          //                        tempkn,
          //                        A(k, k), ldak,
          //                        &info);
        }
        /*
         *  PlasmaUpper
         */
        else {
          CORE_dpotrf(PlasmaUpper, tempkn, A(k, k), ldak, &info);
        }
        CHECK_CUDA(cudaStreamSynchronize(plasma->cuda_stream[rank]));
        end = get_current_time();
        log_event(rank, device, EVENT_COMPUTE, start, end);
        start = get_current_time();
        CHECK_CUDA(cudaMemcpy2DAsync(
            A(k, k), ldak * sizeofTileElement, localDestination.get(),
            ldak * sizeofTileElement, tempkn * sizeofTileElement, tempkn,
            cudaMemcpyDeviceToHost, plasma->cuda_stream[rank]));
        volumeGPU2CPU += ldak * sizeofTileElement * tempkn;
        CHECK_CUDA(cudaMemcpyAsync(&info, infoDevice, sizeof(int),
                                   cudaMemcpyDeviceToHost,
                                   plasma->cuda_stream[rank]));
        volumeGPU2CPU += sizeof(int);
        CHECK_CUDA(cudaStreamSynchronize(plasma->cuda_stream[rank]));
        std::unique_lock<std::mutex> lock{dataTable->mutex};
        std::swap(localDestination, dataTable->diagonal);
        dataTable->diagIdx = k;
        lock.unlock();
        end = get_current_time();
        log_event(rank, device, EVENT_G2C, start, end);
        while (localDestination.use_count() > 1) {
          plasma_yield();
        }
        if (info != 0) {
          plasma_request_fail(sequence, request, info + A.nb * k);
          ss_abort();
        }
        ss_cond_set(k, k, 1);
      } else {
        const auto &tileLocalDestination = getTile(A, k, k);
        const auto sizeofTileElementLocalDestination =
            getSizeofTileElement(tileLocalDestination.dtype);
        if (n == 0) {
          start = get_current_time();
          CHECK_CUDA(cudaMemcpy2DAsync(
              localDestination.get(), ldak * sizeofTileElementLocalDestination,
              A(k, k), ldak * sizeofTileElementLocalDestination,
              tempkn * sizeofTileElementLocalDestination, tempkn,
              cudaMemcpyHostToDevice, plasma->cuda_stream[rank]));
          volumeCPU2GPU += ldak * sizeofTileElementLocalDestination * tempkn;
          CHECK_CUDA(cudaStreamSynchronize(plasma->cuda_stream[rank]));
          end = get_current_time();
          log_event(rank, device, EVENT_C2G, start, end);
        }

        const auto &tileBufferA = getTile(A, k, n);
        const auto sizeofTileElementBufferA =
            getSizeofTileElement(tileBufferA.dtype);
        ss_cond_wait(k, n, 1);
        std::unique_lock<std::mutex> lock{dataTable->mutex};
        auto [bufferA, existA, landedA] = lookupTable(
            dataTable->map, dataTable->mapBound, {k, n}, bufferSizeInByte);
        lock.unlock();
        if (!existA) {
          start = get_current_time();
          CHECK_CUDA(cudaMemcpy2DAsync(
              bufferA.get(), ldan * sizeofTileElementBufferA, A(k, n),
              ldan * sizeofTileElementBufferA, A.nb * sizeofTileElementBufferA,
              tempkn, cudaMemcpyHostToDevice, plasma->cuda_stream[rank]));
          volumeCPU2GPU += ldak * sizeofTileElementBufferA * A.nb;
          CHECK_CUDA(cudaStreamSynchronize(plasma->cuda_stream[rank]));
          *landedA = true;
          end = get_current_time();
          log_event(rank, device, EVENT_C2G, start, end);
        }

        while (!*landedA) {
          plasma_yield();
        }

        start = get_current_time();
        /*
         *  PlasmaLower
         */
        if (uplo == PlasmaLower) {
          CHECK_CUBLAS(mixedKernelsSyrkExRowMajorNTN(
              plasma->cublas_handle[rank], cublasLtMeta, CUBLAS_FILL_MODE_LOWER,
              tempkn, A.nb, none, bufferA.get(),
              (cudaDataType)tileBufferA.dtype, ldan, one,
              localDestination.get(), (cudaDataType)tileLocalDestination.dtype,
              ldak,
              (cublasComputeType_t)getComputeType(tileLocalDestination.dtype),
              CUBLAS_GEMM_DEFAULT, castBufferA,
              (unsigned long)(A.nb) * A.nb * sizeof(double) * 2));
          //                    CORE_dsyrk(
          //                         PlasmaLower, PlasmaNoTrans,
          //                         tempkn, A.nb,
          //                         -1.0, A(k, n), ldak,
          //                          1.0, A(k, k), ldak);
        }
        /*
         *  PlasmaUpper
         */
        else {
          CORE_dsyrk(PlasmaUpper, PlasmaTrans, tempkn, A.nb, -1.0, A(n, k),
                     ldan, 1.0, A(k, k), ldak);
        }
        CHECK_CUDA(cudaStreamSynchronize(plasma->cuda_stream[rank]));
        end = get_current_time();
        log_event(rank, device, EVENT_COMPUTE, start, end);
      }
    } else {
      if (n == k) {
        const auto &tileLocalDestination = getTile(A, m, k);
        const auto sizeofTileElementLocalDestination =
            getSizeofTileElement(tileLocalDestination.dtype);
        if (k == 0) {
          start = get_current_time();
          CHECK_CUDA(cudaMemcpy2DAsync(
              localDestination.get(), ldak * sizeofTileElementLocalDestination,
              A(m, k), ldak * sizeofTileElementLocalDestination,
              A.nb * sizeofTileElementLocalDestination, tempmn,
              cudaMemcpyHostToDevice, plasma->cuda_stream[rank]));
          volumeCPU2GPU += ldak * sizeofTileElementLocalDestination * A.nb;
          CHECK_CUDA(cudaStreamSynchronize(plasma->cuda_stream[rank]));
          end = get_current_time();
          log_event(rank, device, EVENT_C2G, start, end);
        }

        const auto &tileDiagonal = getTile(A, k, k);
        const auto sizeofTileElementDiagonal =
            getSizeofTileElement(tileDiagonal.dtype);
        ss_cond_wait(k, k, 1);
        bool *landed = nullptr;
        std::unique_lock<std::mutex> lock{dataTable->mutex};
        std::shared_ptr<double> diagonal = dataTable->diagonal;
        lock.unlock();
        if (dataTable->diagIdx != k /* will data racing? */) {
          lock.lock();
          auto [bufferA, existA, landedA] = lookupTable(
              dataTable->map, dataTable->mapBound, {k, k}, bufferSizeInByte);
          lock.unlock();
          landed = landedA;
          if (!existA) {
            start = get_current_time();
            CHECK_CUDA(cudaMemcpy2DAsync(
                bufferA.get(), ldak * sizeofTileElementDiagonal, A(k, k),
                ldak * sizeofTileElementDiagonal,
                A.nb * sizeofTileElementDiagonal, A.nb, cudaMemcpyHostToDevice,
                plasma->cuda_stream[rank]));
            volumeCPU2GPU += ldak * sizeofTileElementDiagonal * A.nb;
            CHECK_CUDA(cudaStreamSynchronize(plasma->cuda_stream[rank]));
            *landedA = true;
            end = get_current_time();
            log_event(rank, device, EVENT_C2G, start, end);
          }
          diagonal = std::move(bufferA);
        }
        if (landed != nullptr) {
          while (!*landed) {
            plasma_yield();
          }
        }
        start = get_current_time();
        /*
         *  PlasmaLower
         */
        if (uplo == PlasmaLower) {
          CHECK_CUBLAS(mixedKernelsTrsmExRowMajor(
              plasma->cublas_handle[rank], CUBLAS_SIDE_RIGHT,
              CUBLAS_FILL_MODE_LOWER, CUBLAS_OP_T, CUBLAS_DIAG_NON_UNIT, tempmn,
              A.nb, zone, diagonal.get(), (cudaDataType)tileDiagonal.dtype,
              ldak, localDestination.get(),
              (cudaDataType)tileLocalDestination.dtype, ldak, castBufferA));
          //                    CORE_dtrsm(
          //                        PlasmaRight, PlasmaLower, PlasmaTrans,
          //                        PlasmaNonUnit, tempmn, A.nb, zone, A(k, k),
          //                        ldak,
          //                              A(m, k), ldam);
        }
        /*
         *  PlasmaUpper
         */
        else {
          CORE_dtrsm(PlasmaLeft, PlasmaUpper, PlasmaTrans, PlasmaNonUnit, A.nb,
                     tempmn, zone, A(k, k), ldak, A(k, m), ldak);
        }
        CHECK_CUDA(cudaStreamSynchronize(plasma->cuda_stream[rank]));
        end = get_current_time();
        log_event(rank, device, EVENT_COMPUTE, start, end);
        start = get_current_time();
        CHECK_CUDA(cudaMemcpy2DAsync(
            A(m, k), ldak * sizeofTileElementLocalDestination,
            localDestination.get(), ldak * sizeofTileElementLocalDestination,
            A.nb * sizeofTileElementLocalDestination, tempmn,
            cudaMemcpyDeviceToHost, plasma->cuda_stream[rank]));
        volumeGPU2CPU += ldak * sizeofTileElementLocalDestination * A.nb;
        CHECK_CUDA(cudaStreamSynchronize(plasma->cuda_stream[rank]));
        end = get_current_time();
        log_event(rank, device, EVENT_G2C, start, end);
        ss_cond_set(m, k, 1);
      } else {
        const auto &tileLocalDestination = getTile(A, m, k);
        const auto sizeofTileElementLocalDestination =
            getSizeofTileElement(tileLocalDestination.dtype);
        if (n == 0) {
          start = get_current_time();
          CHECK_CUDA(cudaMemcpy2DAsync(
              localDestination.get(), ldak * sizeofTileElementLocalDestination,
              A(m, k), ldak * sizeofTileElementLocalDestination,
              A.nb * sizeofTileElementLocalDestination, tempmn,
              cudaMemcpyHostToDevice, plasma->cuda_stream[rank]));
          volumeCPU2GPU += ldak * sizeofTileElementLocalDestination * A.nb;
          CHECK_CUDA(cudaStreamSynchronize(plasma->cuda_stream[rank]));
          end = get_current_time();
          log_event(rank, device, EVENT_C2G, start, end);
        }

        const auto &tileBufferB = getTile(A, k, n);
        const auto sizeofTileElementBufferB =
            getSizeofTileElement(tileBufferB.dtype);
        ss_cond_wait(k, n, 1);
        std::unique_lock<std::mutex> lock{dataTable->mutex};
        auto [bufferB, existB, landedB] = lookupTable(
            dataTable->map, dataTable->mapBound, {k, n}, bufferSizeInByte);
        lock.unlock();
        if (!existB) {
          start = get_current_time();
          CHECK_CUDA(cudaMemcpy2DAsync(
              bufferB.get(), ldan * sizeofTileElementBufferB, A(k, n),
              ldan * sizeofTileElementBufferB, A.nb * sizeofTileElementBufferB,
              A.nb, cudaMemcpyHostToDevice, plasma->cuda_stream[rank]));
          volumeCPU2GPU += ldan * sizeofTileElementBufferB * A.nb;
          CHECK_CUDA(cudaStreamSynchronize(plasma->cuda_stream[rank]));
          *landedB = true;
          end = get_current_time();
          log_event(rank, device, EVENT_C2G, start, end);
        }

        const auto &tileBufferC = getTile(A, m, n);
        const auto sizeofTileElementBufferC =
            getSizeofTileElement(tileBufferC.dtype);
        ss_cond_wait(m, n, 1);
        lock.lock();
        auto [bufferC, existC, landedC] = lookupTable(
            dataTable->map, dataTable->mapBound, {m, n}, bufferSizeInByte);
        lock.unlock();
        if (!existC) {
          start = get_current_time();
          CHECK_CUDA(cudaMemcpy2DAsync(
              bufferC.get(), ldan * sizeofTileElementBufferC, A(m, n),
              ldan * sizeofTileElementBufferC, A.nb * sizeofTileElementBufferC,
              tempmn, cudaMemcpyHostToDevice, plasma->cuda_stream[rank]));
          volumeCPU2GPU += ldan * sizeofTileElementBufferC * A.nb;
          CHECK_CUDA(cudaStreamSynchronize(plasma->cuda_stream[rank]));
          *landedC = true;
          end = get_current_time();
          log_event(rank, device, EVENT_C2G, start, end);
        }
        while (!*landedB) {
          plasma_yield();
        }
        while (!*landedC) {
          plasma_yield();
        }
        start = get_current_time();
        /*
         *  PlasmaLower
         */
        if (uplo == PlasmaLower) {
          CHECK_CUBLAS(mixedKernelsGemmExRowMajorNTN(
              plasma->cublas_handle[rank], cublasLtMeta, tempmn, A.nb, A.nb,
              mzone, bufferC.get(), (cudaDataType)tileBufferC.dtype, ldan,
              bufferB.get(), (cudaDataType)tileBufferB.dtype, ldan, zone,
              localDestination.get(), (cudaDataType)tileLocalDestination.dtype,
              ldak,
              (cublasComputeType_t)getComputeType(tileLocalDestination.dtype),
              CUBLAS_GEMM_DEFAULT, castBufferA,
              (unsigned long)(A.nb) * A.nb * sizeof(double) * 2));
          //                    CORE_dgemm(
          //                        PlasmaNoTrans, PlasmaTrans,
          //                        tempmn, A.nb, A.nb,
          //                        mzone, A(m, n), ldam,
          //                               A(k, n), ldak,
          //                         zone, A(m, k), ldam);
        }
        /*
         *  PlasmaUpper
         */
        else {
          CORE_dgemm(PlasmaTrans, PlasmaNoTrans, A.nb, tempmn, A.nb, mzone,
                     A(n, k), ldan, A(n, m), ldan, zone, A(k, m), ldak);
        }
        CHECK_CUDA(cudaStreamSynchronize(plasma->cuda_stream[rank]));
        end = get_current_time();
        log_event(rank, device, EVENT_COMPUTE, start, end);
      }
    }
    n = next_n;
    m = next_m;
    k = next_k;
  }
  CHECK_CUDA(cudaFree(pBuffer));
  CHECK_CUDA(cudaFree(castBufferA));
  CHECK_CUBLAS(
      mixed_kernels::cublasLtFp8RowMajorNTNMeta::destroy(cublasLtMeta));
  //  CHECK_CUDA(cudaFree(buffer));
  ss_finalize();
  if (rank == plasma->rootRankOfDevice[device]) {
    delete dataTable;
  }
#undef volumeCPU2GPU
#undef volumeGPU2CPU
}

void plasma_pdpotrf_gpu_reuse_data_table_all_managed_mixed_precision_dist(
    plasma_context_t *plasma) {
  const int rank = PLASMA_RANK;
  const int device = plasma->rankMapping[rank];
  double start, end;
  cudaSetDevice(device);
  if (rank == 0) {
    plasma->volumeCPU2GPU = 0;
    plasma->volumeGPU2CPU = 0;
  }

  if (rank == plasma->rootRankOfDevice[device]) {
    ((DataTable **)plasma->dataTable)[device] = new DataTable{};
  }

  static_assert(sizeof(std::atomic<decltype(plasma->volumeCPU2GPU)>) ==
                sizeof(decltype(plasma->volumeCPU2GPU)));
  static_assert(sizeof(std::atomic<decltype(plasma->volumeGPU2CPU)>) ==
                sizeof(decltype(plasma->volumeGPU2CPU)));
  PLASMA_enum uplo;
  MixedPrecisionTiledArray A;
  PLASMA_sequence *sequence;
  PLASMA_request *request;

  int k, m, n;
  int next_k;
  int next_m;
  int next_n;
  int ldak, ldam, ldan;
  int *infoDevice;
  int info;
  int tempkn, tempmn;

  double zone = (double)1.0;
  double mzone = (double)-1.0;

  plasma_unpack_args_4(uplo, A, sequence, request);
  if (uplo == PlasmaUpper) {
    printf("Upper is not supported yet\n");
    return;
  }
  if (sequence->status != PLASMA_SUCCESS) return;
#define volumeCPU2GPU                                                \
  (reinterpret_cast<std::atomic<decltype(plasma->volumeCPU2GPU)> &>( \
      plasma->volumeCPU2GPU))
#define volumeGPU2CPU                                                \
  (reinterpret_cast<std::atomic<decltype(plasma->volumeGPU2CPU)> &>( \
      plasma->volumeGPU2CPU))

  k = 0;
  m = PLASMA_RANK;
  while (m >= A.nt) {
    k++;
    m = m - A.nt + k;
  }
  n = 0;

  size_t workspaceInBytes;
  CHECK_CUSOLVER(cusolverDnPotrf_bufferSize(
      plasma->cusolver_handle[rank], plasma->cusolver_params[rank],
      CUBLAS_FILL_MODE_LOWER, A.nb, CUDA_R_64F, NULL, A.nb, CUDA_R_64F,
      &workspaceInBytes));
  //  double *buffer, *bufferA, *bufferB, *bufferC;
  const size_t bufferSizeInByte = (unsigned long)(A.nb) * A.nb * sizeof(double);
  std::shared_ptr<double> localDestination;
  //  cudaMalloc((void **)&buffer, 3 * bufferSizeInByte);
  //  bufferA = buffer;
  //  bufferB = bufferA + (unsigned long)(A.nb) * A.nb;
  //  bufferC = bufferB + (unsigned long)(A.nb) * A.nb;
  void *pBuffer;
  CHECK_CUDA(
      cudaMalloc(&pBuffer, workspaceInBytes +
                               (sizeof(int) - workspaceInBytes % sizeof(int)) +
                               sizeof(int)));
  void *cublasLtMeta;
  mixed_kernels::cublasLtFp8RowMajorNTNMeta::create(&cublasLtMeta, A.nb, A.nb,
                                                    A.nb, A.nb, A.nb, A.nb,
                                                    workspaceInBytes, pBuffer);
  infoDevice = (int *)((long)pBuffer + workspaceInBytes +
                       (sizeof(int) - workspaceInBytes % sizeof(int)));

  double *castBufferA;
  CHECK_CUDA(cudaMalloc(reinterpret_cast<void **>(&castBufferA),
                        (unsigned long)(A.nb) * A.nb * sizeof(double) * 2));

  double one = 1., none = -1.;
  ss_init(A.nt, A.nt, 0);
  // synced in ss_init
  auto dataTable = ((DataTable **)plasma->dataTable)[device];
  if (rank == plasma->rootRankOfDevice[device]) {
    const auto freeBlocks = getMapBound(A.nb);
    // Do not have to sync because the first potrf is done by rank-0
    dataTable->mapBound = max(freeBlocks / plasma->uniqueDevice,
                              min(plasma->uniqueDevice * 3, freeBlocks));
    if (dataTable->mapBound < plasma->uniqueDevice * 3) {
      printf("GPU out of memory, mostly because NB is too large\n");
      ss_abort();
    }
  }
  plasma_barrier(plasma);

  // Initialize localDestination
  {
    next_n = n;
    next_m = m;
    next_k = k;

    next_n++;
    if (next_n > next_k) {
      next_m += PLASMA_SIZE;
      while (next_m >= A.nt && next_k < A.nt) {
        next_k++;
        next_m = next_m - A.nt + next_k;
      }
      next_n = 0;
    }
    if (k == 0) {
      std::unique_lock<std::mutex> lock{dataTable->mutex};
      auto [bufferDestination, _1, _2] = lookupTable(
          dataTable->map, dataTable->mapBound, {m, k}, bufferSizeInByte);
      lock.unlock();
      localDestination = std::move(bufferDestination);
    } else {
      std::unique_lock<std::mutex> lock{dataTable->mutex};
      auto [bufferDestination, _1, _2] =
          lookupTable(dataTable->map, dataTable->mapBound, {next_m, next_k},
                      bufferSizeInByte);
      lock.unlock();
      localDestination = std::move(bufferDestination);
    }
  }

  while (k < A.nt && m < A.nt && !ss_aborted()) {
    next_n = n;
    next_m = m;
    next_k = k;

    next_n++;
    if (next_n > next_k) {
      next_m += PLASMA_SIZE;
      while (next_m >= A.nt && next_k < A.nt) {
        next_k++;
        next_m = next_m - A.nt + next_k;
      }
      next_n = 0;
    }

    tempkn = k == A.nt - 1 ? A.n - k * A.nb : A.nb;
    tempmn = m == A.nt - 1 ? A.n - m * A.nb : A.nb;

    ldak = BLKLDD(A, k);
    ldan = BLKLDD(A, n);
    ldam = BLKLDD(A, m);

    if (m == k) {
      if (n == k) {
        const auto &tile = getTile(A, k, k);
        const auto sizeofTileElement = getSizeofTileElement(tile.dtype);
        if (k == 0) {
          start = get_current_time();
          CHECK_CUDA(cudaMemcpy2DAsync(
              localDestination.get(), ldak * sizeofTileElement, A(k, k),
              ldak * sizeofTileElement, tempkn * sizeofTileElement, tempkn,
              cudaMemcpyHostToDevice, plasma->cuda_stream[rank]));
          volumeCPU2GPU += ldak * sizeofTileElement * tempkn;
          CHECK_CUDA(cudaStreamSynchronize(plasma->cuda_stream[rank]));
          end = get_current_time();
          log_event(rank, device, EVENT_C2G, start, end);
        }
        start = get_current_time();
        /*
         *  PlasmaLower
         */
        if (uplo == PlasmaLower) {
          CHECK_CUSOLVER(cusolverDnPotrfRowMajor(
              plasma->cusolver_handle[rank], plasma->cusolver_params[rank],
              CUBLAS_FILL_MODE_LOWER, tempkn, (cudaDataType)tile.dtype,
              localDestination.get(), ldak, (cudaDataType)tile.dtype, pBuffer,
              workspaceInBytes, infoDevice));
          //                    CORE_dpotrf(
          //                        PlasmaLower,
          //                        tempkn,
          //                        A(k, k), ldak,
          //                        &info);
        }
        /*
         *  PlasmaUpper
         */
        else {
          CORE_dpotrf(PlasmaUpper, tempkn, A(k, k), ldak, &info);
        }
        CHECK_CUDA(cudaStreamSynchronize(plasma->cuda_stream[rank]));
        end = get_current_time();
        log_event(rank, device, EVENT_COMPUTE, start, end);
        start = get_current_time();
        CHECK_CUDA(cudaMemcpy2DAsync(
            A(k, k), ldak * sizeofTileElement, localDestination.get(),
            ldak * sizeofTileElement, tempkn * sizeofTileElement, tempkn,
            cudaMemcpyDeviceToHost, plasma->cuda_stream[rank]));
        volumeGPU2CPU += ldak * sizeofTileElement * tempkn;
        CHECK_CUDA(cudaStreamSynchronize(plasma->cuda_stream[rank]));
        end = get_current_time();
        log_event(rank, device, EVENT_G2C, start, end);
        CHECK_CUDA(cudaMemcpyAsync(&info, infoDevice, sizeof(int),
                                   cudaMemcpyDeviceToHost,
                                   plasma->cuda_stream[rank]));
        volumeGPU2CPU += sizeof(int);
        std::unique_lock<std::mutex> lock{dataTable->mutex};
        dataTable->diagIdx = k;
        dataTable->diagonal = std::move(localDestination);
        setLandedTable(dataTable->map, {k, k}, true);
        auto [bufferNext, existNext, landedNext] =
            lookupTable(dataTable->map, dataTable->mapBound, {next_m, next_k},
                        bufferSizeInByte);
        lock.unlock();
        if (existNext) {
          printf(
              "Unexpected situation happened at %s:%d, please report the bug\n",
              __FILE__, __LINE__);
          ss_abort();
        }
        localDestination = std::move(bufferNext);
        if (info != 0) {
          plasma_request_fail(sequence, request, info + A.nb * k);
          ss_abort();
        }
        CHECK_CUDA(cudaStreamSynchronize(plasma->cuda_stream[rank]));
        ss_cond_set(k, k, 1);
      } else {
        const auto &tileLocalDestination = getTile(A, k, k);
        const auto sizeofTileElementLocalDestination =
            getSizeofTileElement(tileLocalDestination.dtype);
        if (n == 0) {
          start = get_current_time();
          CHECK_CUDA(cudaMemcpy2DAsync(
              localDestination.get(), ldak * sizeofTileElementLocalDestination,
              A(k, k), ldak * sizeofTileElementLocalDestination,
              tempkn * sizeofTileElementLocalDestination, tempkn,
              cudaMemcpyHostToDevice, plasma->cuda_stream[rank]));
          volumeCPU2GPU += ldak * sizeofTileElementLocalDestination * tempkn;
          CHECK_CUDA(cudaStreamSynchronize(plasma->cuda_stream[rank]));
          end = get_current_time();
          log_event(rank, device, EVENT_C2G, start, end);
        }

        const auto &tileBufferA = getTile(A, k, n);
        const auto sizeofTileElementBufferA =
            getSizeofTileElement(tileBufferA.dtype);
        ss_cond_wait(k, n, 1);
        std::unique_lock<std::mutex> lock{dataTable->mutex};
        auto [bufferA, existA, landedA] = lookupTable(
            dataTable->map, dataTable->mapBound, {k, n}, bufferSizeInByte);
        lock.unlock();
        if (!existA) {
          start = get_current_time();
          CHECK_CUDA(cudaMemcpy2DAsync(
              bufferA.get(), ldan * sizeofTileElementBufferA, A(k, n),
              ldan * sizeofTileElementBufferA, A.nb * sizeofTileElementBufferA,
              tempkn, cudaMemcpyHostToDevice, plasma->cuda_stream[rank]));
          volumeCPU2GPU += ldak * sizeofTileElementBufferA * A.nb;
          CHECK_CUDA(cudaStreamSynchronize(plasma->cuda_stream[rank]));
          *landedA = true;
          end = get_current_time();
          log_event(rank, device, EVENT_C2G, start, end);
        }

        while (!*landedA) {
          plasma_yield();
        }
        start = get_current_time();
        /*
         *  PlasmaLower
         */
        if (uplo == PlasmaLower) {
          CHECK_CUBLAS(mixedKernelsSyrkExRowMajorNTN(
              plasma->cublas_handle[rank], cublasLtMeta, CUBLAS_FILL_MODE_LOWER,
              tempkn, A.nb, none, bufferA.get(),
              (cudaDataType)tileBufferA.dtype, ldan, one,
              localDestination.get(), (cudaDataType)tileLocalDestination.dtype,
              ldak,
              (cublasComputeType_t)getComputeType(tileLocalDestination.dtype),
              CUBLAS_GEMM_DEFAULT, castBufferA,
              (unsigned long)(A.nb) * A.nb * sizeof(double) * 2));
          //                    CORE_dsyrk(
          //                         PlasmaLower, PlasmaNoTrans,
          //                         tempkn, A.nb,
          //                         -1.0, A(k, n), ldak,
          //                          1.0, A(k, k), ldak);
        }
        /*
         *  PlasmaUpper
         */
        else {
          CORE_dsyrk(PlasmaUpper, PlasmaTrans, tempkn, A.nb, -1.0, A(n, k),
                     ldan, 1.0, A(k, k), ldak);
        }
        CHECK_CUDA(cudaStreamSynchronize(plasma->cuda_stream[rank]));
        end = get_current_time();
        log_event(rank, device, EVENT_COMPUTE, start, end);
      }
    } else {
      if (n == k) {
        const auto &tileLocalDestination = getTile(A, m, k);
        const auto sizeofTileElementLocalDestination =
            getSizeofTileElement(tileLocalDestination.dtype);
        if (k == 0) {
          start = get_current_time();
          CHECK_CUDA(cudaMemcpy2DAsync(
              localDestination.get(), ldak * sizeofTileElementLocalDestination,
              A(m, k), ldak * sizeofTileElementLocalDestination,
              A.nb * sizeofTileElementLocalDestination, tempmn,
              cudaMemcpyHostToDevice, plasma->cuda_stream[rank]));
          volumeCPU2GPU += ldak * sizeofTileElementLocalDestination * A.nb;
          CHECK_CUDA(cudaStreamSynchronize(plasma->cuda_stream[rank]));
          end = get_current_time();
          log_event(rank, device, EVENT_C2G, start, end);
        }

        const auto &tileDiagonal = getTile(A, k, k);
        const auto sizeofTileElementDiagonal =
            getSizeofTileElement(tileDiagonal.dtype);
        ss_cond_wait(k, k, 1);
        bool *landed = nullptr;
        std::unique_lock<std::mutex> lock{dataTable->mutex};
        std::shared_ptr<double> diagonal = dataTable->diagonal;
        lock.unlock();
        if (dataTable->diagIdx != k /* will data racing? */) {
          lock.lock();
          auto [bufferA, existA, landedA] = lookupTable(
              dataTable->map, dataTable->mapBound, {k, k}, bufferSizeInByte);
          lock.unlock();
          landed = landedA;
          if (!existA) {
            start = get_current_time();
            CHECK_CUDA(cudaMemcpy2DAsync(
                bufferA.get(), ldak * sizeofTileElementDiagonal, A(k, k),
                ldak * sizeofTileElementDiagonal,
                A.nb * sizeofTileElementDiagonal, A.nb, cudaMemcpyHostToDevice,
                plasma->cuda_stream[rank]));
            volumeCPU2GPU += ldak * sizeofTileElementDiagonal * A.nb;
            CHECK_CUDA(cudaStreamSynchronize(plasma->cuda_stream[rank]));
            *landedA = true;
            end = get_current_time();
            log_event(rank, device, EVENT_C2G, start, end);
          }
          diagonal = std::move(bufferA);
        }
        if (landed != nullptr) {
          while (!*landed) {
            plasma_yield();
          }
        }
        start = get_current_time();
        /*
         *  PlasmaLower
         */
        if (uplo == PlasmaLower) {
          CHECK_CUBLAS(mixedKernelsTrsmExRowMajor(
              plasma->cublas_handle[rank], CUBLAS_SIDE_RIGHT,
              CUBLAS_FILL_MODE_LOWER, CUBLAS_OP_T, CUBLAS_DIAG_NON_UNIT, tempmn,
              A.nb, zone, diagonal.get(), (cudaDataType)tileDiagonal.dtype,
              ldak, localDestination.get(),
              (cudaDataType)tileLocalDestination.dtype, ldak, castBufferA));
          //                    CORE_dtrsm(
          //                        PlasmaRight, PlasmaLower, PlasmaTrans,
          //                        PlasmaNonUnit, tempmn, A.nb, zone, A(k, k),
          //                        ldak,
          //                              A(m, k), ldam);
        }
        /*
         *  PlasmaUpper
         */
        else {
          CORE_dtrsm(PlasmaLeft, PlasmaUpper, PlasmaTrans, PlasmaNonUnit, A.nb,
                     tempmn, zone, A(k, k), ldak, A(k, m), ldak);
        }
        CHECK_CUDA(cudaStreamSynchronize(plasma->cuda_stream[rank]));
        end = get_current_time();
        log_event(rank, device, EVENT_COMPUTE, start, end);
        start = get_current_time();
        CHECK_CUDA(cudaMemcpy2DAsync(
            A(m, k), ldak * sizeofTileElementLocalDestination,
            localDestination.get(), ldak * sizeofTileElementLocalDestination,
            A.nb * sizeofTileElementLocalDestination, tempmn,
            cudaMemcpyDeviceToHost, plasma->cuda_stream[rank]));
        volumeGPU2CPU += ldak * sizeofTileElementLocalDestination * A.nb;
        CHECK_CUDA(cudaStreamSynchronize(plasma->cuda_stream[rank]));
        lock.lock();
        setLandedTable(dataTable->map, {m, k}, true);
        auto [bufferNext, existNext, landedNext] =
            lookupTable(dataTable->map, dataTable->mapBound, {next_m, next_k},
                        bufferSizeInByte);
        lock.unlock();
        if (existNext) {
          printf(
              "Unexpected situation happened at %s:%d, please report the bug\n",
              __FILE__, __LINE__);
          ss_abort();
        }
        localDestination = std::move(bufferNext);

        CHECK_CUDA(cudaStreamSynchronize(plasma->cuda_stream[rank]));
        end = get_current_time();
        log_event(rank, device, EVENT_G2C, start, end);
        ss_cond_set(m, k, 1);
      } else {
        const auto &tileLocalDestination = getTile(A, m, k);
        const auto sizeofTileElementLocalDestination =
            getSizeofTileElement(tileLocalDestination.dtype);
        if (n == 0) {
          start = get_current_time();
          CHECK_CUDA(cudaMemcpy2DAsync(
              localDestination.get(), ldak * sizeofTileElementLocalDestination,
              A(m, k), ldak * sizeofTileElementLocalDestination,
              A.nb * sizeofTileElementLocalDestination, tempmn,
              cudaMemcpyHostToDevice, plasma->cuda_stream[rank]));
          volumeCPU2GPU += ldak * sizeofTileElementLocalDestination * A.nb;
          CHECK_CUDA(cudaStreamSynchronize(plasma->cuda_stream[rank]));
          end = get_current_time();
          log_event(rank, device, EVENT_C2G, start, end);
        }

        const auto &tileBufferB = getTile(A, k, n);
        const auto sizeofTileElementBufferB =
            getSizeofTileElement(tileBufferB.dtype);
        ss_cond_wait(k, n, 1);
        std::unique_lock<std::mutex> lock{dataTable->mutex};
        auto [bufferB, existB, landedB] = lookupTable(
            dataTable->map, dataTable->mapBound, {k, n}, bufferSizeInByte);
        lock.unlock();
        if (!existB) {
          start = get_current_time();
          CHECK_CUDA(cudaMemcpy2DAsync(
              bufferB.get(), ldan * sizeofTileElementBufferB, A(k, n),
              ldan * sizeofTileElementBufferB, A.nb * sizeofTileElementBufferB,
              A.nb, cudaMemcpyHostToDevice, plasma->cuda_stream[rank]));
          volumeCPU2GPU += ldan * sizeofTileElementBufferB * A.nb;
          CHECK_CUDA(cudaStreamSynchronize(plasma->cuda_stream[rank]));
          *landedB = true;
          end = get_current_time();
          log_event(rank, device, EVENT_C2G, start, end);
        }

        const auto &tileBufferC = getTile(A, m, n);
        const auto sizeofTileElementBufferC =
            getSizeofTileElement(tileBufferC.dtype);
        ss_cond_wait(m, n, 1);
        lock.lock();
        auto [bufferC, existC, landedC] = lookupTable(
            dataTable->map, dataTable->mapBound, {m, n}, bufferSizeInByte);
        lock.unlock();
        if (!existC) {
          start = get_current_time();
          CHECK_CUDA(cudaMemcpy2DAsync(
              bufferC.get(), ldan * sizeofTileElementBufferC, A(m, n),
              ldan * sizeofTileElementBufferC, A.nb * sizeofTileElementBufferC,
              tempmn, cudaMemcpyHostToDevice, plasma->cuda_stream[rank]));
          volumeCPU2GPU += ldan * sizeofTileElementBufferC * A.nb;
          CHECK_CUDA(cudaStreamSynchronize(plasma->cuda_stream[rank]));
          *landedC = true;
          end = get_current_time();
          log_event(rank, device, EVENT_C2G, start, end);
        }
        while (!*landedB) {
          plasma_yield();
        }
        while (!*landedC) {
          plasma_yield();
        }
        start = get_current_time();
        /*
         *  PlasmaLower
         */
        if (uplo == PlasmaLower) {
          CHECK_CUBLAS(mixedKernelsGemmExRowMajorNTN(
              plasma->cublas_handle[rank], cublasLtMeta, tempmn, A.nb, A.nb,
              mzone, bufferC.get(), (cudaDataType)tileBufferC.dtype, ldan,
              bufferB.get(), (cudaDataType)tileBufferB.dtype, ldan, zone,
              localDestination.get(), (cudaDataType)tileLocalDestination.dtype,
              ldak,
              (cublasComputeType_t)getComputeType(tileLocalDestination.dtype),
              CUBLAS_GEMM_DEFAULT, castBufferA,
              (unsigned long)(A.nb) * A.nb * sizeof(double) * 2));
          //                    CORE_dgemm(
          //                        PlasmaNoTrans, PlasmaTrans,
          //                        tempmn, A.nb, A.nb,
          //                        mzone, A(m, n), ldam,
          //                               A(k, n), ldak,
          //                         zone, A(m, k), ldam);
        }
        /*
         *  PlasmaUpper
         */
        else {
          CORE_dgemm(PlasmaTrans, PlasmaNoTrans, A.nb, tempmn, A.nb, mzone,
                     A(n, k), ldan, A(n, m), ldan, zone, A(k, m), ldak);
        }
        CHECK_CUDA(cudaStreamSynchronize(plasma->cuda_stream[rank]));
        end = get_current_time();
        log_event(rank, device, EVENT_COMPUTE, start, end);
      }
    }
    n = next_n;
    m = next_m;
    k = next_k;
  }
  CHECK_CUDA(cudaFree(pBuffer));
  CHECK_CUDA(cudaFree(castBufferA));
  CHECK_CUBLAS(
      mixed_kernels::cublasLtFp8RowMajorNTNMeta::destroy(cublasLtMeta));
  //  CHECK_CUDA(cudaFree(buffer));
  ss_finalize();
  if (rank == plasma->rootRankOfDevice[device]) {
    delete dataTable;
  }
#undef volumeCPU2GPU
#undef volumeGPU2CPU
}

/***************************************************************************/
/**
 *  Parallel tile Cholesky factorization - dynamic scheduling
 **/
void plasma_pdpotrf_gpu_reuse_data_mixed_precision_dist_quark(
    PLASMA_enum uplo, MixedPrecisionTiledArray A, PLASMA_sequence *sequence,
    PLASMA_request *request) {
  printf("Dynamic scheduling is not supported yet\n");
  return;
}

void plasma_pdpotrf_gpu_reuse_data_table_mixed_precision_dist_quark(
    PLASMA_enum uplo, MixedPrecisionTiledArray A, PLASMA_sequence *sequence,
    PLASMA_request *request) {
  printf("Dynamic scheduling is not supported yet\n");
  return;
}

void plasma_pdpotrf_gpu_reuse_data_table_all_managed_mixed_precision_dist_quark(
    PLASMA_enum uplo, MixedPrecisionTiledArray A, PLASMA_sequence *sequence,
    PLASMA_request *request) {
  printf("Dynamic scheduling is not supported yet\n");
  return;
}

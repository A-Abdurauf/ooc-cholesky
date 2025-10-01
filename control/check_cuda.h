#pragma once
#ifdef PLASMA_WITH_CUDA
#define CHECK_CUSOLVER(statement)                                       \
  do {                                                                  \
    cusolverStatus_t error;                                             \
    error = statement;                                                  \
    if (error != CUSOLVER_STATUS_SUCCESS) {                             \
      printf("cusolver failed with code \"%d\" in %s:%d\n", (int)error, \
             __FILE__, __LINE__);                                       \
    }                                                                   \
  } while (0)

#define CHECK_CUDA(statement)                                                 \
  do {                                                                        \
    cudaError_t error;                                                        \
    error = statement;                                                        \
    if (error != cudaSuccess) {                                               \
      printf("cuda failed with code \"%d\" in %s:%d\n", (int)error, __FILE__, \
             __LINE__);                                                       \
    }                                                                         \
  } while (0)

#define CHECK_CUBLAS(statement)                                       \
  do {                                                                \
    cublasStatus_t error;                                             \
    error = statement;                                                \
    if (error != CUBLAS_STATUS_SUCCESS) {                             \
      printf("cublas failed with code \"%d\" in %s:%d\n", (int)error, \
             __FILE__, __LINE__);                                     \
    }                                                                 \
  } while (0)

#endif

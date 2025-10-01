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
#include "common.h"

#define A(m, n) BLKADDR(A, double, m, n)
/***************************************************************************/
/**
 *  Parallel tile Cholesky factorization - static scheduling
 **/
void plasma_pdpotrf_gpu(plasma_context_t *plasma) {
  PLASMA_enum uplo;
  PLASMA_desc A;
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

  k = 0;
  m = PLASMA_RANK;
  while (m >= A.nt) {
    k++;
    m = m - A.nt + k;
  }
  n = 0;

  const int rank = PLASMA_RANK;
  size_t workspaceInBytes;
  CHECK_CUSOLVER(cusolverDnPotrf_bufferSize(
      plasma->cusolver_handle[rank], plasma->cusolver_params[rank],
      CUBLAS_FILL_MODE_LOWER, A.nb, CUDA_R_64F, NULL, A.nb, CUDA_R_64F,
      &workspaceInBytes));
  double *buffer, *bufferA, *bufferB, *bufferC;
  const size_t bufferSizeInByte = (unsigned long)(A.nb) * A.nb * sizeof(double);
  cudaMalloc((void **)&buffer, 3 * bufferSizeInByte);
  bufferA = buffer;
  bufferB = bufferA + (unsigned long)(A.nb) * A.nb;
  bufferC = bufferB + (unsigned long)(A.nb) * A.nb;
  void *pBuffer;
  CHECK_CUDA(
      cudaMalloc(&pBuffer, workspaceInBytes +
                               (sizeof(int) - workspaceInBytes % sizeof(int)) +
                               sizeof(int)));
  infoDevice = pBuffer + workspaceInBytes +
               (sizeof(int) - workspaceInBytes % sizeof(int));
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
        CHECK_CUDA(cudaMemcpy2DAsync(
            bufferA, ldak * sizeof(double), A(k, k), ldak * sizeof(double),
            tempkn * sizeof(double), tempkn, cudaMemcpyHostToDevice,
            plasma->cuda_stream[rank]));
        /*
         *  PlasmaLower
         */
        if (uplo == PlasmaLower) {
          CHECK_CUSOLVER(cusolverDnPotrf(
              plasma->cusolver_handle[PLASMA_RANK],
              plasma->cusolver_params[PLASMA_RANK], CUBLAS_FILL_MODE_LOWER,
              tempkn, CUDA_R_64F, bufferA, ldak, CUDA_R_64F, pBuffer,
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
        CHECK_CUDA(cudaMemcpy2DAsync(
            A(k, k), ldak * sizeof(double), bufferA, ldak * sizeof(double),
            tempkn * sizeof(double), tempkn, cudaMemcpyDeviceToHost,
            plasma->cuda_stream[rank]));
        CHECK_CUDA(cudaMemcpyAsync(&info, infoDevice, sizeof(int),
                                   cudaMemcpyDeviceToHost,
                                   plasma->cuda_stream[rank]));
        CHECK_CUDA(cudaStreamSynchronize(plasma->cuda_stream[rank]));
        if (info != 0) {
          plasma_request_fail(sequence, request, info + A.nb * k);
          ss_abort();
        }
        ss_cond_set(k, k, 1);
      } else {
        CHECK_CUDA(cudaMemcpy2DAsync(
            bufferB, ldak * sizeof(double), A(k, k), ldak * sizeof(double),
            tempkn * sizeof(double), tempkn, cudaMemcpyHostToDevice,
            plasma->cuda_stream[rank]));
        ss_cond_wait(k, n, 1);
        CHECK_CUDA(cudaMemcpy2DAsync(
            bufferA, ldak * sizeof(double), A(k, n), ldak * sizeof(double),
            tempkn * sizeof(double), A.nb, cudaMemcpyHostToDevice,
            plasma->cuda_stream[rank]));

        /*
         *  PlasmaLower
         */
        if (uplo == PlasmaLower) {
          CHECK_CUBLAS(cublasDsyrk(
              plasma->cublas_handle[rank], CUBLAS_FILL_MODE_LOWER, CUBLAS_OP_N,
              tempkn, A.nb, &none, bufferA, ldak, &one, bufferB, ldak));
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
        CHECK_CUDA(cudaMemcpy2DAsync(
            A(k, k), ldak * sizeof(double), bufferB, ldak * sizeof(double),
            tempkn * sizeof(double), tempkn, cudaMemcpyDeviceToHost,
            plasma->cuda_stream[rank]));
        CHECK_CUDA(cudaStreamSynchronize(plasma->cuda_stream[rank]));
      }
    } else {
      if (n == k) {
        CHECK_CUDA(cudaMemcpy2DAsync(
            bufferB, ldam * sizeof(double), A(m, k), ldam * sizeof(double),
            tempmn * sizeof(double), A.nb, cudaMemcpyHostToDevice,
            plasma->cuda_stream[rank]));
        ss_cond_wait(k, k, 1);
        CHECK_CUDA(cudaMemcpy2DAsync(
            bufferA, ldak * sizeof(double), A(k, k), ldak * sizeof(double),
            A.nb * sizeof(double), A.nb, cudaMemcpyHostToDevice,
            plasma->cuda_stream[rank]));
        /*
         *  PlasmaLower
         */
        if (uplo == PlasmaLower) {
          CHECK_CUBLAS(cublasDtrsm(plasma->cublas_handle[rank],
                                   CUBLAS_SIDE_RIGHT, CUBLAS_FILL_MODE_LOWER,
                                   CUBLAS_OP_T, CUBLAS_DIAG_NON_UNIT, tempmn,
                                   A.nb, &zone, bufferA, ldak, bufferB, ldam));
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
        CHECK_CUDA(cudaMemcpy2DAsync(
            A(m, k), ldam * sizeof(double), bufferB, ldam * sizeof(double),
            tempmn * sizeof(double), A.nb, cudaMemcpyDeviceToHost,
            plasma->cuda_stream[rank]));
        CHECK_CUDA(cudaStreamSynchronize(plasma->cuda_stream[rank]));
        ss_cond_set(m, k, 1);
      } else {
        CHECK_CUDA(cudaMemcpy2DAsync(
            bufferA, ldam * sizeof(double), A(m, k), ldam * sizeof(double),
            tempmn * sizeof(double), A.nb, cudaMemcpyHostToDevice,
            plasma->cuda_stream[rank]));
        ss_cond_wait(k, n, 1);
        CHECK_CUDA(cudaMemcpy2DAsync(
            bufferB, ldak * sizeof(double), A(k, n), ldak * sizeof(double),
            A.nb * sizeof(double), A.nb, cudaMemcpyHostToDevice,
            plasma->cuda_stream[rank]));
        ss_cond_wait(m, n, 1);
        CHECK_CUDA(cudaMemcpy2DAsync(
            bufferC, ldam * sizeof(double), A(m, n), ldam * sizeof(double),
            tempmn * sizeof(double), A.nb, cudaMemcpyHostToDevice,
            plasma->cuda_stream[rank]));
        /*
         *  PlasmaLower
         */
        if (uplo == PlasmaLower) {
          CHECK_CUBLAS(cublasDgemm(plasma->cublas_handle[rank], CUBLAS_OP_N,
                                   CUBLAS_OP_T, tempmn, A.nb, A.nb, &mzone,
                                   bufferC, ldam, bufferB, ldak, &zone, bufferA,
                                   ldam));
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
        CHECK_CUDA(cudaMemcpy2DAsync(
            A(m, k), ldam * sizeof(double), bufferA, ldam * sizeof(double),
            tempmn * sizeof(double), A.nb, cudaMemcpyDeviceToHost,
            plasma->cuda_stream[rank]));
        CHECK_CUDA(cudaStreamSynchronize(plasma->cuda_stream[rank]));
      }
    }
    n = next_n;
    m = next_m;
    k = next_k;
  }
  CHECK_CUDA(cudaFree(pBuffer));
  CHECK_CUDA(cudaFree(buffer));
  ss_finalize();
}

/***************************************************************************/
/**
 *  Parallel tile Cholesky factorization - dynamic scheduling
 **/
void plasma_pdpotrf_gpu_quark(PLASMA_enum uplo, PLASMA_desc A,
                              PLASMA_sequence *sequence,
                              PLASMA_request *request) {
  plasma_context_t *plasma;
  Quark_Task_Flags task_flags = Quark_Task_Flags_Initializer;

  int k, m, n;
  int ldak, ldam, ldan;
  int tempkm, tempmm;

  double zone = (double)1.0;
  double mzone = (double)-1.0;

  plasma = plasma_context_self();
  if (sequence->status != PLASMA_SUCCESS) return;
  QUARK_Task_Flag_Set(&task_flags, TASK_SEQUENCE,
                      (intptr_t)sequence->quark_sequence);
  plasma_profile_by_function(&task_flags, POTRF);

  /*
   *  PlasmaLower
   */
  if (uplo == PlasmaLower) {
    for (k = 0; k < A.mt; k++) {
      tempkm = k == A.mt - 1 ? A.m - k * A.mb : A.mb;
      ldak = BLKLDD(A, k);

      QUARK_CORE_dpotrf(plasma->quark, &task_flags, PlasmaLower, tempkm, A.mb,
                        A(k, k), ldak, sequence, request, A.nb * k);

      for (m = k + 1; m < A.mt; m++) {
        tempmm = m == A.mt - 1 ? A.m - m * A.mb : A.mb;
        ldam = BLKLDD(A, m);
        QUARK_CORE_dtrsm(plasma->quark, &task_flags, PlasmaRight, PlasmaLower,
                         PlasmaTrans, PlasmaNonUnit, tempmm, A.mb, A.mb, zone,
                         A(k, k), ldak, A(m, k), ldam);
      }
      for (m = k + 1; m < A.mt; m++) {
        tempmm = m == A.mt - 1 ? A.m - m * A.mb : A.mb;
        ldam = BLKLDD(A, m);
        QUARK_CORE_dsyrk(plasma->quark, &task_flags, PlasmaLower, PlasmaNoTrans,
                         tempmm, A.mb, A.mb, -1.0, A(m, k), ldam, 1.0, A(m, m),
                         ldam);

        for (n = k + 1; n < m; n++) {
          ldan = BLKLDD(A, n);
          QUARK_CORE_dgemm(plasma->quark, &task_flags, PlasmaNoTrans,
                           PlasmaTrans, tempmm, A.mb, A.mb, A.mb, mzone,
                           A(m, k), ldam, A(n, k), ldan, zone, A(m, n), ldam);
        }
      }
    }
  }
  /*
   *  PlasmaUpper
   */
  else {
    for (k = 0; k < A.nt; k++) {
      tempkm = k == A.nt - 1 ? A.n - k * A.nb : A.nb;
      ldak = BLKLDD(A, k);
      QUARK_CORE_dpotrf(plasma->quark, &task_flags, PlasmaUpper, tempkm, A.mb,
                        A(k, k), ldak, sequence, request, A.nb * k);

      for (m = k + 1; m < A.nt; m++) {
        tempmm = m == A.nt - 1 ? A.n - m * A.nb : A.nb;
        QUARK_CORE_dtrsm(plasma->quark, &task_flags, PlasmaLeft, PlasmaUpper,
                         PlasmaTrans, PlasmaNonUnit, A.nb, tempmm, A.mb, zone,
                         A(k, k), ldak, A(k, m), ldak);
      }
      for (m = k + 1; m < A.nt; m++) {
        tempmm = m == A.nt - 1 ? A.n - m * A.nb : A.nb;
        ldam = BLKLDD(A, m);
        QUARK_CORE_dsyrk(plasma->quark, &task_flags, PlasmaUpper, PlasmaTrans,
                         tempmm, A.mb, A.mb, -1.0, A(k, m), ldak, 1.0, A(m, m),
                         ldam);

        for (n = k + 1; n < m; n++) {
          ldan = BLKLDD(A, n);
          QUARK_CORE_dgemm(plasma->quark, &task_flags, PlasmaTrans,
                           PlasmaNoTrans, A.mb, tempmm, A.mb, A.mb, mzone,
                           A(k, n), ldak, A(k, m), ldak, zone, A(n, m), ldan);
        }
      }
    }
  }
}

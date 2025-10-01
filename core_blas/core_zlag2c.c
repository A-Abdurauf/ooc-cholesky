/**
 *
 * @file core_zlag2c.c
 *
 *  PLASMA core_blas kernel
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @precisions mixed zc -> ds
 *
 **/
#ifdef PLASMA_WITH_MKL
#include <mkl_lapacke.h>
#else
#include <lapacke.h>
#endif
#include "common.h"

/***************************************************************************//**
 *
 * @ingroup CORE_PLASMA_Complex64_t
 *
 *  CORE_zlag2c converts a PLASMA_Complex64_t matrix, A, to a
 *  PLASMA_Complex32_t matrix, B.
 *
 *******************************************************************************
 *
 * @param[in] m
 *          The number of rows of the matrices A and B. m >= 0.
 *
 * @param[in] n
 *          The number of columns of the matrices A and B. n >= 0.
 *
 * @param[in] A
 *          The PLASMA_Complex64_t m-by-n matrix to convert.
 *
 * @param[in] lda
 *          The leading dimension of the array A. lda >= max(1,m).
 *
 * @param[out] B
 *          The PLASMA_Complex32_t m-by-n matrix to convert.
 *
 * @param[in] ldb
 *          The leading dimension of the array B. ldb >= max(1,m).
 *
 * @param[out] info
 *         - 0 on successful exit.
 *         - 1 if an entry of the matrix A is greater than the SINGLE
 *         PRECISION overflow threshold, in this case, the content
 *         of B in exit is unspecified.
 *
 ******************************************************************************/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_zlag2c = PCORE_zlag2c
#define CORE_zlag2c PCORE_zlag2c
#endif
void CORE_zlag2c(int m, int n,
                 const PLASMA_Complex64_t *A, int lda,
                 PLASMA_Complex32_t *B, int ldb, int *info)
{
    *info = LAPACKE_zlag2c_work(LAPACK_COL_MAJOR, m, n, A, lda, B, ldb);
}

/***************************************************************************//**
 *
 * @ingroup CORE_PLASMA_Complex64_t
 *
 *  CORE_PLASMA_clag2z converts a PLASMA_Complex32_t matrix, A, to a
 *  PLASMA_Complex64_t matrix, B.
 *
 *******************************************************************************
 *
 * @param[in] m
 *          The number of rows of the matrices A and B. m >= 0.
 *
 * @param[in] n
 *          The number of columns of the matrices A and B. n >= 0.
 *
 * @param[in] A
 *          The PLASMA_Complex32_t m-by-n matrix to convert.
 *
 * @param[in] lda
 *          The leading dimension of the array A. lda >= max(1,m).
 *
 * @param[out] B
 *          The PLASMA_Complex64_t m-by-n matrix to convert.
 *
 * @param[in] ldb
 *          The leading dimension of the array B. ldb >= max(1,m).
 *
 ******************************************************************************/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_clag2z = PCORE_clag2z
#define CORE_clag2z PCORE_clag2z
#endif
void CORE_clag2z(int m, int n,
                 const PLASMA_Complex32_t *A, int lda,
                 PLASMA_Complex64_t *B, int ldb)
{
    LAPACKE_clag2z_work(LAPACK_COL_MAJOR, m, n, A, lda, B, ldb);
}

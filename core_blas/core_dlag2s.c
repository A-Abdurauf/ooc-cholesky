/**
 *
 * @file core_dlag2s.c
 *
 *  PLASMA core_blas kernel
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @generated ds Fri Apr  1 11:02:29 2016
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
 * @ingroup CORE_double
 *
 *  CORE_dlag2s converts a double matrix, A, to a
 *  float matrix, B.
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
 *          The double m-by-n matrix to convert.
 *
 * @param[in] lda
 *          The leading dimension of the array A. lda >= max(1,m).
 *
 * @param[out] B
 *          The float m-by-n matrix to convert.
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
#pragma weak CORE_dlag2s = PCORE_dlag2s
#define CORE_dlag2s PCORE_dlag2s
#endif
void CORE_dlag2s(int m, int n,
                 const double *A, int lda,
                 float *B, int ldb, int *info)
{
    *info = LAPACKE_dlag2s_work(LAPACK_COL_MAJOR, m, n, A, lda, B, ldb);
}

/***************************************************************************//**
 *
 * @ingroup CORE_double
 *
 *  CORE_PLASMA_slag2d converts a float matrix, A, to a
 *  double matrix, B.
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
 *          The float m-by-n matrix to convert.
 *
 * @param[in] lda
 *          The leading dimension of the array A. lda >= max(1,m).
 *
 * @param[out] B
 *          The double m-by-n matrix to convert.
 *
 * @param[in] ldb
 *          The leading dimension of the array B. ldb >= max(1,m).
 *
 ******************************************************************************/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_slag2d = PCORE_slag2d
#define CORE_slag2d PCORE_slag2d
#endif
void CORE_slag2d(int m, int n,
                 const float *A, int lda,
                 double *B, int ldb)
{
    LAPACKE_slag2d_work(LAPACK_COL_MAJOR, m, n, A, lda, B, ldb);
}

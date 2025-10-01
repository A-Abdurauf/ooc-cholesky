/**
 *
 * @file qwrapper_zgemv.c
 *
 *  PLASMA core_blas quark wrapper
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Mark Gates
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @precisions normal z -> c d s
 *
 **/
#include <cblas.h>
#include "common.h"

/***************************************************************************//**
 *
 **/
void QUARK_CORE_zgemv(Quark *quark, Quark_Task_Flags *task_flags,
                      PLASMA_enum trans, int m, int n,
                      PLASMA_Complex64_t alpha, const PLASMA_Complex64_t *A, int lda,
                                                const PLASMA_Complex64_t *x, int incx,
                      PLASMA_Complex64_t beta,        PLASMA_Complex64_t *y, int incy)
{
    plasma_profile_by_kernel( task_flags, GEMV );

    QUARK_Insert_Task(quark, CORE_zgemv_quark, task_flags,
        sizeof(PLASMA_enum),             &trans,  VALUE,
        sizeof(int),                     &m,      VALUE,
        sizeof(int),                     &n,      VALUE,
        sizeof(PLASMA_Complex64_t),      &alpha,  VALUE,
        sizeof(PLASMA_Complex64_t)*m*n,  A,               INPUT,
        sizeof(int),                     &lda,    VALUE,
        sizeof(PLASMA_Complex64_t)*n,    x,               INPUT,
        sizeof(int),                     &incx,   VALUE,
        sizeof(PLASMA_Complex64_t),      &beta,   VALUE,
        sizeof(PLASMA_Complex64_t)*m,    y,               INOUT,
        sizeof(int),                     &incy,   VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_zgemv_quark = PCORE_zgemv_quark
#define CORE_zgemv_quark PCORE_zgemv_quark
#endif
void CORE_zgemv_quark(Quark *quark)
{
    PLASMA_enum trans;
    int m, n, lda, incx, incy;
    PLASMA_Complex64_t alpha, beta;
    const PLASMA_Complex64_t *A, *x;
    PLASMA_Complex64_t *y;

    quark_unpack_args_11( quark, trans, m, n, alpha, A, lda, x, incx, beta, y, incy );
    cblas_zgemv(
        CblasColMajor,
        (CBLAS_TRANSPOSE)trans,
        m, n,
        CBLAS_SADDR(alpha), A, lda,
                            x, incx,
        CBLAS_SADDR(beta),  y, incy);
}

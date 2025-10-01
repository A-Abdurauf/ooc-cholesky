/**
 *
 * @file qwrapper_sgemv.c
 *
 *  PLASMA core_blas quark wrapper
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Mark Gates
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @generated s Fri Apr  1 11:02:43 2016
 *
 **/
#include <cblas.h>
#include "common.h"

/***************************************************************************//**
 *
 **/
void QUARK_CORE_sgemv(Quark *quark, Quark_Task_Flags *task_flags,
                      PLASMA_enum trans, int m, int n,
                      float alpha, const float *A, int lda,
                                                const float *x, int incx,
                      float beta,        float *y, int incy)
{
    plasma_profile_by_kernel( task_flags, GEMV );

    QUARK_Insert_Task(quark, CORE_sgemv_quark, task_flags,
        sizeof(PLASMA_enum),             &trans,  VALUE,
        sizeof(int),                     &m,      VALUE,
        sizeof(int),                     &n,      VALUE,
        sizeof(float),      &alpha,  VALUE,
        sizeof(float)*m*n,  A,               INPUT,
        sizeof(int),                     &lda,    VALUE,
        sizeof(float)*n,    x,               INPUT,
        sizeof(int),                     &incx,   VALUE,
        sizeof(float),      &beta,   VALUE,
        sizeof(float)*m,    y,               INOUT,
        sizeof(int),                     &incy,   VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_sgemv_quark = PCORE_sgemv_quark
#define CORE_sgemv_quark PCORE_sgemv_quark
#endif
void CORE_sgemv_quark(Quark *quark)
{
    PLASMA_enum trans;
    int m, n, lda, incx, incy;
    float alpha, beta;
    const float *A, *x;
    float *y;

    quark_unpack_args_11( quark, trans, m, n, alpha, A, lda, x, incx, beta, y, incy );
    cblas_sgemv(
        CblasColMajor,
        (CBLAS_TRANSPOSE)trans,
        m, n,
        (alpha), A, lda,
                            x, incx,
        (beta),  y, incy);
}

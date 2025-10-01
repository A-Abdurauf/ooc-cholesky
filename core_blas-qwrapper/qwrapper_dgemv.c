/**
 *
 * @file qwrapper_dgemv.c
 *
 *  PLASMA core_blas quark wrapper
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Mark Gates
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @generated d Fri Apr  1 11:02:43 2016
 *
 **/
#include <cblas.h>
#include "common.h"

/***************************************************************************//**
 *
 **/
void QUARK_CORE_dgemv(Quark *quark, Quark_Task_Flags *task_flags,
                      PLASMA_enum trans, int m, int n,
                      double alpha, const double *A, int lda,
                                                const double *x, int incx,
                      double beta,        double *y, int incy)
{
    plasma_profile_by_kernel( task_flags, GEMV );

    QUARK_Insert_Task(quark, CORE_dgemv_quark, task_flags,
        sizeof(PLASMA_enum),             &trans,  VALUE,
        sizeof(int),                     &m,      VALUE,
        sizeof(int),                     &n,      VALUE,
        sizeof(double),      &alpha,  VALUE,
        sizeof(double)*m*n,  A,               INPUT,
        sizeof(int),                     &lda,    VALUE,
        sizeof(double)*n,    x,               INPUT,
        sizeof(int),                     &incx,   VALUE,
        sizeof(double),      &beta,   VALUE,
        sizeof(double)*m,    y,               INOUT,
        sizeof(int),                     &incy,   VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_dgemv_quark = PCORE_dgemv_quark
#define CORE_dgemv_quark PCORE_dgemv_quark
#endif
void CORE_dgemv_quark(Quark *quark)
{
    PLASMA_enum trans;
    int m, n, lda, incx, incy;
    double alpha, beta;
    const double *A, *x;
    double *y;

    quark_unpack_args_11( quark, trans, m, n, alpha, A, lda, x, incx, beta, y, incy );
    cblas_dgemv(
        CblasColMajor,
        (CBLAS_TRANSPOSE)trans,
        m, n,
        (alpha), A, lda,
                            x, incx,
        (beta),  y, incy);
}

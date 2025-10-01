/**
 *
 * @file qwrapper_dgemv_tile.c
 *
 *  PLASMA core_blas quark wrapper
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Mark Gates
 * @date 2010-11-15
 * @generated d Fri Apr  1 11:02:43 2016
 *
 **/
#include <cblas.h>
#include "common.h"

/***************************************************************************//**
 *
 * Version of zgemv for tile storage, to avoid dependency problem when
 * computations are done within the tile. alpha and beta are passed as
 * pointers so they can depend on runtime values.
 *
 * @param[in] Alock
 *          Pointer to tile owning submatrix A.
 *
 * @param[in] xlock
 *          Pointer to tile owning subvector x.
 *
 * @param[in] ylock
 *          Pointer to tile owning subvector y.
 *
 **/
void QUARK_CORE_dgemv_tile(Quark *quark, Quark_Task_Flags *task_flags,
                           PLASMA_enum trans,
                           int m, int n,
                           const double *alpha, const double *A, int lda,
                                                            const double *x, int incx,
                           const double *beta,        double *y, int incy,
                           const double *Alock,
                           const double *xlock,
                           const double *ylock)
{
    /* Quick return. Bad things happen if sizeof(...)*m*n is zero in QUARK_Insert_Task */
    if ( m == 0 || n == 0 )
        return;

    DAG_SET_PROPERTIES("gemv", "lightslateblue");
    QUARK_Insert_Task(quark, CORE_dgemv_tile_quark, task_flags,
        sizeof(PLASMA_enum),             &trans,  VALUE,
        sizeof(int),                     &m,      VALUE,
        sizeof(int),                     &n,      VALUE,
        sizeof(double),      alpha,           INPUT,
        sizeof(double)*m*n,  A,               NODEP,          /* input; see Alock */
        sizeof(int),                     &lda,    VALUE,
        sizeof(double)*n,    x,               NODEP,          /* input; see xlock */
        sizeof(int),                     &incx,   VALUE,
        sizeof(double),      beta,            INPUT,
        sizeof(double)*m,    y,                       NODEP,  /* inout; see ylock */
        sizeof(int),                     &incy,   VALUE,
        sizeof(double)*m*n,  Alock,           INPUT,
        sizeof(double)*n,    xlock,           INPUT,
        sizeof(double)*m,    ylock,                   INOUT,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_dgemv_tile_quark = PCORE_dgemv_tile_quark
#define CORE_dgemv_tile_quark PCORE_dgemv_tile_quark
#endif
void CORE_dgemv_tile_quark(Quark *quark)
{
    PLASMA_enum trans;
    int m, n, lda, incx, incy;
    const double *alpha, *beta;
    const double *A, *x;
    double *y;

    quark_unpack_args_11( quark, trans, m, n, alpha, A, lda, x, incx, beta, y, incy );
    cblas_dgemv(
        CblasColMajor,
        (CBLAS_TRANSPOSE)trans,
        m, n,
        (*alpha), A, lda,
                             x, incx,
        (*beta),  y, incy );
}

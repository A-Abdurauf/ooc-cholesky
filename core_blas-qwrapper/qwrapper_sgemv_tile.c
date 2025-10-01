/**
 *
 * @file qwrapper_sgemv_tile.c
 *
 *  PLASMA core_blas quark wrapper
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Mark Gates
 * @date 2010-11-15
 * @generated s Fri Apr  1 11:02:43 2016
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
void QUARK_CORE_sgemv_tile(Quark *quark, Quark_Task_Flags *task_flags,
                           PLASMA_enum trans,
                           int m, int n,
                           const float *alpha, const float *A, int lda,
                                                            const float *x, int incx,
                           const float *beta,        float *y, int incy,
                           const float *Alock,
                           const float *xlock,
                           const float *ylock)
{
    /* Quick return. Bad things happen if sizeof(...)*m*n is zero in QUARK_Insert_Task */
    if ( m == 0 || n == 0 )
        return;

    DAG_SET_PROPERTIES("gemv", "lightslateblue");
    QUARK_Insert_Task(quark, CORE_sgemv_tile_quark, task_flags,
        sizeof(PLASMA_enum),             &trans,  VALUE,
        sizeof(int),                     &m,      VALUE,
        sizeof(int),                     &n,      VALUE,
        sizeof(float),      alpha,           INPUT,
        sizeof(float)*m*n,  A,               NODEP,          /* input; see Alock */
        sizeof(int),                     &lda,    VALUE,
        sizeof(float)*n,    x,               NODEP,          /* input; see xlock */
        sizeof(int),                     &incx,   VALUE,
        sizeof(float),      beta,            INPUT,
        sizeof(float)*m,    y,                       NODEP,  /* inout; see ylock */
        sizeof(int),                     &incy,   VALUE,
        sizeof(float)*m*n,  Alock,           INPUT,
        sizeof(float)*n,    xlock,           INPUT,
        sizeof(float)*m,    ylock,                   INOUT,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_sgemv_tile_quark = PCORE_sgemv_tile_quark
#define CORE_sgemv_tile_quark PCORE_sgemv_tile_quark
#endif
void CORE_sgemv_tile_quark(Quark *quark)
{
    PLASMA_enum trans;
    int m, n, lda, incx, incy;
    const float *alpha, *beta;
    const float *A, *x;
    float *y;

    quark_unpack_args_11( quark, trans, m, n, alpha, A, lda, x, incx, beta, y, incy );
    cblas_sgemv(
        CblasColMajor,
        (CBLAS_TRANSPOSE)trans,
        m, n,
        (*alpha), A, lda,
                             x, incx,
        (*beta),  y, incy );
}

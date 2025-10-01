/**
 *
 * @file qwrapper_dpltmg_fiedler.c
 *
 *  PLASMA core_blas quark wrapper
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @generated d Fri Apr  1 11:02:42 2016
 *
 **/
#include "common.h"

void
CORE_dpltmg_fiedler_quark(Quark *quark);

/***************************************************************************//**
 *
 **/
void QUARK_CORE_dpltmg_fiedler(Quark *quark, Quark_Task_Flags *task_flags,
                               int m, int n,
                               const double *X, int incX,
                               const double *Y, int incY,
                                     double *A, int lda)
{
    plasma_profile_by_kernel( task_flags, PLRNT );


    if ( X == Y ) {
        QUARK_Insert_Task(quark, CORE_dpltmg_fiedler_quark, task_flags,
            sizeof(int),                     &m,      VALUE,
            sizeof(int),                     &n,      VALUE,
            sizeof(double)*m,    X,               INPUT,
            sizeof(int),                     &incX,   VALUE,
            sizeof(double*),     &Y,      VALUE,
            sizeof(int),                     &incY,   VALUE,
            sizeof(double)*m*n,  A,               OUTPUT,
            sizeof(int),                     &lda,    VALUE,
            0);
    } else {
        QUARK_Insert_Task(quark, CORE_dpltmg_fiedler_quark, task_flags,
            sizeof(int),                     &m,      VALUE,
            sizeof(int),                     &n,      VALUE,
            sizeof(double)*m,    X,               INPUT,
            sizeof(int),                     &incX,   VALUE,
            sizeof(double)*n,    Y,               INPUT,
            sizeof(int),                     &incY,   VALUE,
            sizeof(double)*m*n,  A,               OUTPUT,
            sizeof(int),                     &lda,    VALUE,
            0);
    }
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_dpltmg_fiedler_quark = PCORE_dpltmg_fiedler_quark
#define CORE_dpltmg_fiedler_quark PCORE_dpltmg_fiedler_quark
#endif
void CORE_dpltmg_fiedler_quark(Quark *quark)
{
    int m, n, lda, incx, incy;
    const double *X, *Y;
    double *A;

    quark_unpack_args_8( quark, m, n, X, incx, Y, incy, A, lda );
    CORE_dpltmg_fiedler( m, n, X, incx, Y, incy, A, lda );
}

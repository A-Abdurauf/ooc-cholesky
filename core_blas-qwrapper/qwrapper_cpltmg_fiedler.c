/**
 *
 * @file qwrapper_cpltmg_fiedler.c
 *
 *  PLASMA core_blas quark wrapper
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @generated c Fri Apr  1 11:02:42 2016
 *
 **/
#include "common.h"

void
CORE_cpltmg_fiedler_quark(Quark *quark);

/***************************************************************************//**
 *
 **/
void QUARK_CORE_cpltmg_fiedler(Quark *quark, Quark_Task_Flags *task_flags,
                               int m, int n,
                               const PLASMA_Complex32_t *X, int incX,
                               const PLASMA_Complex32_t *Y, int incY,
                                     PLASMA_Complex32_t *A, int lda)
{
    plasma_profile_by_kernel( task_flags, PLRNT );


    if ( X == Y ) {
        QUARK_Insert_Task(quark, CORE_cpltmg_fiedler_quark, task_flags,
            sizeof(int),                     &m,      VALUE,
            sizeof(int),                     &n,      VALUE,
            sizeof(PLASMA_Complex32_t)*m,    X,               INPUT,
            sizeof(int),                     &incX,   VALUE,
            sizeof(PLASMA_Complex32_t*),     &Y,      VALUE,
            sizeof(int),                     &incY,   VALUE,
            sizeof(PLASMA_Complex32_t)*m*n,  A,               OUTPUT,
            sizeof(int),                     &lda,    VALUE,
            0);
    } else {
        QUARK_Insert_Task(quark, CORE_cpltmg_fiedler_quark, task_flags,
            sizeof(int),                     &m,      VALUE,
            sizeof(int),                     &n,      VALUE,
            sizeof(PLASMA_Complex32_t)*m,    X,               INPUT,
            sizeof(int),                     &incX,   VALUE,
            sizeof(PLASMA_Complex32_t)*n,    Y,               INPUT,
            sizeof(int),                     &incY,   VALUE,
            sizeof(PLASMA_Complex32_t)*m*n,  A,               OUTPUT,
            sizeof(int),                     &lda,    VALUE,
            0);
    }
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_cpltmg_fiedler_quark = PCORE_cpltmg_fiedler_quark
#define CORE_cpltmg_fiedler_quark PCORE_cpltmg_fiedler_quark
#endif
void CORE_cpltmg_fiedler_quark(Quark *quark)
{
    int m, n, lda, incx, incy;
    const PLASMA_Complex32_t *X, *Y;
    PLASMA_Complex32_t *A;

    quark_unpack_args_8( quark, m, n, X, incx, Y, incy, A, lda );
    CORE_cpltmg_fiedler( m, n, X, incx, Y, incy, A, lda );
}

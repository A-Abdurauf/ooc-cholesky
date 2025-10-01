/**
 *
 * @file qwrapper_zpltmg_fiedler.c
 *
 *  PLASMA core_blas quark wrapper
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @precisions normal z -> c d s
 *
 **/
#include "common.h"

void
CORE_zpltmg_fiedler_quark(Quark *quark);

/***************************************************************************//**
 *
 **/
void QUARK_CORE_zpltmg_fiedler(Quark *quark, Quark_Task_Flags *task_flags,
                               int m, int n,
                               const PLASMA_Complex64_t *X, int incX,
                               const PLASMA_Complex64_t *Y, int incY,
                                     PLASMA_Complex64_t *A, int lda)
{
    plasma_profile_by_kernel( task_flags, PLRNT );


    if ( X == Y ) {
        QUARK_Insert_Task(quark, CORE_zpltmg_fiedler_quark, task_flags,
            sizeof(int),                     &m,      VALUE,
            sizeof(int),                     &n,      VALUE,
            sizeof(PLASMA_Complex64_t)*m,    X,               INPUT,
            sizeof(int),                     &incX,   VALUE,
            sizeof(PLASMA_Complex64_t*),     &Y,      VALUE,
            sizeof(int),                     &incY,   VALUE,
            sizeof(PLASMA_Complex64_t)*m*n,  A,               OUTPUT,
            sizeof(int),                     &lda,    VALUE,
            0);
    } else {
        QUARK_Insert_Task(quark, CORE_zpltmg_fiedler_quark, task_flags,
            sizeof(int),                     &m,      VALUE,
            sizeof(int),                     &n,      VALUE,
            sizeof(PLASMA_Complex64_t)*m,    X,               INPUT,
            sizeof(int),                     &incX,   VALUE,
            sizeof(PLASMA_Complex64_t)*n,    Y,               INPUT,
            sizeof(int),                     &incY,   VALUE,
            sizeof(PLASMA_Complex64_t)*m*n,  A,               OUTPUT,
            sizeof(int),                     &lda,    VALUE,
            0);
    }
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_zpltmg_fiedler_quark = PCORE_zpltmg_fiedler_quark
#define CORE_zpltmg_fiedler_quark PCORE_zpltmg_fiedler_quark
#endif
void CORE_zpltmg_fiedler_quark(Quark *quark)
{
    int m, n, lda, incx, incy;
    const PLASMA_Complex64_t *X, *Y;
    PLASMA_Complex64_t *A;

    quark_unpack_args_8( quark, m, n, X, incx, Y, incy, A, lda );
    CORE_zpltmg_fiedler( m, n, X, incx, Y, incy, A, lda );
}

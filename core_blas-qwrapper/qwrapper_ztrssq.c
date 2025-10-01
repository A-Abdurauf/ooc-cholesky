/**
 *
 * @file qwrapper_ztrssq.c
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
CORE_ztrssq_quark(Quark *quark);
void
CORE_ztrssq_f1_quark(Quark *quark);

/***************************************************************************//**
 *
 **/
void QUARK_CORE_ztrssq_f1( Quark *quark, Quark_Task_Flags *task_flags,
                           PLASMA_enum uplo, PLASMA_enum diag,
                           int m, int n, const PLASMA_Complex64_t *A, int lda,
                           double *scale, double *sumsq,
                           double *fake, int szeF, int paramF )
{
    plasma_profile_by_kernel( task_flags, LASSQ );

    if ( (fake == scale) && (paramF & GATHERV) ) {
        QUARK_Insert_Task(quark, CORE_ztrssq_quark, task_flags,
            sizeof(PLASMA_enum),              &uplo, VALUE,
            sizeof(PLASMA_enum),              &diag, VALUE,
            sizeof(int),                      &m,    VALUE,
            sizeof(int),                      &n,    VALUE,
            sizeof(PLASMA_Complex64_t)*lda*n, A,         INPUT,
            sizeof(int),                      &lda,  VALUE,
            sizeof(double)*1,                 scale,     INOUT | GATHERV,
            sizeof(double)*1,                 sumsq,     INOUT,
            0);
    } else {
        QUARK_Insert_Task(quark, CORE_ztrssq_f1_quark, task_flags,
            sizeof(PLASMA_enum),              &uplo, VALUE,
            sizeof(PLASMA_enum),              &diag, VALUE,
            sizeof(int),                      &m,    VALUE,
            sizeof(int),                      &n,    VALUE,
            sizeof(PLASMA_Complex64_t)*lda*n, A,         INPUT,
            sizeof(int),                      &lda,  VALUE,
            sizeof(double)*1,                 scale,     INOUT,
            sizeof(double)*1,                 sumsq,     INOUT,
            sizeof(double)*szeF,              fake,      paramF,
            0);
    }
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_ztrssq_quark = PCORE_ztrssq_quark
#define CORE_ztrssq_quark PCORE_ztrssq_quark
#endif
void CORE_ztrssq_quark(Quark *quark)
{
    PLASMA_enum uplo, diag;
    int m;
    int n;
    PLASMA_Complex64_t *A;
    int lda;
    double *scale;
    double *sumsq;

    quark_unpack_args_8( quark, uplo, diag, m, n, A, lda, scale, sumsq );
    CORE_ztrssq( uplo, diag, m, n, A, lda, scale, sumsq );
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_ztrssq_f1_quark = PCORE_ztrssq_f1_quark
#define CORE_ztrssq_f1_quark PCORE_ztrssq_f1_quark
#endif
void CORE_ztrssq_f1_quark(Quark *quark)
{
    PLASMA_enum uplo, diag;
    int m;
    int n;
    PLASMA_Complex64_t *A;
    int lda;
    double *scale;
    double *sumsq;
    double *fake;

    quark_unpack_args_9( quark, uplo, diag, m, n, A, lda, scale, sumsq, fake );
    CORE_ztrssq( uplo, diag, m, n, A, lda, scale, sumsq );
}

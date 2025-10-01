/**
 *
 * @file qwrapper_dsyssq.c
 *
 *  PLASMA core_blas quark wrapper
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @generated d Fri Apr  1 11:02:41 2016
 *
 **/
#include "common.h"

/***************************************************************************//**
 *
 **/
void QUARK_CORE_dsyssq_f1( Quark *quark, Quark_Task_Flags *task_flags,
                           PLASMA_enum uplo, int n, const double *A, int lda,
                           double *scale, double *sumsq,
                           double *fake, int szeF, int paramF )
{
    plasma_profile_by_kernel( task_flags, LASSQ );


    if ( (fake == scale) && (paramF & GATHERV) ) {
        QUARK_Insert_Task(quark, CORE_dsyssq_quark, task_flags,
            sizeof(PLASMA_enum),              &uplo, VALUE,
            sizeof(int),                      &n,    VALUE,
            sizeof(double)*lda*n, A,         INPUT,
            sizeof(int),                      &lda,  VALUE,
            sizeof(double)*1,                 scale,     INOUT | paramF,
            sizeof(double)*1,                 sumsq,     INOUT,
            0);
    } else {
        QUARK_Insert_Task(quark, CORE_dsyssq_f1_quark, task_flags,
            sizeof(PLASMA_enum),              &uplo, VALUE,
            sizeof(int),                      &n,    VALUE,
            sizeof(double)*lda*n, A,         INPUT,
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
#pragma weak CORE_dsyssq_quark = PCORE_dsyssq_quark
#define CORE_dsyssq_quark PCORE_dsyssq_quark
#endif
void CORE_dsyssq_quark(Quark *quark)
{
    PLASMA_enum uplo;
    int n;
    double *A;
    int lda;
    double *scale;
    double *sumsq;

    quark_unpack_args_6( quark, uplo, n, A, lda, scale, sumsq );
    CORE_dsyssq( uplo, n, A, lda, scale, sumsq );
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_dsyssq_f1_quark = PCORE_dsyssq_f1_quark
#define CORE_dsyssq_f1_quark PCORE_dsyssq_f1_quark
#endif
void CORE_dsyssq_f1_quark(Quark *quark)
{
    PLASMA_enum uplo;
    int n;
    double *A;
    int lda;
    double *scale;
    double *sumsq;
    double *fake;

    quark_unpack_args_7( quark, uplo, n, A, lda, scale, sumsq, fake );
    CORE_dsyssq( uplo, n, A, lda, scale, sumsq );
}

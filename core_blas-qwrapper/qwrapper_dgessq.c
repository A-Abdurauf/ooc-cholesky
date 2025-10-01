/**
 *
 * @file qwrapper_dgessq.c
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
void QUARK_CORE_dgessq_f1( Quark *quark, Quark_Task_Flags *task_flags,
                           int m, int n, const double *A, int lda,
                           double *scale, double *sumsq,
                           double *fake, int szeF, int paramF )
{
    plasma_profile_by_kernel( task_flags, LASSQ );

    if ( (fake == scale) && (paramF & GATHERV) ) {
        QUARK_Insert_Task(quark, CORE_dgessq_quark, task_flags,
            sizeof(int),                      &m,    VALUE,
            sizeof(int),                      &n,    VALUE,
            sizeof(double)*lda*n, A,         INPUT,
            sizeof(int),                      &lda,  VALUE,
            sizeof(double)*1,                 scale,     INOUT | GATHERV,
            sizeof(double)*1,                 sumsq,     INOUT,
            0);
    } else {
        QUARK_Insert_Task(quark, CORE_dgessq_f1_quark, task_flags,
            sizeof(int),                      &m,    VALUE,
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
#pragma weak CORE_dgessq_quark = PCORE_dgessq_quark
#define CORE_dgessq_quark PCORE_dgessq_quark
#endif
void CORE_dgessq_quark(Quark *quark)
{
    int m;
    int n;
    double *A;
    int lda;
    double *scale;
    double *sumsq;

    quark_unpack_args_6( quark, m, n, A, lda, scale, sumsq );
    CORE_dgessq( m, n, A, lda, scale, sumsq );
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_dgessq_f1_quark = PCORE_dgessq_f1_quark
#define CORE_dgessq_f1_quark PCORE_dgessq_f1_quark
#endif
void CORE_dgessq_f1_quark(Quark *quark)
{
    int m;
    int n;
    double *A;
    int lda;
    double *scale;
    double *sumsq;
    double *fake;

    quark_unpack_args_7( quark, m, n, A, lda, scale, sumsq, fake );
    CORE_dgessq( m, n, A, lda, scale, sumsq );
}

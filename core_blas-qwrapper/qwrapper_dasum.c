/**
 *
 * @file qwrapper_dasum.c
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
void QUARK_CORE_dasum(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum storev, PLASMA_enum uplo, int M, int N,
                       const double *A, int lda, int szeA,
                       double *work, int szeW)
{
    QUARK_Insert_Task(
        quark, CORE_dasum_quark, task_flags,
        sizeof(PLASMA_enum),                &storev,    VALUE,
        sizeof(PLASMA_enum),                &uplo,      VALUE,
        sizeof(int),                        &M,         VALUE,
        sizeof(int),                        &N,         VALUE,
        sizeof(double)*szeA,     A,                 INPUT,
        sizeof(int),                        &lda,       VALUE,
        sizeof(double)*szeW,                 work,              INOUT,
        0);
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_dasum_f1(Quark *quark, Quark_Task_Flags *task_flags,
                          PLASMA_enum storev, PLASMA_enum uplo, int M, int N,
                          const double *A, int lda, int szeA,
                          double *work, int szeW, double *fake, int szeF)
{
    plasma_profile_by_kernel( task_flags, ASUM );

    if ( work == fake ) {
        QUARK_Insert_Task(
            quark, CORE_dasum_quark, task_flags,
            sizeof(PLASMA_enum),                &storev,    VALUE,
            sizeof(PLASMA_enum),                &uplo,      VALUE,
            sizeof(int),                        &M,         VALUE,
            sizeof(int),                        &N,         VALUE,
            sizeof(double)*szeA,     A,                 INPUT,
            sizeof(int),                        &lda,       VALUE,
            sizeof(double)*szeW,                 work,              INOUT | GATHERV,
            0);
    } else {
        QUARK_Insert_Task(
            quark, CORE_dasum_f1_quark, task_flags,
            sizeof(PLASMA_enum),                &storev,    VALUE,
            sizeof(PLASMA_enum),                &uplo,      VALUE,
            sizeof(int),                        &M,         VALUE,
            sizeof(int),                        &N,         VALUE,
            sizeof(double)*szeA,     A,                 INPUT,
            sizeof(int),                        &lda,       VALUE,
            sizeof(double)*szeW,                 work,              INOUT,
            sizeof(double)*szeF,                 fake,              INOUT | GATHERV,
            0);
    }
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_dasum_quark = PCORE_dasum_quark
#define CORE_dasum_quark PCORE_dasum_quark
#endif
void CORE_dasum_quark(Quark *quark)
{
    PLASMA_enum storev;
    PLASMA_enum uplo;
    int M;
    int N;
    double *A;
    int lda;
    double *work;

    quark_unpack_args_7(quark, storev, uplo, M, N, A, lda, work);
    CORE_dasum(storev, uplo, M, N, A, lda, work);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_dasum_f1_quark = PCORE_dasum_f1_quark
#define CORE_dasum_f1_quark PCORE_dasum_f1_quark
#endif
void CORE_dasum_f1_quark(Quark *quark)
{
    PLASMA_enum storev;
    PLASMA_enum uplo;
    int M;
    int N;
    double *A;
    int lda;
    double *work;
    double *fake;

    quark_unpack_args_8(quark, storev, uplo, M, N, A, lda, work, fake);
    CORE_dasum(storev, uplo, M, N, A, lda, work);
}

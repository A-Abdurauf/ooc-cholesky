/**
 *
 * @file qwrapper_dlatro.c
 *
 *  PLASMA core_blas quark wrapper
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Azzam Haidar
 * @date 2010-11-15
 * @generated d Fri Apr  1 11:02:43 2016
 *
 **/
#include "common.h"

/***************************************************************************//**
 *
 **/
void QUARK_CORE_dlatro(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum uplo, PLASMA_enum trans, int m, int n, int nb,
                       const double *A, int lda,
                             double *B, int ldb)
{
    plasma_profile_by_kernel( task_flags, LATRO );

    QUARK_Insert_Task(quark, CORE_dlatro_quark, task_flags,
        sizeof(PLASMA_enum),                &uplo,  VALUE,
        sizeof(PLASMA_enum),                &trans, VALUE,
        sizeof(int),                        &m,     VALUE,
        sizeof(int),                        &n,     VALUE,
        sizeof(double)*nb*nb,    A,             INPUT,
        sizeof(int),                        &lda,   VALUE,
        sizeof(double)*nb*nb,    B,             OUTPUT,
        sizeof(int),                        &ldb,   VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_dlatro_quark = PCORE_dlatro_quark
#define CORE_dlatro_quark PCORE_dlatro_quark
#endif
void CORE_dlatro_quark(Quark *quark)
{
    PLASMA_enum uplo;
    PLASMA_enum trans;
    int M;
    int N;
    const double *A;
    int LDA;
    double *B;
    int LDB;

    quark_unpack_args_8(quark, uplo, trans, M, N, A, LDA, B, LDB);
    CORE_dlatro(uplo, trans, M, N, A, LDA, B, LDB);
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_dlatro_f1(Quark *quark, Quark_Task_Flags *task_flags,
                          PLASMA_enum uplo, PLASMA_enum trans, int m, int n, int nb,
                          const double *A, int lda,
                                double *B, int ldb,
                          double *fake1, int szefake1, int flag1)
{
    plasma_profile_by_kernel( task_flags, LATRO );

    if ( fake1 == B ) {
        QUARK_Insert_Task(quark, CORE_dlatro_quark, task_flags,
            sizeof(PLASMA_enum),                &uplo,  VALUE,
            sizeof(PLASMA_enum),                &trans, VALUE,
            sizeof(int),                        &m,     VALUE,
            sizeof(int),                        &n,     VALUE,
            sizeof(double)*nb*nb,    A,             INPUT,
            sizeof(int),                        &lda,   VALUE,
            sizeof(double)*nb*nb,    B,             OUTPUT,
            sizeof(int),                        &ldb,   VALUE,
            0);
    } else {
        QUARK_Insert_Task(quark, CORE_dlatro_f1_quark, task_flags,
            sizeof(PLASMA_enum),                &uplo,  VALUE,
            sizeof(PLASMA_enum),                &trans, VALUE,
            sizeof(int),                        &m,     VALUE,
            sizeof(int),                        &n,     VALUE,
            sizeof(double)*nb*nb,    A,             INPUT,
            sizeof(int),                        &lda,   VALUE,
            sizeof(double)*nb*nb,    B,             OUTPUT,
            sizeof(int),                        &ldb,   VALUE,
            sizeof(double)*szefake1, fake1,         flag1,
            0);
    }
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_dlatro_f1_quark = PCORE_dlatro_f1_quark
#define CORE_dlatro_f1_quark PCORE_dlatro_f1_quark
#endif
void CORE_dlatro_f1_quark(Quark *quark)
{
    PLASMA_enum uplo;
    PLASMA_enum trans;
    int M;
    int N;
    const double *A;
    int LDA;
    double *B;
    int LDB;
    void *fake1;

    quark_unpack_args_9(quark, uplo, trans, M, N, A, LDA, B, LDB, fake1);
    CORE_dlatro(uplo, trans, M, N, A, LDA, B, LDB);
}


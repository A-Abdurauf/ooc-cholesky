/**
 *
 * @file qwrapper_dlacpy.c
 *
 *  PLASMA core_blas quark wrapper
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Julien Langou
 * @author Henricus Bouwmeester
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @generated d Fri Apr  1 11:02:40 2016
 *
 **/
#ifdef PLASMA_WITH_MKL
#include <mkl_lapacke.h>
#else
#include <lapack.h>
#endif
#include "common.h"

/***************************************************************************//**
 *
 **/
void QUARK_CORE_dlacpy(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum uplo, int m, int n, int nb,
                       const double *A, int lda,
                       double *B, int ldb)
{
    plasma_profile_by_kernel( task_flags, LACPY );

    QUARK_Insert_Task(quark, CORE_dlacpy_quark, task_flags,
        sizeof(PLASMA_enum),                &uplo,  VALUE,
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
#pragma weak CORE_dlacpy_quark = PCORE_dlacpy_quark
#define CORE_dlacpy_quark PCORE_dlacpy_quark
#endif
void CORE_dlacpy_quark(Quark *quark)
{
    PLASMA_enum uplo;
    int M;
    int N;
    const double *A;
    int LDA;
    double *B;
    int LDB;

    quark_unpack_args_7(quark, uplo, M, N, A, LDA, B, LDB);
#ifdef PLASMA_WITH_MKL
    LAPACKE_dlacpy_work(
        LAPACK_COL_MAJOR,
        lapack_const(uplo),
        M, N, A, LDA, B, LDB);
#else
    LAPACK_dlacpy(&lapack_const(uplo), &M, &N, A, &LDA, B, &LDB);
#endif
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_dlacpy_f1(Quark *quark, Quark_Task_Flags *task_flags,
                          PLASMA_enum uplo, int m, int n, int nb,
                          const double *A, int lda,
                          double *B, int ldb,
                          double *fake1, int szefake1, int flag1)
{
    plasma_profile_by_kernel( task_flags, LACPY );

    if ( fake1 == B ) {
        QUARK_Insert_Task(quark, CORE_dlacpy_quark, task_flags,
            sizeof(PLASMA_enum),                &uplo,  VALUE,
            sizeof(int),                        &m,     VALUE,
            sizeof(int),                        &n,     VALUE,
            sizeof(double)*nb*nb,    A,             INPUT,
            sizeof(int),                        &lda,   VALUE,
            sizeof(double)*nb*nb,    B,             OUTPUT | flag1,
            sizeof(int),                        &ldb,   VALUE,
            0);
    }
    else {
        QUARK_Insert_Task(quark, CORE_dlacpy_f1_quark, task_flags,
            sizeof(PLASMA_enum),                &uplo,  VALUE,
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
#pragma weak CORE_dlacpy_f1_quark = PCORE_dlacpy_f1_quark
#define CORE_dlacpy_f1_quark PCORE_dlacpy_f1_quark
#endif
void CORE_dlacpy_f1_quark(Quark *quark)
{
    PLASMA_enum uplo;
    int M;
    int N;
    const double *A;
    int LDA;
    double *B;
    int LDB;
    void *fake1;

    quark_unpack_args_8(quark, uplo, M, N, A, LDA, B, LDB, fake1);
#ifdef PLASMA_WITH_MKL
    LAPACKE_dlacpy_work(
        LAPACK_COL_MAJOR,
        lapack_const(uplo),
        M, N, A, LDA, B, LDB);
#else
    LAPACK_dlacpy(&lapack_const(uplo), &M, &N, A, &LDA, B, &LDB);
#endif
}


/**
 *
 * @file qwrapper_dlantr.c
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
#ifdef PLASMA_WITH_MKL
#include <mkl_lapacke.h>
#else
#include <lapacke.h>
#endif
#include "common.h"

void
CORE_dlantr_quark(Quark *quark);
void
CORE_dlantr_f1_quark(Quark *quark);

/***************************************************************************//**
 *
 **/
void QUARK_CORE_dlantr(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum norm, PLASMA_enum uplo, PLASMA_enum diag, int M, int N,
                       const double *A, int LDA, int szeA,
                       int szeW, double *result)
{
    szeW = max(1, szeW);
    plasma_profile_by_kernel( task_flags, LANGE );

    QUARK_Insert_Task(quark, CORE_dlantr_quark, task_flags,
        sizeof(PLASMA_enum),                &norm,  VALUE,
        sizeof(PLASMA_enum),                &uplo,  VALUE,
        sizeof(PLASMA_enum),                &diag,  VALUE,
        sizeof(int),                        &M,     VALUE,
        sizeof(int),                        &N,     VALUE,
        sizeof(double)*szeA,     A,             INPUT,
        sizeof(int),                        &LDA,   VALUE,
        sizeof(double)*szeW,                 NULL,          SCRATCH,
        sizeof(double),                      result,        OUTPUT,
        0);
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_dlantr_f1(Quark *quark, Quark_Task_Flags *task_flags,
                          PLASMA_enum norm, PLASMA_enum uplo, PLASMA_enum diag, int M, int N,
                          const double *A, int LDA, int szeA,
                          int szeW, double *result,
                          double *fake, int szeF)
{
    szeW = max(1, szeW);
    plasma_profile_by_kernel( task_flags, LANGE );


    if ( result == fake ) {
        QUARK_Insert_Task(quark, CORE_dlantr_quark, task_flags,
            sizeof(PLASMA_enum),                &norm,  VALUE,
            sizeof(PLASMA_enum),                &uplo,  VALUE,
            sizeof(PLASMA_enum),                &diag,  VALUE,
            sizeof(int),                        &M,     VALUE,
            sizeof(int),                        &N,     VALUE,
            sizeof(double)*szeA,     A,             INPUT,
            sizeof(int),                        &LDA,   VALUE,
            sizeof(double)*szeW,                 NULL,          SCRATCH,
            sizeof(double),                      result,        OUTPUT | GATHERV,
            0);
    } else {
        QUARK_Insert_Task(quark, CORE_dlantr_f1_quark, task_flags,
            sizeof(PLASMA_enum),                &norm,  VALUE,
            sizeof(PLASMA_enum),                &uplo,  VALUE,
            sizeof(PLASMA_enum),                &diag,  VALUE,
            sizeof(int),                        &M,     VALUE,
            sizeof(int),                        &N,     VALUE,
            sizeof(double)*szeA,     A,             INPUT,
            sizeof(int),                        &LDA,   VALUE,
            sizeof(double)*szeW,                 NULL,          SCRATCH,
            sizeof(double),                      result,        OUTPUT,
            sizeof(double)*szeF,                 fake,          OUTPUT | GATHERV,
            0);
    }
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_dlantr_quark = PCORE_dlantr_quark
#define CORE_dlantr_quark PCORE_dlantr_quark
#endif
void CORE_dlantr_quark(Quark *quark)
{
    double *normA;
    PLASMA_enum norm, uplo, diag;
    int M;
    int N;
    double *A;
    int LDA;
    double *work;

    quark_unpack_args_9(quark, norm, uplo, diag, M, N, A, LDA, work, normA);
    *normA = LAPACKE_dlantr_work(
        LAPACK_COL_MAJOR,
        lapack_const(norm),
        lapack_const(uplo),
        lapack_const(diag),
        M, N, A, LDA, work);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_dlantr_f1_quark = PCORE_dlantr_f1_quark
#define CORE_dlantr_f1_quark PCORE_dlantr_f1_quark
#endif
void CORE_dlantr_f1_quark(Quark *quark)
{
    double *normA;
    PLASMA_enum norm, uplo, diag;
    int M;
    int N;
    double *A;
    int LDA;
    double *work;
    double *fake;

    quark_unpack_args_10(quark, norm, uplo, diag, M, N, A, LDA, work, normA, fake);
    *normA = LAPACKE_dlantr_work(
        LAPACK_COL_MAJOR,
        lapack_const(norm),
        lapack_const(uplo),
        lapack_const(diag),
        M, N, A, LDA, work);
}


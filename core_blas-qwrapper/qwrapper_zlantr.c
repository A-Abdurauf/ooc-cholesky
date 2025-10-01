/**
 *
 * @file qwrapper_zlantr.c
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
#ifdef PLASMA_WITH_MKL
#include <mkl_lapacke.h>
#else
#include <lapacke.h>
#endif
#include "common.h"

void
CORE_zlantr_quark(Quark *quark);
void
CORE_zlantr_f1_quark(Quark *quark);

/***************************************************************************//**
 *
 **/
void QUARK_CORE_zlantr(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum norm, PLASMA_enum uplo, PLASMA_enum diag, int M, int N,
                       const PLASMA_Complex64_t *A, int LDA, int szeA,
                       int szeW, double *result)
{
    szeW = max(1, szeW);
    plasma_profile_by_kernel( task_flags, LANGE );

    QUARK_Insert_Task(quark, CORE_zlantr_quark, task_flags,
        sizeof(PLASMA_enum),                &norm,  VALUE,
        sizeof(PLASMA_enum),                &uplo,  VALUE,
        sizeof(PLASMA_enum),                &diag,  VALUE,
        sizeof(int),                        &M,     VALUE,
        sizeof(int),                        &N,     VALUE,
        sizeof(PLASMA_Complex64_t)*szeA,     A,             INPUT,
        sizeof(int),                        &LDA,   VALUE,
        sizeof(double)*szeW,                 NULL,          SCRATCH,
        sizeof(double),                      result,        OUTPUT,
        0);
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_zlantr_f1(Quark *quark, Quark_Task_Flags *task_flags,
                          PLASMA_enum norm, PLASMA_enum uplo, PLASMA_enum diag, int M, int N,
                          const PLASMA_Complex64_t *A, int LDA, int szeA,
                          int szeW, double *result,
                          double *fake, int szeF)
{
    szeW = max(1, szeW);
    plasma_profile_by_kernel( task_flags, LANGE );


    if ( result == fake ) {
        QUARK_Insert_Task(quark, CORE_zlantr_quark, task_flags,
            sizeof(PLASMA_enum),                &norm,  VALUE,
            sizeof(PLASMA_enum),                &uplo,  VALUE,
            sizeof(PLASMA_enum),                &diag,  VALUE,
            sizeof(int),                        &M,     VALUE,
            sizeof(int),                        &N,     VALUE,
            sizeof(PLASMA_Complex64_t)*szeA,     A,             INPUT,
            sizeof(int),                        &LDA,   VALUE,
            sizeof(double)*szeW,                 NULL,          SCRATCH,
            sizeof(double),                      result,        OUTPUT | GATHERV,
            0);
    } else {
        QUARK_Insert_Task(quark, CORE_zlantr_f1_quark, task_flags,
            sizeof(PLASMA_enum),                &norm,  VALUE,
            sizeof(PLASMA_enum),                &uplo,  VALUE,
            sizeof(PLASMA_enum),                &diag,  VALUE,
            sizeof(int),                        &M,     VALUE,
            sizeof(int),                        &N,     VALUE,
            sizeof(PLASMA_Complex64_t)*szeA,     A,             INPUT,
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
#pragma weak CORE_zlantr_quark = PCORE_zlantr_quark
#define CORE_zlantr_quark PCORE_zlantr_quark
#endif
void CORE_zlantr_quark(Quark *quark)
{
    double *normA;
    PLASMA_enum norm, uplo, diag;
    int M;
    int N;
    PLASMA_Complex64_t *A;
    int LDA;
    double *work;

    quark_unpack_args_9(quark, norm, uplo, diag, M, N, A, LDA, work, normA);
    *normA = LAPACKE_zlantr_work(
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
#pragma weak CORE_zlantr_f1_quark = PCORE_zlantr_f1_quark
#define CORE_zlantr_f1_quark PCORE_zlantr_f1_quark
#endif
void CORE_zlantr_f1_quark(Quark *quark)
{
    double *normA;
    PLASMA_enum norm, uplo, diag;
    int M;
    int N;
    PLASMA_Complex64_t *A;
    int LDA;
    double *work;
    double *fake;

    quark_unpack_args_10(quark, norm, uplo, diag, M, N, A, LDA, work, normA, fake);
    *normA = LAPACKE_zlantr_work(
        LAPACK_COL_MAJOR,
        lapack_const(norm),
        lapack_const(uplo),
        lapack_const(diag),
        M, N, A, LDA, work);
}


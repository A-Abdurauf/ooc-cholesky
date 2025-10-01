/**
 *
 * @file qwrapper_slansy.c
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
 * @generated s Fri Apr  1 11:02:41 2016
 *
 **/
#ifdef PLASMA_WITH_MKL
#include <mkl_lapacke.h>
#else
#include <lapacke.h>
#endif
#include "common.h"

/***************************************************************************//**
 *
 **/
void QUARK_CORE_slansy(Quark *quark, Quark_Task_Flags *task_flags,
                       int norm, PLASMA_enum uplo, int N,
                       const float *A, int LDA, int szeA,
                       int szeW, float *result)
{
    szeW = max(1, szeW);
    plasma_profile_by_kernel( task_flags, LANSY );

    QUARK_Insert_Task(quark, CORE_slansy_quark, task_flags,
        sizeof(PLASMA_enum),                &norm,  VALUE,
        sizeof(PLASMA_enum),                &uplo,  VALUE,
        sizeof(int),                        &N,     VALUE,
        sizeof(float)*szeA,     A,             INPUT,
        sizeof(int),                        &LDA,   VALUE,
        sizeof(float)*szeW,                 NULL,          SCRATCH,
        sizeof(float),                      result,        OUTPUT,
        0);
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_slansy_f1(Quark *quark, Quark_Task_Flags *task_flags,
                          PLASMA_enum norm, PLASMA_enum uplo, int N,
                          const float *A, int LDA, int szeA,
                          int szeW, float *result,
                          float *fake, int szeF)
{
    szeW = max(1, szeW);
    plasma_profile_by_kernel( task_flags, LANSY );


    if ( result == fake ) {
        QUARK_Insert_Task(quark, CORE_slansy_quark, task_flags,
            sizeof(PLASMA_enum),                &norm,  VALUE,
            sizeof(PLASMA_enum),                &uplo,  VALUE,
            sizeof(int),                        &N,     VALUE,
            sizeof(float)*szeA,     A,             INPUT,
            sizeof(int),                        &LDA,   VALUE,
            sizeof(float)*szeW,                 NULL,          SCRATCH,
            sizeof(float)*szeF,                 result,        OUTPUT | GATHERV,
            0);
    } else {
        QUARK_Insert_Task(quark, CORE_slansy_f1_quark, task_flags,
            sizeof(PLASMA_enum),                &norm,  VALUE,
            sizeof(PLASMA_enum),                &uplo,  VALUE,
            sizeof(int),                        &N,     VALUE,
            sizeof(float)*szeA,     A,             INPUT,
            sizeof(int),                        &LDA,   VALUE,
            sizeof(float)*szeW,                 NULL,          SCRATCH,
            sizeof(float),                      result,        OUTPUT,
            sizeof(float)*szeF,                 fake,          OUTPUT | GATHERV,
            0);
    }
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_slansy_quark = PCORE_slansy_quark
#define CORE_slansy_quark PCORE_slansy_quark
#endif
void CORE_slansy_quark(Quark *quark)
{
    float *normA;
    int norm;
    PLASMA_enum uplo;
    int N;
    float *A;
    int LDA;
    float *work;

    quark_unpack_args_7(quark, norm, uplo, N, A, LDA, work, normA);
    *normA = LAPACKE_slansy_work(
        LAPACK_COL_MAJOR,
        lapack_const(norm), lapack_const(uplo),
        N, A, LDA, work);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_slansy_f1_quark = PCORE_slansy_f1_quark
#define CORE_slansy_f1_quark PCORE_slansy_f1_quark
#endif
void CORE_slansy_f1_quark(Quark *quark)
{
    float *normA;
    int norm;
    PLASMA_enum uplo;
    int N;
    float *A;
    int LDA;
    float *work;
    float *fake;

    quark_unpack_args_8(quark, norm, uplo, N, A, LDA, work, normA, fake);
    *normA = LAPACKE_slansy_work(
        LAPACK_COL_MAJOR,
        lapack_const(norm),
        lapack_const(uplo),
        N, A, LDA, work);
}

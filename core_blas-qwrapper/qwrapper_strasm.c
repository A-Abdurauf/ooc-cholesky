/**
 *
 * @file qwrapper_strasm.c
 *
 *  PLASMA core_blas quark wrapper
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @generated s Fri Apr  1 11:02:41 2016
 *
 **/
#include "common.h"

void
CORE_strasm_quark(Quark *quark);
void
CORE_strasm_f1_quark(Quark *quark);

/***************************************************************************//**
 *
 **/
void QUARK_CORE_strasm(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum storev, PLASMA_enum uplo, PLASMA_enum diag, int M, int N,
                       const float *A, int lda, int szeA,
                       float *work, int szeW)
{
    QUARK_Insert_Task(
        quark, CORE_strasm_quark, task_flags,
        sizeof(PLASMA_enum),                &storev,    VALUE,
        sizeof(PLASMA_enum),                &uplo,      VALUE,
        sizeof(PLASMA_enum),                &diag,      VALUE,
        sizeof(int),                        &M,         VALUE,
        sizeof(int),                        &N,         VALUE,
        sizeof(float)*szeA,     A,                 INPUT,
        sizeof(int),                        &lda,       VALUE,
        sizeof(float)*szeW,                 work,              INOUT,
        0);
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_strasm_f1(Quark *quark, Quark_Task_Flags *task_flags,
                          PLASMA_enum storev, PLASMA_enum uplo, PLASMA_enum diag, int M, int N,
                          const float *A, int lda, int szeA,
                          float *work, int szeW, float *fake, int szeF)
{
    plasma_profile_by_kernel( task_flags, ASUM );

    if ( work == fake ) {
        QUARK_Insert_Task(
            quark, CORE_strasm_quark, task_flags,
            sizeof(PLASMA_enum),                &storev,    VALUE,
            sizeof(PLASMA_enum),                &uplo,      VALUE,
            sizeof(PLASMA_enum),                &diag,      VALUE,
            sizeof(int),                        &M,         VALUE,
            sizeof(int),                        &N,         VALUE,
            sizeof(float)*szeA,     A,                 INPUT,
            sizeof(int),                        &lda,       VALUE,
            sizeof(float)*szeW,                 work,              INOUT | GATHERV,
            0);
    } else {
        QUARK_Insert_Task(
            quark, CORE_strasm_f1_quark, task_flags,
            sizeof(PLASMA_enum),                &storev,    VALUE,
            sizeof(PLASMA_enum),                &uplo,      VALUE,
            sizeof(PLASMA_enum),                &diag,      VALUE,
            sizeof(int),                        &M,         VALUE,
            sizeof(int),                        &N,         VALUE,
            sizeof(float)*szeA,     A,                 INPUT,
            sizeof(int),                        &lda,       VALUE,
            sizeof(float)*szeW,                 work,              INOUT,
            sizeof(float)*szeF,                 fake,              INOUT | GATHERV,
            0);
    }
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_strasm_quark = PCORE_strasm_quark
#define CORE_strasm_quark PCORE_strasm_quark
#endif
void CORE_strasm_quark(Quark *quark)
{
    PLASMA_enum storev;
    PLASMA_enum uplo;
    PLASMA_enum diag;
    int M;
    int N;
    float *A;
    int lda;
    float *work;

    quark_unpack_args_8(quark, storev, uplo, diag, M, N, A, lda, work);
    CORE_strasm(storev, uplo, diag, M, N, A, lda, work);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_strasm_f1_quark = PCORE_strasm_f1_quark
#define CORE_strasm_f1_quark PCORE_strasm_f1_quark
#endif
void CORE_strasm_f1_quark(Quark *quark)
{
    PLASMA_enum storev;
    PLASMA_enum uplo;
    PLASMA_enum diag;
    int M;
    int N;
    float *A;
    int lda;
    float *work;
    float *fake;

    quark_unpack_args_9(quark, storev, uplo, diag, M, N, A, lda, work, fake);
    CORE_strasm(storev, uplo, diag, M, N, A, lda, work);
}

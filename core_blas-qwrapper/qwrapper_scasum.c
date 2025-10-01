/**
 *
 * @file qwrapper_scasum.c
 *
 *  PLASMA core_blas quark wrapper
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @generated c Fri Apr  1 11:02:41 2016
 *
 **/
#include "common.h"

/***************************************************************************//**
 *
 **/
void QUARK_CORE_scasum(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum storev, PLASMA_enum uplo, int M, int N,
                       const PLASMA_Complex32_t *A, int lda, int szeA,
                       float *work, int szeW)
{
    QUARK_Insert_Task(
        quark, CORE_scasum_quark, task_flags,
        sizeof(PLASMA_enum),                &storev,    VALUE,
        sizeof(PLASMA_enum),                &uplo,      VALUE,
        sizeof(int),                        &M,         VALUE,
        sizeof(int),                        &N,         VALUE,
        sizeof(PLASMA_Complex32_t)*szeA,     A,                 INPUT,
        sizeof(int),                        &lda,       VALUE,
        sizeof(float)*szeW,                 work,              INOUT,
        0);
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_scasum_f1(Quark *quark, Quark_Task_Flags *task_flags,
                          PLASMA_enum storev, PLASMA_enum uplo, int M, int N,
                          const PLASMA_Complex32_t *A, int lda, int szeA,
                          float *work, int szeW, float *fake, int szeF)
{
    plasma_profile_by_kernel( task_flags, ASUM );

    if ( work == fake ) {
        QUARK_Insert_Task(
            quark, CORE_scasum_quark, task_flags,
            sizeof(PLASMA_enum),                &storev,    VALUE,
            sizeof(PLASMA_enum),                &uplo,      VALUE,
            sizeof(int),                        &M,         VALUE,
            sizeof(int),                        &N,         VALUE,
            sizeof(PLASMA_Complex32_t)*szeA,     A,                 INPUT,
            sizeof(int),                        &lda,       VALUE,
            sizeof(float)*szeW,                 work,              INOUT | GATHERV,
            0);
    } else {
        QUARK_Insert_Task(
            quark, CORE_scasum_f1_quark, task_flags,
            sizeof(PLASMA_enum),                &storev,    VALUE,
            sizeof(PLASMA_enum),                &uplo,      VALUE,
            sizeof(int),                        &M,         VALUE,
            sizeof(int),                        &N,         VALUE,
            sizeof(PLASMA_Complex32_t)*szeA,     A,                 INPUT,
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
#pragma weak CORE_scasum_quark = PCORE_scasum_quark
#define CORE_scasum_quark PCORE_scasum_quark
#endif
void CORE_scasum_quark(Quark *quark)
{
    PLASMA_enum storev;
    PLASMA_enum uplo;
    int M;
    int N;
    PLASMA_Complex32_t *A;
    int lda;
    float *work;

    quark_unpack_args_7(quark, storev, uplo, M, N, A, lda, work);
    CORE_scasum(storev, uplo, M, N, A, lda, work);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_scasum_f1_quark = PCORE_scasum_f1_quark
#define CORE_scasum_f1_quark PCORE_scasum_f1_quark
#endif
void CORE_scasum_f1_quark(Quark *quark)
{
    PLASMA_enum storev;
    PLASMA_enum uplo;
    int M;
    int N;
    PLASMA_Complex32_t *A;
    int lda;
    float *work;
    float *fake;

    quark_unpack_args_8(quark, storev, uplo, M, N, A, lda, work, fake);
    CORE_scasum(storev, uplo, M, N, A, lda, work);
}

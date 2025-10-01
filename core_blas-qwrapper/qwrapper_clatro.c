/**
 *
 * @file qwrapper_clatro.c
 *
 *  PLASMA core_blas quark wrapper
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Azzam Haidar
 * @date 2010-11-15
 * @generated c Fri Apr  1 11:02:43 2016
 *
 **/
#include "common.h"

/***************************************************************************//**
 *
 **/
void QUARK_CORE_clatro(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum uplo, PLASMA_enum trans, int m, int n, int nb,
                       const PLASMA_Complex32_t *A, int lda,
                             PLASMA_Complex32_t *B, int ldb)
{
    plasma_profile_by_kernel( task_flags, LATRO );

    QUARK_Insert_Task(quark, CORE_clatro_quark, task_flags,
        sizeof(PLASMA_enum),                &uplo,  VALUE,
        sizeof(PLASMA_enum),                &trans, VALUE,
        sizeof(int),                        &m,     VALUE,
        sizeof(int),                        &n,     VALUE,
        sizeof(PLASMA_Complex32_t)*nb*nb,    A,             INPUT,
        sizeof(int),                        &lda,   VALUE,
        sizeof(PLASMA_Complex32_t)*nb*nb,    B,             OUTPUT,
        sizeof(int),                        &ldb,   VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_clatro_quark = PCORE_clatro_quark
#define CORE_clatro_quark PCORE_clatro_quark
#endif
void CORE_clatro_quark(Quark *quark)
{
    PLASMA_enum uplo;
    PLASMA_enum trans;
    int M;
    int N;
    const PLASMA_Complex32_t *A;
    int LDA;
    PLASMA_Complex32_t *B;
    int LDB;

    quark_unpack_args_8(quark, uplo, trans, M, N, A, LDA, B, LDB);
    CORE_clatro(uplo, trans, M, N, A, LDA, B, LDB);
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_clatro_f1(Quark *quark, Quark_Task_Flags *task_flags,
                          PLASMA_enum uplo, PLASMA_enum trans, int m, int n, int nb,
                          const PLASMA_Complex32_t *A, int lda,
                                PLASMA_Complex32_t *B, int ldb,
                          PLASMA_Complex32_t *fake1, int szefake1, int flag1)
{
    plasma_profile_by_kernel( task_flags, LATRO );

    if ( fake1 == B ) {
        QUARK_Insert_Task(quark, CORE_clatro_quark, task_flags,
            sizeof(PLASMA_enum),                &uplo,  VALUE,
            sizeof(PLASMA_enum),                &trans, VALUE,
            sizeof(int),                        &m,     VALUE,
            sizeof(int),                        &n,     VALUE,
            sizeof(PLASMA_Complex32_t)*nb*nb,    A,             INPUT,
            sizeof(int),                        &lda,   VALUE,
            sizeof(PLASMA_Complex32_t)*nb*nb,    B,             OUTPUT,
            sizeof(int),                        &ldb,   VALUE,
            0);
    } else {
        QUARK_Insert_Task(quark, CORE_clatro_f1_quark, task_flags,
            sizeof(PLASMA_enum),                &uplo,  VALUE,
            sizeof(PLASMA_enum),                &trans, VALUE,
            sizeof(int),                        &m,     VALUE,
            sizeof(int),                        &n,     VALUE,
            sizeof(PLASMA_Complex32_t)*nb*nb,    A,             INPUT,
            sizeof(int),                        &lda,   VALUE,
            sizeof(PLASMA_Complex32_t)*nb*nb,    B,             OUTPUT,
            sizeof(int),                        &ldb,   VALUE,
            sizeof(PLASMA_Complex32_t)*szefake1, fake1,         flag1,
            0);
    }
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_clatro_f1_quark = PCORE_clatro_f1_quark
#define CORE_clatro_f1_quark PCORE_clatro_f1_quark
#endif
void CORE_clatro_f1_quark(Quark *quark)
{
    PLASMA_enum uplo;
    PLASMA_enum trans;
    int M;
    int N;
    const PLASMA_Complex32_t *A;
    int LDA;
    PLASMA_Complex32_t *B;
    int LDB;
    void *fake1;

    quark_unpack_args_9(quark, uplo, trans, M, N, A, LDA, B, LDB, fake1);
    CORE_clatro(uplo, trans, M, N, A, LDA, B, LDB);
}


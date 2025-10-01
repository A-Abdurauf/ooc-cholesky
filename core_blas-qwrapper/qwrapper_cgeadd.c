/**
 *
 * @file qwrapper_cgeadd.c
 *
 *  PLASMA core_blas quark wrapper
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @generated c Fri Apr  1 11:02:40 2016
 *
 **/
#include "common.h"

/***************************************************************************//**
 *
 **/
void QUARK_CORE_cgeadd(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum trans, int m, int n, int nb,
                       PLASMA_Complex32_t alpha,
                       const PLASMA_Complex32_t *A, int lda,
                       PLASMA_Complex32_t beta,
                             PLASMA_Complex32_t *B, int ldb)
{
    plasma_profile_by_kernel( task_flags, GEADD );

    QUARK_Insert_Task(quark, CORE_cgeadd_quark, task_flags,
        sizeof(PLASMA_enum),                &trans, VALUE,
        sizeof(int),                        &m,     VALUE,
        sizeof(int),                        &n,     VALUE,
        sizeof(PLASMA_Complex32_t),         &alpha, VALUE,
        sizeof(PLASMA_Complex32_t)*nb*nb,    A,             INPUT,
        sizeof(int),                        &lda,   VALUE,
        sizeof(PLASMA_Complex32_t),         &beta,  VALUE,
        sizeof(PLASMA_Complex32_t)*nb*nb,    B,             INOUT,
        sizeof(int),                        &ldb,   VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_cgeadd_quark = PCORE_cgeadd_quark
#define CORE_cgeadd_quark PCORE_cgeadd_quark
#endif
void CORE_cgeadd_quark(Quark *quark)
{
    PLASMA_enum trans;
    int M;
    int N;
    PLASMA_Complex32_t alpha;
    PLASMA_Complex32_t *A;
    int LDA;
    PLASMA_Complex32_t beta;
    PLASMA_Complex32_t *B;
    int LDB;

    quark_unpack_args_9(quark, trans, M, N, alpha, A, LDA, beta, B, LDB);
    CORE_cgeadd(trans, M, N, alpha, A, LDA, beta, B, LDB);
    return;
}

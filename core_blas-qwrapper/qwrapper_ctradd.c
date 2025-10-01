/**
 *
 * @file qwrapper_ctradd.c
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

void
CORE_ctradd_quark(Quark *quark);

/***************************************************************************//**
 *
 **/
void QUARK_CORE_ctradd(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum uplo, PLASMA_enum trans, int m, int n, int nb,
                       PLASMA_Complex32_t alpha,
                       const PLASMA_Complex32_t *A, int lda,
                       PLASMA_Complex32_t beta,
                             PLASMA_Complex32_t *B, int ldb)
{
    plasma_profile_by_kernel( task_flags, GEADD );

    QUARK_Insert_Task(quark, CORE_ctradd_quark, task_flags,
        sizeof(PLASMA_enum),                &uplo,  VALUE,
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
#pragma weak CORE_ctradd_quark = PCORE_ctradd_quark
#define CORE_ctradd_quark PCORE_ctradd_quark
#endif
void CORE_ctradd_quark(Quark *quark)
{
    PLASMA_enum uplo, trans;
    int M;
    int N;
    PLASMA_Complex32_t alpha;
    PLASMA_Complex32_t *A;
    int LDA;
    PLASMA_Complex32_t beta;
    PLASMA_Complex32_t *B;
    int LDB;

    quark_unpack_args_10(quark, uplo, trans, M, N, alpha, A, LDA, beta, B, LDB);
    CORE_ctradd(uplo, trans, M, N, alpha, A, LDA, beta, B, LDB);
    return;
}

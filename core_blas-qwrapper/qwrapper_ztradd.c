/**
 *
 * @file qwrapper_ztradd.c
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
#include "common.h"

void
CORE_ztradd_quark(Quark *quark);

/***************************************************************************//**
 *
 **/
void QUARK_CORE_ztradd(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum uplo, PLASMA_enum trans, int m, int n, int nb,
                       PLASMA_Complex64_t alpha,
                       const PLASMA_Complex64_t *A, int lda,
                       PLASMA_Complex64_t beta,
                             PLASMA_Complex64_t *B, int ldb)
{
    plasma_profile_by_kernel( task_flags, GEADD );

    QUARK_Insert_Task(quark, CORE_ztradd_quark, task_flags,
        sizeof(PLASMA_enum),                &uplo,  VALUE,
        sizeof(PLASMA_enum),                &trans, VALUE,
        sizeof(int),                        &m,     VALUE,
        sizeof(int),                        &n,     VALUE,
        sizeof(PLASMA_Complex64_t),         &alpha, VALUE,
        sizeof(PLASMA_Complex64_t)*nb*nb,    A,             INPUT,
        sizeof(int),                        &lda,   VALUE,
        sizeof(PLASMA_Complex64_t),         &beta,  VALUE,
        sizeof(PLASMA_Complex64_t)*nb*nb,    B,             INOUT,
        sizeof(int),                        &ldb,   VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_ztradd_quark = PCORE_ztradd_quark
#define CORE_ztradd_quark PCORE_ztradd_quark
#endif
void CORE_ztradd_quark(Quark *quark)
{
    PLASMA_enum uplo, trans;
    int M;
    int N;
    PLASMA_Complex64_t alpha;
    PLASMA_Complex64_t *A;
    int LDA;
    PLASMA_Complex64_t beta;
    PLASMA_Complex64_t *B;
    int LDB;

    quark_unpack_args_10(quark, uplo, trans, M, N, alpha, A, LDA, beta, B, LDB);
    CORE_ztradd(uplo, trans, M, N, alpha, A, LDA, beta, B, LDB);
    return;
}

/**
 *
 * @file qwrapper_stradd.c
 *
 *  PLASMA core_blas quark wrapper
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @generated s Fri Apr  1 11:02:40 2016
 *
 **/
#include "common.h"

void
CORE_stradd_quark(Quark *quark);

/***************************************************************************//**
 *
 **/
void QUARK_CORE_stradd(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum uplo, PLASMA_enum trans, int m, int n, int nb,
                       float alpha,
                       const float *A, int lda,
                       float beta,
                             float *B, int ldb)
{
    plasma_profile_by_kernel( task_flags, GEADD );

    QUARK_Insert_Task(quark, CORE_stradd_quark, task_flags,
        sizeof(PLASMA_enum),                &uplo,  VALUE,
        sizeof(PLASMA_enum),                &trans, VALUE,
        sizeof(int),                        &m,     VALUE,
        sizeof(int),                        &n,     VALUE,
        sizeof(float),         &alpha, VALUE,
        sizeof(float)*nb*nb,    A,             INPUT,
        sizeof(int),                        &lda,   VALUE,
        sizeof(float),         &beta,  VALUE,
        sizeof(float)*nb*nb,    B,             INOUT,
        sizeof(int),                        &ldb,   VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_stradd_quark = PCORE_stradd_quark
#define CORE_stradd_quark PCORE_stradd_quark
#endif
void CORE_stradd_quark(Quark *quark)
{
    PLASMA_enum uplo, trans;
    int M;
    int N;
    float alpha;
    float *A;
    int LDA;
    float beta;
    float *B;
    int LDB;

    quark_unpack_args_10(quark, uplo, trans, M, N, alpha, A, LDA, beta, B, LDB);
    CORE_stradd(uplo, trans, M, N, alpha, A, LDA, beta, B, LDB);
    return;
}

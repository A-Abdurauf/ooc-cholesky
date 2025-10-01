/**
 *
 * @file qwrapper_dtradd.c
 *
 *  PLASMA core_blas quark wrapper
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @generated d Fri Apr  1 11:02:40 2016
 *
 **/
#include "common.h"

void
CORE_dtradd_quark(Quark *quark);

/***************************************************************************//**
 *
 **/
void QUARK_CORE_dtradd(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum uplo, PLASMA_enum trans, int m, int n, int nb,
                       double alpha,
                       const double *A, int lda,
                       double beta,
                             double *B, int ldb)
{
    plasma_profile_by_kernel( task_flags, GEADD );

    QUARK_Insert_Task(quark, CORE_dtradd_quark, task_flags,
        sizeof(PLASMA_enum),                &uplo,  VALUE,
        sizeof(PLASMA_enum),                &trans, VALUE,
        sizeof(int),                        &m,     VALUE,
        sizeof(int),                        &n,     VALUE,
        sizeof(double),         &alpha, VALUE,
        sizeof(double)*nb*nb,    A,             INPUT,
        sizeof(int),                        &lda,   VALUE,
        sizeof(double),         &beta,  VALUE,
        sizeof(double)*nb*nb,    B,             INOUT,
        sizeof(int),                        &ldb,   VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_dtradd_quark = PCORE_dtradd_quark
#define CORE_dtradd_quark PCORE_dtradd_quark
#endif
void CORE_dtradd_quark(Quark *quark)
{
    PLASMA_enum uplo, trans;
    int M;
    int N;
    double alpha;
    double *A;
    int LDA;
    double beta;
    double *B;
    int LDB;

    quark_unpack_args_10(quark, uplo, trans, M, N, alpha, A, LDA, beta, B, LDB);
    CORE_dtradd(uplo, trans, M, N, alpha, A, LDA, beta, B, LDB);
    return;
}

/**
 *
 * @file qwrapper_csymm.c
 *
 *  PLASMA core_blas quark wrapper
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Hatem Ltaief
 * @author Mathieu Faverge
 * @author Jakub Kurzak
 * @date 2010-11-15
 * @generated c Fri Apr  1 11:02:40 2016
 *
 **/
#include "common.h"

/***************************************************************************//**
 *
 **/
void QUARK_CORE_csymm(Quark *quark, Quark_Task_Flags *task_flags,
                      PLASMA_enum side, PLASMA_enum uplo,
                      int m, int n, int nb,
                      PLASMA_Complex32_t alpha, const PLASMA_Complex32_t *A, int lda,
                      const PLASMA_Complex32_t *B, int ldb,
                      PLASMA_Complex32_t beta, PLASMA_Complex32_t *C, int ldc)
{
    plasma_profile_by_kernel( task_flags, SYMM );

    QUARK_Insert_Task(quark, CORE_csymm_quark, task_flags,
        sizeof(PLASMA_enum),                &side,    VALUE,
        sizeof(PLASMA_enum),                &uplo,    VALUE,
        sizeof(int),                        &m,       VALUE,
        sizeof(int),                        &n,       VALUE,
        sizeof(PLASMA_Complex32_t),         &alpha,   VALUE,
        sizeof(PLASMA_Complex32_t)*nb*nb,    A,               INPUT,
        sizeof(int),                        &lda,     VALUE,
        sizeof(PLASMA_Complex32_t)*nb*nb,    B,               INPUT,
        sizeof(int),                        &ldb,     VALUE,
        sizeof(PLASMA_Complex32_t),         &beta,    VALUE,
        sizeof(PLASMA_Complex32_t)*nb*nb,    C,               INOUT,
        sizeof(int),                        &ldc,     VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_csymm_quark = PCORE_csymm_quark
#define CORE_csymm_quark PCORE_csymm_quark
#endif
void CORE_csymm_quark(Quark *quark)
{
    PLASMA_enum side;
    PLASMA_enum uplo;
    int M;
    int N;
    PLASMA_Complex32_t alpha;
    PLASMA_Complex32_t *A;
    int LDA;
    PLASMA_Complex32_t *B;
    int LDB;
    PLASMA_Complex32_t beta;
    PLASMA_Complex32_t *C;
    int LDC;

    quark_unpack_args_12(quark, side, uplo, M, N, alpha, A, LDA, B, LDB, beta, C, LDC);
    cblas_csymm(
        CblasColMajor,
        (CBLAS_SIDE)side, (CBLAS_UPLO)uplo,
        M, N,
        CBLAS_SADDR(alpha), A, LDA,
        B, LDB,
        CBLAS_SADDR(beta), C, LDC);
}

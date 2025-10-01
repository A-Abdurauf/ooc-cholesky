/**
 *
 * @file qwrapper_ztrasm.c
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
CORE_ztrasm_quark(Quark *quark);
void
CORE_ztrasm_f1_quark(Quark *quark);

/***************************************************************************//**
 *
 **/
void QUARK_CORE_ztrasm(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum storev, PLASMA_enum uplo, PLASMA_enum diag, int M, int N,
                       const PLASMA_Complex64_t *A, int lda, int szeA,
                       double *work, int szeW)
{
    QUARK_Insert_Task(
        quark, CORE_ztrasm_quark, task_flags,
        sizeof(PLASMA_enum),                &storev,    VALUE,
        sizeof(PLASMA_enum),                &uplo,      VALUE,
        sizeof(PLASMA_enum),                &diag,      VALUE,
        sizeof(int),                        &M,         VALUE,
        sizeof(int),                        &N,         VALUE,
        sizeof(PLASMA_Complex64_t)*szeA,     A,                 INPUT,
        sizeof(int),                        &lda,       VALUE,
        sizeof(double)*szeW,                 work,              INOUT,
        0);
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_ztrasm_f1(Quark *quark, Quark_Task_Flags *task_flags,
                          PLASMA_enum storev, PLASMA_enum uplo, PLASMA_enum diag, int M, int N,
                          const PLASMA_Complex64_t *A, int lda, int szeA,
                          double *work, int szeW, double *fake, int szeF)
{
    plasma_profile_by_kernel( task_flags, ASUM );

    if ( work == fake ) {
        QUARK_Insert_Task(
            quark, CORE_ztrasm_quark, task_flags,
            sizeof(PLASMA_enum),                &storev,    VALUE,
            sizeof(PLASMA_enum),                &uplo,      VALUE,
            sizeof(PLASMA_enum),                &diag,      VALUE,
            sizeof(int),                        &M,         VALUE,
            sizeof(int),                        &N,         VALUE,
            sizeof(PLASMA_Complex64_t)*szeA,     A,                 INPUT,
            sizeof(int),                        &lda,       VALUE,
            sizeof(double)*szeW,                 work,              INOUT | GATHERV,
            0);
    } else {
        QUARK_Insert_Task(
            quark, CORE_ztrasm_f1_quark, task_flags,
            sizeof(PLASMA_enum),                &storev,    VALUE,
            sizeof(PLASMA_enum),                &uplo,      VALUE,
            sizeof(PLASMA_enum),                &diag,      VALUE,
            sizeof(int),                        &M,         VALUE,
            sizeof(int),                        &N,         VALUE,
            sizeof(PLASMA_Complex64_t)*szeA,     A,                 INPUT,
            sizeof(int),                        &lda,       VALUE,
            sizeof(double)*szeW,                 work,              INOUT,
            sizeof(double)*szeF,                 fake,              INOUT | GATHERV,
            0);
    }
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_ztrasm_quark = PCORE_ztrasm_quark
#define CORE_ztrasm_quark PCORE_ztrasm_quark
#endif
void CORE_ztrasm_quark(Quark *quark)
{
    PLASMA_enum storev;
    PLASMA_enum uplo;
    PLASMA_enum diag;
    int M;
    int N;
    PLASMA_Complex64_t *A;
    int lda;
    double *work;

    quark_unpack_args_8(quark, storev, uplo, diag, M, N, A, lda, work);
    CORE_ztrasm(storev, uplo, diag, M, N, A, lda, work);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_ztrasm_f1_quark = PCORE_ztrasm_f1_quark
#define CORE_ztrasm_f1_quark PCORE_ztrasm_f1_quark
#endif
void CORE_ztrasm_f1_quark(Quark *quark)
{
    PLASMA_enum storev;
    PLASMA_enum uplo;
    PLASMA_enum diag;
    int M;
    int N;
    PLASMA_Complex64_t *A;
    int lda;
    double *work;
    double *fake;

    quark_unpack_args_9(quark, storev, uplo, diag, M, N, A, lda, work, fake);
    CORE_ztrasm(storev, uplo, diag, M, N, A, lda, work);
}

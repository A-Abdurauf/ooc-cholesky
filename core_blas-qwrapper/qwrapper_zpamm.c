/**
 *
 * @file qwrapper_zpamm.c
 *
 *  PLASMA core_blas quark wrapper
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Dulceneia Becker
 * @date 2011-06-14
 * @precisions normal z -> c d s
 *
 **/
#include "common.h"

/***************************************************************************//**
 *
 **/
void
QUARK_CORE_zpamm(Quark *quark, Quark_Task_Flags *task_flags,
                 int op, PLASMA_enum side, int storev,
                 int m, int n, int k, int l,
                 const PLASMA_Complex64_t *A1, int lda1,
                       PLASMA_Complex64_t *A2, int lda2,
                 const PLASMA_Complex64_t *V, int ldv,
                       PLASMA_Complex64_t *W, int ldw)
{
    QUARK_Insert_Task(quark, CORE_zpamm_quark, task_flags,
        sizeof(int),                        &op,      VALUE,
        sizeof(PLASMA_enum),                &side,    VALUE,
        sizeof(PLASMA_enum),                &storev,  VALUE,
        sizeof(int),                        &m,       VALUE,
        sizeof(int),                        &n,       VALUE,
        sizeof(int),                        &k,       VALUE,
        sizeof(int),                        &l,       VALUE,
        sizeof(PLASMA_Complex64_t)*m*k,     A1,           INPUT,
        sizeof(int),                        &lda1,    VALUE,
        sizeof(PLASMA_Complex64_t)*k*n,     A2,           INOUT,
        sizeof(int),                        &lda2,    VALUE,
        sizeof(PLASMA_Complex64_t)*m*n,     V,            INPUT,
        sizeof(int),                        &ldv,     VALUE,
        sizeof(PLASMA_Complex64_t)*m*n,     W,            INOUT,
        sizeof(int),                        &ldw,     VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
void
CORE_zpamm_quark(Quark *quark)
{
    int op;
    PLASMA_enum side;
    int storev;
    int M;
    int N;
    int K;
    int L;
    PLASMA_Complex64_t *A1;
    int LDA1;
    PLASMA_Complex64_t *A2;
    int LDA2;
    PLASMA_Complex64_t *V;
    int LDV;
    PLASMA_Complex64_t *W;
    int LDW;

    quark_unpack_args_15(quark, op, side, storev, M, N, K, L,
            A1, LDA1, A2, LDA2, V, LDV, W, LDW);

    CORE_zpamm( op, side, storev, M, N, K, L, A1, LDA1, A2, LDA2, V, LDV, W, LDW);
}

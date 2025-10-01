/**
 *
 * @file qwrapper_spamm.c
 *
 *  PLASMA core_blas quark wrapper
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Dulceneia Becker
 * @date 2011-06-14
 * @generated s Fri Apr  1 11:02:43 2016
 *
 **/
#include "common.h"

/***************************************************************************//**
 *
 **/
void
QUARK_CORE_spamm(Quark *quark, Quark_Task_Flags *task_flags,
                 int op, PLASMA_enum side, int storev,
                 int m, int n, int k, int l,
                 const float *A1, int lda1,
                       float *A2, int lda2,
                 const float *V, int ldv,
                       float *W, int ldw)
{
    QUARK_Insert_Task(quark, CORE_spamm_quark, task_flags,
        sizeof(int),                        &op,      VALUE,
        sizeof(PLASMA_enum),                &side,    VALUE,
        sizeof(PLASMA_enum),                &storev,  VALUE,
        sizeof(int),                        &m,       VALUE,
        sizeof(int),                        &n,       VALUE,
        sizeof(int),                        &k,       VALUE,
        sizeof(int),                        &l,       VALUE,
        sizeof(float)*m*k,     A1,           INPUT,
        sizeof(int),                        &lda1,    VALUE,
        sizeof(float)*k*n,     A2,           INOUT,
        sizeof(int),                        &lda2,    VALUE,
        sizeof(float)*m*n,     V,            INPUT,
        sizeof(int),                        &ldv,     VALUE,
        sizeof(float)*m*n,     W,            INOUT,
        sizeof(int),                        &ldw,     VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
void
CORE_spamm_quark(Quark *quark)
{
    int op;
    PLASMA_enum side;
    int storev;
    int M;
    int N;
    int K;
    int L;
    float *A1;
    int LDA1;
    float *A2;
    int LDA2;
    float *V;
    int LDV;
    float *W;
    int LDW;

    quark_unpack_args_15(quark, op, side, storev, M, N, K, L,
            A1, LDA1, A2, LDA2, V, LDV, W, LDW);

    CORE_spamm( op, side, storev, M, N, K, L, A1, LDA1, A2, LDA2, V, LDV, W, LDW);
}

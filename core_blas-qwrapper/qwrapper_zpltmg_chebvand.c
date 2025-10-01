/**
 *
 * @file qwrapper_zpltmg_chebvand.c
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
CORE_zpltmg_chebvand_quark(Quark *quark);

/***************************************************************************//**
 *
 **/
void QUARK_CORE_zpltmg_chebvand( Quark *quark, Quark_Task_Flags *task_flags,
                                 int M, int N, PLASMA_Complex64_t *A, int LDA,
                                 int gN, int m0, int n0,
                                 PLASMA_Complex64_t *W )
{
    plasma_profile_by_kernel( task_flags, PLRNT );

    QUARK_Insert_Task(quark, CORE_zpltmg_chebvand_quark, task_flags,
        sizeof(int),                      &M,       VALUE,
        sizeof(int),                      &N,       VALUE,
        sizeof(PLASMA_Complex64_t)*LDA*N, A,            OUTPUT,
        sizeof(int),                      &LDA,     VALUE,
        sizeof(int),                      &gN,      VALUE,
        sizeof(int),                      &m0,      VALUE,
        sizeof(int),                      &n0,      VALUE,
        sizeof(PLASMA_Complex64_t)*N*2,   W,            INOUT,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_zpltmg_chebvand_quark = PCORE_zpltmg_chebvand_quark
#define CORE_zpltmg_chebvand_quark PCORE_zpltmg_chebvand_quark
#endif
void CORE_zpltmg_chebvand_quark(Quark *quark)
{
    int m;
    int n;
    PLASMA_Complex64_t *A;
    int lda;
    int gN;
    int m0;
    int n0;
    PLASMA_Complex64_t *W;

    quark_unpack_args_8( quark, m, n, A, lda, gN, m0, n0, W );
    CORE_zpltmg_chebvand( m, n, A, lda, gN, m0, n0, W );
}

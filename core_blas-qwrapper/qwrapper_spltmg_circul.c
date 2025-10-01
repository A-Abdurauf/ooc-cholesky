/**
 *
 * @file qwrapper_spltmg_circul.c
 *
 *  PLASMA core_blas quark wrapper
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @generated s Fri Apr  1 11:02:42 2016
 *
 **/
#include "common.h"

void
CORE_spltmg_circul_quark(Quark *quark);

/***************************************************************************//**
 *
 **/
void QUARK_CORE_spltmg_circul( Quark *quark, Quark_Task_Flags *task_flags,
                               int M, int N, float *A, int LDA,
                               int gM, int m0, int n0,
                               const float *V )
{
    plasma_profile_by_kernel( task_flags, PLRNT );

    QUARK_Insert_Task(quark, CORE_spltmg_circul_quark, task_flags,
        sizeof(int),                      &M,       VALUE,
        sizeof(int),                      &N,       VALUE,
        sizeof(float)*LDA*N, A,            OUTPUT,
        sizeof(int),                      &LDA,     VALUE,
        sizeof(int),                      &gM,      VALUE,
        sizeof(int),                      &m0,      VALUE,
        sizeof(int),                      &n0,      VALUE,
        sizeof(float)*gM,    V,            INPUT,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_spltmg_circul_quark = PCORE_spltmg_circul_quark
#define CORE_spltmg_circul_quark PCORE_spltmg_circul_quark
#endif
void CORE_spltmg_circul_quark(Quark *quark)
{
    int m;
    int n;
    float *A;
    int lda;
    int gM;
    int m0;
    int n0;
    const float *V;

    quark_unpack_args_8( quark, m, n, A, lda, gM, m0, n0, V );
    CORE_spltmg_circul( m, n, A, lda, gM, m0, n0, V );
}

/**
 *
 * @file qwrapper_dpltmg_circul.c
 *
 *  PLASMA core_blas quark wrapper
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @generated d Fri Apr  1 11:02:42 2016
 *
 **/
#include "common.h"

void
CORE_dpltmg_circul_quark(Quark *quark);

/***************************************************************************//**
 *
 **/
void QUARK_CORE_dpltmg_circul( Quark *quark, Quark_Task_Flags *task_flags,
                               int M, int N, double *A, int LDA,
                               int gM, int m0, int n0,
                               const double *V )
{
    plasma_profile_by_kernel( task_flags, PLRNT );

    QUARK_Insert_Task(quark, CORE_dpltmg_circul_quark, task_flags,
        sizeof(int),                      &M,       VALUE,
        sizeof(int),                      &N,       VALUE,
        sizeof(double)*LDA*N, A,            OUTPUT,
        sizeof(int),                      &LDA,     VALUE,
        sizeof(int),                      &gM,      VALUE,
        sizeof(int),                      &m0,      VALUE,
        sizeof(int),                      &n0,      VALUE,
        sizeof(double)*gM,    V,            INPUT,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_dpltmg_circul_quark = PCORE_dpltmg_circul_quark
#define CORE_dpltmg_circul_quark PCORE_dpltmg_circul_quark
#endif
void CORE_dpltmg_circul_quark(Quark *quark)
{
    int m;
    int n;
    double *A;
    int lda;
    int gM;
    int m0;
    int n0;
    const double *V;

    quark_unpack_args_8( quark, m, n, A, lda, gM, m0, n0, V );
    CORE_dpltmg_circul( m, n, A, lda, gM, m0, n0, V );
}

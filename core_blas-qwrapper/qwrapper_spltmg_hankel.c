/**
 *
 * @file qwrapper_spltmg_hankel.c
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
CORE_spltmg_hankel_quark(Quark *quark);

/***************************************************************************//**
 *
 **/
void QUARK_CORE_spltmg_hankel( Quark *quark, Quark_Task_Flags *task_flags,
                               PLASMA_enum uplo, int M, int N, float *A, int LDA,
                               int m0, int n0, int nb,
                               const float *V1,
                               const float *V2)
{
    plasma_profile_by_kernel( task_flags, PLRNT );

    QUARK_Insert_Task(quark, CORE_spltmg_hankel_quark, task_flags,
        sizeof(PLASMA_enum),              &uplo,    VALUE,
        sizeof(int),                      &M,       VALUE,
        sizeof(int),                      &N,       VALUE,
        sizeof(float)*LDA*N, A,            OUTPUT,
        sizeof(int),                      &LDA,     VALUE,
        sizeof(int),                      &m0,      VALUE,
        sizeof(int),                      &n0,      VALUE,
        sizeof(int),                      &nb,      VALUE,
        sizeof(float)*nb,    V1,           INPUT,
        sizeof(float)*nb,    V2,           INPUT,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_spltmg_hankel_quark = PCORE_spltmg_hankel_quark
#define CORE_spltmg_hankel_quark PCORE_spltmg_hankel_quark
#endif
void CORE_spltmg_hankel_quark(Quark *quark)
{
    PLASMA_enum uplo;
    int m;
    int n;
    float *A;
    int lda;
    int m0;
    int n0;
    int nb;
    const float *V1;
    const float *V2;

    quark_unpack_args_10( quark, uplo, m, n, A, lda, m0, n0, nb, V1, V2 );
    CORE_spltmg_hankel( uplo, m, n, A, lda, m0, n0, nb, V1, V2 );
}

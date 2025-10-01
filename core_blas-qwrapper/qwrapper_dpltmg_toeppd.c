/**
 *
 * @file qwrapper_dpltmg_toeppd.c
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
CORE_dpltmg_toeppd1_quark(Quark *quark);
void
CORE_dpltmg_toeppd2_quark(Quark *quark);

/***************************************************************************//**
 *
 **/
void QUARK_CORE_dpltmg_toeppd1(Quark *quark, Quark_Task_Flags *task_flags,
                               int gM, int m0, int M,
                               double *W,
                               unsigned long long int seed)
{
    plasma_profile_by_kernel( task_flags, PLRNT );

    QUARK_Insert_Task(quark, CORE_dpltmg_toeppd1_quark, task_flags,
        sizeof(int),                     &gM,    VALUE,
        sizeof(int),                     &m0,    VALUE,
        sizeof(int),                     &M,     VALUE,
        sizeof(double)*M*2,   W,        OUTPUT,
        sizeof(unsigned long long int),   &seed, VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_dpltmg_toeppd1_quark = PCORE_dpltmg_toeppd1_quark
#define CORE_dpltmg_toeppd1_quark PCORE_dpltmg_toeppd1_quark
#endif
void CORE_dpltmg_toeppd1_quark(Quark *quark)
{
    int gM;
    int m0;
    int M;
    double *W;
    unsigned long long int seed;

    quark_unpack_args_5( quark, gM, m0, M, W, seed );
    CORE_dpltmg_toeppd1( gM, m0, M, W, seed );
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_dpltmg_toeppd2(Quark *quark, Quark_Task_Flags *task_flags,
                               int M, int N, int K, int m0, int n0,
                               const double *W,
                               double *A, int LDA )
{
    plasma_profile_by_kernel( task_flags, PLRNT );

    QUARK_Insert_Task(quark, CORE_dpltmg_toeppd2_quark, task_flags,
        sizeof(int),                     &M,   VALUE,
        sizeof(int),                     &N,   VALUE,
        sizeof(int),                     &K,   VALUE,
        sizeof(int),                     &m0,  VALUE,
        sizeof(int),                     &n0,  VALUE,
        sizeof(double)*K*2,   W,      INPUT,
        sizeof(double)*LDA*N, A,      INOUT,
        sizeof(int),                     &LDA, VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_dpltmg_toeppd2_quark = PCORE_dpltmg_toeppd2_quark
#define CORE_dpltmg_toeppd2_quark PCORE_dpltmg_toeppd2_quark
#endif
void CORE_dpltmg_toeppd2_quark(Quark *quark)
{
    int M;
    int N;
    int K;
    int m0;
    int n0;
    const double *W;
    double *A;
    int LDA;

    quark_unpack_args_8( quark, M, N, K, m0, n0, W, A, LDA );
    CORE_dpltmg_toeppd2( M, N, K, m0, n0, W, A, LDA );
}

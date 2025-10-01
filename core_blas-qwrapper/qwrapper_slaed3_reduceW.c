/**
 *
 * @file qwrapper_slaed3_reduceW.c
 *
 *  PLASMA core_blas quark wrapper
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Gregoire Pichon
 * @author Azzam Haidar
 * @date 2014-07
 * @generated s Fri Apr  1 11:02:43 2016
 *
 **/
#include "common.h"

void
CORE_slaed3_reduceW_quark(Quark *quark);

void
CORE_slaed3_reduceW_p2_quark(Quark *quark);

/***************************************************************************//**
 *
 **/
void QUARK_CORE_slaed3_reduceW(Quark *quark, Quark_Task_Flags *task_flags,
                               int n, int n1, const int *K, int l,
                               const float *Q, int LDQ,
                               const float *Wred, float *W)
{
    plasma_profile_by_kernel( task_flags, LAED3_REDUCEW );

    QUARK_Insert_Task(quark, CORE_slaed3_reduceW_quark, task_flags,
        sizeof(int),     &n,    VALUE,
        sizeof(int),     &n1,   VALUE,
        sizeof(int),      K,        INOUT,
        sizeof(int),     &l,    VALUE,
        sizeof(float),   Q,        NODEP,
        sizeof(int),     &LDQ,  VALUE,
        sizeof(float),   Wred,     INPUT,
        sizeof(float),   W,        OUTPUT,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_slaed3_reduceW_quark = PCORE_slaed3_reduceW_quark
#define CORE_slaed3_reduceW_quark PCORE_slaed3_reduceW_quark
#endif
void CORE_slaed3_reduceW_quark(Quark *quark)
{
    int n, n1, l;
    const int *K;
    const float *Q;
    int LDQ;
    const float *Wred;
    float *W;

    quark_unpack_args_8(quark, n, n1, K, l,
                        Q, LDQ, Wred, W);

    CORE_slaed3_reduceW(n, n1, *K, l,
                        Q, LDQ, Wred, W);
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_slaed3_reduceW_p2(Quark *quark, Quark_Task_Flags *task_flags,
                                  int n, int n1, const int *K, int l,
                                  float **Q, const int *LDQ,
                                  const float *Wred, float *W)
{
    plasma_profile_by_kernel( task_flags, LAED3_REDUCEW );

    QUARK_Insert_Task(quark, CORE_slaed3_reduceW_p2_quark, task_flags,
        sizeof(int),     &n,    VALUE,
        sizeof(int),     &n1,   VALUE,
        sizeof(int),      K,        INOUT,
        sizeof(int),     &l,    VALUE,
        sizeof(float*),  Q,        INOUT,
        sizeof(int),      LDQ,      NODEP,
        sizeof(float),   Wred,     INPUT,
        sizeof(float),   W,        OUTPUT,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_slaed3_reduceW_p2_quark = PCORE_slaed3_reduceW_p2_quark
#define CORE_slaed3_reduceW_p2_quark PCORE_slaed3_reduceW_p2_quark
#endif
void CORE_slaed3_reduceW_p2_quark(Quark *quark)
{
    int n, n1, l;
    const int *K;
    float **Q;
    int *LDQ;
    const float *Wred;
    float *W;

    quark_unpack_args_8(quark, n, n1, K, l,
                        Q, LDQ, Wred, W);

    CORE_slaed3_reduceW(n, n1, *K, l,
                        *Q, *LDQ, Wred, W);
}

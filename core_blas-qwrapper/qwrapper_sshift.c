/**
 *
 * @file qwrapper_sshift.c
 *
 *  PLASMA InPlaceTransformation module
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 *  This work is the implementation of an inplace transformation
 *  based on the GKK algorithm by Gustavson, Karlsson, Kagstrom
 *  and its fortran implementation.
 *
 * @version 2.8.0
 * @author Mathieu Faverge
 * @date 2010-11-15
 *
 * @generated s Fri Apr  1 11:02:41 2016
 *
 **/
#include "common.h"

/***************************************************************************//**
 *
 **/
void QUARK_CORE_sshiftw(Quark *quark, Quark_Task_Flags *task_flags,
                        int s, int cl, int m, int n, int L, float *A, float *W)
{
    plasma_profile_by_kernel( task_flags, SHIFTW );

    QUARK_Insert_Task(quark, CORE_sshiftw_quark, task_flags,
        sizeof(int),                      &s,   VALUE,
        sizeof(int),                      &cl,  VALUE,
        sizeof(int),                      &m,   VALUE,
        sizeof(int),                      &n,   VALUE,
        sizeof(int),                      &L,   VALUE,
        sizeof(float)*m*n*L, A,        INOUT,
        sizeof(float)*L,     W,        INPUT,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_sshiftw_quark = PCORE_sshiftw_quark
#define CORE_sshiftw_quark PCORE_sshiftw_quark
#endif
void CORE_sshiftw_quark(Quark *quark)
{
    int s;
    int cl;
    int m;
    int n;
    int L;
    float *A;
    float *W;

    quark_unpack_args_7(quark, s, cl, m, n, L, A, W);
    CORE_sshiftw(s, cl, m, n, L, A, W);
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_sshift(Quark *quark, Quark_Task_Flags *task_flags,
                       int s, int m, int n, int L, float *A)
{
    plasma_profile_by_kernel( task_flags, SHIFT );

    QUARK_Insert_Task(quark, CORE_sshift_quark, task_flags,
        sizeof(int),                      &s,    VALUE,
        sizeof(int),                      &m,    VALUE,
        sizeof(int),                      &n,    VALUE,
        sizeof(int),                      &L,    VALUE,
        sizeof(float)*m*n*L, A,        INOUT | GATHERV,
        sizeof(float)*L,     NULL,     SCRATCH,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_sshift_quark = PCORE_sshift_quark
#define CORE_sshift_quark PCORE_sshift_quark
#endif
void CORE_sshift_quark(Quark *quark)
{
    int s;
    int m;
    int n;
    int L;
    float *A;
    float *W;

    quark_unpack_args_6(quark, s, m, n, L, A, W);
    memcpy(W, &(A[s*L]), L*sizeof(float));
    CORE_sshiftw(s, 0, m, n, L, A, W);
}


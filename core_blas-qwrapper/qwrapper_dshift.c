/**
 *
 * @file qwrapper_dshift.c
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
 * @generated d Fri Apr  1 11:02:41 2016
 *
 **/
#include "common.h"

/***************************************************************************//**
 *
 **/
void QUARK_CORE_dshiftw(Quark *quark, Quark_Task_Flags *task_flags,
                        int s, int cl, int m, int n, int L, double *A, double *W)
{
    plasma_profile_by_kernel( task_flags, SHIFTW );

    QUARK_Insert_Task(quark, CORE_dshiftw_quark, task_flags,
        sizeof(int),                      &s,   VALUE,
        sizeof(int),                      &cl,  VALUE,
        sizeof(int),                      &m,   VALUE,
        sizeof(int),                      &n,   VALUE,
        sizeof(int),                      &L,   VALUE,
        sizeof(double)*m*n*L, A,        INOUT,
        sizeof(double)*L,     W,        INPUT,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_dshiftw_quark = PCORE_dshiftw_quark
#define CORE_dshiftw_quark PCORE_dshiftw_quark
#endif
void CORE_dshiftw_quark(Quark *quark)
{
    int s;
    int cl;
    int m;
    int n;
    int L;
    double *A;
    double *W;

    quark_unpack_args_7(quark, s, cl, m, n, L, A, W);
    CORE_dshiftw(s, cl, m, n, L, A, W);
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_dshift(Quark *quark, Quark_Task_Flags *task_flags,
                       int s, int m, int n, int L, double *A)
{
    plasma_profile_by_kernel( task_flags, SHIFT );

    QUARK_Insert_Task(quark, CORE_dshift_quark, task_flags,
        sizeof(int),                      &s,    VALUE,
        sizeof(int),                      &m,    VALUE,
        sizeof(int),                      &n,    VALUE,
        sizeof(int),                      &L,    VALUE,
        sizeof(double)*m*n*L, A,        INOUT | GATHERV,
        sizeof(double)*L,     NULL,     SCRATCH,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_dshift_quark = PCORE_dshift_quark
#define CORE_dshift_quark PCORE_dshift_quark
#endif
void CORE_dshift_quark(Quark *quark)
{
    int s;
    int m;
    int n;
    int L;
    double *A;
    double *W;

    quark_unpack_args_6(quark, s, m, n, L, A, W);
    memcpy(W, &(A[s*L]), L*sizeof(double));
    CORE_dshiftw(s, 0, m, n, L, A, W);
}


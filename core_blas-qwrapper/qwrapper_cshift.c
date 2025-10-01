/**
 *
 * @file qwrapper_cshift.c
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
 * @generated c Fri Apr  1 11:02:41 2016
 *
 **/
#include "common.h"

/***************************************************************************//**
 *
 **/
void QUARK_CORE_cshiftw(Quark *quark, Quark_Task_Flags *task_flags,
                        int s, int cl, int m, int n, int L, PLASMA_Complex32_t *A, PLASMA_Complex32_t *W)
{
    plasma_profile_by_kernel( task_flags, SHIFTW );

    QUARK_Insert_Task(quark, CORE_cshiftw_quark, task_flags,
        sizeof(int),                      &s,   VALUE,
        sizeof(int),                      &cl,  VALUE,
        sizeof(int),                      &m,   VALUE,
        sizeof(int),                      &n,   VALUE,
        sizeof(int),                      &L,   VALUE,
        sizeof(PLASMA_Complex32_t)*m*n*L, A,        INOUT,
        sizeof(PLASMA_Complex32_t)*L,     W,        INPUT,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_cshiftw_quark = PCORE_cshiftw_quark
#define CORE_cshiftw_quark PCORE_cshiftw_quark
#endif
void CORE_cshiftw_quark(Quark *quark)
{
    int s;
    int cl;
    int m;
    int n;
    int L;
    PLASMA_Complex32_t *A;
    PLASMA_Complex32_t *W;

    quark_unpack_args_7(quark, s, cl, m, n, L, A, W);
    CORE_cshiftw(s, cl, m, n, L, A, W);
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_cshift(Quark *quark, Quark_Task_Flags *task_flags,
                       int s, int m, int n, int L, PLASMA_Complex32_t *A)
{
    plasma_profile_by_kernel( task_flags, SHIFT );

    QUARK_Insert_Task(quark, CORE_cshift_quark, task_flags,
        sizeof(int),                      &s,    VALUE,
        sizeof(int),                      &m,    VALUE,
        sizeof(int),                      &n,    VALUE,
        sizeof(int),                      &L,    VALUE,
        sizeof(PLASMA_Complex32_t)*m*n*L, A,        INOUT | GATHERV,
        sizeof(PLASMA_Complex32_t)*L,     NULL,     SCRATCH,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_cshift_quark = PCORE_cshift_quark
#define CORE_cshift_quark PCORE_cshift_quark
#endif
void CORE_cshift_quark(Quark *quark)
{
    int s;
    int m;
    int n;
    int L;
    PLASMA_Complex32_t *A;
    PLASMA_Complex32_t *W;

    quark_unpack_args_6(quark, s, m, n, L, A, W);
    memcpy(W, &(A[s*L]), L*sizeof(PLASMA_Complex32_t));
    CORE_cshiftw(s, 0, m, n, L, A, W);
}


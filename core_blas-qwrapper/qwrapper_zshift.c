/**
 *
 * @file qwrapper_zshift.c
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
 * @precisions normal z -> c d s
 *
 **/
#include "common.h"

/***************************************************************************//**
 *
 **/
void QUARK_CORE_zshiftw(Quark *quark, Quark_Task_Flags *task_flags,
                        int s, int cl, int m, int n, int L, PLASMA_Complex64_t *A, PLASMA_Complex64_t *W)
{
    plasma_profile_by_kernel( task_flags, SHIFTW );

    QUARK_Insert_Task(quark, CORE_zshiftw_quark, task_flags,
        sizeof(int),                      &s,   VALUE,
        sizeof(int),                      &cl,  VALUE,
        sizeof(int),                      &m,   VALUE,
        sizeof(int),                      &n,   VALUE,
        sizeof(int),                      &L,   VALUE,
        sizeof(PLASMA_Complex64_t)*m*n*L, A,        INOUT,
        sizeof(PLASMA_Complex64_t)*L,     W,        INPUT,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_zshiftw_quark = PCORE_zshiftw_quark
#define CORE_zshiftw_quark PCORE_zshiftw_quark
#endif
void CORE_zshiftw_quark(Quark *quark)
{
    int s;
    int cl;
    int m;
    int n;
    int L;
    PLASMA_Complex64_t *A;
    PLASMA_Complex64_t *W;

    quark_unpack_args_7(quark, s, cl, m, n, L, A, W);
    CORE_zshiftw(s, cl, m, n, L, A, W);
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_zshift(Quark *quark, Quark_Task_Flags *task_flags,
                       int s, int m, int n, int L, PLASMA_Complex64_t *A)
{
    plasma_profile_by_kernel( task_flags, SHIFT );

    QUARK_Insert_Task(quark, CORE_zshift_quark, task_flags,
        sizeof(int),                      &s,    VALUE,
        sizeof(int),                      &m,    VALUE,
        sizeof(int),                      &n,    VALUE,
        sizeof(int),                      &L,    VALUE,
        sizeof(PLASMA_Complex64_t)*m*n*L, A,        INOUT | GATHERV,
        sizeof(PLASMA_Complex64_t)*L,     NULL,     SCRATCH,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_zshift_quark = PCORE_zshift_quark
#define CORE_zshift_quark PCORE_zshift_quark
#endif
void CORE_zshift_quark(Quark *quark)
{
    int s;
    int m;
    int n;
    int L;
    PLASMA_Complex64_t *A;
    PLASMA_Complex64_t *W;

    quark_unpack_args_6(quark, s, m, n, L, A, W);
    memcpy(W, &(A[s*L]), L*sizeof(PLASMA_Complex64_t));
    CORE_zshiftw(s, 0, m, n, L, A, W);
}


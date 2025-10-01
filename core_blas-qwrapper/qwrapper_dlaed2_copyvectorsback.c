/**
 *
 * @file qwrapper_dlaed2_copyvectorsback.c
 *
 *  PLASMA core_blas quark wrapper
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Gregoire Pichon
 * @author Azzam Haidar
 * @date 2014-07
 * @precisions normal d -> s
 *
 **/
#include "common.h"

void
CORE_dlaed2_copydef_quark(Quark *quark);

/***************************************************************************//**
 *
 **/
void QUARK_CORE_dlaed2_copydef(Quark *quark, Quark_Task_Flags *task_flags,
                               int n, int n1, const int *K, const int *ctot,
                               double *Q, int LDQ, const double *Q2,
                               /*int large_work, double *W_dep,*/
                               int start, int end)
{
    plasma_profile_by_kernel( task_flags, LAED2_COPYDEF );

    QUARK_Insert_Task(quark, CORE_dlaed2_copydef_quark, task_flags,
        sizeof(int),          &n,     VALUE,
        sizeof(int),          &n1,    VALUE,
        sizeof(int),           K,        INPUT,
        sizeof(int)*4,         ctot,     NODEP,
        sizeof(double)*LDQ*n,  Q,        NODEP,
        sizeof(int),          &LDQ,   VALUE,
        sizeof(double),        Q2,       NODEP,
        sizeof(int),          &start, VALUE,
        sizeof(int),          &end,   VALUE,
        /*
         * Fake dependency to guaranty correct execution of kernels running
         * on a same subset of Q. Note that the kernel does not necessarly
         * work on this subset
         */
        sizeof(double),        Q+start*LDQ,  INOUT,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_dlaed2_copydef_quark = PCORE_dlaed2_copydef_quark
#define CORE_dlaed2_copydef_quark PCORE_dlaed2_copydef_quark
#endif
void CORE_dlaed2_copydef_quark(Quark *quark)
{
    int n, n1;
    const int *K;
    int start, end;
    const int *ctot;
    double *Q;
    int LDQ;
    const double *Q2;
    double *fake;

    quark_unpack_args_10(quark, n, n1, K, ctot,
                         Q, LDQ, Q2, start, end, fake );

    CORE_dlaed2_copydef(n, n1, *K, ctot,
                        Q, LDQ, Q2, start, end);
}

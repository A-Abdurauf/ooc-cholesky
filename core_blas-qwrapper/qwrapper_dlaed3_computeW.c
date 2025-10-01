/**
 *
 * @file qwrapper_dlaed3_computeW.c
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
CORE_dlaed3_compW_p2f1_quark(Quark *quark);

/***************************************************************************//**
 *
 **/
void QUARK_CORE_dlaed3_compW_p2f1(Quark *quark, Quark_Task_Flags *task_flags,
                                  int n, const int *K,
                                  double **Q, const int *LDQ,
                                  const double *DLAMBDA, double *W,
                                  const int *INDX,
                                  int start, int end,
                                  void *fakeQ, int flagfQ,
                                  void *fakeW, int flagfW )
{
    int flagQ = NODEP;
    int flagW = OUTPUT;

    if (fakeQ == Q) {
        flagQ = INOUT;
        fakeQ = NULL;
        flagfQ = NODEP;
    }

    if (fakeW == W) {
        flagW |= flagfW;
        fakeW  = NULL;
        flagfW = NODEP;
    }

    plasma_profile_by_kernel( task_flags, LAED3_COMPW );

    QUARK_Insert_Task(quark, CORE_dlaed3_compW_p2f1_quark, task_flags,
        sizeof(int),      &n,         VALUE,
        sizeof(int*),      K,             INPUT,
        /* Q is NODEP to avoid sequentialization of subproblems */
        sizeof(double**),  Q,             flagQ,
        sizeof(int*),      LDQ,           NODEP,
        sizeof(double*),   DLAMBDA,       NODEP,
        sizeof(double*),   W,             flagW,
        sizeof(int*),      INDX,          NODEP,
        sizeof(int),      &start,     VALUE,
        sizeof(int),      &end,       VALUE,
        /*
         * Fake dependency to guaranty correct execution of kernels running
         * on a same subset of Q Note that the kernel does not necessarly
         * work on this subset
         */
        1, fakeQ, flagfQ,
        1, fakeW, flagfW,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_dlaed3_compW_p2f1_quark = PCORE_dlaed3_compW_p2f1_quark
#define CORE_dlaed3_compW_p2f1_quark PCORE_dlaed3_compW_p2f1_quark
#endif
void CORE_dlaed3_compW_p2f1_quark(Quark *quark)
{
    int n;
    const int *K;
    const double **Q;
    const int *LDQ;
    const double *DLAMBDA;
    double *W;
    const int *INDX;
    int start;
    int end;
    const double *fake1;
    void *fake2;

    quark_unpack_args_11(quark, n, K, Q, LDQ, DLAMBDA, W, INDX, start, end, fake1, fake2);
    CORE_dlaed3_computeW(n, *K, *Q, *LDQ, DLAMBDA, W, INDX, start, end);
}

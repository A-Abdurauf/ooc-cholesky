/**
 *
 * @file qwrapper_dlaed1_pipelined.c
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
CORE_dlaed1_pipelined_quark(Quark *quark);

/***************************************************************************//**
 *
 **/
void QUARK_CORE_dlaed1_pipelined(Quark *quark, Quark_Task_Flags *task_flags,
                                 int n, int n1, const int *K,
                                 const int *INDX, const int *ctot,
                                 double *D, const double *beta,
                                 double *Q, int LDQ, double *Q2,
                                 const double *DLAMBDA, const double *W, double *Wred,
                                 int start, int end)
{
    /* Rk: This kernel is not traced, because each subkernel is already traced */
    plasma_gendag_by_kernel( task_flags, LAED1_PIPELINED );

    QUARK_Insert_Task(quark, CORE_dlaed1_pipelined_quark, task_flags,
        sizeof(int),     &n,       VALUE,
        sizeof(int),     &n1,      VALUE,
        sizeof(int*),     K,           INPUT,
        sizeof(int*),     INDX,    NODEP,
        sizeof(int*),     ctot,    NODEP,
        sizeof(double*),  D,       NODEP,
        sizeof(double*),  beta,    NODEP,
        sizeof(double*),  Q,       NODEP,
        sizeof(int),     &LDQ,     VALUE,
        sizeof(double*),  Q2,      NODEP,
        sizeof(double*),  DLAMBDA, NODEP,
        sizeof(double*),  W,       NODEP,
        sizeof(double*),  Wred,       INOUT,
        sizeof(int),     &start,   VALUE,
        sizeof(int),     &end,     VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
/* Rk: This kernel is not traced, because each subkernel is already traced */
void
CORE_dlaed1_pipelined_quark(Quark *quark)
{
    int n;
    int n1;
    const int *K;
    double *D;
    const double *beta;
    double *Q;
    int LDQ;
    double *Q2;
    const double *DLAMBDA;
    const double *W;
    double *Wred;
    const int *INDX;
    const int *ctot;
    int start;
    int end;

    quark_unpack_args_15(quark, n, n1, K,
                         INDX, ctot,
                         D, beta, Q, LDQ, Q2,
                         DLAMBDA, W, Wred,
                         start, end);

    CORE_dlaed2_compressq(n, n1, INDX, ctot,
                          Q, LDQ, Q2, start, end);

    CORE_dlaed4(n, *K, D, *beta, Q, LDQ,
                DLAMBDA, W, INDX, start, end);

    CORE_dlaed3_computeW(n, *K, Q, LDQ, DLAMBDA, Wred,
                         INDX, start, end);
}

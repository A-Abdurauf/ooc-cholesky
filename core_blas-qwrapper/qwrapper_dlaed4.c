/**
 *
 * @file qwrapper_dlaed4.c
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
CORE_dlaed4_p2f1_quark(Quark *quark);

/***************************************************************************//**
 *
 **/
/*
 * Q and LDQ are passed by pointers, to avoid memory leak in Quark, we don't
 * forward LDQ which is used as K in our case.
 */
void QUARK_CORE_dlaed4_p2f1(Quark *quark, Quark_Task_Flags *task_flags,
                            int n, const int *K,
                            double *D, const double *beta, double **Q, const int *LDQ,
                            const double *DLAMBDA, const double *W, const int *INDX,
                            int start, int end,
                            PLASMA_sequence *sequence, PLASMA_request *request,
                            void *fakeQ, int flagfQ )
{
    int flagQ = NODEP;

    if (fakeQ == Q) {
        flagQ = flagfQ;
        fakeQ = NULL;
        flagfQ = NODEP;
    }

    plasma_profile_by_kernel( task_flags, LAED4 );

    QUARK_Insert_Task(quark, CORE_dlaed4_p2f1_quark, task_flags,
        sizeof(int),              &n,        VALUE,
        sizeof(int*),              K,            INPUT,
        sizeof(double*),           D,            NODEP, /* INOUT */
        sizeof(double*),           beta,         NODEP, /* INPUT */
        sizeof(double**),          Q,            flagQ,
        sizeof(double*),           DLAMBDA,      NODEP,
        sizeof(double*),           W,            NODEP,
        sizeof(int*),              INDX,         NODEP,
        sizeof(int),              &start,    VALUE,
        sizeof(int),              &end,      VALUE,
        sizeof(PLASMA_sequence*), &sequence, VALUE,
        sizeof(PLASMA_request*),  &request,  VALUE,
        /*
         * Fake dependency to guaranty correct execution of kernels running
         * on a same subset of Q Note that the kernel does not necessarly
         * work on this subset
         */
        1, fakeQ, flagfQ,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_dlaed4_p2f1_quark = PCORE_dlaed4_p2f1_quark
#define CORE_dlaed4_p2f1_quark PCORE_dlaed4_p2f1_quark
#endif
void CORE_dlaed4_p2f1_quark(Quark *quark)
{
    int n;
    const int *K;
    double *D;
    const double *beta;
    double **Q;
    const double *DLAMBDA;
    const double *W;
    const int *INDX;
    int start, end;
    PLASMA_sequence *sequence;
    PLASMA_request *request;
    double *fake;
    int info = 0;

    quark_unpack_args_13(quark, n, K, D, beta, Q,
                         DLAMBDA, W, INDX, start, end,
                         sequence, request, fake);

    info = CORE_dlaed4(n, *K, D, *beta, *Q, *K,
                       DLAMBDA, W, INDX, start, end);

    if (info != PLASMA_SUCCESS){
        plasma_sequence_flush(quark, sequence, request, info);
    }
}

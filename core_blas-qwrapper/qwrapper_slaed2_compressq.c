/**
 *
 * @file qwrapper_slaed2_copyvectors.c
 *
 *  PLASMA core_blas quark wrapper
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Grégoire Pichon
 * @author Azzam Haidar
 * @date 2014-07
 * @generated s Fri Apr  1 11:02:43 2016
 *
 **/
#include "common.h"

void
CORE_slaed2_compressq_quark(Quark *quark);

/***************************************************************************//**
 *
 **/
void QUARK_CORE_slaed2_compressq(Quark *quark, Quark_Task_Flags *task_flags,
                                 int n, int n1, int start, int end,
                                 const int *INDX, const int *ctot,
                                 const float *Q, int LDQ,
                                 float *Q2, int *K )
{
    plasma_profile_by_kernel( task_flags, LAED2_COMPRESSQ );


    /* If start == 0, no need for extra fake dependency */
    if ( start == 0 ) {
        QUARK_Insert_Task(quark, CORE_slaed2_compressq_quark, task_flags,
            sizeof(int),           &n,     VALUE,
            sizeof(int),           &n1,    VALUE,
            sizeof(int),           &start, VALUE,
            sizeof(int),           &end,   VALUE,
            sizeof(int)*n,          INDX,    NODEP,
            sizeof(int)*4,          ctot,    NODEP,
            sizeof(float)*LDQ*n,   Q,       INOUT,
            sizeof(int),           &LDQ,   VALUE,
            sizeof(float),         Q2,      INOUT | GATHERV,
            sizeof(int),            K,       INPUT,
            1,                      NULL,    NODEP,
            0);
    }
    else {
        QUARK_Insert_Task(quark, CORE_slaed2_compressq_quark, task_flags,
            sizeof(int),           &n,     VALUE,
            sizeof(int),           &n1,    VALUE,
            sizeof(int),           &start, VALUE,
            sizeof(int),           &end,   VALUE,
            sizeof(int)*n,          INDX,    NODEP,
            sizeof(int)*4,          ctot,    NODEP,
            sizeof(float)*LDQ*n,   Q,       NODEP, /* INOUT | GATHERV: Use nodep here due to submission order */
            sizeof(int),           &LDQ,   VALUE,
            sizeof(float),         Q2,      INOUT | GATHERV,
            sizeof(int),            K,       INPUT,
            /*
             * Fake dependency to guaranty correct execution of kernels running on a same subset of Q
             * Note that the kernel does not necessarly work on this subset
             */
            sizeof(float), Q+LDQ*start, INOUT,
            0);
    }
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_slaed2_compressq_quark = PCORE_slaed2_compressq_quark
#define CORE_slaed2_compressq_quark PCORE_slaed2_compressq_quark
#endif
void CORE_slaed2_compressq_quark(Quark *quark)
{
    int n, n1, start, end;
    const int *INDX;
    const int *ctot;
    const float *Q;
    int LDQ;
    float *Q2;
    int *K;
    float *fake;

    quark_unpack_args_11(quark, n, n1, start, end,
                         INDX, ctot,
                         Q, LDQ, Q2, K, fake );

    CORE_slaed2_compressq(n, n1, INDX, ctot,
                          Q, LDQ, Q2,
                          start, end);
}

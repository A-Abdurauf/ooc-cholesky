/**
 *
 * @file qwrapper_cswap.c
 *
 *  PLASMA core_blas quark wrapper
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Gregoire Pichon
 * @author Azzam Haidar
 * @date 2014-07
 * @generated c Fri Apr  1 11:02:43 2016
 *
 **/
#ifdef PLASMA_WITH_MKL
#include <mkl_lapacke.h>
#else
#include <lapacke.h>
#endif
#include "common.h"

void
CORE_cswap_quark(Quark *quark);

/***************************************************************************//**
 *
 **/
void QUARK_CORE_cswap(Quark *quark, Quark_Task_Flags *task_flags,
                      int m, int n,
                      PLASMA_Complex32_t *Q,
                      int LDQ, PLASMA_Complex32_t *work,
                      int *perm, int begin, int end)
{
    QUARK_Insert_Task(quark, CORE_cswap_quark, task_flags,
        sizeof(int),                &n,           VALUE,
        sizeof(int),                &m,           VALUE,
        sizeof(PLASMA_Complex32_t),  Q,               NODEP,
        sizeof(int),                &LDQ,         VALUE,
        sizeof(PLASMA_Complex32_t),  work,            INPUT,
        sizeof(int*),                perm,            NODEP,
        sizeof(int),                &begin,       VALUE,
        sizeof(int),                &end,         VALUE,
        sizeof(PLASMA_Complex32_t),  Q+begin*LDQ,     OUTPUT,
        0);
    /* in order to force all swap waiting the previous copy.
     * could be done as a independent pointer but to avoid
     * using a freeing task I can use work as the dependencies
     * pointer for barrier
     * */
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_cswap_quark = PCORE_cswap_quark
#define CORE_cswap_quark PCORE_cswap_quark
#endif
void CORE_cswap_quark(Quark *quark)
{
    int m;
    int n;
    PLASMA_Complex32_t *Q;
    int LDQ;
    PLASMA_Complex32_t *work;
    int *perm;
    int begin;
    int end;

    quark_unpack_args_8(quark, m, n, Q, LDQ, work, perm, begin, end);
    CORE_cswap(m, n, Q, LDQ, work, perm, begin, end);
}

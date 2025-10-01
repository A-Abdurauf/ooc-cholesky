/**
 *
 * @file qwrapper_zswap.c
 *
 *  PLASMA core_blas quark wrapper
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Gregoire Pichon
 * @author Azzam Haidar
 * @date 2014-07
 * @precisions normal z -> c d s
 *
 **/
#ifdef PLASMA_WITH_MKL
#include <mkl_lapacke.h>
#else
#include <lapacke.h>
#endif
#include "common.h"

void
CORE_zswap_quark(Quark *quark);

/***************************************************************************//**
 *
 **/
void QUARK_CORE_zswap(Quark *quark, Quark_Task_Flags *task_flags,
                      int m, int n,
                      PLASMA_Complex64_t *Q,
                      int LDQ, PLASMA_Complex64_t *work,
                      int *perm, int begin, int end)
{
    QUARK_Insert_Task(quark, CORE_zswap_quark, task_flags,
        sizeof(int),                &n,           VALUE,
        sizeof(int),                &m,           VALUE,
        sizeof(PLASMA_Complex64_t),  Q,               NODEP,
        sizeof(int),                &LDQ,         VALUE,
        sizeof(PLASMA_Complex64_t),  work,            INPUT,
        sizeof(int*),                perm,            NODEP,
        sizeof(int),                &begin,       VALUE,
        sizeof(int),                &end,         VALUE,
        sizeof(PLASMA_Complex64_t),  Q+begin*LDQ,     OUTPUT,
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
#pragma weak CORE_zswap_quark = PCORE_zswap_quark
#define CORE_zswap_quark PCORE_zswap_quark
#endif
void CORE_zswap_quark(Quark *quark)
{
    int m;
    int n;
    PLASMA_Complex64_t *Q;
    int LDQ;
    PLASMA_Complex64_t *work;
    int *perm;
    int begin;
    int end;

    quark_unpack_args_8(quark, m, n, Q, LDQ, work, perm, begin, end);
    CORE_zswap(m, n, Q, LDQ, work, perm, begin, end);
}

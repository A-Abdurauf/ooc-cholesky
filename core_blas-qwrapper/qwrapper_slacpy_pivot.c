/**
 *
 * @file qwrapper_slacpy_pivot.c
 *
 *  PLASMA core_blas quark wrapper
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Mathieu Faverge
 * @date 2013-02-01
 * @generated s Fri Apr  1 11:02:44 2016
 *
 **/
#include "common.h"

/***************************************************************************//**
 *
 **/
void QUARK_CORE_slacpy_pivot(Quark *quark, Quark_Task_Flags *task_flags,
                             const PLASMA_desc descA,
                             PLASMA_enum direct, int k1, int k2, const int *ipiv,
                             int *rankin, int *rankout,
                             float *A, int lda,
                             int pos, int init)
{
    DAG_SET_PROPERTIES( "CPY_PIV"  , "white"   );
    QUARK_Insert_Task(quark, CORE_slacpy_pivot_quark, task_flags,
        sizeof(PLASMA_desc),                    &descA,         VALUE,
        sizeof(PLASMA_enum),                    &direct,        VALUE,
        sizeof(int),                            &k1,            VALUE,
        sizeof(int),                            &k2,            VALUE,
        sizeof(int)*lda,                         ipiv,                INPUT,
        sizeof(int)*lda,                         rankin,              INOUT,
        sizeof(int)*lda,                         rankout,             OUTPUT | GATHERV,
        sizeof(float)*lda*descA.nb, A,                   INOUT | GATHERV,
        sizeof(int),                            &lda,           VALUE,
        sizeof(int),                            &pos,           VALUE,
        sizeof(int),                            &init,          VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_slacpy_pivot_quark = PCORE_slacpy_pivot_quark
#define CORE_slacpy_pivot_quark PCORE_slacpy_pivot_quark
#endif
void CORE_slacpy_pivot_quark(Quark *quark)
{
    PLASMA_desc descA;
    PLASMA_enum direct;
    float *A;
    int lda, pos, k1, k2;
    int *rankin;
    int *rankout;
    const int *ipiv;
    int init;

    quark_unpack_args_11(quark, descA, direct, k1, k2, ipiv, rankin, rankout, A, lda, pos, init);
    A = (direct == PlasmaRowwise) ? A + pos : A + (pos * lda);
    CORE_slacpy_pivot(descA, direct, k1, k2, ipiv, rankin, rankout+pos, A, lda, init );
}

/**
 *
 * @file qwrapper_cgessm.c
 *
 *  PLASMA core_blas quark wrapper
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Hatem Ltaief
 * @author Mathieu Faverge
 * @author Jakub Kurzak
 * @date 2010-11-15
 * @generated c Fri Apr  1 11:02:40 2016
 *
 **/
#include "common.h"

/***************************************************************************//**
 *
 **/
void QUARK_CORE_cgessm(Quark *quark, Quark_Task_Flags *task_flags,
                       int m, int n, int k, int ib, int nb,
                       const int *IPIV,
                       const PLASMA_Complex32_t *L, int ldl,
                       PLASMA_Complex32_t *A, int lda)
{
    plasma_profile_by_kernel( task_flags, GESSM );

    QUARK_Insert_Task(quark, CORE_cgessm_quark, task_flags,
        sizeof(int),                        &m,     VALUE,
        sizeof(int),                        &n,     VALUE,
        sizeof(int),                        &k,     VALUE,
        sizeof(int),                        &ib,    VALUE,
        sizeof(int)*nb,                      IPIV,          INPUT,
        sizeof(PLASMA_Complex32_t)*nb*nb,    L,             INPUT | QUARK_REGION_L,
        sizeof(int),                        &ldl,   VALUE,
        sizeof(PLASMA_Complex32_t)*nb*nb,    A,             INOUT,
        sizeof(int),                        &lda,   VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_cgessm_quark = PCORE_cgessm_quark
#define CORE_cgessm_quark PCORE_cgessm_quark
#endif
void CORE_cgessm_quark(Quark *quark)
{
    int m;
    int n;
    int k;
    int ib;
    int *IPIV;
    PLASMA_Complex32_t *L;
    int ldl;
    PLASMA_Complex32_t *A;
    int lda;

    quark_unpack_args_9(quark, m, n, k, ib, IPIV, L, ldl, A, lda);
    CORE_cgessm(m, n, k, ib, IPIV, L, ldl, A, lda);
}

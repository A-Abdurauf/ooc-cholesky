/**
 *
 * @file qwrapper_cssssm.c
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
void QUARK_CORE_cssssm(Quark *quark, Quark_Task_Flags *task_flags,
                       int m1, int n1, int m2, int n2, int k, int ib, int nb,
                       PLASMA_Complex32_t *A1, int lda1,
                       PLASMA_Complex32_t *A2, int lda2,
                       const PLASMA_Complex32_t *L1, int ldl1,
                       const PLASMA_Complex32_t *L2, int ldl2,
                       const int *IPIV)
{
    plasma_profile_by_kernel( task_flags, SSSSM );

    QUARK_Insert_Task(quark, CORE_cssssm_quark, task_flags,
        sizeof(int),                        &m1,    VALUE,
        sizeof(int),                        &n1,    VALUE,
        sizeof(int),                        &m2,    VALUE,
        sizeof(int),                        &n2,    VALUE,
        sizeof(int),                        &k,     VALUE,
        sizeof(int),                        &ib,    VALUE,
        sizeof(PLASMA_Complex32_t)*nb*nb,    A1,            INOUT,
        sizeof(int),                        &lda1,  VALUE,
        sizeof(PLASMA_Complex32_t)*nb*nb,    A2,            INOUT | LOCALITY,
        sizeof(int),                        &lda2,  VALUE,
        sizeof(PLASMA_Complex32_t)*ib*nb,    L1,            INPUT,
        sizeof(int),                        &ldl1,  VALUE,
        sizeof(PLASMA_Complex32_t)*ib*nb,    L2,            INPUT,
        sizeof(int),                        &ldl2,  VALUE,
        sizeof(int)*nb,                      IPIV,          INPUT,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_cssssm_quark = PCORE_cssssm_quark
#define CORE_cssssm_quark PCORE_cssssm_quark
#endif
void CORE_cssssm_quark(Quark *quark)
{
    int m1;
    int n1;
    int m2;
    int n2;
    int k;
    int ib;
    PLASMA_Complex32_t *A1;
    int lda1;
    PLASMA_Complex32_t *A2;
    int lda2;
    PLASMA_Complex32_t *L1;
    int ldl1;
    PLASMA_Complex32_t *L2;
    int ldl2;
    int *IPIV;

    quark_unpack_args_15(quark, m1, n1, m2, n2, k, ib, A1, lda1, A2, lda2, L1, ldl1, L2, ldl2, IPIV);
    CORE_cssssm(m1, n1, m2, n2, k, ib, A1, lda1, A2, lda2, L1, ldl1, L2, ldl2, IPIV);
}

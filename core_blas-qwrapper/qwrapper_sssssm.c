/**
 *
 * @file qwrapper_sssssm.c
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
 * @generated s Fri Apr  1 11:02:40 2016
 *
 **/
#include "common.h"

/***************************************************************************//**
 *
 **/
void QUARK_CORE_sssssm(Quark *quark, Quark_Task_Flags *task_flags,
                       int m1, int n1, int m2, int n2, int k, int ib, int nb,
                       float *A1, int lda1,
                       float *A2, int lda2,
                       const float *L1, int ldl1,
                       const float *L2, int ldl2,
                       const int *IPIV)
{
    plasma_profile_by_kernel( task_flags, SSSSM );

    QUARK_Insert_Task(quark, CORE_sssssm_quark, task_flags,
        sizeof(int),                        &m1,    VALUE,
        sizeof(int),                        &n1,    VALUE,
        sizeof(int),                        &m2,    VALUE,
        sizeof(int),                        &n2,    VALUE,
        sizeof(int),                        &k,     VALUE,
        sizeof(int),                        &ib,    VALUE,
        sizeof(float)*nb*nb,    A1,            INOUT,
        sizeof(int),                        &lda1,  VALUE,
        sizeof(float)*nb*nb,    A2,            INOUT | LOCALITY,
        sizeof(int),                        &lda2,  VALUE,
        sizeof(float)*ib*nb,    L1,            INPUT,
        sizeof(int),                        &ldl1,  VALUE,
        sizeof(float)*ib*nb,    L2,            INPUT,
        sizeof(int),                        &ldl2,  VALUE,
        sizeof(int)*nb,                      IPIV,          INPUT,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_sssssm_quark = PCORE_sssssm_quark
#define CORE_sssssm_quark PCORE_sssssm_quark
#endif
void CORE_sssssm_quark(Quark *quark)
{
    int m1;
    int n1;
    int m2;
    int n2;
    int k;
    int ib;
    float *A1;
    int lda1;
    float *A2;
    int lda2;
    float *L1;
    int ldl1;
    float *L2;
    int ldl2;
    int *IPIV;

    quark_unpack_args_15(quark, m1, n1, m2, n2, k, ib, A1, lda1, A2, lda2, L1, ldl1, L2, ldl2, IPIV);
    CORE_sssssm(m1, n1, m2, n2, k, ib, A1, lda1, A2, lda2, L1, ldl1, L2, ldl2, IPIV);
}

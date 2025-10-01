/**
 *
 * @file qwrapper_sttmlq.c
 *
 *  PLASMA core_blas quark wrapper
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Hatem Ltaief
 * @author Mathieu Faverge
 * @author Dulceneia Becker
 * @date 2010-11-15
 * @generated s Fri Apr  1 11:02:40 2016
 *
 **/
#include "common.h"

/***************************************************************************//**
 *
 **/
void QUARK_CORE_sttmlq(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum side, PLASMA_enum trans,
                       int m1, int n1, int m2, int n2, int k, int ib, int nb,
                       float *A1, int lda1,
                       float *A2, int lda2,
                       const float *V, int ldv,
                       const float *T, int ldt)
{
    int ldwork = side == PlasmaLeft ? ib : nb;

    plasma_profile_by_kernel( task_flags, TTMLQ );

    QUARK_Insert_Task(quark, CORE_sttmlq_quark, task_flags,
        sizeof(PLASMA_enum),                &side,  VALUE,
        sizeof(PLASMA_enum),                &trans, VALUE,
        sizeof(int),                        &m1,    VALUE,
        sizeof(int),                        &n1,    VALUE,
        sizeof(int),                        &m2,    VALUE,
        sizeof(int),                        &n2,    VALUE,
        sizeof(int),                        &k,     VALUE,
        sizeof(int),                        &ib,    VALUE,
        sizeof(float)*nb*nb,    A1,            INOUT,
        sizeof(int),                        &lda1,  VALUE,
        sizeof(float)*nb*nb,    A2,            INOUT,
        sizeof(int),                        &lda2,  VALUE,
        sizeof(float)*nb*nb,    V,             INPUT|QUARK_REGION_D|QUARK_REGION_L,
        sizeof(int),                        &ldv,   VALUE,
        sizeof(float)*ib*nb,    T,             INPUT,
        sizeof(int),                        &ldt,   VALUE,
        sizeof(float)*ib*nb,    NULL,          SCRATCH,
        sizeof(int),                        &ldwork,    VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_sttmlq_quark = PCORE_sttmlq_quark
#define CORE_sttmlq_quark PCORE_sttmlq_quark
#endif
void CORE_sttmlq_quark(Quark *quark)
{
    PLASMA_enum side;
    PLASMA_enum trans;
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
    float *V;
    int ldv;
    float *T;
    int ldt;
    float *WORK;
    int ldwork;

    quark_unpack_args_18(quark, side, trans, m1, n1, m2, n2, k, ib,
                         A1, lda1, A2, lda2, V, ldv, T, ldt, WORK, ldwork);
    CORE_sttmlq(side, trans, m1, n1, m2, n2, k, ib, A1, lda1,
                A2, lda2, V, ldv, T, ldt, WORK, ldwork);
}

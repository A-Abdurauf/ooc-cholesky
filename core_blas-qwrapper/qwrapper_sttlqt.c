/**
 *
 * @file qwrapper_sttlqt.c
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
void QUARK_CORE_sttlqt(Quark *quark, Quark_Task_Flags *task_flags,
                       int m, int n, int ib, int nb,
                       float *A1, int lda1,
                       float *A2, int lda2,
                       float *T, int ldt)
{
    plasma_profile_by_kernel( task_flags, TTLQT );

    QUARK_Insert_Task(quark, CORE_sttlqt_quark, task_flags,
        sizeof(int),                        &m,     VALUE,
        sizeof(int),                        &n,     VALUE,
        sizeof(int),                        &ib,    VALUE,
        sizeof(float)*nb*nb,    A1,            INOUT|QUARK_REGION_D|QUARK_REGION_L,
        sizeof(int),                        &lda1,  VALUE,
        sizeof(float)*nb*nb,    A2,            INOUT|QUARK_REGION_D|QUARK_REGION_L|LOCALITY,
        sizeof(int),                        &lda2,  VALUE,
        sizeof(float)*ib*nb,    T,             OUTPUT,
        sizeof(int),                        &ldt,   VALUE,
        sizeof(float)*nb,       NULL,          SCRATCH,
        sizeof(float)*ib*nb,    NULL,          SCRATCH,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_sttlqt_quark = PCORE_sttlqt_quark
#define CORE_sttlqt_quark PCORE_sttlqt_quark
#endif
void CORE_sttlqt_quark(Quark *quark)
{
    int m;
    int n;
    int ib;
    float *A1;
    int lda1;
    float *A2;
    int lda2;
    float *T;
    int ldt;
    float *TAU;
    float *WORK;

    quark_unpack_args_11(quark, m, n, ib, A1, lda1, A2, lda2, T, ldt, TAU, WORK);
    CORE_sttlqt(m, n, ib, A1, lda1, A2, lda2, T, ldt, TAU, WORK);
}

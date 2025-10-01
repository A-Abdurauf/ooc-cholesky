/**
 *
 * @file qwrapper_ctsmlq_corner.c
 *
 *  PLASMA core_blas quark wrapper
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Hatem Ltaief
 * @author Mathieu Faverge
 * @author Azzam Haidar
 * @date 2010-11-15
 * @generated c Fri Apr  1 11:02:42 2016
 *
 **/
#include "common.h"

/***************************************************************************//**
 *
 **/
void QUARK_CORE_ctsmlq_corner(Quark *quark, Quark_Task_Flags *task_flags,
                              int m1, int n1, int m2, int n2, int m3, int n3, int k, int ib, int nb,
                              PLASMA_Complex32_t *A1, int lda1,
                              PLASMA_Complex32_t *A2, int lda2,
                              PLASMA_Complex32_t *A3, int lda3,
                              const PLASMA_Complex32_t *V, int ldv,
                              const PLASMA_Complex32_t *T, int ldt)
{
    int ldwork = nb;

    plasma_profile_by_kernel( task_flags, TSMLQ );

    QUARK_Insert_Task(quark, CORE_ctsmlq_corner_quark, task_flags,
        sizeof(int),                        &m1,    VALUE,
        sizeof(int),                        &n1,    VALUE,
        sizeof(int),                        &m2,    VALUE,
        sizeof(int),                        &n2,    VALUE,
        sizeof(int),                        &m3,    VALUE,
        sizeof(int),                        &n3,    VALUE,
        sizeof(int),                        &k,     VALUE,
        sizeof(int),                        &ib,    VALUE,
        sizeof(int),                        &nb,    VALUE,
        sizeof(PLASMA_Complex32_t)*nb*nb,    A1,            INOUT|QUARK_REGION_D|QUARK_REGION_U,
        sizeof(int),                        &lda1,  VALUE,
        sizeof(PLASMA_Complex32_t)*nb*nb,    A2,            INOUT,
        sizeof(int),                        &lda2,  VALUE,
        sizeof(PLASMA_Complex32_t)*nb*nb,    A3,            INOUT|QUARK_REGION_D|QUARK_REGION_U,
        sizeof(int),                        &lda3,  VALUE,
        sizeof(PLASMA_Complex32_t)*nb*nb,    V,             INPUT,
        sizeof(int),                        &ldv,   VALUE,
        sizeof(PLASMA_Complex32_t)*ib*nb,    T,             INPUT,
        sizeof(int),                        &ldt,   VALUE,
        sizeof(PLASMA_Complex32_t)*4*nb*nb,    NULL,          SCRATCH,
        sizeof(int),                        &ldwork, VALUE,
        0);
}

/***************************************************************************//**
 * This kernel applies right and left transformations as depicted below:
 * |I -VTV'| * | A1  A2| * |I - VT'V'|
 *             | A2' A3 |
 * where A1 and A3 are symmetric matrices.
 * Only the upper part is referenced.
 * This is an adhoc implementation, can be further optimized...
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_ctsmlq_corner_quark = PCORE_ctsmlq_corner_quark
#define CORE_ctsmlq_corner_quark PCORE_ctsmlq_corner_quark
#endif
void CORE_ctsmlq_corner_quark(Quark *quark)
{
    int m1;
    int n1;
    int m2;
    int n2;
    int m3;
    int n3;
    int k;
    int ib;
    int nb;
    PLASMA_Complex32_t *A1;
    int lda1;
    PLASMA_Complex32_t *A2;
    int lda2;
    PLASMA_Complex32_t *A3;
    int lda3;
    PLASMA_Complex32_t *V;
    int ldv;
    PLASMA_Complex32_t *T;
    int ldt;
    PLASMA_Complex32_t *WORK;
    int ldwork;

    quark_unpack_args_21(quark, m1, n1, m2, n2, m3, n3, k, ib, nb,
                         A1, lda1, A2, lda2, A3, lda3, V, ldv, T, ldt, WORK, ldwork);

    CORE_ctsmlq_corner(m1, n1, m2, n2, m3, n3, k, ib, nb,
                       A1, lda1, A2, lda2, A3, lda3, V, ldv, T, ldt, WORK, ldwork);

}

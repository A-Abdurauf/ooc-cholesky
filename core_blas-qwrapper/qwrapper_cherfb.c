/**
 *
 * @file qwrapper_cherfb.c
 *
 *  PLASMA core_blas quark wrapper
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Hatem Ltaief
 * @date 2010-11-15
 * @generated c Fri Apr  1 11:02:42 2016
 *
 **/
#include "common.h"

/***************************************************************************//**
 * This kernel is just a workaround for now... will be deleted eventually
 * and replaced by the one above (Piotr's Task)
 **/
void QUARK_CORE_cherfb(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum uplo,
                       int n, int k, int ib, int nb,
                       const PLASMA_Complex32_t *A, int lda,
                       const PLASMA_Complex32_t *T, int ldt,
                       PLASMA_Complex32_t *C, int ldc)
{
    /* TODO: Understand why A needs to be INOUT and not INPUT */
    plasma_profile_by_kernel( task_flags, HERFB );

    QUARK_Insert_Task(
        quark, CORE_cherfb_quark, task_flags,
        sizeof(PLASMA_enum),                     &uplo,  VALUE,
        sizeof(int),                             &n,     VALUE,
        sizeof(int),                             &k,     VALUE,
        sizeof(int),                             &ib,    VALUE,
        sizeof(int),                             &nb,    VALUE,
        sizeof(PLASMA_Complex32_t)*nb*nb,        A,          uplo == PlasmaUpper ? INOUT|QUARK_REGION_U : INOUT|QUARK_REGION_L,
        sizeof(int),                             &lda,   VALUE,
        sizeof(PLASMA_Complex32_t)*ib*nb,        T,          INPUT,
        sizeof(int),                             &ldt,   VALUE,
        sizeof(PLASMA_Complex32_t)*nb*nb,        C,          uplo == PlasmaUpper ? INOUT|QUARK_REGION_D|QUARK_REGION_U : INOUT|QUARK_REGION_D|QUARK_REGION_L,
        sizeof(int),                             &ldc,   VALUE,
        sizeof(PLASMA_Complex32_t)*2*nb*nb,    NULL,         SCRATCH,
        sizeof(int),                             &nb,    VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_cherfb_quark = PCORE_cherfb_quark
#define CORE_cherfb_quark PCORE_cherfb_quark
#endif
void CORE_cherfb_quark(Quark *quark)
{
    PLASMA_enum uplo;
    int n;
    int k;
    int ib;
    int nb;
    PLASMA_Complex32_t *A;
    int lda;
    PLASMA_Complex32_t *T;
    int ldt;
    PLASMA_Complex32_t *C;
    int ldc;
    PLASMA_Complex32_t *WORK;
    int ldwork;

    quark_unpack_args_13(quark, uplo, n, k, ib, nb, A, lda, T, ldt, C, ldc, WORK, ldwork);
    CORE_cherfb(uplo, n, k, ib, nb, A, lda, T, ldt, C, ldc, WORK, ldwork);
}

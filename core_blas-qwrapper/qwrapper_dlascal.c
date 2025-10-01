/**
 *
 * @file qwrapper_dlascal.c
 *
 *  PLASMA core_blas quark wrapper
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Mathieu Faverge
 * @date 2015-11-05
 * @generated d Fri Apr  1 11:02:43 2016
 *
 **/
#include "common.h"

void
CORE_dlascal_quark(Quark *quark);

/***************************************************************************//**
 *
 **/
void QUARK_CORE_dlascal(Quark *quark, Quark_Task_Flags *task_flags,
                        PLASMA_enum uplo, int m, int n, int nb,
                        double alpha, double *A, int lda)
{
    plasma_profile_by_kernel( task_flags, LASCL );

    QUARK_Insert_Task(quark, CORE_dlascal_quark, task_flags,
                      sizeof(PLASMA_enum),                &uplo,      VALUE,
                      sizeof(int),                        &m,         VALUE,
                      sizeof(int),                        &n,         VALUE,
                      sizeof(double),         &alpha,     VALUE,
                      sizeof(double)*lda*nb,   A,         INOUT,
                      sizeof(int),                        &lda,       VALUE,
                      0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_dlascal_quark = PCORE_dlascal_quark
#define CORE_dlascal_quark PCORE_dlascal_quark
#endif
void CORE_dlascal_quark(Quark *quark)
{
    PLASMA_enum uplo;
    int m, n;
    double alpha;
    double *A;
    int lda;

    quark_unpack_args_6(quark, uplo, m, n, alpha, A, lda);
    CORE_dlascal(uplo, m, n, alpha, A, lda);
}

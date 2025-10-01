/**
 *
 * @file qwrapper_dlaset.c
 *
 *  PLASMA core_blas quark wrapper
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Hatem Ltaief
 * @date 2010-11-15
 * @generated d Fri Apr  1 11:02:42 2016
 *
 **/
#ifdef PLASMA_WITH_MKL
#include <mkl_lapacke.h>
#else
#include <lapacke.h>
#endif
#include "common.h"

/***************************************************************************//**
 *
 **/
void QUARK_CORE_dlaset(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum uplo, int M, int N,
                       double alpha, double beta,
                       double *A, int LDA)
{
    plasma_profile_by_kernel( task_flags, LASET );

    QUARK_Insert_Task(quark, CORE_dlaset_quark, task_flags,
        sizeof(PLASMA_enum),                &uplo,  VALUE,
        sizeof(int),                        &M,     VALUE,
        sizeof(int),                        &N,     VALUE,
        sizeof(double),         &alpha, VALUE,
        sizeof(double),         &beta,  VALUE,
        sizeof(double)*LDA*N,    A,      OUTPUT,
        sizeof(int),                        &LDA,   VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_dlaset_quark = PCORE_dlaset_quark
#define CORE_dlaset_quark PCORE_dlaset_quark
#endif
void CORE_dlaset_quark(Quark *quark)
{
    PLASMA_enum uplo;
    int M;
    int N;
    double alpha;
    double beta;
    double *A;
    int LDA;

    quark_unpack_args_7(quark, uplo, M, N, alpha, beta, A, LDA);
    LAPACKE_dlaset_work(
        LAPACK_COL_MAJOR,
        lapack_const(uplo),
        M, N, alpha, beta, A, LDA);
}

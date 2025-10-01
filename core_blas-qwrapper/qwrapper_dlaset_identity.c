/**
 *
 * @file qwrapper_dlaset_identity.c
 *
 *  PLASMA core_blas quark wrapper
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Gr�goire Pichon
 * @date 2014-08
 * @generated d Fri Apr  1 11:02:43 2016
 *
 **/
#ifdef PLASMA_WITH_MKL
#include <mkl_lapacke.h>
#else
#include <lapacke.h>
#endif
#include "common.h"

void
CORE_dlaset_identity_quark(Quark *quark);
/***************************************************************************//**
 *
 **/
void QUARK_CORE_dlaset_identity(Quark *quark, Quark_Task_Flags *task_flags,
                                int n, int start, int size,
                                double *A)
{
    plasma_profile_by_kernel( task_flags, LASET );

    QUARK_Insert_Task(quark, CORE_dlaset_identity_quark, task_flags,
        sizeof(int),                        &n,     VALUE,
        sizeof(int),                        &start, VALUE,
        sizeof(int),                        &size, VALUE,
        sizeof(double*),        A,     NODEP,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_dlaset_identity_quark = PCORE_dlaset_identity_quark
#define CORE_dlaset_identity_quark PCORE_dlaset_identity_quark
#endif
void CORE_dlaset_identity_quark(Quark *quark)
{
    int n;
    int start;
    int size;
    double *A;

    quark_unpack_args_4(quark, n, start, size, A);

    memset(A+n*start, 0, n*size*sizeof(double));

    int j;
    for (j=start; j<size+start; j++){
        A[n*j+j] = 1.0;
    }
}

/**
 *
 * @file qwrapper_zlaset_identity.c
 *
 *  PLASMA core_blas quark wrapper
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Gr�goire Pichon
 * @date 2014-08
 * @precisions normal z -> c d s
 *
 **/
#ifdef PLASMA_WITH_MKL
#include <mkl_lapacke.h>
#else
#include <lapacke.h>
#endif
#include "common.h"

void
CORE_zlaset_identity_quark(Quark *quark);
/***************************************************************************//**
 *
 **/
void QUARK_CORE_zlaset_identity(Quark *quark, Quark_Task_Flags *task_flags,
                                int n, int start, int size,
                                PLASMA_Complex64_t *A)
{
    plasma_profile_by_kernel( task_flags, LASET );

    QUARK_Insert_Task(quark, CORE_zlaset_identity_quark, task_flags,
        sizeof(int),                        &n,     VALUE,
        sizeof(int),                        &start, VALUE,
        sizeof(int),                        &size, VALUE,
        sizeof(PLASMA_Complex64_t*),        A,     NODEP,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_zlaset_identity_quark = PCORE_zlaset_identity_quark
#define CORE_zlaset_identity_quark PCORE_zlaset_identity_quark
#endif
void CORE_zlaset_identity_quark(Quark *quark)
{
    int n;
    int start;
    int size;
    PLASMA_Complex64_t *A;

    quark_unpack_args_4(quark, n, start, size, A);

    memset(A+n*start, 0, n*size*sizeof(PLASMA_Complex64_t));

    int j;
    for (j=start; j<size+start; j++){
        A[n*j+j] = 1.0;
    }
}

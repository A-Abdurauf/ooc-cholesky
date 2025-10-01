/**
 *
 * @file qwrapper_zlascal.c
 *
 *  PLASMA core_blas quark wrapper
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Mathieu Faverge
 * @date 2015-11-05
 * @precisions normal z -> c d s
 *
 **/
#include "common.h"

void
CORE_zlascal_quark(Quark *quark);

/***************************************************************************//**
 *
 **/
void QUARK_CORE_zlascal(Quark *quark, Quark_Task_Flags *task_flags,
                        PLASMA_enum uplo, int m, int n, int nb,
                        PLASMA_Complex64_t alpha, PLASMA_Complex64_t *A, int lda)
{
    plasma_profile_by_kernel( task_flags, LASCL );

    QUARK_Insert_Task(quark, CORE_zlascal_quark, task_flags,
                      sizeof(PLASMA_enum),                &uplo,      VALUE,
                      sizeof(int),                        &m,         VALUE,
                      sizeof(int),                        &n,         VALUE,
                      sizeof(PLASMA_Complex64_t),         &alpha,     VALUE,
                      sizeof(PLASMA_Complex64_t)*lda*nb,   A,         INOUT,
                      sizeof(int),                        &lda,       VALUE,
                      0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_zlascal_quark = PCORE_zlascal_quark
#define CORE_zlascal_quark PCORE_zlascal_quark
#endif
void CORE_zlascal_quark(Quark *quark)
{
    PLASMA_enum uplo;
    int m, n;
    PLASMA_Complex64_t alpha;
    PLASMA_Complex64_t *A;
    int lda;

    quark_unpack_args_6(quark, uplo, m, n, alpha, A, lda);
    CORE_zlascal(uplo, m, n, alpha, A, lda);
}

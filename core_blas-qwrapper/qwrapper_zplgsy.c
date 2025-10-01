/**
 *
 * @file qwrapper_zplgsy.c
 *
 *  PLASMA core_blas quark wrapper
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Piotr Luszczek
 * @author Pierre Lemarinier
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @precisions normal z -> c d s
 *
 **/
#include "common.h"

/***************************************************************************//**
 *
 **/
void QUARK_CORE_zplgsy( Quark *quark, Quark_Task_Flags *task_flags,
                        PLASMA_Complex64_t bump, int m, int n, PLASMA_Complex64_t *A, int lda,
                        int bigM, int m0, int n0, unsigned long long int seed )
{
    plasma_profile_by_kernel( task_flags, PLGSY );

    QUARK_Insert_Task(quark, CORE_zplgsy_quark, task_flags,
        sizeof(PLASMA_Complex64_t),       &bump, VALUE,
        sizeof(int),                      &m,    VALUE,
        sizeof(int),                      &n,    VALUE,
        sizeof(PLASMA_Complex64_t)*lda*n, A,         OUTPUT,
        sizeof(int),                      &lda,  VALUE,
        sizeof(int),                      &bigM, VALUE,
        sizeof(int),                      &m0,   VALUE,
        sizeof(int),                      &n0,   VALUE,
        sizeof(unsigned long long int),   &seed, VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_zplgsy_quark = PCORE_zplgsy_quark
#define CORE_zplgsy_quark PCORE_zplgsy_quark
#endif
void CORE_zplgsy_quark(Quark *quark)
{
    PLASMA_Complex64_t bump;
    int m;
    int n;
    PLASMA_Complex64_t *A;
    int lda;
    int bigM;
    int m0;
    int n0;
    unsigned long long int seed;

    quark_unpack_args_9( quark, bump, m, n, A, lda, bigM, m0, n0, seed );
    CORE_zplgsy( bump, m, n, A, lda, bigM, m0, n0, seed );
}


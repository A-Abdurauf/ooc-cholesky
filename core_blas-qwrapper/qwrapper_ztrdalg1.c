/**
 *
 * @file qwrapper_ztrdalg1.c
 *
 *  PLASMA core_blas quark wrapper
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Azzam Haidar
 * @date 2013-07-04
 * @precisions normal z -> c d s
 *
 **/
#include "common.h"

/***************************************************************************//**
 *
 **/
void QUARK_CORE_ztrdalg1(Quark *quark, Quark_Task_Flags *task_flags,
                         int n,
                         int nb,
                         PLASMA_Complex64_t *A,
                         int lda,
                         PLASMA_Complex64_t *V,
                         PLASMA_Complex64_t *TAU,
                         int Vblksiz, int wantz,
                         int i, int sweepid, int m, int grsiz,
                         int *PCOL, int *ACOL, int *MCOL)
{
    QUARK_Insert_Task(quark, CORE_ztrdalg1_quark,   task_flags,
        sizeof(int),                      &n, VALUE,
        sizeof(int),                     &nb, VALUE,
        sizeof(PLASMA_Complex64_t),        A,   NODEP,
        sizeof(int),                    &lda, VALUE,
        sizeof(PLASMA_Complex64_t),        V,   NODEP,
        sizeof(PLASMA_Complex64_t),      TAU,   NODEP,
        sizeof(int),                &Vblksiz, VALUE,
        sizeof(int),                  &wantz, VALUE,
        sizeof(int),                      &i, VALUE,
        sizeof(int),                &sweepid, VALUE,
        sizeof(int),                      &m, VALUE,
        sizeof(int),                  &grsiz, VALUE,
        sizeof(PLASMA_Complex64_t)*nb,  NULL,   SCRATCH,
        sizeof(int),                    PCOL, INPUT,
        sizeof(int),                    ACOL, INPUT,
        sizeof(int),                    MCOL, OUTPUT | LOCALITY,
        0);

}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_ztrdalg1_quark = PCORE_ztrdalg1_quark
#define CORE_ztrdalg1_quark PCORE_ztrdalg1_quark
#endif
void CORE_ztrdalg1_quark(Quark *quark)
{
    int n, nb, lda, Vblksiz, wantz, i, sweepid, m, grsiz;
    PLASMA_Complex64_t *A, *V, *TAU, *work;

    quark_unpack_args_13(quark, n, nb, A, lda, V, TAU, Vblksiz, wantz, i, sweepid, m, grsiz, work);
    CORE_ztrdalg1(n, nb, A, lda, V, TAU, Vblksiz, wantz, i, sweepid, m, grsiz, work);
}

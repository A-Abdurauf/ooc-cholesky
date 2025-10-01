/**
 *
 * @file qwrapper_strdalg1.c
 *
 *  PLASMA core_blas quark wrapper
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Azzam Haidar
 * @date 2013-07-04
 * @generated s Fri Apr  1 11:02:44 2016
 *
 **/
#include "common.h"

/***************************************************************************//**
 *
 **/
void QUARK_CORE_strdalg1(Quark *quark, Quark_Task_Flags *task_flags,
                         int n,
                         int nb,
                         float *A,
                         int lda,
                         float *V,
                         float *TAU,
                         int Vblksiz, int wantz,
                         int i, int sweepid, int m, int grsiz,
                         int *PCOL, int *ACOL, int *MCOL)
{
    QUARK_Insert_Task(quark, CORE_strdalg1_quark,   task_flags,
        sizeof(int),                      &n, VALUE,
        sizeof(int),                     &nb, VALUE,
        sizeof(float),        A,   NODEP,
        sizeof(int),                    &lda, VALUE,
        sizeof(float),        V,   NODEP,
        sizeof(float),      TAU,   NODEP,
        sizeof(int),                &Vblksiz, VALUE,
        sizeof(int),                  &wantz, VALUE,
        sizeof(int),                      &i, VALUE,
        sizeof(int),                &sweepid, VALUE,
        sizeof(int),                      &m, VALUE,
        sizeof(int),                  &grsiz, VALUE,
        sizeof(float)*nb,  NULL,   SCRATCH,
        sizeof(int),                    PCOL, INPUT,
        sizeof(int),                    ACOL, INPUT,
        sizeof(int),                    MCOL, OUTPUT | LOCALITY,
        0);

}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_strdalg1_quark = PCORE_strdalg1_quark
#define CORE_strdalg1_quark PCORE_strdalg1_quark
#endif
void CORE_strdalg1_quark(Quark *quark)
{
    int n, nb, lda, Vblksiz, wantz, i, sweepid, m, grsiz;
    float *A, *V, *TAU, *work;

    quark_unpack_args_13(quark, n, nb, A, lda, V, TAU, Vblksiz, wantz, i, sweepid, m, grsiz, work);
    CORE_strdalg1(n, nb, A, lda, V, TAU, Vblksiz, wantz, i, sweepid, m, grsiz, work);
}

/**
 *
 * @file qwrapper_sbrdalg1.c
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
void QUARK_CORE_sbrdalg1(Quark *quark, Quark_Task_Flags *task_flags,
                         PLASMA_enum uplo, int n, int nb,
                         float *A, int lda,
                         float *VQ, float *TAUQ,
                         float *VP, float *TAUP,
                         int Vblksiz, int wantz,
                         int i, int sweepid, int m, int grsiz,
                         int *PCOL, int *ACOL, int *MCOL)
{
    QUARK_Insert_Task(quark, CORE_sbrdalg1_quark, task_flags,
        sizeof(int),                   &uplo, VALUE,
        sizeof(int),                      &n, VALUE,
        sizeof(int),                     &nb, VALUE,
        sizeof(float),        A,    NODEP,
        sizeof(int),                    &lda, VALUE,
        sizeof(float),       VQ,    NODEP,
        sizeof(float),     TAUQ,    NODEP,
        sizeof(float),       VP,    NODEP,
        sizeof(float),     TAUP,    NODEP,
        sizeof(int),                &Vblksiz, VALUE,
        sizeof(int),                  &wantz, VALUE,
        sizeof(int),                      &i, VALUE,
        sizeof(int),                &sweepid, VALUE,
        sizeof(int),                      &m, VALUE,
        sizeof(int),                  &grsiz, VALUE,
        sizeof(float)*nb,  NULL,    SCRATCH,
        sizeof(int),                    PCOL,    INPUT,
        sizeof(int),                    ACOL,    INPUT,
        sizeof(int),                    MCOL,    OUTPUT | LOCALITY,
        0);

}
/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_sbrdalg1_quark = PCORE_sbrdalg1_quark
#define CORE_sbrdalg1_quark PCORE_sbrdalg1_quark
#endif
void CORE_sbrdalg1_quark(Quark *quark)
{
    PLASMA_enum uplo;
    int n, nb, lda, Vblksiz, wantz, i, sweepid, m, grsiz;
    float *A, *VQ, *TAUQ, *VP, *TAUP, *work;

    quark_unpack_args_16(quark, uplo, n, nb, A, lda, VQ, TAUQ, VP, TAUP,
                         Vblksiz, wantz, i, sweepid, m, grsiz, work);
    CORE_sbrdalg1(uplo, n, nb, A, lda, VQ, TAUQ, VP, TAUP,
                  Vblksiz, wantz, i, sweepid, m, grsiz, work);
}

/**
 *
 * @file qwrapper_zbrdalg1.c
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
void QUARK_CORE_zbrdalg1(Quark *quark, Quark_Task_Flags *task_flags,
                         PLASMA_enum uplo, int n, int nb,
                         PLASMA_Complex64_t *A, int lda,
                         PLASMA_Complex64_t *VQ, PLASMA_Complex64_t *TAUQ,
                         PLASMA_Complex64_t *VP, PLASMA_Complex64_t *TAUP,
                         int Vblksiz, int wantz,
                         int i, int sweepid, int m, int grsiz,
                         int *PCOL, int *ACOL, int *MCOL)
{
    QUARK_Insert_Task(quark, CORE_zbrdalg1_quark, task_flags,
        sizeof(int),                   &uplo, VALUE,
        sizeof(int),                      &n, VALUE,
        sizeof(int),                     &nb, VALUE,
        sizeof(PLASMA_Complex64_t),        A,    NODEP,
        sizeof(int),                    &lda, VALUE,
        sizeof(PLASMA_Complex64_t),       VQ,    NODEP,
        sizeof(PLASMA_Complex64_t),     TAUQ,    NODEP,
        sizeof(PLASMA_Complex64_t),       VP,    NODEP,
        sizeof(PLASMA_Complex64_t),     TAUP,    NODEP,
        sizeof(int),                &Vblksiz, VALUE,
        sizeof(int),                  &wantz, VALUE,
        sizeof(int),                      &i, VALUE,
        sizeof(int),                &sweepid, VALUE,
        sizeof(int),                      &m, VALUE,
        sizeof(int),                  &grsiz, VALUE,
        sizeof(PLASMA_Complex64_t)*nb,  NULL,    SCRATCH,
        sizeof(int),                    PCOL,    INPUT,
        sizeof(int),                    ACOL,    INPUT,
        sizeof(int),                    MCOL,    OUTPUT | LOCALITY,
        0);

}
/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_zbrdalg1_quark = PCORE_zbrdalg1_quark
#define CORE_zbrdalg1_quark PCORE_zbrdalg1_quark
#endif
void CORE_zbrdalg1_quark(Quark *quark)
{
    PLASMA_enum uplo;
    int n, nb, lda, Vblksiz, wantz, i, sweepid, m, grsiz;
    PLASMA_Complex64_t *A, *VQ, *TAUQ, *VP, *TAUP, *work;

    quark_unpack_args_16(quark, uplo, n, nb, A, lda, VQ, TAUQ, VP, TAUP,
                         Vblksiz, wantz, i, sweepid, m, grsiz, work);
    CORE_zbrdalg1(uplo, n, nb, A, lda, VQ, TAUQ, VP, TAUP,
                  Vblksiz, wantz, i, sweepid, m, grsiz, work);
}

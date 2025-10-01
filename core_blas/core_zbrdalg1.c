/**
 *
 * @file core_zbrdalg1.c
 *
 *  PLASMA core_blas kernel
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Azzam Haidar
 * @date 2013-07-04
 * @precisions normal z -> c d s
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
 * @ingroup CORE_PLASMA_Complex64_t
 *
 *  CORE_zbrdalg1 is a part of the bidiagonal reduction algorithm (bulgechasing).
 *  It correspond to a local driver of the kernels that should be executed on a
 *  single core.
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *         @arg PlasmaLower:
 *         @arg PlasmaUpper:
 *
 * @param[in] n
 *          The order of the matrix A. n >= 0.
 *
 * @param[in] nb
 *          The size of the Bandwidth of the matrix A,
 *          which correspond to the tile size. nb >= 0.
 *
 * @param[in,out] A
 *          PLASMA_Complex64_t array, dimension (lda,n)
 *          On entry, the (2nb+1)-by-n lower or upper band  general
 *          matrix to be reduced to bidiagonal.
 *          On exit, if uplo = PlasmaUpper, the diagonal and first
 *          superdiagonal of A are overwritten by the corresponding
 *          elements of the bidiagonal matrix B.
 *          if uplo = PlasmaLower the diagonal and first subdiagonal
 *          of A are overwritten by the corresponding elements of the
 *          elements of the bidiagonal matrix B.
 *
 * @param[in] lda
 *          (input) INTEGER
 *          The leading dimension of the array A.  LDA >= max(1,nb+1).
 *
 * @param[out] VQ
 *          PLASMA_Complex64_t array, dimension (n) if wantz=0
 *          or ldv*Vblksiz*blkcnt if wantz>0.
 *          The scalar elementary left reflectors are written in
 *          this array.
 *
 * @param[out] TAUQ
 *          PLASMA_Complex64_t array, dimension (n) if wantz=0
 *          or Vblksiz*Vblksiz*blkcnt if wantz>0.
 *          The scalar factors of the left elementary reflectors
 *          are written in this array.
 *
 * @param[in] VP
 *          PLASMA_Complex64_t array, dimension (n) if wantz=0
 *          or ldv*Vblksiz*blkcnt if wantz>0.
 *          The scalar elementary right reflectors are written in
 *          this array.
 *
 * @param[in] TAUP
 *          PLASMA_Complex64_t array, dimension (n) if wantz=0
 *          or Vblksiz*Vblksiz*blkcnt if wantz>0.
 *          The scalar factors of the right elementary reflectors
 *          are written in this array.
 *
 * @param[in] Vblksiz
 *          Local parameter to Plasma. It correspond to the local bloccking
 *          of the applyQ2 used to apply the orthogonal matrix Q2.
 *
 * @param[in] wantz
 *          integer tobe 0 or 1. if wantz=0 the V and TAU are not stored on
 *          only they are kept for next step then overwritten.
 *
 * @param[in] i
 *          Integer that refer to the current sweep. (outer loop).
 *
 * @param[in] sweepid
 *          Integer that refer to the sweep to chase.(inner loop).
 *
 * @param[in] m
 *          Integer that refer to a sweep step, to ensure order dependencies.
 *
 * @param[in] grsiz
 *          Integer that refer to the size of a group.
 *          group mean the number of kernel that should be executed sequentially
 *          on the same core.
 *          group size is a trade-off between locality (cache reuse) and parallelism.
 *          a small group size increase parallelism while a large group size increase
 *          cache reuse.
 *
 * @param[in] work
 *          Workspace of size nb. Used by the core_zgbtype[123]cb.
 *
 ******************************************************************************/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_zbrdalg1 = PCORE_zbrdalg1
#define CORE_zbrdalg1 PCORE_zbrdalg1
#endif
void CORE_zbrdalg1( PLASMA_enum uplo, int n, int nb,
                    PLASMA_Complex64_t *A, int lda,
                    PLASMA_Complex64_t *VQ, PLASMA_Complex64_t *TAUQ,
                    PLASMA_Complex64_t *VP, PLASMA_Complex64_t *TAUP,
                    int Vblksiz, int wantz,
                    int i, int sweepid, int m, int grsiz,
                    PLASMA_Complex64_t *work)
{
    int    k, shift=3;
    int    myid, colpt, stind, edind, blklastind, stepercol;

    k = shift / grsiz;
    stepercol = (k*grsiz == shift) ? k : k+1;
    for (k = 0; k < grsiz; k++){
        myid = (i-sweepid)*(stepercol*grsiz) +(m-1)*grsiz + k+1;
        if(myid%2 ==0) {
            colpt      = (myid/2) * nb + 1 + sweepid - 1;
            stind      = colpt - nb + 1;
            edind      = min(colpt, n);
            blklastind = colpt;
        } else {
            colpt      = ((myid+1)/2)*nb + 1 +sweepid -1 ;
            stind      = colpt-nb+1;
            edind      = min(colpt,n);
            if( (stind>=edind-1) && (edind==n) )
                blklastind = n;
            else
                blklastind = 0;
        }

        if( myid == 1 ){
           CORE_zgbtype1cb(uplo, n, nb, A, lda, VQ, TAUQ, VP, TAUP, stind-1, edind-1, sweepid-1,  Vblksiz, wantz, work);
        }else if(myid%2 == 0){
           CORE_zgbtype2cb(uplo, n, nb, A, lda, VQ, TAUQ, VP, TAUP, stind-1, edind-1, sweepid-1, Vblksiz, wantz, work);
        }else{ /*if(myid%2 == 1)*/
           CORE_zgbtype3cb(uplo, n, nb, A, lda, VQ, TAUQ, VP, TAUP, stind-1, edind-1, sweepid-1, Vblksiz, wantz, work);
        }

        if(blklastind >= (n-1))  break;
    }
}

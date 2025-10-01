/**
 *
 * @file core_zgeqp3_norms.c
 *
 *  PLASMA core_blas kernel
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Mark Gates
 * @date 2010-11-15
 * @precisions normal z -> c d s
 *
 **/
#ifdef PLASMA_WITH_MKL
#include <mkl_lapacke.h>
#else
#include <lapacke.h>
#endif
#include <math.h>
#include <cblas.h>
#include "common.h"

#define A(m,n) BLKADDR( A, PLASMA_Complex64_t, m, n )

/***************************************************************************//**
 *
 * @ingroup CORE_PLASMA_Complex64_t
 *
 *  CORE_zgeqp3_norms computes the 2-norm of each column of A[ ioff:m, joff:n ]
 *  that is marked with norms2[j] == -1 on entry. Entries that are not marked
 *  are assumed to already contain the correct 2-norm, so that the same routine
 *  can be used for computing the initial norms and for updating bad norms.
 *  The result is stored duplicated in norms1 and norms2.
 *
 *******************************************************************************
 *
 *  @param[in] A
 *          PLASMA descriptor of the matrix A.
 *          On entry, the M-by-N matrix described by the descriptor.
 *
 *  @param[in] ioff
 *          Row offset.
 *
 *  @param[in] joff
 *          Column offset.
 *
 *  @param[in,out] norms1
 *          Vector of size A.n.
 *          On exit, norms1[j] is 2-norm of column j, for j >= joff.
 *
 *  @param[in,out] norms2
 *          Vector of size A.n.
 *          On entry, if norms2[j] == -1, re-compute norm of column j.
 *          On exit, norms2[j] is 2-norm of column j, for j >= joff.
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_zgeqp3_norms = PCORE_zgeqp3_norms
#define CORE_zgeqp3_norms PCORE_zgeqp3_norms
#define CORE_zgessq PCORE_zgessq
int
CORE_zgessq(int M, int N,
            const PLASMA_Complex64_t *A, int LDA,
            double *scale, double *sumsq);
#endif
void CORE_zgeqp3_norms( PLASMA_desc A, int ioff, int joff, double *norms1, double *norms2 )
{
    const PLASMA_Complex64_t *Ai;
    int j, ii, ioff2, len, mb, nb, lda;
    double sumsq, scale;

    if ( A.nt != 1 ) {
        coreblas_error(1, "Illegal value of A.nt");
        return;
    }

    nb = min( A.nb, A.n );
    for( j = joff; j < nb; ++j ) {
        if ( norms2[j] == -1. ) {
            scale = 0.;
            sumsq = 1.;
            ioff2 = ioff;
            for( ii = 0; ii < A.mt; ++ii ) {
                mb = min( A.mb, A.m - ii*A.mb );
                Ai = A(ii,0);
                lda = BLKLDD( A, ii );
                len = mb - ioff2;
                CORE_zgessq( len, 1, Ai + j*lda + ioff2, lda, &scale, &sumsq );
                ioff2 = 0;
            }
            norms2[j] = scale * sqrt( sumsq );
            norms1[j] = norms2[j];
        }
    }
}

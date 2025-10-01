/**
 *
 * @file core_dgeqp3_larfg.c
 *
 *  PLASMA core_blas kernel
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Mark Gates
 * @date 2010-11-15
 * @generated d Fri Apr  1 11:02:31 2016
 *
 **/
#include <math.h>
#include "common.h"

#define A(m,n) BLKADDR( A, double, m, n )

/***************************************************************************//**
 *
 * @ingroup CORE_double
 *
 * CORE_dgeqp3_larfg generates a Householder elementary reflector H, such that
 *
 *     H**T * x = [ beta ]  and  H**T * H = I.
 *                [ 0    ]
 *
 * where alpha and beta are scalars, with beta real, and x is an n element vector.
 * H is reperested in the form
 *
 *     H = I - tau * [ 1 ] * [ 1 v**T ],
 *                   [ v ]
 *
 * where tau is a scalar and v is an (n-1) element vector.
 * If x[1:] = 0 and x[0] is real, then tau = 0 and H = I.
 * Otherwise, 1 <= real(tau) <= 2 and abs(tau-1) <= 1.
 *
 * Here, x = A[ ii*mb + i : m, jj*nb + j ].
 * That is, x is j-th column of the jj-th block-column of A, starting in the
 * i-th row of the ii-th block-row of A, and going to the last row.
 * Note that x spans multiple tiles of A.
 *
 * This DIFFERS from LAPACK in that the 1.0 is stored explicitly in the top
 * element of x and beta is stored separately. (Whereas in LAPACK, the
 * 1.0 is implicit and beta is stored in the top element of x.)
 *
 *******************************************************************************
 *
 * @param[in,out] A
 *         Description of m by n matrix A.
 *         On entry, A[ ii*mb + i : m, jj*nb + j ] is the vector x.
 *         On exit,  A[ ii*mb + i : m, jj*nb + j ] is overwritten with [ 1, v ].
 *
 * @param[in] ii
 *         Index of block row of A to start in.
 *
 * @param[in] jj
 *         Index of block column of A.
 *
 * @param[in] i
 *         Index of row within ii-th block row to start in.
 *
 * @param[in] j
 *         Index of column within jj-th block column.
 *
 * @param[out] tau
 *         The scalar tau.
 *
 * @param[out] beta
 *         The scalar beta.
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_dgeqp3_larfg = PCORE_dgeqp3_larfg
#define CORE_dgeqp3_larfg PCORE_dgeqp3_larfg
#endif
void CORE_dgeqp3_larfg( PLASMA_desc A, int ii, int jj, int i, int j,
                        double *tau, double *beta )
{
    int i2, kk, k, mb, lda;
    double norm2;
    double x0, scale;
    double *Akj, *Aij;

    /* TODO: it might be simpler to pass in global (i,j) indices, */
    /* from which (ii,jj) can be derived. */

    /* compute norm( A[i+1:m,j] )^2 */
    /* todo: use lassq algorithm */
    i2 = i+1;
    norm2 = 0.;
    for( kk = ii; kk < A.mt; ++kk ) {
        mb  = min( A.mb, A.m - kk*A.mb );
        lda = BLKLDD( A, kk );
        Akj = A(kk,jj);
        for( k = i2; k < mb; ++k ) {
            norm2 += Akj[k + j*lda] * ( Akj[k + j*lda] );
        }
        i2 = 0;
    }

    lda = BLKLDD( A, ii );
    Aij = A(ii,jj);
    x0 = Aij[i + j*lda];
    if ( norm2 == 0. && cimag(x0) == 0. ) {
        /* H = I */
        *tau  = 0;
        *beta = x0;
    }
    else {
        /* todo: use lapack's fancy sqrt function */
        /* todo: scale vector as in lapack */
        *beta = sqrt( (x0)*(x0) + cimag(x0)*cimag(x0) + norm2 );
        if ( (x0) >= 0. ) {
            *beta = -(*beta);
        }
        *tau  = (*beta - x0) / (*beta);
        /* todo: use zladiv or dladiv */
        scale = 1. / (x0 - *beta);

        /* x *= scale */
        i2 = i;
        for( kk = ii; kk < A.mt; ++kk ) {
            mb  = min( A.mb, A.m - kk*A.mb );
            lda = BLKLDD( A, kk );
            Akj = A(kk,jj);
            for( k = i2; k < mb; ++k ) {
                Akj[k + j*lda] *= scale;
            }
            i2 = 0;
        }
    }
    /* unlike LAPACK, we explicitly store the 1 in the vector, */
    /* and return beta separately. */
    lda = BLKLDD( A, ii );
    Aij = A(ii,jj);
    Aij[i + j*lda] = 1.;
}

/**
 *
 * @file core_dgeqp3_update.c
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
#include <cblas.h>
#ifdef PLASMA_WITH_MKL
#include <mkl_lapacke.h>
#else
#include <lapacke.h>
#endif
#include "common.h"

#define A(m,n) BLKADDR( A, double, m, n )

/***************************************************************************//**
 *
 * @ingroup CORE_double
 *
 * CORE_dgeqp3_update updates row k of one tile of A
 * and subtracts that row from the column norms.
 *
 *******************************************************************************
 *
 * @param[in] Ajj
 *         Diagonal tile (jj,jj) of A.
 *
 * @param[in] lda1
 *         Leading dimension of Ajj.
 *
 * @param[in,out] Ajk
 *         Tile (jj,kk) of A, kk >= jj.
 *         On exit, updates row joff+k (i.e., as if Q was applied to trailing matrix).
 *
 * @param[in] lda2
 *         Leading dimension of Ajk.
 *
 * @param[in] Fk
 *         Tile kk of F.
 *
 * @param[in] ldf
 *         Leading dimension of Fk.
 *
 * @param[in] joff
 *         Row offset.
 *
 * @param[in] k
 *         Update row joff+k, based on having factored k columns.
 *         (That is, joff columns of this tile were factored in previous panels;
 *          k columns have been factored during this panel.)
 *
 * @param[in] koff
 *         Column to start updating.
 *         For diagonal tile, koff=joff+k+1, else koff=0.
 *
 * @param[in] nb
 *         Number of columns in kk-th block-column of A.
 *
 * @param[in,out] norms1
 *         kk-th block of partial column norms vector, dimension nb.
 *         On exit, norms1[koff:nb] -= Ajk[k, koff:nb ].
 *
 * @param[in,out] norms2
 *         kk-th block of original column norms vector, dimension nb.
 *         Unchanged on exit, except if cancellation is detected for
 *         some column j, sets norm2[j] = -1 and sets info = 1.
 *
 * @param[out] info
 *         Set to true if numerical instability (cancellation) is detected
 *         in updating column norms. dgeqp3 handles this error.
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_dgeqp3_update = PCORE_dgeqp3_update
#define CORE_dgeqp3_update PCORE_dgeqp3_update
#endif
void CORE_dgeqp3_update( const double *Ajj, int lda1,
                         double       *Ajk, int lda2,
                         const double *Fk,  int ldf,
                         int joff, int k, int koff, int nb,
                         double *norms1, double *norms2,
                         int *info )
{
    double temp, temp2;
    double tol3z = sqrt( LAPACKE_dlamch_work('e'));
    const double zone  =  1.0;
    const double mzone = -1.0;
    int j;

    /* update row k of A -- this is vector*matrix */
    /* Ajk[k,j:nb] -= Ajj[k,0:k+1] * Fk[j:nb,0:k+1].T */
    cblas_dgemm( CblasColMajor, CblasNoTrans, CblasTrans, 1, nb-koff, k+1,
                 (mzone), &Ajj[joff+k + joff*lda1], lda1,
                                     &Fk [koff              ], ldf,
                 (zone),  &Ajk[joff+k + koff*lda2], lda2 );

    for( j = koff; j < nb; ++j ) {
        if ( norms1[j] != 0. ) {
            /* NOTE: The following lines follow from the analysis in Lapack Working Note 176. */
            temp = fabs( Ajk[joff+k + j*lda2] ) / norms1[j];
            temp = max( 0., (1. + temp)*(1. - temp) );
            temp2 = norms1[j] / norms2[j];
            temp2 = temp * temp2*temp2;
            norms1[j] = norms1[j]*sqrt( temp );
            if( temp2 <= tol3z ) {
                /* flag numerical problem (i.e., cancellation) in updating norm.
                *  norms1[j] will be re-computed. Above we stored the inaccurate
                *  value anyway to allow comparison with the accurate value, for
                *  easier debugging. */
                norms2[j] = -1;
                *info = 1;
            }
        }
    }
}

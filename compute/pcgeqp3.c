/**
 *
 * @file pcgeqp3.c
 *
 *  PLASMA auxiliary routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Mark Gates
 * @date 2010-11-15
 * @generated c Fri Apr  1 11:03:00 2016
 *
 **/
#include "common.h"

#define A(m,n) BLKADDR(A, PLASMA_Complex32_t, m, n)

static const PLASMA_Complex32_t zone  = (PLASMA_Complex32_t)  1.;
static const PLASMA_Complex32_t mzone = (PLASMA_Complex32_t) -1.;
static const PLASMA_Complex32_t zzero = (PLASMA_Complex32_t)  0.;

#define BLKM(A,k) ((k) == (A).mt-1 ? (A).m - (k)*(A).mb : (A).mb)
#define BLKN(A,k) ((k) == (A).nt-1 ? (A).n - (k)*(A).nb : (A).nb)

/***************************************************************************//**
 *  Parallel tile QR factorization with column pivoting - dynamic scheduling
 **/
void plasma_pcgeqp3_quark( PLASMA_desc A, int *jpvt, PLASMA_Complex32_t *tau,
                           PLASMA_Complex32_t *work, float *rwork,
                           PLASMA_sequence *sequence, PLASMA_request *request )
{
    /* Aij is A(ii,jj) tile; similarly for Ajj, Aik, Ajk.
     * Fk  is F(kk) tile
     * ii, jj, kk are tile indices
     * j is global column index (0 to A.n)
     * k is column index within current panel
     * ioff is offset from beginning of ii-th block row
     * joff, koff are column offsets from beginning of jj-th and kk-th block columns
     * im is # rows in ii-th block row
     * kn is # cols in kk-th block row
     * jb is # cols in current panel -- which is NOT the number of rows
     *       nor cols in A(jj,jj) tile -- see formula
     */
    PLASMA_Complex32_t *Aij, *Ajj, *Aik, *Ajk, *Fk;
    int ii, ioff, im;
    int jj, joff, jb, j;
    int kk, koff, kn, k;
    int nj, jk, ldf, lda, minmn;

    /* TODO: must be static to exist after return. Should be put into workspace.
     * I think aux is k-vector, leaving (n-k) space available in work (? verify). No space if n <= nb.
     * With immediate execution, this shouldn't matter. */
    static PLASMA_Complex32_t beta;
    static int info = 0;

    plasma_context_t *plasma;
    Quark_Task_Flags task_flags = Quark_Task_Flags_Initializer;

    PLASMA_Complex32_t *F = work;
    ldf = A.n;
    PLASMA_Complex32_t *aux = &work[A.nb*ldf];
    float *norms1 = rwork;
    float *norms2 = rwork + A.n;

    plasma = plasma_context_self();
    if (sequence->status != PLASMA_SUCCESS)
        return;
    QUARK_Task_Flag_Set(&task_flags, TASK_SEQUENCE, (intptr_t)sequence->quark_sequence);
    plasma_profile_by_function( &task_flags, GEQRT );

    QUARK_CORE_cgeqp3_init( plasma->quark, &task_flags, A.n, jpvt );

    /* Compute the initial norm of each column of A */
    for( kk = 0; kk < A.nt; kk++ ) {
        kn = BLKN( A, kk );
        /* set norms2[:] = -1. to force cgeqp3_norms to compute all column norms */
        QUARK_CORE_slaset( plasma->quark, &task_flags,
            PlasmaUpperLower, kn, 1, mzone, mzone, &norms2[kk*A.nb], 1 );
        QUARK_CORE_cgeqp3_norms( plasma->quark, &task_flags,
            plasma_desc_submatrix( A, 0, kk*A.nb, A.m, kn ),
            0, 0,
            &norms1[kk*A.nb],
            &norms2[kk*A.nb] );
    }

    minmn = min( A.m, A.n );

    j = 0;
    while( j < minmn ) {
        jj   = (int)( j / A.nb );
        Ajj  = A(jj,jj);
        nj   = A.n - j;
        joff = j % A.nb;
        jb   = min( A.nb - joff, minmn - j );

        /* F = 0, so gemv can accumulate into F */
        for( kk = jj; kk < A.nt; ++kk ) {
            kn = BLKN( A, kk );
            Fk = &F[(kk-jj)*A.nb];
            QUARK_CORE_claset( plasma->quark, &task_flags,
                PlasmaUpperLower, kn, jb,
                zzero, zzero, Fk, ldf );
        }

        k = 0;
        while( k < jb ) {
            jk = j + k;

            /* this executes synchronously, allowing us to check info immediately
             * and break the loop if a bad column norm (due to cancellation) is detected. */
            QUARK_CORE_cgeqp3_pivot( plasma->quark, &task_flags,
                A, F, ldf, jj, joff+k, jpvt, norms1, norms2, &info );
            if ( info != 0 ) {
                break;
            }

            /* update current column */
            if ( k > 0 ) {
                /* A[jk:m,jk] -= A[jk:m,j:jk] * F[k,0:k]^H */
                ioff = joff + k;
                for( ii = jj; ii < A.mt; ++ii ) {
                    im = BLKM( A, ii );
                    lda = BLKLDD( A, ii );
                    Aij = A(ii,jj);
                    /* really gemv, but F needs to be conjfugated, so use gemm */
                    QUARK_CORE_cgemm_tile( plasma->quark, &task_flags,
                        PlasmaNoTrans, PlasmaConjTrans, im-ioff, 1, k, A.nb,
                        &mzone, &Aij[ioff +     joff*lda], lda,
                                &  F[joff+k             ], ldf,
                        &zone,  &Aij[ioff + (joff+k)*lda], lda,
                        Aij, F, Aij );  /* tile dependencies */
                    ioff = 0;
                }
            }

            /* Householder */
            QUARK_CORE_cgeqp3_larfg( plasma->quark, &task_flags,
                A, jj, jj, joff+k, joff+k, &tau[jk], &beta );

            /* compute F */
            if ( k < nj-1 ) {
                /* F[k+1:nj,k] = tau[jk] * A[jk:m,jk+1:n]^T * A[jk:m,jk] */
                /* assumes F = 0, above */
                ioff = joff + k;
                for( ii = jj; ii < A.mt; ++ii ) {
                    im = BLKM( A, ii );
                    lda = BLKLDD( A, ii );
                    Aij = A(ii,jj);
                    koff = joff + k+1;
                    for( kk = jj; kk < A.nt; ++kk ) {
                        kn = BLKN( A, kk );
                        Aik = A(ii,kk);
                        Fk  = &F[(kk-jj)*A.nb];
                        QUARK_CORE_cgemv_tile( plasma->quark, &task_flags,
                            PlasmaConjTrans, im-ioff, kn-koff,
                            &tau[jk], &Aik[ioff +     koff*lda], lda,
                                      &Aij[ioff + (joff+k)*lda], 1,
                            &zone,    & Fk[koff +        k*ldf], 1,
                            Aik, Aij, Fk );  /* tile dependencies */
                        koff = 0;
                    }
                    ioff = 0;
                }
            }

            /* incremental update of F */
            if ( k > 0 ) {
                /* aux = - A[jk:m,j:jk]^T * A[jk:m,jk] */
                ioff = joff + k;
                for( ii = jj; ii < A.mt; ++ii ) {
                    im = BLKM( A, ii );
                    lda = BLKLDD( A, ii );
                    Aij = A(ii,jj);
                    if ( ii == jj ) {
                        /* first gemv, beta=0; later beta=1 */
                        QUARK_CORE_cgemv_tile( plasma->quark, &task_flags,
                            PlasmaConjTrans, im-ioff, k,
                            &mzone, &Aij[ioff +     joff*lda], lda,
                                    &Aij[ioff + (joff+k)*lda], 1,
                            &zzero, aux, 1,
                            Aij, Aij, aux );  /* tile dependencies */
                    }
                    else {
                        QUARK_CORE_cgemv_tile( plasma->quark, &task_flags,
                            PlasmaConjTrans, im-ioff, k,
                            &mzone, &Aij[ioff +     joff*lda], lda,
                                    &Aij[ioff + (joff+k)*lda], 1,
                            &zone,  aux, 1,
                            Aij, Aij, aux );  /* tile dependencies */
                    }
                    ioff = 0;
                }
                /* F[0:nj,k] += tau[jk] * F[0:nj,0:k] * aux */
                /* (compared to lapack, minus moved above into aux) */
                koff = joff;
                for( kk = jj; kk < A.nt; ++kk ) {
                    kn = BLKN( A, kk );
                    Fk = &F[(kk-jj)*A.nb];
                    QUARK_CORE_cgemv_tile( plasma->quark, &task_flags,
                        PlasmaNoTrans, kn-koff, k,
                        &tau[jk], &Fk[koff], ldf,
                                  aux, 1,
                        &zone,    &Fk[koff + k*ldf], 1,
                        Fk, aux, Fk );  /* tile dependencies */
                    koff = 0;
                }
            }

            /* update pivot row and norms */
            /* A[jk,jk+1:n] -= A[jk,j:jk+1] * F[k+1:nj,0:k+1]^T */
            /* norms1[jk+1:n] =  sqrt( norms1[jk+1:n]**2 - A[jk,jk+1:n]**2 ) */
            koff = joff + k+1;
            lda = BLKLDD( A, jj );
            for( kk = jj; kk < A.nt; ++kk ) {
                kn = BLKN( A, kk );
                Ajk = A(jj,kk);
                Fk  = &F[(kk-jj)*A.nb];
                QUARK_CORE_cgeqp3_update( plasma->quark, &task_flags,
                    Ajj, lda, Ajk, lda, Fk, ldf,
                    joff, k, koff, kn,
                    &norms1[kk*A.nb], &norms2[kk*A.nb], &info );
                koff = 0;
            }

            /* save beta (from clarfg) to diagonal */
            QUARK_CORE_csetvar( plasma->quark, &task_flags,
                &beta, &Ajj[joff+k + (joff+k)*lda],
                Ajj );  /* tile dependencies */

            k += 1;
        }

        /* trailing matrix update */
        ioff = joff + k;
        for( ii = jj; ii < A.mt; ++ii ) {
            im = BLKM( A, ii );
            lda = BLKLDD( A, ii );
            koff = joff + k;
            for( kk = jj; kk < A.nt; ++kk ) {
                /* update partial tiles (ii=jj or  kk=jj) if k-loop aborted before jnb (k<jnb), and */
                /* update whole   tiles (ii>jj and kk>jj) */
                if ( k < jb || (ii > jj && kk > jj)) {
                    kn = BLKN( A, kk );
                    Aik = A(ii,kk);          /* Aik  im x kn */
                    Aij = A(ii,jj);          /* Aik  im x jb */
                    Fk  = &F[(kk-jj)*A.nb];  /* Fk^T jb x kn */
                    QUARK_CORE_cgemm_tile( plasma->quark, &task_flags,
                        PlasmaNoTrans, PlasmaConjTrans, im-ioff, kn-koff, jb, A.nb,
                        &mzone, &Aij[ioff + joff*lda], lda,
                                & Fk[koff           ], ldf,
                        &zone,  &Aik[ioff + koff*lda], lda,
                        Aij, Fk, Aik );  /* tile dependencies */
                }
                koff = 0;
            }
            ioff = 0;
        }

        /* re-compute bad column norms */
        if ( info != 0 ) {
            info = 0;
            koff = joff + k;
            for( kk = jj; kk < A.nt; ++kk ) {
                kn = BLKN( A, kk );
                QUARK_CORE_cgeqp3_norms( plasma->quark, &task_flags,
                    plasma_desc_submatrix( A, jj*A.mb, kk*A.nb, A.m - jj*A.mb, kn ),
                    joff+k, koff,
                    &norms1[kk*A.nb],
                    &norms2[kk*A.nb] );
                koff = 0;
            }
        }

        j += k;
    }
}

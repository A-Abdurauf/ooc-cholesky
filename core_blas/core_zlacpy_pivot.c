/**
 *
 * @file core_zlacpy_pivot.c
 *
 *  PLASMA core_blas kernel
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Mathieu Faverge
 * @date 2013-02-01
 * @precisions normal z -> c d s
 *
 **/
#ifdef PLASMA_WITH_MKL
#include <mkl_lapacke.h>
#else
#include <lapacke.h>
#endif
#include "common.h"

#define A(m, n) BLKADDR(descA, PLASMA_Complex64_t, m, n)

/***************************************************************************//**
 *
 * @ingroup CORE_PLASMA_Complex64_t
 *
 *  CORE_zlacpy_pivot extracts the original version of the rows selected by the
 *  ipiv array and copies them into a new buffer.
 *
 *  This kernel is used by tournament pivoting algorithms, to extract the
 *  selected rows from the original matrix that will make it to the next level
 *  of the tournament.
 *
 *******************************************************************************
 *
 *  @param[in] descA
 *          The descriptor of the matrix A in which the kernel will extract the
 *          original rows.
 *
 *  @param[in] direct
 *          @arg PlasmaRowwise: The extracted rows are stored in column major layout.
 *          @arg PlasmaColumnwise: The extracted rows are store in row major layout.
 *
 *  @param[in] k1
 *          The first element of IPIV for which a row interchange will
 *          be done.
 *
 *  @param[in] k2
 *          The last element of IPIV for which a row interchange will
 *          be done.
 *
 *  @param[in] ipiv
 *          The pivot indices; Only the element in position k1 to k2
 *          are accessed. The pivots should be included in the interval 1 to A.m
 *
 *  @param[in,out] rankin
 *          On entry, the global indices relative to the full matrix A
 *          factorized, in the local sub-matrix. If init == 1, rankin is
 *          initialized to A.i, .. A.i+descA.m
 *          On exit, rows are permutted according to ipiv.
 *
 *  @param[out] rankout
 *          On exit, contains the global indices of the first (k2-k1+1) rows.
 *
 *  @param[out] A
 *          An lda-by-descA.n matrix. On exit, A contains the original version
 *          of the rows selected by the pivoting process.
 *
 *  @param[in] lda
 *          The leading dimension of the array A.  lda >= max(1,(k2-k1+1)).
 *
 *  @param[in] init
 *          True if rankin needs to be initialized.
 *          False, if rankin contains already initialized data.
 *
 *******************************************************************************
 *
 * @return
 *         \retval PLASMA_SUCCESS successful exit
 *         \retval <0 if INFO = -k, the k-th argument had an illegal value
 *
 ******************************************************************************/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_zlacpy_pivot = PCORE_zlacpy_pivot
#define CORE_zlacpy_pivot PCORE_zlacpy_pivot
#endif
int CORE_zlacpy_pivot( const PLASMA_desc descA,
                       PLASMA_enum direct, int k1, int k2, const int *ipiv,
                       int *rankin, int *rankout,
                       PLASMA_Complex64_t *A, int lda,
                       int init)
{
    int i, ip, it, ir, ld;
    const int *lpiv;
    int *ro = rankout;

    /* Init rankin if first step */
    if ( init ) {
        int val = descA.i;
        for(i=0; i<descA.m; i++, val++) {
            rankin[i] = val;
        }
    }

    /* Generate the rankout */
    ro = rankout;
    lpiv = ipiv;
    for(i=k1-1; i<k2; i++, ro++, lpiv++) {
        *ro = rankin[ (*lpiv) - 1 ];
        rankin[ (*lpiv) - 1 ] = rankin[ i ];
    }

    /* Extract the rows */
    if (direct == PlasmaRowwise) {
       ro = rankout;
       for(i=k1-1; i<k2; i++, ro++) {
           ip = (*ro) - descA.i;
           it = ip / descA.mb;
           ir = ip % descA.mb;
           ld = BLKLDD(descA, it);
           cblas_zcopy(descA.n, A(it, 0) + ir, ld,
                                A + i,         lda );
       }
    }
    else {
       ro = rankout;
       for(i=k1-1; i<k2; i++, ro++) {
           ip = (*ro) - descA.i;
           it = ip / descA.mb;
           ir = ip % descA.mb;
           ld = BLKLDD(descA, it);
           cblas_zcopy(descA.n, A(it, 0) + ir, ld,
                                A + i*lda,     1 );
       }
    }
    return PLASMA_SUCCESS;
}

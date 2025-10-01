/**
 * @file core_dlaed2_compressq.c
 *
 *  PLASMA computational routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Gregoire Pichon
 * @author Azzam Haidar
 * @date 2014-07
 * @precisions normal d -> s
 *
 **/
#include "common.h"
#include <math.h>
#ifdef PLASMA_WITH_MKL
#include <mkl_lapacke.h>
#else
#include <lapacke.h>
#endif

/**
 ******************************************************************************
 *
 * @ingroup CORE_double
 *
 *  CORE_dlaed2_compressq - Copy a set of vectors from the Q matrix into the
 *  compressed Q2 matrix.
 *  Q2 matrix is oragnized as follow:
 *
 *  Q2 = (  q1_1  mix1  0   def1 )
 *       (   0    mix2  q2  def2 )
 *  and stored in memory by columns in the order (q1_1, mix1, mix2, q2, def),
 *  with q1 and q2 the non deflated values, of the two subproblems, (mix1, mix2)
 *  the mixed eigenvalues and def the deflated eiegnvalues.
 *
 * Rk: Note that here we follow the same order that LAPACK dlaed2 routine, so
 * the copy from Q to Q2 is done such that Q2 is filled in a continuous manner.
 *
 *******************************************************************************
 *
 * @param[in] n
 *          n specifies the dimension of the symmetric tridiagonal matrix
 *
 * @param[in] n1
 *          n1 specifies the location of the last eigenvalue of the first subproblem
 *          min(1, n) <= n1 <= n/2
 *
 * @param[in] start
 *          start specifies the first sorted column index to be considered by
 *          this kernel.
 *
 * @param[in] end
 *          end specifies the last column index to be considered by this kernel.
 *          start <= end <= n.
 *
 * @param[in] INDX
 *          The permutation array used to sort Q into Q2. Array of dimension n,
 *          but only INDX[i] for start <= i < end are referenced.
 *
 * @param[in] ctot
 *          ctot[i] is the number of columns of type i, as defined with INDX.
 *          0- number of colums with non-zero in the upper half only
 *          1- number of dense colums
 *          2- number of colums with non-zero in the lower half only
 *          3- number of deflated columns
 *
 * @param[in] Q
 *          On entry, Q contains the eigenvectors in the uncompressed form.
 *          WARNING: Q is the pointer to the full matrix even if only a subset
 *          of vectors will be extract.
 *
 * @param[in] LDQ
 *          LDQ specifies the leading dimension of Q
 *
 * @param[out] Q2
 *          On exit, columns start to end from Q2 stores the eigenvectors in the
 *          compressed form and sorted by type.
 *          WARNING: Q2 is the pointer to the full matrix even if only columns
 *          start to end are touched.
 *
 *******************************************************************************/

#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_dlaed2_compressq = PCORE_dlaed2_compressq
#define CORE_dlaed2_compressq PCORE_dlaed2_compressq
#endif
void CORE_dlaed2_compressq(int n, int n1, const int *INDX, const int *ctot,
                           const double *Q, int LDQ, double *Q2,
                           int start, int end)
{
    double *Q2_q1, *Q2_q2, *Q2_df;
    int pend = ctot[0]; /* Partial end of for-loops */
    int i, is;
    int n2 = n-n1;

    /* end cannot be larger than n */
    if (end > n)
        end = n;

    Q2_q1 = Q2;
    Q2_q2 = Q2_q1 + (ctot[0] + ctot[1])*n1;
    Q2_df = Q2_q2 + (ctot[1] + ctot[2])*n2;

    /* Vectors of q1 */
    for (i=start; (i<end) && (i<pend); i++) {
        is = INDX[i] * LDQ;
        cblas_dcopy(n1, Q+is, 1, Q2_q1 + n1*i, 1);
    }

    /* Mixed Vectors */
    pend += ctot[1];
    for (; (i<end) && (i<pend); i++) {
        is = INDX[i] * LDQ;
        cblas_dcopy(n1, Q+is,    1, Q2_q1+n1*i,           1);
        cblas_dcopy(n2, Q+is+n1, 1, Q2_q2+n2*(i-ctot[0]), 1);
    }

    /* Vectors of q2 */
    pend += ctot[2];
    for (; (i<end) && (i<pend); i++) {
        is = INDX[i] * LDQ + n1;
        cblas_dcopy(n2, Q+is,    1, Q2_q2+n2*(i-ctot[0]), 1);
    }

    /* Deflated Vectors */
    for (; i<end; i++) {
        is = INDX[i] * LDQ;
        cblas_dcopy(n, Q+is,     1, Q2_df+n *(i-pend),    1);
    }
}


/**
 *****************************************************************************
 *
 * @ingroup CORE_double
 *
 *  CORE_dlaed2_copydef - Copy back a portion of the deflated eigenvectors from
 *  Q2 to Q.
 *
 *******************************************************************************
 *
 * @param[in] n
 *          n specifies the dimension of the symmetric tridiagonal matrix
 *
 * @param[in] n1
 *          n1 specifies the location of the last eigenvalue of the first subproblem
 *          min(1, n) <= n1 <= n/2
 *
 * @param[in] K
 *          K specifies the number of non-deflated eigenvalues
 *
 * @param[in] start
 *          start specifies the first column index to be considered by this kernel
 *          note that this index can be < K, it is supported
 *
 * @param[in] end
 *          end specifies the last column index to be considered by this kernel
 *
 * @param[in] ctot
 *          ctot[i] is the number of columns of type i, as defined with INDX.
 *          0- number of colums with non-zero in the upper half only
 *          1- number of dense colums
 *          2- number of colums with non-zero in the lower half only
 *          3- number of deflated columns
 *
 * @param[out] Q
 *          On exit, Q(start:end-1) contains the updated eigenvectors
 *
 * @param[in] LDQ
 *          LDQ specifies the leading direction of Q
 *
 * @param[in] Q2
 *          On entry, Q2 contains the saved eigenvectors in the compressed form.
 *
 *
 *******************************************************************************/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_dlaed2_copydef = PCORE_dlaed2_copydef
#define CORE_dlaed2_copydef PCORE_dlaed2_copydef
#endif
void CORE_dlaed2_copydef(int n, int n1, int K, const int *ctot,
                         double *Q, int LDQ, const double *Q2,
                         int start, int end )
{
    int size;

    start = max( start, K );
    size  = max( 0, end-start );

    if ( size > 0 ) {
        size_t offset = 0;

        offset  = (size_t)(   n1 ) * (size_t)(ctot[0] + ctot[1]);
        offset += (size_t)((n-n1)) * (size_t)(ctot[1] + ctot[2]);
        offset += (size_t)(   n  ) * (size_t)(start - K);

        LAPACKE_dlacpy_work(LAPACK_COL_MAJOR,
                            lapack_const(PlasmaUpperLower), n, size,
                            Q2+offset, n, Q+LDQ*start, LDQ);
    }
}

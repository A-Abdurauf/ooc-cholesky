/**
 * @file core_dlaed3_merge.c
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
#include <math.h>
#ifdef PLASMA_WITH_MKL
#include <mkl_lapacke.h>
#else
#include <lapacke.h>
#endif
#include "common.h"

/**
 *******************************************************************************
 *
 * @ingroup CORE_double
 *
 *  CORE_dlaed3_merge - Merge the eigenvalues of two subproblems and generates
 *  the index arrays to sort the associated eigenvectors.
 *
 *******************************************************************************
 *
 * @param[in] n
 *          n specifies the dimension of the symmetric tridiagonal matrix
 *
 * @param[in] K
 *          K specifies the number of non-deflated eigenvalues
 *
 * @param[out] D
 *          Array of double of size n.
 *          On exit, stores the sorted eigenvalues.
 *
 * @param[out] INDXQ
 *          Array of integers of size n.
 *          On exit, INDXQ stores the permutation array to sort the
 *          eigenvectors.
 *
 ***************************************************************************/
void
CORE_dlaed3_merge( int n, int K, double *D, int *INDXQ )
{
    int i, n2;

    /*
     * Merge the eigenvalues in D in asscending order and generates INDXQ
     */
    if (K==0){
        for (i=0; i<n; i++){
            INDXQ[i] = i;
        }
    }
    else{
        int ione  = 1;
        int imone = -1;
        n2 = n - K;
        PLASMA_FCALL(dlamrg, DLAMRG)(&K, &n2, D, &ione, &imone, INDXQ);
        for (i=0; i<n; i++){
            INDXQ[i]--; /* -1 to make it in C mode */
        }
    }
}

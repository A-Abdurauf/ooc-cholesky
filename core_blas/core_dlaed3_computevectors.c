/**
 * @file core_dlaed3_computevectors.c
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
 *  CORE_dlaed3_computevectors - Compute the subset(start:end) of eigenvectors
 *  of the modified rank-1 system.
 *
 *******************************************************************************
 *
 * @param[in] K
 *          K specifies the number of non-deflated eigenvalues
 *
 * @param[in] il_nondef
 *          The first eigenvector index when computing only a subset of all eigenpairs
 *
 * @param[in] iu_nondef
 *          The last eigenvector index when computing only a subset of all eigenpairs
 *
 * @param[in,out] Q
 *          Q is an array of dimension (LDQ, end)
 *          Initially the first K columns of the full Q matrix are used as workspace.
 *          On output the Q(start:end) contain the updated eigenvectors.
 *
 * @param[in] LDQ
 *          LDQ specifies the leading dimension of Q
 *
 * @param[in] W
 *          W is the result of the previous computation after reduction. See
 *          core_dlaed3_computeW() and core_dlaed3_reduceW().
 *          The first K elements of this array contain the components of the
 *          deflation-adjusted updating vector. Destroyed on output.
 *
 * @param[in] S
 *          S is local workspace of size K.
 *
 * @param[in] INDXC
 *          INDXC is the permutation used to arrange the columns of the deflated
 *          Q matrix into three groups (see core_dlaed2_computeK()).
 *          The rows of the eigenvectors found by core_dlaed4() must be likewise
 *          permuted before the matrix multiply of core_dlaed3_updatevectors()
 *          can take place.
 *
 * @param[in] start
 *          start specifies the first column index to be considered by this kernel
 *
 * @param[in] end
 *          end specifies the last column index to be considered by this kernel
 *
 ***************************************************************************/

#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_dlaed3_computevectors = PCORE_dlaed3_computevectors
#define CORE_dlaed3_computevectors PCORE_dlaed3_computevectors
#endif
void CORE_dlaed3_computevectors(int K, int il_nondef, int iu_nondef,
                                double *Q, int LDQ, double *W, double *S,
                                const int *INDXC,
                                int start, int end)
{
    double *lQ = NULL;
    int i, j, ind;

    /* Adjust the subset */
    start = max( 0,   max( start, il_nondef ) );
    end   = min( end, min( K,     iu_nondef ) );

    if ( K == 1 ) {
        return;
    }

    lQ = Q + start * LDQ;

    if ( K == 2 ) {
        assert( INDXC[0] == 0 || INDXC[0] == 1 );
        assert( INDXC[1] == 0 || INDXC[1] == 1 );

        for (j=start; j<end; j++, lQ+=LDQ) {
            W[0] = lQ[0];
            W[1] = lQ[1];
            ind  = INDXC[0];
            lQ[0] = W[ind];
            ind  = INDXC[1];
            lQ[1] = W[ind];
        }
        return;
    }

    for (j=start; j<end; j++, lQ+=LDQ){
        double tmp;

        for (i=0; i<K; i++){
            S[i] = W[i] / lQ[i];
        }

        tmp = cblas_dnrm2(K, S, 1);
        for (i=0; i<K; i++){
            ind = INDXC[i];
            lQ[i] = S[ind] / tmp;
        }
    }

    return;
}

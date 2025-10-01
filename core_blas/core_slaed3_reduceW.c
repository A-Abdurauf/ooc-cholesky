/**
 * @file core_slaed3_reduceW.c
 *
 *  PLASMA computational routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Gregoire Pichon
 * @author Azzam Haidar
 * @date 2014-07
 * @generated s Fri Apr  1 11:02:33 2016
 *
 **/
#include "common.h"
#include <math.h>
#include <stdlib.h>

/**
 ******************************************************************************
 *
 * @ingroup CORE_float
 *
 *  CORE_slaed3_reduceW - Computes the reduction of multiple W
 *  l threads were computing Wred(:,j) and this kernel will compute the
 *  Pi( Wred(i,j), j=1..l )
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
 * @param[in] l
 *          l specifies the number of columns of Wred
 *
 * @param[out] Q
 *
 * @param[in] LDQ
 *          LDQ specifies the leading direction of Q
 *
 * @param[in] Wred
 *          Wred[:,j] corresponds to the local W for a previous task
 *
 * @param[out] W
 *          On exit, W(i) = sqrt( Pi( Wred(i,j), j=1..l ) * Q(i,i) )
 *
 *******************************************************************************/

#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_slaed3_reduceW = PCORE_slaed3_reduceW
#define CORE_slaed3_reduceW PCORE_slaed3_reduceW
#endif
void CORE_slaed3_reduceW(int n, int n1, int K, int l,
                         const float *Q, int LDQ,
                         const float *Wred, float *W)
{
    float *S;
    int i, j;

    S = malloc(n*sizeof(float));

    /* multiply the partial W that has been computed by slaed3_computeW
     * to generate the final W */

    if ( K > 2 ) {

        /* Startup reduction */
        cblas_scopy(K, Wred, 1, S, 1);

        for (i=1; i<l; i++){
            for (j=0; j<K; j++){
                S[j] = S[j] * Wred[n*i+j];
            }
        }

        /* Update W according to previous value */
        for (i=0; i<K; i++){
            S[i] = S[i] * Q[LDQ*i+i];

            if (W[i] > 0.0){
                W[i] =  sqrt(-S[i]);
            }
            else {
                W[i] = -sqrt(-S[i]);
            }
        }
    }

    free(S);
}

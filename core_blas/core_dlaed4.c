/**
 * @file core_dlaed4.c
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

/***************************************************************************//**
 *
 * @ingroup CORE_double
 *
 *  CORE_dlaed4 - solve the secular equation for indexes between start and end
 *
 *******************************************************************************
 *
 * @param[in] n
 *          n specifies the dimension of the symmetric tridiagonal matrix
 *
 * @param[in] K
 *          K specifies the number of non-deflated eigenvalues
 *
 * @param[in,out] D
 *          On entry, D contains the eigenvalues of the two submatrices to be merged.
 *          On exit, D contains the updated eigenvalues sorted into increasing order.
 *
 * @param[in] beta
 *          beta_bis[0] specifies the rank-1 approximation that was used for splitting
 *          the problem into two subproblems.
 *
 * @param[in,out] Q
 *          On exit, Q contains the updated eigenvectors
 *
 * @param[in] LDQ
 *          LDQ specifies the leading direction of Q
 *
 * @param[in] D0
 *          On entry, D0 conatins the original sorted eigenvalues.
 *
 * @param[in] Z
 *          Z contains the components of the updating vectors.
 *
 * @param[out] INDX
 *          The permutation used to sort the contents of DLAMBDA into ascending order
 *
 * @param[in] start
 *          start specifies the first column index to be considered by this kernel
 *          0 <= start <= end
 *
 * @param[in] end
 *          end specifies the last column index to be considered by this kernel
 *          start <= end <= n
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *
 *******************************************************************************/

#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_dlaed4 = PCORE_dlaed4
#define CORE_dlaed4 PCORE_dlaed4
#endif
int CORE_dlaed4(int n, int K,
                double *D, double beta,
                double *Q, int LDQ,
                const double *D0, const double *Z,
                const int *INDX,
                int start, int end )
{
    int i, is, id;
    int info;

    end = min(end, n);

    for (i=start; i<end; i++){
        is = INDX[i];
        if (is < K){
            id = is+1;
            PLASMA_FCALL(dlaed4, DLAED4)(&K, &id, D0, Z, Q+LDQ*is, &beta, D+is, &info);
            if (info != 0){
                coreblas_error(info, "numerical error in dlaed4\n");
                return info;
            }
        }
    }
    return PLASMA_SUCCESS;
}

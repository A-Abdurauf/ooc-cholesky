/**
 * @file core_dlaed0_betaapprox.c
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

/**
 *******************************************************************************
 *
 * @ingroup CORE_double
 *
 *  CORE_dlaed0_betaapprox - compute the rank-1 approximation for each
 *  subproblem.  When a larger problem is split into two subproblems, the last
 *  diagonal element of the first subproblem and the first diagonal element of
 *  the second subproblem are set to D[i]-E[i+1]
 *
 *******************************************************************************
 *
 * @param[in] subpbs
 *          specifies the number of merge to be done at the leaf of the tree
 *
 * @param[in] subpbs_info
 *          subpbs_info[i] is the position of the first eigenvalue of the ith
 *          eigenproblem.
 *
 * @param[in,out] D
 *          On entry, D contains the diagonal elements of the two submatrices to
 *          be merged.
 *          On exit, D contains the updated rank-1 approximated diagonal
 *          elements.
 *
 * @param[in] E
 *          E contains the extra-diagonal elements of the two submatrices to be
 *          merged.
 *
 ******************************************************************************/
void
CORE_dlaed0_betaapprox(int subpbs, const int *subpbs_info, double *D, const double *E)
{
    double tmp;
    int i;
    int index;

    for (i=0; i<subpbs; i++){
        index = subpbs_info[i]-1;
        tmp = fabs(E[index]);
        D[index]   -= tmp;
        D[index+1] -= tmp;
    }
}

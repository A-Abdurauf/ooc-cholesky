/**
 *
 * @file core_pivot_update.c
 *
 *  PLASMA core_blas kernel
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Mathieu Faverge
 * @author Ichitaro Yamazaki
 * @date 2010-11-15
 *
 **/
#include "common.h"

/***************************************************************************//**
 *
 *  CORE_pivot_update convert the given ipiv returned by the top level of the
 *  tournament pivoting into a ipiv array for the full matrix.
 *
 *******************************************************************************
 *
 *  @param[in] m
 *          The number of elements in indices array. Correspond to the number of
 *          row of the sub-matrix given to the LU decomposition.  M >= 0.
 *
 *  @param[in] n
 *          The number of element in ipiv array. It also corresponds to the
 *          minimum number of rows and columns of the submatrix given to the LU
 *          decomposition.  0 <= N <= M.
 *
 *  @param[in,out] ipiv
 *          On entry, the IPIV array contains pivot comprised between 1 and M.
 *          On exit, the IPIV array contains pivots for the full matrix
 *          comprised between 1 and the number of rows of the matrix. (Maximum
 *          value stored in indices array).
 *
 *  @param[in,out] indices
 *          On entry, contains the positions of the local rows of the sub-matrix
 *          in the full matrix.
 *          On exit, this array is swapped accordint to the ipiv array.
 *
 *  @param[int] offset
 *          The offset added to every pivot values in the output IPIV array.
 *
 *  @param[int] init
 *          If init == 1, indices array is filled with the values
 *          offset, ... , offset+m
 *
 ******************************************************************************/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_pivot_update = PCORE_pivot_update
#define CORE_pivot_update PCORE_pivot_update
#endif
void CORE_pivot_update(int m, int n, int *ipiv, int *indices,
                       int offset, int init)
{
    int i, piv, ind;

    /* Initialize indices if it is the first level of the tournament */
    if ( init ) {
        for(i=0; i<m; i++) {
            indices[i] = offset+i;
        }
    }

    for(i=0; i<n; i++) {
        piv = ipiv[i]-1;    /* 0 < ... < m  */
        ind = indices[piv]; /* 0 < ... < Am */

        /* Apply the pivot on permutation */
        indices[ piv ] = indices[ i ];
        indices[ i ]   = ind;

        /* While the rows is one that is already swapped,
         * we search its new position */
        while( (ind-offset) < i ) {
            ind = indices[ ind-offset ];
        }
        ipiv[i] = ind+1;
    }
}

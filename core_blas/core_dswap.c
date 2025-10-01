/**
 * @file core_dswap.c
 *
 *  PLASMA computational routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Grégoire Pichon
 * @author Azzam Haidar
 * @date 2014-07
 * @generated d Fri Apr  1 11:02:32 2016
 *
 **/
#include "common.h"
#include <math.h>

/***************************************************************************//**
 *
 * @ingroup CORE_double
 *
 *  CORE_dswap - Extract the eigenvectors in the range [start,end-1]
 *  from work and copy them in order into the Q matrix thanks to the
 *  given permutation.
 *
 *******************************************************************************
 *
 * @param[in] m
 *          m specifies the number of entries in each eigenvector,
 *          i.e. the number of rows of the matrices Q and Work.
 *
 * @param[in] n
 *          n specifies the number of eigenvectors to copy from work
 *          to Q, i.e. the number of columns of the matrices W and
 *          work.
 *
 * @param[out] Q
 *          On entry, matrix of size LDQ -by- n.
 *          On exit, Q will contain the n sorted eigenvectors.
 *
 * @param[in] ldq
 *          ldq specifies the leading dimension of Q. ldq >= max(1,m)
 *
 * @param[in] work
 *          On entry work contains the non-sorted eigenvectors and is
 *          of dimension m-by-n.
 *
 * @param[in] perm
 *          The permutation array used to copy work into Q. On entry,
 *          the i-th eigenvector is stored in the column perm[i] of
 *          work, and is copied to the i-th column of Q.
 *
 * @param[in] start
 *          start specifies the first column index to be considered by
 *          this kernel.
 *
 * @param[in] end
 *          end specifies the last column index to be considered by
 *          this kernel
 *
 ***************************************************************************/

#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_dswap = PCORE_dswap
#define CORE_dswap PCORE_dswap
#endif
void CORE_dswap(int m, int n,
                double *Q, int ldq,
                const double *work,
                const int *perm,
                int start, int end)
{
    int i;
    Q += start*ldq;
    for (i=start; i<end; i++, Q+=ldq){
        cblas_dcopy(m, work+m*perm[i], 1, Q, 1);
    }
}

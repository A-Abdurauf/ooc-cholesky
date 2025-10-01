/**
 *
 * @file core_cpltmg_fiedler.c
 *
 *  PLASMA core_blas kernel
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @generated c Fri Apr  1 11:02:29 2016
 *
 **/
#include <math.h>
#include "common.h"

/***************************************************************************//**
 *
 * @ingroup CORE_PLASMA_Complex32_t
 *
 *  CORE_cpltmg_fiedler is a kernel used in fiedler matrix generation
 *
 *  See http://www.mathworks.fr/fr/help/matlab/ref/gallery.html#f84-999960
 *
 *  Fiedler matrix of size n-by-n is defined throug a random vector c
 *  of size n, such that each element is equal to abs(n(i)-n(j)).
 *
 *  Matrix A has a dominant positive eigenvalue and all the other
 *  eigenvalues are negative.
 *
 *  Explicit formulas for inv(A) and det(A) are given in
 *  [Todd, J., Basic Numerical Mathematics, Vol. 2: Numerical Algebra,
 *  Birkhauser, Basel, and Academic Press, New York, 1977, p. 159] and
 *  attributed to Fiedler. These indicate that inv(A) is tridiagonal
 *  except for nonzero (1,n) and (n,1) elements.
 *
 *******************************************************************************
 *
 * @param[in] M
 *         The number of rows of the tile A to initialize. M >= 0.
 *
 * @param[in] N
 *         The number of columns of the tile A to initialize. N >= 0.
 *
 * @param[in] X
 *          X is a vector of dimension at least: ( 1 + ( M - 1 )*abs( incX ) )
 *          On entry, the vector used to initialize A.
 *
 * @param[in] incX
 *         On entry, incX specifies the increment for the elements of X.
 *         incX != 0.
 *
 * @param[in] Y
 *          Y is a vector of dimension at least: ( 1 + ( N - 1 )*abs( incY ) )
 *          On entry, the vector used to initialize A.
 *
 * @param[in] incY
 *         On entry, incY specifies the increment for the elements of Y.
 *         incY != 0.
 *
 * @param[out] A
 *         On entry, the M-by-N tile to be initialized.
 *         On exit, each element of A is defined by:
 *               A(i,j) = abs( X(i) - Y(j) )
 *
 * @param[in] LDA
 *         The leading dimension of the tile A. LDA >= max(1,M).
 *
 ******************************************************************************/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_cpltmg_fiedler = PCORE_cpltmg_fiedler
#define CORE_cpltmg_fiedler PCORE_cpltmg_fiedler
#endif
void CORE_cpltmg_fiedler( int M, int N,
                          const PLASMA_Complex32_t *X, int incX,
                          const PLASMA_Complex32_t *Y, int incY,
                                PLASMA_Complex32_t *A, int LDA )
{
    const PLASMA_Complex32_t *tmpX;
    int i, j;

    for (j=0; j<N; j++, Y+=incY) {
        tmpX = X;
        for (i=0; i<M; i++, tmpX+=incX, A++) {
            *A = cabsf( *tmpX - *Y );
        }
        A += (LDA - M);
    }
}

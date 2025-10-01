/**
 * @file core_clascal.c
 *
 *  PLASMA computational routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Mathieu Faverge
 * @date 2015-11-05
 * @generated c Fri Apr  1 11:02:32 2016
 *
 **/
#include <cblas.h>
#include "common.h"

/***************************************************************************//**
 *
 * @ingroup CORE_PLASMA_Complex32_t
 *
 *  CORE_clascal scales a two-dimensional matrix A. As opposite to
 *  CORE_clascl(), no checks is performed to prevent under/overflow. This should
 *  have been done at higher level.
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          Specifies the shape of A:
 *          = PlasmaUpperLower: A is a general matrix.
 *          = PlasmaUpper: A is an upper trapezoidal matrix.
 *          = PlasmaLower: A is a lower trapezoidal matrix.
 *
 * @param[in] m is the number of rows of the matrix A. m >= 0
 *
 * @param[in] n is the number of columns of the matrix A. n >= 0
 *
 * @param[in] alpha
 *            The scalar factor.
 *
 * @param[in,out] A is the matrix to be multiplied by alpha
 *
 * @param[in] lda is the leading dimension of the array A. lda >= max(1,m).
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *
 ******************************************************************************/
int
CORE_clascal( PLASMA_enum uplo, int m, int n,
              PLASMA_Complex32_t alpha, PLASMA_Complex32_t *A, int lda )
{
    int i;

    if ( (uplo != PlasmaUpperLower) &&
         (uplo != PlasmaUpper)      &&
         (uplo != PlasmaLower))
    {
        coreblas_error(1, "illegal value of uplo");
        return -1;
    }

    if (m < 0) {
        coreblas_error(2, "Illegal value of m");
        return -2;
    }
    if (n < 0) {
        coreblas_error(3, "Illegal value of n");
        return -3;
    }
    if ( (lda < max(1,m)) && (m > 0) ) {
        coreblas_error(6, "Illegal value of lda");
        return -6;
    }

    switch ( uplo ) {
    case PlasmaUpper:
        for(i=0; i<n; i++) {
            cblas_cscal( min( i+1, m ), CBLAS_SADDR(alpha), A+i*lda, 1 );
        }
        break;

    case PlasmaLower:
        for(i=0; i<n; i++) {
            cblas_cscal( max( m, m-i ), CBLAS_SADDR(alpha), A+i*lda, 1 );
        }
        break;
    default:
        if (m == lda) {
            cblas_cscal( m*n, CBLAS_SADDR(alpha), A, 1 );
        }
        else {
            for(i=0; i<n; i++) {
                cblas_cscal( m, CBLAS_SADDR(alpha), A+i*lda, 1 );
            }
        }
    }

    return PLASMA_SUCCESS;
}

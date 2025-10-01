/**
 *
 * @file core_dlatro.c
 *
 *  PLASMA core_blas kernel
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Azzam Haidar
 * @date 2010-11-15
 * @generated d Fri Apr  1 11:02:31 2016
 *
 **/
#ifdef PLASMA_WITH_MKL
#include <mkl_lapacke.h>
#else
#include <lapacke.h>
#endif
#include "common.h"

/***************************************************************************//**
 *
 * @ingroup CORE_double
 *
 *  CORE_dlatro transposes a m-by-n matrix out of place.
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          Specifies whether the matrix A is upper triangular or lower
 *          triangular:
 *          = PlasmaUpper: the upper triangle of A and the lower triangle of B
 *          are referenced.
 *          = PlasmaLower: the lower triangle of A and the upper triangle of B
 *          are referenced.
 *          = PlasmaUpperLower: All A and B are referenced.
 *
 * @param[in] trans
 *          Specifies whether the matrix A is transposed, not transposed or
 *          ugate transposed:
 *          = PlasmaNoTrans:   B is a copy of A (equivalent to dlacpy);
 *          = PlasmaTrans:     B is the transpose of A;
 *          = PlasmaTrans: B is the ugate transpose of A.
 *
 * @param[in] M
 *         Number of rows of the matrix A and number of columns of the matrix B, if trans == Pasma[Conj]Trans.
 *         Number of rows of the matrix A and the matrix B, if trans == PasmaNoTrans.
 *
 * @param[in] N
 *         Number of columns of the matrix A and number of rows of the matrix B, if trans == Pasma[Conj]Trans.
 *         Number of columns of the matrix A and of the matrix B, if trans == PlasmaNoTrans.
 *
 * @param[in] A
 *         Matrix of size LDA-by-N, if trans == Pasma[Conj]Trans.
 *         Matrix of size LDA-by-M, if trans == PasmaNoTrans.
 *
 * @param[in] LDA
 *         The leading dimension of the array A.
 *         LDA >= max(1,M), if trans == Pasma[Conj]Trans.
 *         LDA >= max(1,N), if trans == PasmaNoTrans.
 *
 * @param[out] B
 *         Matrix of size LDB-by-M, if trans == Pasma[Conj]Trans.
 *         Matrix of size LDB-by-N, if trans == PasmaNoTrans.
 *
 * @param[in] LDB
 *         The leading dimension of the array B.
 *         LDB >= max(1,N), if trans == Pasma[Conj]Trans.
 *         LDB >= max(1,M), if trans == PasmaNoTrans.
 *
 *
 *******************************************************************************
 *
 * @return
 *         \retval PLASMA_SUCCESS successful exit
 *         \retval <0 if INFO = -k, the k-th argument had an illegal value
 *
 ******************************************************************************/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_dlatro = PCORE_dlatro
#define CORE_dlatro PCORE_dlatro
#endif
int CORE_dlatro(PLASMA_enum uplo, PLASMA_enum trans,
                int M, int N,
                const double *A, int LDA,
                      double *B, int LDB)
{
    int i, j;

    /* Check input arguments */
    if ((uplo != PlasmaUpper) && (uplo != PlasmaLower) && (uplo != PlasmaUpperLower) ) {
        coreblas_error(1, "Illegal value of uplo");
        return -1;
    }
    if ((trans != PlasmaTrans) && (trans != PlasmaNoTrans) && (trans != PlasmaTrans) ) {
        coreblas_error(2, "Illegal value of trans");
        return -2;
    }
    if (M < 0) {
        coreblas_error(3, "Illegal value of M");
        return -3;
    }
    if (N < 0) {
        coreblas_error(4, "Illegal value of N");
        return -4;
    }
    if ( (LDA < max(1,M)) && (M > 0) ) {
        coreblas_error(6, "Illegal value of LDA");
        return -6;
    }
    if ( (LDB < max(1,N)) && (N > 0) ) {
        coreblas_error(8, "Illegal value of LDB");
        return -8;
    }

    if (trans == PlasmaNoTrans) {
        CORE_dlacpy(uplo, M, N, A, LDA, B, LDB);
    }
    else {
        if (trans == PlasmaTrans) {
            if(uplo == PlasmaUpper) {
                for(j=0; j<N; j++)
                    for(i=0; i<min(j+1,M); i++)
                        B[j+i*LDB] = (A[i+j*LDA]);
            }
            else if(uplo == PlasmaLower) {
                for(j=0;j<N;j++)
                    for(i=j;i<M;i++)
                        B[j+i*LDB] = (A[i+j*LDA]);
            }
            else {
                for(j=0;j<N;j++)
                    for(i=0;i<M;i++)
                        B[j+i*LDB] = (A[i+j*LDA]);
            }
        }
        else {
            if(uplo==PlasmaUpper) {
                for(j=0;j<N;j++)
                    for(i=0;i<min(j+1,M);i++)
                        B[j+i*LDB] = A[i+j*LDA];
            }
            else if(uplo==PlasmaLower) {
                for(j=0;j<N;j++)
                    for(i=j;i<M;i++)
                        B[j+i*LDB] = A[i+j*LDA];
            }
            else {
                for(j=0;j<N;j++)
                    for(i=0;i<M;i++)
                        B[j+i*LDB] = A[i+j*LDA];
            }
        }
    }

    return PLASMA_SUCCESS;
}

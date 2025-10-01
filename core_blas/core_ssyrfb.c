/**
 *
 * @file core_ssyrfb.c
 *
 *  PLASMA core_blas kernel
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Hatem Ltaief
 * @date 2010-11-15
 * @generated s Fri Apr  1 11:02:30 2016
 *
 **/
#ifdef PLASMA_WITH_MKL
#include <mkl_lapacke.h>
#else
#include <lapacke.h>
#endif
#include "common.h"
#undef COMPLEX
#define REAL

/***************************************************************************//**
 *
 * @ingroup CORE_float
 *
 *  CORE_ssyrfb overwrites the symmetric complex N-by-N tile C with
 *
 *    Q**T*C*Q
 *
 *  where Q is a complex unitary matrix defined as the product of k
 *  elementary reflectors
 *
 *    Q = H(1) H(2) . . . H(k)
 *
 *  as returned by CORE_sgeqrt. Only PlasmaLower supported!
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *         @arg PlasmaLower : the upper part of the symmetric matrix C
 *                            is not referenced.
 *         @arg PlasmaUpper : the lower part of the symmetric matrix C
 *                            is not referenced (not supported).
 *
 * @param[in] n
 *         The number of rows/columns of the tile C.  N >= 0.
 *
 * @param[in] k
 *         The number of elementary reflectors whose product defines
 *         the matrix Q. K >= 0.
 *
 * @param[in] ib
 *         The inner-blocking size.  IB >= 0.
 *
 * @param[in] nb
 *         The blocking size.  NB >= 0.
 *
 * @param[in] A
 *         The i-th column must contain the vector which defines the
 *         elementary reflector H(i), for i = 1,2,...,k, as returned by
 *         CORE_sgeqrt in the first k columns of its array argument A.
 *
 * @param[in] lda
 *         The leading dimension of the array A. LDA >= max(1,N).
 *
 * @param[in] T
 *         The IB-by-K triangular factor T of the block reflector.
 *         T is upper triangular by block (economic storage);
 *         The rest of the array is not referenced.
 *
 * @param[in] ldt
 *         The leading dimension of the array T. LDT >= IB.
 *
 * @param[in,out] C
 *         On entry, the symmetric N-by-N tile C.
 *         On exit, C is overwritten by Q**T*C*Q.
 *
 * @param[in] ldc
 *         The leading dimension of the array C. LDC >= max(1,M).
 *
 * @param[in,out] WORK
 *         On exit, if INFO = 0, WORK(1) returns the optimal LDWORK.
 *
 * @param[in] ldwork
 *         The dimension of the array WORK. LDWORK >= max(1,N);
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *
 ******************************************************************************/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_ssyrfb = PCORE_ssyrfb
#define CORE_ssyrfb PCORE_ssyrfb
#define CORE_sormlq PCORE_sormlq
#define CORE_sormqr PCORE_sormqr
int  CORE_sormlq(PLASMA_enum side, PLASMA_enum trans,
                 int M, int N, int IB, int K,
                 const float *V, int LDV,
                 const float *T, int LDT,
                 float *C, int LDC,
                 float *WORK, int LDWORK);
int  CORE_sormqr(PLASMA_enum side, PLASMA_enum trans,
                 int M, int N, int K, int IB,
                 const float *V, int LDV,
                 const float *T, int LDT,
                 float *C, int LDC,
                 float *WORK, int LDWORK);
#endif
int CORE_ssyrfb( PLASMA_enum uplo, int n,
                 int k, int ib, int nb,
                 const float *A, int lda,
                 const float *T, int ldt,
                 float *C, int ldc,
                 float *WORK, int ldwork )
{
    float tmp;
    int i, j;

    /* Check input arguments */
    if ((uplo != PlasmaUpper) && (uplo != PlasmaLower)) {
        coreblas_error(1, "Illegal value of uplo");
        return -1;
    }
    if (n < 0) {
        coreblas_error(2, "Illegal value of n");
        return -2;
    }
    if (k < 0) {
        coreblas_error(3, "Illegal value of k");
        return -3;
    }
    if (ib < 0) {
        coreblas_error(4, "Illegal value of ib");
        return -4;
    }
    if (nb < 0) {
        coreblas_error(5, "Illegal value of nb");
        return -5;
    }
    if ( (lda < max(1,n)) && (n > 0) ) {
        coreblas_error(7, "Illegal value of lda");
        return -7;
    }
    if ( (ldt < max(1,ib)) && (ib > 0) ) {
        coreblas_error(9, "Illegal value of ldt");
        return -9;
    }
    if ( (ldc < max(1,n)) && (n > 0) ) {
        coreblas_error(11, "Illegal value of ldc");
        return -11;
    }

    if (uplo == PlasmaLower) {
        /* Rebuild the symmetric block: WORK <- C */
        for (j = 0; j < n; j++) {
            *(WORK + j + j * ldwork) =  *(C + ldc*j + j);
            for (i = j+1; i < n; i++){
                tmp = *(C + i + j*ldc);
                *(WORK + i + j * ldwork) = tmp;
                *(WORK + j + i * ldwork) = ( tmp );
            }
        }

        /* Left */
        CORE_sormqr(PlasmaLeft, PlasmaTrans, n, n, k, ib,
                    A, lda, T, ldt, WORK, ldwork, WORK+nb*ldwork, ldwork);
        /* Right */
        CORE_sormqr(PlasmaRight, PlasmaNoTrans, n, n, k, ib,
                    A, lda, T, ldt, WORK, ldwork, WORK+nb*ldwork, ldwork);

        /*
         * Copy back the final result to the lower part of C
         */
        LAPACKE_slacpy_work( LAPACK_COL_MAJOR, lapack_const(PlasmaLower), n, n, WORK, ldwork, C, ldc );
    }
    else {
        /* Rebuild the symmetric block: WORK <- C */
        for (j = 0; j < n; j++) {
            for (i = 0; i < j; i++){
                tmp = *(C + i + j*ldc);
                *(WORK + i + j * ldwork) = tmp;
                *(WORK + j + i * ldwork) = ( tmp );
            }
            *(WORK + j + j * ldwork) =  *(C + ldc*j + j);
        }

        /* Right */
        CORE_sormlq(PlasmaRight, PlasmaTrans, n, n, k, ib,
                    A, lda, T, ldt, WORK, ldwork, WORK+nb*ldwork, ldwork);
        /* Left */
        CORE_sormlq(PlasmaLeft, PlasmaNoTrans, n, n, k, ib,
                    A, lda, T, ldt, WORK, ldwork, WORK+nb*ldwork, ldwork);

        /*
         * Copy back the final result to the upper part of C
         */
        LAPACKE_slacpy_work( LAPACK_COL_MAJOR, lapack_const(PlasmaUpper), n, n, WORK, ldwork, C, ldc );
    }
    return 0;
}

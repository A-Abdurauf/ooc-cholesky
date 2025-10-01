/**
 *
 * @file core_ctstrf.c
 *
 *  PLASMA core_blas kernel
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Hatem Ltaief
 * @author Mathieu Faverge
 * @author Jakub Kurzak
 * @date 2010-11-15
 * @generated c Fri Apr  1 11:02:28 2016
 *
 **/
#include "common.h"
#include <cblas.h>
#include <math.h>

/***************************************************************************//**
 *
 * @ingroup CORE_PLASMA_Complex32_t
 *
 *  CORE_ctstrf computes an LU factorization of a complex matrix formed
 *  by an upper triangular NB-by-N tile U on top of a M-by-N tile A
 *  using partial pivoting with row interchanges.
 *
 *  This is the right-looking Level 2.5 BLAS version of the algorithm.
 *
 *******************************************************************************
 *
 * @param[in] M
 *         The number of rows of the tile A.  M >= 0.
 *
 * @param[in] N
 *         The number of columns of the tile A.  N >= 0.
 *
 * @param[in] IB
 *         The inner-blocking size.  IB >= 0.
 *
 * @param[in] NB
 *
 * @param[in,out] U
 *         On entry, the NB-by-N upper triangular tile.
 *         On exit, the new factor U from the factorization
 *
 * @param[in] LDU
 *         The leading dimension of the array U.  LDU >= max(1,NB).
 *
 * @param[in,out] A
 *         On entry, the M-by-N tile to be factored.
 *         On exit, the factor L from the factorization
 *
 * @param[in] LDA
 *         The leading dimension of the array A.  LDA >= max(1,M).
 *
 * @param[in,out] L
 *         On entry, the IB-by-N lower triangular tile.
 *         On exit, the interchanged rows form the tile A in case of pivoting.
 *
 * @param[in] LDL
 *         The leading dimension of the array L.  LDL >= max(1,IB).
 *
 * @param[out] IPIV
 *         The pivot indices; for 1 <= i <= min(M,N), row i of the
 *         tile U was interchanged with row IPIV(i) of the tile A.
 *
 * @param[in,out] WORK
 *
 * @param[in] LDWORK
 *         The leading dimension of the array WORK.
 *
 * @param[out] INFO
 *
 *******************************************************************************
 *
 * @return
 *         \retval PLASMA_SUCCESS successful exit
 *         \retval <0 if INFO = -k, the k-th argument had an illegal value
 *         \retval >0 if INFO = k, U(k,k) is exactly zero. The factorization
 *              has been completed, but the factor U is exactly
 *              singular, and division by zero will occur if it is used
 *              to solve a system of equations.
 *
 ******************************************************************************/

#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_ctstrf = PCORE_ctstrf
#define CORE_ctstrf PCORE_ctstrf
#define CORE_cssssm PCORE_cssssm
int  CORE_cssssm(int M1, int N1, int M2, int N2, int K, int IB,
                 PLASMA_Complex32_t *A1, int LDA1,
                 PLASMA_Complex32_t *A2, int LDA2,
                 const PLASMA_Complex32_t *L1, int LDL1,
                 const PLASMA_Complex32_t *L2, int LDL2,
                 const int *IPIV);
#endif
int CORE_ctstrf(int M, int N, int IB, int NB,
                PLASMA_Complex32_t *U, int LDU,
                PLASMA_Complex32_t *A, int LDA,
                PLASMA_Complex32_t *L, int LDL,
                int *IPIV,
                PLASMA_Complex32_t *WORK, int LDWORK,
                int *INFO)
{
    static PLASMA_Complex32_t zzero = 0.0;
    static PLASMA_Complex32_t mzone =-1.0;

    PLASMA_Complex32_t alpha;
    int i, j, ii, sb;
    int im, ip;

    /* Check input arguments */
    *INFO = 0;
    if (M < 0) {
        coreblas_error(1, "Illegal value of M");
        return -1;
    }
    if (N < 0) {
        coreblas_error(2, "Illegal value of N");
        return -2;
    }
    if (IB < 0) {
        coreblas_error(3, "Illegal value of IB");
        return -3;
    }
    if ((LDU < max(1,NB)) && (NB > 0)) {
        coreblas_error(6, "Illegal value of LDU");
        return -6;
    }
    if ((LDA < max(1,M)) && (M > 0)) {
        coreblas_error(8, "Illegal value of LDA");
        return -8;
    }
    if ((LDL < max(1,IB)) && (IB > 0)) {
        coreblas_error(10, "Illegal value of LDL");
        return -10;
    }

    /* Quick return */
    if ((M == 0) || (N == 0) || (IB == 0))
        return PLASMA_SUCCESS;

    /* Set L to 0 */
    memset(L, 0, LDL*N*sizeof(PLASMA_Complex32_t));

    ip = 0;
    for (ii = 0; ii < N; ii += IB) {
        sb = min(N-ii, IB);

        for (i = 0; i < sb; i++) {
            im = cblas_icamax(M, &A[LDA*(ii+i)], 1);
            IPIV[ip] = ii+i+1;

            if (cabsf(A[LDA*(ii+i)+im]) > cabsf(U[LDU*(ii+i)+ii+i])) {
                /*
                 * Swap behind.
                 */
                cblas_cswap(i, &L[LDL*ii+i], LDL, &WORK[im], LDWORK );
                /*
                 * Swap ahead.
                 */
                cblas_cswap(sb-i, &U[LDU*(ii+i)+ii+i], LDU, &A[LDA*(ii+i)+im], LDA );
                /*
                 * Set IPIV.
                 */
                IPIV[ip] = NB + im + 1;

                for (j = 0; j < i; j++) {
                    A[LDA*(ii+j)+im] = zzero;
                }
            }

            if ((*INFO == 0) && (cabsf(A[LDA*(ii+i)+im]) == zzero)
                && (cabsf(U[LDU*(ii+i)+ii+i]) == zzero)) {
                *INFO = ii+i+1;
            }

            alpha = ((PLASMA_Complex32_t)1. / U[LDU*(ii+i)+ii+i]);
            cblas_cscal(M, CBLAS_SADDR(alpha), &A[LDA*(ii+i)], 1);
            cblas_ccopy(M, &A[LDA*(ii+i)], 1, &WORK[LDWORK*i], 1);
            cblas_cgeru(
                CblasColMajor, M, sb-i-1,
                CBLAS_SADDR(mzone), &A[LDA*(ii+i)], 1,
                &U[LDU*(ii+i+1)+ii+i], LDU,
                &A[LDA*(ii+i+1)], LDA );
            ip = ip+1;
        }
        /*
         * Apply the subpanel to the rest of the panel.
         */
        if(ii+i < N) {
            for(j = ii; j < ii+sb; j++) {
                if (IPIV[j] <= NB) {
                    IPIV[j] = IPIV[j] - ii;
                }
            }

            CORE_cssssm(
                NB, N-(ii+sb), M, N-(ii+sb), sb, sb,
                &U[LDU*(ii+sb)+ii], LDU,
                &A[LDA*(ii+sb)], LDA,
                &L[LDL*ii], LDL,
                WORK, LDWORK, &IPIV[ii]);

            for(j = ii; j < ii+sb; j++) {
                if (IPIV[j] <= NB) {
                    IPIV[j] = IPIV[j] + ii;
                }
            }
        }
    }
    return PLASMA_SUCCESS;
}

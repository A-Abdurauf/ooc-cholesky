/**
 *
 * @file core_sttqrt.c
 *
 *  PLASMA core_blas kernel
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Hatem Ltaief
 * @author Mathieu Faverge
 * @author Dulceneia Becker
 * @date 2010-11-15
 * @generated s Fri Apr  1 11:02:27 2016
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
 *  CORE_sttqrt computes a QR factorization of a rectangular matrix
 *  formed by coupling a complex N-by-N upper triangular tile A1
 *  on top of a complex M-by-N upper trapezoidal tile A2:
 *
 *    | A1 | = Q * R
 *    | A2 |
 *
 *  The tile Q is represented as a product of elementary reflectors
 *
 *    Q = H(1) H(2) . . . H(k), where k = min(M,N).
 *
 *  Each H(i) has the form
 *
 *    H(i) = I - tau * v * v'
 *
 *  where tau is a complex scalar, and v is a complex vector with
 *  v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A2(1:m,i),
 *  and tau in TAU(i).
 *
 *******************************************************************************
 *
 * @param[in] M
 *         The number of rows of the tile A2.  M >= 0.
 *
 * @param[in] N
 *         The number of columns of the tile A1 and A2.  N >= 0.
 *
 * @param[in] IB
 *         The inner-blocking size.  IB >= 0.
 *
 * @param[in,out] A1
 *         On entry, the N-by-N tile A1.
 *         On exit, the elements on and above the diagonal of the array
 *         contain the N-by-N upper trapezoidal tile R;
 *         the elements below the diagonal are not referenced.
 *
 * @param[in] LDA1
 *         The leading dimension of the array A1.  LDA1 >= max(1,N).
 *
 * @param[in,out] A2
 *         On entry, the M-by-N upper triangular tile A2.
 *         On exit, the elements on and above the diagonal of the array
 *         with the array TAU, represent
 *         the unitary tile Q as a product of elementary reflectors
 *         (see Further Details).
 *
 * @param[in] LDA2
 *         The leading dimension of the array A2.  LDA2 >= max(1,M).
 *
 * @param[out] T
 *         The IB-by-N triangular factor T of the block reflector.
 *         T is upper triangular by block (economic storage);
 *         The rest of the array is not referenced.
 *
 * @param[in] LDT
 *         The leading dimension of the array T. LDT >= IB.
 *
 * @param[out] TAU
 *         The scalar factors of the elementary reflectors (see Further
 *         Details).
 *
 * @param[in,out] WORK
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *
 ******************************************************************************/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_sttqrt = PCORE_sttqrt
#define CORE_sttqrt PCORE_sttqrt
/* Trick to call the version without tracing */
#define CORE_slaset PCORE_slaset
void
CORE_slaset(PLASMA_enum uplo, int n1, int n2,
            float alpha, float beta,
            float *tileA, int ldtilea);
#define CORE_spemv PCORE_spemv
int
CORE_spemv(PLASMA_enum trans, int storev,
           int M, int N, int L,
           float ALPHA,
           const float *A, int LDA,
           const float *X, int INCX,
           float BETA,
           float *Y, int INCY,
           float *WORK);
#endif
int CORE_sttqrt(int M, int N, int IB,
                float *A1, int LDA1,
                float *A2, int LDA2,
                float *T, int LDT,
                float *TAU, float *WORK)
{
    static int                ione  = 1;
    static float zone  = 1.0;
    static float zzero = 0.0;

    float alpha;
    int i, j, l, ii, sb, mi, ni;

    /* Check input arguments */
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
    if ((LDA2 < max(1,M)) && (M > 0)) {
        coreblas_error(7, "Illegal value of LDA2");
        return -7;
    }

    /* Quick return */
    if ((M == 0) || (N == 0) || (IB == 0))
        return PLASMA_SUCCESS;

    /* TODO: Need to check why some cases require
     *  this to not have uninitialized values */
    CORE_slaset( PlasmaUpperLower, IB, N,
                 0., 0., T, LDT);

    for (ii = 0; ii < N; ii += IB) {
        sb = min(N-ii, IB);
        for (i = 0; i < sb; i++) {
            j  = ii + i;
            mi = min( j + 1, M );
            ni = sb-i-1;

            /*
             * Generate elementary reflector H( II*IB+I ) to annihilate
             * A( II*IB+I:mi, II*IB+I ).
             */
            LAPACKE_slarfg_work(
                    mi+1, &A1[LDA1*j+j], &A2[LDA2*j], ione, &TAU[j]);

            if (ni > 0) {
                /*
                 * Apply H( II*IB+I ) to A( II*IB+I:M, II*IB+I+1:II*IB+IB ) from the left.
                 */
                cblas_scopy(
                    ni,
                    &A1[LDA1*(j+1)+j], LDA1,
                    WORK, 1);

#ifdef COMPLEX
                LAPACKE_slacgv_work(ni, WORK, 1);
#endif
                cblas_sgemv(
                    CblasColMajor, (CBLAS_TRANSPOSE)PlasmaTrans,
                    mi, ni,
                    (zone), &A2[LDA2*(j+1)], LDA2,
                                       &A2[LDA2*j],     1,
                    (zone), WORK,            1);
#ifdef COMPLEX
                LAPACKE_slacgv_work(ni, WORK, 1);
#endif
                alpha = -(TAU[j]);
                cblas_saxpy(
                    ni, (alpha),
                    WORK, 1,
                    &A1[LDA1*(j+1)+j], LDA1);
#ifdef COMPLEX
                LAPACKE_slacgv_work(ni, WORK, 1);
#endif
                cblas_sger(
                    CblasColMajor, mi, ni,
                    (alpha), &A2[LDA2*j], 1,
                    WORK, 1,
                    &A2[LDA2*(j+1)], LDA2);
            }

            /*
             * Calculate T
             *
             * T(0:i-1, j) = alpha * A2(0:M-1, ii:j-1)' * A2(0:M-1, j)
             */

            if ( i > 0 ) {

                l = min(i, max(0, M-ii));
                alpha = -(TAU[j]);

                CORE_spemv(
                        PlasmaTrans, PlasmaColumnwise,
                        min(j, M), i, l,
                        alpha, &A2[LDA2*ii], LDA2,
                               &A2[LDA2*j],  1,
                        zzero, &T[LDT*j],    1,
                        WORK);

                /* T(0:i-1, j) = T(0:i-1, ii:j-1) * T(0:i-1, j) */
                cblas_strmv(
                        CblasColMajor, (CBLAS_UPLO)PlasmaUpper,
                        (CBLAS_TRANSPOSE)PlasmaNoTrans,
                        (CBLAS_DIAG)PlasmaNonUnit,
                        i, &T[LDT*ii], LDT,
                           &T[LDT*j], 1);

            }

            T[LDT*j+i] = TAU[j];
        }

        /* Apply Q' to the rest of the matrix to the left  */
        if (N > ii+sb) {
            mi = min(ii+sb, M);
            ni = N-(ii+sb);
            l  = min(sb, max(0, mi-ii));
            CORE_sparfb(
                PlasmaLeft, PlasmaTrans,
                PlasmaForward, PlasmaColumnwise,
                IB, ni, mi, ni, sb, l,             //replaced sb by IB
                &A1[LDA1*(ii+sb)+ii], LDA1,
                &A2[LDA2*(ii+sb)], LDA2,
                &A2[LDA2*ii], LDA2,
                &T[LDT*ii], LDT,
                WORK, sb);
        }
    }
    return PLASMA_SUCCESS;
}

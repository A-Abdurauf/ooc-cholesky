/**
 *
 * @file core_sttlqt.c
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
 *  CORE_sttlqt computes a LQ factorization of a rectangular matrix
 *  formed by coupling side-by-side a complex M-by-M lower triangular tile A1
 *  and a complex M-by-N lower triangular tile A2:
 *
 *    | A1 A2 | = L * Q
 *
 *  The tile Q is represented as a product of elementary reflectors
 *
 *    Q = H(k)' . . . H(2)' H(1)', where k = min(M,N).
 *
 *  Each H(i) has the form
 *
 *    H(i) = I - tau * v * v'
 *
 *  where tau is a complex scalar, and v is a complex vector with
 *  v(1:i-1) = 0 and v(i) = 1; g(v(i+1:n)) is stored on exit in
 *  A2(i,1:n), and tau in TAU(i).
 *
 *******************************************************************************
 *
 * @param[in] M
 *         The number of rows of the tile A1 and A2. M >= 0.
 *         The number of columns of the tile A1.
 *
 * @param[in] N
 *         The number of columns of the tile A2. N >= 0.
 *
 * @param[in] IB
 *         The inner-blocking size.  IB >= 0.
 *
 * @param[in,out] A1
 *         On entry, the M-by-M tile A1.
 *         On exit, the elements on and below the diagonal of the array
 *         contain the M-by-M lower trapezoidal tile L;
 *         the elements above the diagonal are not referenced.
 *
 * @param[in] LDA1
 *         The leading dimension of the array A1.  LDA1 >= max(1,N).
 *
 * @param[in,out] A2
 *         On entry, the M-by-N lower triangular tile A2.
 *         On exit, the elements on and below the diagonal of the array
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
#pragma weak CORE_sttlqt = PCORE_sttlqt
#define CORE_sttlqt PCORE_sttlqt
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
int CORE_sttlqt(int M, int N, int IB,
                float *A1, int LDA1,
                float *A2, int LDA2,
                float *T, int LDT,
                float *TAU, float *WORK)
{
    static float zone  = 1.0;
    static float zzero = 0.0;
#ifdef COMPLEX
    static int                ione  = 1;
#endif

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

    for(ii = 0; ii < M; ii += IB) {
        sb = min(M-ii, IB);
        for(i = 0; i < sb; i++) {
            j  = ii + i;
            mi = sb-i-1;
            ni = min( j + 1, N);
            /*
             * Generate elementary reflector H( II*IB+I ) to annihilate A( II*IB+I, II*IB+I:M ).
             */
#ifdef COMPLEX
            LAPACKE_slacgv_work(ni, &A2[j], LDA2);
            LAPACKE_slacgv_work(ione, &A1[LDA1*j+j], LDA1);
#endif
            LAPACKE_slarfg_work(ni+1, &A1[LDA1*j+j], &A2[j], LDA2, &TAU[j]);

            if (mi > 0) {
                /*
                 * Apply H( j-1 ) to A( j:II+IB-1, j-1:M  ) from the right.
                 */
                cblas_scopy(
                    mi,
                    &A1[LDA1*j+(j+1)], 1,
                    WORK, 1);

                cblas_sgemv(
                    CblasColMajor, (CBLAS_TRANSPOSE)PlasmaNoTrans,
                    mi, ni,
                    (zone), &A2[j+1], LDA2,
                    &A2[j], LDA2,
                    (zone), WORK, 1);

                alpha = -(TAU[j]);
                cblas_saxpy(
                    mi, (alpha),
                    WORK, 1,
                    &A1[LDA1*j+j+1], 1);

                cblas_sger(
                    CblasColMajor, mi, ni,
                    (alpha), WORK, 1,
                    &A2[j], LDA2,
                    &A2[j+1], LDA2);
            }

            /*
             * Calculate T.
             */

            if (i > 0 ) {

                l = min(i, max(0, N-ii));
                alpha = -(TAU[j]);

                CORE_spemv(
                        PlasmaNoTrans, PlasmaRowwise,
                        i , min(j, N), l,
                        alpha, &A2[ii], LDA2,
                        &A2[j], LDA2,
                        zzero, &T[LDT*j], 1,
                        WORK);

                /* T(0:i-1, j) = T(0:i-1, ii:j-1) * T(0:i-1, j) */
                cblas_strmv(
                        CblasColMajor, (CBLAS_UPLO)PlasmaUpper,
                        (CBLAS_TRANSPOSE)PlasmaNoTrans,
                        (CBLAS_DIAG)PlasmaNonUnit,
                        i, &T[LDT*ii], LDT,
                        &T[LDT*j], 1);

            }

#ifdef COMPLEX
            LAPACKE_slacgv_work(ni, &A2[j], LDA2 );
            LAPACKE_slacgv_work(ione, &A1[LDA1*j+j], LDA1 );
#endif

            T[LDT*j+i] = TAU[j];
        }

        /* Apply Q to the rest of the matrix to the right */
        if (M > ii+sb) {
            mi = M-(ii+sb);
            ni = min(ii+sb, N);
            l  = min(sb, max(0, ni-ii));
            CORE_sparfb(
                PlasmaRight, PlasmaNoTrans,
                PlasmaForward, PlasmaRowwise,
                mi, IB, mi, ni, sb, l,
                &A1[LDA1*ii+ii+sb], LDA1,
                &A2[ii+sb], LDA2,
                &A2[ii], LDA2,
                &T[LDT*ii], LDT,
                WORK, M);

        }
    }
    return PLASMA_SUCCESS;
}


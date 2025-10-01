/**
 *
 * @file core_slarfy.c
 *
 *  PLASMA core_blas kernel
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Azzam Haidar
 * @date 2011-05-15
 * @generated s Fri Apr  1 11:02:31 2016
 *
 **/
#include <math.h>
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
 *  CORE_slarfy applies an elementary reflector, or Householder matrix, H,
 *  to a N-by-N hermitian matrix C, from both the left and the right.
 *
 *  H is represented in the form
 *
 *     H = I - tau * v * v'
 *
 *  where  tau  is a scalar and  v  is a vector.
 *
 *  If tau is zero, then H is taken to be the unit matrix.
 *
 *******************************************************************************
 *
 * @param[in] N
 *          The number of rows and columns of the matrix C.  N >= 0.
 *
 * @param[in,out] A
 *          REAL array, dimension (LDA, N)
 *          On entry, the Hermetian matrix A.
 *          On exit, A is overwritten by H * A * H'.
 *
 * @param[in] LDA
 *         The leading dimension of the array A.  LDA >= max(1,N).
 *
 * @param[in] V
 *          The vector V that contains the Householder reflectors.
 *
 * @param[in] TAU
 *          The value tau.
 *
 * @param[out] WORK
 *          Workspace.
 *
 ******************************************************************************/
void
CORE_slarfy(int N,
            float *A, int LDA,
            const float *V,
            const float *TAU,
            float *WORK)
{
    static float zzero =  0.0;
    static float zmone = -1.0;

    int j;
    float dtmp;

    /* Compute dtmp = X'*V */
    /* X = AVtau */
    cblas_ssymv(CblasColMajor, CblasLower,
                N, (*TAU), A, LDA,
                V, 1, (zzero), WORK, 1);

    /* cblas_sdotc_sub(N, WORK, 1, V, 1, &dtmp);*/
    dtmp = 0.;
    for (j = 0; j < N ; j++)
        dtmp = dtmp + (WORK[j]) * V[j];

    /* Compute 1/2 X'*V*t = 1/2*dtmp*tau  */
    dtmp = -dtmp * 0.5 * (*TAU);

   /* Compute W=X-1/2VX'Vt = X - dtmp*V */
    cblas_saxpy(N, (dtmp),
                V, 1, WORK, 1);

    /*
     * Performs the symmetric rank 2 operation
     *    A := alpha*x*y' + alpha*y*x' + A
     */
    cblas_ssyr2(CblasColMajor, CblasLower, N,
                (zmone), WORK, 1,
                                    V,    1,
                                    A,    LDA);

    return;
}
#undef COMPLEX

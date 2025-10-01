/**
 *
 * @file sgecon.c
 *
 *  PLASMA computational routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Jakub Kurzak
 * @author Ichitaro Yamazaki
 * @author Mathieu Faverge
 * @date 2012-10-04
 * @generated s Fri Apr  1 11:02:56 2016
 *
 **/
#ifdef PLASMA_WITH_MKL
#include <mkl_lapacke.h>
#else
#include <lapacke.h>
#endif
#include "common.h"
#include <math.h>
#undef COMPLEX
#define REAL

/***************************************************************************//**
 *
 * @ingroup float
 *
 *  PLASMA_sgecon - estimates the reciprocal of the condition number
 *  of a general complex matrix A, in either the 1-norm or the infinity-norm,
 *  using the LU factorization computed by PLASMA_sgetrf().
 *
 *  An estimate is obtained for norm(inv(A)), and the reciprocal of the condition
 *  number is computed as
 *
 *    \f[ rcond = \frac{1}{\|\|A\-\| \times \-\|A^{-1}\|\|} \f]
 *
 *******************************************************************************
 *
 * @param[in] norm
 *          Specifies whether the 1-norm condition number
 *          or the infinity-norm condition number is required:
 *          = PlasmaOneNorm: One norm
 *          = PlasmaInfNorm: Infinity norm
 *
 * @param[in] N
 *          The order of the matrix A. N >= 0.
 *
 * @param[in] A
 *          The N-by-N matrix A.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,N).
 *
 * @param[in] Anorm
 *          If norm = PlasmaOneNorm, the 1-norm of the original matrix A.
 *          If norm = PlasmaInfNorm, the infinity-norm of the original matrix A.
 *
 * \param[out] rcond
 *          The reciprocal of the condition number of the matrix A,
 *          computed as stated above.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *
 *******************************************************************************
 *
 * @sa PLASMA_sgecon_Tile
 * @sa PLASMA_sgecon_Tile_Async
 * @sa PLASMA_cgecon
 * @sa PLASMA_dgecon
 * @sa PLASMA_sgecon
 *
 ******************************************************************************/
int PLASMA_sgecon(PLASMA_enum norm, int N,
                  float *A, int LDA, float Anorm, float *rcond)
{
    int NB;
    int status;
    plasma_context_t *plasma;
    PLASMA_sequence *sequence = NULL;
    PLASMA_request request = PLASMA_REQUEST_INITIALIZER;
    PLASMA_desc descA;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_sgecon", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    /* Check input arguments */
    if (norm != PlasmaOneNorm && norm != PlasmaInfNorm) {
        plasma_error("PLASMA_sgecon", "illegal value of norm");
        return -1;
    }
    if (N < 0) {
        plasma_error("PLASMA_sgecon", "illegal value of N");
        return -2;
    }
    if (LDA < max(1,N)) {
        plasma_error("PLASMA_sgecon", "illegal value of LDA");
        return -4;
    }
    if (Anorm < 0.) {
        plasma_error("PLASMA_sgecon", "illegal value of Anorm");
        return -5;
    }

    /* Quick return */
    *rcond = (float)0.;
    if (N == 0) {
        *rcond = (float)1.;
        return PLASMA_SUCCESS;
    }
    else if (Anorm == 0.) {
        return PLASMA_SUCCESS;
    }

    /* Tune NB depending on M, N & NRHS; Set NBNB */
    status = plasma_tune(PLASMA_FUNC_SGESV, N, N, 0);
    if (status != PLASMA_SUCCESS) {
        plasma_error("PLASMA_sgecon", "plasma_tune() failed");
        return status;
    }
    /* Set NT */
    NB = PLASMA_NB;
    plasma_sequence_create(plasma, &sequence);

    if (PLASMA_TRANSLATION == PLASMA_OUTOFPLACE) {
        plasma_sooplap2tile(
            descA, A, NB, NB, LDA, N, 0, 0, N, N,
            sequence, &request, plasma_desc_mat_free(&(descA)));
    } else {
        plasma_siplap2tile(
            descA, A, NB, NB, LDA, N, 0, 0, N, N, sequence, &request);
    }

    /* Call the tile interface */
    PLASMA_sgecon_Tile_Async(norm, &descA, Anorm, rcond, sequence, &request);

    if (PLASMA_TRANSLATION == PLASMA_OUTOFPLACE) {
        plasma_dynamic_sync();
        plasma_desc_mat_free(&descA);
    } else {
        plasma_siptile2lap(descA, A, NB, NB, LDA, N, sequence, &request);
        plasma_dynamic_sync();
    }

    status = sequence->status;
    plasma_sequence_destroy(plasma, sequence);
    return status;
}

/***************************************************************************//**
 *
 * @ingroup float_Tile
 *
 *  PLASMA_sgecon_Tile - estimates the reciprocal of the condition number
 *  of a general complex matrix A, in either the 1-norm or the infinity-norm.
 *  Tile equivalent of PLASMA_sgecon().
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] norm
 *          Specifies whether the 1-norm condition number
 *          or the infinity-norm condition number is required:
 *          = PlasmaOneNorm: One norm
 *          = PlasmaInfNorm: Infinity norm
 *
 * @param[in] A
 *          The N-by-N matrix A.
 *
 * @param[in] Anorm
 *          If norm = PlasmaOneNorm, the 1-norm of the original matrix A.
 *          If norm = PlasmaInfNorm, the infinity-norm of the original matrix A.
 *
 * \param[out] rcond
 *          The reciprocal of the condition number of the matrix A,
 *          computed as stated above.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *
 *******************************************************************************
 *
 * @sa PLASMA_sgecon
 * @sa PLASMA_sgecon_Tile_Async
 * @sa PLASMA_cgecon_Tile
 * @sa PLASMA_dgecon_Tile
 * @sa PLASMA_sgecon_Tile
 *
 ******************************************************************************/
int PLASMA_sgecon_Tile(PLASMA_enum norm, PLASMA_desc *A, float Anorm, float *rcond)
{
    plasma_context_t *plasma;
    PLASMA_sequence *sequence = NULL;
    PLASMA_request request = PLASMA_REQUEST_INITIALIZER;
    int status;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_sgecon_Tile", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    plasma_sequence_create(plasma, &sequence);
    PLASMA_sgecon_Tile_Async(norm, A, Anorm, rcond, sequence, &request);
    plasma_dynamic_sync();
    status = sequence->status;
    plasma_sequence_destroy(plasma, sequence);
    return status;
}

/***************************************************************************//**
 *
 * @ingroup float_Tile_Async
 *
 *  PLASMA_sgecon_Tile_Async - estimates the reciprocal of the condition number
 *  of a general complex matrix A, in either the 1-norm or the infinity-norm.
 *  Non-blocking equivalent of PLASMA_sgecon_Tile().
 *  May return before the computation is finished.
 *  Allows for pipelining of operations at runtime.
 *
 *******************************************************************************
 *
 * @param[in] sequence
 *          Identifies the sequence of function calls that this call belongs to
 *          (for completion checks and exception handling purposes).
 *
 * @param[out] request
 *          Identifies this function call (for exception handling purposes).
 *
 *******************************************************************************
 *
 * @sa PLASMA_sgecon
 * @sa PLASMA_sgecon_Tile
 * @sa PLASMA_cgecon_Tile_Async
 * @sa PLASMA_dgecon_Tile_Async
 * @sa PLASMA_sgecon_Tile_Async
 *
 ******************************************************************************/
int PLASMA_sgecon_Tile_Async(PLASMA_enum norm, PLASMA_desc *A, float Anorm, float *rcond,
                             PLASMA_sequence *sequence, PLASMA_request *request)
{
    plasma_context_t *plasma;
    PLASMA_desc descA;
    PLASMA_desc descW;
    float *workN;
    float Ainvnorm;
    int kase, kase1;
    int isave[3], itrs = 0;
    int fallback = PLASMA_FALSE;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_sgecon_Tile_Async", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        plasma_fatal_error("PLASMA_sgecon_Tile_Async", "NULL sequence");
        return PLASMA_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        plasma_fatal_error("PLASMA_sgecon_Tile_Async", "NULL request");
        return PLASMA_ERR_UNALLOCATED;
    }
    /* Check sequence status */
    if (sequence->status == PLASMA_SUCCESS)
        request->status = PLASMA_SUCCESS;
    else
        return plasma_request_fail(sequence, request, PLASMA_ERR_SEQUENCE_FLUSHED);

    /* Check descriptors for correctness */
    if ( A->m != A->n ) {
        plasma_error("PLASMA_sgecon_Tile_Async", "invalid A descriptor (not square)");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }
    if (plasma_desc_check(A) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_sgecon_Tile_Async", "invalid first descriptor");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    } else {
        descA = *A;
    }
    /* Check input arguments */
    if (descA.nb != descA.mb) {
        plasma_error("PLASMA_sgecon_Tile_Async", "only square tiles supported");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }

    /* Quick return */
    *rcond = (float)0.;
    if (descA.m == 0) {
        *rcond = (float)1.;
        return PLASMA_SUCCESS;
    }
    else if (Anorm == 0.) {
        return PLASMA_SUCCESS;
    }

    /* Estimate the norm of inv(A). */
    Ainvnorm = (float)0.;
    if (norm == PlasmaOneNorm)
        kase1 = 1;
    else
        kase1 = 2;
    kase = 0;

#if defined(REAL)
    int *isgn = (int*)plasma_shared_alloc(plasma, descA.m, PlasmaInteger);
#endif

    workN = (float*)plasma_shared_alloc(plasma, descA.m, PlasmaRealFloat);
    plasma_sdesc_alloc( descW, descA.mb, descA.nb,
                        descA.m, 1, 0, 0, descA.m, 1, plasma_desc_mat_free(&(descW)));

    do {
        itrs ++;
#if defined(REAL)
        LAPACKE_slacn2_work( descA.m, workN, descW.mat, isgn, &Ainvnorm, &kase, isave);
#else
        LAPACKE_slacn2_work( descA.m, workN, descW.mat,       &Ainvnorm, &kase, isave);
#endif

#define FALLBACK
#ifdef  FALLBACK
        /*
         * Fall back to LAPACK
         */
        if( isnan(Ainvnorm) || isinf(Ainvnorm) || Ainvnorm > LAPACKE_slamch('O') ) {
            int info;
            float *Atmp = (float*)malloc(descA.m * descA.n * sizeof(float));

            plasma_sooptile2lap( descA, Atmp, descA.mb, descA.nb, descA.m, descA.n, sequence, request);
            plasma_dynamic_sync();
            info = LAPACKE_sgecon(LAPACK_COL_MAJOR, lapack_const(norm), descA.n, Atmp, descA.m, Anorm, rcond);
            free(Atmp);
            fallback = PLASMA_TRUE;
            sequence->status = info;
            kase = 0;
        }
#endif

        if (kase != 0) {
            if (kase == kase1) {
                /* Multiply by inv(L). */
                plasma_parallel_call_9(plasma_pstrsm,
                    PLASMA_enum, PlasmaLeft,
                    PLASMA_enum, PlasmaLower,
                    PLASMA_enum, PlasmaNoTrans,
                    PLASMA_enum, PlasmaUnit,
                    float, 1.0,
                    PLASMA_desc, descA,
                    PLASMA_desc, descW,
                    PLASMA_sequence*, sequence,
                    PLASMA_request*, request);

                /* Multiply by inv(U). */
                plasma_parallel_call_9(plasma_pstrsm,
                    PLASMA_enum, PlasmaLeft,
                    PLASMA_enum, PlasmaUpper,
                    PLASMA_enum, PlasmaNoTrans,
                    PLASMA_enum, PlasmaNonUnit,
                    float, 1.0,
                    PLASMA_desc, descA,
                    PLASMA_desc, descW,
                    PLASMA_sequence*, sequence,
                    PLASMA_request*, request);
            }
            else {
                /* Multiply by inv(U**T). */
                plasma_parallel_call_9(plasma_pstrsm,
                    PLASMA_enum, PlasmaLeft,
                    PLASMA_enum, PlasmaUpper,
                    PLASMA_enum, PlasmaTrans,
                    PLASMA_enum, PlasmaNonUnit,
                    float, 1.0,
                    PLASMA_desc, descA,
                    PLASMA_desc, descW,
                    PLASMA_sequence*, sequence,
                    PLASMA_request*, request);

                /* Multiply by inv(L**T). */
                plasma_parallel_call_9(plasma_pstrsm,
                    PLASMA_enum, PlasmaLeft,
                    PLASMA_enum, PlasmaLower,
                    PLASMA_enum, PlasmaTrans,
                    PLASMA_enum, PlasmaUnit,
                    float, 1.0,
                    PLASMA_desc, descA,
                    PLASMA_desc, descW,
                    PLASMA_sequence*, sequence,
                    PLASMA_request*, request);
            }
        }
        plasma_dynamic_sync();
    }
    while (kase != 0);

    /* Compute the estimate of the reciprocal condition number. */
    if ((Ainvnorm != 0.0) && (fallback == PLASMA_FALSE)) {
        *rcond = ((float)1.0 / Ainvnorm) / Anorm;
    }

#if defined(REAL)
    plasma_shared_free(plasma, isgn);
#endif
    plasma_shared_free(plasma, workN);
    plasma_desc_mat_free(&descW);

    return PLASMA_SUCCESS;
}

/**
 *
 * @file dlaswp.c
 *
 *  PLASMA computational routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @generated d Fri Apr  1 11:02:53 2016
 *
 **/
#include "common.h"

/***************************************************************************//**
 *
 * @ingroup double
 *
 *  PLASMA_dlaswp - performs a series of row interchanges on the matrix A.
 *  One row interchange is initiated for each of rows K1 through K2 of A.
 *
 *******************************************************************************
 *
 * @param[in] N
 *          The order of the matrix A.  N >= 0.
 *
 * @param[in] A
 *          The tile factors L and U from the factorization, computed by PLASMA_dgetrf.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,N).
 *
 * @param[in] K1
 *          The first element of IPIV for which a row interchange will
 *          be done.
 *
 * @param[in] K2
 *          The last element of IPIV for which a row interchange will
 *          be done.
 *
 * @param[in] IPIV
 *          The pivot indices from PLASMA_dgetrf.
 *
 * @param[in] INCX
 *          The increment between successive values of IPIV. If IPIV
 *          is negative, the pivots are applied in reverse order.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *          \return <0 if -i, the i-th argument had an illegal value
 *
 *******************************************************************************
 *
 * @sa PLASMA_dlaswp_Tile
 * @sa PLASMA_dlaswp_Tile_Async
 * @sa PLASMA_claswp
 * @sa PLASMA_dlaswp
 * @sa PLASMA_slaswp
 * @sa PLASMA_dgetrf
 *
 ******************************************************************************/
int PLASMA_dlaswp(int N, double *A, int LDA,
                  int K1, int K2, const int *IPIV, int INCX)
{
    int NB;
    int status;
    plasma_context_t *plasma;
    PLASMA_sequence *sequence = NULL;
    PLASMA_request request = PLASMA_REQUEST_INITIALIZER;
    PLASMA_desc descA;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_dlaswp", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    /* Check input arguments */
    if (N < 0) {
        plasma_error("PLASMA_dlaswp", "illegal value of N");
        return -1;
    }
    if (LDA < max(1, N)) {
        plasma_error("PLASMA_dlaswp", "illegal value of LDA");
        return -3;
    }

    /* Quick return */
    if ( N == 0 )
        return PLASMA_SUCCESS;

    /* Tune NB & IB depending on N & NRHS; Set NBNBSIZE */
    status = plasma_tune(PLASMA_FUNC_DGESV, LDA, N, N);
    if (status != PLASMA_SUCCESS) {
        plasma_error("PLASMA_dlaswp", "plasma_tune() failed");
        return status;
    }

    /* Set NT & NTRHS */
    NB = PLASMA_NB;

    plasma_sequence_create(plasma, &sequence);

    if ( PLASMA_TRANSLATION == PLASMA_OUTOFPLACE ) {
        plasma_dooplap2tile( descA, A, NB, NB, LDA, N, 0, 0, K2, N, sequence, &request,
                             plasma_desc_mat_free(&(descA)) );
    } else {
        plasma_diplap2tile( descA, A, NB, NB, LDA, N, 0, 0, K2, N,
                            sequence, &request);
    }

    /* Call the tile interface */
    PLASMA_dlaswp_Tile_Async(&descA, K1, K2, IPIV, INCX, sequence, &request);

    if ( PLASMA_TRANSLATION == PLASMA_OUTOFPLACE ) {
        plasma_dooptile2lap( descA, A, NB, NB, LDA, N,  sequence, &request);
        plasma_dynamic_sync();
        plasma_desc_mat_free(&descA);
    } else {
        plasma_diptile2lap( descA, A, NB, NB, LDA, N,  sequence, &request);
        plasma_dynamic_sync();
    }

    status = sequence->status;
    plasma_sequence_destroy(plasma, sequence);
    return status;
}

/***************************************************************************//**
 *
 * @ingroup double_Tile
 *
 *  PLASMA_dlaswp_Tile - performs a series of row interchanges on the matrix A.
 *  One row interchange is initiated for each of rows K1 through K2 of A.
 *  Tile equivalent of PLASMA_dlaswp().
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] A
 *          The tile factors L and U from the factorization, computed by PLASMA_dgetrf.
 *
 * @param[in] K1
 *          The first element of IPIV for which a row interchange will
 *          be done.
 *
 * @param[in] K2
 *          The last element of IPIV for which a row interchange will
 *          be done.
 *
 * @param[in] IPIV
 *          The pivot indices from PLASMA_dgetrf.
 *
 * @param[in] INCX
 *          The increment between successive values of IPIV. If IPIV
 *          is negative, the pivots are applied in reverse order.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa PLASMA_dlaswp
 * @sa PLASMA_dlaswp_Tile_Async
 * @sa PLASMA_claswp_Tile
 * @sa PLASMA_dlaswp_Tile
 * @sa PLASMA_slaswp_Tile
 * @sa PLASMA_dgetrf_Tile
 *
 ******************************************************************************/
int PLASMA_dlaswp_Tile(PLASMA_desc *A, int K1, int K2, const int *IPIV, int INCX)
{
    plasma_context_t *plasma;
    PLASMA_sequence *sequence = NULL;
    PLASMA_request request = PLASMA_REQUEST_INITIALIZER;
    int status;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_dlaswp_Tile", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    plasma_sequence_create(plasma, &sequence);
    PLASMA_dlaswp_Tile_Async(A, K1, K2, IPIV, INCX, sequence, &request);
    plasma_dynamic_sync();
    status = sequence->status;
    plasma_sequence_destroy(plasma, sequence);
    return status;
}

/***************************************************************************//**
 *
 * @ingroup double_Tile_Async
 *
 *  PLASMA_dlaswp_Tile_Async - performs a series of row interchanges
 *  on the matrix A.  One row interchange is initiated for each of
 *  rows K1 through K2 of A.
 *  Non-blocking equivalent of PLASMA_dlaswp_Tile().
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
 * @sa PLASMA_dlaswp
 * @sa PLASMA_dlaswp_Tile
 * @sa PLASMA_claswp_Tile_Async
 * @sa PLASMA_dlaswp_Tile_Async
 * @sa PLASMA_slaswp_Tile_Async
 * @sa PLASMA_dgetrf_Tile_Async
 *
 ******************************************************************************/
int PLASMA_dlaswp_Tile_Async(PLASMA_desc *A, int K1, int K2, const int *IPIV, int INCX,
                             PLASMA_sequence *sequence, PLASMA_request *request)
{
    PLASMA_desc descA;
    plasma_context_t *plasma;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_dlaswp_Tile", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        plasma_fatal_error("PLASMA_dlaswp_Tile", "NULL sequence");
        return PLASMA_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        plasma_fatal_error("PLASMA_dlaswp_Tile", "NULL request");
        return PLASMA_ERR_UNALLOCATED;
    }
    /* Check sequence status */
    if (sequence->status == PLASMA_SUCCESS)
        request->status = PLASMA_SUCCESS;
    else
        return plasma_request_fail(sequence, request, PLASMA_ERR_SEQUENCE_FLUSHED);

    /* Check descriptors for correctness */
    if (plasma_desc_check(A) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_dlaswp_Tile", "invalid first descriptor");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    } else {
        descA = *A;
    }

    if ( (K1 != 1) || (K2 != descA.m) ) {
        plasma_error("PLASMA_dlaswp_Tile", "invalid K1 or K2 (1..M is the only interval supported right now)");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }

    plasma_dynamic_call_3(
        plasma_pdbarrier_tl2pnl,
        PLASMA_desc, descA,
        PLASMA_sequence*, sequence,
        PLASMA_request*,  request);

    /* swap */
    plasma_dynamic_call_5(
        plasma_pdlaswp,
        PLASMA_desc, descA,
        const int *,       IPIV,
        int,         INCX,
        PLASMA_sequence*, sequence,
        PLASMA_request*,  request);

    plasma_dynamic_call_3(
        plasma_pdbarrier_pnl2tl,
        PLASMA_desc, descA,
        PLASMA_sequence*, sequence,
        PLASMA_request*,  request);

    return PLASMA_SUCCESS;
}

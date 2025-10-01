/**
 *
 * @file cgetrf_tntpiv.c
 *
 *  PLASMA computational routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Mathieu Faverge
 * @date 2009-11-15
 *
 * @generated c Fri Apr  1 11:02:56 2016
 *
 **/
#include <stdlib.h>
#include "common.h"

/***************************************************************************//**
 *
 * @ingroup PLASMA_Complex32_t
 *
 *  PLASMA_cgetrf_tntpiv - Computes an LU factorization of a general M-by-N
 *  matrix A using the tile LU algorithm with tournament pivoting based on LU
 *  decomposition at each round.
 *
 *******************************************************************************
 *
 * @param[in] M
 *          The number of rows of the matrix A. M >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrix A. N >= 0.
 *
 * @param[in,out] A
 *          On entry, the M-by-N matrix to be factored.
 *          On exit, the tile factors L and U from the factorization.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,M).
 *
 * @param[out] IPIV
 *          The pivot indices that define the permutations.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *          \retval >0 if i, U(i,i) is exactly zero. The factorization has been completed,
 *               but the factor U is exactly singular, and division by zero will occur
 *               if it is used to solve a system of equations.
 *
 *******************************************************************************
 *
 * @sa PLASMA_cgetrf
 * @sa PLASMA_cgetrf_nopiv
 * @sa PLASMA_cgetrf_incpiv
 * @sa PLASMA_cgetrf_tntpiv_Tile
 * @sa PLASMA_cgetrf_tntpiv_Tile_Async
 * @sa PLASMA_cgetrf_tntpiv
 * @sa PLASMA_dgetrf_tntpiv
 * @sa PLASMA_sgetrf_tntpiv
 *
 ******************************************************************************/
int PLASMA_cgetrf_tntpiv(int M, int N,
                         PLASMA_Complex32_t *A, int LDA,
                         int *IPIV)
{
    int NB;
    int status;
    int *Wpivot;
    PLASMA_desc descA, W;
    plasma_context_t *plasma;
    PLASMA_sequence *sequence = NULL;
    PLASMA_request request = PLASMA_REQUEST_INITIALIZER;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_cgetrf", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    /* Check input arguments */
    if (M < 0) {
        plasma_error("PLASMA_cgetrf", "illegal value of M");
        return -1;
    }
    if (N < 0) {
        plasma_error("PLASMA_cgetrf", "illegal value of N");
        return -2;
    }
    if (LDA < max(1, M)) {
        plasma_error("PLASMA_cgetrf", "illegal value of LDA");
        return -4;
    }
    /* Quick return */
    if (min(M, N) == 0)
        return PLASMA_SUCCESS;

    /* Tune NB & IB depending on M, N & NRHS; Set NBNBSIZE */
    status = plasma_tune(PLASMA_FUNC_CGESV, M, N, 0);
    if (status != PLASMA_SUCCESS) {
        plasma_error("PLASMA_cgetrf", "plasma_tune() failed");
        return status;
    }
    /* Set NT & NTRHS */
    NB = PLASMA_NB;

    plasma_sequence_create(plasma, &sequence);

    if ( PLASMA_TRANSLATION == PLASMA_OUTOFPLACE ) {
        plasma_cooplap2tile( descA, A, NB, NB, LDA, N, 0, 0, N, N,
                             sequence, &request,
                             plasma_desc_mat_free(&(descA)) );
    } else {
        plasma_ciplap2tile( descA, A, NB, NB, LDA, N, 0, 0, N, N,
                            sequence, &request);
    }

    /* Call the tile interface */
    PLASMA_Alloc_Workspace_cgetrf_tntpiv_Tile(&descA, &W, &Wpivot);
    PLASMA_cgetrf_tntpiv_Tile_Async(&descA, IPIV, &W, Wpivot, sequence, &request);

    if ( PLASMA_TRANSLATION == PLASMA_OUTOFPLACE ) {
        plasma_cooptile2lap( descA, A, NB, NB, LDA, N,  sequence, &request);
        plasma_dynamic_sync();
        plasma_desc_mat_free(&descA);
    } else {
        plasma_ciptile2lap( descA, A, NB, NB, LDA, N,  sequence, &request);
        plasma_dynamic_sync();
    }

    plasma_desc_mat_free(&W);
    free(Wpivot);

    status = sequence->status;
    plasma_sequence_destroy(plasma, sequence);

    return status;
}

/***************************************************************************//**
 *
 * @ingroup PLASMA_Complex32_t_Tile
 *
 *  PLASMA_cgetrf_tntpiv_Tile - Computes the tile LU factorization of a matrix.
 *  Tile equivalent of PLASMA_cgetrf_tntpiv().
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in,out] A
 *          On entry, the M-by-N matrix to be factored.
 *          On exit, the tile factors L and U from the factorization.
 *
 * @param[out] IPIV
 *          The pivot indices that define the permutations
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *          \retval >0 if i, U(i,i) is exactly zero. The factorization has been completed,
 *               but the factor U is exactly singular, and division by zero will occur
 *               if it is used to solve a system of equations.
 *
 *******************************************************************************
 *
 * @sa PLASMA_cgetrf_tntpiv
 * @sa PLASMA_cgetrf_tntpiv_Tile_Async
 * @sa PLASMA_cgetrf_tntpiv_Tile
 * @sa PLASMA_dgetrf_tntpiv_Tile
 * @sa PLASMA_sgetrf_tntpiv_Tile
 *
 ******************************************************************************/
int PLASMA_cgetrf_tntpiv_Tile(PLASMA_desc *A, int *IPIV)
{
    plasma_context_t *plasma;
    PLASMA_sequence *sequence = NULL;
    PLASMA_request request = PLASMA_REQUEST_INITIALIZER;
    PLASMA_desc W;
    int *Wpivot;
    int status;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_cgetrf_tntpiv_Tile", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    plasma_sequence_create(plasma, &sequence);
    PLASMA_Alloc_Workspace_cgetrf_tntpiv_Tile(A, &W, &Wpivot);
    PLASMA_cgetrf_tntpiv_Tile_Async(A, IPIV, &W, Wpivot, sequence, &request);
    plasma_dynamic_sync();
    plasma_desc_mat_free(&W);
    free(Wpivot);
    status = sequence->status;
    plasma_sequence_destroy(plasma, sequence);
    return status;
}

/***************************************************************************//**
 *
 * @ingroup PLASMA_Complex32_t_Tile_Async
 *
 *  PLASMA_cgetrf_tntpiv_Tile_Async - Computes the tile LU factorization of a
 *  matrix.  Non-blocking equivalent of PLASMA_cgetrf_tntpiv_Tile().
 *  May return before the computation is finished.
 *  Allows for pipelining of operations at runtime.
 *
 *******************************************************************************
 *
 * @param[in,out] A
 *          On entry, the M-by-N matrix to be factored.
 *          On exit, the tile factors L and U from the factorization.
 *
 * @param[out] IPIV
 *          The pivot indices that define the permutations
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
 * @sa PLASMA_cgetrf_tntpiv
 * @sa PLASMA_cgetrf_tntpiv_Tile
 * @sa PLASMA_cgetrf_tntpiv_Tile_Async
 * @sa PLASMA_dgetrf_tntpiv_Tile_Async
 * @sa PLASMA_sgetrf_tntpiv_Tile_Async
 *
 ******************************************************************************/
int PLASMA_cgetrf_tntpiv_Tile_Async(PLASMA_desc *A, int *IPIV,
                                    PLASMA_desc *W, int *Wpivot,
                                    PLASMA_sequence *sequence,
                                    PLASMA_request *request)
{
    PLASMA_desc descA, descW;
    plasma_context_t *plasma;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_cgetrf_tntpiv_Tile", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        plasma_fatal_error("PLASMA_cgetrf_tntpiv_Tile", "NULL sequence");
        return PLASMA_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        plasma_fatal_error("PLASMA_cgetrf_tntpiv_Tile", "NULL request");
        return PLASMA_ERR_UNALLOCATED;
    }
    if (PLASMA_TNT_MODE != PLASMA_TOURNAMENT_LU) {
        plasma_fatal_error("PLASMA_cgetrf_tntpiv_Tile", "Only PLASMA_TOURNAMENT_LU supported");
        return PLASMA_ERR_UNALLOCATED;
    }
    /* Check sequence status */
    if (sequence->status == PLASMA_SUCCESS)
        request->status = PLASMA_SUCCESS;
    else
        return plasma_request_fail(sequence, request, PLASMA_ERR_SEQUENCE_FLUSHED);

    /* Check descriptors for correctness */
    if (plasma_desc_check(A) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_cgetrf_tntpiv_Tile", "invalid first descriptor");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    } else {
        descA = *A;
    }
    /* Check descriptors for correctness */
    if (plasma_desc_check(W) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_cgetrf_tntpiv_Tile_Async", "invalid W descriptor");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    } else {
        descW = *W;
    }
    /* Check input arguments */
    if (descA.nb != descA.mb) {
        plasma_error("PLASMA_cgetrf_tntpiv_Tile", "only square tiles supported");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }

    plasma_dynamic_call_6(
        plasma_pcgetrf_tntpiv,
        PLASMA_desc,      descA,
        int*,             IPIV,
        PLASMA_desc,      descW,
        int*,             Wpivot,
        PLASMA_sequence*, sequence,
        PLASMA_request*,  request);

    return PLASMA_SUCCESS;
}

int PLASMA_Alloc_Workspace_cgetrf_tntpiv_Tile(PLASMA_desc *A,
                                              PLASMA_desc *W,
                                              int **Wpivot)
{
    plasma_context_t *plasma;
    int rndsize;
    int wm, tmp, nblevel = 1;

    plasma  = plasma_context_self();
    rndsize = PLASMA_TNT_SIZE;
    wm = ((A->mt + rndsize - 1) / rndsize )
        * rndsize * A->mb;

    /* Compute the maximum number of level */
    tmp = wm;
    while ( (tmp >> 2) > 1 ){
        nblevel++;
        tmp = tmp >> 2;
    }
    plasma_cdesc_alloc( *W, rndsize*A->mb, A->nb, wm, nblevel*A->nb,
                        0, 0, wm, nblevel*A->nb, plasma_desc_mat_free( W ) );

    (*Wpivot) = (int *)malloc(wm * nblevel * sizeof(int));

    return PLASMA_SUCCESS;
}

/**
 *
 * @file slaset.c
 *
 *  PLASMA computational routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @generated s Fri Apr  1 11:02:55 2016
 *
 **/
#include "common.h"

/***************************************************************************//**
 *
 * @ingroup float
 *
 *  PLASMA_slaset sets all or part of a two-dimensional matrix A.
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          Specifies the part of the matrix A to be set.
 *            = PlasmaUpperLower: All the matrix A
 *            = PlasmaUpper: Upper triangular part is set. The lower
 *            triangle is unchanged.
 *            = PlasmaLower: Lower triangular part is set. The upper
 *            triangle is unchange.
 *
 * @param[in] M
 *          The number of rows of the matrix A. M >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrix A. N >= 0.
 *
 * @param[in] alpha
 *          All the offdiagonal array elements are set to alpha.
 *
 * @param[in] beta
 *          All the diagonal array elements are set to beta.
 *
 * @param[in,out] A
 *          On entry, the m by n matrix A.
 *          On exit, A(i,j) = ALPHA, 1 <= i <= m, 1 <= j <= n, i.ne.j;
 *                   A(i,i) = BETA,  1 <= i <= min(m,n)
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,M).
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *
 *******************************************************************************
 *
 * @sa PLASMA_slaset_Tile
 * @sa PLASMA_slaset_Tile_Async
 * @sa PLASMA_claset
 * @sa PLASMA_slaset
 * @sa PLASMA_slaset
 *
 ******************************************************************************/
int PLASMA_slaset(PLASMA_enum uplo, int M, int N,
                  float alpha, float beta,
                  float *A, int LDA)
{
    int NB;
    int status;
    plasma_context_t *plasma;
    PLASMA_sequence *sequence = NULL;
    PLASMA_request request = PLASMA_REQUEST_INITIALIZER;
    PLASMA_desc descA;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_slaset", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    /* Check input arguments */
    if ( (uplo != PlasmaUpperLower) &&
         (uplo != PlasmaUpper) &&
         (uplo != PlasmaLower) ) {
        plasma_error("PLASMA_slaset", "illegal value of uplo");
        return -1;
    }
    if (M < 0) {
        plasma_error("PLASMA_slaset", "illegal value of M");
        return -2;
    }
    if (N < 0) {
        plasma_error("PLASMA_slaset", "illegal value of N");
        return -3;
    }
    if (LDA < max(1, M)) {
        plasma_error("PLASMA_slaset", "illegal value of LDA");
        return -5;
    }

    /* Quick return */
    if (min(N, M) == 0)
      return (float)0.0;

    /* Tune NB depending on M, N & NRHS; Set NBNB */
    status = plasma_tune(PLASMA_FUNC_SGEMM, M, N, 0);
    if (status != PLASMA_SUCCESS) {
        plasma_error("PLASMA_slaset", "plasma_tune() failed");
        return status;
    }

    /* Set NT */
    NB   = PLASMA_NB;

    plasma_sequence_create(plasma, &sequence);

    if ( PLASMA_TRANSLATION == PLASMA_OUTOFPLACE ) {
        plasma_sooplap2tile( descA, A, NB, NB, LDA, N, 0, 0, M, N, sequence, &request,
                             plasma_desc_mat_free(&(descA)) );
    } else {
        plasma_siplap2tile(  descA, A, NB, NB, LDA, N, 0, 0, M, N,
                            sequence, &request);
    }

    /* Call the tile interface */
    PLASMA_slaset_Tile_Async(uplo, alpha, beta, &descA, sequence, &request);

    if ( PLASMA_TRANSLATION == PLASMA_OUTOFPLACE ) {
        plasma_sooptile2lap( descA, A, NB, NB, LDA, N,  sequence, &request);
        plasma_dynamic_sync();
        plasma_desc_mat_free(&descA);
    } else {
        plasma_siptile2lap( descA, A, NB, NB, LDA, N,  sequence, &request);
        plasma_dynamic_sync();
    }

    plasma_sequence_destroy(plasma, sequence);
    return PLASMA_SUCCESS;
}

/***************************************************************************//**
 *
 * @ingroup float_Tile
 *
 *  PLASMA_slaset_Tile - Tile equivalent of PLASMA_slaset().
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          Specifies the part of the matrix A to be set.
 *            = PlasmaUpperLower: All the matrix A
 *            = PlasmaUpper: Upper triangular part
 *            = PlasmaLower: Lower triangular part
 *
 * @param[in] alpha
 *          All the offdiagonal array elements are set to alpha.
 *
 * @param[in] beta
 *          All the diagonal array elements are set to beta.
 *
 * @param[in,out] A
 *          On entry, the m by n matrix A.
 *          On exit, A(i,j) = ALPHA, 1 <= i <= m, 1 <= j <= n, i.ne.j;
 *                   A(i,i) = BETA,  1 <= i <= min(m,n)
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa PLASMA_slaset
 * @sa PLASMA_slaset_Tile_Async
 * @sa PLASMA_claset_Tile
 * @sa PLASMA_slaset_Tile
 * @sa PLASMA_slaset_Tile
 *
 ******************************************************************************/
int PLASMA_slaset_Tile(PLASMA_enum uplo,
                       float alpha, float beta,
                       PLASMA_desc *A)
{
    plasma_context_t *plasma;
    PLASMA_sequence *sequence = NULL;
    PLASMA_request request = PLASMA_REQUEST_INITIALIZER;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_slaset_Tile", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    plasma_sequence_create(plasma, &sequence);
    PLASMA_slaset_Tile_Async(uplo, alpha, beta, A, sequence, &request);
    plasma_dynamic_sync();
    plasma_sequence_destroy(plasma, sequence);
    return PLASMA_SUCCESS;
}

/***************************************************************************//**
 *
 * @ingroup float_Tile_Async
 *
 *  PLASMA_slaset_Tile_Async - Non-blocking equivalent of PLASMA_slaset_Tile().
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
 * @sa PLASMA_slaset
 * @sa PLASMA_slaset_Tile
 * @sa PLASMA_claset_Tile_Async
 * @sa PLASMA_slaset_Tile_Async
 * @sa PLASMA_slaset_Tile_Async
 *
 ******************************************************************************/
int PLASMA_slaset_Tile_Async(PLASMA_enum uplo,
                             float alpha, float beta,
                             PLASMA_desc *A,
                             PLASMA_sequence *sequence, PLASMA_request *request)
{
    PLASMA_desc descA;
    plasma_context_t *plasma;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_slaset_Tile_Async", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        plasma_fatal_error("PLASMA_slaset_Tile_Async", "NULL sequence");
        return PLASMA_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        plasma_fatal_error("PLASMA_slaset_Tile_Async", "NULL request");
        return PLASMA_ERR_UNALLOCATED;
    }
    /* Check sequence status */
    if (sequence->status == PLASMA_SUCCESS)
        request->status = PLASMA_SUCCESS;
    else
        return plasma_request_fail(sequence, request, PLASMA_ERR_SEQUENCE_FLUSHED);

    /* Check descriptors for correctness */
    if (plasma_desc_check(A) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_slaset_Tile_Async", "invalid descriptor");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    } else {
        descA = *A;
    }
    /* Check input arguments */
    if (descA.nb != descA.mb) {
        plasma_error("PLASMA_slaset_Tile_Async", "only square tiles supported");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }
    /* Check input arguments */
    if ( (uplo != PlasmaUpperLower) &&
         (uplo != PlasmaUpper) &&
         (uplo != PlasmaLower) ) {
        plasma_error("PLASMA_slaset_Tile_Async", "illegal value of uplo");
        return -1;
    }
    /* Quick return */
    if (min(descA.m, descA.n) == 0) {
        return PLASMA_SUCCESS;
    }

    plasma_dynamic_call_6(plasma_pslaset,
        PLASMA_enum, uplo,
        float, alpha,
        float, beta,
        PLASMA_desc, descA,
        PLASMA_sequence*, sequence,
        PLASMA_request*, request);

    return PLASMA_SUCCESS;
}

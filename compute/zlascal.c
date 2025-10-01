/**
 *
 * @file zlascal.c
 *
 *  PLASMA computational routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @precisions normal z -> s d c
 *
 **/
#include "common.h"

/***************************************************************************//**
 *
 * @ingroup PLASMA_Complex64_t
 *
 *  PLASMA_zlascal - Scales a matrix by the scalar alpha as in
 *  ScaLAPACK pzlascal().
 *
 *    \f[ A = \alpha A \f],
 *
 *  alpha is a scalar, and A a general, upper or lower trapezoidal matrix.
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          Specifies the shape of A:
 *          = PlasmaUpperLower: A is a general matrix.
 *          = PlasmaUpper: A is an upper trapezoidal matrix.
 *          = PlasmaLower: A is a lower trapezoidal matrix.
 *
 * @param[in] M
 *          M specifies the number of rows of the matrix A. M >= 0.
 *
 * @param[in] N
 *          N specifies the number of columns of the matrix A. N >= 0.
 *
 * @param[in] alpha
 *          alpha specifies the scalar alpha
 *
 * @param[in,out] A
 *          A is a LDA-by-N matrix.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,M).
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa PLASMA_zlascal_Tile
 * @sa PLASMA_clascal
 * @sa PLASMA_dlascal
 * @sa PLASMA_slascal
 *
 ******************************************************************************/
int PLASMA_zlascal(PLASMA_enum uplo, int M, int N,
                   PLASMA_Complex64_t alpha, PLASMA_Complex64_t *A, int LDA)
{
    int NB;
    int status;
    PLASMA_desc descA;
    plasma_context_t *plasma;
    PLASMA_sequence *sequence = NULL;
    PLASMA_request request = PLASMA_REQUEST_INITIALIZER;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_zlascal", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }

    /* Check input arguments */
    if (uplo != PlasmaUpper && uplo != PlasmaLower && uplo != PlasmaUpperLower) {
        plasma_error("PLASMA_zlascal", "illegal value of uplo");
        return -1;
    }
    if (M < 0) {
        plasma_error("PLASMA_zlascal", "illegal value of M");
        return -2;
    }
    if (N < 0) {
        plasma_error("PLASMA_zlascal", "illegal value of N");
        return -3;
    }
    if (LDA < max(1, M)) {
        plasma_error("PLASMA_zlascal", "illegal value of LDA");
        return -6;
    }

    /* Quick return */
    if (M == 0 || N == 0 ||
        (alpha == (PLASMA_Complex64_t)1.0))
        return PLASMA_SUCCESS;

    /* Tune NB depending on M, N & NRHS; Set NBNBSIZE */
    status = plasma_tune(PLASMA_FUNC_ZGEMM, M, N, 0);
    if (status != PLASMA_SUCCESS) {
        plasma_error("PLASMA_zlascal", "plasma_tune() failed");
        return status;
    }

    /* Set MT & NT & KT */
    NB = PLASMA_NB;

    plasma_sequence_create(plasma, &sequence);

    if ( PLASMA_TRANSLATION == PLASMA_OUTOFPLACE ) {
        plasma_zooplap2tile( descA, A, NB, NB, LDA, N, 0, 0, M, N, sequence, &request,
                             plasma_desc_mat_free(&(descA)) );
    } else {
        plasma_ziplap2tile( descA, A, NB, NB, LDA, N , 0, 0, M, N,
                            sequence, &request);
    }

    /* Call the tile interface */
    PLASMA_zlascal_Tile_Async(
        uplo, alpha, &descA, sequence, &request);

    if ( PLASMA_TRANSLATION == PLASMA_OUTOFPLACE ) {
        plasma_zooptile2lap( descA, A, NB, NB, LDA, N,  sequence, &request);
        plasma_dynamic_sync();
        plasma_desc_mat_free(&descA);
    } else {
        plasma_ziptile2lap( descA, A, NB, NB, LDA, N, sequence, &request);
        plasma_dynamic_sync();
    }

    status = sequence->status;
    plasma_sequence_destroy(plasma, sequence);
    return status;
}

/***************************************************************************//**
 *
 * @ingroup PLASMA_Complex64_t_Tile
 *
 *  PLASMA_zlascal_Tile - Scales a matrix by the scalar alpha as in
 *  ScaLAPACK pzlascal().
 *
 *    \f[ A = \alpha A \f],
 *
 *  alpha is a scalar, and A a general, upper or lower trapezoidal matrix.
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          Specifies the shape of A:
 *          = PlasmaUpperLower: A is a general matrix.
 *          = PlasmaUpper: A is an upper trapezoidal matrix.
 *          = PlasmaLower: A is a lower trapezoidal matrix.
 *
 * @param[in] alpha
 *          alpha specifies the scalar alpha
 *
 * @param[in] A
 *          A is a LDA-by-N matrix.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa PLASMA_zlascal
 * @sa PLASMA_zlascal_Tile_Async
 * @sa PLASMA_clascal_Tile
 * @sa PLASMA_dlascal_Tile
 * @sa PLASMA_slascal_Tile
 *
 ******************************************************************************/
int PLASMA_zlascal_Tile(PLASMA_enum uplo,
                        PLASMA_Complex64_t alpha, PLASMA_desc *A)
{
    plasma_context_t *plasma;
    PLASMA_sequence *sequence = NULL;
    PLASMA_request request = PLASMA_REQUEST_INITIALIZER;
    int status;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_zlascal_Tile", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    plasma_sequence_create(plasma, &sequence);
    PLASMA_zlascal_Tile_Async(uplo, alpha, A, sequence, &request);
    plasma_dynamic_sync();
    status = sequence->status;
    plasma_sequence_destroy(plasma, sequence);
    return status;
}

/***************************************************************************//**
 *
 * @ingroup PLASMA_Complex64_t_Tile_Async
 *
 *  PLASMA_zlascal_Tile_Async - Scales a matrix by the scalar alpha as in
 *  ScaLAPACK pzlascal().
 *  Non-blocking equivalent of PLASMA_zlascal_Tile().
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
 * @sa PLASMA_zlascal
 * @sa PLASMA_zlascal_Tile
 * @sa PLASMA_clascal_Tile_Async
 * @sa PLASMA_dlascal_Tile_Async
 * @sa PLASMA_slascal_Tile_Async
 *
 ******************************************************************************/
int PLASMA_zlascal_Tile_Async(PLASMA_enum uplo,
                              PLASMA_Complex64_t alpha, PLASMA_desc *A,
                              PLASMA_sequence *sequence, PLASMA_request *request)
{
    plasma_context_t *plasma;
    PLASMA_desc descA;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_zlascal_Tile_Async", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        plasma_fatal_error("PLASMA_zlascal_Tile_Async", "NULL sequence");
        return PLASMA_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        plasma_fatal_error("PLASMA_zlascal_Tile_Async", "NULL request");
        return PLASMA_ERR_UNALLOCATED;
    }
    /* Check sequence status */
    if (sequence->status == PLASMA_SUCCESS)
        request->status = PLASMA_SUCCESS;
    else
        return plasma_request_fail(sequence, request, PLASMA_ERR_SEQUENCE_FLUSHED);

    /* Check descriptors for correctness */
    if (plasma_desc_check(A) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_zlascal_Tile_Async", "invalid first descriptor");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    } else {
        descA = *A;
    }
    /* Check input arguments */
    if (uplo != PlasmaUpper && uplo != PlasmaLower && uplo != PlasmaUpperLower) {
        plasma_error("PLASMA_zlascal", "illegal value of uplo");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }

    if ( (descA.i%descA.mb != 0) || (descA.j%descA.nb != 0) ) {
        plasma_error("PLASMA_zlascal", "start indexes have to be multiple of tile size");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }

    /* Quick return */
    if ( (descA.m == 0) || (descA.n == 0) ||
         (alpha == (PLASMA_Complex64_t)1.0) )
        return PLASMA_SUCCESS;

    plasma_dynamic_call_5(plasma_pzlascal,
        PLASMA_enum, uplo,
        PLASMA_Complex64_t, alpha,
        PLASMA_desc, descA,
        PLASMA_sequence*, sequence,
        PLASMA_request*, request);

    return PLASMA_SUCCESS;
}

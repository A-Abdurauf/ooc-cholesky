/**
 *
 * @file ssymm.c
 *
 *  PLASMA computational routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @generated s Fri Apr  1 11:02:54 2016
 *
 **/
#include "common.h"

/***************************************************************************//**
 *
 * @ingroup float
 *
 *  PLASMA_ssymm - Performs one of the matrix-matrix operations
 *
 *     \f[ C = \alpha \times A \times B + \beta \times C \f]
 *
 *  or
 *
 *     \f[ C = \alpha \times B \times A + \beta \times C \f]
 *
 *  where alpha and beta are scalars, A is an symmetric matrix and  B and
 *  C are m by n matrices.
 *
 *******************************************************************************
 *
 * @param[in] side
 *          Specifies whether the symmetric matrix A appears on the
 *          left or right in the operation as follows:
 *          = PlasmaLeft:      \f[ C = \alpha \times A \times B + \beta \times C \f]
 *          = PlasmaRight:     \f[ C = \alpha \times B \times A + \beta \times C \f]
 *
 * @param[in] uplo
 *          Specifies whether the upper or lower triangular part of
 *          the symmetric matrix A is to be referenced as follows:
 *          = PlasmaLower:     Only the lower triangular part of the
 *                             symmetric matrix A is to be referenced.
 *          = PlasmaUpper:     Only the upper triangular part of the
 *                             symmetric matrix A is to be referenced.
 *
 * @param[in] M
 *          Specifies the number of rows of the matrix C. M >= 0.
 *
 * @param[in] N
 *          Specifies the number of columns of the matrix C. N >= 0.
 *
 * @param[in] alpha
 *          Specifies the scalar alpha.
 *
 * @param[in] A
 *          A is a LDA-by-ka matrix, where ka is M when side = PlasmaLeft,
 *          and is N otherwise. Only the uplo triangular part is referenced.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,ka).
 *
 * @param[in] B
 *          B is a LDB-by-N matrix, where the leading M-by-N part of
 *          the array B must contain the matrix B.
 *
 * @param[in] LDB
 *          The leading dimension of the array B. LDB >= max(1,M).
 *
 * @param[in] beta
 *          Specifies the scalar beta.
 *
 * @param[in,out] C
 *          C is a LDC-by-N matrix.
 *          On exit, the array is overwritten by the M by N updated matrix.
 *
 * @param[in] LDC
 *          The leading dimension of the array C. LDC >= max(1,M).
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa PLASMA_ssymm_Tile
 * @sa PLASMA_csymm
 * @sa PLASMA_dsymm
 * @sa PLASMA_ssymm
 *
 ******************************************************************************/
int PLASMA_ssymm(PLASMA_enum side, PLASMA_enum uplo, int M, int N,
                 float alpha, float *A, int LDA,
                                           float *B, int LDB,
                 float beta,  float *C, int LDC)
{
    int NB;
    int Am;
    int status;
    PLASMA_desc descA, descB, descC;
    plasma_context_t *plasma;
    PLASMA_sequence *sequence = NULL;
    PLASMA_request request = PLASMA_REQUEST_INITIALIZER;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_ssymm", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }

    /* Check input arguments */
    if ( (side != PlasmaLeft) && (side != PlasmaRight) ){
        plasma_error("PLASMA_ssymm", "illegal value of side");
        return -1;
    }
    if ((uplo != PlasmaLower) && (uplo != PlasmaUpper)) {
        plasma_error("PLASMA_ssymm", "illegal value of uplo");
        return -2;
    }
    Am = ( side == PlasmaLeft ) ? M : N;
    if (M < 0) {
        plasma_error("PLASMA_ssymm", "illegal value of M");
        return -3;
    }
    if (N < 0) {
        plasma_error("PLASMA_ssymm", "illegal value of N");
        return -4;
    }
    if (LDA < max(1, Am)) {
        plasma_error("PLASMA_ssymm", "illegal value of LDA");
        return -7;
    }
    if (LDB < max(1, M)) {
        plasma_error("PLASMA_ssymm", "illegal value of LDB");
        return -9;
    }
    if (LDC < max(1, M)) {
        plasma_error("PLASMA_ssymm", "illegal value of LDC");
        return -12;
    }

    /* Quick return */
    if (M == 0 || N == 0 ||
        ((alpha == (float)0.0) && beta == (float)1.0))
        return PLASMA_SUCCESS;

    /* Tune NB depending on M, N & NRHS; Set NBNB */
    status = plasma_tune(PLASMA_FUNC_SSYMM, M, N, 0);
    if (status != PLASMA_SUCCESS) {
        plasma_error("PLASMA_ssymm", "plasma_tune() failed");
        return status;
    }

    /* Set MT & NT & KT */
    NB = PLASMA_NB;

    plasma_sequence_create(plasma, &sequence);

    if ( PLASMA_TRANSLATION == PLASMA_OUTOFPLACE ) {
        plasma_sooplap2tile( descA, A, NB, NB, LDA, Am, 0, 0, Am, Am, sequence, &request,
                             plasma_desc_mat_free(&(descA)) );
        plasma_sooplap2tile( descB, B, NB, NB, LDB, N,  0, 0, M,  N, sequence, &request,
                             plasma_desc_mat_free(&(descA)); plasma_desc_mat_free(&(descB)));
        plasma_sooplap2tile( descC, C, NB, NB, LDC, N,  0, 0, M,  N, sequence, &request,
                             plasma_desc_mat_free(&(descA)); plasma_desc_mat_free(&(descB)); plasma_desc_mat_free(&(descC)));
    } else {
        plasma_siplap2tile( descA, A, NB, NB, LDA, Am, 0, 0, Am, Am,
                            sequence, &request);
        plasma_siplap2tile( descB, B, NB, NB, LDB, N,  0, 0, M,  N,
                            sequence, &request);
        plasma_siplap2tile( descC, C, NB, NB, LDC, N,  0, 0, M,  N,
                            sequence, &request);
    }

    /* Call the tile interface */
    PLASMA_ssymm_Tile_Async(
        side, uplo, alpha, &descA, &descB, beta, &descC, sequence, &request);

    if ( PLASMA_TRANSLATION == PLASMA_OUTOFPLACE ) {
        plasma_sooptile2lap( descC, C, NB, NB, LDC, N,  sequence, &request);
        plasma_dynamic_sync();
        plasma_desc_mat_free(&descA);
        plasma_desc_mat_free(&descB);
        plasma_desc_mat_free(&descC);
    } else {
        plasma_siptile2lap( descA, A, NB, NB, LDA, Am,  sequence, &request);
        plasma_siptile2lap( descB, B, NB, NB, LDB, N,  sequence, &request);
        plasma_siptile2lap( descC, C, NB, NB, LDC, N,  sequence, &request);
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
 *  PLASMA_ssymm_Tile - Performs symmetric matrix multiplication.
 *  Tile equivalent of PLASMA_ssymm().
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] side
 *          Specifies whether the symmetric matrix A appears on the
 *          left or right in the operation as follows:
 *          = PlasmaLeft:      \f[ C = \alpha \times A \times B + \beta \times C \f]
 *          = PlasmaRight:     \f[ C = \alpha \times B \times A + \beta \times C \f]
 *
 * @param[in] uplo
 *          Specifies whether the upper or lower triangular part of
 *          the symmetric matrix A is to be referenced as follows:
 *          = PlasmaLower:     Only the lower triangular part of the
 *                             symmetric matrix A is to be referenced.
 *          = PlasmaUpper:     Only the upper triangular part of the
 *                             symmetric matrix A is to be referenced.
 *
 * @param[in] alpha
 *          Specifies the scalar alpha.
 *
 * @param[in] A
 *          A is a LDA-by-ka matrix, where ka is M when side = PlasmaLeft,
 *          and is N otherwise. Only the uplo triangular part is referenced.
 *
 * @param[in] B
 *          B is a LDB-by-N matrix, where the leading M-by-N part of
 *          the array B must contain the matrix B.
 *
 * @param[in] beta
 *          Specifies the scalar beta.
 *
 * @param[in,out] C
 *          C is a LDC-by-N matrix.
 *          On exit, the array is overwritten by the M by N updated matrix.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa PLASMA_ssymm
 * @sa PLASMA_ssymm_Tile_Async
 * @sa PLASMA_csymm_Tile
 * @sa PLASMA_dsymm_Tile
 * @sa PLASMA_ssymm_Tile
 *
 ******************************************************************************/
int PLASMA_ssymm_Tile(PLASMA_enum side, PLASMA_enum uplo,
                      float alpha, PLASMA_desc *A, PLASMA_desc *B,
                      float beta,  PLASMA_desc *C)
{
    plasma_context_t *plasma;
    PLASMA_sequence *sequence = NULL;
    PLASMA_request request = PLASMA_REQUEST_INITIALIZER;
    int status;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_ssymm_Tile", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    plasma_sequence_create(plasma, &sequence);
    PLASMA_ssymm_Tile_Async(side, uplo, alpha, A, B, beta, C, sequence, &request);
    plasma_dynamic_sync();
    status = sequence->status;
    plasma_sequence_destroy(plasma, sequence);
    return status;
}

/***************************************************************************//**
 *
 * @ingroup float_Tile_Async
 *
 *  PLASMA_ssymm_Tile_Async - Performs symmetric matrix multiplication.
 *  Non-blocking equivalent of PLASMA_ssymm_Tile().
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
 * @sa PLASMA_ssymm
 * @sa PLASMA_ssymm_Tile
 * @sa PLASMA_csymm_Tile_Async
 * @sa PLASMA_dsymm_Tile_Async
 * @sa PLASMA_ssymm_Tile_Async
 *
 ******************************************************************************/
int PLASMA_ssymm_Tile_Async(PLASMA_enum side, PLASMA_enum uplo,
                            float alpha, PLASMA_desc *A, PLASMA_desc *B,
                            float beta,  PLASMA_desc *C,
                            PLASMA_sequence *sequence, PLASMA_request *request)
{
    plasma_context_t *plasma;
    PLASMA_desc descA;
    PLASMA_desc descB;
    PLASMA_desc descC;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_ssymm_Tile_Async", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        plasma_fatal_error("PLASMA_ssymm_Tile_Async", "NULL sequence");
        return PLASMA_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        plasma_fatal_error("PLASMA_ssymm_Tile_Async", "NULL request");
        return PLASMA_ERR_UNALLOCATED;
    }
    /* Check sequence status */
    if (sequence->status == PLASMA_SUCCESS)
        request->status = PLASMA_SUCCESS;
    else
        return plasma_request_fail(sequence, request, PLASMA_ERR_SEQUENCE_FLUSHED);

    /* Check descriptors for correctness */
    if (plasma_desc_check(A) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_ssymm_Tile_Async", "invalid first descriptor");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    } else {
        descA = *A;
    }
    if (plasma_desc_check(B) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_ssymm_Tile_Async", "invalid second descriptor");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    } else {
        descB = *B;
    }
    if (plasma_desc_check(C) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_ssymm_Tile_Async", "invalid third descriptor");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    } else {
        descC = *C;
    }
    /* Check input arguments */
    if ( (side != PlasmaLeft) && (side != PlasmaRight) ){
        plasma_error("PLASMA_ssymm_Tile_Async", "illegal value of side");
        return plasma_request_fail(sequence, request, -1);
    }
    if ((uplo != PlasmaLower) && (uplo != PlasmaUpper)) {
        plasma_error("PLASMA_ssymm_Tile_Async", "illegal value of uplo");
        return plasma_request_fail(sequence, request, -2);
    }

    /* Check matrices sizes */
    if ( (descB.m != descC.m) || (descB.n != descC.n) ) {
        plasma_error("PLASMA_ssymm_Tile_Async", "B and C must have the same size");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }
    if ( (descA.m != descA.n) ||
         ( (side == PlasmaLeft)  && (descA.m != descB.m ) ) ||
         ( (side == PlasmaRight) && (descA.m != descB.n ) ) ) {
        plasma_error("PLASMA_ssymm_Tile_Async", "Matrix A must be square of size M or N regarding side.");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }

    /* Check tiles sizes */
    if ( (descB.mb != descC.mb) || (descB.nb != descC.nb) ) {
        plasma_error("PLASMA_ssymm_Tile_Async", "B and C must have the same tile sizes");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }
    if ( (descA.mb != descA.nb) ||
         ( (side == PlasmaLeft)  && (descA.mb != descB.mb ) ) ||
         ( (side == PlasmaRight) && (descA.mb != descB.nb ) ) ) {
        plasma_error("PLASMA_ssymm_Tile_Async", "Matrix A must be square with square tiles wich fits the reagding tile size of B and C");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }

    /* Check submatrix starting point */
    /* Check submatrix starting point */
    if ( (descB.i%descB.mb != descC.i%descC.mb) ||
         (descB.j%descB.nb != descC.j%descC.nb) ) {
        plasma_error("PLASMA_ssymm_Tile_Async", "B and C submatrices doesn't match");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }
    if ( (descA.i%descA.mb != descA.j%descA.nb) ||
         ( (side == PlasmaLeft)  && (descA.i%descA.mb != descB.i%descB.mb ) ) ||
         ( (side == PlasmaRight) && (descA.i%descA.mb != descB.j%descB.nb ) ) ) {
        plasma_error("PLASMA_ssymm_Tile_Async", "Submatrix A must start on diagnonal and match submatrices B and C.");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }

    /* Quick return */
    if (descC.m == 0 || descC.n == 0 ||
        ( (alpha == (float)0.0) && (beta == (float)1.0) ))
        return PLASMA_SUCCESS;

    plasma_parallel_call_9(plasma_pssymm,
        PLASMA_enum, side,
        PLASMA_enum, uplo,
        float, alpha,
        PLASMA_desc, descA,
        PLASMA_desc, descB,
        float, beta,
        PLASMA_desc, descC,
        PLASMA_sequence*, sequence,
        PLASMA_request*, request);

    return PLASMA_SUCCESS;
}

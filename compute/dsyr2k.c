/**
 *
 * @file dsyr2k.c
 *
 *  PLASMA computational routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @generated d Fri Apr  1 11:02:54 2016
 *
 **/
#include "common.h"

/***************************************************************************//**
 *
 * @ingroup double
 *
 *  PLASMA_dsyr2k - Performs one of the symmetric rank 2k operations
 *
 *    \f[ C = \alpha [ op( A ) \times g( op( B )' )] + \alpha [ op( B ) \times g( op( A )' )] + \beta C \f],
 *    or
 *    \f[ C = \alpha [ g( op( A )' ) \times op( B ) ] + \alpha [ g( op( B )' ) \times op( A ) ] + \beta C \f],
 *
 *  where op( X ) is one of
 *
 *    op( X ) = X  or op( X ) = g( X' )
 *
 *  where alpha and beta are real scalars, C is an n-by-n symmetric
 *  matrix and A and B are an n-by-k matrices the first case and k-by-n
 *  matrices in the second case.
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          = PlasmaUpper: Upper triangle of C is stored;
 *          = PlasmaLower: Lower triangle of C is stored.
 *
 * @param[in] trans
 *          Specifies whether the matrix A is transposed or ugate transposed:
 *          = PlasmaNoTrans: \f[ C = \alpha [ op( A ) \times g( op( B )' )] + \alpha [ op( B ) \times g( op( A )' )] + \beta C \f]
 *          = PlasmaTrans: \f[ C = \alpha [ g( op( A )' ) \times op( B ) ] + \alpha [ g( op( B )' ) \times op( A ) ] + \beta C \f]
 *
 * @param[in] N
 *          N specifies the order of the matrix C. N must be at least zero.
 *
 * @param[in] K
 *          K specifies the number of columns of the A and B matrices with trans = PlasmaNoTrans.
 *          K specifies the number of rows of the A and B matrices with trans = PlasmaTrans.
 *
 * @param[in] alpha
 *          alpha specifies the scalar alpha.
 *
 * @param[in] A
 *          A is a LDA-by-ka matrix, where ka is K when trans = PlasmaNoTrans,
 *          and is N otherwise.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA must be at least
 *          max( 1, N ), otherwise LDA must be at least max( 1, K ).
 *
 * @param[in] B
 *          B is a LDB-by-kb matrix, where kb is K when trans = PlasmaNoTrans,
 *          and is N otherwise.
 *
 * @param[in] LDB
 *          The leading dimension of the array B. LDB must be at least
 *          max( 1, N ), otherwise LDB must be at least max( 1, K ).
 *
 * @param[in] beta
 *          beta specifies the scalar beta.
 *
 * @param[in,out] C
 *          C is a LDC-by-N matrix.
 *          On exit, the array uplo part of the matrix is overwritten
 *          by the uplo part of the updated matrix.
 *
 * @param[in] LDC
 *          The leading dimension of the array C. LDC >= max( 1, N ).
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa PLASMA_dsyr2k_Tile
 * @sa PLASMA_dsyr2k_Tile_Async
 * @sa PLASMA_csyr2k
 * @sa PLASMA_dsyr2k
 * @sa PLASMA_ssyr2k
 *
 ******************************************************************************/
int PLASMA_dsyr2k(PLASMA_enum uplo, PLASMA_enum trans, int N, int K,
                  double alpha, double *A, int LDA, double *B, int LDB,
                  double beta,  double *C, int LDC)
{
    int NB;
    int Am, An;
    int status;
    PLASMA_desc descA, descB, descC;
    plasma_context_t *plasma;
    PLASMA_sequence *sequence = NULL;
    PLASMA_request request = PLASMA_REQUEST_INITIALIZER;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_dsyr2k", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }

    /* Check input arguments */
    if ((uplo != PlasmaUpper) && (uplo != PlasmaLower)) {
        plasma_error("PLASMA_dsyr2k", "illegal value of uplo");
        return -1;
    }
    if ((trans != PlasmaNoTrans) && (trans != PlasmaTrans)) {
        plasma_error("PLASMA_dsyr2k", "illegal value of trans");
        return -2;
    }
    if ( trans == PlasmaNoTrans ) {
        Am = N; An = K;
    } else {
        Am = K; An = N;
    }
    if (N < 0) {
        plasma_error("PLASMA_dsyr2k", "illegal value of N");
        return -3;
    }
    if (K < 0) {
        plasma_error("PLASMA_dsyr2k", "illegal value of K");
        return -4;
    }
    if (LDA < max(1, Am)) {
        plasma_error("PLASMA_dsyr2k", "illegal value of LDA");
        return -7;
    }
    if (LDB < max(1, Am)) {
        plasma_error("PLASMA_dsyr2k", "illegal value of LDB");
        return -9;
    }
    if (LDC < max(1, N)) {
        plasma_error("PLASMA_dsyr2k", "illegal value of LDC");
        return -12;
    }

    /* Quick return */
    if (N == 0 ||
        ((alpha == (double)0.0 || K == 0.0) && beta == (double)1.0))
        return PLASMA_SUCCESS;

    /* Tune NB depending on M, N & NRHS; Set NBNB */
    status = plasma_tune(PLASMA_FUNC_DSYRK, N, K, 0);
    if (status != PLASMA_SUCCESS) {
        plasma_error("PLASMA_dsyr2k", "plasma_tune() failed");
        return status;
    }

    /* Set MT & NT & KT */
    NB = PLASMA_NB;

    plasma_sequence_create(plasma, &sequence);

    if ( PLASMA_TRANSLATION == PLASMA_OUTOFPLACE ) {
        plasma_dooplap2tile( descA, A, NB, NB, LDA, An, 0, 0, Am, An, sequence, &request,
                             plasma_desc_mat_free(&(descA)) );
        plasma_dooplap2tile( descB, B, NB, NB, LDB, An, 0, 0, Am, An, sequence, &request,
                             plasma_desc_mat_free(&(descA)); plasma_desc_mat_free(&(descB)));
        plasma_dooplap2tile( descC, C, NB, NB, LDC, N,  0, 0, N,  N,   sequence, &request,
                             plasma_desc_mat_free(&(descA)); plasma_desc_mat_free(&(descB)); plasma_desc_mat_free(&(descC)));
    } else {
        plasma_diplap2tile( descA, A, NB, NB, LDA, An, 0, 0, Am, An,
                            sequence, &request);
        plasma_diplap2tile( descB, B, NB, NB, LDB, An, 0, 0, Am, An,
                            sequence, &request);
        plasma_diplap2tile( descC, C, NB, NB, LDC, N,  0, 0, N,  N,
                            sequence, &request);
    }

    /* Call the tile interface */
    PLASMA_dsyr2k_Tile_Async(uplo, trans, alpha, &descA, &descB, beta, &descC, sequence, &request);

    if ( PLASMA_TRANSLATION == PLASMA_OUTOFPLACE ) {
        plasma_dooptile2lap( descC, C, NB, NB, LDC, N,  sequence, &request);
        plasma_dynamic_sync();
        plasma_desc_mat_free(&descA);
        plasma_desc_mat_free(&descB);
        plasma_desc_mat_free(&descC);
    } else {
        plasma_diptile2lap( descA, A, NB, NB, LDA, An,  sequence, &request);
        plasma_diptile2lap( descB, B, NB, NB, LDB, An,  sequence, &request);
        plasma_diptile2lap( descC, C, NB, NB, LDC, N,  sequence, &request);
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
 *  PLASMA_dsyr2k_Tile - Performs symmetric rank k update.
 *  Tile equivalent of PLASMA_dsyr2k().
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          = PlasmaUpper: Upper triangle of C is stored;
 *          = PlasmaLower: Lower triangle of C is stored.
 *
 * @param[in] trans
 *          Specifies whether the matrix A is transposed or ugate transposed:
 *          = PlasmaNoTrans:   A is not transposed;
 *          = PlasmaTrans: A is ugate transposed.
 *
 * @param[in] alpha
 *          alpha specifies the scalar alpha.
 *
 * @param[in] A
 *          A is a LDA-by-ka matrix, where ka is K when trans = PlasmaNoTrans,
 *          and is N otherwise.
 *
 * @param[in] B
 *          B is a LDB-by-kb matrix, where kb is K when trans = PlasmaNoTrans,
 *          and is N otherwise.
 *
 * @param[in] beta
 *          beta specifies the scalar beta
 *
 * @param[in,out] C
 *          C is a LDC-by-N matrix.
 *          On exit, the array uplo part of the matrix is overwritten
 *          by the uplo part of the updated matrix.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa PLASMA_dsyr2k
 * @sa PLASMA_dsyr2k_Tile_Async
 * @sa PLASMA_csyr2k_Tile
 * @sa PLASMA_dsyr2k_Tile
 * @sa PLASMA_ssyr2k_Tile
 *
 ******************************************************************************/
int PLASMA_dsyr2k_Tile(PLASMA_enum uplo, PLASMA_enum trans,
                       double alpha, PLASMA_desc *A, PLASMA_desc *B,
                       double beta,  PLASMA_desc *C)
{
    plasma_context_t *plasma;
    PLASMA_sequence *sequence = NULL;
    PLASMA_request request = PLASMA_REQUEST_INITIALIZER;
    int status;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_dsyr2k_Tile", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    plasma_sequence_create(plasma, &sequence);
    PLASMA_dsyr2k_Tile_Async(uplo, trans, alpha, A, B, beta, C, sequence, &request);
    plasma_dynamic_sync();
    status = sequence->status;
    plasma_sequence_destroy(plasma, sequence);
    return status;
}

/***************************************************************************//**
 *
 * @ingroup double_Tile_Async
 *
 *  PLASMA_dsyr2k_Tile_Async - Performs symmetric rank-k update.
 *  Non-blocking equivalent of PLASMA_dsyr2k_Tile().
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
 * @sa PLASMA_dsyr2k
 * @sa PLASMA_dsyr2k_Tile
 * @sa PLASMA_csyr2k_Tile_Async
 * @sa PLASMA_dsyr2k_Tile_Async
 * @sa PLASMA_ssyr2k_Tile_Async
 *
 ******************************************************************************/
int PLASMA_dsyr2k_Tile_Async(PLASMA_enum uplo, PLASMA_enum trans,
                             double alpha, PLASMA_desc *A, PLASMA_desc *B,
                             double beta,  PLASMA_desc *C,
                             PLASMA_sequence *sequence, PLASMA_request *request)
{
    plasma_context_t *plasma;
    PLASMA_desc descA;
    PLASMA_desc descB;
    PLASMA_desc descC;
    int N, K;
    int Am, An, Amb;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_dsyr2k_Tile_Async", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        plasma_fatal_error("PLASMA_dsyr2k_Tile_Async", "NULL sequence");
        return PLASMA_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        plasma_fatal_error("PLASMA_dsyr2k_Tile_Async", "NULL request");
        return PLASMA_ERR_UNALLOCATED;
    }
    /* Check sequence status */
    if (sequence->status == PLASMA_SUCCESS)
        request->status = PLASMA_SUCCESS;
    else
        return plasma_request_fail(sequence, request, PLASMA_ERR_SEQUENCE_FLUSHED);

    /* Check descriptors for correctness */
    if (plasma_desc_check(A) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_dsyr2k_Tile_Async", "invalid first descriptor");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    } else {
        descA = *A;
    }
    if (plasma_desc_check(B) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_dsyr2k_Tile_Async", "invalid second descriptor");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    } else {
        descB = *B;
    }
    if (plasma_desc_check(C) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_dsyr2k_Tile_Async", "invalid third descriptor");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    } else {
        descC = *C;
    }
    /* Check input arguments */
    if ((uplo != PlasmaUpper) && (uplo != PlasmaLower)) {
        plasma_error("PLASMA_dsyr2k", "illegal value of uplo");
        return plasma_request_fail(sequence, request, -1);
    }
    if ((trans != PlasmaNoTrans) && (trans != PlasmaTrans)) {
        plasma_error("PLASMA_dsyr2k", "illegal value of trans");
        return plasma_request_fail(sequence, request, -2);
    }

    if ( trans == PlasmaNoTrans ) {
        Am  = descA.m;
        An  = descA.n;
        Amb = descA.mb;
    } else {
        Am  = descA.n;
        An  = descA.m;
        Amb = descA.nb;
    }

    if (descC.mb != descC.nb) {
        plasma_error("PLASMA_dsyr2k_Tile_Async", "only square tiles for C are supported");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }
    if ( (descB.mb != descA.mb) || (descB.nb != descA.nb) || (Amb != descC.mb) ){
        plasma_error("PLASMA_dsyr2k_Tile_Async", "tile sizes have to match");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }
    if (descC.m != descC.n) {
        plasma_error("PLASMA_dsyr2k_Tile_Async", "only square matrix C is supported");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }
    if ( (descB.m != descA.m) || (descB.n != descA.n) || (Am != descC.m) ){
        plasma_error("PLASMA_dsyr2k_Tile_Async", "sizes of matrices have to match");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }

    N = descC.m;
    K = An;

    /* Quick return */
    if ( N == 0 ||
        ((alpha == (double)0.0 || K == 0) && beta == (double)1.0))
        return PLASMA_SUCCESS;

    plasma_parallel_call_9(plasma_pdsyr2k,
        PLASMA_enum, uplo,
        PLASMA_enum, trans,
        double, alpha,
        PLASMA_desc, descA,
        PLASMA_desc, descB,
        double, beta,
        PLASMA_desc, descC,
        PLASMA_sequence*, sequence,
        PLASMA_request*, request);

    return PLASMA_SUCCESS;
}

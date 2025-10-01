/**
 *
 * @file stile.c
 *
 *  PLASMA auxiliary routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Jakub Kurzak
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @generated s Fri Apr  1 11:03:01 2016
 *
 **/
#include "common.h"
#include "auxiliary.h"
#include "tile.h"

/**
 *******************************************************************************
 *
 * @ingroup float
 *
 *  PLASMA_sLapack_to_Tile - Conversion from LAPACK layout to tile layout.
 *  It automatically infers the layout translation mode base on the values of
 *  Af77 and A->mat pointers.
 *  If (A->mat == Af77 || A->mat == NULL || Af77 == NULL),
 *  PLASMA_TRANSLATION_MODE is set to PLASMA_INPLACE to do inplace layout
 *  translation of the matrix.
 *  If (A->mat != Af77 && A->mat != NULL && Af77 != NULL),
 *  PLASMA_TRANSLATION_MODE is set to PLASMA_OUTOFPLACE to do out of place
 *  layout translation of the matrix. Both A->mat and Af77 have to be allocated
 *  to store the full matrix.
 *
 *******************************************************************************
 *
 * @param[in] Af77
 *          LAPACK matrix.
 *
 * @param[in] LDA
 *          The leading dimension of the matrix Af77.
 *
 * @param[in,out] A
 *          Descriptor of the PLASMA matrix in tile layout.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa PLASMA_sLapack_to_Tile_Async
 * @sa PLASMA_sTile_to_Lapack
 * @sa PLASMA_cLapack_to_Tile
 * @sa PLASMA_dLapack_to_Tile
 * @sa PLASMA_sLapack_to_Tile
 *
 ******************************************************************************/
int PLASMA_sLapack_to_Tile(float *Af77, int LDA, PLASMA_desc *A)
{
    PLASMA_desc descA = *A;
    plasma_context_t *plasma;
    PLASMA_sequence *sequence = NULL;
    PLASMA_request request;
    int status;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_sLapack_to_Tile", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    /* Check descriptor for correctness */
    if (plasma_desc_check(&descA) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_sLapack_to_Tile", "invalid descriptor");
        return PLASMA_ERR_ILLEGAL_VALUE;
    }
    plasma_sequence_create(plasma, &sequence);
    PLASMA_sLapack_to_Tile_Async(Af77, LDA, A, sequence, &request);
    plasma_dynamic_sync();
    status = sequence->status;
    plasma_sequence_destroy(plasma, sequence);
    return status;
}

/**
 *******************************************************************************
 *
 * @ingroup float_Tile_Async
 *
 *  PLASMA_sLapack_to_Tile_Async - Conversion from LAPACK layout to tile layout.
 *  Non-blocking equivalent of PLASMA_sLapack_to_Tile().  May return before the
 *  computation is finished.
 *  Allows for pipelining of operations ar runtime.
 *
 *******************************************************************************
 *
 * @param[in] Af77
 *          LAPACK matrix.
 *
 * @param[in] LDA
 *          The leading dimension of the matrix Af77.
 *
 * @param[in,out] A
 *          Descriptor of the PLASMA matrix in tile layout.
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
 * @sa PLASMA_sTile_to_Lapack_Async
 * @sa PLASMA_sLapack_to_Tile
 * @sa PLASMA_cLapack_to_Tile_Async
 * @sa PLASMA_dLapack_to_Tile_Async
 * @sa PLASMA_sLapack_to_Tile_Async
 *
 ******************************************************************************/
int PLASMA_sLapack_to_Tile_Async(float *Af77, int LDA, PLASMA_desc *A,
                                 PLASMA_sequence *sequence, PLASMA_request *request)
{
    PLASMA_desc descA = *A;
    plasma_context_t *plasma;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_sLapack_to_Tile", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    /* Check descriptor for correctness */
    if (plasma_desc_check(&descA) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_sLapack_to_Tile", "invalid descriptor");
        return PLASMA_ERR_ILLEGAL_VALUE;
    }
    if ( A->lm != LDA ){
        plasma_error("PLASMA_sLapack_to_Tile_Async",
                     "The leading dimension of the output matrix must be equal to the full matrix A->lm");
        return PLASMA_ERR_ILLEGAL_VALUE;
    }

    if ((A->mat == Af77) || (A->mat == NULL) || (Af77 == NULL)) {
        if ( A->mat == NULL )
            A->mat = Af77;

        PLASMA_sgecfi_Async( A->lm, A->ln, A->mat,
                             PlasmaCM, A->lm, 1,
                             PlasmaCCRB, A->mb, A->nb,
                             sequence, request);
    } else {
        plasma_parallel_call_5(
            plasma_pslapack_to_tile,
            float*, Af77,
            int,                 LDA,
            PLASMA_desc,         descA,
            PLASMA_sequence*,    sequence,
            PLASMA_request*,     request);
    }
    return PLASMA_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @ingroup float
 *
 *  PLASMA_Tile_to_Lapack - Conversion from tile layout to LAPACK layout.
 *  It automatically infers the layout translation mode base on the values of
 *  Af77 and A->mat pointers.
 *  If (A->mat == Af77 || A->mat == NULL || Af77 == NULL),
 *  PLASMA_TRANSLATION_MODE is set to PLASMA_INPLACE to do inplace layout
 *  translation of the matrix.
 *  If (A->mat != Af77 && A->mat != NULL && Af77 != NULL),
 *  PLASMA_TRANSLATION_MODE is set to PLASMA_OUTOFPLACE to do out of place
 *  layout translation of the matrix. Both A->mat and Af77 have to be allocated
 *  to store the full matrix.
 *
 *******************************************************************************
 *
 * @param[in] A
 *          Descriptor of the PLASMA matrix in tile layout.
 *
 * @param[in,out] Af77
 *          LAPACK matrix.
 *
 * @param[in] LDA
 *          The leading dimension of the matrix Af77.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa PLASMA_sTile_to_Lapack_Async
 * @sa PLASMA_sLapack_to_Tile
 * @sa PLASMA_cTile_to_Lapack
 * @sa PLASMA_dTile_to_Lapack
 * @sa PLASMA_sTile_to_Lapack
 *
******************************************************************************/
int PLASMA_sTile_to_Lapack(PLASMA_desc *A, float *Af77, int LDA)
{
    PLASMA_desc descA = *A;
    plasma_context_t *plasma;
    PLASMA_sequence *sequence = NULL;
    PLASMA_request request;
    int status;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_sTile_to_Lapack", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    /* Check descriptor for correctness */
    if (plasma_desc_check(&descA) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_sTile_to_Lapack", "invalid descriptor");
        return PLASMA_ERR_ILLEGAL_VALUE;
    }
    plasma_sequence_create(plasma, &sequence);
    PLASMA_sTile_to_Lapack_Async(A, Af77, LDA, sequence, &request);
    plasma_dynamic_sync();
    status = sequence->status;
    plasma_sequence_destroy(plasma, sequence);
    return status;
}

/**
 *******************************************************************************
 *
 * @ingroup float_Tile_Async
 *
 *  PLASMA_sTile_to_Lapack_Async - Conversion from LAPACK layout to tile layout.
 *  Non-blocking equivalent of PLASMA_sTile_to_Lapack().  May return before the
 *  computation is finished.
 *  Allows for pipelining of operations ar runtime.
 *
 *******************************************************************************
 *
 * @param[in] A
 *          Descriptor of the PLASMA matrix in tile layout.
 *
 * @param[in,out] Af77
 *          LAPACK matrix.
 *
 * @param[in] LDA
 *          The leading dimension of the matrix Af77.
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
 * @sa PLASMA_sLapack_to_Tile_Async
 * @sa PLASMA_sTile_to_Lapack
 * @sa PLASMA_cTile_to_Lapack_Async
 * @sa PLASMA_dTile_to_Lapack_Async
 * @sa PLASMA_sTile_to_Lapack_Async
 *
 ******************************************************************************/
int PLASMA_sTile_to_Lapack_Async(PLASMA_desc *A, float *Af77, int LDA,
                                 PLASMA_sequence *sequence, PLASMA_request *request)
{
    PLASMA_desc descA = *A;
    plasma_context_t *plasma;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_sTile_to_Lapack", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    /* Check descriptor for correctness */
    if (plasma_desc_check(&descA) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_sTile_to_Lapack", "invalid descriptor");
        return PLASMA_ERR_ILLEGAL_VALUE;
    }
    if ( A->lm != LDA ){
        plasma_error("PLASMA_sTile_to_Lapack_Async",
                     "The leading dimension of the output matrix must be equal to the full matrix A->lm");
        return PLASMA_ERR_ILLEGAL_VALUE;
    }

    if ((A->mat == Af77) || (A->mat == NULL) || (Af77 == NULL)) {
        if ( A->mat == NULL )
            A->mat = Af77;

        PLASMA_sgecfi_Async( A->lm, A->ln, A->mat,
                             PlasmaCCRB, A->mb, A->nb,
                             PlasmaCM, A->lm, 1,
                             sequence, request);
    } else {
        plasma_static_call_5(
            plasma_pstile_to_lapack,
            PLASMA_desc,         descA,
            float*, Af77,
            int,                 LDA,
            PLASMA_sequence*,    sequence,
            PLASMA_request*,     request);
    }
    return PLASMA_SUCCESS;
}

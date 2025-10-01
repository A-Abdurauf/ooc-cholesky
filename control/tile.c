/**
 *
 * @file tile.c
 *
 *  PLASMA auxiliary routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Jakub Kurzak
 * @date 2010-11-15
 *
 **/
#include "common.h"
#include "auxiliary.h"
#include "tile.h"

/***************************************************************************//**
 *
 * @ingroup Auxiliary
 *
 *  PLASMA_Lapack_to_Tile - Conversion from LAPACK layout to tile layout.
 *
 *******************************************************************************
 *
 * @param[in] Af77
 *          LAPACK matrix.
 *
 * @param[in] LDA
 *          The leading dimension of the matrix Af77.
 *
 * @param[out] A
 *          Descriptor of the PLASMA matrix in tile layout.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *
 ******************************************************************************/
int PLASMA_Lapack_to_Tile(void *Af77, int LDA, PLASMA_desc *A)
{
    PLASMA_desc descA = *A;
    plasma_context_t *plasma;
    PLASMA_sequence *sequence = NULL;
    PLASMA_request request;
    int status;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_Lapack_to_Tile", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    /* Check descriptor for correctness */
    if (plasma_desc_check(&descA) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_Lapack_to_Tile", "invalid descriptor");
        return PLASMA_ERR_ILLEGAL_VALUE;
    }

    plasma_sequence_create(plasma, &sequence);
    switch( A->dtyp ) {
    case PlasmaRealFloat:
        PLASMA_sLapack_to_Tile_Async( Af77, LDA, A, sequence, &request );
        break;
    case PlasmaRealDouble:
        PLASMA_dLapack_to_Tile_Async( Af77, LDA, A, sequence, &request );
        break;
    case PlasmaComplexFloat:
        PLASMA_cLapack_to_Tile_Async( Af77, LDA, A, sequence, &request );
        break;
    case PlasmaComplexDouble:
        PLASMA_zLapack_to_Tile_Async( Af77, LDA, A, sequence, &request );
        break;
    default:
        plasma_error("PLASMA_Lapack_to_Tile", "Type unknown");
    }
    plasma_dynamic_sync();
    status = sequence->status;
    plasma_sequence_destroy(plasma, sequence);
    return status;
}

/***************************************************************************//**
 *
 * @ingroup Auxiliary
 *
 *  PLASMA_Tile_to_Lapack - Conversion from tile layout to LAPACK layout.
 *
 *******************************************************************************
 *
 * @param[out] A
 *          Descriptor of the PLASMA matrix in tile layout.
 *
 * @param[in] Af77
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
 ******************************************************************************/
int PLASMA_Tile_to_Lapack(PLASMA_desc *A, void *Af77, int LDA)
{
    PLASMA_desc descA = *A;
    plasma_context_t *plasma;
    PLASMA_sequence *sequence = NULL;
    PLASMA_request request;
    int status;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_Tile_to_Lapack", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    /* Check descriptor for correctness */
    if (plasma_desc_check(&descA) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_Tile_to_Lapack", "invalid descriptor");
        return PLASMA_ERR_ILLEGAL_VALUE;
    }

    plasma_sequence_create(plasma, &sequence);
    switch( A->dtyp ) {
    case PlasmaRealFloat:
        PLASMA_sTile_to_Lapack_Async( A, Af77, LDA, sequence, &request );
        break;
    case PlasmaRealDouble:
        PLASMA_dTile_to_Lapack_Async( A, Af77, LDA, sequence, &request );
        break;
    case PlasmaComplexFloat:
        PLASMA_cTile_to_Lapack_Async( A, Af77, LDA, sequence, &request );
        break;
    case PlasmaComplexDouble:
        PLASMA_zTile_to_Lapack_Async( A, Af77, LDA, sequence, &request );
        break;
    default:
        plasma_error("PLASMA_Tile_to_Lapack", "Type unknown");
    }
    plasma_dynamic_sync();
    status = sequence->status;
    plasma_sequence_destroy(plasma, sequence);
    return status;
}

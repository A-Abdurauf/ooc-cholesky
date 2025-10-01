/**
 *
 * @file plasma_f77.c
 *
 *  PLASMA auxiliary routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Bilel Hadri
 * @date 2010-11-15
 *
 **/
#include <stdlib.h>
#include "common.h"

#define PLASMA_INIT                PLASMA_FNAME(init               , INIT               )
#define PLASMA_FINALIZE            PLASMA_FNAME(finalize           , FINALIZE           )
#define PLASMA_ENABLE              PLASMA_FNAME(enable             , ENABLE             )
#define PLASMA_DISABLE             PLASMA_FNAME(disable            , DISABLE            )
#define PLASMA_SET                 PLASMA_FNAME(set                , SET                )
#define PLASMA_GET                 PLASMA_FNAME(get                , GET                )
#define PLASMA_DEALLOC_HANDLE      PLASMA_FNAME(dealloc_handle     , DEALLOC_HANDLE     )
#define PLASMA_DEALLOC_HANDLE_TILE PLASMA_FNAME(dealloc_handle_tile, DEALLOC_HANDLE_TILE)
#define PLASMA_VERSION             PLASMA_FNAME(version            , VERSION            )
#define PLASMA_DESC_CREATE         PLASMA_FNAME(desc_create        , DESC_CREATE        )
#define PLASMA_DESC_DESTROY        PLASMA_FNAME(desc_destroy       , DESC_DESTROY       )
#define PLASMA_LAPACK_TO_TILE      PLASMA_FNAME(lapack_to_tile     , LAPACK_TO_TILE     )
#define PLASMA_TILE_TO_LAPACK      PLASMA_FNAME(tile_to_lapack     , TILE_TO_LAPACK     )

#ifdef __cplusplus
extern "C" {
#endif

/***************************************************************************//**
 *  FORTRAN API - auxiliary function prototypes
 **/
void PLASMA_INIT(int *CORES, int *INFO)
{   *INFO = PLASMA_Init(*CORES); }

void PLASMA_FINALIZE(int *INFO)
{   *INFO = PLASMA_Finalize(); }

void PLASMA_ENABLE(PLASMA_enum *lever, int *INFO)
{   *INFO = PLASMA_Enable(*lever); }

void PLASMA_DISABLE(PLASMA_enum *lever, int *INFO)
{   *INFO = PLASMA_Disable(*lever); }

void PLASMA_SET(PLASMA_enum *param, int *value, int *INFO)
{   *INFO = PLASMA_Set(*param, *value); }

void PLASMA_GET(PLASMA_enum *param, int *value, int *INFO)
{   *INFO = PLASMA_Get(*param, value); }

void PLASMA_DEALLOC_HANDLE(intptr_t *sp, int *INFO)
{
    free((void *)(*sp));
    *INFO = PLASMA_SUCCESS;
}

void PLASMA_DEALLOC_HANDLE_TILE(PLASMA_desc **sp, int *INFO)
{
    PLASMA_Dealloc_Handle_Tile( sp );
    *INFO = PLASMA_SUCCESS;
}

void PLASMA_VERSION(int *VER_MAJOR, int *VER_MINOR, int *VER_MICRO, int *INFO)
{
    *VER_MAJOR = PLASMA_VERSION_MAJOR;
    *VER_MINOR = PLASMA_VERSION_MINOR;
    *VER_MICRO = PLASMA_VERSION_MICRO;
    *INFO = PLASMA_SUCCESS;
}

/***************************************************************************//**
 *  FORTRAN API - descriptor allocation and deallocation
 **/
void PLASMA_DESC_CREATE(PLASMA_desc **desc, void *mat, PLASMA_enum *dtyp, int *mb, int *nb, int *bsiz, int *lm, int *ln, int *i, int *j, int *m, int *n, int *INFO)
{   *INFO = PLASMA_Desc_Create(desc, mat, *dtyp, *mb, *nb, *bsiz, *lm, *ln, *i, *j, *m, *n); }

void PLASMA_DESC_DESTROY(PLASMA_desc **desc, int *INFO)
{   *INFO = PLASMA_Desc_Destroy(desc); }

/***************************************************************************//**
 *  FORTRAN API - conversion from LAPACK F77 matrix layout to tile layout
 **/
void PLASMA_LAPACK_TO_TILE(intptr_t *Af77, int *LDA, PLASMA_desc **A, int *INFO)
{   *INFO = PLASMA_Lapack_to_Tile( (void *)Af77, *LDA, *A); }

void PLASMA_TILE_TO_LAPACK(PLASMA_desc **A, intptr_t *Af77, int *LDA, int *INFO)
{   *INFO = PLASMA_Tile_to_Lapack( *A, (void *)Af77, *LDA); }

#ifdef __cplusplus
}
#endif

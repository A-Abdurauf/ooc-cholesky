/**
 *
 * @file workspace_d.c
 *
 *  PLASMA auxiliary routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Jakub Kurzak
 * @author Hatem Ltaief
 * @author Azzam Haidar
 * @date 2010-11-15
 * @generated d Fri Apr  1 11:03:01 2016
 *
 **/
#include "common.h"
#include "workspace.h"

/***************************************************************************//**
 *
 * @ingroup Auxiliary
 *
 *  PLASMA_Alloc_Workspace_dgebrd - Allocates workspace for PLASMA_dgebrd or PLASMA_dgebrd_Tile routine.
 *
 *******************************************************************************
 *
 * @param[in] M
 *          The number of rows of the matrix A. M >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrix A.  N >= 0.
 *
 * @param[out] descT
 *          On exit, workspace handle for storage of the extra T factors required by the tile BRD.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *
 ******************************************************************************/
int PLASMA_Alloc_Workspace_dgebrd(int M, int N, PLASMA_desc **descT) {
    return plasma_alloc_ibnb_tile(M, N, PLASMA_FUNC_DGESVD, PlasmaRealDouble, descT); }

/***************************************************************************//**
 *
 * @ingroup Auxiliary
 *
 *  PLASMA_Alloc_Workspace_dgeev - Allocates workspace for PLASMA_dgeev or PLASMA_dgeev_Tile routine.
 *
 *******************************************************************************
 *
 * @param[in] N
 *          The order of the matrix A.  N >= 0.
 *
 * @param[out] descT
 *          On exit, workspace handle for storage of the extra T factors required by the tile Hessenberg.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *
 ******************************************************************************/
int PLASMA_Alloc_Workspace_dgeev(int N, PLASMA_desc **descT) {
    return plasma_alloc_ibnb_tile(N, N, PLASMA_FUNC_DGEEV, PlasmaRealDouble, descT); }

/***************************************************************************//**
 *
 * @ingroup Auxiliary
 *
 *  PLASMA_Alloc_Workspace_dgehrd - Allocates workspace for PLASMA_dgehrd or PLASMA_dgehrd_Tile routine.
 *
 *******************************************************************************
 *
 * @param[in] N
 *          The number of columns of the matrix A.  N >= 0.
 *
 * @param[out] descT
 *          On exit, workspace handle for storage of the extra T factors required by the tile Hessenberg.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *
 ******************************************************************************/
int PLASMA_Alloc_Workspace_dgehrd(int N, PLASMA_desc **descT) {
    return plasma_alloc_ibnb_tile(N, N, PLASMA_FUNC_DGEHRD, PlasmaRealDouble, descT); }

/***************************************************************************//**
 *
 * @ingroup Auxiliary
 *
 *  PLASMA_Alloc_Workspace_dgels - Allocates workspace for PLASMA_dgels or PLASMA_dgels_Tile routine.
 *
 *******************************************************************************
 *
 * @param[in] M
 *          The number of rows of the matrix A. M >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrix A.  N >= 0.
 *
 * @param[out] descT
 *          On exit, workspace handle for storage of the extra T factors required by the tile QR
 *          or the tile LQ factorization.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *
 ******************************************************************************/
int PLASMA_Alloc_Workspace_dgels(int M, int N, PLASMA_desc **descT) {
    return plasma_alloc_ibnb_tile(M, N, PLASMA_FUNC_DGELS, PlasmaRealDouble, descT); }

/***************************************************************************//**
@@ -192,33 +120,6 @@
 *
 * @ingroup Auxiliary
 *
 *  PLASMA_Alloc_Workspace_dgels_Tile - Allocates tile workspace for PLASMA_dgels_Tile routine.
 *
 *******************************************************************************
 *
 * @param[in] M
 *          The number of rows of the matrix A. M >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrix A.  N >= 0.
 *
 * @param[out] descT
 *          On exit, PLASMA_desc on workspace handle for storage of the extra T factors required by the tile QR
 *          or the tile LQ factorization.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *
 ******************************************************************************/
int PLASMA_Alloc_Workspace_dgels_Tile(int M, int N, PLASMA_desc **descT) {
    return plasma_alloc_ibnb_tile(M, N, PLASMA_FUNC_DGELS, PlasmaRealDouble, descT); }

/***************************************************************************//**
 *
 * @ingroup Auxiliary
 *
 *  PLASMA_Alloc_Workspace_dgeqrf - Allocates workspace for PLASMA_dgeqrf or PLASMA_dgeqrf_Tile routine.
 *
 *******************************************************************************
 *
 * @param[in] M
 *          The number of rows of the matrix A. M >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrix A.  N >= 0.
 *
 * @param[out] descT
 *          On exit, workspace handle for storage of the extra T factors required by the tile QR
 *          factorization.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *
 ******************************************************************************/
int PLASMA_Alloc_Workspace_dgeqrf(int M, int N, PLASMA_desc **descT) {
    return plasma_alloc_ibnb_tile(M, N, PLASMA_FUNC_DGELS, PlasmaRealDouble, descT); }

/***************************************************************************//**
 *
 * @ingroup Auxiliary
 *
 *  PLASMA_Alloc_Workspace_dgeqrf_Tile - Allocates tile workspace for PLASMA_dgels_Tile routine.
 *
 *******************************************************************************
 *
 * @param[in] M
 *          The number of rows of the matrix A. M >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrix A.  N >= 0.
 *
 * @param[out] descT
 *          On exit, PLASMA_desc on workspace handle for storage of the extra T factors required by the tile QR
 *          or the tile LQ factorization.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *
 ******************************************************************************/
int PLASMA_Alloc_Workspace_dgeqrf_Tile(int M, int N, PLASMA_desc **descT) {
    return plasma_alloc_ibnb_tile(M, N, PLASMA_FUNC_DGELS, PlasmaRealDouble, descT); }

/***************************************************************************//**
 *
 * @ingroup Auxiliary
 *
 *  PLASMA_Alloc_Workspace_dgelqf - Allocates workspace for PLASMA_dgelqf or PLASMA_dgelqf_Tile routines.
 *
 *******************************************************************************
 *
 * @param[in] M
 *          The number of rows of the matrix A. M >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrix A.  N >= 0.
 *
 * @param[out] descT
 *          On exit, workspace handle for storage of the extra T factors required by the tile LQ
 *          factorization.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *
 ******************************************************************************/
int PLASMA_Alloc_Workspace_dgelqf(int M, int N, PLASMA_desc **descT) {
    return plasma_alloc_ibnb_tile(M, N, PLASMA_FUNC_DGELS, PlasmaRealDouble, descT); }

/***************************************************************************//**
 *
 * @ingroup Auxiliary
 *
 *  PLASMA_Alloc_Workspace_dgelqf_Tile - Allocates tile workspace for PLASMA_dgels_Tile routine.
 *
 *******************************************************************************
 *
 * @param[in] M
 *          The number of rows of the matrix A. M >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrix A.  N >= 0.
 *
 * @param[out] descT
 *          On exit, PLASMA_desc on workspace handle for storage of the extra T factors required by the tile QR
 *          or the tile LQ factorization.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *
 ******************************************************************************/
int PLASMA_Alloc_Workspace_dgelqf_Tile(int M, int N, PLASMA_desc **descT) {
    return plasma_alloc_ibnb_tile(M, N, PLASMA_FUNC_DGELS, PlasmaRealDouble, descT); }

/***************************************************************************//**
 *
 * @ingroup Auxiliary
 *
 *  PLASMA_Alloc_Workspace_dgesdd - Allocates workspace for PLASMA_dgesdd or PLASMA_dgesdd_Tile routine.
 *
 *******************************************************************************
 *
 * @param[in] M
 *          The number of rows of the matrix A. M >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrix A.  N >= 0.
 *
 * @param[out] descT
 *          On exit, workspace handle for storage of the extra T factors required by the tile BRD.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *
 ******************************************************************************/
int PLASMA_Alloc_Workspace_dgesdd(int M, int N, PLASMA_desc **descT) {
    return plasma_alloc_ibnb_tile(M, N, PLASMA_FUNC_DGESVD, PlasmaRealDouble, descT); }

/***************************************************************************//**
 *
 * @ingroup Auxiliary
 *
 *  PLASMA_Alloc_Workspace_dgesv - Allocates workspace for PLASMA_dgesv or PLASMA_dgesv_Tile routines.
 *
 *******************************************************************************
 *
 * @param[in] N
 *          The number of linear equations, i.e., the order of the matrix A. N >= 0.
 *
 * @param[out] descL
 *          On exit, workspace handle for storage of the extra L factors required by the tile LU
 *          factorization.
 *
 * @param[out] IPIV
 *          On exit, workspace handle for storage of pivot indexes required by the tile LU
 *          factorization (not equivalent to LAPACK).
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *
 ******************************************************************************/
int PLASMA_Alloc_Workspace_dgesv_incpiv(int N, PLASMA_desc **descL, int **IPIV) {
    int status = plasma_alloc_ibnb_tile(N, N, PLASMA_FUNC_DGESV, PlasmaRealDouble, descL);
    if (status != PLASMA_SUCCESS)
        return status;
    return plasma_alloc_ipiv(N, N, PLASMA_FUNC_DGESV, (void**)IPIV); }

/***************************************************************************//**
 *
 * @ingroup Auxiliary
 *
 *  PLASMA_Alloc_Workspace_dgesv_Tile - Allocates workspace for PLASMA_dgesv_Tile routines.
 *
 *******************************************************************************
 *
 * @param[in] N
 *          The number of linear equations, i.e., the order of the matrix A. N >= 0.
 *
 * @param[out] descL
 *          On exit, PLASMA descriptor on workspace handle for storage of the extra
 *          L factors required by the tile LU factorization.
 *
 * @param[out] IPIV
 *          On exit, workspace handle for storage of pivot indexes required by the tile LU
 *          factorization (not equivalent to LAPACK).
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *
 ******************************************************************************/
int PLASMA_Alloc_Workspace_dgesv_incpiv_Tile(int N, PLASMA_desc **descL, int **IPIV)
{
    int status = plasma_alloc_ibnb_tile(N, N, PLASMA_FUNC_DGESV, PlasmaRealDouble, descL);
    if (status != PLASMA_SUCCESS)
        return status;
    return plasma_alloc_ipiv(N, N, PLASMA_FUNC_DGESV, (void **)IPIV);
}
/***************************************************************************//**
 *
 * @ingroup Auxiliary
 *
 *  PLASMA_Alloc_Workspace_dgesvd - Allocates workspace for PLASMA_dgesvd or PLASMA_dgesvd_Tile routine.
 *
 *******************************************************************************
 *
 * @param[in] M
 *          The number of rows of the matrix A. M >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrix A.  N >= 0.
 *
 * @param[out] descT
 *          On exit, workspace handle for storage of the extra T factors required by the tile BRD.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *
 ******************************************************************************/
int PLASMA_Alloc_Workspace_dgesvd(int M, int N, PLASMA_desc **descT) {
    return plasma_alloc_ibnb_tile(M, N, PLASMA_FUNC_DGESVD, PlasmaRealDouble, descT); }

/***************************************************************************//**
 * @ingroup Auxiliary
 *
 *  PLASMA_Alloc_Workspace_dgetrf_incpiv - Allocates workspace for
 *  PLASMA_dgetrf_incpiv or PLASMA_dgetrf_incpiv_Tile or
 *  PLASMA_dgetrf_incpiv_Tile_Async routines.
 *
 *******************************************************************************
 *
 * @param[in] M
 *          The number of rows of the matrix A. M >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrix A.  N >= 0.
 *
 * @param[out] descL
 *          On exit, workspace handle for storage of the extra L factors required by the tile LU
 *          factorization.
 *
 * @param[out] IPIV
 *          On exit, workspace handle for storage of pivot indexes required by the tile LU
 *          factorization (not equivalent to LAPACK).
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *
 ******************************************************************************
 *
 * @sa PLASMA_dgetrf_incpiv
 * @sa PLASMA_dgetrf_incpiv_Tile
 * @sa PLASMA_dgetrf_incpiv_Tile_Async
 *
 ******************************************************************************/
int PLASMA_Alloc_Workspace_dgetrf_incpiv(int M, int N, PLASMA_desc **descL, int **IPIV) {
    int status = plasma_alloc_ibnb_tile(M, N, PLASMA_FUNC_DGESV, PlasmaRealDouble, descL);
    if (status != PLASMA_SUCCESS)
        return status;
    return plasma_alloc_ipiv(M, N, PLASMA_FUNC_DGESV, (void**)IPIV); }

/***************************************************************************//**
 *
 * @ingroup Auxiliary
 *
 *  PLASMA_Alloc_Workspace_dgetrf_incpiv_Tile - Allocates workspace for
 *  PLASMA_dgesv_incpiv_Tile or PLASMA_dgesv_incpiv_Tile_Async routines.
 *
 *******************************************************************************
 *
 * @param[in] N
 *          The number of linear equations, i.e., the order of the matrix A. N >= 0.
 *
 * @param[out] descL
 *          On exit, PLASMA descriptor on workspace handle for storage of the extra
 *          L factors required by the tile LU factorization.
 *
 * @param[out] IPIV
 *          On exit, workspace handle for storage of pivot indexes required by the tile LU
 *          factorization (not equivalent to LAPACK).
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *
 ******************************************************************************/
int PLASMA_Alloc_Workspace_dgetrf_incpiv_Tile(int N, PLASMA_desc **descL, int **IPIV)
{
    int status = plasma_alloc_ibnb_tile(N, N, PLASMA_FUNC_DGESV, PlasmaRealDouble, descL);
    if (status != PLASMA_SUCCESS)
        return status;
    return plasma_alloc_ipiv(N, N, PLASMA_FUNC_DGESV, (void **)IPIV);
}
/***************************************************************************//**
 *
 * @ingroup Auxiliary
 *
 *  PLASMA_Alloc_Workspace_dsyev - Allocates workspace for PLASMA_dsyev or PLASMA_dsyev_Tile routine.
 *
 *******************************************************************************
 *
 * @param[in] M
 *          The number of rows of the matrix A. M >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrix A.  N >= 0.
 *
 * @param[out] descT
 *          On exit, workspace handle for storage of the extra T factors required by the tile TRD.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *
 ******************************************************************************/
int PLASMA_Alloc_Workspace_dsyev(int M, int N, PLASMA_desc **descT) {
    return plasma_alloc_ibnb_tile(M, N, PLASMA_FUNC_DSYEV, PlasmaRealDouble, descT); }

/***************************************************************************//**
 *
 * @ingroup Auxiliary
 *
 *  PLASMA_Alloc_Workspace_dsyevd - Allocates workspace for PLASMA_dsyevd or PLASMA_dsyevd_Tile routine.
 *
 *******************************************************************************
 *
 * @param[in] M
 *          The number of rows of the matrix A. M >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrix A.  N >= 0.
 *
 * @param[out] descT
 *          On exit, workspace handle for storage of the extra T factors required by the tile TRD.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *
 ******************************************************************************/
int PLASMA_Alloc_Workspace_dsyevd(int M, int N, PLASMA_desc **descT) {
    return plasma_alloc_ibnb_tile(M, N, PLASMA_FUNC_DSYEVD, PlasmaRealDouble, descT); }

/***************************************************************************//**
 *
 * @ingroup Auxiliary
 *
 *  PLASMA_Alloc_Workspace_dsyevr - Allocates workspace for PLASMA_dsyevr or PLASMA_dsyevr_Tile routine.
 *
 *******************************************************************************
 *
 * @param[in] M
 *          The number of rows of the matrix A. M >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrix A.  N >= 0.
 *
 * @param[out] descT
 *          On exit, workspace handle for storage of the extra T factors required by the tile TRD.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *
 ******************************************************************************/
int PLASMA_Alloc_Workspace_dsyevr(int M, int N, PLASMA_desc **descT) {
    return plasma_alloc_ibnb_tile(M, N, PLASMA_FUNC_DSYEVD, PlasmaRealDouble, descT); }


/***************************************************************************//**
 *
 * @ingroup Auxiliary
 *
 *  PLASMA_Alloc_Workspace_dsygv - Allocates workspace for PLASMA_dsygv or PLASMA_dsygv_Tile routine.
 *
 *******************************************************************************
 *
 * @param[in] M
 *          The number of rows of the matrix A. M >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrix A.  N >= 0.
 *
 * @param[out] descT
 *          On exit, workspace handle for storage of the extra T factors required by the tile TRD.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *
 ******************************************************************************/
int PLASMA_Alloc_Workspace_dsygv(int M, int N, PLASMA_desc **descT) {
    return plasma_alloc_ibnb_tile(M, N, PLASMA_FUNC_DSYGV, PlasmaRealDouble, descT); }

/***************************************************************************//**
 *
 * @ingroup Auxiliary
 *
 *  PLASMA_Alloc_Workspace_dsygvd - Allocates workspace for PLASMA_dsygvd or PLASMA_dsygvd_Tile routine.
 *
 *******************************************************************************
 *
 * @param[in] M
 *          The number of rows of the matrix A. M >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrix A.  N >= 0.
 *
 * @param[out] descT
 *          On exit, workspace handle for storage of the extra T factors required by the tile TRD.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *
 ******************************************************************************/
int PLASMA_Alloc_Workspace_dsygvd(int M, int N, PLASMA_desc **descT) {
    return plasma_alloc_ibnb_tile(M, N, PLASMA_FUNC_DSYGVD, PlasmaRealDouble, descT); }

 /***************************************************************************//**
 *
 * @ingroup Auxiliary
 *
 *  PLASMA_Alloc_Workspace_dsytrd - Allocates workspace for PLASMA_dsytrd or PLASMA_dsytrd_Tile routine.
 *
 *******************************************************************************
 *
 * @param[in] M
 *          The number of rows of the matrix A. M >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrix A.  N >= 0.
 *
 * @param[out] descT
 *          On exit, workspace handle for storage of the extra T factors required by the tile TRD.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *
 ******************************************************************************/
int PLASMA_Alloc_Workspace_dsytrd(int M, int N, PLASMA_desc **descT) {
    return plasma_alloc_ibnb_tile(M, N, PLASMA_FUNC_DSYTRD, PlasmaRealDouble, descT); }

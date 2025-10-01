/**
 *
 * @file workspace_s.c
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
 * @generated s Fri Apr  1 11:03:01 2016
 *
 **/
#include "common.h"
#include "workspace.h"

/***************************************************************************//**
 *
 * @ingroup Auxiliary
 *
 *  PLASMA_Alloc_Workspace_sgebrd - Allocates workspace for PLASMA_sgebrd or PLASMA_sgebrd_Tile routine.
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
int PLASMA_Alloc_Workspace_sgebrd(int M, int N, PLASMA_desc **descT) {
    return plasma_alloc_ibnb_tile(M, N, PLASMA_FUNC_SGESVD, PlasmaRealFloat, descT); }

/***************************************************************************//**
 *
 * @ingroup Auxiliary
 *
 *  PLASMA_Alloc_Workspace_sgeev - Allocates workspace for PLASMA_sgeev or PLASMA_sgeev_Tile routine.
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
int PLASMA_Alloc_Workspace_sgeev(int N, PLASMA_desc **descT) {
    return plasma_alloc_ibnb_tile(N, N, PLASMA_FUNC_SGEEV, PlasmaRealFloat, descT); }

/***************************************************************************//**
 *
 * @ingroup Auxiliary
 *
 *  PLASMA_Alloc_Workspace_sgehrd - Allocates workspace for PLASMA_sgehrd or PLASMA_sgehrd_Tile routine.
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
int PLASMA_Alloc_Workspace_sgehrd(int N, PLASMA_desc **descT) {
    return plasma_alloc_ibnb_tile(N, N, PLASMA_FUNC_SGEHRD, PlasmaRealFloat, descT); }

/***************************************************************************//**
 *
 * @ingroup Auxiliary
 *
 *  PLASMA_Alloc_Workspace_sgels - Allocates workspace for PLASMA_sgels or PLASMA_sgels_Tile routine.
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
int PLASMA_Alloc_Workspace_sgels(int M, int N, PLASMA_desc **descT) {
    return plasma_alloc_ibnb_tile(M, N, PLASMA_FUNC_SGELS, PlasmaRealFloat, descT); }

/***************************************************************************//**
@@ -192,33 +120,6 @@
 *
 * @ingroup Auxiliary
 *
 *  PLASMA_Alloc_Workspace_sgels_Tile - Allocates tile workspace for PLASMA_sgels_Tile routine.
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
int PLASMA_Alloc_Workspace_sgels_Tile(int M, int N, PLASMA_desc **descT) {
    return plasma_alloc_ibnb_tile(M, N, PLASMA_FUNC_SGELS, PlasmaRealFloat, descT); }

/***************************************************************************//**
 *
 * @ingroup Auxiliary
 *
 *  PLASMA_Alloc_Workspace_sgeqrf - Allocates workspace for PLASMA_sgeqrf or PLASMA_sgeqrf_Tile routine.
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
int PLASMA_Alloc_Workspace_sgeqrf(int M, int N, PLASMA_desc **descT) {
    return plasma_alloc_ibnb_tile(M, N, PLASMA_FUNC_SGELS, PlasmaRealFloat, descT); }

/***************************************************************************//**
 *
 * @ingroup Auxiliary
 *
 *  PLASMA_Alloc_Workspace_sgeqrf_Tile - Allocates tile workspace for PLASMA_sgels_Tile routine.
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
int PLASMA_Alloc_Workspace_sgeqrf_Tile(int M, int N, PLASMA_desc **descT) {
    return plasma_alloc_ibnb_tile(M, N, PLASMA_FUNC_SGELS, PlasmaRealFloat, descT); }

/***************************************************************************//**
 *
 * @ingroup Auxiliary
 *
 *  PLASMA_Alloc_Workspace_sgelqf - Allocates workspace for PLASMA_sgelqf or PLASMA_sgelqf_Tile routines.
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
int PLASMA_Alloc_Workspace_sgelqf(int M, int N, PLASMA_desc **descT) {
    return plasma_alloc_ibnb_tile(M, N, PLASMA_FUNC_SGELS, PlasmaRealFloat, descT); }

/***************************************************************************//**
 *
 * @ingroup Auxiliary
 *
 *  PLASMA_Alloc_Workspace_sgelqf_Tile - Allocates tile workspace for PLASMA_sgels_Tile routine.
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
int PLASMA_Alloc_Workspace_sgelqf_Tile(int M, int N, PLASMA_desc **descT) {
    return plasma_alloc_ibnb_tile(M, N, PLASMA_FUNC_SGELS, PlasmaRealFloat, descT); }

/***************************************************************************//**
 *
 * @ingroup Auxiliary
 *
 *  PLASMA_Alloc_Workspace_sgesdd - Allocates workspace for PLASMA_sgesdd or PLASMA_sgesdd_Tile routine.
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
int PLASMA_Alloc_Workspace_sgesdd(int M, int N, PLASMA_desc **descT) {
    return plasma_alloc_ibnb_tile(M, N, PLASMA_FUNC_SGESVD, PlasmaRealFloat, descT); }

/***************************************************************************//**
 *
 * @ingroup Auxiliary
 *
 *  PLASMA_Alloc_Workspace_sgesv - Allocates workspace for PLASMA_sgesv or PLASMA_sgesv_Tile routines.
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
int PLASMA_Alloc_Workspace_sgesv_incpiv(int N, PLASMA_desc **descL, int **IPIV) {
    int status = plasma_alloc_ibnb_tile(N, N, PLASMA_FUNC_SGESV, PlasmaRealFloat, descL);
    if (status != PLASMA_SUCCESS)
        return status;
    return plasma_alloc_ipiv(N, N, PLASMA_FUNC_SGESV, (void**)IPIV); }

/***************************************************************************//**
 *
 * @ingroup Auxiliary
 *
 *  PLASMA_Alloc_Workspace_sgesv_Tile - Allocates workspace for PLASMA_sgesv_Tile routines.
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
int PLASMA_Alloc_Workspace_sgesv_incpiv_Tile(int N, PLASMA_desc **descL, int **IPIV)
{
    int status = plasma_alloc_ibnb_tile(N, N, PLASMA_FUNC_SGESV, PlasmaRealFloat, descL);
    if (status != PLASMA_SUCCESS)
        return status;
    return plasma_alloc_ipiv(N, N, PLASMA_FUNC_SGESV, (void **)IPIV);
}
/***************************************************************************//**
 *
 * @ingroup Auxiliary
 *
 *  PLASMA_Alloc_Workspace_sgesvd - Allocates workspace for PLASMA_sgesvd or PLASMA_sgesvd_Tile routine.
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
int PLASMA_Alloc_Workspace_sgesvd(int M, int N, PLASMA_desc **descT) {
    return plasma_alloc_ibnb_tile(M, N, PLASMA_FUNC_SGESVD, PlasmaRealFloat, descT); }

/***************************************************************************//**
 * @ingroup Auxiliary
 *
 *  PLASMA_Alloc_Workspace_sgetrf_incpiv - Allocates workspace for
 *  PLASMA_sgetrf_incpiv or PLASMA_sgetrf_incpiv_Tile or
 *  PLASMA_sgetrf_incpiv_Tile_Async routines.
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
 * @sa PLASMA_sgetrf_incpiv
 * @sa PLASMA_sgetrf_incpiv_Tile
 * @sa PLASMA_sgetrf_incpiv_Tile_Async
 *
 ******************************************************************************/
int PLASMA_Alloc_Workspace_sgetrf_incpiv(int M, int N, PLASMA_desc **descL, int **IPIV) {
    int status = plasma_alloc_ibnb_tile(M, N, PLASMA_FUNC_SGESV, PlasmaRealFloat, descL);
    if (status != PLASMA_SUCCESS)
        return status;
    return plasma_alloc_ipiv(M, N, PLASMA_FUNC_SGESV, (void**)IPIV); }

/***************************************************************************//**
 *
 * @ingroup Auxiliary
 *
 *  PLASMA_Alloc_Workspace_sgetrf_incpiv_Tile - Allocates workspace for
 *  PLASMA_sgesv_incpiv_Tile or PLASMA_sgesv_incpiv_Tile_Async routines.
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
int PLASMA_Alloc_Workspace_sgetrf_incpiv_Tile(int N, PLASMA_desc **descL, int **IPIV)
{
    int status = plasma_alloc_ibnb_tile(N, N, PLASMA_FUNC_SGESV, PlasmaRealFloat, descL);
    if (status != PLASMA_SUCCESS)
        return status;
    return plasma_alloc_ipiv(N, N, PLASMA_FUNC_SGESV, (void **)IPIV);
}
/***************************************************************************//**
 *
 * @ingroup Auxiliary
 *
 *  PLASMA_Alloc_Workspace_ssyev - Allocates workspace for PLASMA_ssyev or PLASMA_ssyev_Tile routine.
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
int PLASMA_Alloc_Workspace_ssyev(int M, int N, PLASMA_desc **descT) {
    return plasma_alloc_ibnb_tile(M, N, PLASMA_FUNC_SSYEV, PlasmaRealFloat, descT); }

/***************************************************************************//**
 *
 * @ingroup Auxiliary
 *
 *  PLASMA_Alloc_Workspace_ssyevd - Allocates workspace for PLASMA_ssyevd or PLASMA_ssyevd_Tile routine.
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
int PLASMA_Alloc_Workspace_ssyevd(int M, int N, PLASMA_desc **descT) {
    return plasma_alloc_ibnb_tile(M, N, PLASMA_FUNC_SSYEVD, PlasmaRealFloat, descT); }

/***************************************************************************//**
 *
 * @ingroup Auxiliary
 *
 *  PLASMA_Alloc_Workspace_ssyevr - Allocates workspace for PLASMA_ssyevr or PLASMA_ssyevr_Tile routine.
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
int PLASMA_Alloc_Workspace_ssyevr(int M, int N, PLASMA_desc **descT) {
    return plasma_alloc_ibnb_tile(M, N, PLASMA_FUNC_SSYEVD, PlasmaRealFloat, descT); }


/***************************************************************************//**
 *
 * @ingroup Auxiliary
 *
 *  PLASMA_Alloc_Workspace_ssygv - Allocates workspace for PLASMA_ssygv or PLASMA_ssygv_Tile routine.
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
int PLASMA_Alloc_Workspace_ssygv(int M, int N, PLASMA_desc **descT) {
    return plasma_alloc_ibnb_tile(M, N, PLASMA_FUNC_SSYGV, PlasmaRealFloat, descT); }

/***************************************************************************//**
 *
 * @ingroup Auxiliary
 *
 *  PLASMA_Alloc_Workspace_ssygvd - Allocates workspace for PLASMA_ssygvd or PLASMA_ssygvd_Tile routine.
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
int PLASMA_Alloc_Workspace_ssygvd(int M, int N, PLASMA_desc **descT) {
    return plasma_alloc_ibnb_tile(M, N, PLASMA_FUNC_SSYGVD, PlasmaRealFloat, descT); }

 /***************************************************************************//**
 *
 * @ingroup Auxiliary
 *
 *  PLASMA_Alloc_Workspace_ssytrd - Allocates workspace for PLASMA_ssytrd or PLASMA_ssytrd_Tile routine.
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
int PLASMA_Alloc_Workspace_ssytrd(int M, int N, PLASMA_desc **descT) {
    return plasma_alloc_ibnb_tile(M, N, PLASMA_FUNC_SSYTRD, PlasmaRealFloat, descT); }

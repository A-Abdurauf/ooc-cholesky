/**
 *
 * @file workspace_c.c
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
 * @generated c Fri Apr  1 11:03:01 2016
 *
 **/
#include "common.h"
#include "workspace.h"

/***************************************************************************//**
 *
 * @ingroup Auxiliary
 *
 *  PLASMA_Alloc_Workspace_cgebrd - Allocates workspace for PLASMA_cgebrd or PLASMA_cgebrd_Tile routine.
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
int PLASMA_Alloc_Workspace_cgebrd(int M, int N, PLASMA_desc **descT) {
    return plasma_alloc_ibnb_tile(M, N, PLASMA_FUNC_CGESVD, PlasmaComplexFloat, descT); }

/***************************************************************************//**
 *
 * @ingroup Auxiliary
 *
 *  PLASMA_Alloc_Workspace_cgeev - Allocates workspace for PLASMA_cgeev or PLASMA_cgeev_Tile routine.
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
int PLASMA_Alloc_Workspace_cgeev(int N, PLASMA_desc **descT) {
    return plasma_alloc_ibnb_tile(N, N, PLASMA_FUNC_CGEEV, PlasmaComplexFloat, descT); }

/***************************************************************************//**
 *
 * @ingroup Auxiliary
 *
 *  PLASMA_Alloc_Workspace_cgehrd - Allocates workspace for PLASMA_cgehrd or PLASMA_cgehrd_Tile routine.
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
int PLASMA_Alloc_Workspace_cgehrd(int N, PLASMA_desc **descT) {
    return plasma_alloc_ibnb_tile(N, N, PLASMA_FUNC_CGEHRD, PlasmaComplexFloat, descT); }

/***************************************************************************//**
 *
 * @ingroup Auxiliary
 *
 *  PLASMA_Alloc_Workspace_cgels - Allocates workspace for PLASMA_cgels or PLASMA_cgels_Tile routine.
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
int PLASMA_Alloc_Workspace_cgels(int M, int N, PLASMA_desc **descT) {
    return plasma_alloc_ibnb_tile(M, N, PLASMA_FUNC_CGELS, PlasmaComplexFloat, descT); }

/***************************************************************************//**
@@ -192,33 +120,6 @@
 *
 * @ingroup Auxiliary
 *
 *  PLASMA_Alloc_Workspace_cgels_Tile - Allocates tile workspace for PLASMA_cgels_Tile routine.
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
int PLASMA_Alloc_Workspace_cgels_Tile(int M, int N, PLASMA_desc **descT) {
    return plasma_alloc_ibnb_tile(M, N, PLASMA_FUNC_CGELS, PlasmaComplexFloat, descT); }

/***************************************************************************//**
 *
 * @ingroup Auxiliary
 *
 *  PLASMA_Alloc_Workspace_cgeqrf - Allocates workspace for PLASMA_cgeqrf or PLASMA_cgeqrf_Tile routine.
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
int PLASMA_Alloc_Workspace_cgeqrf(int M, int N, PLASMA_desc **descT) {
    return plasma_alloc_ibnb_tile(M, N, PLASMA_FUNC_CGELS, PlasmaComplexFloat, descT); }

/***************************************************************************//**
 *
 * @ingroup Auxiliary
 *
 *  PLASMA_Alloc_Workspace_cgeqrf_Tile - Allocates tile workspace for PLASMA_cgels_Tile routine.
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
int PLASMA_Alloc_Workspace_cgeqrf_Tile(int M, int N, PLASMA_desc **descT) {
    return plasma_alloc_ibnb_tile(M, N, PLASMA_FUNC_CGELS, PlasmaComplexFloat, descT); }

/***************************************************************************//**
 *
 * @ingroup Auxiliary
 *
 *  PLASMA_Alloc_Workspace_cgelqf - Allocates workspace for PLASMA_cgelqf or PLASMA_cgelqf_Tile routines.
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
int PLASMA_Alloc_Workspace_cgelqf(int M, int N, PLASMA_desc **descT) {
    return plasma_alloc_ibnb_tile(M, N, PLASMA_FUNC_CGELS, PlasmaComplexFloat, descT); }

/***************************************************************************//**
 *
 * @ingroup Auxiliary
 *
 *  PLASMA_Alloc_Workspace_cgelqf_Tile - Allocates tile workspace for PLASMA_cgels_Tile routine.
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
int PLASMA_Alloc_Workspace_cgelqf_Tile(int M, int N, PLASMA_desc **descT) {
    return plasma_alloc_ibnb_tile(M, N, PLASMA_FUNC_CGELS, PlasmaComplexFloat, descT); }

/***************************************************************************//**
 *
 * @ingroup Auxiliary
 *
 *  PLASMA_Alloc_Workspace_cgesdd - Allocates workspace for PLASMA_cgesdd or PLASMA_cgesdd_Tile routine.
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
int PLASMA_Alloc_Workspace_cgesdd(int M, int N, PLASMA_desc **descT) {
    return plasma_alloc_ibnb_tile(M, N, PLASMA_FUNC_CGESVD, PlasmaComplexFloat, descT); }

/***************************************************************************//**
 *
 * @ingroup Auxiliary
 *
 *  PLASMA_Alloc_Workspace_cgesv - Allocates workspace for PLASMA_cgesv or PLASMA_cgesv_Tile routines.
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
int PLASMA_Alloc_Workspace_cgesv_incpiv(int N, PLASMA_desc **descL, int **IPIV) {
    int status = plasma_alloc_ibnb_tile(N, N, PLASMA_FUNC_CGESV, PlasmaComplexFloat, descL);
    if (status != PLASMA_SUCCESS)
        return status;
    return plasma_alloc_ipiv(N, N, PLASMA_FUNC_CGESV, (void**)IPIV); }

/***************************************************************************//**
 *
 * @ingroup Auxiliary
 *
 *  PLASMA_Alloc_Workspace_cgesv_Tile - Allocates workspace for PLASMA_cgesv_Tile routines.
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
int PLASMA_Alloc_Workspace_cgesv_incpiv_Tile(int N, PLASMA_desc **descL, int **IPIV)
{
    int status = plasma_alloc_ibnb_tile(N, N, PLASMA_FUNC_CGESV, PlasmaComplexFloat, descL);
    if (status != PLASMA_SUCCESS)
        return status;
    return plasma_alloc_ipiv(N, N, PLASMA_FUNC_CGESV, (void **)IPIV);
}
/***************************************************************************//**
 *
 * @ingroup Auxiliary
 *
 *  PLASMA_Alloc_Workspace_cgesvd - Allocates workspace for PLASMA_cgesvd or PLASMA_cgesvd_Tile routine.
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
int PLASMA_Alloc_Workspace_cgesvd(int M, int N, PLASMA_desc **descT) {
    return plasma_alloc_ibnb_tile(M, N, PLASMA_FUNC_CGESVD, PlasmaComplexFloat, descT); }

/***************************************************************************//**
 * @ingroup Auxiliary
 *
 *  PLASMA_Alloc_Workspace_cgetrf_incpiv - Allocates workspace for
 *  PLASMA_cgetrf_incpiv or PLASMA_cgetrf_incpiv_Tile or
 *  PLASMA_cgetrf_incpiv_Tile_Async routines.
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
 * @sa PLASMA_cgetrf_incpiv
 * @sa PLASMA_cgetrf_incpiv_Tile
 * @sa PLASMA_cgetrf_incpiv_Tile_Async
 *
 ******************************************************************************/
int PLASMA_Alloc_Workspace_cgetrf_incpiv(int M, int N, PLASMA_desc **descL, int **IPIV) {
    int status = plasma_alloc_ibnb_tile(M, N, PLASMA_FUNC_CGESV, PlasmaComplexFloat, descL);
    if (status != PLASMA_SUCCESS)
        return status;
    return plasma_alloc_ipiv(M, N, PLASMA_FUNC_CGESV, (void**)IPIV); }

/***************************************************************************//**
 *
 * @ingroup Auxiliary
 *
 *  PLASMA_Alloc_Workspace_cgetrf_incpiv_Tile - Allocates workspace for
 *  PLASMA_cgesv_incpiv_Tile or PLASMA_cgesv_incpiv_Tile_Async routines.
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
int PLASMA_Alloc_Workspace_cgetrf_incpiv_Tile(int N, PLASMA_desc **descL, int **IPIV)
{
    int status = plasma_alloc_ibnb_tile(N, N, PLASMA_FUNC_CGESV, PlasmaComplexFloat, descL);
    if (status != PLASMA_SUCCESS)
        return status;
    return plasma_alloc_ipiv(N, N, PLASMA_FUNC_CGESV, (void **)IPIV);
}
/***************************************************************************//**
 *
 * @ingroup Auxiliary
 *
 *  PLASMA_Alloc_Workspace_cheev - Allocates workspace for PLASMA_cheev or PLASMA_cheev_Tile routine.
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
int PLASMA_Alloc_Workspace_cheev(int M, int N, PLASMA_desc **descT) {
    return plasma_alloc_ibnb_tile(M, N, PLASMA_FUNC_CHEEV, PlasmaComplexFloat, descT); }

/***************************************************************************//**
 *
 * @ingroup Auxiliary
 *
 *  PLASMA_Alloc_Workspace_cheevd - Allocates workspace for PLASMA_cheevd or PLASMA_cheevd_Tile routine.
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
int PLASMA_Alloc_Workspace_cheevd(int M, int N, PLASMA_desc **descT) {
    return plasma_alloc_ibnb_tile(M, N, PLASMA_FUNC_CHEEVD, PlasmaComplexFloat, descT); }

/***************************************************************************//**
 *
 * @ingroup Auxiliary
 *
 *  PLASMA_Alloc_Workspace_cheevr - Allocates workspace for PLASMA_cheevr or PLASMA_cheevr_Tile routine.
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
int PLASMA_Alloc_Workspace_cheevr(int M, int N, PLASMA_desc **descT) {
    return plasma_alloc_ibnb_tile(M, N, PLASMA_FUNC_CHEEVD, PlasmaComplexFloat, descT); }


/***************************************************************************//**
 *
 * @ingroup Auxiliary
 *
 *  PLASMA_Alloc_Workspace_chegv - Allocates workspace for PLASMA_chegv or PLASMA_chegv_Tile routine.
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
int PLASMA_Alloc_Workspace_chegv(int M, int N, PLASMA_desc **descT) {
    return plasma_alloc_ibnb_tile(M, N, PLASMA_FUNC_CHEGV, PlasmaComplexFloat, descT); }

/***************************************************************************//**
 *
 * @ingroup Auxiliary
 *
 *  PLASMA_Alloc_Workspace_chegvd - Allocates workspace for PLASMA_chegvd or PLASMA_chegvd_Tile routine.
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
int PLASMA_Alloc_Workspace_chegvd(int M, int N, PLASMA_desc **descT) {
    return plasma_alloc_ibnb_tile(M, N, PLASMA_FUNC_CHEGVD, PlasmaComplexFloat, descT); }

 /***************************************************************************//**
 *
 * @ingroup Auxiliary
 *
 *  PLASMA_Alloc_Workspace_chetrd - Allocates workspace for PLASMA_chetrd or PLASMA_chetrd_Tile routine.
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
int PLASMA_Alloc_Workspace_chetrd(int M, int N, PLASMA_desc **descT) {
    return plasma_alloc_ibnb_tile(M, N, PLASMA_FUNC_CHETRD, PlasmaComplexFloat, descT); }

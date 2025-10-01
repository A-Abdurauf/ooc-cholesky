/**
 *
 * @file core_dsygst.c
 *
 *  PLASMA core_blas kernel
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Hatem Ltaief
 * @date 2010-11-15
 * @generated d Fri Apr  1 11:02:30 2016
 *
 **/
#ifdef PLASMA_WITH_MKL
#include <mkl_lapacke.h>
#else
#include <lapacke.h>
#endif
#include "common.h"

/***************************************************************************//**
 *
 * @ingroup CORE_double
 *
 *  CORE_dsygst - reduces a complex Hermitian-definite generalized
 *  eigenproblem to standard form.
 *  If PlasmaItype == 1, the problem is A*x = lambda*B*x, and A is
 *  overwritten by inv(U**T)*A*inv(U) or inv(L)*A*inv(L**T)
 *  If PlasmaItype == 2 or 3, the problem is A*B*x = lambda*x or B*A*x
 *  = lambda*x, and A is overwritten by U*A*U**T or L**T*A*L.  B must
 *  have been previously factorized as U**T*U or L*L**T by
 *  CORE_dpotrf.
 *
 *******************************************************************************
 *
 * @param[in] itype
 *          Intended usage:
 *          = 1: A*x=(lambda)*B*x
 *          = 2: A*Bx=(lambda)*x
 *          = 3: B*A*x=(lambda)*x
 *
 * @param[in] uplo
 *          Specifies whether the matrix A is upper triangular or
 *          lower triangular:
 *          = PlasmaUpper: Upper triangle of A is stored;
 *          = PlasmaLower: Lower triangle of A is stored.
 *
 * @param[in] N
 *          The order of the matrices A and B. N >= 0.
 *
 * @param[in,out] A
 *          On entry, the symmetric (or Hermitian) matrix A.
 *          If uplo = PlasmaUpper, the leading N-by-N upper triangular
 *          part of A contains the upper triangular part of the matrix
 *          A, and the strictly lower triangular part of A is not
 *          referenced.
 *          If uplo = PlasmaLower, the leading N-by-N lower triangular
 *          part of A contains the lower triangular part of the matrix
 *          A, and the strictly upper triangular part of A is not
 *          referenced.
 *          On exit, if return value == 0, the transformed matrix,
 *          stored in the same format as A.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,N).
 *
 * @param[in,out] B
 *          On entry, the triangular factor from the Cholesky
 *          factorization of B, as returned by PLASMA_DPOTRF.
 *
 * @param[in] LDB
 *          The leading dimension of the array B. LDB >= max(1,N).
 *
 * @param[out] INFO
 *          - 0 on successful exit
 *          - <0 if -i, the i-th argument had an illegal value
 *
 ******************************************************************************/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_dsygst = PCORE_dsygst
#define CORE_dsygst PCORE_dsygst
#endif
void CORE_dsygst(int itype, PLASMA_enum uplo, int N,
                 double *A, int LDA,
                 double *B, int LDB, int *INFO)
{
    *INFO = LAPACKE_dsygst_work(
        LAPACK_COL_MAJOR,
        itype,
        lapack_const(uplo),
        N, A, LDA, B, LDB );
}

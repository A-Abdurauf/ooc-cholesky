/**
 * @file core_dstedc.c
 *
 *  PLASMA computational routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Gr�goire Pichon
 * @author Azzam Haidar
 * @date 2014-07
 * @generated d Fri Apr  1 11:02:33 2016
 *
 **/
#include "common.h"
#ifdef PLASMA_WITH_MKL
#include <mkl_lapacke.h>
#else
#include <lapacke.h>
#endif

#define REAL
#undef COMPLEX

/**
 *******************************************************************************
 *
 * @ingroup CORE_double
 *
 *  CORE_dstedc - solves the symmetric tridiagonal eigensystem using Divide &
 *  Conquer
 *
 *******************************************************************************
 *
 * @param[in] compz
 *          = PlasmaNoVec: computes eigenvalues only.
 *          = PlasmaVec: computes eigenpairs of the original symmetric
 *            matrix. On entry, Z must contain the orthogonal matrix used to
 *            reduce the original matrix to tridiagonal form.
 *          = PlasmaIVec: computes eigenpairs of the tridiagonal matrix. Z is
 *            initialized to the Identity Matrix.
 *
 * @param[in] n
 *          n specifies the order of the matrix. N >= 0
 *
 * @param[in,out] D
 *          On entry, D contains the diagonal elements of the tridiagonal matrix.
 *          On exit, D contains the eigenvalues sorted into increasing order.
 *
 * @param[in] E
 *          On entry, E contains the extra-diagonal elements of the tridiagonal
 *          matrix.
 *          On exit, E is destroyed.
 *
 * @param[in,out] Z
 *          On entry, Z has to be set to the Identity matrix.
 *          On exit, Z contains the eigenvectors.
 *
 * @param[in] LDZ
 *          LDZ specifies the leading direction of Z
 *
 * @param[in] WORK
 *          plasma_complex64_t workspace
 *
 * @param[in] LWORK
 *          Size of plasma_complex64_t workspace
 *
 * @param[in] RWORK
 *           workspace
 *
 * @param[in] LRWORK
 *          Size of double workspace
 *
 * @param[in] IWORK
 *          Integer workspace
 *
 * @param[in] LIWORK
 *          Size of integer workspace
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *
 *******************************************************************************/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_dstedc = PCORE_dstedc
#define CORE_dstedc PCORE_dstedc
#endif
int CORE_dstedc( PLASMA_enum compz, int n,
                 double *D, double *E,
                 double *Z, int LDZ,
                 double *WORK, int LWORK,
#ifdef COMPLEX
                 double *RWORK, int LRWORK,
#endif
                 int *IWORK, int LIWORK )
{
    int info;

    if (WORK == NULL){
        info = LAPACKE_dstedc( LAPACK_COL_MAJOR, lapack_const(compz),
                               n, D, E, Z, LDZ );
    } else {
        info = LAPACKE_dstedc_work( LAPACK_COL_MAJOR,
                                    lapack_const(compz),
                                    n, D, E, Z, LDZ,
                                    WORK, LWORK,
#ifdef COMPLEX
                                    RWORK, LRWORK,
#endif
                                    IWORK, LIWORK );
    }

    assert(!info);
    return info;
}

/**
 * @file core_slaed3_updatevectors.c
 *
 *  PLASMA computational routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Gregoire Pichon
 * @author Azzam Haidar
 * @date 2014-07
 * @generated s Fri Apr  1 11:02:33 2016
 *
 **/
#ifdef PLASMA_WITH_MKL
#include <mkl_lapacke.h>
#else
#include <lapacke.h>
#endif
#include <math.h>
#include "common.h"

/**
 *******************************************************************************
 *
 * @ingroup CORE_float
 *
 *  CORE_slaed3_updatevectors -
 *
 *******************************************************************************
 *
 * @param[in] op
 *          Type of operation to apply
 *          = PlasmaLaed3Update1: Apply the GEMM on the first subproblem
 *          = PlasmaLaed3Update2: Apply the GEMM on the second subproblem and
 *          merge the eigenvalues.
 *          = PlasmaLaed3UpdateAll: Apply all operations in one call.
 *
 * @param[in] wsmode
 *          Specifies the amount of extra workspace used for the computations.
 *          = 0 : means that no extra workspace is available (WORK = NULL).
 *            S = Q2 + n12 * n1 + n23 * (n-n1) is used as the workspace to
 *            sequentially apply the updates
 *          = 1 : means that extra workspace has been allocated in
 *            CORE_slaed3_computevectors() to apply the updates in parallel. If
 *            op = PlasmaLaed3Update1, the first subproblem is updated. If op =
 *            PlasmaLaed3Update2, the second subproblem is updated.
 *          = 3 : means that extra workspace has been allocated in after K has
 *            been computed in CORE_slaed2_computeK(). Allows to apply both
 *            updates in parallel as previously.
 *
 * @param[in] n
 *          n specifies the dimension of the symmetric tridiagonal matrix
 *
 * @param[in] n1
 *          n1 specifies the location of the last eigenvalue of the first subproblem
 *          min(1, n) <= n1 <= n/2
 *
 * @param[in] K
 *          K specifies the number of non-deflated eigenvalues
 *
 * @param[in] il_nondef
 *          The first eigenvector index when computing only a subset of all eigenpairs
 *
 * @param[in] iu_nondef
 *          The last eigenvector index when computing only a subset of all eigenpairs
 *
 * @param[out] Q
 *          The current eigenvectors
 *
 * @param[in] ldq
 *          LDQ specifies the leading direction of Q
 *
 * @param[in,out] Q2
 *          The updated eigenvectors
 *
 * @param[in] ctot
 *          ctot indicates the number of vectors of each type in the Q matrix.
 *          0- non-zero in the upper half only
 *          1- dense
 *          2- non-zero in the lower half only
 *          3- deflated
 *
 * @param[in] W
 *          W is the place were will be stored the updated eigenvectors
 *
 * @param[in] start
 *          start specifies the first column index to be considered by this
 *          kernel
 *
 * @param[in] end
 *          end specifies the last column index to be considered by this kernel
 *
 ***************************************************************************/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_slaed3_updatevectors = PCORE_slaed3_updatevectors
#define CORE_slaed3_updatevectors PCORE_slaed3_updatevectors
#endif
void CORE_slaed3_updatevectors( int op, int wsmode,
                                int n, int n1, int K, int il_nondef, int iu_nondef,
                                float *Q, int ldq, float *Q2,
                                const int *ctot, float *W,
                                int start, int end)
{
    float *last_wsvector;
    float *S;
    int n2, n12, n23, size, sze1;
    int shift, LDS;

    /* Adjust the subset */
    start = max( 0,   max( start, il_nondef ) );
    end   = min( end, min( K,     iu_nondef ) );
    size  = max( 0, end-start );
    sze1  = size;

    /* Quick return */
    if ( K == 0 ){
        return;
    }

    if( start > K )
        return;

    /*
     * Update
     */
    n2  = n-n1;
    n12 = ctot[0] + ctot[1];
    n23 = ctot[1] + ctot[2];

    /* In sequentiel use n23 then n12 in spite of LDS */
    LDS   = wsmode == 0 ? max(n12, n23) : K;
    shift = wsmode == 3 ? K * start : 0;

    /*****************************************
     *               GEMM 2
     *****************************************
     */

    /*
     * If wsmode = 0, then the n-by-n workspace used gor the computation is
     * insufficient to store the deflated eigenvectors and the data before
     * applying the updates.
     * In LAPACK slaed3 routines, it is possible to do it without extra
     * workspace because this one is already allocaet with an extra space of
     * size n used in other auxiliary routines for the eigenvalue problem. In
     * our case, as we used the Q for the final user as the workspace, we cannot
     * provide this extra space of size n. We handle it by allocating here and
     * splitting the update for the last column.
     */
    if (wsmode == 0 && end == n){
        last_wsvector = malloc(n*sizeof(float));
        sze1--;
    }
    else {
        last_wsvector = NULL;
    }

    if( op & PlasmaLaed3Update2 ){
        float *lQ  = Q + ctot[0] + ldq * start;
        float *lQ2 = Q2 + n1 * n12;

        /* copy to S */
        if( wsmode == 0 ) {
            S = Q2 + n12 * n1 + n23 *n2 + LDS * start;

            LAPACKE_slacpy_work(LAPACK_COL_MAJOR,
                                lapack_const(PlasmaUpperLower), n23, sze1,
                                lQ, ldq, S, LDS);

            if (last_wsvector){
                LAPACKE_slacpy_work(LAPACK_COL_MAJOR,
                                    lapack_const(PlasmaUpperLower), n23, 1,
                                    lQ + ldq * sze1, ldq, last_wsvector, 1);
            }
        }
        else {
            S = W + shift + ctot[0];
        }

        /* do GEMM 2 */
        lQ = Q + n1 + start * ldq;
        if (n23 != 0){
            cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                        n2, sze1, n23,
                        1., lQ2, n2,
                            S,   LDS,
                        0., lQ,  ldq);

            if (last_wsvector){
                cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                            n2, 1, n23,
                            1., lQ2,           n2,
                                last_wsvector, n23,
                            0., lQ+sze1*ldq,   ldq);
            }
        }
        else{
            LAPACKE_slaset_work(LAPACK_COL_MAJOR, 'A', n2, size,
                                0., 0., lQ, ldq);
        }
    }


    /*****************************************
     *               GEMM 1
     *****************************************
     */
    if( op & PlasmaLaed3Update1 ){
        float *lQ = Q + ldq * start;

        /* copy to S */
        if( wsmode == 0 ){
            S = Q2 + n12 * n1 + n23 * n2 + LDS * start;

            LAPACKE_slacpy_work(LAPACK_COL_MAJOR,
                                lapack_const(PlasmaUpperLower), n12, sze1,
                                lQ, ldq, S, LDS);

            if (last_wsvector){
                LAPACKE_slacpy_work(LAPACK_COL_MAJOR,
                                    lapack_const(PlasmaUpperLower), n12, 1,
                                    lQ+ldq*sze1, ldq, last_wsvector, 1);
            }
        }
        else {
            S = W + shift;
        }

        /* do GEMM 1 */
        if (n12 != 0){
            cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                        n1, sze1, n12,
                        1., Q2,          n1,
                            S,           LDS,
                        0., lQ, ldq);

           if (last_wsvector){
                cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                            n1, 1, n12,
                            1., Q2,            n1,
                                last_wsvector, n12,
                            0., lQ+ldq*sze1,   ldq);
            }
        }
        else{
            LAPACKE_slaset_work(LAPACK_COL_MAJOR, 'A', n1, size,
                                0., 0., lQ, ldq);
        }
    }

    if (last_wsvector){
        free(last_wsvector);
    }
}

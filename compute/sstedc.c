/**
 * @file sstedc.c
 *
 *  PLASMA computational routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Gregoire Pichon
 * @author Azzam Haidar
 * @date 2014-07
 * @generated s Fri Apr  1 11:02:56 2016
 *
 **/
#include <math.h>
#ifdef PLASMA_WITH_MKL
#include <mkl_lapacke.h>
#else
#include <lapacke.h>
#endif
#include "common.h"

#undef COMPLEX
#define REAL

/**
 *****************************************************************************
 *
 * @ingroup float_Tile
 *
 *  PLASMA_sstedc - Computes all eigenpairs of a symmetric tridiagonal matrix
 *
 *******************************************************************************
 *
 * @param[in] jobz
 *          Intended usage:
 *          = PlasmaIVec: computes eigenpairs of the symmetric tridiagonal matrix
 *          = PlasmaVec: computes eigenpairs of the original matrix (not supported now)
 *
 * @param[in] n
 *          n specifies the dimension of the original matrix
 *
 * @param[in,out] D
 *          On entry, D contains the diagonal elements of the tridiagonal matrix
 *          On exit, D contains the eigenvalues
 *
 * @param[in] E
 *          On entry, E contains the extra-diagonal elements of the tridiagonal matrix
 *
 * @param[out] Z
 *          On exit, if jobz = PlasmaVec and info = 0, the eigenvectors.
 *
 * @param[in] LDZ
 *          The leading dimention of the eigenvectors matrix Z. LDZ >= max(1,N).
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *
 *******************************************************************************
 *
 * @sa PLASMA_sstedc
 * @sa PLASMA_sstedc_Async
 * @sa PLASMA_cstedc
 * @sa PLASMA_sstedc
 * @sa PLASMA_sstedc
 *
 ******************************************************************************/
int PLASMA_sstedc(PLASMA_enum jobz, int n,
                  float *D, float *E,
                  float *Z, int LDZ)
{
    plasma_context_t *plasma;
    PLASMA_sequence *sequence = NULL;
    PLASMA_request request = PLASMA_REQUEST_INITIALIZER;
    int status;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_sstedc", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    plasma_sequence_create(plasma, &sequence);
    PLASMA_sstedc_Async(jobz, n, D, E, Z, LDZ, sequence, &request);
    plasma_dynamic_sync();
    status = sequence->status;
    plasma_sequence_destroy(plasma, sequence);
    return status;
}

/**
 *****************************************************************************
 *
 * @ingroup float_Tile_Async
 *
 *  PLASMA_sstedc_Async - Computes all eigenpairs of a symmetric tridiagonal matrix
 *
 *******************************************************************************
 *
 * @param[in] jobz
 *          Intended usage:
 *          = PlasmaIVec: computes eigenpairs of the symmetric tridiagonal matrix
 *          = PlasmaVec: computes eigenpairs of the original matrix (not supported now)
 *
 * @param[in] n
 *          n specifies the dimension of the original matrix
 *
 * @param[in,out] D
 *          On entry, D contains the diagonal elements of the tridiagonal matrix
 *          On exit, D contains the eigenvalues
 *
 * @param[in] E
 *          On entry, E contains the extra-diagonal elements of the tridiagonal matrix
 *
 * @param[out] Z
 *          On exit, if jobz = PlasmaVec and info = 0, the eigenvectors.
 *
 * @param[in] LDZ
 *          The leading dimention of the eigenvectors matrix Z. LDZ >= max(1,N).
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
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *
 *******************************************************************************/
#ifdef REAL
int PLASMA_sstedc_Async(PLASMA_enum jobz, int n,
                        float *D, float *E,
                        float *Z, int LDZ,
                        PLASMA_sequence *sequence, PLASMA_request *request)
{
    int info = 0;
    int SMLSIZ;
    plasma_context_t *plasma;

    /* Variables for sorting eigenpairs */
    int act_perm = 0;

    plasma = plasma_context_self();

    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_sstedc_Async", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        plasma_fatal_error("PLASMA_sstedc_Async", "NULL sequence");
        return PLASMA_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        plasma_fatal_error("PLASMA_sstedc_Async", "NULL request");
        return PLASMA_ERR_UNALLOCATED;
    }
    /* Check sequence status */
    if (sequence->status == PLASMA_SUCCESS)
        request->status = PLASMA_SUCCESS;
    else
        return plasma_request_fail(sequence, request, PLASMA_ERR_SEQUENCE_FLUSHED);

    SMLSIZ = plasma->ev_smlsze;

    /* Check input arguments */
    if (jobz != PlasmaNoVec && jobz != PlasmaIvec) {
        plasma_error("PLASMA_sstedc_Async", "illegal value of jobz");
        return -1;
    }

    if (n < 0){
        plasma_error("PLASMA_sstedc_Async", "illegal value of n");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }

    if (D == NULL){
        plasma_error("PLASMA_sstedc_Async", "illegal value of D");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }

    if (E == NULL){
        plasma_error("PLASMA_sstedc_Async", "illegal value of E");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }

    if (Z == NULL && (jobz == PlasmaVec) ){
        plasma_error("PLASMA_sstedc_Async", "illegal value of Z");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }

    if (LDZ < max(1, n)){
        plasma_error("PLASMA_sstedc_Async", "illegal value of LDZ");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }
    if ((n < SMLSIZ) || (jobz == PlasmaNoVec)) {
        if (n < SMLSIZ) LAPACKE_slaset_work( LAPACK_COL_MAJOR, lapack_const(PlasmaUpperLower), n, n, 0.0, 1.0, Z, LDZ);
        info = LAPACKE_sstedc( LAPACK_COL_MAJOR, lapack_const(jobz),
                               n, D, E, Z, LDZ);
        if (info != 0){
            plasma_error("PLASMA_sstedc_Async", "LAPACKE sstedc failed");
            return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
        }
        else{
            return PLASMA_SUCCESS;
        }
    }

    plasma_time_t timestedc, timeonesolve, timeallsolve, timeswap, timesort;
    timestedc = plasma_gettime(plasma);

    int work_pos   = 0;
    int work_pos2  = 0;
    float *WORK   = malloc(n*n*sizeof(float));
    float *WORK2  = malloc((4*n)*sizeof(float));
    int *IWORK     = malloc((5*n)*sizeof(int));
    int *localdata = malloc(n*sizeof(int));
    int LDWORK = LDZ;

    float eps = LAPACKE_slamch_work('e');
    float tiny;
    int current = 0;
    int nsubsml = 0, nsubdlaed0 = 0, nsub1 = 0;

    int il = 0;
    int iu = n;
    float vl = 0., vu = 0.;
    char range = 'A';

    int start, size;
    int i = 0;

    memset(localdata, 0, n*sizeof(int));

    /* Set WORK to identity */
    plasma_dynamic_call_5( plasma_pslaset_identity,
                           int, n,
                           float*, WORK,
                           int, n,
                           PLASMA_sequence*, sequence,
                           PLASMA_request*, request);
    plasma_dynamic_sync();

    timeallsolve = plasma_gettime(plasma);

#if defined(ENABLE_DEBUG2)
    plasma_dynamic_sync();
    printf("start sstedc\n");
#endif


    /******************************************************
     * loop over the main or possible-subproblem and solve
     * ****************************************************
     * */
    i = 0;
    while (i<n-1){
        timeonesolve = plasma_gettime(plasma);

        tiny = eps*sqrt(fabs(D[i]))*sqrt(fabs(D[i+1]));

        if ((fabs(E[i]) <= tiny) || (i==n-2)){
            start = current;
            size = i-current+1;

            /* The last index (n-1) is not treated in the loop  */
            if (i == n-2){
                size++;
            }

            if (size == 1){
                nsub1 += 1;
            }

            else if (size < SMLSIZ){
                nsubsml += 1;
                plasma_dynamic_call_8( plasma_psstedc,
                                       PLASMA_enum, jobz,
                                       int,     size,
                                       float*, D+start,
                                       float*, E+start,
                                       float*, WORK+n*start+start,
                                       int,     n,
                                       PLASMA_sequence*, sequence,
                                       PLASMA_request*, request);
            }
            else{

#if defined(ENABLE_DEBUG2)
                plasma_dynamic_sync();
                printf("  start  solving subproblems of size %5d\n",size);
#endif


                nsubdlaed0 += 1;

                plasma_dynamic_call_21( plasma_pslaed0,
                                        int, 2, /* int or Plasma_enum ??? */
                                        char, range,
                                        int, size,
                                        int, size,
                                        float*, D+start,
                                        float*, E+start,
                                        int, il,
                                        int, iu,
                                        float, vl,
                                        float, vu,
                                        float*, WORK+n*start+start,
                                        int, n,
                                        float*, NULL,
                                        int, size,
                                        float*, Z+LDZ*start+start,
                                        float*, WORK2+work_pos2,
                                        int, LDWORK,
                                        int*, IWORK+5*start,
                                        int*, localdata+start,
                                        PLASMA_sequence*, sequence,
                                        PLASMA_request*, request);

                work_pos2 += 4*size;
                work_pos  += size*size;
            }

            current += size;

            timeonesolve = plasma_gettime(plasma)-timeonesolve;
            plasma_printtime("  Finish solving subproblems of size %5d timing= %lf \n",size, timeonesolve);
        }    /* End solving one independant subproblem */

        i++;
    } /* End While */

    timeallsolve = plasma_gettime(plasma)-timeallsolve;
    plasma_printtime("  Finish all solve dlaed0  nbsub %5d     timing= %lf \n", nsubdlaed0, timeallsolve);

    /* Wait for the end of each independant subproblem */
    plasma_dynamic_sync();

    timesort = plasma_gettime(plasma);

    /* Create the permutation to sort D eigenvalues into increasing order */
    CORE_slapst(PlasmaIncreasingOrder, n, D, IWORK);
    memcpy(WORK2, D, n*sizeof(float));
    for (i=0; i<n; i++){
        if (IWORK[i] != i){
            act_perm = 1;
            D[i] = WORK2[IWORK[i]];
        }
    }

    timeswap = plasma_gettime(plasma);

    act_perm = 1;                 /* always copy back to Q */
    if (act_perm == 1){
        plasma_dynamic_call_7(plasma_psswaps,
                              int, n,
                              int*, IWORK,
                              float*, Z,
                              int, LDZ,
                              float*, WORK,
                              PLASMA_sequence*, sequence,
                              PLASMA_request*, request);
    }

    plasma_dynamic_sync();

    timesort  = timeswap-timesort;
    timeswap  = plasma_gettime(plasma)-timeswap;
    timestedc = plasma_gettime(plasma)-timestedc;
    plasma_printtime("  Finish sort                              timing= %lf \n", timesort);
    plasma_printtime("  Finish swap                              timing= %lf \n", timeswap);
    plasma_printtime("  Finish sstedc with nsub_dlaed0 %5d     nsub_smlsiz %5d   nsub1 %5d timing= %lf \n",nsubdlaed0, nsubsml, nsub1, timestedc);

    free(localdata);
    free(WORK);
    free(WORK2);
    free(IWORK);
    return info;
}





#else  /* COMPLEX */
int PLASMA_sstedc_Async(PLASMA_enum jobz, int n,
                        float *D, float *E,
                        float *Z, int LDZ,
                        PLASMA_sequence *sequence, PLASMA_request *request)
{
    plasma_context_t *plasma;
    plasma = plasma_context_self();

    float *Q = malloc(n*n*sizeof(float));
    memset(Q, 0, n*n*sizeof(float));

    PLASMA_sstedc_Async(jobz, n, D, E,
                        Q, n,
                        sequence, request);

    plasma_dynamic_call_8(plasma_pslag2c,
                          int,     n,
                          int,     n,
                          float*,             Q,
                          int,                 n,
                          float*, Z,
                          int,                 LDZ,
                          PLASMA_sequence*, sequence,
                          PLASMA_request*, request);

    // TODO: need to be removed by using dependencies on Q
    QUARK_Barrier(plasma->quark);
    free(Q);

    return 0;
}

#endif

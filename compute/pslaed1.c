/**
 * @file pslaed1.c
 *
 *  PLASMA computational routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Gregoire Pichon
 * @author Azzam Haidar
 * @date 2014-07
 * @generated s Fri Apr  1 11:03:00 2016
 *
 **/
#include "common.h"
#include <math.h>
#include <stdlib.h>

/***************************************************************************//**
 *  Parallel routine slaed1 using divide and conquer for finding eigenvalue of
 *  tridiagonal matrix dynamic scheduling
 **/

/**
 ******************************************************************************
 *
 * plasma_pslaed1_quark - Merge two subproblems after splitting
 *
 *******************************************************************************
 *
 * @param[in] n
 *          n specifies the dimension of the original matrix
 *
 * @param[in] n1
 *          n1 specifies the location of the last eigenvalue of the first
 *          subproblem
 *          min(1, n) <= n1 <= n/2
 *
 * @param[in,out] D
 *          D refers to the eigenvalues
 *
 * @param[in,out] Q
 *          Q refers to the eigenvectors
 *
 * @param[in,out] INDXQ
 *          On entry, the permutation which separately sorts the two
 *          subproblems in D into ascending order.
 *          On exit, the permutation which will reintegrate the
 *          subproblems back into sorted order,
 *          i.e. D( INDXQ( I = 1, N ) ) will be in ascending order.
 *
 * @param[in,out] beta
 *          The subdiagonal entry used to create the rank-1 modification.
 *
 * @param[in,out] work
 *          real workspace
 *
 * @param[in,out] iwork
 *          integer workspace
 *
 * @param[in,out] K1
 *          K1 points to the localdata of the first suproblem
 *          K1[0] corresponds to the number of non-deflated eigenvalues within
 *          last merge
 *
 * @param[in] K2
 *          K2 points to the localdata of the second suproblem
 *
 * @param[in] sequence
 *          Identifies the sequence of function calls that this call belongs to
 *          (for completion checks and exception handling purposes).
 *
 * @param[out] request
 *          Identifies this function call (for exception handling purposes).
 *
 *******************************************************************************/
void plasma_pslaed1_quark(char range,
                          int n, int n1,
                          float *D,
                          int il, int iu,
                          float vl, float vu,
                          float *Q, int LDQ,
                          int *INDXQ, float *beta,
                          float *work, float *work2, int *iwork,
                          int *K1, int *K2,
                          int last_merge,
                          PLASMA_sequence *sequence, PLASMA_request *request)
{
    plasma_context_t *plasma;
    Quark_Task_Flags task_flags = Quark_Task_Flags_Initializer;
    int i, nb, start, end, task_size;
    int wsmode;
    float *Z=NULL;
    float *DLAMBDA=NULL;
    float *W=NULL;
    float *Q2=NULL;
    int *INDX;
    int *INDXC;
    int *INDXP;
    int *COLTYP;
    int nb_tasks;
    float **Qcurr;

    plasma = plasma_context_self();
    if (sequence->status != PLASMA_SUCCESS)
        return;
    QUARK_Task_Flag_Set(&task_flags, TASK_SEQUENCE, (intptr_t)sequence->quark_sequence);

    nb     = plasma->ev_tasknb; //DnC_TASK_NB;
    wsmode = plasma->ev_wsmode;
    if (plasma->world_size == 1) nb = n;

#if defined(ENABLE_DEBUG3)
    QUARK_Barrier(plasma->quark);
    printf("      start  slaed1 %d \n",n);
#endif

    Z = work2;
    DLAMBDA = Z + n;
    W = DLAMBDA + n;
    Q2 = work;

    INDX   = iwork;
    INDXC  = INDX + n;
    COLTYP = INDXC + n;
    INDXP  = COLTYP + n;

    nb_tasks = plasma_ceildiv(n,nb);
    /* TODO: free S. Can be done with sort in GEMM or SCRATCH in reduceW */
    Qcurr = malloc(nb_tasks * sizeof(float*));
    memset( Qcurr, 0, nb_tasks * sizeof(float*));

    /* for the reduction, a NB_TASKS*n matrix will be used */
    /* each compute_W task will compute the local reduction */
    /* the reduce_W task will merge results */
    /* W_red is free at the end of reduction */
    float *W_red = malloc(nb_tasks*n*sizeof(float));

    /* Parameters when computing a subset of eigenpairs */
    int *il_nondef = NULL;
    int *iu_nondef = NULL;

    /*
     * Computes K, the number of non-deflated eigenvalues
     * The two subproblems have to be computed before
     */
    QUARK_CORE_slaed2_computeK(plasma->quark, &task_flags,
                               /* K handles the dependency on the subproblem 1 */
                               K1, n, n1,
                               beta, D, Q, LDQ,
                               Z, DLAMBDA, W,
                               INDX, INDXC, INDXP, INDXQ,
                               COLTYP,
                               Qcurr, wsmode,
                               /* K2 handles the dependency on the subproblem 2 */
                               K2);

    /*
     * Now that the merge is prepared, let's merge the two subproblems in
     * parallel
     *
     * When wsmode == 3, We can perfom in parallel laed2_compressq and laed4, so
     * we split the three operations.
     * Otherwise, we submit the pipelined kernel that combines both operation in
     * sequence.
     */
    if( wsmode == 3 ) {
        for(i=0; i<n; i+=nb){
            task_size = min(nb, n-i);

            QUARK_CORE_slaed2_compressq(plasma->quark, &task_flags,
                                        n, n1, i, i+task_size,
                                        INDX, COLTYP,
                                        Q, LDQ, Q2, K1);

            QUARK_CORE_slaed4_p2f1(plasma->quark, &task_flags,
                                   n, K1, D, beta, Qcurr, K1,
                                   DLAMBDA, W, INDX,
                                   i, i+task_size,
                                   sequence, request,
                                   /* Fake dependency for next operation on the panel */
                                   Qcurr + (i/nb), INOUT);

            QUARK_CORE_slaed3_compW_p2f1(plasma->quark, &task_flags,
                                         n, K1, Qcurr, K1,
                                         DLAMBDA, W_red+n*(i/nb), INDX,
                                         i, i+task_size,
                                         /* Fake dependency for next operation on the panel */
                                         Qcurr + (i/nb), INOUT,
                                         /* Dependency for the barrier represented by the reduce operation */
                                         W_red, OUTPUT | GATHERV);
        }

        QUARK_CORE_slaed3_reduceW_p2(plasma->quark, &task_flags,
                                     n, n1, K1, nb_tasks,
                                     Qcurr, K1, W_red, W);
    }
    else {
        /*
         * No need to submit independently the 3 kernels, things will be done in
         * a pipeline, we just create parallelism by panel
         */
        for(i=0; i<n; i+=nb){
            task_size = min(nb, n-i);

            /* K1 is used as dependency before the barrier */
            QUARK_CORE_slaed1_pipelined(plasma->quark, &task_flags,
                                        n, n1, K1, INDX, COLTYP,
                                        D, beta, Q, LDQ, Q2,
                                        DLAMBDA, W, W_red + n*(i/nb),
                                        i, i+task_size);
        }

        /* BARRIER */
        QUARK_CORE_slaed3_reduceW(plasma->quark, &task_flags,
                                  n, n1, K1, nb_tasks,
                                  Q, LDQ, W_red, W);
    }

    /* Free W_red */
    QUARK_CORE_free( plasma->quark, &task_flags, W_red, sizeof(float) );

    /* Submit the full copy of the deflated vectors back to Q */
    for(start=0; start<n; start+=nb){
        end = min(start + nb, n);
        QUARK_CORE_slaed2_copydef(plasma->quark, &task_flags,
                                  n, n1, K1, COLTYP,
                                  Q, LDQ, Q2,
                                  start, end);
    }

    if(wsmode == 0){
        QUARK_CORE_sDC_fakedep(plasma->quark, &task_flags,
                               nb_tasks, nb, Q, LDQ, W);
    }

    /*
     * The problems have been merged together, we know need to compute the
     * eigenvectors of the merged problem and update the solution
     * If wsmode = 0, we have no workspaces to split the updates in parallel, so
     * we submit the pipelined laed3 kernel,
     * Otherwise we can submit different kernels for each operation as in the
     * previous loop.
     */
    for(start=0, i=0; start<n; start+=nb, i++){
        end = min(start + nb, n);
        task_size = end - start;

        if( wsmode == 0 ){
            QUARK_CORE_slaed3_pipelined(plasma->quark, &task_flags,
                                        n, n1, K1, il_nondef, iu_nondef,
                                        D, Q, LDQ, Q2,
                                        INDXC, INDXQ, COLTYP, W,
                                        start, end);
        }
        else {
            /* Compute eigenvectors of the rank-1 modification */
            QUARK_CORE_slaed3_computevectors(plasma->quark, &task_flags,
                                             wsmode, n, K1, il_nondef, iu_nondef,
                                             Q, LDQ, W, INDXC,
                                             Qcurr,
                                             Qcurr + i,
                                             start, end);

            if (wsmode == 1) {
                /* Allocate ws and copy part of Q into it for update */
                QUARK_CORE_slaed3_wscopy(plasma->quark, &task_flags,
                                         K1, il_nondef, iu_nondef,
                                         Q, LDQ,
                                         Qcurr + i,
                                         start, end);
            }

            QUARK_CORE_slaed3_updatevectors(
                plasma->quark, &task_flags,
                PlasmaLaed3Update1, wsmode, n, n1, K1, il_nondef, iu_nondef,
                D, Q, LDQ, Q2,
                INDXQ, COLTYP,
                (wsmode == 1 ? Qcurr + i : Qcurr),
                start, end, Qcurr + i);

            QUARK_CORE_slaed3_updatevectors(
                plasma->quark, &task_flags,
                PlasmaLaed3Update2, wsmode, n, n1, K1, il_nondef, iu_nondef,
                D, Q, LDQ, Q2,
                INDXQ, COLTYP,
                (wsmode == 1 ? Qcurr + i : Qcurr),
                start, end, Qcurr + i);

            if(wsmode == 1){
                QUARK_CORE_slaed3_freebigwork(plasma->quark, &task_flags, K1, 1, Qcurr + i);
            }
        }
    }

    if(wsmode == 1){
        QUARK_CORE_slaed3_freebigwork(plasma->quark, &task_flags, K1, 5, Qcurr );
    }
    if(wsmode == 3){
        QUARK_CORE_slaed3_freebigwork(plasma->quark, &task_flags, K1, 3, Qcurr );
    }
}

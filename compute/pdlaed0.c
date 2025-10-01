/**
 * @file pdlaed0.c
 *
 *  PLASMA computational routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Gr�goire Pichon
 * @author Azzam Haidar
 * @date 2014-07
 * @precisions normal d -> s
 *
 **/
#ifdef PLASMA_WITH_MKL
#include <mkl_lapacke.h>
#else
#include <lapacke.h>
#endif
#include <math.h>
#include "common.h"

static double done = 1.0;

/***************************************************************************//**
 *  Parallel routine dlaed0 using divide and conquer for finding eigenvalue of tridiagonal matrix
 *  dynamic scheduling
 **/

/**
 ******************************************************************************
 *
 * plasma_pdlaed0_quark - Apply divide and conquer algorithm to an independent subproblem
 *
 *******************************************************************************
 *
 * @param[in] icompq
 *          Intended usage:
 *          = 0: compute eigenvalues only (not supported now)
 *          = 1: compute eigenpairs of the original dense symmetric matrix (not supported now)
 *          = 2: compute eigenpairs of the symmetric tridiagonal matrix
 *
 * @param[in] qsiz
 *          qsiz specifies the dimension of the orthogonal matrix used to reduce
 *          the full matrix to tridiagonal form. Used only if icompq=1
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
 * @param[in,out] Q
 *          In case icompq = 2:
 *          On entry, Q is the identity matrix
 *          On exit, Q will contain the eigenvectors of the tridiagonal matrix
 *
 * @param[in] LDQ
 *          The leading dimention of the eigenvectors matrix Z. LDZ >= max(1,N).
 *
 * @param[in,out] qstore
 *          used only is icompq = 1
 *
 * @param[in] ldqs
 *          used only is icompq = 1
 *
 * @param[in,out] WORK
 *          real workspace
 *
 * @param[in,out] IWORK
 *          integer workspace
 *
 * @param[in,out] localdata
 *         localdata is used for containing informations about all supbroblems
 *
 * @param[in] sequence
 *          Identifies the sequence of function calls that this call belongs to
 *          (for completion checks and exception handling purposes).
 *
 * @param[out] request
 *          Identifies this function call (for exception handling purposes).
 *
 *******************************************************************************/
void plasma_pdlaed0_quark(int icompq, char range,
                          int qsiz, int n,
                          double *D, double *E,
                          int il, int iu,
                          double vl, double vu,
                          double *Q, int LDQ, double *qstore,
                          int ldqs, double *WORK, double *WORK2, int LDWORK,
                          int *IWORK, int *localdata,
                          PLASMA_sequence *sequence, PLASMA_request *request)
{
    plasma_context_t *plasma;
    Quark_Task_Flags task_flags = Quark_Task_Flags_Initializer;
    double *scale;

    int i, j, k;
    int SMLSIZ;
    int subpbs, spm1, submat, matsiz;
    int *subpbs_info, *subpbs_info2, *testiwork;
    int INDXQ;
    int K_sub1, K_sub2;
    int iwork_pos, work_pos;
    int msd2;
    int *INDXQ_ptr;
    double *beta_t;

    plasma = plasma_context_self();
    if (sequence->status != PLASMA_SUCCESS)
        return;
    QUARK_Task_Flag_Set(&task_flags, TASK_SEQUENCE, (intptr_t)sequence->quark_sequence);


    /* First steps of dlaed0: norms computation and scaling */
    scale = (double*)malloc(sizeof(double));
    QUARK_CORE_dlaed0_lascl(plasma->quark, &task_flags,
                            n, scale, D, E);

    /*
     * Determine the size and placement of the submatrices, and save in the
     * leading elements of IWORK.
     */
    SMLSIZ = plasma->ev_smlsze;
    IWORK[0] = n;
    subpbs = 1;
    while (IWORK[subpbs-1] > SMLSIZ){
        for (i=subpbs-1; i>-1; i--){
            IWORK[2*i+1] = (IWORK[i]+1)/2;
            IWORK[2*i  ] =  IWORK[i]   /2;
        }
        subpbs *= 2;
    }

    testiwork   = malloc(subpbs * sizeof(int));
    subpbs_info = malloc(subpbs * sizeof(int));
    memcpy(subpbs_info, IWORK, subpbs * sizeof(int));

    testiwork[0] = 1;
    for (i=1; i<subpbs; i++){
        subpbs_info[i] += subpbs_info[i-1];
        testiwork[i] = i+1;
    }

    /*
     * Divide the matrix into SUBPBS submatrices of size at most SMLSIZ+1 using
     * rank-1 modifications (cuts).
     *
     * As it will be modifed by the next sequential part, subproblems
     * informations is backup and then freed by the kernel
     */
    subpbs_info2 = malloc(subpbs * sizeof(int));
    memcpy(subpbs_info2, subpbs_info, subpbs * sizeof(int));
    QUARK_CORE_dlaed0_betaapprox( plasma->quark, &task_flags,
                                  subpbs-1, subpbs_info2, D, E);

    /* Free subpbs displacements backup */
    QUARK_CORE_free( plasma->quark, &task_flags, subpbs_info2, sizeof(int) );

    INDXQ = 4*n;

#if defined(ENABLE_TIMER) && defined(ENABLE_DEBUG2)
    PLASMA_Double_t timedlaed0=0.0;
    timedlaed0 = PLASMA_Wtime();
#endif

#if defined(ENABLE_DEBUG2)
    QUARK_Barrier(plasma->quark);
    printf("  start dlaed0 %d\n",n);
#endif

    /*
     * Solve each submatrix eigenproblem at the bottom of
     * the divide and conquer tree.
     */
    for (i=0; i<subpbs; i++){
        if (i == 0){
            submat = 0;
            matsiz = subpbs_info[0];
        }
        else{
            submat = subpbs_info[i-1];
            matsiz = subpbs_info[i] - subpbs_info[i-1];
        }

        if (icompq == 2) {
            QUARK_CORE_dstedc_f2(plasma->quark, &task_flags,
                                 PlasmaIvec, matsiz,
                                 D+submat, E+submat,
                                 Q+LDQ*submat+submat, LDQ,
                                 /*  */
                                 localdata+i, sizeof(int), INOUT,
                                  /*
                                   * Fake dependency on D to ensure that all
                                   * supproblems have been split by laed0_beta
                                   */
                                 D, sizeof(double)*n, INPUT);
        } else {
            /* Not implemented for now */
            assert(0);
        }

        k = 0;
        for (j=submat; j<subpbs_info[i]; j++, k++){
            IWORK[INDXQ+j] = k;
        }
    }

    while (subpbs > 1){
        spm1  = subpbs - 1;
        work_pos  = 0;
        iwork_pos = 0;
        for (i=0; i<spm1; i+=2){
            if (i == 0){
                K_sub1 = 0;
                K_sub2 = testiwork[1]/2;
                submat = 0;
                matsiz = subpbs_info[1];
                msd2 = subpbs_info[0];
            }
            else{
                K_sub1 = testiwork[i-1];
                K_sub2 = K_sub1 + (testiwork[i+1]-testiwork[i-1])/2;
                submat = subpbs_info[i-1];
                matsiz = subpbs_info[i+1] - subpbs_info[i-1];
                msd2 = matsiz/2;
            }

            INDXQ_ptr = IWORK+INDXQ+submat;
            beta_t = E+submat+msd2-1;

            int last_merge = 0;
            if (subpbs == 2){
                last_merge = 1;
            }

            plasma_pdlaed1_quark(range,
                                 matsiz, msd2,
                                 D+submat,
                                 il, iu,
                                 vl, vu,
                                 Q+LDQ*submat+submat, LDQ,
                                 INDXQ_ptr, beta_t,
                                 WORK+LDWORK*submat+submat, WORK2+work_pos,
                                 IWORK+iwork_pos,  //4*submat
                                 localdata+K_sub1, localdata+K_sub2,
                                 last_merge,
                                 sequence, request);
            work_pos  += 3*matsiz;
            iwork_pos += 4*matsiz;
            subpbs_info[i/2] = subpbs_info[i+1];
            testiwork[i/2] = testiwork[i+1];
        } // loop over the subproblems to be merged
        subpbs = subpbs/2;
#if defined(ENABLE_DEBUG2)
        QUARK_Barrier(plasma->quark);
#endif
    } // while over the level

    /* Scale back the eigenvalues */
    /* Rk: We can use scale because previous dependencies guaranty that scale is initialized */
    QUARK_CORE_dlascl_p2f1( plasma->quark, &task_flags,
                            PlasmaGeneral, 0, 0, &done, scale,
                            n, 1, D, n,
                            /* Fake data to ensure dependencies on the last merge for D */
                            localdata, sizeof(int), INOUT );

    /* Free scale */
    QUARK_CORE_free( plasma->quark, &task_flags, scale, sizeof(double) );

#if defined(ENABLE_TIMER) && defined(ENABLE_DEBUG2)
    timedlaed0   = PLASMA_Wtime()-timedlaed0;
    printf("  Finish dlaed0                            timing= %lf \n",timedlaed0);
#endif

    free(testiwork);
    free(subpbs_info);
}

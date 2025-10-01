/**
 *
 * @file psgetrf_tntpiv.c
 *
 *  PLASMA auxiliary routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 *  LU decomposition with tournament pivoting.
 *  Two algorithms are avialable to select the pivots at each step:
 *    - LU decomposition
 *    - rank revealing QR. (A strong revealing QR would improve the stability,
 *      but is not yet available in LAPACK)
 *
 * @version 2.8.0
 * @author Mathieu Faverge
 * @author Ichitaro Yamazaki
 * @date 2013-02-01
 *
 * @generated s Fri Apr  1 11:03:00 2016
 *
 **/
#include "common.h"

#define A(m,n)    (BLKADDR(A, float, m, n))
#define Acpy(__m) (BLKADDR(W, float, (__m), 0))
#define IPIV(k)    (&(IPIV[(int64_t)A.mb*(int64_t)(k)]))
#define RANK(__m, __k) (&(Wi[(int64_t)W.mb*(int64_t)(__m) + (int64_t)W.lm*(__k)]))

void plasma_psgetrf_rectil_panel_quark(plasma_context_t *plasma,
                                       int *panel_thread_count,
                                       CORE_sgetrf_data_t *data, PLASMA_desc A, int *IPIV,
                                       Quark_Task_Flags *task_flags,
                                       PLASMA_sequence *sequence, PLASMA_request *request)
{
    while ( ((*panel_thread_count * 4 * A.mb) > A.m)
            && (*panel_thread_count > 1) ) {
        (*panel_thread_count)--;
        QUARK_Task_Flag_Set(task_flags, TASK_THREAD_COUNT, *panel_thread_count );
    }

    QUARK_CORE_sgetrf_rectil(
        plasma->quark, task_flags,
        data, A, A(0, 0), A.mb*A.nb, IPIV,
        sequence, request, 1, A.i,
        *panel_thread_count );
}

void plasma_psgetrf_tntpiv_panel_quark(plasma_context_t *plasma,
                                       PLASMA_desc A, int *IPIV,
                                       PLASMA_desc W, int *Wi,
                                       Quark_Task_Flags *task_flags,
                                       PLASMA_sequence *sequence, PLASMA_request *request)
{
    int tempkm, tempmm, tempr;
    int tempm, nexti, pos;
    int ldak, ldam;
    int round, round_size = PLASMA_TNT_SIZE;
    int prev_round_size = 1;
    int curr_round_size = round_size;
    int next_round_size = curr_round_size * round_size;
    int m, i;

    tempkm = min(A.m, A.mb);
    ldak = BLKLDD(A, 0);

    /* Create a first copy of the panel */
    for (m = 0; m < A.mt; m++) {
        tempmm = m == A.mt-1 ? A.m - m*A.mb : A.mb;
        ldam = BLKLDD(A, m);
        i   = m / round_size;
        pos = m % round_size;

        QUARK_CORE_slacpy_f1(
            plasma->quark, task_flags,
            PlasmaUpperLower,
            tempmm, A.n, A.mb,
            A(m, 0), ldam,
            Acpy(i) + pos*A.mb, W.mb,
            Acpy(i), W.mb*W.nb, OUTPUT | GATHERV );
    }

    round = 0;
    while ( ((A.mt - 1) / curr_round_size) > 0 ) {
        /* Let's submit all the factorizations */
        for (m=0, i=0; m<A.mt; m+=curr_round_size, i++) {

            if ( (m+curr_round_size) < A.mt ) {
                tempm = round_size * A.mb;
                tempr = curr_round_size * A.mb;
            } else {
                tempm = ((A.mt-1-m) / prev_round_size) * A.mb;
                if ( ( (A.mt-1-m) % prev_round_size == 0 )
                     && (A.m % A.mb != 0) ) {
                    tempm += A.m % A.mb;
                } else {
                    tempm += A.mb;
                }
                tempr = A.m - m * A.mb;
            }

            QUARK_CORE_sgetrf(
                plasma->quark, task_flags,
                tempm, A.n, A.mb,
                Acpy( i ), W.mb,
                IPIV( i ),
                sequence, request,
                0, A.i );

            tempm = min(tempm, A.n);
            nexti = i / round_size;
            pos   = i % round_size;

#ifdef RESET_DATA_TO_ZERO
            if ( pos == 0 ) {
                assert( nexti <= i );
                QUARK_CORE_slaset(
                    plasma->quark, task_flags,
                    PlasmaUpperLower, W.mb, W.nb, 0., 0., Acpy( nexti ), W.mb );
            }
#endif
            QUARK_CORE_slacpy_pivot(
                plasma->quark, task_flags,
                plasma_desc_submatrix(A, m*A.mb, 0, tempr, A.n),
                PlasmaRowwise, 1, tempm, IPIV( i ),
                RANK( i, round ), RANK( nexti, round+1 ),
                Acpy( nexti ), W.mb,
                pos*A.mb, prev_round_size==1 );
        }

        round++;
        prev_round_size = curr_round_size;
        curr_round_size = next_round_size;
        next_round_size = curr_round_size * round_size;
    }

    /* Last factorization */
    tempm = ((A.mt-1) / prev_round_size) * A.mb;
    if ( ( (A.mt-1) % prev_round_size == 0 )
         && (A.m % A.mb != 0) )
        tempm += A.m % A.mb;
    else
        tempm += A.mb;

    QUARK_CORE_sgetrf(
        plasma->quark, task_flags,
        tempm, A.n, A.mb,
        Acpy(0), W.mb,
        IPIV(0),
        sequence, request,
        1, A.i );

    QUARK_CORE_pivot_update(
        plasma->quark, task_flags,
        tempm, min(tempm, A.n),
        IPIV( 0 ), RANK( 0, round ),
        A.i, (int)(prev_round_size == 1));

    /* Finish to factorize the panel */
    QUARK_CORE_slaswp_ontile(
        plasma->quark, task_flags,
        A, A(0, 0), 1, min(tempkm, A.n), IPIV(0), 1, A(0, 0) );

    /* Copy back the factorization result, once it has been swapped */
    QUARK_CORE_slacpy(
        plasma->quark, task_flags,
        PlasmaUpperLower,
        tempkm, A.n, A.mb,
        Acpy(0), W.mb,
        A(0, 0), ldak);

    /* Apply TRSM on the panel
     * Using A(k,k) ensures that the panel is swapped */
    for (m=1; m<A.mt; m++) {
        tempmm = m == A.mt-1 ? A.m-m*A.mb : A.mb;
        ldam = BLKLDD(A, m);

        QUARK_CORE_strsm(
            plasma->quark, task_flags,
            PlasmaRight, PlasmaUpper, PlasmaNoTrans, PlasmaNonUnit,
            tempmm, A.n, A.mb,
            1., A(0, 0), ldak,
                A(m, 0), ldam);
    }
}

void plasma_psgetrf_tntpiv_qrrr_panel_quark(plasma_context_t *plasma,
                                            PLASMA_desc A, int *IPIV,
                                            PLASMA_desc W, int *Wi,
                                            Quark_Task_Flags *task_flags,
                                            PLASMA_sequence *sequence, PLASMA_request *request)
{
    int tempkm, tempmm, tempr;
    int tempm, nexti, pos;
    int ldak, ldam;
    int round, round_size = PLASMA_TNT_SIZE;
    int prev_round_size = 1;
    int curr_round_size = round_size;
    int next_round_size = curr_round_size * round_size;
    int m, i, copy;

    tempkm = min(A.m, A.mb);
    ldak = BLKLDD(A, 0);

    if ( ((A.mt - 1) / curr_round_size) > 0 )
        copy = PlasmaTrans;
    else
        copy = PlasmaNoTrans;

    /* Create a first copy of the panel */
    for (m = 0; m < A.mt; m++) {
        tempmm = m == A.mt-1 ? A.m - m*A.mb : A.mb;
        ldam = BLKLDD(A, m);
        i   = m / round_size;
        pos = m % round_size;

        if (copy == PlasmaNoTrans) {
            QUARK_CORE_slatro_f1(
                plasma->quark, task_flags,
                PlasmaUpperLower, copy,
                tempmm, A.n, A.mb,
                A(m, 0), ldam,
                Acpy(i) + pos*A.mb, W.mb,
                Acpy(i), W.mb*W.nb, OUTPUT | GATHERV );
        } else {
            QUARK_CORE_slatro_f1(
                plasma->quark, task_flags,
                PlasmaUpperLower, copy,
                tempmm, A.n, A.mb,
                A(m, 0), ldam,
                Acpy(i) + pos*(A.mb*A.nb), W.nb,
                Acpy(i), W.mb*W.nb, OUTPUT | GATHERV );
        }
    }

    round = 0;
    while ( ((A.mt - 1) / curr_round_size) > 0 ) {

        /* Let's submit all the factorizations */
        for (m=0, i=0; m<A.mt; m+=curr_round_size, i++) {

            if ( (m+curr_round_size) < A.mt ) {
                tempm = round_size * A.mb;
                tempr = curr_round_size * A.mb;
            } else {
                tempm = ((A.mt-1-m) / prev_round_size) * A.mb;
                if ( ( (A.mt-1-m) % prev_round_size == 0 )
                     && (A.m % A.mb != 0) ) {
                    tempm += A.m % A.mb;
                } else {
                    tempm += A.mb;
                }
                tempr = A.m - m * A.mb;
            }

            QUARK_CORE_sgeqp3_tntpiv(
                plasma->quark, task_flags,
                A.n, tempm, A.mb,
                Acpy( i ), W.nb,
                IPIV( i ),
                sequence, request,
                0, A.i );

            tempm = min(tempm, A.n);
            nexti = i / round_size;
            pos   = i % round_size;

            if ( ((A.mt - 1) / next_round_size) > 0 )
                copy = PlasmaTrans;
            else
                copy = PlasmaNoTrans;

            /* If it's the last round, we don't transpose */
            if ( copy == PlasmaNoTrans ) {
                QUARK_CORE_slacpy_pivot(
                    plasma->quark, task_flags,
                    plasma_desc_submatrix(A, m*A.mb, 0, tempr, A.n),
                    PlasmaRowwise, 1, tempm, IPIV( i ),
                    RANK( i, round ), RANK( nexti, round+1 ),
                    Acpy( nexti ), W.mb,
                    pos*A.mb, prev_round_size==1 );
            } else {
                QUARK_CORE_slacpy_pivot(
                    plasma->quark, task_flags,
                    plasma_desc_submatrix(A, m*A.mb, 0, tempr, A.n),
                    PlasmaColumnwise, 1, tempm, IPIV( i ),
                    RANK( i, round ), RANK( nexti, round+1 ),
                    Acpy( nexti ), W.nb,
                    pos*A.mb, prev_round_size==1 );
            }
        }

        round++;
        prev_round_size = curr_round_size;
        curr_round_size = next_round_size;
        next_round_size = curr_round_size * round_size;
    }

    /* Last factorization */
    tempm = ((A.mt-1) / prev_round_size) * A.mb;
    if ( ( (A.mt-1) % prev_round_size == 0 )
         && (A.m % A.mb != 0) )
        tempm += A.m % A.mb;
    else
        tempm += A.mb;

    QUARK_CORE_sgetrf(
        plasma->quark, task_flags,
        tempm, A.n, A.mb,
        Acpy(0), W.mb,
        IPIV(0),
        sequence, request,
        1, A.i );

    QUARK_CORE_pivot_update(
        plasma->quark, task_flags,
        tempm, min(tempm, A.n),
        IPIV( 0 ), RANK( 0, round ),
        A.i, (int)(prev_round_size == 1));

    /* Finish to factorize the panel */
    QUARK_CORE_slaswp_ontile(
        plasma->quark, task_flags,
        A, A(0, 0), 1, min(tempkm, A.n), IPIV(0), 1, A(0, 0) );

    /* Copy back the factorization result, once it has been swapped */
    QUARK_CORE_slacpy(
        plasma->quark, task_flags,
        PlasmaUpperLower,
        tempkm, A.n, A.mb,
        Acpy(0), W.mb,
        A(0, 0), ldak);

    /* Apply TRSM on the panel
     * Using A(k,k) ensures that the panel is swapped */
    for (m=1; m<A.mt; m++) {
        tempmm = m == A.mt-1 ? A.m-m*A.mb : A.mb;
        ldam = BLKLDD(A, m);

        QUARK_CORE_strsm(
            plasma->quark, task_flags,
            PlasmaRight, PlasmaUpper, PlasmaNoTrans, PlasmaNonUnit,
            tempmm, A.n, A.mb,
            1., A(0, 0), ldak,
                A(m, 0), ldam);
    }
}
/*
 * W is a workspace in tile layout with tile of size max_round_size*A.mb -by- A.nb
 * Wi is an integer workspace to store the rank of the lines involved in each round.
 * This workspace has tiles of size max_round_size*A.mb -by- 1
 */

/***************************************************************************//**
 *  Parallel tile LU factorization - dynamic scheduling - Right looking
 **/
void plasma_psgetrf_tntpiv_quark(PLASMA_desc A, int *IPIV,
                                 PLASMA_desc W, int *Wi,
                                 PLASMA_sequence *sequence, PLASMA_request *request)
{
    int k, m, n, minmnt;
    plasma_context_t *plasma;
    int tempkm, tempkn, tempmm, tempnn;
    int tempm, tempk;
    int ldak, ldam;
    Quark_Task_Flags task_flags = Quark_Task_Flags_Initializer;

    float zone  = (float)1.0;
    float mzone = (float)-1.0;
    void * fakedep;

    void (*panel)(plasma_context_t *, PLASMA_desc, int *,
                  PLASMA_desc, int *, Quark_Task_Flags *,
                  PLASMA_sequence *, PLASMA_request *) = NULL;

    plasma = plasma_context_self();
    if (sequence->status != PLASMA_SUCCESS)
        return;
    QUARK_Task_Flag_Set(&task_flags, TASK_SEQUENCE, (intptr_t)sequence->quark_sequence);
    plasma_profile_by_function( &task_flags, GETRF );

    if (PLASMA_TNT_MODE == PLASMA_TOURNAMENT_LU)
        panel = plasma_psgetrf_tntpiv_panel_quark;
    else
        panel = plasma_psgetrf_tntpiv_qrrr_panel_quark;

    minmnt = min(A.mt, A.nt);
    for (k = 0; k < minmnt; k++)
    {
        tempk  = k * A.mb;
        tempm  = A.m - tempk;
        tempkm = k == A.mt-1 ? A.m-k*A.mb : A.mb;
        tempkn = k == A.nt-1 ? A.n-k*A.nb : A.nb;
        ldak = BLKLDD(A, k);

        QUARK_Task_Flag_Set(&task_flags, TASK_PRIORITY, QUARK_TASK_MAX_PRIORITY - k );

        /*
         * Factorize the panel
         */
        panel(plasma,
              plasma_desc_submatrix(A, tempk, tempk, tempm, tempkn),
              IPIV(k), W, Wi, &task_flags,
              sequence, request);

        /*
         * Update the trailing submatrix
         */
        fakedep = (void *)(intptr_t)(k+1);
        for (n = k+1; n < A.nt; n++)
        {
            QUARK_Task_Flag_Set(&task_flags, TASK_PRIORITY, QUARK_TASK_MAX_PRIORITY - n );
            /*
             * Apply row interchange after the panel (work on the panel)
             */
            tempnn = n == A.nt-1 ? A.n-n*A.nb : A.nb;
            QUARK_CORE_sswptr_ontile(
                plasma->quark, &task_flags,
                plasma_desc_submatrix(A, tempk, n*A.nb, tempm, tempnn),
                A(k, n), 1, tempkm, IPIV(k), 1,
                A(k, k), ldak);

            m = k+1;
            if ( m < A.mt ) {
                tempmm = m == A.mt-1 ? A.m-m*A.mb : A.mb;
                ldam = BLKLDD(A, m);

                QUARK_CORE_sgemm2(
                    plasma->quark, &task_flags,
                    PlasmaNoTrans, PlasmaNoTrans,
                    tempmm, tempnn, A.nb, A.mb,
                    mzone, A(m, k), ldam,
                           A(k, n), ldak,
                    zone,  A(m, n), ldam);

                for (m = k+2; m < A.mt; m++)
                {
                    tempmm = m == A.mt-1 ? A.m-m*A.mb : A.mb;
                    ldam = BLKLDD(A, m);

                    QUARK_CORE_sgemm_f2(
                        plasma->quark, &task_flags,
                        PlasmaNoTrans, PlasmaNoTrans,
                        tempmm, tempnn, A.nb, A.mb,
                        mzone, A(m, k), ldam,
                               A(k, n), ldak,
                        zone,  A(m, n), ldam,
                        /* Dependency on next swapa (gemm need to be done before) */
                        A(k+1, n), A.mb*A.nb, INOUT | GATHERV,
                        /* Dependency on next swapb (gemm need to use panel k before it has to be swaped */
                        fakedep,   1,         INPUT );
                }
            }
        }
    }

    QUARK_Task_Flag_Set(&task_flags, TASK_PRIORITY, QUARK_TASK_MIN_PRIORITY );
    for (k = 0; k < min(A.mt, A.nt); k++)
    {
        int mintmp;
        tempk  = k * A.mb;
        tempm  = A.m - tempk;
        tempkm = k == A.mt-1 ? tempm : A.mb;
        tempkn = k == A.nt-1 ? A.n - k * A.nb : A.nb;
        mintmp = min(tempkm, tempkn);
        ldak = BLKLDD(A, k);

        /*
         * Apply row interchange behind the panel (work on the panel)
         */
        fakedep = (void*)(intptr_t)k;
        for (n = 0; n < k; n++)
        {
            tempnn = n == A.nt-1 ? A.n-n*A.nb : A.nb;
            QUARK_CORE_slaswp_ontile_f2(
                plasma->quark, &task_flags,
                plasma_desc_submatrix(A, tempk, n*A.nb, tempm, tempnn),
                A(k, n), 1, mintmp, IPIV(k), 1,
                /* Dependency on previous swapb */
                A(k-1,n), A.lm*A.nb, INPUT,
                /* Dependency on all GEMM from previous step */
                fakedep,  1,         INOUT | GATHERV );
        }
    }
}

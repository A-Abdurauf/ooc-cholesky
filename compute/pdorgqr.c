/**
 *
 * @file pdorgqr.c
 *
 *  PLASMA auxiliary routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Hatem Ltaief
 * @author Jakub Kurzak
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @generated d Fri Apr  1 11:02:58 2016
 *
 **/
#include "common.h"

#define A(m,n) BLKADDR(A, double, m, n)
#define Q(m,n) BLKADDR(Q, double, m, n)
#define T(m,n) BLKADDR(T, double, m, n)
/***************************************************************************//**
 *  Parallel construction of Q using tile V (application to identity) - dynamic scheduling
 **/
void plasma_pdorgqr_quark(PLASMA_desc A, PLASMA_desc Q, PLASMA_desc T,
                          PLASMA_sequence *sequence, PLASMA_request *request)
{
    plasma_context_t *plasma;
    Quark_Task_Flags task_flags = Quark_Task_Flags_Initializer;

    int k, m, n;
    int ldak, ldqk, ldam, ldqm;
    int tempmm, tempnn, tempkmin, tempkm;
    int tempAkm, tempAkn;
    int ib, minmnt;

    plasma = plasma_context_self();
    if (sequence->status != PLASMA_SUCCESS)
        return;
    QUARK_Task_Flag_Set(&task_flags, TASK_SEQUENCE, (intptr_t)sequence->quark_sequence);
    plasma_profile_by_function( &task_flags, UNMQR );

    ib = PLASMA_IB;
    minmnt = min(A.mt, A.nt);
    for (k = minmnt-1; k >= 0; k--) {
        tempAkm  = k == A.mt-1 ? A.m-k*A.mb : A.mb;
        tempAkn  = k == A.nt-1 ? A.n-k*A.nb : A.nb;
        tempkmin = min( tempAkn, tempAkm );
        tempkm   = k == Q.mt-1 ? Q.m-k*Q.mb : Q.mb;
        ldak = BLKLDD(A, k);
        ldqk = BLKLDD(Q, k);
        for (m = Q.mt - 1; m > k; m--) {
            tempmm = m == Q.mt-1 ? Q.m-m*Q.mb : Q.mb;
            ldam = BLKLDD(A, m);
            ldqm = BLKLDD(Q, m);
            for (n = k; n < Q.nt; n++) {
                tempnn = n == Q.nt-1 ? Q.n-n*Q.nb : Q.nb;
                QUARK_CORE_dtsmqr(
                    plasma->quark, &task_flags,
                    PlasmaLeft, PlasmaNoTrans,
                    Q.mb, tempnn, tempmm, tempnn, tempAkn, ib, T.nb,
                    Q(k, n), ldqk,
                    Q(m, n), ldqm,
                    A(m, k), ldam,
                    T(m, k), T.mb);
            }
        }
        for (n = k; n < Q.nt; n++) {
            tempnn = n == Q.nt-1 ? Q.n-n*Q.nb : Q.nb;
            QUARK_CORE_dormqr(
                plasma->quark, &task_flags,
                PlasmaLeft, PlasmaNoTrans,
                tempkm, tempnn, tempkmin, ib, T.nb,
                A(k, k), ldak,
                T(k, k), T.mb,
                Q(k, n), ldqk);
        }
    }
}

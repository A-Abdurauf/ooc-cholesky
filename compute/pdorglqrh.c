/**
 *
 * @file pdorglqrh.c
 *
 *  PLASMA auxiliary routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Dulceneia Becker
 * @date 2011-05-24
 * @generated d Fri Apr  1 11:02:58 2016
 *
 **/
#include "common.h"

#define A(m,n)  BLKADDR(A, double, (m), (n))
#define Q(m,n)  BLKADDR(Q, double, (m), (n))
#define T(m,n)  BLKADDR(T, double, (m), (n))
#define T2(m,n) BLKADDR(T, double, (m), (n)+(A.nt))
/***************************************************************************//**
 *  Parallel construction of Q using tile V (application to identity;
 *  reduction Householder) - dynamic scheduling
 **/
void plasma_pdorglqrh_quark(PLASMA_desc A, PLASMA_desc Q,
        PLASMA_desc T, int BS,
        PLASMA_sequence *sequence, PLASMA_request *request)
{
    plasma_context_t *plasma;
    Quark_Task_Flags task_flags = Quark_Task_Flags_Initializer;

    int k, m, n;
    int K, N, RD, lastRD;
    int ldak;
    int ldqm;
    int tempkm, tempkmin, tempNn, tempnn, tempmm, tempNRDn;
    int ib;

    plasma = plasma_context_self();
    if (sequence->status != PLASMA_SUCCESS)
        return;
    QUARK_Task_Flag_Set(&task_flags, TASK_SEQUENCE, (intptr_t)sequence->quark_sequence);
    plasma_profile_by_function( &task_flags, UNMLQ );

    ib = PLASMA_IB;
    K = min(A.mt, A.nt);

    for (k = K-1; k >= 0; k--) {
        tempkm = k == A.mt-1 ? A.m-k*A.mb : A.mb;
        ldak = BLKLDD(A, k);
        lastRD = 0;
        for (RD = BS; RD < A.nt-k; RD *= 2)
            lastRD = RD;
        for (RD = lastRD; RD >= BS; RD /= 2) {
            for (N = k; N+RD < A.nt; N += 2*RD) {
                tempNRDn = N+RD == A.nt-1 ? A.n-(N+RD)*A.nb : A.nb;
                for (m = k; m < Q.mt; m++) {
                    tempmm = m == Q.mt-1 ? Q.m-m*Q.mb : Q.mb;
                    ldqm   = BLKLDD(Q, m   );
                    QUARK_CORE_dttmlq(
                                      plasma->quark, &task_flags,
                                      PlasmaRight, PlasmaNoTrans,
                                      tempmm, Q.nb, tempmm, tempNRDn,
                                      tempkm, ib, T.nb,
                                      Q (m, N   ), ldqm,
                                      Q (m, N+RD), ldqm,
                                      A (k, N+RD), ldak,
                                      T2(k, N+RD), T.mb);
                }
            }
        }
        for (N = k; N < A.nt; N += BS) {
            tempNn = N == A.nt-1 ? A.n-N*A.nb : A.nb;
            tempkmin = min(tempkm, tempNn);
            for (n = min(N+BS, A.nt)-1; n > N; n--) {
                tempnn = n == Q.nt-1 ? Q.n-n*Q.nb : Q.nb;

                for (m = k; m < Q.mt; m++) {
                    tempmm = m == Q.mt-1 ? Q.m-m*Q.mb : Q.mb;
                    ldqm = BLKLDD(Q, m);
                    QUARK_CORE_dtsmlq(
                                      plasma->quark, &task_flags,
                                      PlasmaRight, PlasmaNoTrans,
                                      tempmm, Q.nb, tempmm, tempnn,
                                      tempkm, ib, T.nb,
                                      Q(m, N), ldqm,
                                      Q(m, n), ldqm,
                                      A(k, n), ldak,
                                      T(k, n), T.mb);
                }
            }
            for (m = k; m < Q.mt; m++) {
                tempmm = m == Q.mt-1 ? Q.m-m*Q.mb : Q.mb;
                ldqm = BLKLDD(Q, m);
                QUARK_CORE_dormlq(
                                  plasma->quark, &task_flags,
                                  PlasmaRight, PlasmaNoTrans,
                                  tempmm, tempNn,
                                  tempkmin, ib, T.nb,
                                  A(k, N), ldak,
                                  T(k, N), T.mb,
                                  Q(m, N), ldqm);
            }
        }
    }
}

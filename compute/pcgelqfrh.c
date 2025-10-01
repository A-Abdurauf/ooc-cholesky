/**
 *
 * @file pcgelqfrh.c
 *
 *  PLASMA auxiliary routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Jakub Kurzak
 * @author Hatem Ltaief
 * @author Dulceneia Becker
 * @date 2010-11-15
 * @generated c Fri Apr  1 11:02:57 2016
 *
 **/
#include "common.h"

#define A(m,n)  BLKADDR(A, PLASMA_Complex32_t, (m), (n))
#define T(m,n)  BLKADDR(T, PLASMA_Complex32_t, (m), (n))
#define T2(m,n) BLKADDR(T, PLASMA_Complex32_t, (m), (n)+A.nt)
/***************************************************************************//**
 *  Parallel tile LQ factorization (reduction Householder) - dynamic scheduling
 **/
void plasma_pcgelqfrh_quark(PLASMA_desc A, PLASMA_desc T, int BS,
                            PLASMA_sequence *sequence, PLASMA_request *request)
{
    plasma_context_t *plasma;
    Quark_Task_Flags task_flags = Quark_Task_Flags_Initializer;

    int k, m, n;
    int N, RD;
    int ldak, ldam;
    int tempkmin, tempkm, tempNn, tempnn, tempmm, tempNRDn;
    int ib;

    plasma = plasma_context_self();
    if (sequence->status != PLASMA_SUCCESS)
        return;
    QUARK_Task_Flag_Set(&task_flags, TASK_SEQUENCE, (intptr_t)sequence->quark_sequence);
    plasma_profile_by_function( &task_flags, GELQT );

    ib = PLASMA_IB;

    for (k = 0; k < min(A.mt, A.nt); k++) {
        tempkm = k == A.mt-1 ? A.m-k*A.mb : A.mb;
        ldak = BLKLDD(A, k);
        for (N = k; N < A.nt; N += BS) {
            tempNn = N == A.nt-1 ? A.n-N*A.nb : A.nb;
            tempkmin = min(tempkm, tempNn);
            QUARK_CORE_cgelqt(
                plasma->quark, &task_flags,
                tempkm, tempNn, ib, T.nb,
                A(k, N), ldak,
                T(k, N), T.mb);
            for (m = k+1; m < A.mt; m++) {
                tempmm = m == A.mt-1 ? A.m-m*A.mb : A.mb;
                ldam = BLKLDD(A, m);
                QUARK_CORE_cunmlq(
                    plasma->quark, &task_flags,
                    PlasmaRight, PlasmaConjTrans,
                    tempmm, tempNn, tempkmin, ib, T.nb,
                    A(k, N), ldak,
                    T(k, N), T.mb,
                    A(m, N), ldam);
            }
            for (n = N+1; n < min(N+BS, A.nt); n++) {
                tempnn = n == A.nt-1 ? A.n-n*A.nb : A.nb;
                QUARK_CORE_ctslqt(
                    plasma->quark, &task_flags,
                    tempkm, tempnn, ib, T.nb,
                    A(k, N), ldak,
                    A(k, n), ldak,
                    T(k, n), T.mb);

                for (m = k+1; m < A.mt; m++) {
                    tempmm = m == A.mt-1 ? A.m-m*A.mb : A.mb;
                    ldam = BLKLDD(A, m);
                    QUARK_CORE_ctsmlq(
                        plasma->quark, &task_flags,
                        PlasmaRight, PlasmaConjTrans,
                        tempmm, A.nb, tempmm, tempnn, tempkm, ib, T.nb,
                        A(m, N), ldam,
                        A(m, n), ldam,
                        A(k, n), ldak,
                        T(k, n), T.mb);
                }
            }
        }
        for (RD = BS; RD < A.nt-k; RD *= 2) {
            for (N = k; N+RD < A.nt; N += 2*RD) {
                tempNRDn = N+RD == A.nt-1 ? A.n-(N+RD)*A.nb : A.nb;
                QUARK_CORE_cttlqt(
                    plasma->quark, &task_flags,
                    tempkm, tempNRDn, ib, T.nb,
                    A (k, N   ), ldak,
                    A (k, N+RD), ldak,
                    T2(k, N+RD), T.mb);

                for (m = k+1; m < A.mt; m++) {
                    tempmm = m == A.mt-1 ? A.m-m*A.mb : A.mb;
                    ldam   = BLKLDD(A, m );
                    QUARK_CORE_cttmlq(
                        plasma->quark, &task_flags,
                        PlasmaRight, PlasmaConjTrans,
                        tempmm, A.nb, tempmm, tempNRDn, tempkm, ib, T.nb,
                        A (m, N   ), ldam,
                        A (m, N+RD), ldam,
                        A (k, N+RD), ldak,
                        T2(k, N+RD), T.mb);
                }
            }
        }
    }
}

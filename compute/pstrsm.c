/**
 *
 * @file pstrsm.c
 *
 *  PLASMA auxiliary routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Jakub Kurzak
 * @author Hatem Ltaief
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @generated s Fri Apr  1 11:02:57 2016
 *
 **/
#include "common.h"

#define A(m,n) BLKADDR(A, float, m, n)
#define B(m,n) BLKADDR(B, float, m, n)
/***************************************************************************//**
 *  Parallel tile triangular solve - static scheduling
 **/
void plasma_pstrsm(plasma_context_t *plasma)
{
    PLASMA_enum side;
    PLASMA_enum uplo;
    PLASMA_enum trans;
    PLASMA_enum diag;
    float alpha;
    PLASMA_desc A;
    PLASMA_desc B;
    PLASMA_sequence *sequence;
    PLASMA_request *request;

    int k, m, n;
    int next_k;
    int next_m;
    int next_n;
    int ldam, ldan, ldak, ldbm, ldbk;
    int tempkm, tempkn, tempmm, tempnn;

    float zone  = (float) 1.0;
    float mzone = (float)-1.0;
    float lalpha;
    float minvalpha;

    plasma_unpack_args_9(side, uplo, trans, diag, alpha, A, B, sequence, request);
    minvalpha = mzone / alpha;
    if (sequence->status != PLASMA_SUCCESS)
        return;
    ss_init(B.mt, B.nt, -1);
    /*
     *  PlasmaLeft
     */
    if (side == PlasmaLeft) {
        k = 0;
        m = PLASMA_RANK;
        while (m >= B.mt) {
            k++;
            m = m - B.mt + k;
        }
        n = 0;

        while (k < B.mt && m < B.mt) {
            next_k = k;
            next_m = m;
            next_n = n;

            next_n++;
            if (next_n >= B.nt) {
                next_m += PLASMA_SIZE;
                while (next_m >= B.mt && next_k < B.mt) {
                    next_k++;
                    next_m = next_m - B.mt + next_k;
                }
                next_n = 0;
            }

            tempnn = n == B.nt-1 ? B.n-n*B.nb : B.nb;
            tempmm = m == B.mt-1 ? B.m-m*B.mb : B.mb;

            lalpha = k == 0 ? alpha : zone;
            if (m == k) {
                ss_cond_wait(m, n, k-1);
                /*
                 *  PlasmaLeft / PlasmaLower / PlasmaNoTrans
                 *  PlasmaLeft / PlasmaUpper / Plasma[Conj]Trans
                 */
                if ((uplo == PlasmaLower && trans == PlasmaNoTrans)
                    || (uplo == PlasmaUpper && trans != PlasmaNoTrans)) {
                    tempkm = k == B.mt-1 ? B.m-k*B.mb : B.mb;
                    ldak = BLKLDD(A, k);
                    ldbk = BLKLDD(B, k);
                    CORE_strsm(
                        side, uplo, trans, diag,
                        tempkm, tempnn,
                        lalpha, A(k, k), ldak,
                                B(k, n), ldbk);
                }
                /*
                 *  PlasmaLeft / PlasmaLower / Plasma[Cojn]Trans
                 *  PlasmaLeft / PlasmaUpper / PlasmaNoTrans
                 */
                else {
                    tempkm = k == 0 ? B.m-(B.mt-1)*B.mb : B.mb;
                    ldak = BLKLDD(A, B.mt-1-k);
                    ldbk = BLKLDD(B, B.mt-1-k);
                    CORE_strsm(
                        side, uplo, trans, diag,
                        tempkm, tempnn,
                        lalpha, A(B.mt-1-k, B.mt-1-k), ldak,
                                B(B.mt-1-k, n       ), ldbk);
                }
                ss_cond_set(k, n, k);
            }
            else {
                ss_cond_wait(k, n, k);
                ss_cond_wait(m, n, k-1);
                /*
                 *  PlasmaRight / PlasmaLower / PlasmaNoTrans
                 */
                if (uplo == PlasmaLower) {
                    if (trans == PlasmaNoTrans) {
                        ldam = BLKLDD(A, m);
                        ldbk = BLKLDD(B, k);
                        ldbm = BLKLDD(B, m);
                        CORE_sgemm(
                            PlasmaNoTrans, PlasmaNoTrans,
                            tempmm, tempnn, B.mb,
                            mzone,  A(m, k), ldam,
                                    B(k, n), ldbk,
                            lalpha, B(m, n), ldbm);
                    }
                    /*
                     *  PlasmaRight / PlasmaLower / Plasma[Conj]Trans
                     */
                    else {
                        tempkm = k == 0 ? A.m-(A.mt-1)*A.mb : A.mb;
                        ldak = BLKLDD(A, B.mt-1-k);
                        ldbk = BLKLDD(B, B.mt-1-k);
                        ldbm = BLKLDD(B, B.mt-1-m);
                        CORE_sgemm(
                            trans, PlasmaNoTrans,
                            B.mb, tempnn, tempkm,
                            mzone,  A(A.mt-1-k, A.mt-1-m), ldak,
                                    B(B.mt-1-k, n       ), ldbk,
                            lalpha, B(B.mt-1-m, n       ), ldbm);
                    }
                }
                else {
                    /*
                     *  PlasmaRight / PlasmaUpper / PlasmaNoTrans
                     */
                    if (trans == PlasmaNoTrans) {
                        tempkm = k == 0 ? A.m-(A.mt-1)*A.mb : A.mb;
                        ldam = BLKLDD(A, B.mt-1-m);
                        ldbk = BLKLDD(B, B.mt-1-k);
                        ldbm = BLKLDD(B, B.mt-1-m);
                        CORE_sgemm(
                            PlasmaNoTrans, PlasmaNoTrans,
                            B.mb, tempnn, tempkm,
                            mzone,  A(A.mt-1-m, A.mt-1-k), ldam,
                                    B(B.mt-1-k, n       ), ldbk,
                            lalpha, B(B.mt-1-m, n       ), ldbm);
                    }
                    /*
                     *  PlasmaRight / PlasmaUpper / Plasma[Conj]Trans
                     */
                    else {
                        ldak = BLKLDD(A, k);
                        ldbk = BLKLDD(B, k);
                        ldbm = BLKLDD(B, m);
                        CORE_sgemm(
                            trans, PlasmaNoTrans,
                            tempmm, tempnn, B.mb,
                            mzone,  A(k, m), ldak,
                                    B(k, n), ldbk,
                            lalpha, B(m, n), ldbm);
                    }
                }
                ss_cond_set(m, n, k);
            }
            n = next_n;
            m = next_m;
            k = next_k;
        }
    }
    /*
     *  PlasmaRight
     */
    else {
        k = 0;
        n = PLASMA_RANK;
        while (n >= B.nt) {
            k++;
            n = n - B.nt + k;
        }
        m = 0;

        while (k < B.nt && n < B.nt) {
            next_n = n;
            next_m = m;
            next_k = k;

            next_m++;
            if (next_m >= B.mt) {
                next_n += PLASMA_SIZE;
                while (next_n >= B.nt && next_k < B.nt) {
                    next_k++;
                    next_n = next_n - B.nt + next_k;
                }
                next_m = 0;
            }

            tempmm = m == B.mt-1 ? B.m-m*B.mb : B.mb;
            tempnn = n == B.nt-1 ? B.n-n*B.nb : B.nb;

            lalpha = k == 0 ? alpha : zone;
            if (n == k) {
                ss_cond_wait(m, n, k-1);
                /*
                 *  PlasmaRight / PlasmaLower / PlasmaNoTrans
                 */
                if (uplo == PlasmaLower) {
                    if (trans == PlasmaNoTrans) {
                        tempkn = k == 0 ? B.n-(B.nt-1)*B.nb : B.nb;
                        ldak = BLKLDD(A, B.nt-1-k);
                        ldbm = BLKLDD(B, m);
                        CORE_strsm(
                            side, uplo, trans, diag,
                            tempmm, tempkn,
                            lalpha, A(B.nt-1-k, B.nt-1-k), ldak,
                                    B(m,        B.nt-1-k), ldbm);
                    }
                    /*
                     *  PlasmaRight / PlasmaLower / Plasma[Conj]Trans
                     */
                    else {
                        tempkn = k == B.nt-1 ? B.n-k*B.nb : B.nb;
                        ldak = BLKLDD(A, k);
                        ldbm = BLKLDD(B, m);
                        CORE_strsm(
                            side, uplo, trans, diag,
                            tempmm, tempkn,
                            alpha, A(k, k), ldak,
                                   B(m, k), ldbm);
                    }
                }
                else {
                    /*
                     *  PlasmaRight / PlasmaUpper / PlasmaNoTrans
                     */
                    if (trans == PlasmaNoTrans) {
                        tempkn = k == B.nt-1 ? B.n-k*B.nb : B.nb;
                        ldak = BLKLDD(A, k);
                        ldbm = BLKLDD(B, m);
                        CORE_strsm(
                            side, uplo, trans, diag,
                            tempmm, tempkn,
                            lalpha, A(k, k), ldak,
                                    B(m, k), ldbm);
                    }
                    /*
                     *  PlasmaRight / PlasmaUpper / Plasma[Conj]Trans
                     */
                    else {
                        tempkn = k == 0 ? B.n-(B.nt-1)*B.nb : B.nb;
                        ldak = BLKLDD(A, B.nt-1-k);
                        ldbm = BLKLDD(B, m);
                        CORE_strsm(
                            side, uplo, trans, diag,
                            tempmm, tempkn,
                            alpha, A(B.nt-1-k, B.nt-1-k), ldak,
                                   B(m,        B.nt-1-k), ldbm);
                    }
                }
                ss_cond_set(m, k, k);
            }
            else {
                ss_cond_wait(m, k, k);
                ss_cond_wait(m, n, k-1);
                /*
                 *  PlasmaRight / PlasmaLower / PlasmaNoTrans
                 */
                if (uplo == PlasmaLower) {
                    if (trans == PlasmaNoTrans) {
                        tempkn = k == 0 ? B.n-(B.nt-1)*B.nb : B.nb;
                        ldak = BLKLDD(A, B.nt-1-k);
                        ldbm = BLKLDD(B, m);
                        CORE_sgemm(
                            PlasmaNoTrans, PlasmaNoTrans,
                            tempmm, B.mb, tempkn,
                            mzone,  B(m,        B.nt-1-k), ldbm,
                                    A(B.nt-1-k, B.nt-1-n), ldak,
                            lalpha, B(m,        B.nt-1-n), ldbm);
                    }
                    /*
                     *  PlasmaRight / PlasmaLower / Plasma[Conj]Trans
                     */
                    else {
                        ldan = BLKLDD(A, n);
                        ldbm = BLKLDD(B, m);
                        CORE_sgemm(
                            PlasmaNoTrans, trans,
                            tempmm, tempnn, B.mb,
                            minvalpha, B(m, k), ldbm,
                                       A(n, k), ldan,
                            zone,      B(m, n), ldbm);
                    }
                }
                else {
                    /*
                     *  PlasmaRight / PlasmaUpper / PlasmaNoTrans
                     */
                    if (trans == PlasmaNoTrans) {
                        ldak = BLKLDD(A, k);
                        ldbm = BLKLDD(B, m);
                        CORE_sgemm(
                            PlasmaNoTrans, PlasmaNoTrans,
                            tempmm, tempnn, B.mb,
                            mzone,  B(m, k), ldbm,
                                    A(k, n), ldak,
                            lalpha, B(m, n), ldbm);
                    }
                    /*
                     *  PlasmaRight / PlasmaUpper / Plasma[Conj]Trans
                     */
                    else {
                        tempkn = k == 0 ? B.n-(B.nt-1)*B.nb : B.nb;
                        ldan = BLKLDD(A, n);
                        ldbm = BLKLDD(B, m);
                        CORE_sgemm(
                            PlasmaNoTrans, trans,
                            tempmm, B.nb, tempkn,
                            minvalpha, B(m,        B.nt-1-k), ldbm,
                                       A(B.nt-1-n, B.nt-1-k), ldan,
                            zone,      B(m,        B.nt-1-n), ldbm);
                    }
                }
                ss_cond_set(m, n, k);
            }
            n = next_n;
            m = next_m;
            k = next_k;
        }
    }
    ss_finalize();
}

/***************************************************************************//**
 *  Parallel tile triangular solve - dynamic scheduling
 **/
void plasma_pstrsm_quark(PLASMA_enum side, PLASMA_enum uplo, PLASMA_enum trans, PLASMA_enum diag,
                         float alpha, PLASMA_desc A, PLASMA_desc B,
                         PLASMA_sequence *sequence, PLASMA_request *request)
{
    plasma_context_t *plasma;
    Quark_Task_Flags task_flags = Quark_Task_Flags_Initializer;

    int k, m, n;
    int ldak, ldam, ldan, ldbk, ldbm;
    int tempkm, tempkn, tempmm, tempnn;

    float zone       = (float) 1.0;
    float mzone      = (float)-1.0;
    float minvalpha  = (float)-1.0 / alpha;
    float lalpha;

    plasma = plasma_context_self();
    if (sequence->status != PLASMA_SUCCESS)
        return;
    QUARK_Task_Flag_Set(&task_flags, TASK_SEQUENCE, (intptr_t)sequence->quark_sequence);
    plasma_profile_by_function( &task_flags, TRSM );

    /*
     *  PlasmaLeft / PlasmaUpper / PlasmaNoTrans
     */
    if (side == PlasmaLeft) {
        if (uplo == PlasmaUpper) {
            if (trans == PlasmaNoTrans) {
                for (k = 0; k < B.mt; k++) {
                    tempkm = k == 0 ? B.m-(B.mt-1)*B.mb : B.mb;
                    ldak = BLKLDD(A, B.mt-1-k);
                    ldbk = BLKLDD(B, B.mt-1-k);
                    lalpha = k == 0 ? alpha : zone;
                    for (n = 0; n < B.nt; n++) {
                        tempnn = n == B.nt-1 ? B.n-n*B.nb : B.nb;
                        QUARK_CORE_strsm(
                            plasma->quark, &task_flags,
                            side, uplo, trans, diag,
                            tempkm, tempnn, A.mb,
                            lalpha, A(B.mt-1-k, B.mt-1-k), ldak,  /* lda * tempkm */
                                    B(B.mt-1-k,        n), ldbk); /* ldb * tempnn */
                    }
                    for (m = k+1; m < B.mt; m++) {
                        ldam = BLKLDD(A, B.mt-1-m);
                        ldbm = BLKLDD(B, B.mt-1-m);
                        for (n = 0; n < B.nt; n++) {
                            tempnn = n == B.nt-1 ? B.n-n*B.nb : B.nb;
                            QUARK_CORE_sgemm(
                                plasma->quark, &task_flags,
                                PlasmaNoTrans, PlasmaNoTrans,
                                B.mb, tempnn, tempkm, A.mb,
                                mzone,  A(B.mt-1-m, B.mt-1-k), ldam,
                                        B(B.mt-1-k, n       ), ldbk,
                                lalpha, B(B.mt-1-m, n       ), ldbm);
                        }
                    }
                }
            }
            /*
             *  PlasmaLeft / PlasmaUpper / Plasma[Conj]Trans
             */
            else {
                for (k = 0; k < B.mt; k++) {
                    tempkm = k == B.mt-1 ? B.m-k*B.mb : B.mb;
                    ldak = BLKLDD(A, k);
                    ldbk = BLKLDD(B, k);
                    lalpha = k == 0 ? alpha : zone;
                    for (n = 0; n < B.nt; n++) {
                        tempnn = n == B.nt-1 ? B.n-n*B.nb : B.nb;
                        QUARK_CORE_strsm(
                            plasma->quark, &task_flags,
                            side, uplo, trans, diag,
                            tempkm, tempnn, A.mb,
                            lalpha, A(k, k), ldak,
                                    B(k, n), ldbk);
                    }
                    for (m = k+1; m < B.mt; m++) {
                        tempmm = m == B.mt-1 ? B.m-m*B.mb : B.mb;
                        ldbm = BLKLDD(B, m);
                        for (n = 0; n < B.nt; n++) {
                            tempnn = n == B.nt-1 ? B.n-n*B.nb : B.nb;
                            QUARK_CORE_sgemm(
                                plasma->quark, &task_flags,
                                trans, PlasmaNoTrans,
                                tempmm, tempnn, B.mb, A.mb,
                                mzone,  A(k, m), ldak,
                                        B(k, n), ldbk,
                                lalpha, B(m, n), ldbm);
                        }
                    }
                }
            }
        }
        /*
         *  PlasmaLeft / PlasmaLower / PlasmaNoTrans
         */
        else {
            if (trans == PlasmaNoTrans) {
                for (k = 0; k < B.mt; k++) {
                    tempkm = k == B.mt-1 ? B.m-k*B.mb : B.mb;
                    ldak = BLKLDD(A, k);
                    ldbk = BLKLDD(B, k);
                    lalpha = k == 0 ? alpha : zone;
                    for (n = 0; n < B.nt; n++) {
                        tempnn = n == B.nt-1 ? B.n-n*B.nb : B.nb;
                        QUARK_CORE_strsm(
                            plasma->quark, &task_flags,
                            side, uplo, trans, diag,
                            tempkm, tempnn, A.mb,
                            lalpha, A(k, k), ldak,
                                    B(k, n), ldbk);
                    }
                    for (m = k+1; m < B.mt; m++) {
                        tempmm = m == B.mt-1 ? B.m-m*B.mb : B.mb;
                        ldam = BLKLDD(A, m);
                        ldbm = BLKLDD(B, m);
                        for (n = 0; n < B.nt; n++) {
                            tempnn = n == B.nt-1 ? B.n-n*B.nb : B.nb;
                            QUARK_CORE_sgemm(
                                plasma->quark, &task_flags,
                                PlasmaNoTrans, PlasmaNoTrans,
                                tempmm, tempnn, B.mb, A.mb,
                                mzone,  A(m, k), ldam,
                                        B(k, n), ldbk,
                                lalpha, B(m, n), ldbm);
                        }
                    }
                }
            }
            /*
             *  PlasmaLeft / PlasmaLower / Plasma[Conj]Trans
             */
            else {
                for (k = 0; k < B.mt; k++) {
                    tempkm = k == 0 ? B.m-(B.mt-1)*B.mb : B.mb;
                    ldak = BLKLDD(A, B.mt-1-k);
                    ldbk = BLKLDD(B, B.mt-1-k);
                    lalpha = k == 0 ? alpha : zone;
                    for (n = 0; n < B.nt; n++) {
                        tempnn = n == B.nt-1 ? B.n-n*B.nb : B.nb;
                        QUARK_CORE_strsm(
                            plasma->quark, &task_flags,
                            side, uplo, trans, diag,
                            tempkm, tempnn, A.mb,
                            lalpha, A(B.mt-1-k, B.mt-1-k), ldak,
                                    B(B.mt-1-k,        n), ldbk);
                    }
                    for (m = k+1; m < B.mt; m++) {
                        tempmm = m == B.mt-1 ? B.m-m*B.mb : B.mb;
                        ldbm = BLKLDD(B, B.mt-1-m);
                        for (n = 0; n < B.nt; n++) {
                            tempnn = n == B.nt-1 ? B.n-n*B.nb : B.nb;
                            QUARK_CORE_sgemm(
                                plasma->quark, &task_flags,
                                trans, PlasmaNoTrans,
                                B.mb, tempnn, tempkm, A.mb,
                                mzone,  A(B.mt-1-k, B.mt-1-m), ldak,
                                        B(B.mt-1-k, n       ), ldbk,
                                lalpha, B(B.mt-1-m, n       ), ldbm);
                        }
                    }
                }
            }
        }
    }
    /*
     *  PlasmaRight / PlasmaUpper / PlasmaNoTrans
     */
    else {
        if (uplo == PlasmaUpper) {
            if (trans == PlasmaNoTrans) {
                for (k = 0; k < B.nt; k++) {
                    tempkn = k == B.nt-1 ? B.n-k*B.nb : B.nb;
                    ldak = BLKLDD(A, k);
                    lalpha = k == 0 ? alpha : zone;
                    for (m = 0; m < B.mt; m++) {
                        tempmm = m == B.mt-1 ? B.m-m*B.mb : B.mb;
                        ldbm = BLKLDD(B, m);
                        QUARK_CORE_strsm(
                            plasma->quark, &task_flags,
                            side, uplo, trans, diag,
                            tempmm, tempkn, A.mb,
                            lalpha, A(k, k), ldak,  /* lda * tempkn */
                                    B(m, k), ldbm); /* ldb * tempkn */
                    }
                    for (m = 0; m < B.mt; m++) {
                        tempmm = m == B.mt-1 ? B.m-m*B.mb : B.mb;
                        ldbm = BLKLDD(B, m);
                        for (n = k+1; n < B.nt; n++) {
                            tempnn = n == B.nt-1 ? B.n-n*B.nb : B.nb;
                            QUARK_CORE_sgemm(
                                plasma->quark, &task_flags,
                                PlasmaNoTrans, PlasmaNoTrans,
                                tempmm, tempnn, B.mb, A.mb,
                                mzone,  B(m, k), ldbm,  /* ldb * B.mb   */
                                        A(k, n), ldak,  /* lda * tempnn */
                                lalpha, B(m, n), ldbm); /* ldb * tempnn */
                        }
                    }
                }
            }
            /*
             *  PlasmaRight / PlasmaUpper / Plasma[Conj]Trans
             */
            else {
                for (k = 0; k < B.nt; k++) {
                    tempkn = k == 0 ? B.n-(B.nt-1)*B.nb : B.nb;
                    ldak = BLKLDD(A, B.nt-1-k);
                    for (m = 0; m < B.mt; m++) {
                        tempmm = m == B.mt-1 ? B.m-m*B.mb : B.mb;
                        ldbm = BLKLDD(B, m);
                        QUARK_CORE_strsm(
                            plasma->quark, &task_flags,
                            side, uplo, trans, diag,
                            tempmm, tempkn, A.mb,
                            alpha, A(B.nt-1-k, B.nt-1-k), ldak,  /* lda * tempkn */
                                   B(       m, B.nt-1-k), ldbm); /* ldb * tempkn */

                        for (n = k+1; n < B.nt; n++) {
                            ldan = BLKLDD(A, B.nt-1-n);
                            QUARK_CORE_sgemm(
                                plasma->quark, &task_flags,
                                PlasmaNoTrans, trans,
                                tempmm, B.nb, tempkn, A.mb,
                                minvalpha, B(m,        B.nt-1-k), ldbm,  /* ldb  * tempkn */
                                           A(B.nt-1-n, B.nt-1-k), ldan, /* A.mb * tempkn (Never last row) */
                                zone,      B(m,        B.nt-1-n), ldbm); /* ldb  * B.nb   */
                        }
                    }
                }
            }
        }
        /*
         *  PlasmaRight / PlasmaLower / PlasmaNoTrans
         */
        else {
            if (trans == PlasmaNoTrans) {
                for (k = 0; k < B.nt; k++) {
                    tempkn = k == 0 ? B.n-(B.nt-1)*B.nb : B.nb;
                    ldak = BLKLDD(A, B.nt-1-k);
                    lalpha = k == 0 ? alpha : zone;
                    for (m = 0; m < B.mt; m++) {
                        tempmm = m == B.mt-1 ? B.m-m*B.mb : B.mb;
                        ldbm = BLKLDD(B, m);
                        QUARK_CORE_strsm(
                            plasma->quark, &task_flags,
                            side, uplo, trans, diag,
                            tempmm, tempkn, A.mb,
                            lalpha, A(B.nt-1-k, B.nt-1-k), ldak,  /* lda * tempkn */
                                    B(       m, B.nt-1-k), ldbm); /* ldb * tempkn */

                        for (n = k+1; n < B.nt; n++) {
                            QUARK_CORE_sgemm(
                                plasma->quark, &task_flags,
                                PlasmaNoTrans, PlasmaNoTrans,
                                tempmm, B.nb, tempkn, A.mb,
                                mzone,  B(m,        B.nt-1-k), ldbm,  /* ldb * tempkn */
                                        A(B.nt-1-k, B.nt-1-n), ldak,  /* lda * B.nb   */
                                lalpha, B(m,        B.nt-1-n), ldbm); /* ldb * B.nb   */
                        }
                    }
                }
            }
            /*
             *  PlasmaRight / PlasmaLower / Plasma[Conj]Trans
             */
            else {
                for (k = 0; k < B.nt; k++) {
                    tempkn = k == B.nt-1 ? B.n-k*B.nb : B.nb;
                    ldak = BLKLDD(A, k);
                    for (m = 0; m < B.mt; m++) {
                        tempmm = m == B.mt-1 ? B.m-m*B.mb : B.mb;
                        ldbm = BLKLDD(B, m);
                        QUARK_CORE_strsm(
                            plasma->quark, &task_flags,
                            side, uplo, trans, diag,
                            tempmm, tempkn, A.mb,
                            alpha, A(k, k), ldak,  /* lda * tempkn */
                                   B(m, k), ldbm); /* ldb * tempkn */

                        for (n = k+1; n < B.nt; n++) {
                            tempnn = n == B.nt-1 ? B.n-n*B.nb : B.nb;
                            ldan = BLKLDD(A, n);
                            QUARK_CORE_sgemm(
                                plasma->quark, &task_flags,
                                PlasmaNoTrans, trans,
                                tempmm, tempnn, B.mb, A.mb,
                                minvalpha, B(m, k), ldbm,  /* ldb  * tempkn */
                                           A(n, k), ldan, /* ldan * tempkn */
                                zone,      B(m, n), ldbm); /* ldb  * tempnn */
                        }
                    }
                }
            }
        }
    }
}

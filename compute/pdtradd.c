/**
 *
 * @file pdtradd.c
 *
 *  PLASMA auxiliary routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Emmanuel Agullo
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @generated d Fri Apr  1 11:02:58 2016
 *
 **/
#include "common.h"

#define A(m,n) BLKADDR(A, double, m, n)
#define B(m,n) BLKADDR(B, double, m, n)
/***************************************************************************//**
                                                                              *
                                                                              **/
void plasma_pdtradd(plasma_context_t *plasma)
{
    PLASMA_enum uplo, trans;
    double alpha, beta;
    PLASMA_desc A;
    PLASMA_desc B;
    PLASMA_sequence *sequence;
    PLASMA_request *request;

    int tempmm, tempnn;
    int m, n;
    int next_m;
    int next_n;
    int ldan, ldam, ldbm;

    plasma_unpack_args_8(uplo, trans, alpha, A, beta, B, sequence, request);
    if (sequence->status != PLASMA_SUCCESS)
        return;

    n = 0;
    m = PLASMA_RANK;
    switch(uplo){
    case PlasmaLower:
        while (m >= B.mt && n < B.nt) {
            n++;
            m = m-B.mt+n;
        }

        while (n < B.nt) {
            next_m = m;
            next_n = n;

            next_m += PLASMA_SIZE;
            while (next_m >= B.mt && next_n < B.nt) {
                next_n++;
                next_m = next_m-B.mt+next_n;
            }

            tempmm = m == B.mt-1 ? B.m-B.mb*m : B.nb;
            tempnn = n == B.nt-1 ? B.n-B.nb*n : B.nb;
            ldbm = BLKLDD(B, m);
            if (m == n) {
                ldam = BLKLDD(A, m);
                CORE_dtradd(uplo, trans, tempmm, tempnn,
                            alpha, A(m, n), ldam,
                            beta,  B(m, n), ldbm);
            }
            else {
                if (trans == PlasmaNoTrans) {
                    ldam = BLKLDD(A, m);
                    CORE_dgeadd(trans, tempmm, tempnn,
                                alpha, A(m, n), ldam,
                                beta,  B(m, n), ldbm);
                }
                else {
                    ldan = BLKLDD(A, n);
                    CORE_dgeadd(trans, tempmm, tempnn,
                                alpha, A(n, m), ldan,
                                beta,  B(m, n), ldbm);
                }
            }
            m = next_m;
            n = next_n;
        }
        break;

    case PlasmaUpper:
        while ((m >= min(n+1, B.mt)) && (n < B.nt)) {
            n++;
            m = m - min(n, B.mt);
        }

        while (n < B.nt) {
            next_m = m;
            next_n = n;

            next_m += PLASMA_SIZE;
            while ((next_m >= min(next_n+1, B.mt)) && (next_n < B.nt)) {
                next_n++;
                next_m = next_m - min(next_n, B.mt);
            }

            tempmm = m == B.mt-1 ? B.m-B.mb*m : B.nb;
            tempnn = n == B.nt-1 ? B.n-B.nb*n : B.nb;
            ldbm = BLKLDD(B, m);
            if (m == n) {
                ldam = BLKLDD(A, m);
                CORE_dtradd(uplo, trans, tempmm, tempnn,
                            alpha, A(m, n), ldam,
                            beta,  B(m, n), ldbm);
            }
            else {
                if (trans == PlasmaNoTrans) {
                    ldam = BLKLDD(A, m);
                    CORE_dgeadd(trans, tempmm, tempnn,
                                alpha, A(m, n), ldam,
                                beta,  B(m, n), ldbm);
                }
                else {
                    ldan = BLKLDD(A, n);
                    CORE_dgeadd(trans, tempmm, tempnn,
                                alpha, A(n, m), ldan,
                                beta,  B(m, n), ldbm);
                }
            }
            m = next_m;
            n = next_n;
        }
        break;

    case PlasmaUpperLower:
    default:
        while (m >= B.mt && n < B.nt) {
            n++;
            m = m-B.mt;
        }

        while (n < B.nt) {
            next_m = m;
            next_n = n;

            next_m += PLASMA_SIZE;
            while (next_m >= B.mt && next_n < B.nt) {
                next_n++;
                next_m = next_m-B.mt;
            }

            tempmm = m == B.mt-1 ? B.m-B.mb*m : B.nb;
            tempnn = n == B.nt-1 ? B.n-B.nb*n : B.nb;
            ldbm = BLKLDD(B, m);
            if (trans == PlasmaNoTrans) {
                ldam = BLKLDD(A, m);
                CORE_dgeadd(trans, tempmm, tempnn,
                            alpha, A(m, n), ldam,
                            beta,  B(m, n), ldbm);
            }
            else {
                ldan = BLKLDD(A, n);
                CORE_dgeadd(trans, tempmm, tempnn,
                            alpha, A(n, m), ldan,
                            beta,  B(m, n), ldbm);
            }
            m = next_m;
            n = next_n;
        }
    }
}

/***************************************************************************//**
                                                                              *
                                                                              **/
void plasma_pdtradd_quark(PLASMA_enum uplo, PLASMA_enum trans,
                          double alpha, PLASMA_desc A,
                          double beta,  PLASMA_desc B,
                          PLASMA_sequence *sequence, PLASMA_request *request)
{
    plasma_context_t *plasma;
    Quark_Task_Flags task_flags = Quark_Task_Flags_Initializer;

    int tempmm, tempnn, tempmn, tempnm;
    int m, n;
    int ldam, ldan, ldbm, ldbn;

    plasma = plasma_context_self();
    if (sequence->status != PLASMA_SUCCESS)
        return;
    QUARK_Task_Flag_Set(&task_flags, TASK_SEQUENCE, (intptr_t)sequence->quark_sequence);
    plasma_profile_by_function( &task_flags, GEADD );

    switch(uplo){
    case PlasmaLower:
        if (trans == PlasmaNoTrans) {
            for (n = 0; n < min(B.mt,B.nt); n++) {
                tempnm = n == B.mt-1 ? B.m-n*B.mb : B.mb;
                tempnn = n == B.nt-1 ? B.n-n*B.nb : B.nb;
                ldan = BLKLDD(A, n);
                ldbn = BLKLDD(B, n);

                QUARK_CORE_dtradd(
                    plasma->quark, &task_flags,
                    uplo, trans, tempnm, tempnn, B.mb,
                    alpha, A(n, n), ldan,
                    beta,  B(n, n), ldbn);

                for (m = n+1; m < B.mt; m++) {
                    tempmm = m == B.mt-1 ? B.m-B.mb*m : B.nb;
                    ldam = BLKLDD(A, m);
                    ldbm = BLKLDD(B, m);

                    QUARK_CORE_dgeadd(
                        plasma->quark, &task_flags,
                        trans, tempmm, tempnn, B.mb,
                        alpha, A(m, n), ldam,
                        beta,  B(m, n), ldbm);
                }
            }
        }
        else {
            for (n = 0; n < min(B.mt,B.nt); n++) {
                tempnm = n == B.mt-1 ? B.m-n*B.mb : B.mb;
                tempnn = n == B.nt-1 ? B.n-n*B.nb : B.nb;
                ldan = BLKLDD(A, n);
                ldbn = BLKLDD(B, n);

                QUARK_CORE_dtradd(
                    plasma->quark, &task_flags,
                    uplo, trans, tempnm, tempnn, B.mb,
                    alpha, A(n, n), ldan,
                    beta,  B(n, n), ldbn);

                for (m = n+1; m < B.mt; m++) {
                    tempmm = m == B.mt-1 ? B.m-B.mb*m : B.nb;
                    ldbm = BLKLDD(B, m);

                    QUARK_CORE_dgeadd(
                        plasma->quark, &task_flags,
                        trans, tempmm, tempnn, B.mb,
                        alpha, A(n, m), ldan,
                        beta,  B(m, n), ldbm);
                }
            }
        }
        break;
    case PlasmaUpper:
        if (trans == PlasmaNoTrans) {
            for (m = 0; m < min(B.mt,B.nt); m++) {
                tempmm = m == B.mt-1 ? B.m-B.mb*m : B.nb;
                tempmn = m == B.nt-1 ? B.n-m*B.nb : B.nb;
                ldam = BLKLDD(A, m);
                ldbm = BLKLDD(B, m);

                QUARK_CORE_dtradd(
                    plasma->quark, &task_flags,
                    uplo, trans, tempmm, tempmn, B.mb,
                    alpha, A(m, m), ldam,
                    beta,  B(m, m), ldbm);

                for (n = m+1; n < B.nt; n++) {
                    tempnn = n == B.nt-1 ? B.n-n*B.nb : B.nb;

                    QUARK_CORE_dgeadd(
                        plasma->quark, &task_flags,
                        trans, tempmm, tempnn, B.mb,
                        alpha, A(m, n), ldam,
                        beta,  B(m, n), ldbm);
                }
            }
        }
        else {
            for (m = 0; m < min(B.mt,B.nt); m++) {
                tempmm = m == B.mt-1 ? B.m-B.mb*m : B.nb;
                tempmn = m == B.nt-1 ? B.n-m*B.nb : B.nb;
                ldam = BLKLDD(A, m);
                ldbm = BLKLDD(B, m);

                QUARK_CORE_dtradd(
                    plasma->quark, &task_flags,
                    uplo, trans, tempmm, tempmn, B.mb,
                    alpha, A(m, m), ldam,
                    beta,  B(m, m), ldbm);

                for (n = m+1; n < B.nt; n++) {
                    tempnn = n == B.nt-1 ? B.n-n*B.nb : B.nb;
                    ldan = BLKLDD(A, n);

                    QUARK_CORE_dgeadd(
                        plasma->quark, &task_flags,
                        trans, tempmm, tempnn, B.mb,
                        alpha, A(n, m), ldan,
                        beta,  B(m, n), ldbm);
                }
            }
        }
        break;
    case PlasmaUpperLower:
    default:
        if (trans == PlasmaNoTrans) {
            for (m = 0; m < B.mt; m++) {
                tempmm = m == B.mt-1 ? B.m-B.mb*m : B.nb;
                ldam = BLKLDD(A, m);
                ldbm = BLKLDD(B, m);

                for (n = 0; n < B.nt; n++) {
                    tempnn = n == B.nt-1 ? B.n-n*B.nb : B.nb;

                    QUARK_CORE_dgeadd(
                        plasma->quark, &task_flags,
                        trans, tempmm, tempnn, B.mb,
                        alpha, A(m, n), ldam,
                        beta,  B(m, n), ldbm);
                }
            }
        }
        else {
            for (m = 0; m < B.mt; m++) {
                tempmm = m == B.mt-1 ? B.m-B.mb*m : B.nb;
                ldam = BLKLDD(A, m);
                ldbm = BLKLDD(B, m);

                for (n = 0; n < B.nt; n++) {
                    tempnn = n == B.nt-1 ? B.n-n*B.nb : B.nb;
                    ldan = BLKLDD(A, n);

                    QUARK_CORE_dgeadd(
                        plasma->quark, &task_flags,
                        trans, tempmm, tempnn, B.mb,
                        alpha, A(n, m), ldan,
                        beta,  B(m, n), ldbm);
                }
            }
        }
    }
}

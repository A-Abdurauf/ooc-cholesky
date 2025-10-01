/**
 *
 * @file pctrtri.c
 *
 *  PLASMA auxiliary routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Julien Langou
 * @author Henricus Bouwmeester
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @generated c Fri Apr  1 11:02:57 2016
 *
 **/
#include "common.h"

#define A(m,n) BLKADDR(A, PLASMA_Complex32_t, m, n)
/***************************************************************************//**
 *  Parallel tile triangular matrix inverse - dynamic scheduling
 **/
void plasma_pctrtri_quark(PLASMA_enum uplo, PLASMA_enum diag, PLASMA_desc A,
                          PLASMA_sequence *sequence, PLASMA_request *request)
{
    plasma_context_t *plasma;
    Quark_Task_Flags task_flags = Quark_Task_Flags_Initializer;

    int k, m, n;
    int ldam, ldak;
    int tempkn, tempkm, tempmm, tempnn;

    PLASMA_Complex32_t zone  = (PLASMA_Complex32_t) 1.0;
    PLASMA_Complex32_t mzone = (PLASMA_Complex32_t)-1.0;

    plasma = plasma_context_self();
    if (sequence->status != PLASMA_SUCCESS)
        return;
    QUARK_Task_Flag_Set(&task_flags, TASK_SEQUENCE, (intptr_t)sequence->quark_sequence);
    plasma_profile_by_function( &task_flags, TRTRI );

    /*
     *  PlasmaLower
     */
    if (uplo == PlasmaLower) {
        for (k = 0; k < A.nt; k++) {
            tempkn = k == A.nt-1 ? A.n-k*A.nb : A.nb;
            ldak = BLKLDD(A, k);
            for (m = k+1; m < A.mt; m++) {
                tempmm = m == A.mt-1 ? A.m-m*A.mb : A.mb;
                ldam = BLKLDD(A, m);
                QUARK_CORE_ctrsm(
                    plasma->quark, &task_flags,
                    PlasmaRight, uplo, PlasmaNoTrans, diag,
                    tempmm, tempkn, A.mb,
                    mzone, A(k, k), ldak,
                           A(m, k), ldam);
            }
            for (m = k+1; m < A.mt; m++) {
                tempmm = m == A.mt-1 ? A.m-m*A.mb : A.mb;
                ldam = BLKLDD(A, m);
                for (n = 0; n < k; n++) {
                    QUARK_CORE_cgemm(
                        plasma->quark, &task_flags,
                        PlasmaNoTrans, PlasmaNoTrans,
                        tempmm, A.nb, tempkn, A.mb,
                        zone, A(m, k), ldam,
                              A(k, n), ldak,
                        zone, A(m, n), ldam);
                }
            }
            for (n = 0; n < k; n++) {
                QUARK_CORE_ctrsm(
                    plasma->quark, &task_flags,
                    PlasmaLeft, uplo, PlasmaNoTrans, diag,
                    tempkn, A.nb, A.mb,
                    zone, A(k, k), ldak,
                          A(k, n), ldak);
            }
            QUARK_CORE_ctrtri(
                plasma->quark, &task_flags,
                uplo, diag,
                tempkn, A.mb,
                A(k, k), ldak,
                sequence, request, A.nb*k);
        }
    }
    /*
     *  PlasmaUpper
     */
    else {
        for (k = 0; k < A.mt; k++) {
            tempkm = k == A.mt-1 ? A.m-k*A.mb : A.mb;
            ldak = BLKLDD(A, k);
            for (n = k+1; n < A.nt; n++) {
                tempnn = n == A.nt-1 ? A.n-n*A.nb : A.nb;
                QUARK_CORE_ctrsm(
                    plasma->quark, &task_flags,
                    PlasmaLeft, uplo, PlasmaNoTrans, diag,
                    tempkm, tempnn, A.mb,
                    mzone, A(k, k), ldak,
                           A(k, n), ldak);
            }
            for (m = 0; m < k; m++) {
                ldam = BLKLDD(A, m);
                for (n = k+1; n < A.nt; n++) {
                    tempnn = n == A.nt-1 ? A.n-n*A.nb : A.nb;
                    QUARK_CORE_cgemm(
                        plasma->quark, &task_flags,
                        PlasmaNoTrans, PlasmaNoTrans,
                        A.mb, tempnn, tempkm, A.mb,
                        zone, A(m, k), ldam,
                              A(k, n), ldak,
                        zone, A(m, n), ldam);
                }
                QUARK_CORE_ctrsm(
                    plasma->quark, &task_flags,
                    PlasmaRight, uplo, PlasmaNoTrans, diag,
                    A.mb, tempkm, A.mb,
                    zone, A(k, k), ldak,
                          A(m, k), ldam);
            }
            QUARK_CORE_ctrtri(
                plasma->quark, &task_flags,
                uplo, diag,
                tempkm, A.mb,
                A(k, k), ldak,
                sequence, request, A.mb*k);
        }
    }
}

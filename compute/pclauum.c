/**
 *
 * @file pclauum.c
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
 *  Parallel UU' or L'L operation - dynamic scheduling
 **/
void plasma_pclauum_quark(PLASMA_enum uplo, PLASMA_desc A,
                          PLASMA_sequence *sequence, PLASMA_request *request)
{
    plasma_context_t *plasma;
    Quark_Task_Flags task_flags = Quark_Task_Flags_Initializer;

    int k, m, n;
    int ldak, ldam, ldan;
    int tempkm, tempkn;

    plasma = plasma_context_self();
    if (sequence->status != PLASMA_SUCCESS)
        return;
    QUARK_Task_Flag_Set(&task_flags, TASK_SEQUENCE, (intptr_t)sequence->quark_sequence);
    plasma_profile_by_function( &task_flags, LAUUM );

    /*
     *  PlasmaLower
     */
    if (uplo == PlasmaLower) {
        for (k = 0; k < A.mt; k++) {
            tempkm = k == A.mt-1 ? A.m-k*A.mb : A.mb;
            ldak = BLKLDD(A, k);
            for(n = 0; n < k; n++) {
                ldan = BLKLDD(A, n);
                QUARK_CORE_cherk(
                    plasma->quark, &task_flags,
                    uplo, PlasmaConjTrans,
                    A.mb, tempkm, A.mb,
                    1.0, A(k, n), ldak,
                    1.0, A(n, n), ldan);

                for(m = n+1; m < k; m++) {
                    ldam = BLKLDD(A, m);
                    QUARK_CORE_cgemm(
                        plasma->quark, &task_flags,
                        PlasmaConjTrans, PlasmaNoTrans,
                        A.mb, A.nb, tempkm, A.mb,
                        1.0, A(k, m), ldak,
                             A(k, n), ldak,
                        1.0, A(m, n), ldam);
                }
            }
            for (n = 0; n < k; n++) {
                QUARK_CORE_ctrmm(
                    plasma->quark, &task_flags,
                    PlasmaLeft, uplo, PlasmaConjTrans, PlasmaNonUnit,
                    tempkm, A.nb, A.mb,
                    1.0, A(k, k), ldak,
                         A(k, n), ldak);
            }
            QUARK_CORE_clauum(
                plasma->quark, &task_flags,
                uplo, tempkm, A.mb,
                A(k, k), ldak);
        }
    }
    /*
     *  PlasmaUpper
     */
    else {
        for (k = 0; k < A.mt; k++) {
            tempkn = k == A.nt-1 ? A.n-k*A.nb : A.nb;
            ldak = BLKLDD(A, k);

            for (m = 0; m < k; m++) {
                ldam = BLKLDD(A, m);
                QUARK_CORE_cherk(
                    plasma->quark, &task_flags,
                    uplo, PlasmaNoTrans,
                    A.mb, tempkn, A.mb,
                    1.0, A(m, k), ldam,
                    1.0, A(m, m), ldam);

                for (n = m+1; n < k; n++){
                    ldan = BLKLDD(A, n);
                    QUARK_CORE_cgemm(
                        plasma->quark, &task_flags,
                        PlasmaNoTrans, PlasmaConjTrans,
                        A.mb, A.nb, tempkn, A.mb,
                        1.0, A(m, k), ldam,
                             A(n, k), ldan,
                        1.0, A(m, n), ldam);
                }
            }
            for (m = 0; m < k; m++) {
                ldam = BLKLDD(A, m);
                QUARK_CORE_ctrmm(
                    plasma->quark, &task_flags,
                    PlasmaRight, uplo, PlasmaConjTrans, PlasmaNonUnit,
                    A.mb, tempkn, A.mb,
                    1.0, A(k, k), ldak,
                         A(m, k), ldam);
            }
            QUARK_CORE_clauum(
                plasma->quark, &task_flags,
                uplo, tempkn, A.mb,
                A(k, k), ldak);
        }
    }
}

/**
 *
 * @file pslascal.c
 *
 *  PLASMA auxiliary routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @generated s Fri Apr  1 11:02:58 2016
 *
 **/
#include "common.h"

#define A(m,n) BLKADDR(A, float, m, n)
/***************************************************************************//**
 *  Parallel scale of a  matrix A.
 **/
void plasma_pslascal_quark(PLASMA_enum uplo, float alpha, PLASMA_desc A,
                           PLASMA_sequence *sequence, PLASMA_request *request)
{
    plasma_context_t *plasma;
    Quark_Task_Flags task_flags = Quark_Task_Flags_Initializer;

    int tempmm, tempnn, tempmn, tempnm;
    int m, n;
    int ldam, ldan;
    int minmnt = min(A.mt, A.nt);

    plasma = plasma_context_self();
    if (sequence->status != PLASMA_SUCCESS)
        return;

    QUARK_Task_Flag_Set(&task_flags, TASK_SEQUENCE, (intptr_t)sequence->quark_sequence);
    plasma_profile_by_function( &task_flags, LASCL );

    switch(uplo) {
    case PlasmaLower:
        for (n = 0; n < minmnt; n++) {
            tempnm = n == A.mt-1 ? A.m-n*A.mb : A.mb;
            tempnn = n == A.nt-1 ? A.n-n*A.nb : A.nb;
            ldan = BLKLDD(A, n);

            QUARK_CORE_slascal(
                plasma->quark, &task_flags,
                PlasmaUpper, tempnm, tempnn, A.mb,
                alpha, A(n, n), ldan);

            for (m = n+1; m < A.mt; m++) {
                tempmm = m == A.mt-1 ? A.m-A.mb*m : A.nb;
                ldam = BLKLDD(A, m);

                QUARK_CORE_slascal(
                    plasma->quark, &task_flags,
                    PlasmaUpperLower, tempmm, tempnn, A.mb,
                    alpha, A(m, n), ldam);
            }
        }
        break;

    case PlasmaUpper:
        for (m = 0; m < minmnt; m++) {
            tempmm = m == A.mt-1 ? A.m-A.mb*m : A.nb;
            tempmn = m == A.nt-1 ? A.n-m*A.nb : A.nb;
            ldam = BLKLDD(A, m);

            QUARK_CORE_slascal(
                plasma->quark, &task_flags,
                PlasmaUpper, tempmm, tempmn, A.mb,
                alpha, A(m, m), ldam);

            for (n = m+1; n < A.nt; n++) {
                tempnn = n == A.nt-1 ? A.n-n*A.nb : A.nb;

                QUARK_CORE_slascal(
                    plasma->quark, &task_flags,
                    PlasmaUpperLower, tempmm, tempnn, A.mb,
                    alpha, A(m, n), ldam);
            }
        }
        break;

    case PlasmaUpperLower:
    default:
        for (m = 0; m < A.mt; m++) {
            tempmm = m == A.mt-1 ? A.m-A.mb*m : A.nb;
            ldam = BLKLDD(A, m);

            for (n = 0; n < A.nt; n++) {
                tempnn = n == A.nt-1 ? A.n-n*A.nb : A.nb;

                QUARK_CORE_slascal(
                    plasma->quark, &task_flags,
                    PlasmaUpperLower, tempmm, tempnn, A.mb,
                    alpha, A(m, n), ldam);
            }
        }
    }
}

/**
 *
 * @file pdlaswp.c
 *
 *  PLASMA auxiliary routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @generated d Fri Apr  1 11:02:59 2016
 *
 **/
#include "common.h"

#define B(m, n) BLKADDR(B, double, m, n)
#define IPIV(k) &(IPIV[(int64_t)B.mb*(int64_t)(k)])

/***************************************************************************//**
 *  Parallel tile row interchanges - dynamic scheduling
 **/
void plasma_pdlaswp_quark(PLASMA_desc B, const int *IPIV, int inc,
                          PLASMA_sequence *sequence, PLASMA_request *request)
{
    plasma_context_t *plasma;
    Quark_Task_Flags task_flags = Quark_Task_Flags_Initializer;

    int m, n;
    int tempi, tempm, tempmm, tempnn;

    plasma = plasma_context_self();
    if (sequence->status != PLASMA_SUCCESS)
        return;
    QUARK_Task_Flag_Set(&task_flags, TASK_SEQUENCE, (intptr_t)sequence->quark_sequence);
    plasma_profile_by_function( &task_flags, LASWP );

    if ( inc > 0 )
    {
        for (m = 0; m < B.mt; m++) {
            tempi = m * B.mb;
            tempm = B.m - tempi;
            tempmm = m == B.mt-1 ? tempm : B.mb;

            for (n = 0; n < B.nt; n++) {
                tempnn = n == B.nt-1 ? B.n - n * B.nb : B.nb;

                QUARK_CORE_dlaswp_ontile(
                    plasma->quark, &task_flags,
                    plasma_desc_submatrix(B, tempi, n*B.nb, tempm, tempnn),
                    B(m, n), 1, tempmm, IPIV(m), inc, B(B.mt-1, n) );
            }
        }
    }
    else
    {
        for (m = B.mt-1; m > -1; m--) {
            tempi = m * B.mb;
            tempm = B.m - tempi;
            tempmm = m == B.mt-1 ? tempm : B.mb;

            for (n = 0; n < B.nt; n++) {
                tempnn = n == B.nt-1 ? B.n - n * B.nb : B.nb;

                QUARK_CORE_dlaswp_ontile(
                    plasma->quark, &task_flags,
                    plasma_desc_submatrix(B, tempi, n*B.nb, tempm, tempnn),
                    B(m, n), 1, tempmm, IPIV(m), inc, B(0, n) );
            }
        }
    }
}

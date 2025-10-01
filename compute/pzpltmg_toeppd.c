/**
 *
 * @file pzpltmg_toeppd.c
 *
 *  PLASMA auxiliary routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @precisions normal z -> s d c
 *
 **/
#include <stdlib.h>
#include "common.h"

#define A(m, n) BLKADDR(A, PLASMA_Complex64_t, m, n)

/***************************************************************************//**
 *  Parallel tile Toeppd matrix generation - dynamic scheduling
 **/
void plasma_pzpltmg_toeppd_quark( PLASMA_desc A, unsigned long long int seed,
                                  PLASMA_sequence *sequence, PLASMA_request *request )
{
    plasma_context_t *plasma;
    Quark_Task_Flags task_flags = Quark_Task_Flags_Initializer;

    PLASMA_Complex64_t **work;
    int m, n, k;
    int ldam;
    int tempm0, tempn0, tempk0;
    int tempmm, tempnn, tempkm;

    plasma = plasma_context_self();
    if (sequence->status != PLASMA_SUCCESS)
        return;
    QUARK_Task_Flag_Set(&task_flags, TASK_SEQUENCE, (intptr_t)sequence->quark_sequence);
    plasma_profile_by_function( &task_flags, PLRNT );

    work = (PLASMA_Complex64_t**) malloc ( A.mt * sizeof( PLASMA_Complex64_t* ) );

    /* Generate a random vector of size A.m */
    for (k = 0; k < A.mt; k++) {
        tempk0 = k * A.mb;
        tempkm = k == A.mt-1 ? A.m-k*A.mb : A.mb;

        /* Allocate temporary vector and initialize it randomly */
        work[k] = (PLASMA_Complex64_t*)plasma_shared_alloc(plasma, 2*tempkm, PlasmaComplexDouble);

        QUARK_CORE_zpltmg_toeppd1(
            plasma->quark, &task_flags,
            A.m, tempk0+1, tempkm, work[k], seed );
    }

    /* Compute the toeppd Matrix tile by tile */
    for (m = 0; m < A.mt; m++) {
        tempm0 = m * A.mb;
        tempmm = m == A.mt-1 ? A.m-m*A.mb : A.mb;
        ldam = BLKLDD(A, m);

        for (n = 0; n < A.nt; n++) {
            tempnn = n == A.nt-1 ? A.n-n*A.nb : A.nb;
            tempn0 = n * A.nb;

            QUARK_CORE_zlaset(
                plasma->quark, &task_flags,
                PlasmaUpperLower, tempmm, tempnn, 0., 0.,
                A(m, n), ldam);

            for (k = 0; k < A.mt; k++) {
                tempkm = k == A.mt-1 ? A.m-k*A.mb : A.mb;

                QUARK_CORE_zpltmg_toeppd2(
                    plasma->quark, &task_flags,
                    tempmm, tempnn, tempkm,
                    tempm0, tempn0, work[k],
                    A(m, n), ldam );
            }
        }
    }

    /* Submit the workspace free */
    for (k = 0; k < A.mt; k++) {
        tempk0 = k * A.mb;
        tempkm = k == A.mt-1 ? A.m-k*A.mb : A.mb;

        QUARK_CORE_free(plasma->quark, &task_flags,
                        work[k], (2*tempkm)*sizeof(PLASMA_Complex64_t));
    }

    /* We can free work and loose all pointers because they are already saved by Quark */
    free(work);
}

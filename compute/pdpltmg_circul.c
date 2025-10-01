/**
 *
 * @file pdpltmg_circul.c
 *
 *  PLASMA auxiliary routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Ichitaro Yamazaki
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @generated d Fri Apr  1 11:02:59 2016
 *
 **/
#include <stdlib.h>
#include "common.h"

#define A(m,n) BLKADDR(A, double, m, n)

/***************************************************************************//**
 *  Parallel tile circulant matrix generation -- Dynamic scheduling
 *
 *  See http://www.mathworks.fr/fr/help/matlab/ref/gallery.html#f84-999880
 *
 *  Circulant matrix
 *
 *  A circulant matrix has the property that each row is obtained from the
 *  previous one by cyclically permuting the entries one step forward. It is a
 *  special Toeplitz matrix in which the diagonals "wrap around."
 *
 *  The eigensystem of C (n-by-n) is known explicitly: If t is an nth root of
 *  unity, then the inner product of v and w = [1 t t2 ... t(n – 1)] is an
 *  eigenvalue of C and w(n:-1:1) is an eigenvector, where v is the first column of
 *  C.
 *
 */
void plasma_pdpltmg_circul_quark( PLASMA_desc A, unsigned long long int seed,
                                 PLASMA_sequence *sequence, PLASMA_request *request )
{
    plasma_context_t *plasma;
    Quark_Task_Flags task_flags = Quark_Task_Flags_Initializer;
    double *V;
    int m, n, ldam;
    int tempm0, tempn0;
    int tempmm, tempnn;

    plasma = plasma_context_self();
    if (sequence->status != PLASMA_SUCCESS)
        return;
    QUARK_Task_Flag_Set(&task_flags, TASK_SEQUENCE, (intptr_t)sequence->quark_sequence);
    plasma_profile_by_function( &task_flags, PLRNT );

    /*
     * Allocate and initialize the first row of the circulant matrix
     */
    V = (double*) plasma_shared_alloc( plasma, A.m, PlasmaRealDouble );

    QUARK_CORE_dplrnt( plasma->quark, &task_flags,
                       A.m, 1, V, A.m, A.m, 0, 0, seed );

    /* Initialize the full matrix */
    for (m = 0; m < A.mt; m++) {
        tempm0 = m * A.mb;
        tempmm = m == A.mt-1 ? A.m - tempm0 : A.mb;
        ldam = BLKLDD(A, m);

        for (n = 0; n < A.nt; n++) {
            tempn0 = n * A.nb;
            tempnn = n == A.nt-1 ? A.n - tempn0 : A.nb;

            QUARK_CORE_dpltmg_circul(
                plasma->quark, &task_flags,
                tempmm, tempnn, A(m, n), ldam,
                A.m, tempm0, tempn0, V );
        }
    }

    QUARK_CORE_free(plasma->quark, &task_flags,
                    V, (A.m)*sizeof(double));
}

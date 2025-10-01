/**
 *
 * @file pzpltmg_hankel.c
 *
 *  PLASMA auxiliary routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Ichitaro Yamazaki
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @precisions normal z -> s d c
 *
 **/
#include <stdlib.h>
#include "common.h"

#define A(m,n) BLKADDR(A, PLASMA_Complex64_t, m, n)

/***************************************************************************//**
 *  Parallel tile Hankel matrix generation -- Dynamic scheduling
 *
 *  See http://en.wikipedia.org/wiki/Hankel_matrix
 *
 *  Hankel matrix
 *
 *  In linear algebra, a Hankel matrix (or catalecticant matrix), named after
 *  Hermann Hankel, is a square matrix with constant skew-diagonals (positive
 *  sloping diagonals), e.g.:
 *
 *  \begin{bmatrix}
 *  a & b & c & d & e \\
 *  b & c & d & e & f \\
 *  c & d & e & f & g \\
 *  d & e & f & g & h \\
 *  e & f & g & h & i \\
 *  \end{bmatrix}.
 *
 *  A(i,j) = A(i-1,j+1)
 *
 */
void plasma_pzpltmg_hankel_quark( PLASMA_desc A, unsigned long long int seed,
                                 PLASMA_sequence *sequence, PLASMA_request *request )
{
    plasma_context_t *plasma;
    Quark_Task_Flags task_flags = Quark_Task_Flags_Initializer;
    PLASMA_Complex64_t **V;
    int sizeMT, sizeM;
    int m, n, k;
    int ldam;
    int tempm0, tempn0, tempk0;
    int tempmm, tempnn;

    plasma = plasma_context_self();
    if (sequence->status != PLASMA_SUCCESS)
        return;
    QUARK_Task_Flag_Set(&task_flags, TASK_SEQUENCE, (intptr_t)sequence->quark_sequence);
    plasma_profile_by_function( &task_flags, PLRNT );

    /* Generate a random vector of size A.m+A.n-1 */
    sizeMT = A.mt      + A.nt;
    sizeM  = A.mt*A.mb + A.n - 1;

    V = (PLASMA_Complex64_t**) malloc ( sizeMT * sizeof( PLASMA_Complex64_t* ) );
    for (k = 0; k < sizeMT; k++) {
        tempk0 = k * A.mb;

        /* Allocate temporary vector and initialize it randomly */
        V[k] = (PLASMA_Complex64_t*)plasma_shared_alloc(plasma, A.mb, PlasmaComplexDouble);

        QUARK_CORE_zplrnt(
            plasma->quark, &task_flags,
            A.mb, 1, V[k], A.mb,
            sizeM, tempk0+1, 0, seed );
    }

    /* Initialize the full matrix */
    for (m = 0; m < A.mt; m++) {
        tempm0 = m * A.mb;
        tempmm = m == A.mt-1 ? A.m - tempm0 : A.mb;
        ldam = BLKLDD(A, m);

        for (n = 0; n < A.nt; n++) {
            tempn0 = n * A.nb;
            tempnn = n == A.nt-1 ? A.n - tempn0 : A.nb;

            QUARK_CORE_zpltmg_hankel(
                plasma->quark, &task_flags,
                PlasmaUpperLower, tempmm, tempnn, A(m, n), ldam,
                tempm0, tempn0, A.mb, V[m+n], V[m+n+1]);
        }
    }

    /* Submit the workspace free */
    for (k = 0; k < sizeMT; k++) {
        QUARK_CORE_free(plasma->quark, &task_flags,
                        V[k], A.mb*sizeof(PLASMA_Complex64_t));
    }

    /* We can free work and loose all pointers because they are already saved by Quark */
    free(V);
}

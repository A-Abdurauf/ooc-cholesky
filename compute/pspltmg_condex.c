/**
 *
 * @file pspltmg_condex.c
 *
 *  PLASMA auxiliary routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @generated s Fri Apr  1 11:02:59 2016
 *
 **/
#include "common.h"

#define A(m, n) BLKADDR(A, float, m, n)

/***************************************************************************//**
 *  Parallel tile Condex matrix generation
 *  sequential call to either static or dynamic scheduling
 *
 *  See http://www.mathworks.fr/fr/help/matlab/ref/gallery.html#f84-999898
 *  gallery('condex',n,4,100)
 *
 *  Returns a "counter-example" matrix to a condition estimator. It has order n
 *  and scalar parameter theta (default 100).
 *
 * LAPACK (RCOND): It is the inverse of this matrix that is a counter-example.
 *
 **/
void plasma_spltmg_condex( PLASMA_desc A,
                           PLASMA_sequence *sequence, PLASMA_request *request )
{
    plasma_context_t *plasma;
    Quark_Task_Flags task_flags = Quark_Task_Flags_Initializer;
    PLASMA_desc descQ;
    float *Q;
    float theta = 100.;

    plasma = plasma_context_self();
    if (sequence->status != PLASMA_SUCCESS)
        return;
    QUARK_Task_Flag_Set(&task_flags, TASK_SEQUENCE, (intptr_t)sequence->quark_sequence);
    plasma_profile_by_function( &task_flags, PLRNT );

    /*
     * Allocate the descriptor
     * Exploit the fact that allocated buffer has the same
     * size in PlasmaCM, and PlasmaCCRB formats
     */
    Q = (float*) plasma_shared_alloc( plasma, A.m * 3, PlasmaRealFloat );
    descQ = plasma_desc_init(
        PlasmaRealFloat, A.mb, A.nb, A.mb*A.nb,
        A.m, 3, 0, 0, A.m, 3 );
    descQ.mat = Q;

    /*
     * Direct call to the kernel, this matrix generation and
     * is not on the critical path, or at least it should not.
     */
    CORE_spltmg_condexq( descQ.m, A.n, descQ.mat, descQ.m );

    /* Submit kernels for layout conversion (Mostly to keep track of dependencies) */
    PLASMA_sgecfi_Async( descQ.m, 3, descQ.mat,
                         PlasmaCM,   descQ.m,  3,
                         PlasmaCCRB, descQ.mb, descQ.nb,
                         sequence, request );

    /* Init A to (1. + theta) * I */
    plasma_dynamic_call_6(plasma_pslaset,
                          PLASMA_enum,        PlasmaUpperLower,
                          float, 0.,
                          float, 1. + theta,
                          PLASMA_desc,        A,
                          PLASMA_sequence*,   sequence,
                          PLASMA_request*,    request);

    /* Computes A = A - theta Q * Q' */
    plasma_parallel_call_9(plasma_psgemm,
                           PLASMA_enum,        PlasmaNoTrans,
                           PLASMA_enum,        PlasmaTrans,
                           float, (-theta),
                           PLASMA_desc,        descQ,
                           PLASMA_desc,        descQ,
                           float, 1.,
                           PLASMA_desc,        A,
                           PLASMA_sequence*,   sequence,
                           PLASMA_request*,    request);

    plasma_dynamic_sync();
    plasma_shared_free( plasma, Q );
}

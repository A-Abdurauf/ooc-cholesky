/**
 *
 * @file pzpltmg_house.c
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
#ifdef PLASMA_WITH_MKL
#include <mkl_lapacke.h>
#else
#include <lapacke.h>
#endif
#include "common.h"

/***************************************************************************//**
 *  Parallel tile Householder matrix generation
 *  sequential call to either static or dynamic scheduling
 *
 * See http://www.mathworks.fr/fr/help/matlab/ref/gallery.html#f84-999993
 *
 * Householder matrix
 *
 * Generates a random column vector of size M, and returns the housholder matrix
 * H = eye(n,n) - beta*v*v' that satisfies the relationship
 *
 * H*x = -sign(x(1))*norm(x)*e1
 *
 * where e1 is the first column of eye(n,n). Note that if x is complex, then
 * sign(x) exp(i*arg(x)) (which equals x./abs(x) when x is nonzero).
 *
 */
void plasma_zpltmg_house( PLASMA_desc A, unsigned long long int seed,
                          PLASMA_sequence *sequence, PLASMA_request *request )
{
    plasma_context_t *plasma;
    Quark_Task_Flags task_flags = Quark_Task_Flags_Initializer;
    PLASMA_desc descV;
    PLASMA_Complex64_t *V;
    PLASMA_Complex64_t tau;

    plasma = plasma_context_self();
    if (sequence->status != PLASMA_SUCCESS)
        return;
    QUARK_Task_Flag_Set(&task_flags, TASK_SEQUENCE, (intptr_t)sequence->quark_sequence);
    plasma_profile_by_function( &task_flags, PLRNT );

    /*
     * Submit A initialization, such that other thread have work while main
     * thread is initializing V
     */
    plasma_dynamic_call_6(plasma_pzlaset,
                          PLASMA_enum,        PlasmaUpperLower,
                          PLASMA_Complex64_t, 0.,
                          PLASMA_Complex64_t, 1.,
                          PLASMA_desc,        A,
                          PLASMA_sequence*,   sequence,
                          PLASMA_request*,    request);

    /*
     * Allocate the descriptor of the Householder reflector.
     * Exploit the fact that for one column, the PlasmaCM and PlasmaCCRB storage
     * are identical
     */
    V = (PLASMA_Complex64_t*) plasma_shared_alloc( plasma, A.m, PlasmaComplexDouble );
    descV = plasma_desc_init(
        PlasmaComplexDouble, A.mb, A.nb, A.mb*A.nb,
        A.m, 1, 0, 0, A.m, 1 );
    descV.mat = V;

    /* Initialize Householder vector */
    {
        CORE_zplrnt( A.m, 1, V, A.m, A.m, 0, 0, seed );
        LAPACKE_zlarfg_work( A.m, V, V+1, 1, &tau );
        V[0] = 1.;
    }

    /* Computes A = A - tau * V * V' (gerc == gemm) */
    plasma_parallel_call_9(plasma_pzgemm,
                           PLASMA_enum,        PlasmaNoTrans,
                           PLASMA_enum,        PlasmaConjTrans,
                           PLASMA_Complex64_t, (-tau),
                           PLASMA_desc,        descV,
                           PLASMA_desc,        descV,
                           PLASMA_Complex64_t, 1.,
                           PLASMA_desc,        A,
                           PLASMA_sequence*,   sequence,
                           PLASMA_request*,    request);

    plasma_dynamic_sync();
    plasma_shared_free( plasma, V );
}

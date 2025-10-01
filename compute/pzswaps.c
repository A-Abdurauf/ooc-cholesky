/**
 * @file pzswaps.c
 *
 *  PLASMA computational routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Gregoire Pichon
 * @author Azzam Haidar
 * @date 2014-07
 * @precisions normal z -> s d c
 *
 **/
#ifdef PLASMA_WITH_MKL
#include <mkl_lapacke.h>
#else
#include <lapacke.h>
#endif
#include <math.h>
#include "common.h"

/***************************************************************************//**
 *
 * @ingroup CORE_PLASMA_Complex64_t
 *
 *  plasma_pzswaps_quark - swap vectors in Z according to the permutation
 *  given starts by copying Z in work and then copying-back according
 *  the permutation given
 *
 *******************************************************************************
 *
 * @param[in] n
 *          n specifies the dimension of the matrix
 *
 * @param[in] IWORK
 *          The permutation used to sort Z
 *
 * @param[in,out] Z
 *          On entry, the non-sorted eigenvectors
 *          On exit, the sorted eigenvectors
 *
 * @param[in] LDZ
 *          LDZ specifies the leading dimension of Z
 *
 * @param[out] work
 *          Space used to copy Z
 *
 * @param[in] sequence
 *          Identifies the sequence of function calls that this call belongs to
 *          (for completion checks and exception handling purposes).
 *
 * @param[out] request
 *          Identifies this function call (for exception handling purposes).
 *
***************************************************************************/
void plasma_pzswaps_quark(int n, int *IWORK,
                          PLASMA_Complex64_t *Z, int LDZ,
                          PLASMA_Complex64_t *work,
                          PLASMA_sequence *sequence, PLASMA_request *request)
{
    plasma_context_t *plasma;
    plasma = plasma_context_self();

    Quark_Task_Flags task_flags = Quark_Task_Flags_Initializer;
    if (sequence->status != PLASMA_SUCCESS)
         return;
    QUARK_Task_Flag_Set(&task_flags, TASK_SEQUENCE, (intptr_t)sequence->quark_sequence);
    plasma_profile_by_function( &task_flags, LASWP );

    int i;
    int nb = plasma->ev_tasknb;
    int task_size;

    for(i=0; i<n; i+=nb){
        task_size = min(nb, n-i);
        QUARK_CORE_zswap(plasma->quark, &task_flags,
                         n, n, Z, LDZ, work, IWORK,
                         i, i+task_size);
    }
}

/**
 *
 * @file pzlaset_identity.c
 *
 *  PLASMA auxiliary routines
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
#include "common.h"
/***************************************************************************//**
 *  Parallel initialization a 1-D to identity
 **/
void plasma_pzlaset_identity_quark(int n,
                                   PLASMA_Complex64_t *A, int lda,
                                   PLASMA_sequence *sequence, PLASMA_request *request)
{
    plasma_context_t *plasma;
    Quark_Task_Flags task_flags = Quark_Task_Flags_Initializer;

    int i;

    plasma = plasma_context_self();
    if (sequence->status != PLASMA_SUCCESS)
        return;

    QUARK_Task_Flag_Set(&task_flags, TASK_SEQUENCE, (intptr_t)sequence->quark_sequence);
    plasma_profile_by_function( &task_flags, LASET );

    int nb = plasma->ev_tasknb;
    int task_size;

    for (i=0; i<n; i+=nb){
        task_size = min(nb, n-i);
        QUARK_CORE_zlaset_identity(plasma->quark, &task_flags,
                                   n, i, task_size,
                                   A);
    }
}

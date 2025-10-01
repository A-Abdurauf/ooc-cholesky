/**
 *
 * @file pclaset_identity.c
 *
 *  PLASMA auxiliary routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Gregoire Pichon
 * @author Azzam Haidar
 * @date 2014-07
 * @generated c Fri Apr  1 11:03:00 2016
 *
 **/
#include "common.h"
/***************************************************************************//**
 *  Parallel initialization a 1-D to identity
 **/
void plasma_pclaset_identity_quark(int n,
                                   PLASMA_Complex32_t *A, int lda,
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
        QUARK_CORE_claset_identity(plasma->quark, &task_flags,
                                   n, i, task_size,
                                   A);
    }
}

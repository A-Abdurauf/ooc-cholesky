/**
 *
 * @file qwrapper_sDC_fakedep.c
 *
 *  PLASMA core_blas quark wrapper
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Gregoire Pichon
 * @author Azzam Haidar
 * @date 2014-07
 * @generated s Fri Apr  1 11:02:43 2016
 *
 **/
#include "common.h"

/***************************************************************************//**
 *
 **/
void
CORE_sDC_fakedep_quark(Quark *quark)
{
    int nb;
    quark_unpack_args_1(quark, nb);
    return;
}

/***************************************************************************//**
 *
 **/
void
QUARK_CORE_sDC_fakedep(Quark *quark, Quark_Task_Flags *task_flags,
                       int nb_tasks, int nb, float *Q, int LDQ, float *W)
{
    //plasma_profile_by_kernel( task_flags, FAKE_DC );

    /* DEP: W_red+LDQ*i for i=0..nb_tasks */
    Quark_Task *task;
    int i;

    task = QUARK_Task_Init( quark, CORE_sDC_fakedep_quark,   task_flags);
    QUARK_Task_Pack_Arg(quark, task, sizeof(int),    &nb, VALUE );
    QUARK_Task_Pack_Arg(quark, task, sizeof(float),  W,      OUTPUT );
    for (i=0; i<nb_tasks; i++) {
        QUARK_Task_Pack_Arg(quark, task, sizeof(float), Q+LDQ*i*nb, INPUT );
    }
    QUARK_Insert_Task_Packed(quark, task);
}

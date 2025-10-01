/**
 *
 * @file qwrapper_slaed3_freebigwork.c
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
#include <stdlib.h>
#include "common.h"

void
CORE_slaed3_freebigwork_quark(Quark *quark);

/***************************************************************************//**
 *
 **/
void QUARK_CORE_slaed3_freebigwork(Quark *quark, Quark_Task_Flags *task_flags,
                                   int *K_bis, int wsmode, float **WORK)
{
    plasma_gendag_by_kernel( task_flags, LAED3_FREE );

    QUARK_Insert_Task(quark, CORE_slaed3_freebigwork_quark, task_flags,
        sizeof(int),      &wsmode,  VALUE,
        sizeof(float**),  WORK,      INOUT,
        sizeof(int),       K_bis, ((wsmode == 3 || wsmode == 5) ? INOUT : INPUT),
        0);
}
/***************************************************************************//**
 *
 **/
void
CORE_slaed3_freebigwork_quark(Quark *quark)
{
    int wsmode;
    float **WORK;
    int *K;

    quark_unpack_args_3(quark, wsmode, WORK, K);

    if((wsmode == 1) || (wsmode == 3)){
        free(*WORK);
        *WORK = NULL;
    }
    if((wsmode == 3) || (wsmode == 5)){
        free(WORK);
        WORK=NULL;
    }
}

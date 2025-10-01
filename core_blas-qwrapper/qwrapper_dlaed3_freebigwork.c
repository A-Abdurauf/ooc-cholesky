/**
 *
 * @file qwrapper_dlaed3_freebigwork.c
 *
 *  PLASMA core_blas quark wrapper
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Gregoire Pichon
 * @author Azzam Haidar
 * @date 2014-07
 * @precisions normal d -> s
 *
 **/
#include <stdlib.h>
#include "common.h"

void
CORE_dlaed3_freebigwork_quark(Quark *quark);

/***************************************************************************//**
 *
 **/
void QUARK_CORE_dlaed3_freebigwork(Quark *quark, Quark_Task_Flags *task_flags,
                                   int *K_bis, int wsmode, double **WORK)
{
    plasma_gendag_by_kernel( task_flags, LAED3_FREE );

    QUARK_Insert_Task(quark, CORE_dlaed3_freebigwork_quark, task_flags,
        sizeof(int),      &wsmode,  VALUE,
        sizeof(double**),  WORK,      INOUT,
        sizeof(int),       K_bis, ((wsmode == 3 || wsmode == 5) ? INOUT : INPUT),
        0);
}
/***************************************************************************//**
 *
 **/
void
CORE_dlaed3_freebigwork_quark(Quark *quark)
{
    int wsmode;
    double **WORK;
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

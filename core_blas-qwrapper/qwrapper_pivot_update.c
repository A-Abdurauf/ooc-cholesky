/**
 *
 * @file qwrapper_pivot_update.c
 *
 *  PLASMA core_blas quark wrapper
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Mathieu Faverge
 * @author Ichitaro Yamazaki
 * @date 2010-11-15
 *
 **/
#include "common.h"

/***************************************************************************//**
 *
 **/
void CORE_pivot_update_quark(Quark *quark)
{
    int m, n, offset, init;
    int *indices;
    int *ipiv;

    quark_unpack_args_6(quark, m, n, ipiv, indices, offset, init);
    CORE_pivot_update(m, n, ipiv, indices, offset, init);
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_pivot_update(Quark *quark, Quark_Task_Flags *task_flags,
                             int m, int n, int *ipiv, int *indices,
                             int offset, int init)
{
    DAG_SET_PROPERTIES( "PIV_UP"  , "white"   );
    QUARK_Insert_Task(quark, CORE_pivot_update_quark, task_flags,
        sizeof(int),   &m,       VALUE,
        sizeof(int),   &n,       VALUE,
        sizeof(int)*n,  ipiv,        INOUT,
        sizeof(int)*m,  indices,     INOUT,
        sizeof(int),   &offset,  VALUE,
        sizeof(int),   &init,    VALUE,
        0);
}



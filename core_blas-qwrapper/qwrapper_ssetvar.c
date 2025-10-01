/**
 *
 * @file qwrapper_ssetvar.c
 *
 *  PLASMA core_blas quark wrapper
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Mark Gates
 * @date 2010-11-15
 * @generated s Fri Apr  1 11:02:43 2016
 *
 **/
#include "common.h"

/***************************************************************************//**
 *
 * QUARK_CORE_ssetvar sets a single variable, x := alpha.
 * Since x can be in the middle of a tile, we need to depend on the whole tile,
 * so add xlock argument.
 *
 *******************************************************************************
 *
 * @param[in] alpha
 *         Scalar to set x to, passed by pointer so it can depend on runtime value.
 *
 * @param[out] x
 *         On exit, x = alpha.
 *
 * @param[out] xlock
 *         Pointer to tile owning output variable x.
 *
 **/
void QUARK_CORE_ssetvar(Quark *quark, Quark_Task_Flags *task_flags,
                        const float *alpha, float *x,
                        float *xlock)
{
    DAG_SET_PROPERTIES("setvar", "orange");
    QUARK_Insert_Task(quark, CORE_ssetvar_quark, task_flags,
        sizeof(float),  alpha,  INPUT,
        sizeof(float),  x,              NODEP,  /* INOUT; see xlock */
        sizeof(float),  xlock,          INOUT,  /* tile dependency containing x */
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_ssetvar_quark = PCORE_ssetvar_quark
#define CORE_ssetvar_quark PCORE_ssetvar_quark
#endif
void CORE_ssetvar_quark(Quark *quark)
{
    const float *alpha;
    float *x;

    quark_unpack_args_2( quark, alpha, x );
    *x = *alpha;
}

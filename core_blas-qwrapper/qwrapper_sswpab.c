/**
 *
 * @file qwrapper_sswpab.c
 *
 *  PLASMA InPlaceTransformation module
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 *  This work is the implementation of an inplace transformation
 *  based on the GKK algorithm by Gustavson, Karlsson, Kagstrom
 *  and its fortran implementation.
 *
 * @version 2.8.0
 * @author Mathieu Faverge
 * @date 2010-11-15
 *
 * @generated s Fri Apr  1 11:02:41 2016
 *
 **/
#include "common.h"

/***************************************************************************//**
 *
 **/
void QUARK_CORE_sswpab(Quark *quark, Quark_Task_Flags *task_flags,
                       int i, int n1, int n2,
                       float *A, int szeA)
{
    plasma_profile_by_kernel( task_flags, SWPAB );

    QUARK_Insert_Task(
        quark, CORE_sswpab_quark, task_flags,
        sizeof(int),                           &i,   VALUE,
        sizeof(int),                           &n1,  VALUE,
        sizeof(int),                           &n2,  VALUE,
        sizeof(float)*szeA,       A,            INOUT,
        sizeof(float)*min(n1,n2), NULL,         SCRATCH,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_sswpab_quark = PCORE_sswpab_quark
#define CORE_sswpab_quark PCORE_sswpab_quark
#endif
void CORE_sswpab_quark(Quark *quark)
{
    int i;
    int n1;
    int n2;
    float *A;
    float *work;

    quark_unpack_args_5(quark, i, n1, n2, A, work);
    CORE_sswpab( i, n1, n2, A, work);
}

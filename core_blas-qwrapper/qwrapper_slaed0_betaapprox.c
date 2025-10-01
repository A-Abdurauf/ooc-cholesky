/**
 *
 * @file qwrapper_slaed0_betaapprox.c
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

void
CORE_slaed0_lascl_quark(Quark *quark);

/***************************************************************************//**
 *
 **/
void QUARK_CORE_slaed0_lascl( Quark *quark, Quark_Task_Flags *task_flags,
                              int n, float *scale, float *D, float *E)
{
    plasma_profile_by_kernel( task_flags, LASCL );

    QUARK_Insert_Task(quark, CORE_slaed0_lascl_quark, task_flags,
        sizeof(int),          &n,     VALUE,
        sizeof(float),        scale,     INOUT,
        sizeof(float)*n,      D,         INOUT,
        sizeof(float)*(n-1),  E,         INOUT,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_slaed0_lascl_quark = PCORE_slaed0_lascl_quark
#define CORE_slaed0_lascl_quark PCORE_slaed0_lascl_quark
#endif
void CORE_slaed0_lascl_quark(Quark *quark)
{
    int n;
    float *scale, *D, *E;

    quark_unpack_args_4(quark, n, scale, D, E);

    *scale = PLASMA_FCALL(slanst, SLANST)( &(lapack_const(PlasmaMaxNorm)), &n, D, E);

    CORE_slascl(PlasmaGeneral, 0, 0, *scale, 1., n,   1, D, n  );
    CORE_slascl(PlasmaGeneral, 0, 0, *scale, 1., n-1, 1, E, n-1);
}

void
CORE_slaed0_betaapprox_quark(Quark *quark);

/***************************************************************************//**
 *
 **/
void QUARK_CORE_slaed0_betaapprox(Quark *quark, Quark_Task_Flags *task_flags,
                                  int subpbs, const int *subpbs_displ,
                                  float *D, const float *E)
{
    plasma_gendag_by_kernel( task_flags, LAED0_BETAAPPROX );

    QUARK_Insert_Task(
        quark, CORE_slaed0_betaapprox_quark, task_flags,
        sizeof(int),     &subpbs,      VALUE,
        sizeof(int*),    subpbs_displ, INPUT,
        sizeof(float*), D,            INOUT,
        sizeof(float*), E,            INPUT,
        0);
}

/***************************************************************************//**
 *
 **/
void
CORE_slaed0_betaapprox_quark(Quark *quark)
{
    int subpbs;
    int *subpbs_displ;
    float *D;
    const float *E;

    quark_unpack_args_4(quark, subpbs, subpbs_displ, D, E);
    CORE_slaed0_betaapprox(subpbs, subpbs_displ, D, E);
}

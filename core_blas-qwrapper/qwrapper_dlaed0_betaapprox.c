/**
 *
 * @file qwrapper_dlaed0_betaapprox.c
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
#include "common.h"

void
CORE_dlaed0_lascl_quark(Quark *quark);

/***************************************************************************//**
 *
 **/
void QUARK_CORE_dlaed0_lascl( Quark *quark, Quark_Task_Flags *task_flags,
                              int n, double *scale, double *D, double *E)
{
    plasma_profile_by_kernel( task_flags, LASCL );

    QUARK_Insert_Task(quark, CORE_dlaed0_lascl_quark, task_flags,
        sizeof(int),          &n,     VALUE,
        sizeof(double),        scale,     INOUT,
        sizeof(double)*n,      D,         INOUT,
        sizeof(double)*(n-1),  E,         INOUT,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_dlaed0_lascl_quark = PCORE_dlaed0_lascl_quark
#define CORE_dlaed0_lascl_quark PCORE_dlaed0_lascl_quark
#endif
void CORE_dlaed0_lascl_quark(Quark *quark)
{
    int n;
    double *scale, *D, *E;

    quark_unpack_args_4(quark, n, scale, D, E);

    *scale = PLASMA_FCALL(dlanst, DLANST)( &(lapack_const(PlasmaMaxNorm)), &n, D, E);

    CORE_dlascl(PlasmaGeneral, 0, 0, *scale, 1., n,   1, D, n  );
    CORE_dlascl(PlasmaGeneral, 0, 0, *scale, 1., n-1, 1, E, n-1);
}

void
CORE_dlaed0_betaapprox_quark(Quark *quark);

/***************************************************************************//**
 *
 **/
void QUARK_CORE_dlaed0_betaapprox(Quark *quark, Quark_Task_Flags *task_flags,
                                  int subpbs, const int *subpbs_displ,
                                  double *D, const double *E)
{
    plasma_gendag_by_kernel( task_flags, LAED0_BETAAPPROX );

    QUARK_Insert_Task(
        quark, CORE_dlaed0_betaapprox_quark, task_flags,
        sizeof(int),     &subpbs,      VALUE,
        sizeof(int*),    subpbs_displ, INPUT,
        sizeof(double*), D,            INOUT,
        sizeof(double*), E,            INPUT,
        0);
}

/***************************************************************************//**
 *
 **/
void
CORE_dlaed0_betaapprox_quark(Quark *quark)
{
    int subpbs;
    int *subpbs_displ;
    double *D;
    const double *E;

    quark_unpack_args_4(quark, subpbs, subpbs_displ, D, E);
    CORE_dlaed0_betaapprox(subpbs, subpbs_displ, D, E);
}

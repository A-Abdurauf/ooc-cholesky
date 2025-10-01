/**
 *
 * @file qwrapper_dsteqr.c
 *
 *  PLASMA core_blas quark wrapper
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Grégoire Pichon
 * @date 2014-07
 * @generated d Fri Apr  1 11:02:43 2016
 *
 **/
#include "common.h"

void
CORE_dsteqr_quark(Quark *quark);

/***************************************************************************//**
 *
 **/
void QUARK_CORE_dsteqr(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum compz, int n,
                       double *D, double *E,
                       double *Z, int ldz)
{
    plasma_profile_by_kernel( task_flags, STEQR );

    QUARK_Insert_Task(quark, CORE_dsteqr_quark, task_flags,
        sizeof(PLASMA_enum),              &compz,    VALUE,
        sizeof(int),                      &n,        VALUE,
        sizeof(double)*n,                  D,        INOUT,
        sizeof(double)*(n-1),              E,        INOUT,
        sizeof(double)*ldz*n,  Z,        INOUT,
        sizeof(int),                      &ldz,      VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_dsteqr_quark = PCORE_dsteqr_quark
#define CORE_dsteqr_quark PCORE_dsteqr_quark
#endif
void CORE_dsteqr_quark(Quark *quark)
{
    PLASMA_enum compz;
    int n;
    double *D;
    double *E;
    double *Z;
    int ldz;

    quark_unpack_args_6(quark, compz, n, D, E, Z, ldz);
    CORE_dsteqr(compz, n, D, E, Z, ldz, NULL);
}

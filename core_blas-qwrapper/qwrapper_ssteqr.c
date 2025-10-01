/**
 *
 * @file qwrapper_ssteqr.c
 *
 *  PLASMA core_blas quark wrapper
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Grégoire Pichon
 * @date 2014-07
 * @generated s Fri Apr  1 11:02:43 2016
 *
 **/
#include "common.h"

void
CORE_ssteqr_quark(Quark *quark);

/***************************************************************************//**
 *
 **/
void QUARK_CORE_ssteqr(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum compz, int n,
                       float *D, float *E,
                       float *Z, int ldz)
{
    plasma_profile_by_kernel( task_flags, STEQR );

    QUARK_Insert_Task(quark, CORE_ssteqr_quark, task_flags,
        sizeof(PLASMA_enum),              &compz,    VALUE,
        sizeof(int),                      &n,        VALUE,
        sizeof(float)*n,                  D,        INOUT,
        sizeof(float)*(n-1),              E,        INOUT,
        sizeof(float)*ldz*n,  Z,        INOUT,
        sizeof(int),                      &ldz,      VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_ssteqr_quark = PCORE_ssteqr_quark
#define CORE_ssteqr_quark PCORE_ssteqr_quark
#endif
void CORE_ssteqr_quark(Quark *quark)
{
    PLASMA_enum compz;
    int n;
    float *D;
    float *E;
    float *Z;
    int ldz;

    quark_unpack_args_6(quark, compz, n, D, E, Z, ldz);
    CORE_ssteqr(compz, n, D, E, Z, ldz, NULL);
}

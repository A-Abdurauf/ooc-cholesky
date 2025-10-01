/**
 *
 * @file qwrapper_zsteqr.c
 *
 *  PLASMA core_blas quark wrapper
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Grégoire Pichon
 * @date 2014-07
 * @precisions normal z -> c d s
 *
 **/
#include "common.h"

void
CORE_zsteqr_quark(Quark *quark);

/***************************************************************************//**
 *
 **/
void QUARK_CORE_zsteqr(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum compz, int n,
                       double *D, double *E,
                       PLASMA_Complex64_t *Z, int ldz)
{
    plasma_profile_by_kernel( task_flags, STEQR );

    QUARK_Insert_Task(quark, CORE_zsteqr_quark, task_flags,
        sizeof(PLASMA_enum),              &compz,    VALUE,
        sizeof(int),                      &n,        VALUE,
        sizeof(double)*n,                  D,        INOUT,
        sizeof(double)*(n-1),              E,        INOUT,
        sizeof(PLASMA_Complex64_t)*ldz*n,  Z,        INOUT,
        sizeof(int),                      &ldz,      VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_zsteqr_quark = PCORE_zsteqr_quark
#define CORE_zsteqr_quark PCORE_zsteqr_quark
#endif
void CORE_zsteqr_quark(Quark *quark)
{
    PLASMA_enum compz;
    int n;
    double *D;
    double *E;
    PLASMA_Complex64_t *Z;
    int ldz;

    quark_unpack_args_6(quark, compz, n, D, E, Z, ldz);
    CORE_zsteqr(compz, n, D, E, Z, ldz, NULL);
}

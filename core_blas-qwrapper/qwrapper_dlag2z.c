/**
 *
 * @file qwrapper_dlag2z.c
 *
 *  PLASMA core_blas quark wrapper
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Gregoire Pichon
 * @author Azzam Haidar
 * @date 2014-07
 * @precisions normal z -> c
 *
 **/
#ifdef PLASMA_WITH_MKL
#include <mkl_lapacke.h>
#else
#include <lapacke.h>
#endif
#include "common.h"

void
CORE_dlag2z_quark(Quark *quark);

/***************************************************************************//**
 *
 **/
void QUARK_CORE_dlag2z(Quark *quark, Quark_Task_Flags *task_flags,
                       int m, int n,
                       const double *R, int LDR,
                       PLASMA_Complex64_t *Z, int LDZ )
{
    QUARK_Insert_Task(quark, CORE_dlag2z_quark, task_flags,
        sizeof(int),                        &m,       VALUE,
        sizeof(int),                        &n,       VALUE,
        sizeof(double*),                    R,        INPUT,
        sizeof(int),                        &LDR,     VALUE,
        sizeof(PLASMA_Complex64_t*),        Z,        OUTPUT,
        sizeof(int),                        &LDZ,     VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_dlag2z_quark = PCORE_dlag2z_quark
#define CORE_dlag2z_quark PCORE_dlag2z_quark
#endif
void CORE_dlag2z_quark(Quark *quark)
{
    int m;
    int n;
    const double *R;
    int LDR;
    PLASMA_Complex64_t *Z;
    int LDZ;

    quark_unpack_args_6(quark, m, n, R, LDR, Z, LDZ);

    CORE_dlag2z(m, n, R, LDR, Z, LDZ);
}

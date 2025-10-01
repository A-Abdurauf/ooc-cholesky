/**
 *
 * @file qwrapper_slag2c.c
 *
 *  PLASMA core_blas quark wrapper
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Gregoire Pichon
 * @author Azzam Haidar
 * @date 2014-07
 * @generated c Fri Apr  1 11:02:43 2016
 *
 **/
#ifdef PLASMA_WITH_MKL
#include <mkl_lapacke.h>
#else
#include <lapacke.h>
#endif
#include "common.h"

void
CORE_slag2c_quark(Quark *quark);

/***************************************************************************//**
 *
 **/
void QUARK_CORE_slag2c(Quark *quark, Quark_Task_Flags *task_flags,
                       int m, int n,
                       const float *R, int LDR,
                       PLASMA_Complex32_t *Z, int LDZ )
{
    QUARK_Insert_Task(quark, CORE_slag2c_quark, task_flags,
        sizeof(int),                        &m,       VALUE,
        sizeof(int),                        &n,       VALUE,
        sizeof(float*),                    R,        INPUT,
        sizeof(int),                        &LDR,     VALUE,
        sizeof(PLASMA_Complex32_t*),        Z,        OUTPUT,
        sizeof(int),                        &LDZ,     VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_slag2c_quark = PCORE_slag2c_quark
#define CORE_slag2c_quark PCORE_slag2c_quark
#endif
void CORE_slag2c_quark(Quark *quark)
{
    int m;
    int n;
    const float *R;
    int LDR;
    PLASMA_Complex32_t *Z;
    int LDZ;

    quark_unpack_args_6(quark, m, n, R, LDR, Z, LDZ);

    CORE_slag2c(m, n, R, LDR, Z, LDZ);
}

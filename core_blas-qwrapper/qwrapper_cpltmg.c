/**
 *
 * @file qwrapper_cpltmg.c
 *
 *  PLASMA core_blas quark wrapper
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @generated c Fri Apr  1 11:02:42 2016
 *
 **/
#include "common.h"

/***************************************************************************//**
 *
 **/
void QUARK_CORE_cpltmg( Quark *quark, Quark_Task_Flags *task_flags,
                         PLASMA_enum mtxtype, int m, int n, PLASMA_Complex32_t *A, int lda,
                         int gM, int gN, int m0, int n0, unsigned long long int seed )
{
    plasma_profile_by_kernel( task_flags, PLRNT );

    QUARK_Insert_Task(quark, CORE_cpltmg_quark, task_flags,
        sizeof(int),                      &mtxtype, VALUE,
        sizeof(int),                      &m,       VALUE,
        sizeof(int),                      &n,       VALUE,
        sizeof(PLASMA_Complex32_t)*lda*n, A,            OUTPUT,
        sizeof(int),                      &lda,     VALUE,
        sizeof(int),                      &gM,      VALUE,
        sizeof(int),                      &gN,      VALUE,
        sizeof(int),                      &m0,      VALUE,
        sizeof(int),                      &n0,      VALUE,
        sizeof(unsigned long long int),   &seed,    VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_cpltmg_quark = PCORE_cpltmg_quark
#define CORE_cpltmg_quark PCORE_cpltmg_quark
#endif
void CORE_cpltmg_quark(Quark *quark)
{
    int mtxtype;
    int m;
    int n;
    PLASMA_Complex32_t *A;
    int lda;
    int gM;
    int gN;
    int m0;
    int n0;
    unsigned long long int seed;

    quark_unpack_args_10( quark, mtxtype, m, n, A, lda, gM, gN, m0, n0, seed );
    CORE_cpltmg( mtxtype, m, n, A, lda, gM, gN, m0, n0, seed );
}

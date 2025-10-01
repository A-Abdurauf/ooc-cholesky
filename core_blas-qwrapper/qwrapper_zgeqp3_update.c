/**
 *
 * @file qwrapper_zgeqp3_update.c
 *
 *  PLASMA core_blas quark wrapper
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Mark Gates
 * @date 2010-11-15
 * @precisions normal z -> c d s
 *
 **/
#include "common.h"

/***************************************************************************//**
 *
 **/
void QUARK_CORE_zgeqp3_update( Quark *quark, Quark_Task_Flags *task_flags,
                               PLASMA_Complex64_t *Ajj, int lda1,
                               PLASMA_Complex64_t *Ajk, int lda2,
                               PLASMA_Complex64_t *Fk,  int ldf,
                               int joff, int k, int koff, int nb,
                               double *norms1, double *norms2,
                               int *info )
{
    DAG_SET_PROPERTIES("update", "magenta");
    QUARK_Insert_Task(
        quark, CORE_zgeqp3_update_quark, task_flags,
        sizeof(PLASMA_Complex64_t)*nb*nb,  Ajj,             INPUT,
        sizeof(int),                       &lda1,   VALUE,
        sizeof(PLASMA_Complex64_t)*nb*nb,  Ajk,                     INOUT,
        sizeof(int),                       &lda2,   VALUE,
        sizeof(PLASMA_Complex64_t)*nb*nb,  Fk,              INPUT,
        sizeof(int),                       &ldf,    VALUE,
        sizeof(int),                       &joff,   VALUE,
        sizeof(int),                       &k,      VALUE,
        sizeof(int),                       &koff,   VALUE,
        sizeof(int),                       &nb,     VALUE,
        sizeof(double)*nb,                 norms1,                  INOUT,
        sizeof(double)*nb,                 norms2,                  NODEP,  /* INOUT, but implied by norms1 */
        sizeof(int),                       info,                    OUTPUT,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_zgeqp3_update_quark = PCORE_zgeqp3_update_quark
#define CORE_zgeqp3_update_quark PCORE_zgeqp3_update_quark
#endif
void CORE_zgeqp3_update_quark( Quark *quark )
{
    const PLASMA_Complex64_t *Ajj, *Fk;
    PLASMA_Complex64_t *Ajk;
    int lda1, lda2, ldf, joff, k, koff, nb;
    double *norms1, *norms2;
    int *info;

    quark_unpack_args_13( quark, Ajj, lda1, Ajk, lda2, Fk, ldf, joff, k, koff, nb, norms1, norms2, info );
    CORE_zgeqp3_update(          Ajj, lda1, Ajk, lda2, Fk, ldf, joff, k, koff, nb, norms1, norms2, info );
}

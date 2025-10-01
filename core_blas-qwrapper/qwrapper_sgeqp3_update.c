/**
 *
 * @file qwrapper_sgeqp3_update.c
 *
 *  PLASMA core_blas quark wrapper
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Mark Gates
 * @date 2010-11-15
 * @generated s Fri Apr  1 11:02:44 2016
 *
 **/
#include "common.h"

/***************************************************************************//**
 *
 **/
void QUARK_CORE_sgeqp3_update( Quark *quark, Quark_Task_Flags *task_flags,
                               float *Ajj, int lda1,
                               float *Ajk, int lda2,
                               float *Fk,  int ldf,
                               int joff, int k, int koff, int nb,
                               float *norms1, float *norms2,
                               int *info )
{
    DAG_SET_PROPERTIES("update", "magenta");
    QUARK_Insert_Task(
        quark, CORE_sgeqp3_update_quark, task_flags,
        sizeof(float)*nb*nb,  Ajj,             INPUT,
        sizeof(int),                       &lda1,   VALUE,
        sizeof(float)*nb*nb,  Ajk,                     INOUT,
        sizeof(int),                       &lda2,   VALUE,
        sizeof(float)*nb*nb,  Fk,              INPUT,
        sizeof(int),                       &ldf,    VALUE,
        sizeof(int),                       &joff,   VALUE,
        sizeof(int),                       &k,      VALUE,
        sizeof(int),                       &koff,   VALUE,
        sizeof(int),                       &nb,     VALUE,
        sizeof(float)*nb,                 norms1,                  INOUT,
        sizeof(float)*nb,                 norms2,                  NODEP,  /* INOUT, but implied by norms1 */
        sizeof(int),                       info,                    OUTPUT,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_sgeqp3_update_quark = PCORE_sgeqp3_update_quark
#define CORE_sgeqp3_update_quark PCORE_sgeqp3_update_quark
#endif
void CORE_sgeqp3_update_quark( Quark *quark )
{
    const float *Ajj, *Fk;
    float *Ajk;
    int lda1, lda2, ldf, joff, k, koff, nb;
    float *norms1, *norms2;
    int *info;

    quark_unpack_args_13( quark, Ajj, lda1, Ajk, lda2, Fk, ldf, joff, k, koff, nb, norms1, norms2, info );
    CORE_sgeqp3_update(          Ajj, lda1, Ajk, lda2, Fk, ldf, joff, k, koff, nb, norms1, norms2, info );
}

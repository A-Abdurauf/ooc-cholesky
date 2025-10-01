/**
 *
 * @file qwrapper_clascl.c
 *
 *  PLASMA core_blas quark wrapper
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Grégoire Pichon
 * @author Azzam Haidar
 * @date 2014-07
 * @generated c Fri Apr  1 11:02:43 2016
 *
 **/
#include "common.h"

void
CORE_clascl_quark(Quark *quark);

void
CORE_clascl_p2f1_quark(Quark *quark);

/***************************************************************************//**
 *
 **/
void QUARK_CORE_clascl(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum type, int kl, int ku, float cfrom, float cto,
                       int m, int n, PLASMA_Complex32_t *A, int lda)
{
    plasma_profile_by_kernel( task_flags, LASCL );

    QUARK_Insert_Task(quark, CORE_clascl_quark, task_flags,
                      sizeof(PLASMA_enum),                &type,      VALUE,
                      sizeof(int),                        &kl,        VALUE,
                      sizeof(int),                        &ku,        VALUE,
                      sizeof(float),                     &cfrom,     VALUE,
                      sizeof(float),                     &cto,       VALUE,
                      sizeof(int),                        &m,         VALUE,
                      sizeof(int),                        &n,         VALUE,
                      sizeof(PLASMA_Complex32_t)*lda*n,    A,         INOUT,
                      sizeof(int),                        &lda,       VALUE,
                      0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_clascl_quark = PCORE_clascl_quark
#define CORE_clascl_quark PCORE_clascl_quark
#endif
void CORE_clascl_quark(Quark *quark)
{
    PLASMA_enum type;
    int kl, ku;
    float cfrom, cto;
    int m, n;
    PLASMA_Complex32_t *A;
    int lda;

    quark_unpack_args_9(quark, type, kl, ku, cfrom, cto, m, n, A, lda);
    CORE_clascl(type, kl, ku, cfrom, cto, m, n, A, lda);
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_clascl_p2f1(Quark *quark, Quark_Task_Flags *task_flags,
                            PLASMA_enum type, int kl, int ku, float *cfrom, float *cto,
                            int m, int n, PLASMA_Complex32_t *A, int lda,
                            void *fake, int szefake, int flag)
{
    plasma_profile_by_kernel( task_flags, LASCL );

    QUARK_Insert_Task(quark, CORE_clascl_p2f1_quark, task_flags,
                      sizeof(PLASMA_enum),                &type,      VALUE,
                      sizeof(int),                        &kl,        VALUE,
                      sizeof(int),                        &ku,        VALUE,
                      sizeof(float),                      cfrom,     INPUT,
                      sizeof(float),                      cto,       INPUT,
                      sizeof(int),                        &m,         VALUE,
                      sizeof(int),                        &n,         VALUE,
                      sizeof(PLASMA_Complex32_t)*lda*n,    A,         INOUT,
                      sizeof(int),                        &lda,       VALUE,
                      szefake,                            fake,       flag,
                      0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_clascl_p2f1_quark = PCORE_clascl_p2f1_quark
#define CORE_clascl_p2f1_quark PCORE_clascl_p2f1_quark
#endif
void CORE_clascl_p2f1_quark(Quark *quark)
{
    PLASMA_enum type;
    int kl, ku;
    float *cfrom, *cto;
    int m, n;
    PLASMA_Complex32_t *A;
    int lda;
    void *fake;

    quark_unpack_args_10(quark, type, kl, ku, cfrom, cto, m, n, A, lda, fake);
    CORE_clascl(type, kl, ku, *cfrom, *cto, m, n, A, lda);
}


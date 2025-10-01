/**
 *
 * @file qwrapper_zlaswp.c
 *
 *  PLASMA core_blas quark wrapper
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @precisions normal z -> c d s
 *
 **/
#ifdef PLASMA_WITH_MKL
#include <mkl_lapacke.h>
#else
#include <lapacke.h>
#endif
#include "common.h"

/***************************************************************************//**
 *
 **/
void QUARK_CORE_zlaswp(Quark *quark, Quark_Task_Flags *task_flags,
                       int n, PLASMA_Complex64_t *A, int lda,
                       int i1,  int i2, const int *ipiv, int inc)
{
    plasma_profile_by_kernel( task_flags, LASWP );

    QUARK_Insert_Task(
        quark, CORE_zlaswp_quark, task_flags,
        sizeof(int),                      &n,    VALUE,
        sizeof(PLASMA_Complex64_t)*lda*n,  A,        INOUT | LOCALITY,
        sizeof(int),                      &lda,  VALUE,
        sizeof(int),                      &i1,   VALUE,
        sizeof(int),                      &i2,   VALUE,
        sizeof(int)*n,                     ipiv,     INPUT,
        sizeof(int),                      &inc,  VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_zlaswp_quark = PCORE_zlaswp_quark
#define CORE_zlaswp_quark PCORE_zlaswp_quark
#endif
void CORE_zlaswp_quark(Quark *quark)
{
    int n, lda, i1, i2, inc;
    int *ipiv;
    PLASMA_Complex64_t *A;

    quark_unpack_args_7(quark, n, A, lda, i1, i2, ipiv, inc);
    LAPACKE_zlaswp_work(LAPACK_COL_MAJOR, n, A, lda, i1, i2, ipiv, inc );
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_zlaswp_f2(Quark *quark, Quark_Task_Flags *task_flags,
                          int n, PLASMA_Complex64_t *A, int lda,
                          int i1,  int i2, const int *ipiv, int inc,
                          PLASMA_Complex64_t *fake1, int szefake1, int flag1,
                          PLASMA_Complex64_t *fake2, int szefake2, int flag2)
{
    plasma_profile_by_kernel( task_flags, LASWP );

    QUARK_Insert_Task(
        quark, CORE_zlaswp_f2_quark, task_flags,
        sizeof(int),                        &n,     VALUE,
        sizeof(PLASMA_Complex64_t)*lda*n,    A,         INOUT | LOCALITY,
        sizeof(int),                        &lda,   VALUE,
        sizeof(int),                        &i1,    VALUE,
        sizeof(int),                        &i2,    VALUE,
        sizeof(int)*n,                       ipiv,      INPUT,
        sizeof(int),                        &inc,   VALUE,
        sizeof(PLASMA_Complex64_t)*szefake1, fake1,     flag1,
        sizeof(PLASMA_Complex64_t)*szefake2, fake2,     flag2,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_zlaswp_f2_quark = PCORE_zlaswp_f2_quark
#define CORE_zlaswp_f2_quark PCORE_zlaswp_f2_quark
#endif
void CORE_zlaswp_f2_quark(Quark* quark)
{
    int n, lda, i1, i2, inc;
    int *ipiv;
    PLASMA_Complex64_t *A;
    void *fake1, *fake2;

    quark_unpack_args_9(quark, n, A, lda, i1, i2, ipiv, inc, fake1, fake2);
    LAPACKE_zlaswp_work(LAPACK_COL_MAJOR, n, A, lda, i1, i2, ipiv, inc );
}


/***************************************************************************//**
 *
 **/
void QUARK_CORE_zlaswp_ontile(Quark *quark, Quark_Task_Flags *task_flags,
                              PLASMA_desc descA, PLASMA_Complex64_t *Aij,
                              int i1,  int i2, const int *ipiv, int inc, PLASMA_Complex64_t *fakepanel)
{
    plasma_profile_by_kernel( task_flags, LASWP );

    if (fakepanel == Aij) {
        QUARK_Insert_Task(
            quark, CORE_zlaswp_ontile_quark, task_flags,
            sizeof(PLASMA_desc),              &descA,     VALUE,
            sizeof(PLASMA_Complex64_t)*1,      Aij,           INOUT | LOCALITY,
            sizeof(int),                      &i1,        VALUE,
            sizeof(int),                      &i2,        VALUE,
            sizeof(int)*(i2-i1+1)*abs(inc),   ipiv,           INPUT,
            sizeof(int),                      &inc,       VALUE,
            sizeof(PLASMA_Complex64_t)*1,      fakepanel,     SCRATCH,
            0);
    } else {
        QUARK_Insert_Task(
            quark, CORE_zlaswp_ontile_quark, task_flags,
            sizeof(PLASMA_desc),              &descA,     VALUE,
            sizeof(PLASMA_Complex64_t)*1,      Aij,           INOUT | LOCALITY,
            sizeof(int),                      &i1,        VALUE,
            sizeof(int),                      &i2,        VALUE,
            sizeof(int)*(i2-i1+1)*abs(inc),   ipiv,           INPUT,
            sizeof(int),                      &inc,       VALUE,
            sizeof(PLASMA_Complex64_t)*1,      fakepanel,     INOUT,
            0);
    }
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_zlaswp_ontile_quark = PCORE_zlaswp_ontile_quark
#define CORE_zlaswp_ontile_quark PCORE_zlaswp_ontile_quark
#endif
void CORE_zlaswp_ontile_quark(Quark *quark)
{
    int i1, i2, inc;
    int *ipiv;
    PLASMA_Complex64_t *A, *fake;
    PLASMA_desc descA;

    quark_unpack_args_7(quark, descA, A, i1, i2, ipiv, inc, fake);
    CORE_zlaswp_ontile(descA, i1, i2, ipiv, inc);
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_zlaswp_ontile_f2(Quark *quark, Quark_Task_Flags *task_flags,
                                 PLASMA_desc descA, PLASMA_Complex64_t *Aij,
                                 int i1,  int i2, const int *ipiv, int inc,
                                 PLASMA_Complex64_t *fake1, int szefake1, int flag1,
                                 PLASMA_Complex64_t *fake2, int szefake2, int flag2)
{
    plasma_profile_by_kernel( task_flags, LASWP );

    QUARK_Insert_Task(
        quark, CORE_zlaswp_ontile_f2_quark, task_flags,
        sizeof(PLASMA_desc),                &descA, VALUE,
        sizeof(PLASMA_Complex64_t)*1,        Aij,       INOUT | LOCALITY,
        sizeof(int),                        &i1,    VALUE,
        sizeof(int),                        &i2,    VALUE,
        sizeof(int)*(i2-i1+1)*abs(inc),      ipiv,      INPUT,
        sizeof(int),                        &inc,   VALUE,
        sizeof(PLASMA_Complex64_t)*szefake1, fake1, flag1,
        sizeof(PLASMA_Complex64_t)*szefake2, fake2, flag2,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_zlaswp_ontile_f2_quark = PCORE_zlaswp_ontile_f2_quark
#define CORE_zlaswp_ontile_f2_quark PCORE_zlaswp_ontile_f2_quark
#endif
void CORE_zlaswp_ontile_f2_quark(Quark *quark)
{
    int i1, i2, inc;
    int *ipiv;
    PLASMA_Complex64_t *A;
    PLASMA_desc descA;
    void *fake1, *fake2;

    quark_unpack_args_8(quark, descA, A, i1, i2, ipiv, inc, fake1, fake2);
    CORE_zlaswp_ontile(descA, i1, i2, ipiv, inc);
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_zswptr_ontile(Quark *quark, Quark_Task_Flags *task_flags,
                              PLASMA_desc descA, PLASMA_Complex64_t *Aij,
                              int i1,  int i2, const int *ipiv, int inc,
                              const PLASMA_Complex64_t *Akk, int ldak)
{
    plasma_profile_by_kernel( task_flags, TRSM );

    QUARK_Insert_Task(
        quark, CORE_zswptr_ontile_quark, task_flags,
        sizeof(PLASMA_desc),              &descA, VALUE,
        sizeof(PLASMA_Complex64_t)*1,      Aij,       INOUT | LOCALITY,
        sizeof(int),                      &i1,    VALUE,
        sizeof(int),                      &i2,    VALUE,
        sizeof(int)*(i2-i1+1)*abs(inc),    ipiv,      INPUT,
        sizeof(int),                      &inc,   VALUE,
        sizeof(PLASMA_Complex64_t)*ldak,   Akk,       INPUT,
        sizeof(int),                      &ldak,  VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_zswptr_ontile_quark = PCORE_zswptr_ontile_quark
#define CORE_zswptr_ontile_quark PCORE_zswptr_ontile_quark
#endif
void CORE_zswptr_ontile_quark(Quark *quark)
{
    int i1, i2, inc, ldak;
    int *ipiv;
    PLASMA_Complex64_t *A, *Akk;
    PLASMA_desc descA;

    quark_unpack_args_8(quark, descA, A, i1, i2, ipiv, inc, Akk, ldak);
    CORE_zswptr_ontile(descA, i1, i2, ipiv, inc, Akk, ldak);
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_zlaswpc_ontile(Quark *quark, Quark_Task_Flags *task_flags,
                              PLASMA_desc descA, PLASMA_Complex64_t *Aij,
                              int i1,  int i2, const int *ipiv, int inc, PLASMA_Complex64_t *fakepanel)
{
    plasma_profile_by_kernel( task_flags, LASWP );

    if (fakepanel == Aij) {
        QUARK_Insert_Task(
            quark, CORE_zlaswpc_ontile_quark, task_flags,
            sizeof(PLASMA_desc),              &descA,     VALUE,
            sizeof(PLASMA_Complex64_t)*1,      Aij,           INOUT | LOCALITY,
            sizeof(int),                      &i1,        VALUE,
            sizeof(int),                      &i2,        VALUE,
            sizeof(int)*(i2-i1+1)*abs(inc),   ipiv,           INPUT,
            sizeof(int),                      &inc,       VALUE,
            sizeof(PLASMA_Complex64_t)*1,      fakepanel,     SCRATCH,
            0);
    } else {
        QUARK_Insert_Task(
            quark, CORE_zlaswpc_ontile_quark, task_flags,
            sizeof(PLASMA_desc),              &descA,     VALUE,
            sizeof(PLASMA_Complex64_t)*1,      Aij,           INOUT | LOCALITY,
            sizeof(int),                      &i1,        VALUE,
            sizeof(int),                      &i2,        VALUE,
            sizeof(int)*(i2-i1+1)*abs(inc),   ipiv,           INPUT,
            sizeof(int),                      &inc,       VALUE,
            sizeof(PLASMA_Complex64_t)*1,      fakepanel,     INOUT,
            0);
    }
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_zlaswpc_ontile_quark = PCORE_zlaswpc_ontile_quark
#define CORE_zlaswpc_ontile_quark PCORE_zlaswpc_ontile_quark
#endif
void CORE_zlaswpc_ontile_quark(Quark *quark)
{
    int i1, i2, inc;
    int *ipiv;
    PLASMA_Complex64_t *A, *fake;
    PLASMA_desc descA;

    quark_unpack_args_7(quark, descA, A, i1, i2, ipiv, inc, fake);
    CORE_zlaswpc_ontile(descA, i1, i2, ipiv, inc);
}


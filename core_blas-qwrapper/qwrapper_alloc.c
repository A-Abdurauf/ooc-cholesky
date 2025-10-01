/**
 *
 * @file qwrapper_alloc.c
 *
 *  PLASMA core_blas quark wrapper
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Mathieu Faverge
 * @date 2010-11-15
 *
 **/
#include <stdlib.h>
#include "common.h"

/***************************************************************************//**
 *
 **/
void
CORE_free_quark(Quark *quark)
{
    void *A;

    quark_unpack_args_1(quark, A);
    if (A != NULL)
        free(A);
}

/***************************************************************************//**
 *
 **/
void
QUARK_CORE_free(Quark *quark, Quark_Task_Flags *task_flags, void *A, int szeA)
{
    QUARK_Insert_Task(
        quark, CORE_free_quark, task_flags,
        szeA, A, INOUT,
        0);
}

void
CORE_foo_quark(Quark *quark) {
    void *A;
    quark_unpack_args_1(quark, A);
}

void
CORE_foo2_quark(Quark *quark) {
    void *A, *B;
    quark_unpack_args_2(quark, A, B);
}


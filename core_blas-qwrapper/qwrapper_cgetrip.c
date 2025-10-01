/**
 *
 * @file qwrapper_cgetrip.c
 *
 *  PLASMA InPlaceTransformation module
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 *  This work is the implementation of an inplace transformation
 *  based on the GKK algorithm by Gustavson, Karlsson, Kagstrom
 *  and its fortran implementation.
 *
 * @version 2.8.0
 * @author Mathieu Faverge
 * @date 2010-11-15
 *
 * @generated c Fri Apr  1 11:02:41 2016
 *
 **/
#include "common.h"

/***************************************************************************//**
 *
 **/
void QUARK_CORE_cgetrip(Quark *quark, Quark_Task_Flags *task_flags,
                        int m, int n, PLASMA_Complex32_t *A, int szeA)
{
    plasma_profile_by_kernel( task_flags, GETRIP );

    QUARK_Insert_Task(quark, CORE_cgetrip_quark, task_flags,
        sizeof(int),                     &m,   VALUE,
        sizeof(int),                     &n,   VALUE,
        sizeof(PLASMA_Complex32_t)*szeA, A,        INOUT,
        sizeof(PLASMA_Complex32_t)*szeA, NULL,     SCRATCH,
        0);
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_cgetrip_f1(Quark *quark, Quark_Task_Flags *task_flags,
                           int m, int n,
                           PLASMA_Complex32_t *A,    int szeA,
                           PLASMA_Complex32_t *fake, int szeF, int paramF)
{
    plasma_profile_by_kernel( task_flags, GETRIP );

    if ( (fake == A) && (paramF & GATHERV) ) {
        QUARK_Insert_Task(
            quark, CORE_cgetrip_quark, task_flags,
            sizeof(int),                     &m,   VALUE,
            sizeof(int),                     &n,   VALUE,
            sizeof(PLASMA_Complex32_t)*szeA, A,        INOUT | paramF,
            sizeof(PLASMA_Complex32_t)*szeA, NULL,     SCRATCH,
            0);
    } else {
        QUARK_Insert_Task(
            quark, CORE_cgetrip_f1_quark, task_flags,
            sizeof(int),                     &m,   VALUE,
            sizeof(int),                     &n,   VALUE,
            sizeof(PLASMA_Complex32_t)*szeA, A,        INOUT,
            sizeof(PLASMA_Complex32_t)*szeA, NULL,     SCRATCH,
            sizeof(PLASMA_Complex32_t)*szeF, fake,     paramF,
            0);
    }
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_cgetrip_f2(Quark *quark, Quark_Task_Flags *task_flags,
                           int m, int n,
                           PLASMA_Complex32_t *A,    int szeA,
                           PLASMA_Complex32_t *fake1, int szeF1, int paramF1,
                           PLASMA_Complex32_t *fake2, int szeF2, int paramF2)
{
    plasma_profile_by_kernel( task_flags, GETRIP );

    if ( (fake2 == A) && (paramF2 & GATHERV) ) {
        QUARK_Insert_Task(
            quark, CORE_cgetrip_f1_quark, task_flags,
            sizeof(int),                     &m,   VALUE,
            sizeof(int),                     &n,   VALUE,
            sizeof(PLASMA_Complex32_t)*szeA, A,        INOUT | paramF2,
            sizeof(PLASMA_Complex32_t)*szeA, NULL,     SCRATCH,
            sizeof(PLASMA_Complex32_t)*szeF1, fake1,     paramF1,
            0);
    } else if ( (fake1 == A) && (paramF1 & GATHERV) ) {
        QUARK_Insert_Task(
            quark, CORE_cgetrip_f1_quark, task_flags,
            sizeof(int),                     &m,   VALUE,
            sizeof(int),                     &n,   VALUE,
            sizeof(PLASMA_Complex32_t)*szeA, A,        INOUT | paramF1,
            sizeof(PLASMA_Complex32_t)*szeA, NULL,     SCRATCH,
            sizeof(PLASMA_Complex32_t)*szeF2, fake2,     paramF2,
            0);
    } else {
        QUARK_Insert_Task(
            quark, CORE_cgetrip_f2_quark, task_flags,
            sizeof(int),                     &m,   VALUE,
            sizeof(int),                     &n,   VALUE,
            sizeof(PLASMA_Complex32_t)*szeA, A,        INOUT,
            sizeof(PLASMA_Complex32_t)*szeA, NULL,     SCRATCH,
            sizeof(PLASMA_Complex32_t)*szeF1, fake1,     paramF1,
            sizeof(PLASMA_Complex32_t)*szeF2, fake2,     paramF2,
            0);
    }
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_cgetrip_quark = PCORE_cgetrip_quark
#define CORE_cgetrip_quark PCORE_cgetrip_quark
#endif
void CORE_cgetrip_quark(Quark *quark)
{
    int m;
    int n;
    PLASMA_Complex32_t *A;
    PLASMA_Complex32_t *W;

    quark_unpack_args_4(quark, m, n, A, W);
    CORE_cgetrip(m, n, A, W);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_cgetrip_f1_quark = PCORE_cgetrip_f1_quark
#define CORE_cgetrip_f1_quark PCORE_cgetrip_f1_quark
#endif
void CORE_cgetrip_f1_quark(Quark *quark)
{
    int m;
    int n;
    PLASMA_Complex32_t *A;
    PLASMA_Complex32_t *W;
    PLASMA_Complex32_t *fake;

    quark_unpack_args_5(quark, m, n, A, W, fake);
    CORE_cgetrip(m, n, A, W);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_cgetrip_f2_quark = PCORE_cgetrip_f2_quark
#define CORE_cgetrip_f2_quark PCORE_cgetrip_f2_quark
#endif
void CORE_cgetrip_f2_quark(Quark *quark)
{
    int m;
    int n;
    PLASMA_Complex32_t *A;
    PLASMA_Complex32_t *W;
    PLASMA_Complex32_t *fake1;
    PLASMA_Complex32_t *fake2;

    quark_unpack_args_6(quark, m, n, A, W, fake1, fake2);
    CORE_cgetrip(m, n, A, W);
}

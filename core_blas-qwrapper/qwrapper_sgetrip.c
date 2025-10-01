/**
 *
 * @file qwrapper_sgetrip.c
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
 * @generated s Fri Apr  1 11:02:41 2016
 *
 **/
#include "common.h"

/***************************************************************************//**
 *
 **/
void QUARK_CORE_sgetrip(Quark *quark, Quark_Task_Flags *task_flags,
                        int m, int n, float *A, int szeA)
{
    plasma_profile_by_kernel( task_flags, GETRIP );

    QUARK_Insert_Task(quark, CORE_sgetrip_quark, task_flags,
        sizeof(int),                     &m,   VALUE,
        sizeof(int),                     &n,   VALUE,
        sizeof(float)*szeA, A,        INOUT,
        sizeof(float)*szeA, NULL,     SCRATCH,
        0);
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_sgetrip_f1(Quark *quark, Quark_Task_Flags *task_flags,
                           int m, int n,
                           float *A,    int szeA,
                           float *fake, int szeF, int paramF)
{
    plasma_profile_by_kernel( task_flags, GETRIP );

    if ( (fake == A) && (paramF & GATHERV) ) {
        QUARK_Insert_Task(
            quark, CORE_sgetrip_quark, task_flags,
            sizeof(int),                     &m,   VALUE,
            sizeof(int),                     &n,   VALUE,
            sizeof(float)*szeA, A,        INOUT | paramF,
            sizeof(float)*szeA, NULL,     SCRATCH,
            0);
    } else {
        QUARK_Insert_Task(
            quark, CORE_sgetrip_f1_quark, task_flags,
            sizeof(int),                     &m,   VALUE,
            sizeof(int),                     &n,   VALUE,
            sizeof(float)*szeA, A,        INOUT,
            sizeof(float)*szeA, NULL,     SCRATCH,
            sizeof(float)*szeF, fake,     paramF,
            0);
    }
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_sgetrip_f2(Quark *quark, Quark_Task_Flags *task_flags,
                           int m, int n,
                           float *A,    int szeA,
                           float *fake1, int szeF1, int paramF1,
                           float *fake2, int szeF2, int paramF2)
{
    plasma_profile_by_kernel( task_flags, GETRIP );

    if ( (fake2 == A) && (paramF2 & GATHERV) ) {
        QUARK_Insert_Task(
            quark, CORE_sgetrip_f1_quark, task_flags,
            sizeof(int),                     &m,   VALUE,
            sizeof(int),                     &n,   VALUE,
            sizeof(float)*szeA, A,        INOUT | paramF2,
            sizeof(float)*szeA, NULL,     SCRATCH,
            sizeof(float)*szeF1, fake1,     paramF1,
            0);
    } else if ( (fake1 == A) && (paramF1 & GATHERV) ) {
        QUARK_Insert_Task(
            quark, CORE_sgetrip_f1_quark, task_flags,
            sizeof(int),                     &m,   VALUE,
            sizeof(int),                     &n,   VALUE,
            sizeof(float)*szeA, A,        INOUT | paramF1,
            sizeof(float)*szeA, NULL,     SCRATCH,
            sizeof(float)*szeF2, fake2,     paramF2,
            0);
    } else {
        QUARK_Insert_Task(
            quark, CORE_sgetrip_f2_quark, task_flags,
            sizeof(int),                     &m,   VALUE,
            sizeof(int),                     &n,   VALUE,
            sizeof(float)*szeA, A,        INOUT,
            sizeof(float)*szeA, NULL,     SCRATCH,
            sizeof(float)*szeF1, fake1,     paramF1,
            sizeof(float)*szeF2, fake2,     paramF2,
            0);
    }
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_sgetrip_quark = PCORE_sgetrip_quark
#define CORE_sgetrip_quark PCORE_sgetrip_quark
#endif
void CORE_sgetrip_quark(Quark *quark)
{
    int m;
    int n;
    float *A;
    float *W;

    quark_unpack_args_4(quark, m, n, A, W);
    CORE_sgetrip(m, n, A, W);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_sgetrip_f1_quark = PCORE_sgetrip_f1_quark
#define CORE_sgetrip_f1_quark PCORE_sgetrip_f1_quark
#endif
void CORE_sgetrip_f1_quark(Quark *quark)
{
    int m;
    int n;
    float *A;
    float *W;
    float *fake;

    quark_unpack_args_5(quark, m, n, A, W, fake);
    CORE_sgetrip(m, n, A, W);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_sgetrip_f2_quark = PCORE_sgetrip_f2_quark
#define CORE_sgetrip_f2_quark PCORE_sgetrip_f2_quark
#endif
void CORE_sgetrip_f2_quark(Quark *quark)
{
    int m;
    int n;
    float *A;
    float *W;
    float *fake1;
    float *fake2;

    quark_unpack_args_6(quark, m, n, A, W, fake1, fake2);
    CORE_sgetrip(m, n, A, W);
}

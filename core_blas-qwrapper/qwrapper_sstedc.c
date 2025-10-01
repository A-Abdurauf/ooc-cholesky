/**
 *
 * @file qwrapper_sstedc.c
 *
 *  PLASMA core_blas quark wrapper
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Gregoire Pichon
 * @author Azzam Haidar
 * @date 2014-07
 * @generated s Fri Apr  1 11:02:43 2016
 *
 **/
#include "common.h"

#define REAL
#undef COMPLEX

void
CORE_sstedc_quark(Quark *quark);
void
CORE_sstedc_f2_quark(Quark *quark);

/***************************************************************************//**
 *
 **/
void QUARK_CORE_sstedc(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum compz, int n,
                       float *D, float *E,
                       float *Z, int ldz)
{
    plasma_profile_by_kernel( task_flags, STEDC );

    QUARK_Insert_Task(quark, CORE_sstedc_quark, task_flags,
        sizeof(PLASMA_enum),              &compz,    VALUE,
        sizeof(int),                      &n,        VALUE,
        sizeof(float)*n,                  D,        INOUT,
        sizeof(float)*(n-1),              E,        INOUT,
        sizeof(float)*ldz*n,  Z,        INOUT,
        sizeof(int),                      &ldz,      VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_sstedc_quark = PCORE_sstedc_quark
#define CORE_sstedc_quark PCORE_sstedc_quark
#endif
void CORE_sstedc_quark(Quark *quark)
{
    PLASMA_enum compz;
    int n;
    float *D;
    float *E;
    float *Z;
    int ldz;

    quark_unpack_args_6(quark, compz, n, D, E, Z, ldz);
    CORE_sstedc(compz, n, D, E, Z, ldz,
                NULL, -1,
#ifdef COMPLEX
                NULL, -1,
#endif
                NULL, -1);
}


/***************************************************************************//**
 *
 **/
void QUARK_CORE_sstedc_f2(Quark *quark, Quark_Task_Flags *task_flags,
                          PLASMA_enum compz, int n,
                          float *D, float *E,
                          float *Z, int ldz,
                          void *fake1, int szefake1, int flag1,
                          void *fake2, int szefake2, int flag2)
{
    plasma_profile_by_kernel( task_flags, STEDC );


    if ( D == fake2 ) {
        QUARK_Insert_Task(quark, CORE_sstedc_f2_quark, task_flags,
            sizeof(PLASMA_enum),              &compz,    VALUE,
            sizeof(int),                      &n,        VALUE,
            sizeof(float)*n,                  D,        INPUT,
            sizeof(float)*(n-1),              E,        NODEP,
            sizeof(float)*ldz*n,  Z,        NODEP,
            sizeof(int),                      &ldz,      VALUE,
            szefake1,                          fake1,    flag1,
            1,                                 NULL,     NODEP,
            0);
    }
    else {
        QUARK_Insert_Task(quark, CORE_sstedc_f2_quark, task_flags,
            sizeof(PLASMA_enum),              &compz,    VALUE,
            sizeof(int),                      &n,        VALUE,
            sizeof(float)*n,                  D,        NODEP,
            sizeof(float)*(n-1),              E,        NODEP,
            sizeof(float)*ldz*n,  Z,        NODEP,
            sizeof(int),                      &ldz,      VALUE,
            szefake1,                          fake1,    flag1,
            szefake2,                          fake2,    flag2,
            0);
    }
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_sstedc_f2_quark = PCORE_sstedc_f2_quark
#define CORE_sstedc_f2_quark PCORE_sstedc_f2_quark
#endif
void CORE_sstedc_f2_quark(Quark *quark)
{
    PLASMA_enum compz;
    int n;
    float *D;
    float *E;
    float *Z;
    int ldz;
    void *fake1, *fake2;

    quark_unpack_args_8(quark, compz, n, D, E, Z, ldz, fake1, fake2);
    CORE_sstedc(compz, n, D, E, Z, ldz,
                NULL, -1,
#ifdef COMPLEX
                NULL, -1,
#endif
                NULL, -1);
}

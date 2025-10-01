/**
 *
 * @file qwrapper_zstedc.c
 *
 *  PLASMA core_blas quark wrapper
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Gregoire Pichon
 * @author Azzam Haidar
 * @date 2014-07
 * @precisions normal z -> c d s
 *
 **/
#include "common.h"

#define COMPLEX
#undef REAL

void
CORE_zstedc_quark(Quark *quark);
void
CORE_zstedc_f2_quark(Quark *quark);

/***************************************************************************//**
 *
 **/
void QUARK_CORE_zstedc(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum compz, int n,
                       double *D, double *E,
                       PLASMA_Complex64_t *Z, int ldz)
{
    plasma_profile_by_kernel( task_flags, STEDC );

    QUARK_Insert_Task(quark, CORE_zstedc_quark, task_flags,
        sizeof(PLASMA_enum),              &compz,    VALUE,
        sizeof(int),                      &n,        VALUE,
        sizeof(double)*n,                  D,        INOUT,
        sizeof(double)*(n-1),              E,        INOUT,
        sizeof(PLASMA_Complex64_t)*ldz*n,  Z,        INOUT,
        sizeof(int),                      &ldz,      VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_zstedc_quark = PCORE_zstedc_quark
#define CORE_zstedc_quark PCORE_zstedc_quark
#endif
void CORE_zstedc_quark(Quark *quark)
{
    PLASMA_enum compz;
    int n;
    double *D;
    double *E;
    PLASMA_Complex64_t *Z;
    int ldz;

    quark_unpack_args_6(quark, compz, n, D, E, Z, ldz);
    CORE_zstedc(compz, n, D, E, Z, ldz,
                NULL, -1,
#ifdef COMPLEX
                NULL, -1,
#endif
                NULL, -1);
}


/***************************************************************************//**
 *
 **/
void QUARK_CORE_zstedc_f2(Quark *quark, Quark_Task_Flags *task_flags,
                          PLASMA_enum compz, int n,
                          double *D, double *E,
                          PLASMA_Complex64_t *Z, int ldz,
                          void *fake1, int szefake1, int flag1,
                          void *fake2, int szefake2, int flag2)
{
    plasma_profile_by_kernel( task_flags, STEDC );


    if ( D == fake2 ) {
        QUARK_Insert_Task(quark, CORE_zstedc_f2_quark, task_flags,
            sizeof(PLASMA_enum),              &compz,    VALUE,
            sizeof(int),                      &n,        VALUE,
            sizeof(double)*n,                  D,        INPUT,
            sizeof(double)*(n-1),              E,        NODEP,
            sizeof(PLASMA_Complex64_t)*ldz*n,  Z,        NODEP,
            sizeof(int),                      &ldz,      VALUE,
            szefake1,                          fake1,    flag1,
            1,                                 NULL,     NODEP,
            0);
    }
    else {
        QUARK_Insert_Task(quark, CORE_zstedc_f2_quark, task_flags,
            sizeof(PLASMA_enum),              &compz,    VALUE,
            sizeof(int),                      &n,        VALUE,
            sizeof(double)*n,                  D,        NODEP,
            sizeof(double)*(n-1),              E,        NODEP,
            sizeof(PLASMA_Complex64_t)*ldz*n,  Z,        NODEP,
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
#pragma weak CORE_zstedc_f2_quark = PCORE_zstedc_f2_quark
#define CORE_zstedc_f2_quark PCORE_zstedc_f2_quark
#endif
void CORE_zstedc_f2_quark(Quark *quark)
{
    PLASMA_enum compz;
    int n;
    double *D;
    double *E;
    PLASMA_Complex64_t *Z;
    int ldz;
    void *fake1, *fake2;

    quark_unpack_args_8(quark, compz, n, D, E, Z, ldz, fake1, fake2);
    CORE_zstedc(compz, n, D, E, Z, ldz,
                NULL, -1,
#ifdef COMPLEX
                NULL, -1,
#endif
                NULL, -1);
}

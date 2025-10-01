/**
 *
 * @file qwrapper_ctstrf.c
 *
 *  PLASMA core_blas quark wrapper
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Hatem Ltaief
 * @author Mathieu Faverge
 * @author Jakub Kurzak
 * @date 2010-11-15
 * @generated c Fri Apr  1 11:02:40 2016
 *
 **/
#include "common.h"

/***************************************************************************//**
 *
 **/
void QUARK_CORE_ctstrf(Quark *quark, Quark_Task_Flags *task_flags,
                       int m, int n, int ib, int nb,
                       PLASMA_Complex32_t *U, int ldu,
                       PLASMA_Complex32_t *A, int lda,
                       PLASMA_Complex32_t *L, int ldl,
                       int *IPIV,
                       PLASMA_sequence *sequence, PLASMA_request *request,
                       PLASMA_bool check_info, int iinfo)
{
    plasma_profile_by_kernel( task_flags, TSTRF );

    QUARK_Insert_Task(quark, CORE_ctstrf_quark, task_flags,
        sizeof(int),                        &m,             VALUE,
        sizeof(int),                        &n,             VALUE,
        sizeof(int),                        &ib,            VALUE,
        sizeof(int),                        &nb,            VALUE,
        sizeof(PLASMA_Complex32_t)*nb*nb,    U,                     INOUT | QUARK_REGION_D | QUARK_REGION_U,
        sizeof(int),                        &ldu,           VALUE,
        sizeof(PLASMA_Complex32_t)*nb*nb,    A,                     INOUT | LOCALITY,
        sizeof(int),                        &lda,           VALUE,
        sizeof(PLASMA_Complex32_t)*ib*nb,    L,                     OUTPUT,
        sizeof(int),                        &ldl,           VALUE,
        sizeof(int)*nb,                      IPIV,                  OUTPUT,
        sizeof(PLASMA_Complex32_t)*ib*nb,    NULL,                  SCRATCH,
        sizeof(int),                        &nb,            VALUE,
        sizeof(PLASMA_sequence*),           &sequence,      VALUE,
        sizeof(PLASMA_request*),            &request,       VALUE,
        sizeof(PLASMA_bool),                &check_info,    VALUE,
        sizeof(int),                        &iinfo,         VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_ctstrf_quark = PCORE_ctstrf_quark
#define CORE_ctstrf_quark PCORE_ctstrf_quark
#endif
void CORE_ctstrf_quark(Quark *quark)
{
    int m;
    int n;
    int ib;
    int nb;
    PLASMA_Complex32_t *U;
    int ldu;
    PLASMA_Complex32_t *A;
    int lda;
    PLASMA_Complex32_t *L;
    int ldl;
    int *IPIV;
    PLASMA_Complex32_t *WORK;
    int ldwork;
    PLASMA_sequence *sequence;
    PLASMA_request *request;
    PLASMA_bool check_info;
    int iinfo;

    int info;

    quark_unpack_args_17(quark, m, n, ib, nb, U, ldu, A, lda, L, ldl, IPIV, WORK, ldwork, sequence, request, check_info, iinfo);
    CORE_ctstrf(m, n, ib, nb, U, ldu, A, lda, L, ldl, IPIV, WORK, ldwork, &info);
    if (info != PLASMA_SUCCESS && check_info)
        plasma_sequence_flush(quark, sequence, request, iinfo + info);
}

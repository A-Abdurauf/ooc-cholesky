/**
 *
 * @file qwrapper_sgeqp3_tntpiv.c
 *
 *  PLASMA core_blas quark wrapper
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @generated s Fri Apr  1 11:02:44 2016
 *
 **/
#include "common.h"

/***************************************************************************//**
 *
 **/
void QUARK_CORE_sgeqp3_tntpiv(Quark *quark, Quark_Task_Flags *task_flags,
                              int m, int n, int nb,
                              float *A, int lda,
                              int *IPIV,
                              PLASMA_sequence *sequence, PLASMA_request *request,
                              PLASMA_bool check_info, int iinfo)
{
    plasma_profile_by_kernel( task_flags, GEQRT );

    QUARK_Insert_Task(quark, CORE_sgeqp3_tntpiv_quark, task_flags,
        sizeof(int),                        &m,             VALUE,
        sizeof(int),                        &n,             VALUE,
        sizeof(float)*nb*nb,    A,                     INOUT | LOCALITY,
        sizeof(int),                        &lda,           VALUE,
        sizeof(int)*nb,                      IPIV,                  OUTPUT,
        sizeof(float)*min(m,n),    NULL,               SCRATCH,
        sizeof(int)*n,                          NULL,               SCRATCH,
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
#pragma weak CORE_sgeqp3_tntpiv_quark = PCORE_sgeqp3_tntpiv_quark
#define CORE_sgeqp3_tntpiv_quark PCORE_sgeqp3_tntpiv_quark
#endif
void CORE_sgeqp3_tntpiv_quark(Quark *quark)
{
    int m;
    int n;
    float *A;
    int lda;
    int *IPIV;
    int *iwork;
    float *tau;
    PLASMA_sequence *sequence;
    PLASMA_request *request;
    PLASMA_bool check_info;
    int iinfo;
    int info;

    quark_unpack_args_11(quark, m, n, A, lda, IPIV, tau, iwork,
                         sequence, request, check_info, iinfo);
    info = CORE_sgeqp3_tntpiv(m, n, A, lda, IPIV, tau, iwork );
    if (info != PLASMA_SUCCESS) {
        if (check_info || info < 0)
            plasma_sequence_flush(quark, sequence, request, iinfo+info);
        else{
            int i;
            for(i=info-1; i<min(m,n); i++)
                IPIV[i] = i+1;
        }
    }
}

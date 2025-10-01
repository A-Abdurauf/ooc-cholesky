/**
 *
 * @file qwrapper_zgetrf_nopiv.c
 *
 *  PLASMA core_blas quark wrapper
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Omar Zenati
 * @author Mathieu Faverge
 * @date 2013-02-01
 * @precisions normal z -> c d s
 *
 **/
#include "common.h"

/***************************************************************************//**
 *
 **/
void QUARK_CORE_zgetrf_nopiv(Quark *quark, Quark_Task_Flags *task_flags,
                             int m, int n, int ib, int nb,
                             PLASMA_Complex64_t *A, int lda,
                             PLASMA_sequence *sequence, PLASMA_request *request,
                             int iinfo)
{
    plasma_profile_by_kernel( task_flags, GETRF );

    QUARK_Insert_Task(quark, CORE_zgetrf_nopiv_quark, task_flags,
        sizeof(int),                        &m,             VALUE,
        sizeof(int),                        &n,             VALUE,
        sizeof(int),                        &ib,            VALUE,
        sizeof(PLASMA_Complex64_t)*nb*nb,    A,                     INOUT | LOCALITY,
        sizeof(int),                        &lda,           VALUE,
        sizeof(PLASMA_sequence*),           &sequence,      VALUE,
        sizeof(PLASMA_request*),            &request,       VALUE,
        sizeof(int),                        &iinfo,         VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_zgetrf_nopiv_quark = PCORE_zgetrf_nopiv_quark
#define CORE_zgetrf_nopiv_quark PCORE_zgetrf_nopiv_quark
#endif
void CORE_zgetrf_nopiv_quark(Quark *quark)
{
    int m;
    int n;
    int ib;
    PLASMA_Complex64_t *A;
    int lda;
    PLASMA_sequence *sequence;
    PLASMA_request *request;
    int iinfo;
    int info = 0;

    quark_unpack_args_8(quark, m, n, ib, A, lda, sequence, request, iinfo );
    info = CORE_zgetrf_nopiv(m, n, ib, A, lda );
    if (info != PLASMA_SUCCESS) {
        plasma_sequence_flush(quark, sequence, request, iinfo+info);
    }
}

/**
 *
 * @file qwrapper_cgetrf_rectil.c
 *
 *  PLASMA core_blas quark wrapper
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Hatem Ltaief
 * @author Mathieu Faverge
 * @author Piotr Luszczek
 * @date 2009-11-15
 *
 * @generated c Fri Apr  1 11:02:42 2016
 *
 **/
#include "common.h"

/***************************************************************************//**
 *
 **/
void QUARK_CORE_cgetrf_rectil(Quark *quark, Quark_Task_Flags *task_flags,
                              CORE_cgetrf_data_t *data,
                              PLASMA_desc A, PLASMA_Complex32_t *Amn, int size,
                              int *IPIV,
                              PLASMA_sequence *sequence, PLASMA_request *request,
                              PLASMA_bool check_info, int iinfo,
                              int nbthread)
{
    plasma_profile_by_kernel( task_flags, GETRF );

    QUARK_Insert_Task(quark, CORE_cgetrf_rectil_quark, task_flags,
        sizeof(void*),                       data,              INOUT,
        sizeof(PLASMA_desc),                &A,             VALUE,
        sizeof(PLASMA_Complex32_t)*size,     Amn,               INOUT,
        sizeof(int)*A.n,                     IPIV,              OUTPUT,
        sizeof(PLASMA_sequence*),           &sequence,      VALUE,
        sizeof(PLASMA_request*),            &request,       VALUE,
        sizeof(PLASMA_bool),                &check_info,    VALUE,
        sizeof(int),                        &iinfo,         VALUE,
        sizeof(int),                        &nbthread,      VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_cgetrf_rectil_quark = PCORE_cgetrf_rectil_quark
#define CORE_cgetrf_rectil_quark PCORE_cgetrf_rectil_quark
#endif
void CORE_cgetrf_rectil_quark(Quark* quark)
{
    CORE_cgetrf_data_t *data;
    PLASMA_desc A;
    PLASMA_Complex32_t *Amn;
    int *IPIV;
    PLASMA_sequence *sequence;
    PLASMA_request *request;
    PLASMA_bool check_info;
    int iinfo;

    int info[3];
    int maxthreads;

    quark_unpack_args_9(quark, data, A, Amn, IPIV, sequence, request,
                        check_info, iinfo, maxthreads );

    info[1] = QUARK_Get_RankInTask(quark);
    info[2] = maxthreads;

    CORE_cgetrf_rectil( data, A, IPIV, info );
    if (info[1] == 0 && info[0] != PLASMA_SUCCESS && check_info)
        plasma_sequence_flush(quark, sequence, request, iinfo + info[0] );
}

/**
 *
 * @file qwrapper_chegst.c
 *
 *  PLASMA core_blas quark wrapper
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Hatem Ltaief
 * @date 2010-11-15
 * @generated c Fri Apr  1 11:02:42 2016
 *
 **/
#ifdef PLASMA_WITH_MKL
#include <mkl_lapacke.h>
#else
#include <lapacke.h>
#endif
#include "common.h"

#undef REAL
#define COMPLEX

/***************************************************************************//**
 *
 **/
void QUARK_CORE_chegst(Quark *quark, Quark_Task_Flags *task_flags,
                       int itype, PLASMA_enum uplo, int n,
                       PLASMA_Complex32_t *A, int lda,
                       PLASMA_Complex32_t *B, int ldb,
                       PLASMA_sequence *sequence, PLASMA_request *request,
                       int iinfo)
{
    QUARK_Insert_Task(quark, CORE_chegst_quark, task_flags,
        sizeof(int),                        &itype,      VALUE,
        sizeof(PLASMA_enum),                &uplo,      VALUE,
        sizeof(int),                        &n,         VALUE,
        sizeof(PLASMA_Complex32_t)*lda*n,    A,             INOUT,
        sizeof(int),                        &lda,       VALUE,
#ifdef COMPLEX
        sizeof(PLASMA_Complex32_t)*ldb*n,    B,             INOUT,
#else
        sizeof(PLASMA_Complex32_t)*ldb*n,    B,             INPUT,
#endif
        sizeof(int),                        &ldb,       VALUE,
        sizeof(PLASMA_sequence*),           &sequence,  VALUE,
        sizeof(PLASMA_request*),            &request,   VALUE,
        sizeof(int),                        &iinfo,     VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_chegst_quark = PCORE_chegst_quark
#define CORE_chegst_quark PCORE_chegst_quark
#endif
void CORE_chegst_quark(Quark *quark)
{
    int itype;
    PLASMA_enum uplo;
    int n;
    PLASMA_Complex32_t *A;
    int lda;
    PLASMA_Complex32_t *B;
    int ldb;
    PLASMA_sequence *sequence;
    PLASMA_request *request;
    int iinfo;

    int info;

    quark_unpack_args_10(quark, itype, uplo, n, A, lda, B, ldb, sequence, request, iinfo);
    info = LAPACKE_chegst_work(
        LAPACK_COL_MAJOR,
        itype,
        lapack_const(uplo),
        n, A, lda, B, ldb);
    if (sequence->status == PLASMA_SUCCESS && info != 0)
      plasma_sequence_flush(quark, sequence, request, iinfo+info);
}

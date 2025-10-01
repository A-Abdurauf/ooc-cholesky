/**
 *
 * @file qwrapper_sgemm.c
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
 * @generated s Fri Apr  1 11:02:40 2016
 *
 **/
#include <cblas.h>
#include "common.h"

/***************************************************************************//**
 *
 **/
void QUARK_CORE_sgemm(Quark *quark, Quark_Task_Flags *task_flags,
                      PLASMA_enum transA, int transB,
                      int m, int n, int k, int nb,
                      float alpha, const float *A, int lda,
                                                const float *B, int ldb,
                      float beta, float *C, int ldc)
{
    plasma_profile_by_kernel( task_flags, GEMM );

    QUARK_Insert_Task(quark, CORE_sgemm_quark, task_flags,
        sizeof(PLASMA_enum),                &transA,    VALUE,
        sizeof(PLASMA_enum),                &transB,    VALUE,
        sizeof(int),                        &m,         VALUE,
        sizeof(int),                        &n,         VALUE,
        sizeof(int),                        &k,         VALUE,
        sizeof(float),         &alpha,     VALUE,
        sizeof(float)*nb*nb,    A,                 INPUT,
        sizeof(int),                        &lda,       VALUE,
        sizeof(float)*nb*nb,    B,                 INPUT,
        sizeof(int),                        &ldb,       VALUE,
        sizeof(float),         &beta,      VALUE,
        sizeof(float)*nb*nb,    C,                 INOUT,
        sizeof(int),                        &ldc,       VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_sgemm2( Quark *quark, Quark_Task_Flags *task_flags,
                        PLASMA_enum transA, int transB,
                        int m, int n, int k, int nb,
                        float alpha, const float *A, int lda,
                        const float *B, int ldb,
                        float beta, float *C, int ldc)
{
    plasma_profile_by_kernel( task_flags, GEMM );

    QUARK_Insert_Task(quark, CORE_sgemm_quark, task_flags,
        sizeof(PLASMA_enum),                &transA,    VALUE,
        sizeof(PLASMA_enum),                &transB,    VALUE,
        sizeof(int),                        &m,         VALUE,
        sizeof(int),                        &n,         VALUE,
        sizeof(int),                        &k,         VALUE,
        sizeof(float),         &alpha,     VALUE,
        sizeof(float)*nb*nb,    A,                 INPUT,
        sizeof(int),                        &lda,       VALUE,
        sizeof(float)*nb*nb,    B,                 INPUT,
        sizeof(int),                        &ldb,       VALUE,
        sizeof(float),         &beta,      VALUE,
        sizeof(float)*nb*nb,    C,                 INOUT | LOCALITY | GATHERV,
        sizeof(int),                        &ldc,       VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_sgemm_quark = PCORE_sgemm_quark
#define CORE_sgemm_quark PCORE_sgemm_quark
#endif
void CORE_sgemm_quark(Quark *quark)
{
    PLASMA_enum transA;
    PLASMA_enum transB;
    int m;
    int n;
    int k;
    float alpha;
    float *A;
    int lda;
    float *B;
    int ldb;
    float beta;
    float *C;
    int ldc;

    quark_unpack_args_13(quark, transA, transB, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
    cblas_sgemm(
        CblasColMajor,
        (CBLAS_TRANSPOSE)transA, (CBLAS_TRANSPOSE)transB,
        m, n, k,
        (alpha), A, lda,
        B, ldb,
        (beta), C, ldc);
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_sgemm_f2(Quark *quark, Quark_Task_Flags *task_flags,
                         PLASMA_enum transA, int transB,
                         int m, int n, int k, int nb,
                         float alpha, const float *A, int lda,
                                                   const float *B, int ldb,
                         float beta, float *C, int ldc,
                         float *fake1, int szefake1, int flag1,
                         float *fake2, int szefake2, int flag2)
{
    plasma_profile_by_kernel( task_flags, GEMM );

    QUARK_Insert_Task(quark, CORE_sgemm_f2_quark, task_flags,
        sizeof(PLASMA_enum),                &transA,    VALUE,
        sizeof(PLASMA_enum),                &transB,    VALUE,
        sizeof(int),                        &m,         VALUE,
        sizeof(int),                        &n,         VALUE,
        sizeof(int),                        &k,         VALUE,
        sizeof(float),         &alpha,     VALUE,
        sizeof(float)*nb*nb,    A,                 INPUT,
        sizeof(int),                        &lda,       VALUE,
        sizeof(float)*nb*nb,    B,                 INPUT,
        sizeof(int),                        &ldb,       VALUE,
        sizeof(float),         &beta,      VALUE,
        sizeof(float)*nb*nb,    C,                 INOUT | LOCALITY,
        sizeof(int),                        &ldc,       VALUE,
        sizeof(float)*szefake1, fake1,             flag1,
        sizeof(float)*szefake2, fake2,             flag2,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_sgemm_f2_quark = PCORE_sgemm_f2_quark
#define CORE_sgemm_f2_quark PCORE_sgemm_f2_quark
#endif
void CORE_sgemm_f2_quark(Quark* quark)
{
    PLASMA_enum transA;
    PLASMA_enum transB;
    int M;
    int N;
    int K;
    float alpha;
    float *A;
    int LDA;
    float *B;
    int LDB;
    float beta;
    float *C;
    int LDC;
    void *fake1, *fake2;

    quark_unpack_args_15(quark, transA, transB, M, N, K, alpha,
                         A, LDA, B, LDB, beta, C, LDC, fake1, fake2);
    cblas_sgemm(
        CblasColMajor,
        (CBLAS_TRANSPOSE)transA, (CBLAS_TRANSPOSE)transB,
        M, N, K,
        (alpha), A, LDA,
        B, LDB,
        (beta), C, LDC);
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_sgemm_p2(Quark *quark, Quark_Task_Flags *task_flags,
                         PLASMA_enum transA, int transB,
                         int m, int n, int k, int nb,
                         float alpha, const float *A, int lda,
                         const float **B, int ldb,
                         float beta, float *C, int ldc)
{
    plasma_profile_by_kernel( task_flags, GEMM );

    QUARK_Insert_Task(quark, CORE_sgemm_p2_quark, task_flags,
        sizeof(PLASMA_enum),                &transA,    VALUE,
        sizeof(PLASMA_enum),                &transB,    VALUE,
        sizeof(int),                        &m,         VALUE,
        sizeof(int),                        &n,         VALUE,
        sizeof(int),                        &k,         VALUE,
        sizeof(float),         &alpha,     VALUE,
        sizeof(float)*lda*nb,   A,                 INPUT,
        sizeof(int),                        &lda,       VALUE,
        sizeof(float*),         B,                 INPUT,
        sizeof(int),                        &ldb,       VALUE,
        sizeof(float),         &beta,      VALUE,
        sizeof(float)*ldc*nb,    C,                 INOUT | LOCALITY,
        sizeof(int),                        &ldc,       VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_sgemm_p2_quark = PCORE_sgemm_p2_quark
#define CORE_sgemm_p2_quark PCORE_sgemm_p2_quark
#endif
void CORE_sgemm_p2_quark(Quark* quark)
{
    PLASMA_enum transA;
    PLASMA_enum transB;
    int M;
    int N;
    int K;
    float alpha;
    float *A;
    int LDA;
    float **B;
    int LDB;
    float beta;
    float *C;
    int LDC;

    quark_unpack_args_13(quark, transA, transB, M, N, K, alpha,
                         A, LDA, B, LDB, beta, C, LDC);
    cblas_sgemm(
        CblasColMajor,
        (CBLAS_TRANSPOSE)transA, (CBLAS_TRANSPOSE)transB,
        M, N, K,
        (alpha), A, LDA,
        *B, LDB,
        (beta), C, LDC);
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_sgemm_p3(Quark *quark, Quark_Task_Flags *task_flags,
                           PLASMA_enum transA, int transB,
                           int m, int n, int k, int nb,
                           float alpha, const float *A, int lda,
                           const float *B, int ldb,
                           float beta, float **C, int ldc)
{
    plasma_profile_by_kernel( task_flags, GEMM );

    QUARK_Insert_Task(quark, CORE_sgemm_p3_quark, task_flags,
        sizeof(PLASMA_enum),                &transA,    VALUE,
        sizeof(PLASMA_enum),                &transB,    VALUE,
        sizeof(int),                        &m,         VALUE,
        sizeof(int),                        &n,         VALUE,
        sizeof(int),                        &k,         VALUE,
        sizeof(float),         &alpha,     VALUE,
        sizeof(float)*lda*nb,   A,                 INPUT,
        sizeof(int),                        &lda,       VALUE,
        sizeof(float)*ldb*nb,   B,                 INPUT,
        sizeof(int),                        &ldb,       VALUE,
        sizeof(float),         &beta,      VALUE,
        sizeof(float*),         C,                 INOUT | LOCALITY,
        sizeof(int),                        &ldc,       VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_sgemm_p3_quark = PCORE_sgemm_p3_quark
#define CORE_sgemm_p3_quark PCORE_sgemm_p3_quark
#endif
void CORE_sgemm_p3_quark(Quark* quark)
{
    PLASMA_enum transA;
    PLASMA_enum transB;
    int M;
    int N;
    int K;
    float alpha;
    float *A;
    int LDA;
    float *B;
    int LDB;
    float beta;
    float **C;
    int LDC;

    quark_unpack_args_13(quark, transA, transB, M, N, K, alpha,
                         A, LDA, B, LDB, beta, C, LDC);
    cblas_sgemm(
        CblasColMajor,
        (CBLAS_TRANSPOSE)transA, (CBLAS_TRANSPOSE)transB,
        M, N, K,
        (alpha), A, LDA,
        B, LDB,
        (beta), *C, LDC);
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_sgemm_p2f1(Quark *quark, Quark_Task_Flags *task_flags,
                           PLASMA_enum transA, int transB,
                           int m, int n, int k, int nb,
                           float alpha, const float *A, int lda,
                           const float **B, int ldb,
                           float beta, float *C, int ldc,
                           float *fake1, int szefake1, int flag1)
{
    plasma_profile_by_kernel( task_flags, GEMM );

    QUARK_Insert_Task(quark, CORE_sgemm_p2f1_quark, task_flags,
        sizeof(PLASMA_enum),                &transA,    VALUE,
        sizeof(PLASMA_enum),                &transB,    VALUE,
        sizeof(int),                        &m,         VALUE,
        sizeof(int),                        &n,         VALUE,
        sizeof(int),                        &k,         VALUE,
        sizeof(float),         &alpha,     VALUE,
        sizeof(float)*lda*nb,   A,                 INPUT,
        sizeof(int),                        &lda,       VALUE,
        sizeof(float*),         B,                 INPUT,
        sizeof(int),                        &ldb,       VALUE,
        sizeof(float),         &beta,      VALUE,
        sizeof(float)*ldc*nb,    C,                 INOUT | LOCALITY,
        sizeof(int),                        &ldc,       VALUE,
        sizeof(float)*szefake1, fake1,             flag1,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_sgemm_p2f1_quark = PCORE_sgemm_p2f1_quark
#define CORE_sgemm_p2f1_quark PCORE_sgemm_p2f1_quark
#endif
void CORE_sgemm_p2f1_quark(Quark* quark)
{
    PLASMA_enum transA;
    PLASMA_enum transB;
    int M;
    int N;
    int K;
    float alpha;
    float *A;
    int LDA;
    float **B;
    int LDB;
    float beta;
    float *C;
    int LDC;
    void *fake1;

    quark_unpack_args_14(quark, transA, transB, M, N, K, alpha,
                         A, LDA, B, LDB, beta, C, LDC, fake1);
    cblas_sgemm(
        CblasColMajor,
        (CBLAS_TRANSPOSE)transA, (CBLAS_TRANSPOSE)transB,
        M, N, K,
        (alpha), A, LDA,
        *B, LDB,
        (beta), C, LDC);
}

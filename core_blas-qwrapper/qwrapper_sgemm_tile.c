/**
 *
 * @file qwrapper_sgemm_tile.c
 *
 *  PLASMA core_blas quark wrapper
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Mark Gates
 * @date 2010-11-15
 * @generated s Fri Apr  1 11:02:43 2016
 *
 **/
#include <cblas.h>
#include "common.h"

/***************************************************************************//**
 *
 * Version of sgemm for tile storage, to avoid dependency problem when
 * computations are done within the tile. alpha and beta are passed as
 * pointers so they can depend on runtime values.
 *
 * @param[in] Alock
 *          Pointer to tile owning submatrix A.
 *
 * @param[in] Block
 *          Pointer to tile owning submatrix B.
 *
 * @param[in] Clock
 *          Pointer to tile owning submatrix C.
 *
 **/
void QUARK_CORE_sgemm_tile(Quark *quark, Quark_Task_Flags *task_flags,
                           PLASMA_enum transA, PLASMA_enum transB,
                           int m, int n, int k, int nb,
                           const float *alpha, const float *A, int lda,
                                                            const float *B, int ldb,
                           const float *beta,        float *C, int ldc,
                           const float *Alock,
                           const float *Block,
                           const float *Clock)
{
    plasma_profile_by_kernel( task_flags, GEMM );

    QUARK_Insert_Task(quark, CORE_sgemm_tile_quark, task_flags,
        sizeof(PLASMA_enum),              &transA, VALUE,
        sizeof(PLASMA_enum),              &transB, VALUE,
        sizeof(int),                      &m,      VALUE,
        sizeof(int),                      &n,      VALUE,
        sizeof(int),                      &k,      VALUE,
        sizeof(float),       alpha,           INPUT,
        sizeof(float)*nb*nb, A,               NODEP,          /* input; see Alock */
        sizeof(int),                      &lda,    VALUE,
        sizeof(float)*nb*nb, B,               NODEP,          /* input; see Block */
        sizeof(int),                      &ldb,    VALUE,
        sizeof(float),       beta,            INPUT,
        sizeof(float)*nb*nb, C,                       NODEP,  /* inout; see Clock */
        sizeof(int),                      &ldc,    VALUE,
        sizeof(float)*nb*nb, Alock,           INPUT,
        sizeof(float)*nb,    Block,           INPUT,
        sizeof(float)*nb,    Clock,                   INOUT,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_sgemm_tile_quark = PCORE_sgemm_tile_quark
#define CORE_sgemm_tile_quark PCORE_sgemm_tile_quark
#endif
void CORE_sgemm_tile_quark(Quark *quark)
{
    PLASMA_enum transA, transB;
    int m, n, k, lda, ldb, ldc;
    const float *alpha, *beta;
    const float *A, *B;
    float *C;

    quark_unpack_args_13( quark, transA, transB, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc );
    cblas_sgemm(
        CblasColMajor,
        (CBLAS_TRANSPOSE)transA, (CBLAS_TRANSPOSE)transB,
        m, n, k,
        (*alpha), A, lda,
                             B, ldb,
        (*beta),  C, ldc );
}

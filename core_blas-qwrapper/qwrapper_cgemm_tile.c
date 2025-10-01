/**
 *
 * @file qwrapper_cgemm_tile.c
 *
 *  PLASMA core_blas quark wrapper
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Mark Gates
 * @date 2010-11-15
 * @generated c Fri Apr  1 11:02:43 2016
 *
 **/
#include <cblas.h>
#include "common.h"

/***************************************************************************//**
 *
 * Version of cgemm for tile storage, to avoid dependency problem when
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
void QUARK_CORE_cgemm_tile(Quark *quark, Quark_Task_Flags *task_flags,
                           PLASMA_enum transA, PLASMA_enum transB,
                           int m, int n, int k, int nb,
                           const PLASMA_Complex32_t *alpha, const PLASMA_Complex32_t *A, int lda,
                                                            const PLASMA_Complex32_t *B, int ldb,
                           const PLASMA_Complex32_t *beta,        PLASMA_Complex32_t *C, int ldc,
                           const PLASMA_Complex32_t *Alock,
                           const PLASMA_Complex32_t *Block,
                           const PLASMA_Complex32_t *Clock)
{
    plasma_profile_by_kernel( task_flags, GEMM );

    QUARK_Insert_Task(quark, CORE_cgemm_tile_quark, task_flags,
        sizeof(PLASMA_enum),              &transA, VALUE,
        sizeof(PLASMA_enum),              &transB, VALUE,
        sizeof(int),                      &m,      VALUE,
        sizeof(int),                      &n,      VALUE,
        sizeof(int),                      &k,      VALUE,
        sizeof(PLASMA_Complex32_t),       alpha,           INPUT,
        sizeof(PLASMA_Complex32_t)*nb*nb, A,               NODEP,          /* input; see Alock */
        sizeof(int),                      &lda,    VALUE,
        sizeof(PLASMA_Complex32_t)*nb*nb, B,               NODEP,          /* input; see Block */
        sizeof(int),                      &ldb,    VALUE,
        sizeof(PLASMA_Complex32_t),       beta,            INPUT,
        sizeof(PLASMA_Complex32_t)*nb*nb, C,                       NODEP,  /* inout; see Clock */
        sizeof(int),                      &ldc,    VALUE,
        sizeof(PLASMA_Complex32_t)*nb*nb, Alock,           INPUT,
        sizeof(PLASMA_Complex32_t)*nb,    Block,           INPUT,
        sizeof(PLASMA_Complex32_t)*nb,    Clock,                   INOUT,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_cgemm_tile_quark = PCORE_cgemm_tile_quark
#define CORE_cgemm_tile_quark PCORE_cgemm_tile_quark
#endif
void CORE_cgemm_tile_quark(Quark *quark)
{
    PLASMA_enum transA, transB;
    int m, n, k, lda, ldb, ldc;
    const PLASMA_Complex32_t *alpha, *beta;
    const PLASMA_Complex32_t *A, *B;
    PLASMA_Complex32_t *C;

    quark_unpack_args_13( quark, transA, transB, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc );
    cblas_cgemm(
        CblasColMajor,
        (CBLAS_TRANSPOSE)transA, (CBLAS_TRANSPOSE)transB,
        m, n, k,
        CBLAS_SADDR(*alpha), A, lda,
                             B, ldb,
        CBLAS_SADDR(*beta),  C, ldc );
}

/**
 *
 * @file core_dgeqp3_init.c
 *
 *  PLASMA core_blas kernel
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Mark Gates
 * @date 2010-11-15
 * @generated d Fri Apr  1 11:02:31 2016
 *
 **/
#ifdef PLASMA_WITH_MKL
#include <mkl_lapacke.h>
#else
#include <lapacke.h>
#endif
#include "common.h"

/***************************************************************************//**
 *
 * @ingroup CORE_double
 *
 *  CORE_dgeqp3_init initializes jpvt to [1, ..., n].
 *  Uses 1-based indexing for Fortran compatability.
 *
 *******************************************************************************
 *
 *  @param[in] n
 *          Size of vector.
 *
 *  @param[in,out] jpvt
 *          Vector of size n.
 *          On exit, jpvt[i] = i+1.
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_dgeqp3_init = PCORE_dgeqp3_init
#define CORE_dgeqp3_init PCORE_dgeqp3_init
#endif
void CORE_dgeqp3_init( int n, int *jpvt )
{
    int j;
    for( j = 0; j < n; ++j ) {
        jpvt[j] = j+1;
    }
}

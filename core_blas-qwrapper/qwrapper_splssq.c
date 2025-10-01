/**
 *
 * @file qwrapper_splssq.c
 *
 *  PLASMA core_blas quark wrapper
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @generated s Fri Apr  1 11:02:41 2016
 *
 **/
#include <math.h>
#include "common.h"

/*****************************************************************************
 *
 * @ingroup CORE_float
 *
 *  QUARK_CORE_splssq returns: scl * sqrt(ssq)
 *
 * with scl and ssq such that
 *
 *    ( scl**2 )*ssq = sum( A( 2*i )**2 * A( 2*i+1 ) )
 *                      i
 *
 * The values of A(2*i+1) are assumed to be at least unity.
 * The values of A(2*i) are assumed to be non-negative and scl is
 *
 *    scl = max( A( 2*i ) ),
 *           i
 *
 * The routine makes only one pass through the matrix A.
 *
 *******************************************************************************
 *
 *  @param[in] M
 *          The number of couple (scale, sumsq) in the matrix A.
 *
 *  @param[in] A
 *          The 2-by-M matrix.
 *
 *  @param[out] result
 *          On exit, result contains scl * sqrt( ssq )
 *
 */
void QUARK_CORE_splssq( Quark *quark, Quark_Task_Flags *task_flags,
                        int m, const float *A, float *result )
{
    plasma_profile_by_kernel( task_flags, LASSQ );

    QUARK_Insert_Task(quark, CORE_splssq_quark, task_flags,
        sizeof(int),       &m,      VALUE,
        sizeof(float)*m,   A,          INOUT,
        sizeof(float)*1,   result,     OUTPUT,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_splssq_quark = PCORE_splssq_quark
#define CORE_splssq_quark PCORE_splssq_quark
#endif
void CORE_splssq_quark(Quark *quark)
{
    int m, i;
    float *A;
    float *result;

    quark_unpack_args_3( quark, m, A, result );

    for( i=1; i<m; i++ ) {
        if( A[0] < A[2*i] ) {
            A[1] = A[2*i+1] + (A[1]     * (( A[0] / A[2*i] ) * ( A[0] / A[2*i] )));
            A[0] = A[2*i];
        } else {
            A[1] = A[1]     + (A[2*i+1] * (( A[2*i] / A[0] ) * ( A[2*i] / A[0] )));
        }
    }

    *result = A[0] * sqrt( A[1] );
}

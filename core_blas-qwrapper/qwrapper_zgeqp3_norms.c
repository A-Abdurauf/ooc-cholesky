/**
 *
 * @file qwrapper_zgeqp3_norms.c
 *
 *  PLASMA core_blas quark wrapper
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Mark Gates
 * @date 2010-11-15
 * @precisions normal z -> c d s
 *
 **/
#include "common.h"

#define A(m,n) BLKADDR( A, PLASMA_Complex64_t, m, n )

/***************************************************************************//**
 *
 **/
void QUARK_CORE_zgeqp3_norms( Quark *quark, Quark_Task_Flags *task_flags,
                              PLASMA_desc A,
                              int ioff, int joff,
                              double *norms1, double *norms2 )
{
    Quark_Task *task;
    int ii, jj;

    DAG_SET_PROPERTIES("norms", "brown");
    task = QUARK_Task_Init( quark, CORE_zgeqp3_norms_quark, task_flags );

    QUARK_Task_Pack_Arg( quark, task, sizeof(PLASMA_desc),  &A,      VALUE          );
    QUARK_Task_Pack_Arg( quark, task, sizeof(int),          &ioff,   VALUE          );
    QUARK_Task_Pack_Arg( quark, task, sizeof(int),          &joff,   VALUE          );
    QUARK_Task_Pack_Arg( quark, task, sizeof(double)*A.nb,  norms1,          INOUT  );
    QUARK_Task_Pack_Arg( quark, task, sizeof(double)*A.nb,  norms2,          INOUT  );

    /* depends on block column */
    /* note: A.nt must be 1. checked in CORE_zgeqp3_norms. */
    for( jj = 0; jj < A.nt; ++jj ) {
        for( ii = 0; ii < A.mt; ++ii ) {
            QUARK_Task_Pack_Arg( quark, task, sizeof(PLASMA_Complex64_t)*A.mb*A.nb, A(ii,jj), INPUT );
        }
    }

    QUARK_Insert_Task_Packed( quark, task );
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_zgeqp3_norms_quark = PCORE_zgeqp3_norms_quark
#define CORE_zgeqp3_norms_quark PCORE_zgeqp3_norms_quark
#endif
void CORE_zgeqp3_norms_quark( Quark *quark )
{
    PLASMA_desc A;
    int ioff, joff;
    double *norms1, *norms2;

    quark_unpack_args_5( quark, A, ioff, joff, norms1, norms2 );
    CORE_zgeqp3_norms(          A, ioff, joff, norms1, norms2 );
}

/**
 *
 * @file qwrapper_dgeqp3_larfg.c
 *
 *  PLASMA core_blas quark wrapper
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Mark Gates
 * @date 2010-11-15
 * @generated d Fri Apr  1 11:02:43 2016
 *
 **/
#include "common.h"

#define A(m,n) BLKADDR( A, double, m, n )

/***************************************************************************//**
 *
 **/
void QUARK_CORE_dgeqp3_larfg(Quark *quark, Quark_Task_Flags *task_flags,
                             PLASMA_desc A, int ii, int jj, int i, int j,
                             double *tau, double *beta )
{
    Quark_Task *task;
    int kk;

    DAG_SET_PROPERTIES("larfg", "red");
    task = QUARK_Task_Init( quark, CORE_dgeqp3_larfg_quark, task_flags );

    QUARK_Task_Pack_Arg( quark, task, sizeof(PLASMA_desc),        &A,    VALUE          );
    QUARK_Task_Pack_Arg( quark, task, sizeof(int),                &ii,   VALUE          );
    QUARK_Task_Pack_Arg( quark, task, sizeof(int),                &jj,   VALUE          );
    QUARK_Task_Pack_Arg( quark, task, sizeof(int),                &i,    VALUE          );
    QUARK_Task_Pack_Arg( quark, task, sizeof(int),                &j,    VALUE          );
    QUARK_Task_Pack_Arg( quark, task, sizeof(double), tau,           OUTPUT );
    QUARK_Task_Pack_Arg( quark, task, sizeof(double), beta,          OUTPUT );

    /* depends on block column */
    for( kk = ii; kk < A.mt; ++kk ) {
        QUARK_Task_Pack_Arg( quark, task, sizeof(double)*A.nb*A.nb, A(kk,jj), INOUT );
    }

    QUARK_Insert_Task_Packed( quark, task );
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_dgeqp3_larfg_quark = PCORE_dgeqp3_larfg_quark
#define CORE_dgeqp3_larfg_quark PCORE_dgeqp3_larfg_quark
#endif
void CORE_dgeqp3_larfg_quark(Quark *quark)
{
    PLASMA_desc A;
    int ii, jj, i, j;
    double *tau, *beta;

    quark_unpack_args_7( quark, A, ii, jj, i, j, tau, beta );
    CORE_dgeqp3_larfg(          A, ii, jj, i, j, tau, beta );
}

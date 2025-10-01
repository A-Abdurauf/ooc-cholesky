/**
 *
 * @file qwrapper_zgeqp3_pivot.c
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
 * Note: this executes synchronously.
 **/

void QUARK_CORE_zgeqp3_pivot( Quark *quark, Quark_Task_Flags *task_flags,
                              PLASMA_desc A,
                              PLASMA_Complex64_t *F, int ldf,
                              int jj, int k, int *jpvt,
                              double *norms1, double *norms2, int *info )
{
    Quark_Task *task;
    int ii, jj2;

    DAG_SET_PROPERTIES("pivot", "lightblue");
    task = QUARK_Task_Init( quark, CORE_zgeqp3_pivot_quark, task_flags );

    QUARK_Task_Pack_Arg( quark, task, sizeof(PLASMA_desc),                  &A,      VALUE          );
    QUARK_Task_Pack_Arg( quark, task, sizeof(PLASMA_Complex64_t)*A.nb*A.nb, F,               INOUT  );
    QUARK_Task_Pack_Arg( quark, task, sizeof(int),                          &ldf,    VALUE          );
    QUARK_Task_Pack_Arg( quark, task, sizeof(int),                          &jj,     VALUE          );
    QUARK_Task_Pack_Arg( quark, task, sizeof(int),                          &k,      VALUE          );
    QUARK_Task_Pack_Arg( quark, task, sizeof(double)*A.n,                   jpvt,            INOUT  );
    QUARK_Task_Pack_Arg( quark, task, sizeof(double)*A.nb,                  norms1,          INOUT  );
    QUARK_Task_Pack_Arg( quark, task, sizeof(double)*A.nb,                  norms2,          NODEP  );  /* INOUT, but implied by norms1 */
    QUARK_Task_Pack_Arg( quark, task, sizeof(int),                          info,            NODEP  );  /* IN,    but implied by norms1 */

    /* depends on whole matrix to right: A[:,jj:] */
    /* TODO is this implicit in norms1 dependency? */
    for( jj2 = jj; jj2 < A.nt; ++jj2 ) {
        for( ii = 0; ii < A.mt; ++ii ) {
            QUARK_Task_Pack_Arg( quark, task, sizeof(PLASMA_Complex64_t)*A.nb*A.nb, A(ii,jj), INOUT );
        }
    }

    /* depends on whole F matrix */
    /* TODO is this implicit in norms1 dependency? */
    /* TODO use Fk = F[(kk-jj)*A.nb] for kk = jj to A.nt ? */
    for( ii = 1; ii < A.nt; ++ii ) {
        QUARK_Task_Pack_Arg( quark, task, sizeof(PLASMA_Complex64_t)*A.nb*A.nb, &F[ii*A.nb], INOUT );
    }

    /* depends on whole norms1 vector (and implicitly, norms2 vector) */
    for( ii = 1; ii < A.nt; ++ii ) {
        QUARK_Task_Pack_Arg( quark, task, sizeof(PLASMA_Complex64_t)*A.nb, &norms1[ii*A.nb], INOUT );
    }

    /* run synchronously -- as soon as dependencies have been fulfilled. */
    QUARK_Execute_Task_Packed( quark, task );
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_zgeqp3_pivot_quark = PCORE_zgeqp3_pivot_quark
#define CORE_zgeqp3_pivot_quark PCORE_zgeqp3_pivot_quark
#endif
void CORE_zgeqp3_pivot_quark( Quark *quark )
{
    PLASMA_desc A;
    PLASMA_Complex64_t *F;
    int ldf, jj, k;
    int *jpvt;
    double *norms1, *norms2;
    int *info;

    quark_unpack_args_9( quark, A, F, ldf, jj, k, jpvt, norms1, norms2, info );
    CORE_zgeqp3_pivot(          A, F, ldf, jj, k, jpvt, norms1, norms2, info );
}

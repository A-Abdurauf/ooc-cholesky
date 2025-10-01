/**
 *
 * @file qwrapper_dlaed2_computeK.c
 *
 *  PLASMA core_blas quark wrapper
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Gregoire Pichon
 * @author Azzam Haidar
 * @date 2014-07
 * @precisions normal d -> s
 *
 **/
#include <stdlib.h>
#include "common.h"

void
CORE_dlaed2_computeK_quark(Quark *quark);

/***************************************************************************//**
 *
 **/
void QUARK_CORE_dlaed2_computeK(Quark *quark, Quark_Task_Flags *task_flags,
                                int *K1, int n, int n1,
                                double *beta, double *D, double *Q, int LDQ,
                                double *Z, double *DLAMBDA, double *W,
                                int *INDX, int *INDXC, int *INDXP, int *INDXQ,
                                int *COLTYP,
                                double **Qmerge, int wsmode,
                                int *K2)
{
    plasma_profile_by_kernel( task_flags, LAED2_COMPUTEK );

    QUARK_Insert_Task(quark, CORE_dlaed2_computeK_quark, task_flags,
                      sizeof(int),           K1,          INOUT,
                      sizeof(int),          &n,         VALUE,
                      sizeof(int),          &n1,        VALUE,
                      sizeof(double),        beta,        NODEP, /* VALUE by ptr */
                      sizeof(double)*n,      D,           NODEP, /* INOUT  */
                      sizeof(double)*LDQ*n,  Q,           NODEP, /* INOUT  */
                      sizeof(int),          &LDQ,       VALUE,
                      sizeof(double)*n,      Z,           NODEP, /* OUTPUT */
                      sizeof(double)*n,      DLAMBDA,     NODEP, /* OUTPUT */
                      sizeof(double)*n,      W,           NODEP, /* OUTPUT */
                      sizeof(int)*n,         INDX,        NODEP, /* OUTPUT */
                      sizeof(int)*n,         INDXC,       NODEP, /* OUTPUT */
                      sizeof(int)*n,         INDXP,       NODEP, /* OUTPUT */
                      sizeof(int)*n,         INDXQ,       NODEP, /* INOUT  */
                      sizeof(int)*n,         COLTYP,      NODEP, /* OUTPUT */
                      sizeof(double*),       Qmerge,      NODEP, /* OUTPUT */
                      sizeof(int),          &wsmode, VALUE,
                      sizeof(int),           K2,          INOUT,
                      0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_dlaed2_computeK_quark = PCORE_dlaed2_computeK_quark
#define CORE_dlaed2_computeK_quark PCORE_dlaed2_computeK_quark
#endif
void CORE_dlaed2_computeK_quark(Quark *quark)
{
    int n;
    int n1;
    int *K1, *K2;
    double *D;
    double *beta;
    double *Q;
    int LDQ;
    double *Z;
    double *DLAMBDA;
    double *W;
    int *INDX;
    int *INDXC;
    int *INDXP;
    int *INDXQ;
    int *COLTYP;
    double **Qmerge;
    int wsmode;

    quark_unpack_args_18(quark, K1, n, n1,
               beta, D, Q, LDQ,
               Z, DLAMBDA, W,
               INDX, INDXC, INDXP, INDXQ,
               COLTYP,
               Qmerge, wsmode, K2);


    CORE_dlaed2_computeK(K1, n, n1,
               beta, D, Q, LDQ,
               Z, DLAMBDA, W,
               INDX, INDXC, INDXP, INDXQ,
               COLTYP);

    /* If workspace mode is equal to 3, we prepare Qmerge space for the merge step */
    *Qmerge = NULL;
    if( wsmode == 3 ) {
        size_t size = (*K1) * (*K1) * sizeof(double);
        if (size > 0) {
            *Qmerge = malloc( size );
        }
    }
}

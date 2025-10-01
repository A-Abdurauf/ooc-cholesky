/**
 *
 * @file qwrapper_slaed2_computeK.c
 *
 *  PLASMA core_blas quark wrapper
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Gregoire Pichon
 * @author Azzam Haidar
 * @date 2014-07
 * @generated s Fri Apr  1 11:02:43 2016
 *
 **/
#include <stdlib.h>
#include "common.h"

void
CORE_slaed2_computeK_quark(Quark *quark);

/***************************************************************************//**
 *
 **/
void QUARK_CORE_slaed2_computeK(Quark *quark, Quark_Task_Flags *task_flags,
                                int *K1, int n, int n1,
                                float *beta, float *D, float *Q, int LDQ,
                                float *Z, float *DLAMBDA, float *W,
                                int *INDX, int *INDXC, int *INDXP, int *INDXQ,
                                int *COLTYP,
                                float **Qmerge, int wsmode,
                                int *K2)
{
    plasma_profile_by_kernel( task_flags, LAED2_COMPUTEK );

    QUARK_Insert_Task(quark, CORE_slaed2_computeK_quark, task_flags,
                      sizeof(int),           K1,          INOUT,
                      sizeof(int),          &n,         VALUE,
                      sizeof(int),          &n1,        VALUE,
                      sizeof(float),        beta,        NODEP, /* VALUE by ptr */
                      sizeof(float)*n,      D,           NODEP, /* INOUT  */
                      sizeof(float)*LDQ*n,  Q,           NODEP, /* INOUT  */
                      sizeof(int),          &LDQ,       VALUE,
                      sizeof(float)*n,      Z,           NODEP, /* OUTPUT */
                      sizeof(float)*n,      DLAMBDA,     NODEP, /* OUTPUT */
                      sizeof(float)*n,      W,           NODEP, /* OUTPUT */
                      sizeof(int)*n,         INDX,        NODEP, /* OUTPUT */
                      sizeof(int)*n,         INDXC,       NODEP, /* OUTPUT */
                      sizeof(int)*n,         INDXP,       NODEP, /* OUTPUT */
                      sizeof(int)*n,         INDXQ,       NODEP, /* INOUT  */
                      sizeof(int)*n,         COLTYP,      NODEP, /* OUTPUT */
                      sizeof(float*),       Qmerge,      NODEP, /* OUTPUT */
                      sizeof(int),          &wsmode, VALUE,
                      sizeof(int),           K2,          INOUT,
                      0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_slaed2_computeK_quark = PCORE_slaed2_computeK_quark
#define CORE_slaed2_computeK_quark PCORE_slaed2_computeK_quark
#endif
void CORE_slaed2_computeK_quark(Quark *quark)
{
    int n;
    int n1;
    int *K1, *K2;
    float *D;
    float *beta;
    float *Q;
    int LDQ;
    float *Z;
    float *DLAMBDA;
    float *W;
    int *INDX;
    int *INDXC;
    int *INDXP;
    int *INDXQ;
    int *COLTYP;
    float **Qmerge;
    int wsmode;

    quark_unpack_args_18(quark, K1, n, n1,
               beta, D, Q, LDQ,
               Z, DLAMBDA, W,
               INDX, INDXC, INDXP, INDXQ,
               COLTYP,
               Qmerge, wsmode, K2);


    CORE_slaed2_computeK(K1, n, n1,
               beta, D, Q, LDQ,
               Z, DLAMBDA, W,
               INDX, INDXC, INDXP, INDXQ,
               COLTYP);

    /* If workspace mode is equal to 3, we prepare Qmerge space for the merge step */
    *Qmerge = NULL;
    if( wsmode == 3 ) {
        size_t size = (*K1) * (*K1) * sizeof(float);
        if (size > 0) {
            *Qmerge = malloc( size );
        }
    }
}

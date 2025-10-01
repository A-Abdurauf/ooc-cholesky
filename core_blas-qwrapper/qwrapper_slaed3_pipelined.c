/**
 *
 * @file qwrapper_slaed3_pipelined.c
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
#include "common.h"

void
CORE_slaed3_pipelined_quark(Quark *quark)
{
    int n, n1;
    const int *K;
    const int *il_nondef;
    const int *iu_nondef;
    float *D, *Q;
    int LDQ;
    float *Q2;
    float *W, *S;
    const int *INDXC;
    int *INDXQ, *COLTYP;
    int start;
    int end;
    int il, iu;

    quark_unpack_args_16(
        quark, n, n1, K, il_nondef, iu_nondef,
        D, Q, LDQ, Q2,
        INDXC, INDXQ, COLTYP, W, S,
        start, end);

    il = ( il_nondef == NULL ) ? 0 : *il_nondef;
    iu = ( iu_nondef == NULL ) ? n : *iu_nondef;

    CORE_slaed3_computevectors(*K, il, iu,
                               Q, LDQ, W, S, INDXC,
                               start, end);

    if (start == 0)
    {
        CORE_slaed3_merge(n, *K, D, INDXQ);
    }

    CORE_slaed3_updatevectors(
        PlasmaLaed3UpdateAll, 0, n, n1, *K, il, iu,
        Q, LDQ, Q2, COLTYP, NULL, start, end);
}

/***************************************************************************//**
                                                                              *
                                                                              **/
void QUARK_CORE_slaed3_pipelined(Quark *quark, Quark_Task_Flags *task_flags,
                                 int n, int n1, int *K, int *il_nondef, int *iu_nondef,
                                 float *D, float *Q, int LDQ, float *Q2,
                                 int *INDXC, int *INDXQ, int *COLTYP, float *W,
                                 int start, int end)
{
    plasma_gendag_by_kernel( task_flags, LAED3_PIPELINED );

    QUARK_Insert_Task(quark, CORE_slaed3_pipelined_quark, task_flags,
        sizeof(int),        &n,         VALUE,
        sizeof(int),        &n1,        VALUE,
        sizeof(int),         K,             INPUT,
        sizeof(int),         il_nondef,     NODEP,
        sizeof(int),         iu_nondef,     NODEP,
        sizeof(float),      D,             NODEP,
        sizeof(float),      Q,             NODEP,
        sizeof(int),        &LDQ,       VALUE,
        sizeof(float),      Q2,            NODEP,
        sizeof(int),         INDXC,         NODEP,
        sizeof(int),         INDXQ,         NODEP,
        sizeof(int),         COLTYP,        NODEP,
        sizeof(float),      W,             INPUT,
        sizeof(float)*LDQ,  NULL,          SCRATCH,//replace LDQ by n.should be n
        sizeof(int),        &start,     VALUE,
        sizeof(int),        &end,       VALUE,
        0);
}

/**
 *
 * @file qwrapper_dlaed3_updatevectors.c
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
#include "common.h"

void
CORE_dlaed3_updatevectors_quark(Quark *quark);

/***************************************************************************//**
 *
 **/
void QUARK_CORE_dlaed3_updatevectors(Quark *quark, Quark_Task_Flags *task_flags,
                                     int oper, int wsmode, int n, int n1, int *K,
                                     int *il_nondef, int *iu_nondef,
                                     double *D, double *Q, int LDQ, double *Q2,
                                     int *INDXQ, int *COLTYP, double **WORK,
                                     int start, int end, double **WORKDEP)
{
    plasma_profile_by_kernel( task_flags, LAED3_UPDATEVECTORS );

    QUARK_Insert_Task(quark, CORE_dlaed3_updatevectors_quark, task_flags,
        sizeof(int),    &oper,        VALUE,
        sizeof(int),    &wsmode,      VALUE,
        sizeof(int),    &n,           VALUE,
        sizeof(int),    &n1,          VALUE,
        sizeof(int),     K,               INPUT,
        sizeof(int),     il_nondef,       NODEP,
        sizeof(int),     iu_nondef,       NODEP,
        sizeof(double),  D,               NODEP,
        sizeof(double),  Q,               NODEP,
        sizeof(int),    &LDQ,         VALUE,
        sizeof(double),  Q2,              NODEP,
        sizeof(int),     INDXQ,           NODEP,
        sizeof(int),     COLTYP,          NODEP,
        sizeof(double),  WORK,            NODEP,
        sizeof(int),    &start,       VALUE,
        sizeof(int),    &end,         VALUE,
        sizeof(double),  Q+start*LDQ, ( wsmode > 0 ? NODEP : INPUT ),
        sizeof(double),  WORKDEP,     ( wsmode > 0 ? INPUT : NODEP ),
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_dlaed3_updatevectors_quark = PCORE_dlaed3_updatevectors_quark
#define CORE_dlaed3_updatevectors_quark PCORE_dlaed3_updatevectors_quark
#endif
void CORE_dlaed3_updatevectors_quark(Quark *quark)
{
    int oper, wsmode, n, n1;
    const int *K;
    int *il_nondef, *iu_nondef;
    int il, iu;
    double *D;
    double *Q;
    int LDQ;
    double *Q2;
    int *INDXQ, *COLTYP;
    double **WORK;
    int start;
    int end;
    void *fake1, *fake2;

    quark_unpack_args_18(
        quark, oper, wsmode, n, n1, K, il_nondef, iu_nondef,
        D, Q, LDQ, Q2,
        INDXQ, COLTYP, WORK,
        start, end, fake1, fake2 );

    if ((start == 0) &&
        (oper & PlasmaLaed3Update2))
    {
        CORE_dlaed3_merge(n, *K, D, INDXQ);
    }

    il = ( il_nondef == NULL ) ? 0 : *il_nondef;
    iu = ( iu_nondef == NULL ) ? n : *iu_nondef;
    CORE_dlaed3_updatevectors(oper, wsmode, n, n1, *K, il, iu,
                              Q, LDQ, Q2, COLTYP, *WORK, start, end);
}

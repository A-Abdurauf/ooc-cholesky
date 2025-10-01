/**
 *
 * @file cgecfi2.c
 *
 *  PLASMA InPlaceTransformation module
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 *  This work is the implementation of an inplace transformation
 *  based on the GKK algorithm by Gustavson, Karlsson, Kagstrom
 *  and its fortran implementation.
 *
 * @version 2.8.0
 * @author Mathieu Faverge
 * @date 2010-11-15
 *
 * @generated c Fri Apr  1 11:03:00 2016
 *
 **/

#include <sys/types.h>
#include "common.h"
#include "cgecfi2.h"

#define PLASMA_pcgetmi2(idep, odep, storev, m, n, mb, nb, A)     \
    plasma_parallel_call_10(                                     \
        plasma_pcgetmi2,                                         \
        PLASMA_enum,  (idep),                                    \
        PLASMA_enum,  (odep),                                    \
        PLASMA_enum,  (storev),                                  \
        int,  (m),                                               \
        int,  (n),                                               \
        int,  (mb),                                              \
        int,  (nb),                                              \
        PLASMA_Complex32_t*, (A),                                \
        PLASMA_sequence*, sequence,                              \
        PLASMA_request*, request);

#define PLASMA_cshift(m, n, mb, nb, A)                          \
    plasma_cshift(plasma, (m), (n), (A),                        \
                  ( (n) / (nb) ), ( (m) / (mb) ), (nb), (mb),   \
                  sequence, request);

#define PLASMA_cshiftr(m, n, mb, nb, A)                         \
    plasma_cshift(plasma, (m), (n), (A),                        \
                  ( (n) / (nb) ), (nb), ( (m) / (mb) ), (mb),   \
                  sequence, request);

/** ****************************************************************************
 *
 * @ingroup InPlaceTransformation
 *
 *  ipt_ccm2ccrb converts a matrix from CM format to CCRB format
 *
 *******************************************************************************
 *
 * @param[in] plasma
 *          Plasma context to which this call belong to.
 *
 * @param[in] m
 *         Number of rows of matrix A
 *
 * @param[in] n
 *         Number of columns of matrix A
 *
 * @param[in,out] A
 *         Matrix of size m*n.
 *
 * @param[in] mb
 *         Number of rows of each block
 *
 * @param[in] nb
 *         Number of columns of each block
 *
 * @param[in] sequence
 *          Identifies the sequence of function calls that this call belongs to
 *          (for completion checks and exception handling purposes).
 *
 * @param[out] request
 *          Identifies this function call (for exception handling purposes).
 *
 ******************************************************************************/


/** ****************************************************************************
 *
 *  Self-contained functions
 *
 *******************************************************************************/
/*
 * Shift inside panels
 */
int ipt_ccm2ccrb(plasma_context_t *plasma, int m, int n, PLASMA_Complex32_t *A, int mb, int nb,
                 PLASMA_sequence *sequence, PLASMA_request *request)
{
    if( (m == 0) || (n == 0) )
        return PLASMA_SUCCESS;

    PLASMA_cshift(m, n, mb, nb, A);
    ipt_cpanel2tile(plasma, m, n, A, mb, nb, sequence, request);

    return PLASMA_SUCCESS;
}

int ipt_cccrb2cm(plasma_context_t *plasma, int m, int n, PLASMA_Complex32_t *A, int mb, int nb,
                 PLASMA_sequence *sequence, PLASMA_request *request)
{
    if( (m == 0) || (n == 0) )
        return PLASMA_SUCCESS;

    ipt_ctile2panel(plasma, m, n, A, mb, nb, sequence, request);
    PLASMA_cshiftr(m, n, mb, nb, A);

    return PLASMA_SUCCESS;
}

/*
 * Transpose each tile
 */
int ipt_cccrb2crrb(plasma_context_t *plasma, PLASMA_enum idep, PLASMA_enum odep, int m, int n, PLASMA_Complex32_t *A, int mb, int nb,
                   PLASMA_sequence *sequence, PLASMA_request *request)
{
    if( (m == 0) || (n == 0) )
        return PLASMA_SUCCESS;

    PLASMA_pcgetmi2(idep, odep, PlasmaColumnwise, m, n, mb, nb, A);

    return PLASMA_SUCCESS;
}

int ipt_ccrrb2ccrb(plasma_context_t *plasma, PLASMA_enum idep, PLASMA_enum odep, int m, int n, PLASMA_Complex32_t *A, int mb, int nb,
                   PLASMA_sequence *sequence, PLASMA_request *request)
{
    if( (m == 0) || (n == 0) )
        return PLASMA_SUCCESS;

    PLASMA_pcgetmi2(idep, odep, PlasmaRowwise, n, m, nb, mb, A);

    return PLASMA_SUCCESS;
}

int ipt_crcrb2rrrb(plasma_context_t *plasma, PLASMA_enum idep, PLASMA_enum odep, int m, int n, PLASMA_Complex32_t *A, int mb, int nb,
                   PLASMA_sequence *sequence, PLASMA_request *request)
{
    if( (m == 0) || (n == 0) )
        return PLASMA_SUCCESS;

    PLASMA_pcgetmi2(idep, odep, PlasmaRowwise, m, n, mb, nb, A);

    return PLASMA_SUCCESS;
}

int ipt_crrrb2rcrb(plasma_context_t *plasma, PLASMA_enum idep, PLASMA_enum odep, int m, int n, PLASMA_Complex32_t *A, int mb, int nb,
                   PLASMA_sequence *sequence, PLASMA_request *request)
{
    if( (m == 0) || (n == 0) )
        return PLASMA_SUCCESS;

    PLASMA_pcgetmi2(idep, odep, PlasmaColumnwise, n, m, nb, mb, A);

    return PLASMA_SUCCESS;
}

/*
 * Transpose all tiles
 */
int ipt_cccrb2rcrb(plasma_context_t *plasma, int m, int n, PLASMA_Complex32_t *A, int mb, int nb,
                   PLASMA_sequence *sequence, PLASMA_request *request)
{
    int M_, N_;

    if( (m == 0) || (n == 0) )
        return PLASMA_SUCCESS;

    M_ = m / mb;
    N_ = n / nb;

    /* quick return */
    if( (M_ < 2) || (N_ < 2) ) {
        return PLASMA_SUCCESS;
    }

    plasma_cshift(plasma, m, n, A, 1, ( m / mb ), ( n / nb ), (mb*nb),
                  sequence, request);

    return PLASMA_SUCCESS;
}

/** ****************************************************************************
 *
 *  Composition of 2 sub-routines
 *
 *******************************************************************************/
/*
 * Shift inside panels + Transpose all tiles
 */
int ipt_ccm2rcrb(plasma_context_t *plasma, int m, int n, PLASMA_Complex32_t *A, int mb, int nb,
                 PLASMA_sequence *sequence, PLASMA_request *request)
{
    if( (m == 0) || (n == 0) )
        return PLASMA_SUCCESS;
    //ipt_ccm2ccrb(  plasma, m, n, A, mb, nb, sequence, request);
    PLASMA_cshift(m, n, mb, nb, A);
    ipt_cpanel2all(plasma, m, n, A, mb, nb, sequence, request);
    ipt_cccrb2rcrb(plasma, m, n, A, mb, nb, sequence, request);

    return PLASMA_SUCCESS;
}

int ipt_crcrb2cm(plasma_context_t *plasma, int m, int n, PLASMA_Complex32_t *A, int mb, int nb,
                 PLASMA_sequence *sequence, PLASMA_request *request)
{
    if( (m == 0) || (n == 0) )
        return PLASMA_SUCCESS;
    ipt_crcrb2ccrb(plasma, m, n, A, mb, nb, sequence, request);
    ipt_call2panel(plasma, m, n, A, mb, nb, sequence, request);
    //ipt_cccrb2cm(  plasma, m, n, A, mb, nb, sequence, request);
    PLASMA_cshiftr(m, n, mb, nb, A);

    return PLASMA_SUCCESS;
}

/*
 * Transpose each tile  + Transpose all tiles
 */
int ipt_cccrb2rrrb(plasma_context_t *plasma, int m, int n, PLASMA_Complex32_t *A, int mb, int nb,
                   PLASMA_sequence *sequence, PLASMA_request *request)
{
    if( (m == 0) || (n == 0) )
        return PLASMA_SUCCESS;
    ipt_cccrb2rcrb(plasma, m, n, A, mb, nb, sequence, request);
    ipt_crcrb2rrrb(plasma, PlasmaIPT_All, PlasmaIPT_NoDep, m, n, A, mb, nb, sequence, request);

    return PLASMA_SUCCESS;
}

int ipt_crrrb2ccrb(plasma_context_t *plasma, int m, int n, PLASMA_Complex32_t *A, int mb, int nb,
                 PLASMA_sequence *sequence, PLASMA_request *request)
{
    if( (m == 0) || (n == 0) )
        return PLASMA_SUCCESS;
    ipt_crrrb2rcrb(plasma, PlasmaIPT_NoDep, PlasmaIPT_All, m, n, A, mb, nb, sequence, request);
    ipt_crcrb2ccrb(plasma, m, n, A, mb, nb, sequence, request);

    return PLASMA_SUCCESS;
}

int ipt_crcrb2crrb(plasma_context_t *plasma, int m, int n, PLASMA_Complex32_t *A, int mb, int nb,
                 PLASMA_sequence *sequence, PLASMA_request *request)
{
    if( (m == 0) || (n == 0) )
        return PLASMA_SUCCESS;
    ipt_crcrb2ccrb(plasma, m, n, A, mb, nb, sequence, request);
    ipt_cccrb2crrb(plasma, PlasmaIPT_All, PlasmaIPT_NoDep, m, n, A, mb, nb, sequence, request);

    return PLASMA_SUCCESS;
}

int ipt_ccrrb2rcrb(plasma_context_t *plasma, int m, int n, PLASMA_Complex32_t *A, int mb, int nb,
                   PLASMA_sequence *sequence, PLASMA_request *request)
{
    if( (m == 0) || (n == 0) )
        return PLASMA_SUCCESS;
    ipt_ccrrb2ccrb(plasma, PlasmaIPT_NoDep, PlasmaIPT_All, m, n, A, mb, nb, sequence, request);
    ipt_cccrb2rcrb(plasma, m, n, A, mb, nb, sequence, request);

    return PLASMA_SUCCESS;
}

/*
 * Transpose each tile  + Shift inside panels
 */
int ipt_ccm2crrb(plasma_context_t *plasma, int m, int n, PLASMA_Complex32_t *A, int mb, int nb,
                 PLASMA_sequence *sequence, PLASMA_request *request)
{
    if( (m == 0) || (n == 0) )
        return PLASMA_SUCCESS;
    //ipt_ccm2ccrb(  plasma, m, n, A, mb, nb, sequence, request);
    PLASMA_cshift(m, n, mb, nb, A);
    ipt_cccrb2crrb(plasma, PlasmaIPT_Panel, PlasmaIPT_NoDep, m, n, A, mb, nb, sequence, request);

    return PLASMA_SUCCESS;
}

int ipt_ccrrb2cm(plasma_context_t *plasma, int m, int n, PLASMA_Complex32_t *A, int mb, int nb,
                 PLASMA_sequence *sequence, PLASMA_request *request)
{
    if( (m == 0) || (n == 0) )
        return PLASMA_SUCCESS;
    ipt_ccrrb2ccrb(plasma, PlasmaIPT_NoDep, PlasmaIPT_Panel, m, n, A, mb, nb, sequence, request);
    //ipt_cccrb2cm(  plasma, m, n, A, mb, nb, sequence, request);
    PLASMA_cshiftr(m, n, mb, nb, A);

    return PLASMA_SUCCESS;
}

int ipt_crm2rcrb(plasma_context_t *plasma, int m, int n, PLASMA_Complex32_t *A, int mb, int nb,
                 PLASMA_sequence *sequence, PLASMA_request *request)
{
    if( (m == 0) || (n == 0) )
        return PLASMA_SUCCESS;
    ipt_crm2rrrb(  plasma, m, n, A, mb, nb, sequence, request);
    ipt_crrrb2rcrb(plasma, PlasmaIPT_Panel, PlasmaIPT_NoDep, m, n, A, mb, nb, sequence, request);

    return PLASMA_SUCCESS;
}

int ipt_crcrb2rm(plasma_context_t *plasma, int m, int n, PLASMA_Complex32_t *A, int mb, int nb,
                 PLASMA_sequence *sequence, PLASMA_request *request)
{
    if( (m == 0) || (n == 0) )
        return PLASMA_SUCCESS;
    ipt_crcrb2rrrb(plasma, PlasmaIPT_NoDep, PlasmaIPT_Panel, m, n, A, mb, nb, sequence, request);
    ipt_crrrb2rm(  plasma, m, n, A, mb, nb, sequence, request);

    return PLASMA_SUCCESS;
}

/** ****************************************************************************
 *
 *  Composition of 3 sub-routines
 *
 *******************************************************************************/
/*
 * Shift inside panels + Transpose all tiles + Transpose inside each tile
 */
int ipt_ccm2rrrb(plasma_context_t *plasma, int m, int n, PLASMA_Complex32_t *A, int mb, int nb,
                 PLASMA_sequence *sequence, PLASMA_request *request)
{
    if( (m == 0) || (n == 0) )
        return PLASMA_SUCCESS;
    //ipt_ccm2ccrb(  plasma, m, n, A, mb, nb, sequence, request);
    PLASMA_cshift(m, n, mb, nb, A);
    ipt_cccrb2crrb(plasma, PlasmaIPT_Panel, PlasmaIPT_All, m, n, A, mb, nb, sequence, request);
    ipt_ccrrb2rrrb(plasma, m, n, A, mb, nb, sequence, request);

    return PLASMA_SUCCESS;
}

int ipt_crrrb2cm(plasma_context_t *plasma, int m, int n, PLASMA_Complex32_t *A, int mb, int nb,
                 PLASMA_sequence *sequence, PLASMA_request *request)
{
    if( (m == 0) || (n == 0) )
        return PLASMA_SUCCESS;
    ipt_crrrb2crrb(plasma, m, n, A, mb, nb, sequence, request);
    ipt_ccrrb2ccrb(plasma, PlasmaIPT_All, PlasmaIPT_Panel, m, n, A, mb, nb, sequence, request);
    //ipt_cccrb2cm(  plasma, m, n, A, mb, nb, sequence, request);
    PLASMA_cshiftr(m, n, mb, nb, A);

    return PLASMA_SUCCESS;
}

int ipt_cccrb2rm(plasma_context_t *plasma, int m, int n, PLASMA_Complex32_t *A, int mb, int nb,
                 PLASMA_sequence *sequence, PLASMA_request *request)
{
    if( (m == 0) || (n == 0) )
        return PLASMA_SUCCESS;
    ipt_cccrb2rcrb(plasma, m, n, A, mb, nb, sequence, request);
    ipt_crcrb2rrrb(plasma, PlasmaIPT_All, PlasmaIPT_Panel, m, n, A, mb, nb, sequence, request);
    ipt_crrrb2rm(  plasma, m, n, A, mb, nb, sequence, request);

    return PLASMA_SUCCESS;
}

int ipt_crm2ccrb(plasma_context_t *plasma, int m, int n, PLASMA_Complex32_t *A, int mb, int nb,
                 PLASMA_sequence *sequence, PLASMA_request *request)
{
    if( (m == 0) || (n == 0) )
        return PLASMA_SUCCESS;
    ipt_crm2rrrb(  plasma, m, n, A, mb, nb, sequence, request);
    ipt_crrrb2rcrb(plasma, PlasmaIPT_Panel, PlasmaIPT_All, m, n, A, mb, nb, sequence, request);
    ipt_crcrb2ccrb(plasma, m, n, A, mb, nb, sequence, request);

    return PLASMA_SUCCESS;
}

/** ****************************************************************************
 *
 *  Composition of 4 sub-routines
 *
 *******************************************************************************/
/*
 * Shift inside panels + Transpose all tiles
 *    + Transpose inside each tile + Shift inside panels
 */
int ipt_ccm2rm(plasma_context_t *plasma, int m, int n, PLASMA_Complex32_t *A, int mb, int nb,
               PLASMA_sequence *sequence, PLASMA_request *request)
{
    if( (m == 0) || (n == 0) )
        return PLASMA_SUCCESS;
    //ipt_ccm2ccrb(  plasma, m, n, A, mb, nb, sequence, request);
    PLASMA_cshift(m, n, mb, nb, A);
    ipt_cpanel2all(plasma, m, n, A, mb, nb, sequence, request);
    ipt_cccrb2rcrb(plasma, m, n, A, mb, nb, sequence, request);
    ipt_crcrb2rrrb(plasma, PlasmaIPT_All, PlasmaIPT_Panel, m, n, A, mb, nb, sequence, request);
    ipt_crrrb2rm(  plasma, m, n, A, mb, nb, sequence, request);

    return PLASMA_SUCCESS;
}

int ipt_crm2cm(plasma_context_t *plasma, int m, int n, PLASMA_Complex32_t *A, int mb, int nb,
                 PLASMA_sequence *sequence, PLASMA_request *request)
{
    if( (m == 0) || (n == 0) )
        return PLASMA_SUCCESS;
    ipt_crm2rrrb(  plasma, m, n, A, mb, nb, sequence, request);
    ipt_crrrb2rcrb(plasma, PlasmaIPT_Panel, PlasmaIPT_All, m, n, A, mb, nb, sequence, request);
    ipt_crcrb2ccrb(plasma, m, n, A, mb, nb, sequence, request);
    ipt_call2panel(plasma, m, n, A, mb, nb, sequence, request);
    //ipt_cccrb2cm(  plasma, m, n, A, mb, nb, sequence, request);
    PLASMA_cshiftr(m, n, mb, nb, A);

    return PLASMA_SUCCESS;
}


/** ****************************************************************************
 *
 *  Barriers
 *
 *******************************************************************************/

int ipt_ctile2panel( plasma_context_t *plasma, int m, int n, PLASMA_Complex32_t *A, int mb, int nb,
                     PLASMA_sequence *sequence, PLASMA_request *request)
{
    if( PLASMA_SCHEDULING != PLASMA_DYNAMIC_SCHEDULING )
        return PLASMA_SUCCESS;

    PLASMA_Complex32_t *Al;
    int i,j;
    int M_ = m / mb;
    int N_ = n / nb;
    int bsiz = mb*nb;
    int psiz = m*nb;
    Quark_Task_Flags task_flags = Quark_Task_Flags_Initializer;
    QUARK_Task_Flag_Set(&task_flags, TASK_SEQUENCE, (intptr_t)sequence->quark_sequence);

    plasma_dynamic_spawn();
    for(j=0; j<N_; j++) {
        Al = &(A[psiz*j]);

        for(i=1; i<M_; i++) {

#ifdef TRACE_IPT
            char str[30];
            sprintf(str, "Foo2 C2RI %d", i*m*nb);
#endif
            QUARK_Insert_Task(plasma->quark, CORE_foo2_quark, &task_flags,
                              sizeof(PLASMA_Complex32_t)*psiz,  Al,           INOUT | GATHERV,
                              sizeof(PLASMA_Complex32_t)*bsiz, &(Al[i*bsiz]), INOUT,
#ifdef TRACE_IPT
                              30, str,   VALUE | TASKLABEL,
                              4, "red",  VALUE | TASKCOLOR,
#endif
                              0);
        }
    }

    return PLASMA_SUCCESS;
}

int ipt_cpanel2tile( plasma_context_t *plasma, int m, int n, PLASMA_Complex32_t *A, int mb, int nb,
                     PLASMA_sequence *sequence, PLASMA_request *request )
{
    if( PLASMA_SCHEDULING != PLASMA_DYNAMIC_SCHEDULING )
        return PLASMA_SUCCESS;

    PLASMA_Complex32_t *Al;
    int i,j;
    int M_ = m / mb;
    int N_ = n / nb;
    int bsiz = mb*nb;
    int psiz = m*nb;
    Quark_Task_Flags task_flags = Quark_Task_Flags_Initializer;
    QUARK_Task_Flag_Set(&task_flags, TASK_SEQUENCE, (intptr_t)sequence->quark_sequence);

    plasma_dynamic_spawn();
    for(j=0; j<N_; j++) {
        Al = &(A[psiz*j]);

        for(i=1; i<M_; i++) {

#ifdef TRACE_IPT
            char str[30];
            sprintf(str, "Foo2 C2RI %d", i*m*nb);
#endif
            QUARK_Insert_Task(plasma->quark, CORE_foo2_quark, &task_flags,
                              sizeof(PLASMA_Complex32_t)*psiz,  Al,           INPUT,
                              sizeof(PLASMA_Complex32_t)*bsiz, &(Al[i*bsiz]), INOUT,
#ifdef TRACE_IPT
                              30, str,   VALUE | TASKLABEL,
                              4, "red",  VALUE | TASKCOLOR,
#endif
                              0);
        }
    }

    return PLASMA_SUCCESS;
}

int ipt_cpanel2all(plasma_context_t *plasma, int m, int n, PLASMA_Complex32_t *A, int mb, int nb,
                   PLASMA_sequence *sequence, PLASMA_request *request)
{
    if (PLASMA_SCHEDULING != PLASMA_DYNAMIC_SCHEDULING)
        return PLASMA_SUCCESS;

    int i;
    int N_ = n / nb;
    Quark_Task_Flags task_flags = Quark_Task_Flags_Initializer;
    QUARK_Task_Flag_Set(&task_flags, TASK_SEQUENCE, (intptr_t)sequence->quark_sequence);

    if ( N_ > 1 ) {
        plasma_dynamic_spawn();
        for(i=1; i<N_; i++) {
#ifdef TRACE_IPT
            char str[30];
            sprintf(str, "Foo2 C2RI %d", i*m*nb);
#endif
            QUARK_Insert_Task(plasma->quark, CORE_foo2_quark, &task_flags,
                              sizeof(PLASMA_Complex32_t)*m*n,  A,            INOUT | GATHERV,
                              sizeof(PLASMA_Complex32_t)*m*nb, &(A[i*m*nb]), INPUT,
#ifdef TRACE_IPT
                              30, str,   VALUE | TASKLABEL,
                              4, "red",  VALUE | TASKCOLOR,
#endif
                              0);
        }
    }

    return PLASMA_SUCCESS;
}

int ipt_call2panel(plasma_context_t *plasma, int m, int n, PLASMA_Complex32_t *A, int mb, int nb,
                   PLASMA_sequence *sequence, PLASMA_request *request)
{
    if (PLASMA_SCHEDULING != PLASMA_DYNAMIC_SCHEDULING)
        return PLASMA_SUCCESS;

    int i;
    int N_ = n / nb;
    Quark_Task_Flags task_flags = Quark_Task_Flags_Initializer;
    QUARK_Task_Flag_Set(&task_flags, TASK_SEQUENCE, (intptr_t)sequence->quark_sequence);

    if ( N_ > 1 ) {
        plasma_dynamic_spawn();
        for(i=1; i<N_; i++) {
#ifdef TRACE_IPT
            char str[30];
            sprintf(str, "Foo2 C2RI %d", i*m*nb);
#endif
            QUARK_Insert_Task(plasma->quark, CORE_foo2_quark, &task_flags,
                              sizeof(PLASMA_Complex32_t)*m*n,  A,            INPUT,
                              sizeof(PLASMA_Complex32_t)*m*nb, &(A[i*m*nb]), INOUT,
#ifdef TRACE_IPT
                              30, str,   VALUE | TASKLABEL,
                              4, "red",  VALUE | TASKCOLOR,
#endif
                              0);
        }
    }
    return PLASMA_SUCCESS;
}

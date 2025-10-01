/**
 *
 * @file sgecfi2.c
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
 * @generated s Fri Apr  1 11:03:00 2016
 *
 **/

#include <sys/types.h>
#include "common.h"
#include "sgecfi2.h"

#define PLASMA_psgetmi2(idep, odep, storev, m, n, mb, nb, A)     \
    plasma_parallel_call_10(                                     \
        plasma_psgetmi2,                                         \
        PLASMA_enum,  (idep),                                    \
        PLASMA_enum,  (odep),                                    \
        PLASMA_enum,  (storev),                                  \
        int,  (m),                                               \
        int,  (n),                                               \
        int,  (mb),                                              \
        int,  (nb),                                              \
        float*, (A),                                \
        PLASMA_sequence*, sequence,                              \
        PLASMA_request*, request);

#define PLASMA_sshift(m, n, mb, nb, A)                          \
    plasma_sshift(plasma, (m), (n), (A),                        \
                  ( (n) / (nb) ), ( (m) / (mb) ), (nb), (mb),   \
                  sequence, request);

#define PLASMA_sshiftr(m, n, mb, nb, A)                         \
    plasma_sshift(plasma, (m), (n), (A),                        \
                  ( (n) / (nb) ), (nb), ( (m) / (mb) ), (mb),   \
                  sequence, request);

/** ****************************************************************************
 *
 * @ingroup InPlaceTransformation
 *
 *  ipt_scm2ccrb converts a matrix from CM format to CCRB format
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
int ipt_scm2ccrb(plasma_context_t *plasma, int m, int n, float *A, int mb, int nb,
                 PLASMA_sequence *sequence, PLASMA_request *request)
{
    if( (m == 0) || (n == 0) )
        return PLASMA_SUCCESS;

    PLASMA_sshift(m, n, mb, nb, A);
    ipt_spanel2tile(plasma, m, n, A, mb, nb, sequence, request);

    return PLASMA_SUCCESS;
}

int ipt_sccrb2cm(plasma_context_t *plasma, int m, int n, float *A, int mb, int nb,
                 PLASMA_sequence *sequence, PLASMA_request *request)
{
    if( (m == 0) || (n == 0) )
        return PLASMA_SUCCESS;

    ipt_stile2panel(plasma, m, n, A, mb, nb, sequence, request);
    PLASMA_sshiftr(m, n, mb, nb, A);

    return PLASMA_SUCCESS;
}

/*
 * Transpose each tile
 */
int ipt_sccrb2crrb(plasma_context_t *plasma, PLASMA_enum idep, PLASMA_enum odep, int m, int n, float *A, int mb, int nb,
                   PLASMA_sequence *sequence, PLASMA_request *request)
{
    if( (m == 0) || (n == 0) )
        return PLASMA_SUCCESS;

    PLASMA_psgetmi2(idep, odep, PlasmaColumnwise, m, n, mb, nb, A);

    return PLASMA_SUCCESS;
}

int ipt_scrrb2ccrb(plasma_context_t *plasma, PLASMA_enum idep, PLASMA_enum odep, int m, int n, float *A, int mb, int nb,
                   PLASMA_sequence *sequence, PLASMA_request *request)
{
    if( (m == 0) || (n == 0) )
        return PLASMA_SUCCESS;

    PLASMA_psgetmi2(idep, odep, PlasmaRowwise, n, m, nb, mb, A);

    return PLASMA_SUCCESS;
}

int ipt_srcrb2rrrb(plasma_context_t *plasma, PLASMA_enum idep, PLASMA_enum odep, int m, int n, float *A, int mb, int nb,
                   PLASMA_sequence *sequence, PLASMA_request *request)
{
    if( (m == 0) || (n == 0) )
        return PLASMA_SUCCESS;

    PLASMA_psgetmi2(idep, odep, PlasmaRowwise, m, n, mb, nb, A);

    return PLASMA_SUCCESS;
}

int ipt_srrrb2rcrb(plasma_context_t *plasma, PLASMA_enum idep, PLASMA_enum odep, int m, int n, float *A, int mb, int nb,
                   PLASMA_sequence *sequence, PLASMA_request *request)
{
    if( (m == 0) || (n == 0) )
        return PLASMA_SUCCESS;

    PLASMA_psgetmi2(idep, odep, PlasmaColumnwise, n, m, nb, mb, A);

    return PLASMA_SUCCESS;
}

/*
 * Transpose all tiles
 */
int ipt_sccrb2rcrb(plasma_context_t *plasma, int m, int n, float *A, int mb, int nb,
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

    plasma_sshift(plasma, m, n, A, 1, ( m / mb ), ( n / nb ), (mb*nb),
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
int ipt_scm2rcrb(plasma_context_t *plasma, int m, int n, float *A, int mb, int nb,
                 PLASMA_sequence *sequence, PLASMA_request *request)
{
    if( (m == 0) || (n == 0) )
        return PLASMA_SUCCESS;
    //ipt_scm2ccrb(  plasma, m, n, A, mb, nb, sequence, request);
    PLASMA_sshift(m, n, mb, nb, A);
    ipt_spanel2all(plasma, m, n, A, mb, nb, sequence, request);
    ipt_sccrb2rcrb(plasma, m, n, A, mb, nb, sequence, request);

    return PLASMA_SUCCESS;
}

int ipt_srcrb2cm(plasma_context_t *plasma, int m, int n, float *A, int mb, int nb,
                 PLASMA_sequence *sequence, PLASMA_request *request)
{
    if( (m == 0) || (n == 0) )
        return PLASMA_SUCCESS;
    ipt_srcrb2ccrb(plasma, m, n, A, mb, nb, sequence, request);
    ipt_sall2panel(plasma, m, n, A, mb, nb, sequence, request);
    //ipt_sccrb2cm(  plasma, m, n, A, mb, nb, sequence, request);
    PLASMA_sshiftr(m, n, mb, nb, A);

    return PLASMA_SUCCESS;
}

/*
 * Transpose each tile  + Transpose all tiles
 */
int ipt_sccrb2rrrb(plasma_context_t *plasma, int m, int n, float *A, int mb, int nb,
                   PLASMA_sequence *sequence, PLASMA_request *request)
{
    if( (m == 0) || (n == 0) )
        return PLASMA_SUCCESS;
    ipt_sccrb2rcrb(plasma, m, n, A, mb, nb, sequence, request);
    ipt_srcrb2rrrb(plasma, PlasmaIPT_All, PlasmaIPT_NoDep, m, n, A, mb, nb, sequence, request);

    return PLASMA_SUCCESS;
}

int ipt_srrrb2ccrb(plasma_context_t *plasma, int m, int n, float *A, int mb, int nb,
                 PLASMA_sequence *sequence, PLASMA_request *request)
{
    if( (m == 0) || (n == 0) )
        return PLASMA_SUCCESS;
    ipt_srrrb2rcrb(plasma, PlasmaIPT_NoDep, PlasmaIPT_All, m, n, A, mb, nb, sequence, request);
    ipt_srcrb2ccrb(plasma, m, n, A, mb, nb, sequence, request);

    return PLASMA_SUCCESS;
}

int ipt_srcrb2crrb(plasma_context_t *plasma, int m, int n, float *A, int mb, int nb,
                 PLASMA_sequence *sequence, PLASMA_request *request)
{
    if( (m == 0) || (n == 0) )
        return PLASMA_SUCCESS;
    ipt_srcrb2ccrb(plasma, m, n, A, mb, nb, sequence, request);
    ipt_sccrb2crrb(plasma, PlasmaIPT_All, PlasmaIPT_NoDep, m, n, A, mb, nb, sequence, request);

    return PLASMA_SUCCESS;
}

int ipt_scrrb2rcrb(plasma_context_t *plasma, int m, int n, float *A, int mb, int nb,
                   PLASMA_sequence *sequence, PLASMA_request *request)
{
    if( (m == 0) || (n == 0) )
        return PLASMA_SUCCESS;
    ipt_scrrb2ccrb(plasma, PlasmaIPT_NoDep, PlasmaIPT_All, m, n, A, mb, nb, sequence, request);
    ipt_sccrb2rcrb(plasma, m, n, A, mb, nb, sequence, request);

    return PLASMA_SUCCESS;
}

/*
 * Transpose each tile  + Shift inside panels
 */
int ipt_scm2crrb(plasma_context_t *plasma, int m, int n, float *A, int mb, int nb,
                 PLASMA_sequence *sequence, PLASMA_request *request)
{
    if( (m == 0) || (n == 0) )
        return PLASMA_SUCCESS;
    //ipt_scm2ccrb(  plasma, m, n, A, mb, nb, sequence, request);
    PLASMA_sshift(m, n, mb, nb, A);
    ipt_sccrb2crrb(plasma, PlasmaIPT_Panel, PlasmaIPT_NoDep, m, n, A, mb, nb, sequence, request);

    return PLASMA_SUCCESS;
}

int ipt_scrrb2cm(plasma_context_t *plasma, int m, int n, float *A, int mb, int nb,
                 PLASMA_sequence *sequence, PLASMA_request *request)
{
    if( (m == 0) || (n == 0) )
        return PLASMA_SUCCESS;
    ipt_scrrb2ccrb(plasma, PlasmaIPT_NoDep, PlasmaIPT_Panel, m, n, A, mb, nb, sequence, request);
    //ipt_sccrb2cm(  plasma, m, n, A, mb, nb, sequence, request);
    PLASMA_sshiftr(m, n, mb, nb, A);

    return PLASMA_SUCCESS;
}

int ipt_srm2rcrb(plasma_context_t *plasma, int m, int n, float *A, int mb, int nb,
                 PLASMA_sequence *sequence, PLASMA_request *request)
{
    if( (m == 0) || (n == 0) )
        return PLASMA_SUCCESS;
    ipt_srm2rrrb(  plasma, m, n, A, mb, nb, sequence, request);
    ipt_srrrb2rcrb(plasma, PlasmaIPT_Panel, PlasmaIPT_NoDep, m, n, A, mb, nb, sequence, request);

    return PLASMA_SUCCESS;
}

int ipt_srcrb2rm(plasma_context_t *plasma, int m, int n, float *A, int mb, int nb,
                 PLASMA_sequence *sequence, PLASMA_request *request)
{
    if( (m == 0) || (n == 0) )
        return PLASMA_SUCCESS;
    ipt_srcrb2rrrb(plasma, PlasmaIPT_NoDep, PlasmaIPT_Panel, m, n, A, mb, nb, sequence, request);
    ipt_srrrb2rm(  plasma, m, n, A, mb, nb, sequence, request);

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
int ipt_scm2rrrb(plasma_context_t *plasma, int m, int n, float *A, int mb, int nb,
                 PLASMA_sequence *sequence, PLASMA_request *request)
{
    if( (m == 0) || (n == 0) )
        return PLASMA_SUCCESS;
    //ipt_scm2ccrb(  plasma, m, n, A, mb, nb, sequence, request);
    PLASMA_sshift(m, n, mb, nb, A);
    ipt_sccrb2crrb(plasma, PlasmaIPT_Panel, PlasmaIPT_All, m, n, A, mb, nb, sequence, request);
    ipt_scrrb2rrrb(plasma, m, n, A, mb, nb, sequence, request);

    return PLASMA_SUCCESS;
}

int ipt_srrrb2cm(plasma_context_t *plasma, int m, int n, float *A, int mb, int nb,
                 PLASMA_sequence *sequence, PLASMA_request *request)
{
    if( (m == 0) || (n == 0) )
        return PLASMA_SUCCESS;
    ipt_srrrb2crrb(plasma, m, n, A, mb, nb, sequence, request);
    ipt_scrrb2ccrb(plasma, PlasmaIPT_All, PlasmaIPT_Panel, m, n, A, mb, nb, sequence, request);
    //ipt_sccrb2cm(  plasma, m, n, A, mb, nb, sequence, request);
    PLASMA_sshiftr(m, n, mb, nb, A);

    return PLASMA_SUCCESS;
}

int ipt_sccrb2rm(plasma_context_t *plasma, int m, int n, float *A, int mb, int nb,
                 PLASMA_sequence *sequence, PLASMA_request *request)
{
    if( (m == 0) || (n == 0) )
        return PLASMA_SUCCESS;
    ipt_sccrb2rcrb(plasma, m, n, A, mb, nb, sequence, request);
    ipt_srcrb2rrrb(plasma, PlasmaIPT_All, PlasmaIPT_Panel, m, n, A, mb, nb, sequence, request);
    ipt_srrrb2rm(  plasma, m, n, A, mb, nb, sequence, request);

    return PLASMA_SUCCESS;
}

int ipt_srm2ccrb(plasma_context_t *plasma, int m, int n, float *A, int mb, int nb,
                 PLASMA_sequence *sequence, PLASMA_request *request)
{
    if( (m == 0) || (n == 0) )
        return PLASMA_SUCCESS;
    ipt_srm2rrrb(  plasma, m, n, A, mb, nb, sequence, request);
    ipt_srrrb2rcrb(plasma, PlasmaIPT_Panel, PlasmaIPT_All, m, n, A, mb, nb, sequence, request);
    ipt_srcrb2ccrb(plasma, m, n, A, mb, nb, sequence, request);

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
int ipt_scm2rm(plasma_context_t *plasma, int m, int n, float *A, int mb, int nb,
               PLASMA_sequence *sequence, PLASMA_request *request)
{
    if( (m == 0) || (n == 0) )
        return PLASMA_SUCCESS;
    //ipt_scm2ccrb(  plasma, m, n, A, mb, nb, sequence, request);
    PLASMA_sshift(m, n, mb, nb, A);
    ipt_spanel2all(plasma, m, n, A, mb, nb, sequence, request);
    ipt_sccrb2rcrb(plasma, m, n, A, mb, nb, sequence, request);
    ipt_srcrb2rrrb(plasma, PlasmaIPT_All, PlasmaIPT_Panel, m, n, A, mb, nb, sequence, request);
    ipt_srrrb2rm(  plasma, m, n, A, mb, nb, sequence, request);

    return PLASMA_SUCCESS;
}

int ipt_srm2cm(plasma_context_t *plasma, int m, int n, float *A, int mb, int nb,
                 PLASMA_sequence *sequence, PLASMA_request *request)
{
    if( (m == 0) || (n == 0) )
        return PLASMA_SUCCESS;
    ipt_srm2rrrb(  plasma, m, n, A, mb, nb, sequence, request);
    ipt_srrrb2rcrb(plasma, PlasmaIPT_Panel, PlasmaIPT_All, m, n, A, mb, nb, sequence, request);
    ipt_srcrb2ccrb(plasma, m, n, A, mb, nb, sequence, request);
    ipt_sall2panel(plasma, m, n, A, mb, nb, sequence, request);
    //ipt_sccrb2cm(  plasma, m, n, A, mb, nb, sequence, request);
    PLASMA_sshiftr(m, n, mb, nb, A);

    return PLASMA_SUCCESS;
}


/** ****************************************************************************
 *
 *  Barriers
 *
 *******************************************************************************/

int ipt_stile2panel( plasma_context_t *plasma, int m, int n, float *A, int mb, int nb,
                     PLASMA_sequence *sequence, PLASMA_request *request)
{
    if( PLASMA_SCHEDULING != PLASMA_DYNAMIC_SCHEDULING )
        return PLASMA_SUCCESS;

    float *Al;
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
                              sizeof(float)*psiz,  Al,           INOUT | GATHERV,
                              sizeof(float)*bsiz, &(Al[i*bsiz]), INOUT,
#ifdef TRACE_IPT
                              30, str,   VALUE | TASKLABEL,
                              4, "red",  VALUE | TASKCOLOR,
#endif
                              0);
        }
    }

    return PLASMA_SUCCESS;
}

int ipt_spanel2tile( plasma_context_t *plasma, int m, int n, float *A, int mb, int nb,
                     PLASMA_sequence *sequence, PLASMA_request *request )
{
    if( PLASMA_SCHEDULING != PLASMA_DYNAMIC_SCHEDULING )
        return PLASMA_SUCCESS;

    float *Al;
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
                              sizeof(float)*psiz,  Al,           INPUT,
                              sizeof(float)*bsiz, &(Al[i*bsiz]), INOUT,
#ifdef TRACE_IPT
                              30, str,   VALUE | TASKLABEL,
                              4, "red",  VALUE | TASKCOLOR,
#endif
                              0);
        }
    }

    return PLASMA_SUCCESS;
}

int ipt_spanel2all(plasma_context_t *plasma, int m, int n, float *A, int mb, int nb,
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
                              sizeof(float)*m*n,  A,            INOUT | GATHERV,
                              sizeof(float)*m*nb, &(A[i*m*nb]), INPUT,
#ifdef TRACE_IPT
                              30, str,   VALUE | TASKLABEL,
                              4, "red",  VALUE | TASKCOLOR,
#endif
                              0);
        }
    }

    return PLASMA_SUCCESS;
}

int ipt_sall2panel(plasma_context_t *plasma, int m, int n, float *A, int mb, int nb,
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
                              sizeof(float)*m*n,  A,            INPUT,
                              sizeof(float)*m*nb, &(A[i*m*nb]), INOUT,
#ifdef TRACE_IPT
                              30, str,   VALUE | TASKLABEL,
                              4, "red",  VALUE | TASKCOLOR,
#endif
                              0);
        }
    }
    return PLASMA_SUCCESS;
}

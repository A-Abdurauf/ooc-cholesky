/**
 * @file pslag2c.c
 *
 *  PLASMA computational routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Grégoire Pichon
 * @author Azzam Haidar
 * @date 2014-07
 * @generated c Fri Apr  1 11:03:00 2016
 *
 **/
#include "common.h"

/***************************************************************************//**
 *
 * @ingroup CORE_PLASMA_Complex32_t
 *
 *  PLASMA_pslag2c - copy the float real precision matrix Q in the float complex
 *  precision matrix Z
 *
 *******************************************************************************
 *
 * @param[in] n
 *          n specifies the dimension of the matrices
 *
 * @param[out] Z
 *          On exit, the eigenvectors stored in complex precison
 *
 * @param[in] LDZ
 *          LDZ specifies the leading dimension of Z
 *
 * @param[in] Q
 *          On entry, the eigenvectors stored in real precision
 *
 * @param[in] sequence
 *          Identifies the sequence of function calls that this call belongs to
 *          (for completion checks and exception handling purposes).
 *
 * @param[out] request
 *          Identifies this function call (for exception handling purposes).
 *
***************************************************************************/
void plasma_pslag2c_quark(int m, int n,
                          const float *R, int ldr,
                          PLASMA_Complex32_t *Z, int ldz,
                          PLASMA_sequence *sequence, PLASMA_request *request)
{
    plasma_context_t *plasma;
    Quark_Task_Flags task_flags = Quark_Task_Flags_Initializer;
    int i, nb, tempnn;

    plasma = plasma_context_self();
    if (sequence->status != PLASMA_SUCCESS)
        return;
    QUARK_Task_Flag_Set(&task_flags, TASK_SEQUENCE, (intptr_t)sequence->quark_sequence);

    nb = plasma->ev_tasknb;
    for(i=0; i<n; i+=nb){
        tempnn = min(nb, n-i);
        QUARK_CORE_slag2c(plasma->quark, &task_flags,
                          m, tempnn, R+i*ldr, ldr, Z+i*ldz, ldz);
    }
}

/**
 * @file psstedc.c
 *
 *  PLASMA computational routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Gregoire Pichon
 * @date 2014-07
 * @generated s Fri Apr  1 11:03:00 2016
 *
 **/
#include "common.h"

#undef COMPLEX
#define REAL

/**
 ******************************************************************************
 *
 * @ingroup CORE_float
 *
 *  plasma_psstedc - solves the symmetric tridiagonal eigensystem using LAPACK's
 *  divide and conquer
 *
 *******************************************************************************
 *
 * @param[in] compz
 *          = PlasmaNoVec: computes eigenvalues only.
 *          = PlasmaIVec: computes eigenpairs of the symmetric tridiagonal matrix
 *          = PlasmaVec: computes eigenpairs of the original matrix
 *
 * @param[in] n
 *          n specifies the dimension of the symmetric tridiagonal matrix
 *
 * @param[in,out] D
 *          On entry, D contains the diagonal elements of the tridiagonal matrix.
 *          On exit, D contains the eigenvalues sorted into increasing order.
 *
 * @param[in] E
 *          On entry, E contains the extra-diagonal elements of the tridiagonal matrix
 *
 * @param[in,out] Z
 *          On entry, Z must be set to 0.
 *          On exit, Z contains the eigenvectors.
 *
 * @param[in] LDZ
 *          LDZ specifies the leading direction of Z
 *
 * @param[in] sequence
 *          Identifies the sequence of function calls that this call belongs to
 *          (for completion checks and exception handling purposes).
 *
 * @param[out] request
 *          Identifies this function call (for exception handling purposes).
 *
 *******************************************************************************/

void plasma_psstedc_quark(PLASMA_enum compz, int n, float *D, float *E,
                          float *Z, int LDZ,
                          PLASMA_sequence *sequence, PLASMA_request *request)
{
    plasma_context_t *plasma;
    plasma = plasma_context_self();

    Quark_Task_Flags task_flags = Quark_Task_Flags_Initializer;

    if (sequence->status != PLASMA_SUCCESS){
        return;
    }

    QUARK_Task_Flag_Set(&task_flags, TASK_SEQUENCE, (intptr_t)sequence->quark_sequence);

    QUARK_CORE_sstedc(plasma->quark, &task_flags,
                      compz, n, D, E, Z, LDZ);
}

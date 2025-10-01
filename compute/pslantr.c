/**
 *
 * @file pslantr.c
 *
 *  PLASMA auxiliary routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @generated s Fri Apr  1 11:02:57 2016
 *
 **/
#include <stdlib.h>
#include <math.h>
#include "common.h"

#define A(m, n, i, j, ldt)  (BLKADDR(A, float, m, n)+((j)*(ldt)+(i)))

/***************************************************************************//**
 *
 **/
void plasma_pslantr_quark(PLASMA_enum norm, PLASMA_enum uplo, PLASMA_enum diag,
                          PLASMA_desc A, float *work, float *result,
                          PLASMA_sequence *sequence, PLASMA_request *request)
{
    plasma_context_t *plasma;
    Quark_Task_Flags task_flags = Quark_Task_Flags_Initializer;

    float* lwork;
    int X, X1, X2, Y, Y1, Y2;
    int ldam, ldan;
    int m, n, k, minMNT;
    int szeW, pos;

    plasma = plasma_context_self();
    if (sequence->status != PLASMA_SUCCESS)
        return;
    QUARK_Task_Flag_Set(&task_flags, TASK_SEQUENCE, (intptr_t)sequence->quark_sequence);
    plasma_profile_by_function( &task_flags, LANGE );

    minMNT = min( A.mt, A.nt );

    *result = 0.0;
    switch ( norm ) {
    /*
     *  PlasmaMaxNorm
     */
    case PlasmaMaxNorm:
        szeW = minMNT*(minMNT+1)/2;
        if ( uplo == PlasmaLower )
            szeW += max((A.mt - A.nt), 0) * A.nt;
        else
            szeW += max((A.nt - A.mt), 0) * A.mt;

        pos = 0;
        lwork = (float *)plasma_shared_alloc(plasma, szeW, PlasmaRealDouble);
        memset(lwork, 0, szeW*sizeof(float));

        /*
         *  PlasmaLower
         */
        if (uplo == PlasmaLower) {
            for(n = 0; n < minMNT; n++) {
                X1 = n == 0      ?  A.i       %A.mb   : 0;
                X2 = n == A.mt-1 ? (A.i+A.m-1)%A.mb+1 : A.mb;
                X = X2 - X1;

                Y1 = n == 0      ?  A.j       %A.nb   : 0;
                Y2 = n == A.nt-1 ? (A.j+A.n-1)%A.nb+1 : A.nb;
                Y = Y2 - Y1;

                ldan = BLKLDD(A, n);
                QUARK_CORE_slantr_f1(
                    plasma->quark, &task_flags,
                    PlasmaMaxNorm, uplo, diag, X, Y,
                    A(n, n, X1, Y1, ldan), ldan, ldan*Y,
                    0, lwork + pos,
                    result, 1);
                pos++;

                for(m = n+1; m < A.mt; m++) {
                    X = m == A.mt-1 ? (A.i+A.m-1)%A.mb+1 : A.mb;
                    ldam = BLKLDD(A, m);

                    QUARK_CORE_slange_f1(
                        plasma->quark, &task_flags,
                        PlasmaMaxNorm, X, Y,
                        A(m, n, 0, Y1, ldam), ldam, ldam*Y,
                        0, lwork + pos,
                        result, 1);
                    pos++;
                }
            }
        }
        /*
         *  PlasmaUpper
         */
        else {
            for(m = 0; m < minMNT; m++) {
                X1 = m == 0      ?  A.i       %A.mb   : 0;
                X2 = m == A.mt-1 ? (A.i+A.m-1)%A.mb+1 : A.mb;
                X = X2 - X1;

                Y1 = m == 0      ?  A.j       %A.nb   : 0;
                Y2 = m == A.nt-1 ? (A.j+A.n-1)%A.nb+1 : A.nb;
                Y = Y2 - Y1;

                ldam = BLKLDD(A, m);
                QUARK_CORE_slantr_f1(
                    plasma->quark, &task_flags,
                    PlasmaMaxNorm, uplo, diag, X, Y,
                    A(m, m, X1, Y1, ldam), ldam, ldam*Y,
                    0, lwork + pos,
                    result, 1);
                pos++;

                for(n = m+1; n < A.nt; n++) {
                    Y = n == A.nt-1 ? (A.j+A.n-1)%A.nb+1 : A.nb;
                    QUARK_CORE_slange_f1(
                        plasma->quark, &task_flags,
                        PlasmaMaxNorm, X, Y,
                        A(m, n, X1, 0, ldam), ldam, ldam*Y,
                        0, lwork + pos,
                        result, 1);
                    pos++;
                }
            }
        }
        QUARK_CORE_slange(
            plasma->quark, &task_flags,
            PlasmaMaxNorm, szeW, 1,
            lwork, 1, szeW,
            0, result);

        QUARK_CORE_free(
            plasma->quark, &task_flags,
            lwork, szeW*sizeof(float));
        break;
    /*
     *  PlasmaOneNorm
     */
    case PlasmaOneNorm:
        lwork = (float *)plasma_shared_alloc(plasma, A.n, PlasmaRealDouble);
        memset(lwork, 0, A.n*sizeof(float));

        /*
         *  PlasmaUpper
         */
        if (uplo == PlasmaUpper) {
            for(m = 0; m < minMNT; m++) {
                X1 = m == 0      ?  A.i       %A.mb   : 0;
                X2 = m == A.mt-1 ? (A.i+A.m-1)%A.mb+1 : A.mb;
                X = X2 - X1;

                Y1 = m == 0      ?  A.j       %A.nb   : 0;
                Y2 = m == A.nt-1 ? (A.j+A.n-1)%A.nb+1 : A.nb;
                Y = Y2 - Y1;

                ldam = BLKLDD(A, m);

                QUARK_CORE_strasm_f1(
                    plasma->quark, &task_flags,
                    PlasmaColumnwise, uplo, diag, X, Y,
                    A(m, m, X1, Y1, ldam), ldam, ldam*Y,
                    lwork + m*A.nb, A.nb,
                    result, 1);

                for(n = m+1; n < A.nt; n++) {
                    Y = n == A.nt-1 ? (A.j+A.n-1)%A.nb+1 : A.nb;

                    QUARK_CORE_sasum_f1(
                        plasma->quark, &task_flags,
                        PlasmaColumnwise, PlasmaUpperLower, X, Y,
                        A(m, n, X1, 0, A.mb), ldam, ldam*Y,
                        lwork + n*A.nb, A.nb,
                        result, 1);
                }
            }
        }
        /*
         *  PlasmaLower
         */
        else {
            for(n = 0; n < minMNT; n++) {
                X1 = n == 0      ?  A.i       %A.mb   : 0;
                X2 = n == A.mt-1 ? (A.i+A.m-1)%A.mb+1 : A.mb;
                X = X2 - X1;

                Y1 = n == 0      ?  A.j       %A.nb   : 0;
                Y2 = n == A.nt-1 ? (A.j+A.n-1)%A.nb+1 : A.nb;
                Y = Y2 - Y1;

                ldan = BLKLDD(A, n);
                QUARK_CORE_strasm_f1(
                    plasma->quark, &task_flags,
                    PlasmaColumnwise, uplo, diag, X, Y,
                    A(n, n, X1, Y1, ldan), ldan, ldan*Y,
                    lwork + n*A.nb, A.nb,
                    result, 1);

                for(m = n+1; m < A.mt; m++) {
                    X = m == A.mt-1 ? (A.i+A.m-1)%A.mb+1 : A.mb;
                    ldam = BLKLDD(A, m);

                    QUARK_CORE_sasum_f1(
                        plasma->quark, &task_flags,
                        PlasmaColumnwise, PlasmaUpperLower, X, Y,
                        A(m, n, 0, Y1, ldam), ldam, ldam*Y,
                        lwork + n*A.nb, A.nb,
                        result, 1);
                }
            }
        }
        QUARK_CORE_slange(
            plasma->quark, &task_flags,
            PlasmaMaxNorm, A.n, 1,
            lwork, 1, A.n,
            0, result);

        QUARK_CORE_free(
            plasma->quark, &task_flags,
            lwork, A.n*sizeof(float) );
        break;

    /*
     *  PlasmaInfNorm
     */
    case PlasmaInfNorm:
        lwork = (float *)plasma_shared_alloc(plasma, A.m, PlasmaRealDouble);
        memset(lwork, 0, A.m*sizeof(float));

        /*
         *  PlasmaLower
         */
        if (uplo == PlasmaLower) {
            for(n = 0; n < minMNT; n++) {
                X1 = n == 0      ?  A.i       %A.mb   : 0;
                X2 = n == A.mt-1 ? (A.i+A.m-1)%A.mb+1 : A.mb;
                X = X2 - X1;

                Y1 = n == 0      ?  A.j       %A.nb   : 0;
                Y2 = n == A.nt-1 ? (A.j+A.n-1)%A.nb+1 : A.nb;
                Y = Y2 - Y1;

                ldan = BLKLDD(A, n);
                QUARK_CORE_strasm_f1(
                    plasma->quark, &task_flags,
                    PlasmaRowwise, uplo, diag, X, Y,
                    A(n, n, X1, Y1, ldan), ldan, ldan*X,
                    lwork + n*A.mb, A.mb,
                    result, 1);

                for(m = n+1; m < A.mt; m++) {
                    X = m == A.mt-1 ? (A.i+A.m-1)%A.mb+1 : A.mb;
                    ldam = BLKLDD(A, m);

                    QUARK_CORE_sasum_f1(
                        plasma->quark, &task_flags,
                        PlasmaRowwise, PlasmaUpperLower,
                        X, Y,
                        A(m, n, 0, Y1, ldam), ldam, ldam*Y,
                        lwork + m*A.mb, A.mb,
                        result, 1);
                }
            }
        }
        /*
         *  PlasmaUpper
         */
        else {
            for(m = 0; m < minMNT; m++) {
                X1 = m == 0      ?  A.i       %A.mb   : 0;
                X2 = m == A.mt-1 ? (A.i+A.m-1)%A.mb+1 : A.mb;
                X = X2 - X1;

                Y1 = m == 0      ?  A.j       %A.nb   : 0;
                Y2 = m == A.nt-1 ? (A.j+A.n-1)%A.nb+1 : A.nb;
                Y = Y2 - Y1;

                ldam = BLKLDD(A, m);
                QUARK_CORE_strasm_f1(
                    plasma->quark, &task_flags,
                    PlasmaRowwise, uplo, diag, X, Y,
                    A(m, m, X1, Y1, ldam), ldam, ldam*X,
                    lwork + m*A.mb, A.mb,
                    result, 1);

                for(n = m+1; n < A.nt; n++) {
                    Y = n == A.nt-1 ? (A.j+A.n-1)%A.nb+1 : A.nb;
                    QUARK_CORE_sasum_f1(
                        plasma->quark, &task_flags,
                        PlasmaRowwise, PlasmaUpperLower,
                        X, Y,
                        A(m, n, X1, 0, ldam), ldam, ldam*Y,
                        lwork + m*A.mb, A.mb,
                        result, 1);
                }
            }
        }
        QUARK_CORE_slange(
            plasma->quark, &task_flags,
            PlasmaMaxNorm, A.m, 1,
            lwork, 1, A.m,
            0, result);

        QUARK_CORE_free(
            plasma->quark, &task_flags,
            lwork, A.m*sizeof(float));
        break;
    /*
     *  PlasmaFrobeniusNorm
     */
    case PlasmaFrobeniusNorm:
        szeW = 2*(PLASMA_SIZE);
        lwork = (float*)plasma_shared_alloc(plasma, szeW, PlasmaRealDouble);

        for(m = 0; m < PLASMA_SIZE; m++) {
            lwork[2*m  ] = 0.;
            lwork[2*m+1] = 1.;
        }

        k = 0;
        /*
         *  PlasmaLower
         */
        if (uplo == PlasmaLower) {
            for(n = 0; n < minMNT; n++) {
                X1 = n == 0      ?  A.i       %A.mb   : 0;
                X2 = n == A.mt-1 ? (A.i+A.m-1)%A.mb+1 : A.mb;
                X = X2 - X1;

                Y1 = n == 0      ?  A.j       %A.nb   : 0;
                Y2 = n == A.nt-1 ? (A.j+A.n-1)%A.nb+1 : A.nb;
                Y = Y2 - Y1;

                ldan = BLKLDD(A, n);

                QUARK_CORE_strssq_f1(
                    plasma->quark, &task_flags,
                    uplo, diag, X, Y,
                    A(n, n, X1, Y1, ldan), ldan,
                    lwork + 2*k,
                    lwork + 2*k + 1,
                    result, 1, INOUT | GATHERV );
                k++;
                k = k % PLASMA_SIZE;

                for(m = n+1; m < A.mt; m++) {
                    X = m == A.mt-1 ? (A.i+A.m-1)%A.mb+1 : A.mb;
                    ldam = BLKLDD(A, m);

                    QUARK_CORE_sgessq_f1(
                        plasma->quark, &task_flags,
                        X, Y,
                        A(m, n, 0, Y1, ldam), ldam,
                        lwork + 2*k,
                        lwork + 2*k + 1,
                        result, 1, INOUT | GATHERV );
                    k++;
                    k = k % PLASMA_SIZE;
                }
            }
        }
        /*
         *  PlasmaUpper
         */
        else {
            for(m = 0; m < minMNT; m++) {
                X1 = m == 0      ?  A.i       %A.mb   : 0;
                X2 = m == A.mt-1 ? (A.i+A.m-1)%A.mb+1 : A.mb;
                X = X2 - X1;

                Y1 = m == 0      ?  A.j       %A.nb   : 0;
                Y2 = m == A.nt-1 ? (A.j+A.n-1)%A.nb+1 : A.nb;
                Y = Y2 - Y1;

                ldam = BLKLDD(A, m);

                QUARK_CORE_strssq_f1(
                    plasma->quark, &task_flags,
                    uplo, diag, X, Y,
                    A(m, m, X1, Y1, ldam), ldam,
                    lwork + 2*k,
                    lwork + 2*k + 1,
                    result, 1, INOUT | GATHERV );
                k++;
                k = k % PLASMA_SIZE;

                for(n = m+1; n < A.nt; n++) {
                    Y = n == A.nt-1 ? (A.j+A.n-1)%A.nb+1 : A.nb;

                    QUARK_CORE_sgessq_f1(
                        plasma->quark, &task_flags,
                        X, Y,
                        A(m, n, X1, 0, ldam), ldam,
                        lwork + 2*k,
                        lwork + 2*k + 1,
                        result, 1, INOUT | GATHERV );

                    k++;
                    k = k % PLASMA_SIZE;
                }
            }
        }
        QUARK_CORE_splssq(
            plasma->quark, &task_flags,
            PLASMA_SIZE, lwork, result );

        QUARK_CORE_free(
            plasma->quark, &task_flags,
            lwork, szeW*sizeof(float));
    default:;
    }
}

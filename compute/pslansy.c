/**
 *
 * @file pslansy.c
 *
 *  PLASMA auxiliary routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Emmanuel Agullo
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
void plasma_pslansy(plasma_context_t *plasma)
{
    PLASMA_enum norm;
    PLASMA_enum uplo;
    PLASMA_desc A;
    float *work;
    float *result;
    PLASMA_sequence *sequence;
    PLASMA_request *request;

    int m, n;
    int next_m;
    int next_n;
    int ldam, ldan;
    int step, lrank;
    int X, X1, X2, Y, Y1, Y2;

    float* lwork;
    float normtmp, normtmp2;
    float *scale, *sumsq;
    float scale2, sumsq2;

    plasma_unpack_args_7(norm, uplo, A, work, result, sequence, request);
    *result = 0.0;

    if (PLASMA_RANK == 0) {
        if ( norm == PlasmaFrobeniusNorm ) {
            memset(work, 0, 2*PLASMA_SIZE*sizeof(float));
        } else {
            memset(work, 0,   PLASMA_SIZE*sizeof(float));
        }
    }
    ss_init(PLASMA_SIZE, 1, 0);

    switch (norm) {
    /*
     *  PlasmaMaxNorm
     */
    case PlasmaMaxNorm:
        n = 0;
        m = PLASMA_RANK;
        while (m >= A.mt && n < A.nt) {
            n++;
            m = m-A.mt+n;
        }

        while (n < A.nt) {
            next_m = m;
            next_n = n;

            next_m += PLASMA_SIZE;
            while (next_m >= A.mt && next_n < A.nt) {
                next_n++;
                next_m = next_m-A.mt+next_n;
            }

            if (m == n) {
                X1 = m == 0      ?  A.i       %A.mb   : 0;
                X2 = m == A.mt-1 ? (A.i+A.m-1)%A.mb+1 : A.mb;
                X = X2 - X1;

                ldam = BLKLDD(A, m);
                CORE_slansy(PlasmaMaxNorm, uplo, X, A(m, n, X1, X1, ldam), ldam, NULL, &normtmp);
            }
            else {
                /*
                 *  PlasmaLower
                 */
                if (uplo == PlasmaLower) {
                    X1 = m == 0      ?  A.i       %A.mb   : 0;
                    X2 = m == A.mt-1 ? (A.i+A.m-1)%A.mb+1 : A.mb;
                    X = X2 - X1;

                    Y1 = n == 0      ?  A.j       %A.nb   : 0;
                    Y2 = n == A.nt-1 ? (A.j+A.n-1)%A.nb+1 : A.nb;
                    Y = Y2 - Y1;

                    ldam = BLKLDD(A, m);
                    CORE_slange(PlasmaMaxNorm, X, Y, A(m, n, X1, Y1, ldam), ldam, NULL, &normtmp);
                }
                /*
                 *  PlasmaUpper
                 */
                else {
                    X1 = n == 0      ?  A.i       %A.mb   : 0;
                    X2 = n == A.mt-1 ? (A.i+A.m-1)%A.mb+1 : A.mb;
                    X = X2 - X1;

                    Y1 = m == 0      ?  A.j       %A.nb   : 0;
                    Y2 = m == A.nt-1 ? (A.j+A.n-1)%A.nb+1 : A.nb;
                    Y = Y2 - Y1;

                    ldan = BLKLDD(A, n);
                    CORE_slange(PlasmaMaxNorm, X, Y, A(n, m, X1, Y1, ldan), ldan, NULL, &normtmp);
                }
            }
            if ( normtmp > work[PLASMA_RANK] )
                work[PLASMA_RANK] = normtmp;

            m = next_m;
            n = next_n;
        }
        ss_cond_set(PLASMA_RANK, 0, 1);
        break;
    /*
     *  PlasmaOneNorm / PlasmaInfNorm
     */
    case PlasmaOneNorm:
    case PlasmaInfNorm:
        m = PLASMA_RANK;
        normtmp2 = 0.0;
        lwork = (float*)plasma_private_alloc(plasma, A.mb, PlasmaRealDouble);

        while (m < A.mt) {
            X1 = m == 0      ?  A.i       %A.mb   : 0;
            X2 = m == A.mt-1 ? (A.i+A.m-1)%A.mb+1 : A.mb;
            X = X2 - X1;

            ldam = BLKLDD(A, m);
            memset(lwork, 0, A.mb*sizeof(float));
            /*
             *  PlasmaLower
             */
            if (uplo == PlasmaLower) {
                for (n = 0; n < m; n++) {
                    Y1 = n == 0 ? A.j%A.nb : 0;
                    Y = A.nb - Y1;
                    CORE_sasum(PlasmaRowwise, PlasmaUpperLower, X, Y, A(m, n, X1, Y1, ldam), ldam, lwork);
                }
                CORE_sasum(PlasmaRowwise, uplo, X, X, A(m, m, X1, X1, ldam), ldam, lwork);

                for (n = m+1; n < A.mt; n++) {
                    Y = n == A.nt-1 ? (A.j+A.n-1)%A.nb+1 : A.nb;
                    ldan = BLKLDD(A, n);
                    CORE_sasum(PlasmaColumnwise, PlasmaUpperLower, Y, X, A(n, m, 0, X1, ldan), ldan, lwork);
                }
            }
            /*
             *  PlasmaUpper
             */
            else {
                for (n = 0; n < m; n++) {
                    Y1 = n == 0 ?  A.j%A.nb : 0;
                    Y = A.nb - Y1;
                    CORE_sasum(PlasmaColumnwise, PlasmaUpperLower, Y, X, A(n, m, Y1, X1, A.nb), A.nb, lwork);
                }
                CORE_sasum(PlasmaRowwise, uplo, X, X, A(m, m, X1, X1, ldam), ldam, lwork);

                for (n = m+1; n < A.mt; n++) {
                    Y = n == A.nt-1 ? (A.j+A.n-1)%A.nb+1 : A.nb;
                    CORE_sasum(PlasmaRowwise, PlasmaUpperLower, X, Y, A(m, n, X1, 0, ldam), ldam, lwork);
                }
            }
            CORE_slange(PlasmaMaxNorm, X, 1, lwork, 1, NULL, &normtmp);
            if ( normtmp > normtmp2 )
                normtmp2 = normtmp;

            m += PLASMA_SIZE;
        }
        work[PLASMA_RANK] = normtmp2;
        ss_cond_set(PLASMA_RANK, 0, 1);
        plasma_private_free(plasma, lwork);
        break;
    /*
     *  PlasmaFrobeniusNorm
     */
    case PlasmaFrobeniusNorm:
        scale = work + 2 * PLASMA_RANK;
        sumsq = work + 2 * PLASMA_RANK + 1;

        *scale = 0.;
        *sumsq = 1.;

        n = 0;
        m = PLASMA_RANK;
        while (m >= A.mt && n < A.nt) {
            n++;
            m = m-A.mt+n;
        }

        while (n < A.nt) {
            next_m = m;
            next_n = n;

            next_m += PLASMA_SIZE;
            while (next_m >= A.mt && next_n < A.nt) {
                next_n++;
                next_m = next_m-A.mt+next_n;
            }

            if (m == n) {
                X1 = m == 0      ?  A.i       %A.mb   : 0;
                X2 = m == A.mt-1 ? (A.i+A.m-1)%A.mb+1 : A.mb;
                X = X2 - X1;

                ldam = BLKLDD(A, m);
                CORE_ssyssq( uplo, X, A(m, n, X1, X1, ldam), ldam, scale, sumsq );
            }
            else {
                scale2 = 0.;
                sumsq2 = 1.;

                /*
                 *  PlasmaLower
                 */
                if (uplo == PlasmaLower) {
                    X1 = m == 0      ?  A.i       %A.mb   : 0;
                    X2 = m == A.mt-1 ? (A.i+A.m-1)%A.mb+1 : A.mb;
                    X = X2 - X1;

                    Y1 = n == 0      ?  A.j       %A.nb   : 0;
                    Y2 = n == A.nt-1 ? (A.j+A.n-1)%A.nb+1 : A.nb;
                    Y = Y2 - Y1;

                    ldam = BLKLDD(A, m);
                    CORE_sgessq( X, Y, A(m, n, X1, Y1, ldam), ldam, &scale2, &sumsq2 );
                }
                /*
                 *  PlasmaUpper
                 */
                else {
                    X1 = n == 0      ?  A.i       %A.mb   : 0;
                    X2 = n == A.mt-1 ? (A.i+A.m-1)%A.mb+1 : A.mb;
                    X = X2 - X1;

                    Y1 = m == 0      ?  A.j       %A.nb   : 0;
                    Y2 = m == A.nt-1 ? (A.j+A.n-1)%A.nb+1 : A.nb;
                    Y = Y2 - Y1;

                    ldan = BLKLDD(A, n);
                    CORE_sgessq( X, Y, A(n, m, X1, Y1, ldan), ldan, &scale2, &sumsq2 );
                }
                sumsq2 *= 2.;

                if ( scale2 != 0. ){
                    if ( *scale < scale2 ) {
                        *sumsq = sumsq2 + (*sumsq) * ( *scale / scale2 ) * ( *scale / scale2 );
                        *scale = scale2;
                    } else {
                        *sumsq = *sumsq + sumsq2 * ( scale2 / *scale ) *  ( scale2 / *scale );
                    }
                }
            }

            m = next_m;
            n = next_n;
        }
        ss_cond_set(PLASMA_RANK, 0, 1);
        break;
    default:;
    }

    if (norm != PlasmaFrobeniusNorm) {
        step = 1;
        lrank = PLASMA_RANK;
        while ( (lrank%2 == 0) && (PLASMA_RANK+step < PLASMA_SIZE) ) {
            ss_cond_wait(PLASMA_RANK+step, 0, step);
            work[PLASMA_RANK] = max(work[PLASMA_RANK], work[PLASMA_RANK+step]);
            lrank = lrank >> 1;
            step  = step << 1;
            ss_cond_set(PLASMA_RANK, 0, step);
        }
        if (PLASMA_RANK > 0) {
            while( lrank != 0 ) {
                if (lrank%2 == 1) {
                    ss_cond_set(PLASMA_RANK, 0, step);
                    lrank = 0;
                } else {
                    lrank = lrank >> 1;
                    step  = step << 1;
                    ss_cond_set(PLASMA_RANK, 0, step);
                }
            }
        }

        if (PLASMA_RANK == 0)
            *result = work[0];
    }
    else {
        step = 1;
        lrank = PLASMA_RANK;
        while ( (lrank%2 == 0) && (PLASMA_RANK+step < PLASMA_SIZE) ) {
            float scale1, scale2;
            float sumsq1, sumsq2;

            ss_cond_wait(PLASMA_RANK+step, 0, step);

            scale1 = work[ 2 * PLASMA_RANK ];
            sumsq1 = work[ 2 * PLASMA_RANK + 1 ];
            scale2 = work[ 2 * (PLASMA_RANK+step) ];
            sumsq2 = work[ 2 * (PLASMA_RANK+step) + 1 ];

            if ( scale2 != 0. ){
                if( scale1 < scale2 ) {
                    work[2 * PLASMA_RANK+1] = sumsq2 + (sumsq1 * (( scale1 / scale2 ) * ( scale1 / scale2 )));
                    work[2 * PLASMA_RANK  ] = scale2;
                } else {
                    work[2 * PLASMA_RANK+1] = sumsq1 + (sumsq2 * (( scale2 / scale1 ) * ( scale2 / scale1 )));
                }
            }
            lrank = lrank >> 1;
            step  = step << 1;
            ss_cond_set(PLASMA_RANK, 0, step);
        }
        if (PLASMA_RANK > 0) {
            while( lrank != 0 ) {
                if (lrank%2 == 1) {
                    ss_cond_set(PLASMA_RANK, 0, step);
                    lrank = 0;
                } else {
                    lrank = lrank >> 1;
                    step  = step << 1;
                    ss_cond_set(PLASMA_RANK, 0, step);
                }
            }
        }

        if (PLASMA_RANK == 0)
            *result = work[0] * sqrt( work[1] );
    }
    ss_finalize();
}

/***************************************************************************//**
 *
 **/
void plasma_pslansy_quark(PLASMA_enum norm, PLASMA_enum uplo, PLASMA_desc A, float *work, float *result,
                          PLASMA_sequence *sequence, PLASMA_request *request)
{
    plasma_context_t *plasma;
    Quark_Task_Flags task_flags = Quark_Task_Flags_Initializer;

    float* lwork;
    int X, X1, X2, Y, Y1;
    int ldam;
    int m, n, k;
    int szeW, pos;

    plasma = plasma_context_self();
    if (sequence->status != PLASMA_SUCCESS)
        return;
    QUARK_Task_Flag_Set(&task_flags, TASK_SEQUENCE, (intptr_t)sequence->quark_sequence);
    plasma_profile_by_function( &task_flags, LANSY );

    *result = 0.0;
    switch ( norm ) {
    /*
     *  PlasmaMaxNorm
     */
    case PlasmaMaxNorm:
        szeW = (A.mt*(A.mt+1))/2;
        pos = 0;
        lwork = (float *)plasma_shared_alloc(plasma, szeW, PlasmaRealDouble);
        memset(lwork, 0, szeW*sizeof(float));
        for(m = 0; m < A.mt; m++) {
            X1 = m == 0      ?  A.i       %A.mb   : 0;
            X2 = m == A.mt-1 ? (A.i+A.m-1)%A.mb+1 : A.mb;
            X = X2 - X1;
            ldam = BLKLDD(A, m);
            QUARK_CORE_slansy_f1(
                plasma->quark, &task_flags,
                PlasmaMaxNorm, uplo, X,
                A(m, m, X1, X1, ldam), ldam, ldam*X,
                0, lwork + pos,
                result, 1);
            pos++;
            /*
             *  PlasmaLower
             */
            if (uplo == PlasmaLower) {
                for(n = 0; n < m; n++) {
                    Y1 = n == 0 ? A.j%A.nb : 0;
                    Y  = A.nb - Y1;
                    QUARK_CORE_slange_f1(
                        plasma->quark, &task_flags,
                        PlasmaMaxNorm, X, Y,
                        A(m, n, X1, Y1, ldam), ldam, ldam*Y,
                        0, lwork + pos,
                        result, 1);
                    pos++;
                }
            }
            /*
             *  PlasmaUpper
             */
            else {
                for(n = m+1; n < A.mt; n++) {
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
     *  PlasmaOneNorm / PlasmaInfNorm
     */
    case PlasmaOneNorm:
    case PlasmaInfNorm:
        lwork = (float *)plasma_shared_alloc(plasma, A.m, PlasmaRealDouble);
        memset(lwork, 0, A.m*sizeof(float));
        for(m = 0; m < A.mt; m++) {
            X1 = m == 0      ?  A.i       %A.mb   : 0;
            X2 = m == A.mt-1 ? (A.i+A.m-1)%A.mb+1 : A.mb;
            X = X2 - X1;
            ldam = BLKLDD(A, m);
            QUARK_CORE_sasum_f1(
                plasma->quark, &task_flags,
                PlasmaRowwise, uplo, X, X,
                A(m, m, X1, X1, ldam), ldam, ldam*X,
                lwork + m*A.mb, A.mb,
                result, 1);
            /*
             *  PlasmaLower
             */
            if (uplo == PlasmaLower) {
                for(n = 0; n < m; n++) {
                    Y1 = n == 0 ? A.j%A.nb : 0;
                    Y = A.nb - Y1;
                    QUARK_CORE_sasum_f1(
                        plasma->quark, &task_flags,
                        PlasmaRowwise, PlasmaUpperLower, X, Y,
                        A(m, n, X1, Y1, ldam), ldam, ldam*Y,
                        lwork + m*A.mb, A.mb,
                        result, 1);

                    QUARK_CORE_sasum_f1(
                        plasma->quark, &task_flags,
                        PlasmaColumnwise, PlasmaUpperLower, X, Y,
                        A(m, n, X1, Y1, ldam), ldam, ldam*Y,
                        lwork + n*A.mb, A.mb,
                        result, 1);
                }
            }
            /*
             *  PlasmaUpper
             */
            else {
                for(n = m+1; n < A.mt; n++) {
                    Y = n == A.nt-1 ? (A.j+A.n-1)%A.nb+1 : A.nb;
                    QUARK_CORE_sasum_f1(
                        plasma->quark, &task_flags,
                        PlasmaRowwise, PlasmaUpperLower, X, Y,
                        A(m, n, X1, 0, ldam), ldam, ldam*Y,
                        lwork + m*A.mb, A.mb,
                        result, 1);

                    QUARK_CORE_sasum_f1(
                        plasma->quark, &task_flags,
                        PlasmaColumnwise, PlasmaUpperLower, X, Y,
                        A(m, n, X1, 0, ldam), ldam, ldam*Y,
                        lwork + n*A.mb, A.mb,
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
        for(m = 0; m < A.mt; m++) {
            X1 = m == 0      ?  A.i       %A.mb   : 0;
            X2 = m == A.mt-1 ? (A.i+A.m-1)%A.mb+1 : A.mb;
            X = X2 - X1;
            ldam = BLKLDD(A, m);

            QUARK_CORE_ssyssq_f1(
                plasma->quark, &task_flags,
                uplo, X,
                A(m, m, X1, X1, ldam), ldam,
                lwork + 2*k,
                lwork + 2*k + 1,
                result, 1, INOUT | GATHERV );
            k++;
            k = k % PLASMA_SIZE;

            /*
             *  PlasmaLower
             */
            if (uplo == PlasmaLower) {
                for(n = 0; n < m; n++) {
                    Y1 = n == 0 ? A.j%A.nb : 0;
                    Y = A.nb - Y1;

                    QUARK_CORE_sgessq_f1(
                        plasma->quark, &task_flags,
                        X, Y,
                        A(m, n, X1, Y1, ldam), ldam,
                        lwork + 2*k,
                        lwork + 2*k + 1,
                        result, 1, INOUT | GATHERV );

                    QUARK_CORE_sgessq_f1(
                        plasma->quark, &task_flags,
                        X, Y,
                        A(m, n, X1, Y1, ldam), ldam,
                        lwork + 2*k,
                        lwork + 2*k + 1,
                        result, 1, INOUT | GATHERV );
                    k++;
                    k = k % PLASMA_SIZE;

                }
            }
            /*
             *  PlasmaUpper
             */
            else {
                for(n = m+1; n < A.mt; n++) {
                    Y = n == A.nt-1 ? (A.j+A.n-1)%A.nb+1 : A.nb;

                    QUARK_CORE_sgessq_f1(
                        plasma->quark, &task_flags,
                        X, Y,
                        A(m, n, X1, 0, ldam), ldam,
                        lwork + 2*k,
                        lwork + 2*k + 1,
                        result, 1, INOUT | GATHERV );

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

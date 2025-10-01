/**
 *
 * @file core_cgetrf_reclap.c
 *
 *  PLASMA core_blas kernel
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Hatem Ltaief
 * @author Mathieu Faverge
 * @author Piotr Luszczek
 * @date 2009-11-15
 *
 * @generated c Fri Apr  1 11:02:30 2016
 *
 **/

#include <math.h>
#include <cblas.h>
#ifdef PLASMA_WITH_MKL
#include <mkl_lapacke.h>
#else
#include <lapacke.h>
#endif
#include "common.h"

struct CORE_cgetrf_data_s {
    volatile PLASMA_Complex32_t *CORE_camax;
    volatile int                *CORE_cstep;
};

static inline void
CORE_cgetrf_reclap_update(CORE_cgetrf_data_t *data,
                          int M, int column, int n1, int n2,
                          PLASMA_Complex32_t *A, int LDA, int *IPIV,
                          int thidx, int thcnt);
static inline void
CORE_cgetrf_reclap_rec(CORE_cgetrf_data_t *data,
                       int M, int N,
                       PLASMA_Complex32_t *A, int LDA,
                       int *IPIV, int *info,
                       int thidx, int thcnt, int column);

/***************************************************************************//**
 *
 * @ingroup CORE_PLASMA_Complex32_t
 *
 *  CORE_cgetrf_reclap computes a LU factorization of a general M-by-N
 *  matrix A stored in CCRB layout using partial pivoting with row
 *  interchanges.
 *
 *  The factorization has the form
 *
 *    A = P * L * U
 *
 *  where P is a permutation matrix, L is lower triangular with unit
 *  diagonal elements (lower trapezoidal if m > n), and U is upper
 *  triangular (upper trapezoidal if m < n).
 *
 *  This is the recursive version of the algorithm applied on column
 *  major layout.
 *
 *  WARNINGS:
 *     - The function CORE_cgetrf_reclap_init has to be called prior
 *     to any call to this function.
 *     - You cannot call this kernel on different matrices at the same
 *     time.
 *     - The matrix A cannot be more than one tile wide.
 *     - The number of threads calling this function has to be excatly
 *     the number defined by info[2] with each one of them a different
 *     index between 0 included and info[2] excluded.
 *
 *******************************************************************************
 *
 *  @param[in] data
 *          Common data structure to all threads initialized by
 *          CORE_cgetrf_reclap_init() that contains information for thread
 *          synchronization and maximum reductions. All threads working on a
 *          given matrix A must provide the same data structure.
 *
 *  @param[in] M
 *          The number of rows of the matrix A.  M >= 0.
 *
 *  @param[in] N
 *          The number of columns of the matrix A.  N >= 0.
 *
 *  @param[in,out] A
 *          On entry, the M-by-N matrix to be factorized.
 *          On exit, the factors L and U from the factorization
 *          A = P*L*U; the unit diagonal elements of L are not stored.
 *
 *  @param[in] LDA
 *          The leading dimension of the array A.  LDA >= max(1,M).
 *
 *  @param[out] IPIV
 *          The pivot indices; for 1 <= i <= min(M,N), row i of the
 *          matrix was interchanged with row IPIV(i).
 *          1 <= IPIV[i] <= M.
 *
 *  @param[in,out] info
 *          Array of 3 integers
 *          - info[0], see returned value
 *          - info[1], is the thread index 0 <= info[0] < info[2]
 *          - info[2], on entry is the number of threads trying to
 *                     participate to the factorization,
 *                     on exit is the real number of threads used to
 *                     perform the factorization.
 *    Info[2] threads, and exactly info[2], have to call this function
 *    to avoid dead lock.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *          \retval -k, the k-th argument had an illegal value
 *          \retval k if U(k,k) is exactly zero. The factorization
 *                  has been completed, but the factor U is exactly
 *                  singular, and division by zero will occur if it is used
 *                  to solve a system of equations.
 *
 */
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_cgetrf_reclap = PCORE_cgetrf_reclap
#define CORE_cgetrf_reclap PCORE_cgetrf_reclap
#endif
int CORE_cgetrf_reclap(CORE_cgetrf_data_t *data,
                       int M, int N,
                       PLASMA_Complex32_t *A, int LDA,
                       int *IPIV, int *info)
{
    int thidx = info[1];
    int thcnt = min( info[2], M / N );
    int minMN = min(M, N);

    info[0] = 0;
    info[2] = thcnt;

    if( M < 0 ) {
        coreblas_error(1, "illegal value of M");
        return -1;
    }
    if( N < 0 ) {
        coreblas_error(2, "illegal value of N");
        return -2;
    }
    if( LDA < max(1, M) ) {
        coreblas_error(5, "illegal value of LDA");
        return -5;
    }

    /*
     * Quick return
     */
    if ( (M == 0) || (N == 0) || (thidx >= thcnt) ){
      return PLASMA_SUCCESS;
    }

    *info = 0;
    CORE_cgetrf_reclap_rec( data, M, minMN, A, LDA, IPIV, info,
                            thidx, thcnt, 0 );

    if ( N > minMN ) {
        CORE_cgetrf_reclap_update(data, M, 0, minMN, N-minMN,
                                  A, LDA, IPIV,
                                  thidx, thcnt);
    }

    return info[0];
}

/*******************************************************************
 *   Additional routines
 */
static float sfmin = 0.;

CORE_cgetrf_data_t *
CORE_cgetrf_reclap_init(int nbthrd)
{
    int i;
    CORE_cgetrf_data_t *data;

    data = (CORE_cgetrf_data_t*)malloc( nbthrd * (sizeof(PLASMA_Complex32_t)+sizeof(int))
                                        + 2 * sizeof(void*) );
    data->CORE_camax = (PLASMA_Complex32_t*)((char*)data + 2 * sizeof(void*));
    data->CORE_cstep = (int*)((char*)data + 2 * sizeof(void*) + nbthrd * sizeof(PLASMA_Complex32_t));

    for (i = 0; i < nbthrd; ++i) {
        data->CORE_camax[i] = 0.;
        data->CORE_cstep[i] = -1;
    }
    if ( sfmin == 0. ) {
        sfmin = LAPACKE_slamch_work('S');
    }
    return data;
}

static inline void
psplit(int n, int pidx, int pcnt, int *poff_p, int *psiz_p)
{
    int q = n / pcnt, r = n % pcnt;

    if (pidx < r) {
        q++;
        *psiz_p = q;
        *poff_p = pidx * q;
    } else {
        *psiz_p = q;
        *poff_p = r * (q + 1) + (pidx - r) * q;
    }
}

static inline void
CORE_camax1_thread(CORE_cgetrf_data_t *data,
                   PLASMA_Complex32_t localamx,
                   int thidx, int thcnt, int *thwinner,
                   PLASMA_Complex32_t *globalamx,
                   int pividx, int *ipiv)
{
    volatile PLASMA_Complex32_t *CORE_camax = data->CORE_camax;
    volatile int                *CORE_cstep = data->CORE_cstep;

    if (thidx == 0) {
        int i, j = 0;
        PLASMA_Complex32_t curval = localamx, tmp;
        float curamx = cabsf(localamx);

        /* make sure everybody filled in their value */
        for (i = 1; i < thcnt; ++i) {
            while (CORE_cstep[i] == -1) { /* wait for thread i to store its value */
            }
        }

        /* better not fuse the loop above and below to make sure data is sync'd */

        for (i = 1; i < thcnt; ++i) {
            tmp = CORE_camax[i];
            if (cabsf(tmp) > curamx) {
                curamx = cabsf(tmp);
                curval = tmp;
                j = i;
            }
        }

        if (0 == j)
            ipiv[0] = pividx;

        /* make sure everybody knows the amax value */
        for (i = 1; i < thcnt; ++i)
            CORE_camax[ i ] = curval;

        CORE_cstep[0] = -j - 2; /* set the index of the winning thread */

        *thwinner = j;
        *globalamx = curval;

        for (i = 1; i < thcnt; ++i)
            CORE_cstep[i] = -3;

        /* make sure everybody read the max value */
        for (i = 1; i < thcnt; ++i) {
            while (CORE_cstep[i] != -1) {
            }
        }

        CORE_cstep[0] = -1;
    } else {
        CORE_camax[thidx] = localamx;
        CORE_cstep[thidx] = -2;  /* announce to thread 0 that local amax was stored */
        while (CORE_cstep[0] == -1) { /* wait for thread 0 to finish calculating the global amax */
        }
        while (CORE_cstep[thidx] != -3) { /* wait for thread 0 to store amax */
        }
        *thwinner  = -CORE_cstep[0] - 2;
        *globalamx = CORE_camax[thidx];     /* read the amax from the location adjacent to the one in the above loop */
        CORE_cstep[thidx] = -1;  /* signal thread 0 that this thread is done reading */

        if (thidx == *thwinner)
            ipiv[0] = pividx;

        while (CORE_cstep[0] != -1) { /* wait for thread 0 to finish */
        }
    }
}

static inline void
CORE_cbarrier_thread(CORE_cgetrf_data_t *data,
                     int thidx, int thcnt)
{
    int idum1, idum2;
    PLASMA_Complex32_t ddum2 = 0.;
    /* it's probably faster to implement a dedicated barrier */
    CORE_camax1_thread( data, 1.0, thidx, thcnt, &idum1, &ddum2, 0, &idum2 );
}

static inline void
CORE_claswap1(int ncol, PLASMA_Complex32_t *a, int lda,
              int idxStart, int idxMax, const int *piv)
{
    int i, j;
    PLASMA_Complex32_t tmp;

    for (j = 0; j < ncol; ++j) {
        for (i = idxStart; i < idxMax; ++i) {
            tmp = a[j*lda + piv[i] - 1];
            a[j*lda + piv[i] - 1] = a[i + j*lda];
            a[i + j*lda] = tmp;
        }
    }
}

static inline void
CORE_cgetrf_reclap_update(CORE_cgetrf_data_t *data,
                          int M, int column, int n1, int n2,
                          PLASMA_Complex32_t *A, int LDA, int *IPIV,
                          int thidx, int thcnt)
{
    static PLASMA_Complex32_t posone =  1.0;
    static PLASMA_Complex32_t negone = -1.0;
    PLASMA_Complex32_t *Atop  = A    + column*LDA;
    PLASMA_Complex32_t *Atop2 = Atop + n1    *LDA;
    int coff, ccnt, lm, loff;

    CORE_cbarrier_thread( data, thidx, thcnt );

    psplit( n2, thidx, thcnt, &coff, &ccnt );

    if (ccnt > 0) {
        CORE_claswap1( ccnt, Atop2 + coff*LDA, LDA, column, n1 + column, IPIV ); /* swap to the right */

        cblas_ctrsm( CblasColMajor, CblasLeft, CblasLower, CblasNoTrans, CblasUnit,
                     n1, ccnt, CBLAS_SADDR(posone), Atop + column, LDA, Atop2 + coff*LDA + column, LDA );
    }

    /* __sync_synchronize(); */ /* hopefully we will not need memory fences */

    /* need to wait for pivoting and triangular solve to finish */
    CORE_cbarrier_thread( data, thidx, thcnt );

    psplit( M, thidx, thcnt, &loff, &lm );
    if (thidx == 0) {
        loff = column + n1;
        lm  -= column + n1;
    };

    cblas_cgemm( CblasColMajor, CblasNoTrans, CblasNoTrans, lm, n2, n1,
                 CBLAS_SADDR(negone), Atop+loff, LDA, Atop2 + column, LDA, CBLAS_SADDR(posone), Atop2+loff, LDA );
}

static void
CORE_cgetrf_reclap_rec(CORE_cgetrf_data_t *data, int M, int N,
                       PLASMA_Complex32_t *A, int LDA,
                       int *IPIV, int *info,
                       int thidx, int thcnt, int column)
{
    int jp, n1, n2, lm, loff;
    PLASMA_Complex32_t tmp1, tmp2, tmp3;
    PLASMA_Complex32_t *Atop = A + column*LDA;

    /* Assumption: N = min( M, N ); */
    if (N > 1) {
        int coff, ccnt;

        n1 = N / 2;
        n2 = N - n1;

        CORE_cgetrf_reclap_rec( data, M, n1, A, LDA, IPIV, info,
                                thidx, thcnt, column );
        if ( *info != 0 )
            return;

        CORE_cgetrf_reclap_update(data, M, column, n1, n2,
                                  A, LDA, IPIV,
                                  thidx, thcnt);

        CORE_cgetrf_reclap_rec( data, M, n2, A, LDA, IPIV, info,
                                thidx, thcnt, column + n1 );
        if ( *info != 0 )
            return;

        psplit( n1, thidx, thcnt, &coff, &ccnt );

        if (ccnt > 0) {
            CORE_claswap1( ccnt, Atop+coff*LDA, LDA, n1 + column, N + column, IPIV ); /* swap to the left */
        }

    } else {
        int thrd;

        CORE_cbarrier_thread( data, thidx, thcnt );

        psplit( M, thidx, thcnt, &loff, &lm );

        if (thidx == 0) {
            loff = column;
            lm -= column;
        }

        tmp2 = Atop[column]; /* all threads read the pivot element in case they need it */

        jp = cblas_icamax( lm, Atop + loff, 1 );
        tmp1 = Atop[loff + jp];

        CORE_camax1_thread( data, tmp1, thidx, thcnt, &thrd,
                            &tmp3, loff + jp + 1, IPIV + column );

        Atop[column] = tmp3; /* all threads set the pivot element: no need for synchronization */

        if ( tmp3 != 0.0 ) {
            if ( cabsf(tmp3) >= sfmin ) {
                PLASMA_Complex32_t tmp = (PLASMA_Complex32_t)1.0 / tmp3;
                n1 = (thidx == 0) ? 1 : 0;
                cblas_cscal( lm - n1, CBLAS_SADDR(tmp), Atop + loff + n1, 1 );
            } else {
                int i;
                PLASMA_Complex32_t *Atop2;
                n1 = (thidx == 0) ? 1 : 0;
                Atop2 = Atop + loff + n1;

                for( i=0; i < lm-n1; i++, Atop2++)
                    *Atop2 = *Atop2 / tmp3;
            }

            if (thrd == thidx) { /* the thread that owns the best pivot */
              if (loff + jp != column) /* if there is a need to exchange the pivot */
                Atop[loff + jp] = tmp2 / tmp3;
            }

        } else {
            *info = column + 1;
            return;
        }

        CORE_cbarrier_thread( data, thidx, thcnt );
    }
}

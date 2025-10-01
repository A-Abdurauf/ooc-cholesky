/**
 * @file core_dlaed2_computeK.c
 *
 *  PLASMA computational routines
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
#include <math.h>
#ifdef PLASMA_WITH_MKL
#include <mkl_lapacke.h>
#else
#include <lapacke.h>
#endif
#include "common.h"

static int ione = 1;

/***************************************************************************//**
 *
 * @ingroup CORE_double
 *
 *  CORE_dlaed2_computeK - Computes the number of deflated eigenvalues
 *  of a symmetric tridiagonal matrix.
 *  A eigenvalue is deflated if it's close enough to another eigenvalue
 *  or if the vector Z contains a tiny entry.
 *  This operation is useful when two subproblems, which were previously solved,
 *  are merge into a larger problem.
 *
 *******************************************************************************
 *
 * @param[out] Kptr
 *          On exit, Kptr stores the number of non-deflated eigenvalues
 *
 * @param[in] n
 *          n specifies the dimension of the symmetric tridiagonal matrix
 *
 * @param[in] n1
 *          n1 specifies the location of the last eigenvalue of the first subproblem
 *          min(1, n) <= n1 <= n/2
 *
 * @param[in,out] betaptr
 *          betaptr[0] specifies the rank-1 approximation that was used for splitting
 *          the problem into two subproblems.
 *          On exit, it contains the updated beta for LAED3.
 *
 * @param[in,out] D
 *          On entry, D contains the eigenvalues of the two submatrices to be merged.
 *          On exit, D contains the updated eigenvalues sorted into increasing order.
 *
 * @param[in,out] Q
 *          On entry, Q contains the eigenvectors of each subproblem.
 *          On exit, Q contains the updated eigenvectors, which can be modified through
 *          givens rotation (with close eigenvalues).
 *
 * @param[in] LDQ
 *          LDQ specifies the leading direction of Q
 *
 * @param[out] Z
 *          Z will contain the updating vectors, composed of the last row of the first
 *          subproblem and the first row of the second subproblem.
 *
 * @param[out] DLAMBDA
 *          DLAMBDA will contain a copy of the first K eigenvalues and then an update of
 *          those eigenvalues for usage in LAED3. In LAPACK, the call of dlamc3 is done at
 *          the beginning of LAED3.
 *
 * @param[out] W
 *          W will contain the fisrt K values of the final Z vector, for future use in LAED3.
 *
 * @param[out] INDX
 *          The permutation used to sort the contents of DLAMBDA into ascending order
 *
 * @param[out] INDXC
 *          The permutation used to arrange the columns of the deflated Q matrix
 *          1- non-zero elements only at and above n1
 *          2- non-zero elements only below n1
 *          3- dense
 *
 * @param[out] INDXP
 *          INDXP is the permutation used to place deflated eigenvalues at the end of the array
 *          INDXP(0:K-1) points to the non-deflated D values
 *          INDXP(K-1:N-1) points to the deflated D values
 *
 * @param[in,out] INDXQ
 *          The permutation which separetely sorts each subprobem in D into ascending order.
 *          In the second half of the permutation, elements will have n1 added to their value.
 *
 * @param[out] COLTYP
 *          COLTYP indicate the type of the Q matrix columns during the process.
 *          1- non-zero in the upper half only
 *          2- dense
 *          3- non-zero in the lower half only
 *          4- deflated
 *
 *******************************************************************************/

#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_dlaed2_computeK = PCORE_dlaed2_computeK
#define CORE_dlaed2_computeK PCORE_dlaed2_computeK
#endif
void CORE_dlaed2_computeK(int *Kptr, int n, int n1,
                          double *betaptr, double *D,
                          double *Q, int LDQ,
                          double *Z, double *DLAMBDA, double *W,
                          int *INDX, int *INDXC, int *INDXP, int *INDXQ,
                          int *COLTYP)
{
    double maxz, tol, t, tau, C, S;
    double beta = betaptr[0];
    int K, K2, i, j, ct, ni;
    int n2 = n - n1;
    int pi = -1;
    int go = 1;

    int pos[4], ctot[4];

    /*
     * Form the z-vector which consists of the last row of Q_1 and the
     * first row of Q_2.
     */
    cblas_dcopy(n1,   Q+n1-1,      LDQ, Z,    1);
    cblas_dcopy(n-n1, Q+LDQ*n1+n1, LDQ, Z+n1, 1);

    /*
     * Normalize z so that norm(z) = 1.  Since z is the concatenation of
     * two normalized vectors, norm2(z) = sqrt(2).
     */
    {
        double scale = 1. / sqrt(2.);
        cblas_dscal(n1, scale, Z, 1);
        if (beta < 0.0) {
            cblas_dscal(n2, -scale, Z+n1, 1);
        } else {
            cblas_dscal(n2,  scale, Z+n1, 1);
        }
    }

    /* beta = abs( norm(z)**2 * beta ) : cancel Z normalization */
    beta = fabs(2. * beta);

    /*
     * Sort eigenvalues into increasing order using INDXQ
     */
    /* Increment the second part of INDXQ */
    for (i=n1; i<n; i++){
        INDXQ[i] += n1;
    }
    /* Copy D in Dlambda modulo permutation */
    for (i=0; i<n; i++){
        DLAMBDA[i] = D[INDXQ[i]];
    }
    /* Compute the merge of the two sorted subsets */
    PLASMA_FCALL(dlamrg, DLARMG)(&n1, &n2, DLAMBDA, &ione, &ione, INDXC);
    /* Compute the final indx array */
    {
        /*
         * INDXC has indices in Fortran mode, so we use a temporary pointer to
         * previous element
         */
        int *idxq = INDXQ-1;
        for (i=0; i<n; i++){
            INDX[i] = idxq[INDXC[i]];
        }
    }

    /* Calculate the allowable deflation tolerance */
    {
        double maxd;
        double eps = LAPACKE_dlamch_work( 'e' );
        int imax = cblas_idamax( n, Z, 1 );
        int jmax = cblas_idamax( n, D, 1 );

        maxz = fabs( Z[imax] );
        maxd = fabs( D[jmax] );
        tol = 8. * eps * (( maxd > maxz ) ? maxd : maxz );
    }

    memset( ctot, 0, 4 * sizeof(int) );
    memset( pos,  0, 4 * sizeof(int) );

    /* If beta is small enough, just reorganize Q */
    /* Contrary to Lapack dlaed2, eigen vectors are later sorted and copy back to Q */
    if ( (beta*maxz) <= tol ) {
        K = 0;
        /* Change the number of type 3 */
        ctot[3] = n;
        for (i=0; i<n; i++){
            DLAMBDA[i] = D[ INDX[i] ];
        }
        cblas_dcopy(n, DLAMBDA, 1, D, 1);
        goto end;
    }

    /* Otherwise search multiple eigenvalues */
    /* 1 for T1, 2 for mixte, 3 for T2 and 4 for deflated */
    for (i=0; i<n1; i++){
        COLTYP[i] = 1;
    }
    for (i=n1; i<n; i++){
        COLTYP[i] = 3;
    }

    K = 0;
    K2 = n;  /* index of last deflated eigenvalue */

    /* ni is the permuted index of i and pi the next value */
    i = 0;
    while (i<n && go == 1){
        ni = INDX[i];

        if (beta*fabs(Z[ni]) <= tol){
            K2--;
            COLTYP[ni] = 4;
            INDXP[K2] = ni;
            i++;
        }

        else{
            pi = ni;
            go = 0;
        }
    }

    while (i < n-1){
        i++;
        ni = INDX[i];

        if (beta*fabs(Z[ni]) <= tol){
            K2--;
            COLTYP[ni] = 4;
            INDXP[K2] = ni;
        }
        else{
            S = Z[pi];
            C = Z[ni];
            /*
             * Find sqrt(a**2+b**2) without overflow or destructive underflow.
             */
            tau = PLASMA_FCALL(dlapy2, DLAPY2)(&C, &S);
            t = D[ni] - D[pi];
            C = C/tau;
            S = -S/tau;

            /* Deflation is possible */
            if (fabs(t*C*S) <= tol){
                Z[ni] = tau;
                Z[pi] = 0.0;

                /* Mixte eigenvector */
                if (COLTYP[ni] != COLTYP[pi]){
                    COLTYP[ni] = 2;
                }
                COLTYP[pi] = 4;

                cblas_drot(n, Q+LDQ*pi, 1, Q+LDQ*ni, 1, C, S);

                t = D[pi]*C*C + D[ni]*S*S;
                D[ni] = D[pi]*S*S + D[ni]*C*C;
                D[pi] = t;
                K2--;

                j = 1;
                go = 1;
                while (j+K2 < n && go == 1){
                    if (D[pi] < D[INDXP[j+K2]]){
                        INDXP[K2+j-1] = INDXP[K2+j];
                        INDXP[K2+j] = pi;
                        j++;
                    }
                    else{
                        INDXP[K2+j-1] = pi;
                        go = 0;
                    }
                }
                if (go == 1){
                    INDXP[K2+j-1] = pi;
                }
                pi = ni;
            }

            /* No deflation */
            else{
                DLAMBDA[K] = D[pi];
                W[K] = Z[pi];
                INDXP[K] = pi;
                K++;
                pi = ni;
            }
        }
    }

    DLAMBDA[K] = D[pi];
    W[K] = Z[pi];
    INDXP[K] = pi;
    K++;

    /* For accuracy before solving the secular equation */
    for (i=0; i<K; i++){
        DLAMBDA[i] = PLASMA_FCALL(dlamc3, DLACM3)(&DLAMBDA[i], &DLAMBDA[i]) - DLAMBDA[i];
    }

    /* Different kind of eigenvalues: from T1, mixte, from T2 and deflated */
    for (i=0; i<4; i++){
        ctot[i] = 0;
    }
    for (i=0; i<n; i++){
        ct = COLTYP[i];
        ctot[ct-1]++;
    }

    /* Positions of eigenvectors in matrix */
    pos[0] = 0;
    pos[1] = ctot[0];
    pos[2] = pos[1] + ctot[1];
    pos[3] = pos[2] + ctot[2];
    K = n - ctot[3];

    /* Fill out INDXC for Q2 columns */
    {
        int is;
        for (i=0; i<K; i++){
            is = INDXP[i];
            ct = COLTYP[is];
            INDX[pos[ct-1]] = is;
            INDXC[pos[ct-1]] = i;
            pos[ct-1]++;
        }
        for (i=K; i<n; i++){
            is = INDXP[i];
            ct = COLTYP[is];
            INDX[pos[ct-1]] = is;
            INDXC[pos[ct-1]] = i;
            pos[ct-1]++;
        }
    }

    for (i=0; i<n; i++){
        Z[i] = D[ INDX[i] ];
    }

    cblas_dcopy(n-K, Z+K, 1, D+K, 1);

    /*
     *     Modify values DLAMDA(i) to make sure all DLAMDA(i)-DLAMDA(j) can
     *     be computed with high relative accuracy (barring over/underflow).
     *     This is a problem on machines without a guard digit in
     *     add/subtract (Cray XMP, Cray YMP, Cray C 90 and Cray 2).
     *     The following code replaces DLAMDA(I) by 2*DLAMDA(I)-DLAMDA(I),
     *     which on any of these machines zeros out the bottommost
     *     bit of DLAMDA(I) if it is 1; this makes the subsequent
     *     subtractions DLAMDA(I)-DLAMDA(J) unproblematic when cancellation
     *     occurs. On binary machines with a guard digit (almost all
     *     machines) it does not change DLAMDA(I) at all. On hexadecimal
     *     and decimal machines with a guard digit, it slightly
     *     changes the bottommost bits of DLAMDA(I). It does not account
     *     for hexadecimal or decimal machines without guard digits
     *     (we know of none). We use a subroutine call to compute
     *     2*DLAMBDA(I) to prevent optimizing compilers from eliminating
     *     this code.
     */
    /* For accuracy. Previously was in DLAED3 */
    for (i=0; i<K; i++){
        DLAMBDA[i] = PLASMA_FCALL(dlamc3, DLAMC3)(&DLAMBDA[i], &DLAMBDA[i]) - DLAMBDA[i];
    }

 end:

    /* Save ctot */
    for (i=0; i<4; i++){
        COLTYP[i] = ctot[i];
    }

    *Kptr = K;
    *betaptr = beta;
}

/**
 *
 * @file example_strsm.c
 *
 *  PLASMA testing routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @brief Example for solving a system of linear equations using 
 *        Cholesky factorization and PLASMA_strsm routine
 *
 * @version 2.8.0
 * @author Bilel Hadri
 * @date 2010-11-15
 * @generated s Fri Apr  1 11:03:08 2016
 *
 **/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <plasma.h>
#include <cblas.h>
#ifdef PLASMA_WITH_MKL
#include <mkl_lapacke.h>
#else
#include <lapacke.h>
#endif
#include <core_blas.h>

int check_solution(int, int, float*, int, float*, float*, int);

int IONE=1;
int ISEED[4] = {0,0,0,1};   /* initial seed for slarnv() */

int main ()
{

    int cores = 2;
    int N     = 10 ;
    int LDA   = 10 ;
    int NRHS  = 5 ;
    int LDB   = 10 ;
    int info_solution;

    float *A1   = (float *)malloc(LDA*N*sizeof(float));
    float *A2   = (float *)malloc(LDA*N*sizeof(float));
    float *B1   = (float *)malloc(LDB*NRHS*sizeof(float));
    float *B2   = (float *)malloc(LDB*NRHS*sizeof(float));

    /* Check if unable to allocate memory */
    if ((!A1)||(!A2)||(!B1)||(!B2)){
        printf("Out of Memory \n ");
        return EXIT_SUCCESS;
    }

    /* Plasma Initialize */
    PLASMA_Init(cores);
    printf("-- PLASMA is initialized to run on %d cores. \n",cores);

    /* Initialize A1 and A2 for Symmetric Positive Matrix */
    PLASMA_splgsy( (float)N, N, A1, LDA, 51 );
    PLASMA_slacpy( PlasmaUpperLower, N, N, A1, LDA, A2, LDA );

    /* Initialize B1 and B2 */
    PLASMA_splrnt( N, NRHS, B1, LDB, 371 );
    PLASMA_slacpy( PlasmaUpperLower, N, NRHS, B1, LDB, B2, LDB );

    /* PLASMA routines */
    PLASMA_spotrf(PlasmaLower, N, A2, LDA);
    PLASMA_strsm(PlasmaLeft, PlasmaLower, PlasmaNoTrans, PlasmaNonUnit, 
                 N, NRHS, (float)1.0, A2, LDA, B2, LDB);
    PLASMA_strsm(PlasmaLeft, PlasmaLower, PlasmaTrans, PlasmaNonUnit,
                 N, NRHS, (float)1.0, A2, LDA, B2, LDB);

    /* Check the solution */
    info_solution = check_solution(N, NRHS, A1, LDA, B1, B2, LDB);

    if ( info_solution != 0 )
       printf("-- Error in STRSM example ! \n");
    else
       printf("-- Run of STRSM example successful ! \n");

    free(A1); free(A2); free(B1); free(B2);

    PLASMA_Finalize();

    return EXIT_SUCCESS;
}


/*------------------------------------------------------------------------
 *  Check the accuracy of the solution of the linear system
 */

int check_solution(int N, int NRHS, float *A1, int LDA, float *B1, float *B2, int LDB)
{
    int info_solution;
    float Rnorm, Anorm, Xnorm, Bnorm;
    float alpha, beta;
    float *work = (float *)malloc(N*sizeof(float));
    float eps;

    eps = LAPACKE_slamch_work('e');

    alpha = 1.0;
    beta  = -1.0;

    Xnorm = LAPACKE_slange_work(LAPACK_COL_MAJOR, lapack_const(PlasmaInfNorm), N, NRHS, B2, LDB, work);
    Anorm = LAPACKE_slange_work(LAPACK_COL_MAJOR, lapack_const(PlasmaInfNorm), N, N, A1, LDA, work);
    Bnorm = LAPACKE_slange_work(LAPACK_COL_MAJOR, lapack_const(PlasmaInfNorm), N, NRHS, B1, LDB, work);

    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, N, NRHS, N, (alpha), A1, LDA, B2, LDB, (beta), B1, LDB);
    Rnorm = LAPACKE_slange_work(LAPACK_COL_MAJOR, lapack_const(PlasmaInfNorm), N, NRHS, B1, LDB, work);

    printf("============\n");
    printf("Checking the Residual of the solution \n");
    printf("-- ||Ax-B||_oo/((||A||_oo||x||_oo+||B||_oo).N.eps) = %e \n",Rnorm/((Anorm*Xnorm+Bnorm)*N*eps));

    if (Rnorm/((Anorm*Xnorm+Bnorm)*N*eps) > 10.0){
        printf("-- The solution is suspicious ! \n");
        info_solution = 1;
     }
    else{
        printf("-- The solution is CORRECT ! \n");
        info_solution = 0;
    }

    free(work);

    return info_solution;
}

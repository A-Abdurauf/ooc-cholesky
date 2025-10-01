      PROGRAM EXAMPLE_DPOSV_F
*
*********************************************************************
*     PLASMA example routine (version 2.8.0)
*     Author: Bilel Hadri
*     Release Date: November, 15th 2010
*     PLASMA is a software package provided by Univ. of Tennessee,
*     Univ. of California Berkeley and Univ. of Colorado Denver.
*     @generated d Fri Apr  1 11:03:08 2016
*********************************************************************
*
      IMPLICIT NONE
*
      INCLUDE "plasmaf.h"
*
*     Purpose
*     =======
*
*     FORTRAN EXAMPLE FOR PLASMA_DPOSV
*     Example for solving a system of linear equations using Cholesky
*     factorization
*
*     =====================================================================
*
*     .. Parameters ..
      INTEGER           CORES, N, NRHS
      PARAMETER         ( CORES = 2 )
      PARAMETER         ( N = 355 )
      PARAMETER         ( NRHS = 5 )
      DOUBLE PRECISION        DONE
      PARAMETER         ( DONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER*8   descA, descB
      DOUBLE PRECISION  A1( N, N ), B1( N, NRHS )
      DOUBLE PRECISION  A2( N, N ), B2( N, NRHS )
      DOUBLE PRECISION  AT( N, N ), BT( N, NRHS )
      DOUBLE PRECISION  ALPHA
      DOUBLE PRECISION  RWORK ( N )
      INTEGER           INFO
      DOUBLE PRECISION  XNORM, ANORM, BNORM, RNORM, RESULT, EPS
      DOUBLE PRECISION  DLAMCH, DLANGE
*     ..
*     .. External Subroutines ..
      EXTERNAL          DLARNV, ZLAGHE, DLAMCH, DLANGE
      EXTERNAL          PLASMA_INIT, PLASMA_DPOSV, PLASMA_FINALIZE
      EXTERNAL          DGEMM
*
*     Initialize Plasma
*
      CALL PLASMA_INIT( CORES, INFO )
      WRITE(*,*) "-- PLASMA is initialized on", CORES, "cores."

      CALL PLASMA_DISABLE( PLASMA_AUTOTUNING, INFO )
      CALL PLASMA_SET( PLASMA_TILE_SIZE, 200, INFO )
      CALL PLASMA_SET( PLASMA_INNER_BLOCK_SIZE, 32, INFO )
*
*     Initialization of the matrix A1
*
      ALPHA = N * 1.
      CALL PLASMA_DPLGSY( ALPHA, N, A1, N, 51, INFO );
      A2(:,:) = A1(:,:)
*
*     Initialization of the RHS
*
      CALL PLASMA_DPLRNT( N, NRHS, B1, N, 371, INFO )
      B2(:,:) = B1(:,:)
*
*     Perform the Cholesky solve
*
      CALL PLASMA_DESC_CREATE( descA, AT, PlasmaRealDouble,
     $     200, 200, 40000, N, N, 0, 0, N, N, INFO );

      CALL PLASMA_DESC_CREATE( descB, BT, PlasmaRealDouble,
     $     200, 200, 40000, N, NRHS, 0, 0, N, NRHS, INFO );

      CALL PLASMA_LAPACK_TO_TILE( A2, N, descA, INFO )
      CALL PLASMA_LAPACK_TO_TILE( B2, N, descB, INFO )

      CALL PLASMA_DPOTRF_TILE( PlasmaUpper, descA, INFO )
      CALL PLASMA_DTRSM_TILE( PlasmaLeft, PlasmaUpper,
     $     PlasmaTrans, PlasmaNonUnit,
     $     DONE, descA, descB, INFO )
      CALL PLASMA_DTRSM_TILE( PlasmaLeft, PlasmaUpper,
     $     PlasmaNoTrans, PlasmaNonUnit,
     $     DONE, descA, descB, INFO )

      CALL PLASMA_TILE_TO_LAPACK( descB, B2, N, INFO )

      CALL PLASMA_DESC_DESTROY( descA, INFO )
      CALL PLASMA_DESC_DESTROY( descB, INFO )
*
*     Check the solution
*
      XNORM = DLANGE('I', N, NRHS, B2, N, RWORK)
      ANORM = DLANGE('I', N, N,    A1, N, RWORK)
      BNORM = DLANGE('I', N, NRHS, B1, N, RWORK)

      CALL DGEMM('No transpose','No transpose', N, NRHS, N, DONE,
     $     A1, N, B2, N, -DONE, B1, N)

      RNORM = DLANGE('I', N, NRHS, B1, N, RWORK)

      EPS= DLAMCH('Epsilon')
      RESULT = RNORM / ((ANORM * XNORM + BNORM) * N * EPS)

      WRITE(*,*) '============'
      WRITE(*,*) 'Checking the Residual of the solution '
      WRITE(*,*) '-- ||Ax-B||_oo/((||A||_oo||x||_oo+||B||_oo).N.eps)=',
     $     RESULT

      IF (( INFO .ne. 0 ) .OR. (RESULT > 60.0)) THEN
         WRITE(*,*) "-- Error in DPOSV example !"
      ELSE
         WRITE(*,*) "-- Run of DPOSV example successful !"
      ENDIF
*
*     Finalize Plasma
*
      CALL PLASMA_FINALIZE( INFO )
*
*     End of EXAMPLE_DPOSV.
*
      END PROGRAM EXAMPLE_DPOSV_F

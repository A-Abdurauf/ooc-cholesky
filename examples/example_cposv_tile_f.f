      PROGRAM EXAMPLE_CPOSV_F
*
*********************************************************************
*     PLASMA example routine (version 2.8.0)
*     Author: Bilel Hadri
*     Release Date: November, 15th 2010
*     PLASMA is a software package provided by Univ. of Tennessee,
*     Univ. of California Berkeley and Univ. of Colorado Denver.
*     @generated c Fri Apr  1 11:03:08 2016
*********************************************************************
*
      IMPLICIT NONE
*
      INCLUDE "plasmaf.h"
*
*     Purpose
*     =======
*
*     FORTRAN EXAMPLE FOR PLASMA_CPOSV
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
      COMPLEX*8        ZONE
      PARAMETER         ( ZONE = ( 1.0, 0.0 ) )
*     ..
*     .. Local Scalars ..
      INTEGER*8   descA, descB
      COMPLEX*8  A1( N, N ), B1( N, NRHS )
      COMPLEX*8  A2( N, N ), B2( N, NRHS )
      COMPLEX*8  AT( N, N ), BT( N, NRHS )
      REAL  ALPHA
      REAL  RWORK ( N )
      INTEGER           INFO
      REAL  XNORM, ANORM, BNORM, RNORM, RESULT, EPS
      REAL  SLAMCH, CLANGE
*     ..
*     .. External Subroutines ..
      EXTERNAL          DLARNV, ZLAGHE, SLAMCH, CLANGE
      EXTERNAL          PLASMA_INIT, PLASMA_CPOSV, PLASMA_FINALIZE
      EXTERNAL          CGEMM
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
      CALL PLASMA_CPLGHE( ALPHA, N, A1, N, 51, INFO );
      A2(:,:) = A1(:,:)
*
*     Initialization of the RHS
*
      CALL PLASMA_CPLRNT( N, NRHS, B1, N, 371, INFO )
      B2(:,:) = B1(:,:)
*
*     Perform the Cholesky solve
*
      CALL PLASMA_DESC_CREATE( descA, AT, PlasmaComplexFloat,
     $     200, 200, 40000, N, N, 0, 0, N, N, INFO );

      CALL PLASMA_DESC_CREATE( descB, BT, PlasmaComplexFloat,
     $     200, 200, 40000, N, NRHS, 0, 0, N, NRHS, INFO );

      CALL PLASMA_LAPACK_TO_TILE( A2, N, descA, INFO )
      CALL PLASMA_LAPACK_TO_TILE( B2, N, descB, INFO )

      CALL PLASMA_CPOTRF_TILE( PlasmaUpper, descA, INFO )
      CALL PLASMA_CTRSM_TILE( PlasmaLeft, PlasmaUpper,
     $     PlasmaConjTrans, PlasmaNonUnit,
     $     ZONE, descA, descB, INFO )
      CALL PLASMA_CTRSM_TILE( PlasmaLeft, PlasmaUpper,
     $     PlasmaNoTrans, PlasmaNonUnit,
     $     ZONE, descA, descB, INFO )

      CALL PLASMA_TILE_TO_LAPACK( descB, B2, N, INFO )

      CALL PLASMA_DESC_DESTROY( descA, INFO )
      CALL PLASMA_DESC_DESTROY( descB, INFO )
*
*     Check the solution
*
      XNORM = CLANGE('I', N, NRHS, B2, N, RWORK)
      ANORM = CLANGE('I', N, N,    A1, N, RWORK)
      BNORM = CLANGE('I', N, NRHS, B1, N, RWORK)

      CALL CGEMM('No transpose','No transpose', N, NRHS, N, ZONE,
     $     A1, N, B2, N, -ZONE, B1, N)

      RNORM = CLANGE('I', N, NRHS, B1, N, RWORK)

      EPS= SLAMCH('Epsilon')
      RESULT = RNORM / ((ANORM * XNORM + BNORM) * N * EPS)

      WRITE(*,*) '============'
      WRITE(*,*) 'Checking the Residual of the solution '
      WRITE(*,*) '-- ||Ax-B||_oo/((||A||_oo||x||_oo+||B||_oo).N.eps)=',
     $     RESULT

      IF (( INFO .ne. 0 ) .OR. (RESULT > 60.0)) THEN
         WRITE(*,*) "-- Error in CPOSV example !"
      ELSE
         WRITE(*,*) "-- Run of CPOSV example successful !"
      ENDIF
*
*     Finalize Plasma
*
      CALL PLASMA_FINALIZE( INFO )
*
*     End of EXAMPLE_CPOSV.
*
      END PROGRAM EXAMPLE_CPOSV_F

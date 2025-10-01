      PROGRAM EXAMPLE_ZPOSV_F
*
*********************************************************************
*     PLASMA example routine (version 2.8.0)
*     Author: Bilel Hadri
*     Release Date: November, 15th 2010
*     PLASMA is a software package provided by Univ. of Tennessee,
*     Univ. of California Berkeley and Univ. of Colorado Denver.
*     @precisions normal z -> c d s
*********************************************************************
*
      IMPLICIT NONE
*
      INCLUDE "plasmaf.h"
*
*     Purpose
*     =======
*
*     FORTRAN EXAMPLE FOR PLASMA_ZPOSV
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
      COMPLEX*16        ZONE
      PARAMETER         ( ZONE = ( 1.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      COMPLEX*16  A1( N, N ), B1( N, NRHS )
      COMPLEX*16  A2( N, N ), B2( N, NRHS )
      DOUBLE PRECISION  ALPHA
      DOUBLE PRECISION  RWORK ( N )
      INTEGER           INFO
      DOUBLE PRECISION  XNORM, ANORM, BNORM, RNORM, RESULT, EPS
      DOUBLE PRECISION  DLAMCH, ZLANGE
*     ..
*     .. External Subroutines ..
      EXTERNAL          DLARNV, ZLAGHE, DLAMCH, ZLANGE
      EXTERNAL          PLASMA_INIT, PLASMA_ZPOSV, PLASMA_FINALIZE
      EXTERNAL          ZGEMM
*
*     Initialize Plasma
*
      CALL PLASMA_INIT( CORES, INFO )
      WRITE(*,*) "-- PLASMA is initialized on", CORES, "cores."
*
*     Initialization of the matrix A1
*
      ALPHA = N * 1.
      CALL PLASMA_ZPLGHE( ALPHA, N, A1, N, 51, INFO );
      A2(:,:) = A1(:,:)
*
*     Initialization of the RHS
*
      CALL PLASMA_ZPLRNT( N, NRHS, B1, N, 371, INFO )
      B2(:,:)=B1(:,:)
*
*     Perform the Cholesky solve
*
      CALL PLASMA_ZPOSV( PlasmaUpper, N, NRHS, A2, N, B2, N, INFO )
*
*     Check the solution
*
      XNORM = ZLANGE('I',N, NRHS, B2, N, RWORK)
      ANORM = ZLANGE('I',N, N, A1, N, RWORK)
      BNORM = ZLANGE('I',N, NRHS, B1, N, RWORK)

      CALL ZGEMM('No transpose','No transpose', N, NRHS, N, ZONE,
     $     A1, N, B2, N, -ZONE, B1, N)

      RNORM = ZLANGE('I',N, NRHS, B1, N, RWORK)

      EPS= DLAMCH('Epsilon')
      RESULT = RNORM / ((ANORM * XNORM + BNORM) * N * EPS)

      WRITE(*,*) '============'
      WRITE(*,*) 'Checking the Residual of the solution '
      WRITE(*,*) '-- ||Ax-B||_oo/((||A||_oo||x||_oo+||B||_oo).N.eps)=',
     $     RESULT

      IF (( INFO .ne. 0 ) .OR. (RESULT > 60.0)) THEN
         WRITE(*,*) "-- Error in ZPOSV example !"
      ELSE
         WRITE(*,*) "-- Run of ZPOSV example successful !"
      ENDIF
*
*     Finalize Plasma
*
      CALL PLASMA_FINALIZE( INFO )
*
*     End of EXAMPLE_ZPOSV.
*
      END PROGRAM EXAMPLE_ZPOSV_F

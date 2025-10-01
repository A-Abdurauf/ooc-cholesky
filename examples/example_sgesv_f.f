        PROGRAM EXAMPLE_SGESV_F
*
*********************************************************************
*     PLASMA example routine (version 2.8.0)
*     Author: Bilel Hadri
*     Release Date: November, 15th 2010
*     PLASMA is a software package provided by Univ. of Tennessee,
*     Univ. of California Berkeley and Univ. of Colorado Denver.
*     @generated s Fri Apr  1 11:03:08 2016
*********************************************************************
*
        IMPLICIT NONE
*
        INCLUDE "plasmaf.h"
*
*  Purpose
*  =======
*
*   FORTRAN EXAMPLE FOR PLASMA_SGESV
*   Example for solving a system of linear equations using LU
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER           CORES, N, NRHS
      PARAMETER         ( CORES = 2 )
      PARAMETER         ( N = 10 )
      PARAMETER         ( NRHS = 5 )
      REAL        DONE
      PARAMETER         ( DONE = 1.0 )
*     ..
*     .. Local Scalars ..
      REAL        A1( N, N ), B1( N, NRHS )
      REAL        A2( N, N ), B2( N, NRHS )
      REAL  RWORK( N )
      INTEGER           IPIV( N )
      INTEGER           I, INFO
      INTEGER           ISEED( 4 )
      REAL  XNORM, ANORM, BNORM, RNORM, RESULT, EPS
      REAL  SLAMCH, SLANGE
*     ..
*     .. External Subroutines ..
      EXTERNAL          SLARNV, SLAMCH, SLANGE
      EXTERNAL          PLASMA_INIT
      EXTERNAL          PLASMA_SGESV, PLASMA_FINALIZE
      EXTERNAL          PLASMA_DEALLOC_HANDLE
      EXTERNAL          SGEMM
*     ..
*     .. Executable Statements ..
*
      DO  I = 1, 3
          ISEED( I ) = 0
      ENDDO
      ISEED( 4 ) = 1
*
*     Initialize Plasma
*
      CALL PLASMA_INIT( CORES, INFO )
      WRITE(*,*) "-- PLASMA is initialized on", CORES, "cores."
*
*     Initialization of the matrix A1
*
      CALL SLARNV( 1, ISEED, N*N, A1 )
      A2(:,:)=A1(:,:)
*
*     Initialization of the RHS
*
      CALL SLARNV( 1, ISEED, N*NRHS, B1 )
      B2(:,:)=B1(:,:)
*
*     PLASMA SGESV
*
      CALL PLASMA_SGESV( N, NRHS, A2, N, IPIV, B2, N, INFO )
*
*     Check the solution
*
      ANORM = SLANGE('I', N, N,    A1, N, RWORK)
      XNORM = SLANGE('I', N, NRHS, B2, N, RWORK)
      BNORM = SLANGE('I', N, NRHS, B1, N, RWORK)

      CALL SGEMM('No transpose', 'No transpose', N, NRHS, N, DONE,
     $     A1, N, B2, N, -DONE, B1, N)

      RNORM = SLANGE('I', N, NRHS, B1, N, RWORK)

      EPS = SLAMCH('Epsilon')
      RESULT = RNORM / ((ANORM * XNORM + BNORM) * N * EPS)

      WRITE(*,*) '============'
      WRITE(*,*) 'Checking the Residual of the solution '
      WRITE(*,*) '-- ||Ax-B||_oo/((||A||_oo||x||_oo+||B||_oo).N.eps)=',
     $        RESULT

      IF (( INFO .ne. 0 ) .OR. (RESULT > 60.0)) THEN
          WRITE(*,*) "-- Error in SGESV example !"
      ELSE
          WRITE(*,*) "-- Run of SGESV example successful !"
      ENDIF
*
*     Finalize Plasma
*
      CALL PLASMA_FINALIZE( INFO )
*
*     End of EXAMPLE_SGESV.
*
      END PROGRAM EXAMPLE_SGESV_F

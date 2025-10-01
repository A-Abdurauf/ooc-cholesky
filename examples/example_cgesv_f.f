        PROGRAM EXAMPLE_CGESV_F
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
*  Purpose
*  =======
*
*   FORTRAN EXAMPLE FOR PLASMA_CGESV
*   Example for solving a system of linear equations using LU
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER           CORES, N, NRHS
      PARAMETER         ( CORES = 2 )
      PARAMETER         ( N = 10 )
      PARAMETER         ( NRHS = 5 )
      COMPLEX*8        ZONE
      PARAMETER         ( ZONE = ( 1.0, 0.0 ) )
*     ..
*     .. Local Scalars ..
      COMPLEX*8        A1( N, N ), B1( N, NRHS )
      COMPLEX*8        A2( N, N ), B2( N, NRHS )
      REAL  RWORK( N )
      INTEGER           IPIV( N )
      INTEGER           I, INFO
      INTEGER           ISEED( 4 )
      REAL  XNORM, ANORM, BNORM, RNORM, RESULT, EPS
      REAL  SLAMCH, CLANGE
*     ..
*     .. External Subroutines ..
      EXTERNAL          CLARNV, SLAMCH, CLANGE
      EXTERNAL          PLASMA_INIT
      EXTERNAL          PLASMA_CGESV, PLASMA_FINALIZE
      EXTERNAL          PLASMA_DEALLOC_HANDLE
      EXTERNAL          CGEMM
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
      CALL CLARNV( 1, ISEED, N*N, A1 )
      A2(:,:)=A1(:,:)
*
*     Initialization of the RHS
*
      CALL CLARNV( 1, ISEED, N*NRHS, B1 )
      B2(:,:)=B1(:,:)
*
*     PLASMA CGESV
*
      CALL PLASMA_CGESV( N, NRHS, A2, N, IPIV, B2, N, INFO )
*
*     Check the solution
*
      ANORM = CLANGE('I', N, N,    A1, N, RWORK)
      XNORM = CLANGE('I', N, NRHS, B2, N, RWORK)
      BNORM = CLANGE('I', N, NRHS, B1, N, RWORK)

      CALL CGEMM('No transpose', 'No transpose', N, NRHS, N, ZONE,
     $     A1, N, B2, N, -ZONE, B1, N)

      RNORM = CLANGE('I', N, NRHS, B1, N, RWORK)

      EPS = SLAMCH('Epsilon')
      RESULT = RNORM / ((ANORM * XNORM + BNORM) * N * EPS)

      WRITE(*,*) '============'
      WRITE(*,*) 'Checking the Residual of the solution '
      WRITE(*,*) '-- ||Ax-B||_oo/((||A||_oo||x||_oo+||B||_oo).N.eps)=',
     $        RESULT

      IF (( INFO .ne. 0 ) .OR. (RESULT > 60.0)) THEN
          WRITE(*,*) "-- Error in CGESV example !"
      ELSE
          WRITE(*,*) "-- Run of CGESV example successful !"
      ENDIF
*
*     Finalize Plasma
*
      CALL PLASMA_FINALIZE( INFO )
*
*     End of EXAMPLE_CGESV.
*
      END PROGRAM EXAMPLE_CGESV_F

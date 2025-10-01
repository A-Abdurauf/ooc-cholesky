        PROGRAM EXAMPLE_DGESV_F
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
*  Purpose
*  =======
*
*   FORTRAN EXAMPLE FOR PLASMA_DGESV
*   Example for solving a system of linear equations using LU
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER           CORES, N, NRHS
      PARAMETER         ( CORES = 2 )
      PARAMETER         ( N = 10 )
      PARAMETER         ( NRHS = 5 )
      DOUBLE PRECISION        DONE
      PARAMETER         ( DONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION        A1( N, N ), B1( N, NRHS )
      DOUBLE PRECISION        A2( N, N ), B2( N, NRHS )
      DOUBLE PRECISION  RWORK( N )
      INTEGER           IPIV( N )
      INTEGER           I, INFO
      INTEGER           ISEED( 4 )
      DOUBLE PRECISION  XNORM, ANORM, BNORM, RNORM, RESULT, EPS
      DOUBLE PRECISION  DLAMCH, DLANGE
*     ..
*     .. External Subroutines ..
      EXTERNAL          DLARNV, DLAMCH, DLANGE
      EXTERNAL          PLASMA_INIT
      EXTERNAL          PLASMA_DGESV, PLASMA_FINALIZE
      EXTERNAL          PLASMA_DEALLOC_HANDLE
      EXTERNAL          DGEMM
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
      CALL DLARNV( 1, ISEED, N*N, A1 )
      A2(:,:)=A1(:,:)
*
*     Initialization of the RHS
*
      CALL DLARNV( 1, ISEED, N*NRHS, B1 )
      B2(:,:)=B1(:,:)
*
*     PLASMA DGESV
*
      CALL PLASMA_DGESV( N, NRHS, A2, N, IPIV, B2, N, INFO )
*
*     Check the solution
*
      ANORM = DLANGE('I', N, N,    A1, N, RWORK)
      XNORM = DLANGE('I', N, NRHS, B2, N, RWORK)
      BNORM = DLANGE('I', N, NRHS, B1, N, RWORK)

      CALL DGEMM('No transpose', 'No transpose', N, NRHS, N, DONE,
     $     A1, N, B2, N, -DONE, B1, N)

      RNORM = DLANGE('I', N, NRHS, B1, N, RWORK)

      EPS = DLAMCH('Epsilon')
      RESULT = RNORM / ((ANORM * XNORM + BNORM) * N * EPS)

      WRITE(*,*) '============'
      WRITE(*,*) 'Checking the Residual of the solution '
      WRITE(*,*) '-- ||Ax-B||_oo/((||A||_oo||x||_oo+||B||_oo).N.eps)=',
     $        RESULT

      IF (( INFO .ne. 0 ) .OR. (RESULT > 60.0)) THEN
          WRITE(*,*) "-- Error in DGESV example !"
      ELSE
          WRITE(*,*) "-- Run of DGESV example successful !"
      ENDIF
*
*     Finalize Plasma
*
      CALL PLASMA_FINALIZE( INFO )
*
*     End of EXAMPLE_DGESV.
*
      END PROGRAM EXAMPLE_DGESV_F

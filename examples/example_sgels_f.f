      PROGRAM EXAMPLE_SGELS_F
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
*     Purpose
*     =======
*
*     FORTRAN EXAMPLE FOR PLASMA_SGELS
*     Example for solving a system of linear equations using QR factorization
*
*     =====================================================================
*
*     .. Parameters ..
      INTEGER           CORES, M, N, NRHS, LDA, LDB
      PARAMETER         ( CORES = 2 )
      PARAMETER         ( M = 20 )
      PARAMETER         ( N = 15 )
      PARAMETER         ( NRHS = 5 )
      PARAMETER         ( LDA = 20 )
      PARAMETER         ( LDB = 20 )
      REAL        DONE
      PARAMETER         ( DONE = 1.0 )
*     ..
*     .. Local Scalars ..
      REAL  A1( LDA, N ), B1( LDB, NRHS )
      REAL  A2( LDA, N ), B2( LDB, NRHS )
      REAL  RISU( MAX(M,N), NRHS )
      REAL  RWORK( MAX(M,N) )
      INTEGER           HT( 2 )
      REAL  XNORM, ANORM, BNORM, RNORM, RESULT, EPS
      INTEGER           I, INFO
      INTEGER           ISEED( 4 )
*     ..
*     .. External Subroutines ..
      REAL  SLAMCH, SLANGE
      EXTERNAL          SLARNV, SLAMCH, SLANGE
      EXTERNAL          PLASMA_INIT, PLASMA_ALLOC_WORKSPACE_SGELS
      EXTERNAL          PLASMA_SGELS, PLASMA_FINALIZE
      EXTERNAL          PLASMA_DEALLOC_HANDLE
      EXTERNAL          SGEMM
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
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
*     Initialization of the matrix
*
      CALL SLARNV( 1, ISEED, LDA*N, A1 )
      A2(:,:) = A1(:,:)
*
*     Initialization of the RHS
*
      CALL SLARNV( 1, ISEED, LDB*NRHS, B1 )
      B2(:,:) = B1(:,:)

      RISU(:,:) = 0.
*
*     Allocate T
*
      CALL PLASMA_ALLOC_WORKSPACE_SGELS( M, N, HT, INFO )
*
*     Perform the QR solve
*
      CALL PLASMA_SGELS( PlasmaNoTrans, M, N, NRHS,
     &     A2, LDA, HT, B2, LDB, INFO )
*
*     Check the solution
*
      ANORM = SLANGE('I', M, N,    A1, LDA, RWORK)
      XNORM = SLANGE('I', M, NRHS, B2, LDB, RWORK)
      BNORM = SLANGE('I', M, NRHS, B1, LDB, RWORK)

      CALL SGEMM('No transpose','No transpose', M, NRHS, N, DONE,
     $     A1, LDA, B2, LDB, -DONE, B1, LDB)

      IF ( M >=N ) THEN
         CALL SGEMM('ConjTranspose','No transpose', N, NRHS, M, DONE,
     $        A1, LDA, B1, LDB, -DONE, RISU, M)
         RNORM = SLANGE('I', M, NRHS, RISU, M, RWORK)
      ELSE
         CALL SGEMM('ConjTranspose','No transpose', N, NRHS, M, DONE,
     $        A1, LDA, B1, LDB, -DONE, RISU, N)
         RNORM = SLANGE('I', N, NRHS, RISU, N, RWORK)
      ENDIF

      EPS= SLAMCH('Epsilon')
      RESULT = RNORM / ((ANORM * XNORM + BNORM) * N * EPS)

      WRITE(*,*) '============'
      WRITE(*,*) 'Checking the Residual of the solution '
      WRITE(*,*) '-- ||Ax-B||_oo/((||A||_oo||x||_oo+||B||_oo).N.eps)=',
     $     RESULT

      IF (( INFO .ne. 0 ) .OR. (RESULT > 60.0)) THEN
         WRITE(*,*) "-- Error in SGELS example !"
      ELSE
         WRITE(*,*) "-- Run of SGELS example successful !"
      ENDIF
*
*     Deallocate T
*
      CALL PLASMA_DEALLOC_HANDLE_TILE( HT, INFO )
*
*     Finalize Plasma
*
      CALL PLASMA_FINALIZE( INFO )
*
*     End of EXAMPLE_SGELS.
*
      END PROGRAM EXAMPLE_SGELS_F

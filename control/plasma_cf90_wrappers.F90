!
!     Copyright © 2011 The Numerical Algorithms Group Ltd. All rights reserved.
!   
!     Redistribution and use in source and binary forms, with or without
!     modification, are permitted provided that the following conditions are
!     met:
!     - Redistributions of source code must retain the above copyright notice,
!       this list of conditions, and the following disclaimer.
!     - Redistributions in binary form must reproduce the above copyright
!       notice, this list of conditions and the following disclaimer listed in
!       this license in the documentation and/or other materials provided with
!       the distribution.
!     - Neither the name of the copyright holders nor the names of its
!       contributors may be used to endorse or promote products derived from
!       this software without specific prior written permission.
!     
!     This software is provided by the copyright holders and contributors "as
!     is" and any express or implied warranties, including, but not limited
!     to, the implied warranties of merchantability and fitness for a
!     particular purpose are disclaimed. in no event shall the copyright owner
!     or contributors be liable for any direct, indirect, incidental, special,
!     exemplary, or consequential damages (including, but not limited to,
!     procurement of substitute goods or services; loss of use, data, or
!     profits; or business interruption) however caused and on any theory of
!     liability, whether in contract, strict liability, or tort (including
!     negligence or otherwise) arising in any way out of the use of this
!     software, even if advised of the possibility of such damage.
!
!
!
! @file plasma_cf90_wrappers.F90
!
!  PLASMA fortran wrapper for BLAS and LAPACK subroutines.
!  PLASMA is a software package provided by Univ. of Tennessee,
!  Univ. of California Berkeley and Univ. of Colorado Denver
!
! @version 2.8.0
! @author Numerical Algorithm Group
! @author Mathieu Faverge
! @date 2011-09-15
! @generated c Fri Apr  1 11:03:02 2016
!
!
! Wrappers to PLASMA functions are provided for the following BLAS
! subroutines since the PLASMA and BLAS interfaces match exactly:
! CGEMM  PLASMA_cgemm
! CHEMM  PLASMA_chemm
! CHER2K PLASMA_cher2k
! CHERK  PLASMA_cherk
! CSYMM  PLASMA_csymm
! CSYR2K PLASMA_csyr2k
! CSYRK  PLASMA_csyrk
! CTRMM  PLASMA_ctrmm
! CTRSM  PLASMA_ctrsm
!
! Wrappers to PLASMA functions are provided for the following LAPACK
! subroutines since the PLASMA and LAPACK interfaces match exactly:
! CGESV  PLASMA_cgesv
! CGETRF PLASMA_cgetrf
! CGETRS PLASMA_cgetrs
! CHEGST PLASMA_chegst
! CLASWP PLASMA_claswp
! CLAUUM PLASMA_clauum
! CPOSV  PLASMA_cposv
! CPOTRF PLASMA_cpotrf
! CPOTRI PLASMA_cpotri
! CPOTRS PLASMA_cpotrs
! CTRTRI PLASMA_ctrtri
! CLACPY PLASMA_clacpy
! CLASET PLASMA_claset
#define PRECISION_c

      subroutine plasma_wrap_CGESV(N,NRHS,A,LDA,IPIV,B,LDB,INFO)
            use iso_c_binding
            use plasma
            implicit none
            integer, parameter :: wp = kind(0.0d0)
            integer, intent(in) :: LDA
            integer, intent(in) :: LDB
            integer, intent(in) :: N
            integer, intent(in) :: NRHS
            integer, intent(out) :: INFO
            integer, intent(out), target :: IPIV(*)
            complex, intent(inout), target :: A(LDA,*)
            complex, intent(inout), target :: B(LDB,*)
            if (.not. plasma_initialized) call plasma_init(24,INFO)
            ! write(*,*) " Calling PLASMA_CGESV"
            call PLASMA_CGESV(N,NRHS,A,LDA,IPIV,B,LDB,INFO)
      end subroutine plasma_wrap_CGESV

      subroutine plasma_wrap_CGETRF(M,N,A,LDA,IPIV,INFO)
            use iso_c_binding
            use plasma
            implicit none
            integer, parameter :: wp = kind(0.0d0)
            integer, intent(in) :: LDA
            integer, intent(in) :: M
            integer, intent(in) :: N
            integer, intent(out) :: INFO
            integer, intent(out), target :: IPIV(*)
            complex, intent(inout), target :: A(LDA,*)
            if (.not. plasma_initialized) call plasma_init(24,INFO)
            ! write(*,*) " Calling PLASMA_CGETRF"
            call PLASMA_CGETRF(M,N,A,LDA,IPIV,INFO)
      end subroutine plasma_wrap_CGETRF

      subroutine plasma_wrap_CGETRS(TRANS,N,NRHS,A,LDA,IPIV,B,LDB,INFO)
            use iso_c_binding
            use plasma
            implicit none
            integer, parameter :: wp = kind(0.0d0)
            integer, intent(in) :: LDA
            integer, intent(in) :: LDB
            integer, intent(in) :: N
            integer, intent(in) :: NRHS
            integer, intent(in), target :: IPIV(*)
            integer, intent(out) :: INFO
            character, intent(in) :: TRANS
            complex, intent(inout), target :: A(LDA,*)
            complex, intent(inout), target :: B(LDB,*)
            integer :: local_TRANS
            if(TRANS=='N' .or. TRANS=='n')then
               local_TRANS = PlasmaNoTrans
            else if(TRANS=='T' .or. TRANS=='t')then
               local_TRANS = PlasmaTrans
            else if(TRANS=='C' .or. TRANS=='c')then
               local_TRANS = PlasmaConjTrans
            else
               local_TRANS = -1
            end if
            if (.not. plasma_initialized) call plasma_init(24,INFO)
            ! write(*,*) " Calling PLASMA_CGETRS"
            call PLASMA_CGETRS(local_TRANS,N,NRHS,A,LDA,IPIV,B,LDB,INFO)
      end subroutine plasma_wrap_CGETRS

      subroutine plasma_wrap_CHEGST(ITYPE,UPLO,N,A,LDA,B,LDB,INFO)
            use iso_c_binding
            use plasma
            implicit none
            integer, parameter :: wp = kind(0.0d0)
            integer, intent(in) :: ITYPE
            integer, intent(in) :: LDA
            integer, intent(in) :: LDB
            integer, intent(in) :: N
            integer, intent(out) :: INFO
            character, intent(in) :: UPLO
            complex, intent(inout), target :: B(LDB,*)
            complex, intent(inout), target :: A(LDA,*)
            integer :: local_UPLO
            if(UPLO=='U' .or. UPLO=='u')then
               local_UPLO = PlasmaUpper
            else if(UPLO=='L' .or. UPLO=='l')then
               local_UPLO = PlasmaLower
            else
               local_UPLO = -1
            end if
            if (.not. plasma_initialized) call plasma_init(24,INFO)
            ! write(*,*) " Calling PLASMA_CHEGST"
            call PLASMA_CHEGST(ITYPE,local_UPLO,N,A,LDA,B,LDB,INFO)
      end subroutine plasma_wrap_CHEGST

      subroutine plasma_wrap_CLASWP(N,A,LDA,K1,K2,IPIV,INCX)
            use iso_c_binding
            use plasma
            implicit none
            integer, parameter :: wp = kind(0.0d0)
            integer, intent(in) :: INCX
            integer, intent(in) :: K1
            integer, intent(in) :: K2
            integer, intent(in) :: LDA
            integer, intent(in) :: N
            integer, intent(in), target :: IPIV(*)
            complex, intent(inout), target :: A(LDA,*)
            integer :: local_ret
            if (.not. plasma_initialized) call plasma_init(24,local_ret)
            ! write(*,*) " Calling PLASMA_CLASWP"
            call PLASMA_CLASWP(N,A,LDA,K1,K2,IPIV,INCX,local_ret)
      end subroutine plasma_wrap_CLASWP

      subroutine plasma_wrap_CLAUUM(UPLO,N,A,LDA,INFO)
            use iso_c_binding
            use plasma
            implicit none
            integer, parameter :: wp = kind(0.0d0)
            integer, intent(in) :: LDA
            integer, intent(in) :: N
            integer, intent(out) :: INFO
            character, intent(in) :: UPLO
            complex, intent(inout), target :: A(LDA,*)
            integer :: local_UPLO
            if(UPLO=='U' .or. UPLO=='u')then
               local_UPLO = PlasmaUpper
            else if(UPLO=='L' .or. UPLO=='l')then
               local_UPLO = PlasmaLower
            else
               local_UPLO = -1
            end if
            if (.not. plasma_initialized) call plasma_init(24,INFO)
            ! write(*,*) " Calling PLASMA_CLAUUM"
            call PLASMA_CLAUUM(local_UPLO,N,A,LDA,INFO)
      end subroutine plasma_wrap_CLAUUM

      subroutine plasma_wrap_CPOSV(UPLO,N,NRHS,A,LDA,B,LDB,INFO)
            use iso_c_binding
            use plasma
            implicit none
            integer, parameter :: wp = kind(0.0d0)
            integer, intent(in) :: LDA
            integer, intent(in) :: LDB
            integer, intent(in) :: N
            integer, intent(in) :: NRHS
            integer, intent(out) :: INFO
            character, intent(in) :: UPLO
            complex, intent(inout), target :: A(LDA,*)
            complex, intent(inout), target :: B(LDB,*)
            integer :: local_UPLO
            if(UPLO=='U' .or. UPLO=='u')then
               local_UPLO = PlasmaUpper
            else if(UPLO=='L' .or. UPLO=='l')then
               local_UPLO = PlasmaLower
            else
               local_UPLO = -1
            end if
            if (.not. plasma_initialized) call plasma_init(24,INFO)
            ! write(*,*) " Calling PLASMA_CPOSV"
            call PLASMA_CPOSV(local_UPLO,N,NRHS,A,LDA,B,LDB,INFO)
      end subroutine plasma_wrap_CPOSV

      subroutine plasma_wrap_CPOTRF(UPLO,N,A,LDA,INFO)
            use iso_c_binding
            use plasma
            implicit none
            integer, parameter :: wp = kind(0.0d0)
            integer, intent(in) :: LDA
            integer, intent(in) :: N
            integer, intent(out) :: INFO
            character, intent(in) :: UPLO
            complex, intent(inout), target :: A(LDA,*)
            integer :: local_UPLO
            if(UPLO=='U' .or. UPLO=='u')then
               local_UPLO = PlasmaUpper
            else if(UPLO=='L' .or. UPLO=='l')then
               local_UPLO = PlasmaLower
            else
               local_UPLO = -1
            end if
            if (.not. plasma_initialized) call plasma_init(24,INFO)
            ! write(*,*) " Calling PLASMA_CPOTRF"
            call PLASMA_CPOTRF(local_UPLO,N,A,LDA,INFO)
      end subroutine plasma_wrap_CPOTRF

      subroutine plasma_wrap_CPOTRI(UPLO,N,A,LDA,INFO)
            use iso_c_binding
            use plasma
            implicit none
            integer, parameter :: wp = kind(0.0d0)
            integer, intent(in) :: LDA
            integer, intent(in) :: N
            integer, intent(out) :: INFO
            character, intent(in) :: UPLO
            complex, intent(inout), target :: A(LDA,*)
            integer :: local_UPLO
            if(UPLO=='U' .or. UPLO=='u')then
               local_UPLO = PlasmaUpper
            else if(UPLO=='L' .or. UPLO=='l')then
               local_UPLO = PlasmaLower
            else
               local_UPLO = -1
            end if
            if (.not. plasma_initialized) call plasma_init(24,INFO)
            ! write(*,*) " Calling PLASMA_CPOTRI"
            call PLASMA_CPOTRI(local_UPLO,N,A,LDA,INFO)
      end subroutine plasma_wrap_CPOTRI

      subroutine plasma_wrap_CPOTRS(UPLO,N,NRHS,A,LDA,B,LDB,INFO)
            use iso_c_binding
            use plasma
            implicit none
            integer, parameter :: wp = kind(0.0d0)
            integer, intent(in) :: LDA
            integer, intent(in) :: LDB
            integer, intent(in) :: N
            integer, intent(in) :: NRHS
            integer, intent(out) :: INFO
            character, intent(in) :: UPLO
            complex, intent(inout), target :: A(LDA,*)
            complex, intent(inout), target :: B(LDB,*)
            integer :: local_UPLO
            if(UPLO=='U' .or. UPLO=='u')then
               local_UPLO = PlasmaUpper
            else if(UPLO=='L' .or. UPLO=='l')then
               local_UPLO = PlasmaLower
            else
               local_UPLO = -1
            end if
            if (.not. plasma_initialized) call plasma_init(24,INFO)
            ! write(*,*) " Calling PLASMA_CPOTRS"
            call PLASMA_CPOTRS(local_UPLO,N,NRHS,A,LDA,B,LDB,INFO)
      end subroutine plasma_wrap_CPOTRS

      subroutine plasma_wrap_CTRTRI(UPLO,DIAG,N,A,LDA,INFO)
            use iso_c_binding
            use plasma
            implicit none
            integer, parameter :: wp = kind(0.0d0)
            integer, intent(in) :: LDA
            integer, intent(in) :: N
            integer, intent(out) :: INFO
            character, intent(in) :: DIAG
            character, intent(in) :: UPLO
            complex, intent(inout), target :: A(LDA,*)
            integer :: local_DIAG
            integer :: local_UPLO
            if(UPLO=='U' .or. UPLO=='u')then
               local_UPLO = PlasmaUpper
            else if(UPLO=='L' .or. UPLO=='l')then
               local_UPLO = PlasmaLower
            else
               local_UPLO = -1
            end if
            if(DIAG=='U' .or. DIAG=='u')then
               local_DIAG = PlasmaUnit
            else if(DIAG=='N' .or. DIAG=='n')then
               local_DIAG = PlasmaNonUnit
            else
               local_DIAG = -1
            end if
            if (.not. plasma_initialized) call plasma_init(24,INFO)
            ! write(*,*) " Calling PLASMA_CTRTRI"
            call PLASMA_CTRTRI(local_UPLO,local_DIAG,N,A,LDA,INFO)
      end subroutine plasma_wrap_CTRTRI

      subroutine plasma_wrap_CGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
            use iso_c_binding
            use plasma
            implicit none
            integer, parameter :: wp = kind(0.0d0)
            integer, intent(in) :: K
            integer, intent(in) :: LDA
            integer, intent(in) :: LDB
            integer, intent(in) :: LDC
            integer, intent(in) :: M
            integer, intent(in) :: N
            character, intent(in) :: TRANSA
            character, intent(in) :: TRANSB
            complex, intent(in) :: ALPHA
            complex, intent(in) :: BETA
            complex, intent(inout), target :: A(LDA,*)
            complex, intent(inout), target :: B(LDB,*)
            complex, intent(inout), target :: C(LDC,*)
            integer :: local_TRANSA
            integer :: local_TRANSB
            integer :: local_ret
            if(TRANSA=='N' .or. TRANSA=='n')then
               local_TRANSA = PlasmaNoTrans
            else if(TRANSA=='T' .or. TRANSA=='t')then
               local_TRANSA = PlasmaTrans
            else if(TRANSA=='C' .or. TRANSA=='c')then
               local_TRANSA = PlasmaConjTrans
            else
               local_TRANSA = -1
            end if
            if(TRANSB=='N' .or. TRANSB=='n')then
               local_TRANSB = PlasmaNoTrans
            else if(TRANSB=='T' .or. TRANSB=='t')then
               local_TRANSB = PlasmaTrans
            else if(TRANSB=='C' .or. TRANSB=='c')then
               local_TRANSB = PlasmaConjTrans
            else
               local_TRANSB = -1
            end if
            if (.not. plasma_initialized) call plasma_init(24,local_ret)
            ! write(*,*) " Calling PLASMA_CGEMM"
            call PLASMA_CGEMM(local_TRANSA,local_TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC,local_ret)
      end subroutine plasma_wrap_CGEMM

#if defined(PRECISION_z) || defined(PRECISION_c)
      subroutine plasma_wrap_CHEMM(SIDE,UPLO,M,N,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
            use iso_c_binding
            use plasma
            implicit none
            integer, parameter :: wp = kind(0.0d0)
            integer, intent(in) :: LDA
            integer, intent(in) :: LDB
            integer, intent(in) :: LDC
            integer, intent(in) :: M
            integer, intent(in) :: N
            character, intent(in) :: SIDE
            character, intent(in) :: UPLO
            complex, intent(in) :: ALPHA
            complex, intent(in) :: BETA
            complex, intent(inout), target :: A(LDA,*)
            complex, intent(inout), target :: B(LDB,*)
            complex, intent(inout), target :: C(LDC,*)
            integer :: local_SIDE
            integer :: local_UPLO
            integer :: local_ret
            if(SIDE=='L' .or. SIDE=='l')then
               local_SIDE = PlasmaLeft
            else if(SIDE=='R' .or. SIDE=='r')then
               local_SIDE = PlasmaRight
            else
               local_SIDE = -1
            end if
            if(UPLO=='U' .or. UPLO=='u')then
               local_UPLO = PlasmaUpper
            else if(UPLO=='L' .or. UPLO=='l')then
               local_UPLO = PlasmaLower
            else
               local_UPLO = -1
            end if
            if (.not. plasma_initialized) call plasma_init(24,local_ret)
            ! write(*,*) " Calling PLASMA_CHEMM"
            call PLASMA_CHEMM(local_SIDE,local_UPLO,M,N,ALPHA,A,LDA,B,LDB,BETA,C,LDC,local_ret)
      end subroutine plasma_wrap_CHEMM

      subroutine plasma_wrap_CHER2K(UPLO,TRANS,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
            use iso_c_binding
            use plasma
            implicit none
            integer, parameter :: wp = kind(0.0d0)
            integer, intent(in) :: K
            integer, intent(in) :: LDA
            integer, intent(in) :: LDB
            integer, intent(in) :: LDC
            integer, intent(in) :: N
            character, intent(in) :: TRANS
            character, intent(in) :: UPLO
            complex, intent(in) :: ALPHA
            complex, intent(inout), target :: A(LDA,*)
            complex, intent(inout), target :: B(LDB,*)
            complex, intent(inout), target :: C(LDC,*)
            real, intent(in) :: BETA
            integer :: local_TRANS
            integer :: local_UPLO
            integer :: local_ret
            if(UPLO=='U' .or. UPLO=='u')then
               local_UPLO = PlasmaUpper
            else if(UPLO=='L' .or. UPLO=='l')then
               local_UPLO = PlasmaLower
            else
               local_UPLO = -1
            end if
            if(TRANS=='N' .or. TRANS=='n')then
               local_TRANS = PlasmaNoTrans
            else if(TRANS=='T' .or. TRANS=='t')then
               local_TRANS = PlasmaTrans
            else if(TRANS=='C' .or. TRANS=='c')then
               local_TRANS = PlasmaConjTrans
            else
               local_TRANS = -1
            end if
            if (.not. plasma_initialized) call plasma_init(24,local_ret)
            ! write(*,*) " Calling PLASMA_CHER2K"
            call PLASMA_CHER2K(local_UPLO,local_TRANS,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC,local_ret)
      end subroutine plasma_wrap_CHER2K

      subroutine plasma_wrap_CHERK(UPLO,TRANS,N,K,ALPHA,A,LDA,BETA,C,LDC)
            use iso_c_binding
            use plasma
            implicit none
            integer, parameter :: wp = kind(0.0d0)
            integer, intent(in) :: K
            integer, intent(in) :: LDA
            integer, intent(in) :: LDC
            integer, intent(in) :: N
            character, intent(in) :: TRANS
            character, intent(in) :: UPLO
            complex, intent(inout), target :: A(LDA,*)
            complex, intent(inout), target :: C(LDC,*)
            real, intent(in) :: ALPHA
            real, intent(in) :: BETA
            integer :: local_TRANS
            integer :: local_UPLO
            integer :: local_ret
            if(UPLO=='U' .or. UPLO=='u')then
               local_UPLO = PlasmaUpper
            else if(UPLO=='L' .or. UPLO=='l')then
               local_UPLO = PlasmaLower
            else
               local_UPLO = -1
            end if
            if(TRANS=='N' .or. TRANS=='n')then
               local_TRANS = PlasmaNoTrans
            else if(TRANS=='T' .or. TRANS=='t')then
               local_TRANS = PlasmaTrans
            else if(TRANS=='C' .or. TRANS=='c')then
               local_TRANS = PlasmaConjTrans
            else
               local_TRANS = -1
            end if
            if (.not. plasma_initialized) call plasma_init(24,local_ret)
            ! write(*,*) " Calling PLASMA_CHERK"
            call PLASMA_CHERK(local_UPLO,local_TRANS,N,K,ALPHA,A,LDA,BETA,C,LDC,local_ret)
      end subroutine plasma_wrap_CHERK
#endif

      subroutine plasma_wrap_CSYMM(SIDE,UPLO,M,N,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
            use iso_c_binding
            use plasma
            implicit none
            integer, parameter :: wp = kind(0.0d0)
            integer, intent(in) :: LDA
            integer, intent(in) :: LDB
            integer, intent(in) :: LDC
            integer, intent(in) :: M
            integer, intent(in) :: N
            character, intent(in) :: SIDE
            character, intent(in) :: UPLO
            complex, intent(in) :: ALPHA
            complex, intent(in) :: BETA
            complex, intent(inout), target :: A(LDA,*)
            complex, intent(inout), target :: B(LDB,*)
            complex, intent(inout), target :: C(LDC,*)
            integer :: local_SIDE
            integer :: local_UPLO
            integer :: local_ret
            if(SIDE=='L' .or. SIDE=='l')then
               local_SIDE = PlasmaLeft
            else if(SIDE=='R' .or. SIDE=='r')then
               local_SIDE = PlasmaRight
            else
               local_SIDE = -1
            end if
            if(UPLO=='U' .or. UPLO=='u')then
               local_UPLO = PlasmaUpper
            else if(UPLO=='L' .or. UPLO=='l')then
               local_UPLO = PlasmaLower
            else
               local_UPLO = -1
            end if
            if (.not. plasma_initialized) call plasma_init(24,local_ret)
            ! write(*,*) " Calling PLASMA_CSYMM"
            call PLASMA_CSYMM(local_SIDE,local_UPLO,M,N,ALPHA,A,LDA,B,LDB,BETA,C,LDC,local_ret)
      end subroutine plasma_wrap_CSYMM

      subroutine plasma_wrap_CSYR2K(UPLO,TRANS,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
            use iso_c_binding
            use plasma
            implicit none
            integer, parameter :: wp = kind(0.0d0)
            integer, intent(in) :: K
            integer, intent(in) :: LDA
            integer, intent(in) :: LDB
            integer, intent(in) :: LDC
            integer, intent(in) :: N
            character, intent(in) :: TRANS
            character, intent(in) :: UPLO
            complex, intent(in) :: ALPHA
            complex, intent(in) :: BETA
            complex, intent(inout), target :: A(LDA,*)
            complex, intent(inout), target :: B(LDB,*)
            complex, intent(inout), target :: C(LDC,*)
            integer :: local_TRANS
            integer :: local_UPLO
            integer :: local_ret
            if(UPLO=='U' .or. UPLO=='u')then
               local_UPLO = PlasmaUpper
            else if(UPLO=='L' .or. UPLO=='l')then
               local_UPLO = PlasmaLower
            else
               local_UPLO = -1
            end if
            if(TRANS=='N' .or. TRANS=='n')then
               local_TRANS = PlasmaNoTrans
            else if(TRANS=='T' .or. TRANS=='t')then
               local_TRANS = PlasmaTrans
            else if(TRANS=='C' .or. TRANS=='c')then
               local_TRANS = PlasmaConjTrans
            else
               local_TRANS = -1
            end if
            if (.not. plasma_initialized) call plasma_init(24,local_ret)
            ! write(*,*) " Calling PLASMA_CSYR2K"
            call PLASMA_CSYR2K(local_UPLO,local_TRANS,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC,local_ret)
      end subroutine plasma_wrap_CSYR2K

      subroutine plasma_wrap_CSYRK(UPLO,TRANS,N,K,ALPHA,A,LDA,BETA,C,LDC)
            use iso_c_binding
            use plasma
            implicit none
            integer, parameter :: wp = kind(0.0d0)
            integer, intent(in) :: K
            integer, intent(in) :: LDA
            integer, intent(in) :: LDC
            integer, intent(in) :: N
            character, intent(in) :: TRANS
            character, intent(in) :: UPLO
            complex, intent(in) :: ALPHA
            complex, intent(in) :: BETA
            complex, intent(inout), target :: A(LDA,*)
            complex, intent(inout), target :: C(LDC,*)
            integer :: local_TRANS
            integer :: local_UPLO
            integer :: local_ret
            if(UPLO=='U' .or. UPLO=='u')then
               local_UPLO = PlasmaUpper
            else if(UPLO=='L' .or. UPLO=='l')then
               local_UPLO = PlasmaLower
            else
               local_UPLO = -1
            end if
            if(TRANS=='N' .or. TRANS=='n')then
               local_TRANS = PlasmaNoTrans
            else if(TRANS=='T' .or. TRANS=='t')then
               local_TRANS = PlasmaTrans
            else if(TRANS=='C' .or. TRANS=='c')then
               local_TRANS = PlasmaConjTrans
            else
               local_TRANS = -1
            end if
            if (.not. plasma_initialized) call plasma_init(24,local_ret)
            ! write(*,*) " Calling PLASMA_CSYRK"
            call PLASMA_CSYRK(local_UPLO,local_TRANS,N,K,ALPHA,A,LDA,BETA,C,LDC,local_ret)
      end subroutine plasma_wrap_CSYRK

      subroutine plasma_wrap_CTRMM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
            use iso_c_binding
            use plasma
            implicit none
            integer, parameter :: wp = kind(0.0d0)
            integer, intent(in) :: LDA
            integer, intent(in) :: LDB
            integer, intent(in) :: M
            integer, intent(in) :: N
            character, intent(in) :: DIAG
            character, intent(in) :: SIDE
            character, intent(in) :: TRANSA
            character, intent(in) :: UPLO
            complex, intent(in) :: ALPHA
            complex, intent(inout), target :: A(LDA,*)
            complex, intent(inout), target :: B(LDB,*)
            integer :: local_DIAG
            integer :: local_SIDE
            integer :: local_TRANSA
            integer :: local_UPLO
            integer :: local_ret
            if(SIDE=='L' .or. SIDE=='l')then
               local_SIDE = PlasmaLeft
            else if(SIDE=='R' .or. SIDE=='r')then
               local_SIDE = PlasmaRight
            else
               local_SIDE = -1
            end if
            if(UPLO=='U' .or. UPLO=='u')then
               local_UPLO = PlasmaUpper
            else if(UPLO=='L' .or. UPLO=='l')then
               local_UPLO = PlasmaLower
            else
               local_UPLO = -1
            end if
            if(TRANSA=='N' .or. TRANSA=='n')then
               local_TRANSA = PlasmaNoTrans
            else if(TRANSA=='T' .or. TRANSA=='t')then
               local_TRANSA = PlasmaTrans
            else if(TRANSA=='C' .or. TRANSA=='c')then
               local_TRANSA = PlasmaConjTrans
            else
               local_TRANSA = -1
            end if
            if(DIAG=='U' .or. DIAG=='u')then
               local_DIAG = PlasmaUnit
            else if(DIAG=='N' .or. DIAG=='n')then
               local_DIAG = PlasmaNonUnit
            else
               local_DIAG = -1
            end if
            if (.not. plasma_initialized) call plasma_init(24,local_ret)
            ! write(*,*) " Calling PLASMA_CTRMM"
            call PLASMA_CTRMM(local_SIDE,local_UPLO,local_TRANSA,local_DIAG,M,N,ALPHA,A,LDA,B,LDB,local_ret)
      end subroutine plasma_wrap_CTRMM

      subroutine plasma_wrap_CTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
            use iso_c_binding
            use plasma
            implicit none
            integer, parameter :: wp = kind(0.0d0)
            integer, intent(in) :: LDA
            integer, intent(in) :: LDB
            integer, intent(in) :: M
            integer, intent(in) :: N
            character, intent(in) :: DIAG
            character, intent(in) :: SIDE
            character, intent(in) :: TRANSA
            character, intent(in) :: UPLO
            complex, intent(in) :: ALPHA
            complex, intent(inout), target :: A(LDA,*)
            complex, intent(inout), target :: B(LDB,*)
            integer :: local_DIAG
            integer :: local_SIDE
            integer :: local_TRANSA
            integer :: local_UPLO
            integer :: local_ret
            if(SIDE=='L' .or. SIDE=='l')then
               local_SIDE = PlasmaLeft
            else if(SIDE=='R' .or. SIDE=='r')then
               local_SIDE = PlasmaRight
            else
               local_SIDE = -1
            end if
            if(UPLO=='U' .or. UPLO=='u')then
               local_UPLO = PlasmaUpper
            else if(UPLO=='L' .or. UPLO=='l')then
               local_UPLO = PlasmaLower
            else
               local_UPLO = -1
            end if
            if(TRANSA=='N' .or. TRANSA=='n')then
               local_TRANSA = PlasmaNoTrans
            else if(TRANSA=='T' .or. TRANSA=='t')then
               local_TRANSA = PlasmaTrans
            else if(TRANSA=='C' .or. TRANSA=='c')then
               local_TRANSA = PlasmaConjTrans
            else
               local_TRANSA = -1
            end if
            if(DIAG=='U' .or. DIAG=='u')then
               local_DIAG = PlasmaUnit
            else if(DIAG=='N' .or. DIAG=='n')then
               local_DIAG = PlasmaNonUnit
            else
               local_DIAG = -1
            end if
            if (.not. plasma_initialized) call plasma_init(24,local_ret)
            ! write(*,*) " Calling PLASMA_CTRSM"
            call PLASMA_CTRSM(local_SIDE,local_UPLO,local_TRANSA,local_DIAG,M,N,ALPHA,A,LDA,B,LDB,local_ret)
      end subroutine plasma_wrap_CTRSM

      subroutine plasma_wrap_CLACPY(UPLO,M,N,A,LDA,B,LDB)
            use iso_c_binding
            use plasma
            implicit none
            integer, parameter :: wp = kind(0.0d0)
            integer, intent(in) :: LDA
            integer, intent(in) :: LDB
            integer, intent(in) :: M
            integer, intent(in) :: N
            character, intent(in) :: UPLO
            complex, intent(inout), target :: A(LDA,*)
            complex, intent(out), target :: B(LDB,*)
            integer :: local_UPLO
            integer :: local_ret
            if(UPLO=='U' .or. UPLO=='u')then
               local_UPLO = PlasmaUpper
            else if(UPLO=='L' .or. UPLO=='l')then
               local_UPLO = PlasmaLower
            else
               local_UPLO = -1
            end if
            if (.not. plasma_initialized) call plasma_init(24,local_ret)
            ! write(*,*) " Calling PLASMA_CLACPY"
            call PLASMA_CLACPY(local_UPLO,M,N,A,LDA,B,LDB,local_ret)
      end subroutine plasma_wrap_CLACPY

      subroutine plasma_wrap_CLASET(UPLO,M,N,ALPHA,BETA,A,LDA)
            use iso_c_binding
            use plasma
            implicit none
            integer, parameter :: wp = kind(0.0d0)
            integer, intent(in) :: LDA
            integer, intent(in) :: M
            integer, intent(in) :: N
            character, intent(in) :: UPLO
            complex, intent(in) :: ALPHA
            complex, intent(in) :: BETA
            complex, intent(inout), target :: A(LDA,*)
            integer :: local_UPLO
            integer :: local_ret
            if(UPLO=='U' .or. UPLO=='u')then
               local_UPLO = PlasmaUpper
            else if(UPLO=='L' .or. UPLO=='l')then
               local_UPLO = PlasmaLower
            else
               local_UPLO = -1
            end if
            if (.not. plasma_initialized) call plasma_init(24,local_ret)
            ! write(*,*) " Calling PLASMA_CLASET"
            call PLASMA_CLASET(local_UPLO,M,N,ALPHA,BETA,A,LDA,local_ret)
      end subroutine plasma_wrap_CLASET

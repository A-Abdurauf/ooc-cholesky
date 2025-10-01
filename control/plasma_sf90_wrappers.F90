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
! @file plasma_sf90_wrappers.F90
!
!  PLASMA fortran wrapper for BLAS and LAPACK subroutines.
!  PLASMA is a software package provided by Univ. of Tennessee,
!  Univ. of California Berkeley and Univ. of Colorado Denver
!
! @version 2.8.0
! @author Numerical Algorithm Group
! @author Mathieu Faverge
! @date 2011-09-15
! @generated s Fri Apr  1 11:03:02 2016
!
!
! Wrappers to PLASMA functions are provided for the following BLAS
! subroutines since the PLASMA and BLAS interfaces match exactly:
! SGEMM  PLASMA_sgemm
! SSYMM  PLASMA_ssymm
! SSYR2K PLASMA_ssyr2k
! SSYRK  PLASMA_ssyrk
! SSYMM  PLASMA_ssymm
! SSYR2K PLASMA_ssyr2k
! SSYRK  PLASMA_ssyrk
! STRMM  PLASMA_strmm
! STRSM  PLASMA_strsm
!
! Wrappers to PLASMA functions are provided for the following LAPACK
! subroutines since the PLASMA and LAPACK interfaces match exactly:
! SGESV  PLASMA_sgesv
! SGETRF PLASMA_sgetrf
! SGETRS PLASMA_sgetrs
! SSYGST PLASMA_ssygst
! SLASWP PLASMA_slaswp
! SLAUUM PLASMA_slauum
! SPOSV  PLASMA_sposv
! SPOTRF PLASMA_spotrf
! SPOTRI PLASMA_spotri
! SPOTRS PLASMA_spotrs
! STRTRI PLASMA_strtri
! SLACPY PLASMA_slacpy
! SLASET PLASMA_slaset
#define PRECISION_s

      subroutine plasma_wrap_SGESV(N,NRHS,A,LDA,IPIV,B,LDB,INFO)
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
            real, intent(inout), target :: A(LDA,*)
            real, intent(inout), target :: B(LDB,*)
            if (.not. plasma_initialized) call plasma_init(24,INFO)
            ! write(*,*) " Calling PLASMA_SGESV"
            call PLASMA_SGESV(N,NRHS,A,LDA,IPIV,B,LDB,INFO)
      end subroutine plasma_wrap_SGESV

      subroutine plasma_wrap_SGETRF(M,N,A,LDA,IPIV,INFO)
            use iso_c_binding
            use plasma
            implicit none
            integer, parameter :: wp = kind(0.0d0)
            integer, intent(in) :: LDA
            integer, intent(in) :: M
            integer, intent(in) :: N
            integer, intent(out) :: INFO
            integer, intent(out), target :: IPIV(*)
            real, intent(inout), target :: A(LDA,*)
            if (.not. plasma_initialized) call plasma_init(24,INFO)
            ! write(*,*) " Calling PLASMA_SGETRF"
            call PLASMA_SGETRF(M,N,A,LDA,IPIV,INFO)
      end subroutine plasma_wrap_SGETRF

      subroutine plasma_wrap_SGETRS(TRANS,N,NRHS,A,LDA,IPIV,B,LDB,INFO)
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
            real, intent(inout), target :: A(LDA,*)
            real, intent(inout), target :: B(LDB,*)
            integer :: local_TRANS
            if(TRANS=='N' .or. TRANS=='n')then
               local_TRANS = PlasmaNoTrans
            else if(TRANS=='T' .or. TRANS=='t')then
               local_TRANS = PlasmaTrans
            else if(TRANS=='C' .or. TRANS=='c')then
               local_TRANS = PlasmaTrans
            else
               local_TRANS = -1
            end if
            if (.not. plasma_initialized) call plasma_init(24,INFO)
            ! write(*,*) " Calling PLASMA_SGETRS"
            call PLASMA_SGETRS(local_TRANS,N,NRHS,A,LDA,IPIV,B,LDB,INFO)
      end subroutine plasma_wrap_SGETRS

      subroutine plasma_wrap_SSYGST(ITYPE,UPLO,N,A,LDA,B,LDB,INFO)
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
            real, intent(inout), target :: B(LDB,*)
            real, intent(inout), target :: A(LDA,*)
            integer :: local_UPLO
            if(UPLO=='U' .or. UPLO=='u')then
               local_UPLO = PlasmaUpper
            else if(UPLO=='L' .or. UPLO=='l')then
               local_UPLO = PlasmaLower
            else
               local_UPLO = -1
            end if
            if (.not. plasma_initialized) call plasma_init(24,INFO)
            ! write(*,*) " Calling PLASMA_SSYGST"
            call PLASMA_SSYGST(ITYPE,local_UPLO,N,A,LDA,B,LDB,INFO)
      end subroutine plasma_wrap_SSYGST

      subroutine plasma_wrap_SLASWP(N,A,LDA,K1,K2,IPIV,INCX)
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
            real, intent(inout), target :: A(LDA,*)
            integer :: local_ret
            if (.not. plasma_initialized) call plasma_init(24,local_ret)
            ! write(*,*) " Calling PLASMA_SLASWP"
            call PLASMA_SLASWP(N,A,LDA,K1,K2,IPIV,INCX,local_ret)
      end subroutine plasma_wrap_SLASWP

      subroutine plasma_wrap_SLAUUM(UPLO,N,A,LDA,INFO)
            use iso_c_binding
            use plasma
            implicit none
            integer, parameter :: wp = kind(0.0d0)
            integer, intent(in) :: LDA
            integer, intent(in) :: N
            integer, intent(out) :: INFO
            character, intent(in) :: UPLO
            real, intent(inout), target :: A(LDA,*)
            integer :: local_UPLO
            if(UPLO=='U' .or. UPLO=='u')then
               local_UPLO = PlasmaUpper
            else if(UPLO=='L' .or. UPLO=='l')then
               local_UPLO = PlasmaLower
            else
               local_UPLO = -1
            end if
            if (.not. plasma_initialized) call plasma_init(24,INFO)
            ! write(*,*) " Calling PLASMA_SLAUUM"
            call PLASMA_SLAUUM(local_UPLO,N,A,LDA,INFO)
      end subroutine plasma_wrap_SLAUUM

      subroutine plasma_wrap_SPOSV(UPLO,N,NRHS,A,LDA,B,LDB,INFO)
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
            real, intent(inout), target :: A(LDA,*)
            real, intent(inout), target :: B(LDB,*)
            integer :: local_UPLO
            if(UPLO=='U' .or. UPLO=='u')then
               local_UPLO = PlasmaUpper
            else if(UPLO=='L' .or. UPLO=='l')then
               local_UPLO = PlasmaLower
            else
               local_UPLO = -1
            end if
            if (.not. plasma_initialized) call plasma_init(24,INFO)
            ! write(*,*) " Calling PLASMA_SPOSV"
            call PLASMA_SPOSV(local_UPLO,N,NRHS,A,LDA,B,LDB,INFO)
      end subroutine plasma_wrap_SPOSV

      subroutine plasma_wrap_SPOTRF(UPLO,N,A,LDA,INFO)
            use iso_c_binding
            use plasma
            implicit none
            integer, parameter :: wp = kind(0.0d0)
            integer, intent(in) :: LDA
            integer, intent(in) :: N
            integer, intent(out) :: INFO
            character, intent(in) :: UPLO
            real, intent(inout), target :: A(LDA,*)
            integer :: local_UPLO
            if(UPLO=='U' .or. UPLO=='u')then
               local_UPLO = PlasmaUpper
            else if(UPLO=='L' .or. UPLO=='l')then
               local_UPLO = PlasmaLower
            else
               local_UPLO = -1
            end if
            if (.not. plasma_initialized) call plasma_init(24,INFO)
            ! write(*,*) " Calling PLASMA_SPOTRF"
            call PLASMA_SPOTRF(local_UPLO,N,A,LDA,INFO)
      end subroutine plasma_wrap_SPOTRF

      subroutine plasma_wrap_SPOTRI(UPLO,N,A,LDA,INFO)
            use iso_c_binding
            use plasma
            implicit none
            integer, parameter :: wp = kind(0.0d0)
            integer, intent(in) :: LDA
            integer, intent(in) :: N
            integer, intent(out) :: INFO
            character, intent(in) :: UPLO
            real, intent(inout), target :: A(LDA,*)
            integer :: local_UPLO
            if(UPLO=='U' .or. UPLO=='u')then
               local_UPLO = PlasmaUpper
            else if(UPLO=='L' .or. UPLO=='l')then
               local_UPLO = PlasmaLower
            else
               local_UPLO = -1
            end if
            if (.not. plasma_initialized) call plasma_init(24,INFO)
            ! write(*,*) " Calling PLASMA_SPOTRI"
            call PLASMA_SPOTRI(local_UPLO,N,A,LDA,INFO)
      end subroutine plasma_wrap_SPOTRI

      subroutine plasma_wrap_SPOTRS(UPLO,N,NRHS,A,LDA,B,LDB,INFO)
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
            real, intent(inout), target :: A(LDA,*)
            real, intent(inout), target :: B(LDB,*)
            integer :: local_UPLO
            if(UPLO=='U' .or. UPLO=='u')then
               local_UPLO = PlasmaUpper
            else if(UPLO=='L' .or. UPLO=='l')then
               local_UPLO = PlasmaLower
            else
               local_UPLO = -1
            end if
            if (.not. plasma_initialized) call plasma_init(24,INFO)
            ! write(*,*) " Calling PLASMA_SPOTRS"
            call PLASMA_SPOTRS(local_UPLO,N,NRHS,A,LDA,B,LDB,INFO)
      end subroutine plasma_wrap_SPOTRS

      subroutine plasma_wrap_STRTRI(UPLO,DIAG,N,A,LDA,INFO)
            use iso_c_binding
            use plasma
            implicit none
            integer, parameter :: wp = kind(0.0d0)
            integer, intent(in) :: LDA
            integer, intent(in) :: N
            integer, intent(out) :: INFO
            character, intent(in) :: DIAG
            character, intent(in) :: UPLO
            real, intent(inout), target :: A(LDA,*)
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
            ! write(*,*) " Calling PLASMA_STRTRI"
            call PLASMA_STRTRI(local_UPLO,local_DIAG,N,A,LDA,INFO)
      end subroutine plasma_wrap_STRTRI

      subroutine plasma_wrap_SGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
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
            real, intent(in) :: ALPHA
            real, intent(in) :: BETA
            real, intent(inout), target :: A(LDA,*)
            real, intent(inout), target :: B(LDB,*)
            real, intent(inout), target :: C(LDC,*)
            integer :: local_TRANSA
            integer :: local_TRANSB
            integer :: local_ret
            if(TRANSA=='N' .or. TRANSA=='n')then
               local_TRANSA = PlasmaNoTrans
            else if(TRANSA=='T' .or. TRANSA=='t')then
               local_TRANSA = PlasmaTrans
            else if(TRANSA=='C' .or. TRANSA=='c')then
               local_TRANSA = PlasmaTrans
            else
               local_TRANSA = -1
            end if
            if(TRANSB=='N' .or. TRANSB=='n')then
               local_TRANSB = PlasmaNoTrans
            else if(TRANSB=='T' .or. TRANSB=='t')then
               local_TRANSB = PlasmaTrans
            else if(TRANSB=='C' .or. TRANSB=='c')then
               local_TRANSB = PlasmaTrans
            else
               local_TRANSB = -1
            end if
            if (.not. plasma_initialized) call plasma_init(24,local_ret)
            ! write(*,*) " Calling PLASMA_SGEMM"
            call PLASMA_SGEMM(local_TRANSA,local_TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC,local_ret)
      end subroutine plasma_wrap_SGEMM

#if defined(PRECISION_z) || defined(PRECISION_c)
      subroutine plasma_wrap_SSYMM(SIDE,UPLO,M,N,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
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
            real, intent(in) :: ALPHA
            real, intent(in) :: BETA
            real, intent(inout), target :: A(LDA,*)
            real, intent(inout), target :: B(LDB,*)
            real, intent(inout), target :: C(LDC,*)
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
            ! write(*,*) " Calling PLASMA_SSYMM"
            call PLASMA_SSYMM(local_SIDE,local_UPLO,M,N,ALPHA,A,LDA,B,LDB,BETA,C,LDC,local_ret)
      end subroutine plasma_wrap_SSYMM

      subroutine plasma_wrap_SSYR2K(UPLO,TRANS,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
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
            real, intent(in) :: ALPHA
            real, intent(inout), target :: A(LDA,*)
            real, intent(inout), target :: B(LDB,*)
            real, intent(inout), target :: C(LDC,*)
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
               local_TRANS = PlasmaTrans
            else
               local_TRANS = -1
            end if
            if (.not. plasma_initialized) call plasma_init(24,local_ret)
            ! write(*,*) " Calling PLASMA_SSYR2K"
            call PLASMA_SSYR2K(local_UPLO,local_TRANS,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC,local_ret)
      end subroutine plasma_wrap_SSYR2K

      subroutine plasma_wrap_SSYRK(UPLO,TRANS,N,K,ALPHA,A,LDA,BETA,C,LDC)
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
            real, intent(inout), target :: A(LDA,*)
            real, intent(inout), target :: C(LDC,*)
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
               local_TRANS = PlasmaTrans
            else
               local_TRANS = -1
            end if
            if (.not. plasma_initialized) call plasma_init(24,local_ret)
            ! write(*,*) " Calling PLASMA_SSYRK"
            call PLASMA_SSYRK(local_UPLO,local_TRANS,N,K,ALPHA,A,LDA,BETA,C,LDC,local_ret)
      end subroutine plasma_wrap_SSYRK
#endif

      subroutine plasma_wrap_SSYMM(SIDE,UPLO,M,N,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
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
            real, intent(in) :: ALPHA
            real, intent(in) :: BETA
            real, intent(inout), target :: A(LDA,*)
            real, intent(inout), target :: B(LDB,*)
            real, intent(inout), target :: C(LDC,*)
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
            ! write(*,*) " Calling PLASMA_SSYMM"
            call PLASMA_SSYMM(local_SIDE,local_UPLO,M,N,ALPHA,A,LDA,B,LDB,BETA,C,LDC,local_ret)
      end subroutine plasma_wrap_SSYMM

      subroutine plasma_wrap_SSYR2K(UPLO,TRANS,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
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
            real, intent(in) :: ALPHA
            real, intent(in) :: BETA
            real, intent(inout), target :: A(LDA,*)
            real, intent(inout), target :: B(LDB,*)
            real, intent(inout), target :: C(LDC,*)
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
               local_TRANS = PlasmaTrans
            else
               local_TRANS = -1
            end if
            if (.not. plasma_initialized) call plasma_init(24,local_ret)
            ! write(*,*) " Calling PLASMA_SSYR2K"
            call PLASMA_SSYR2K(local_UPLO,local_TRANS,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC,local_ret)
      end subroutine plasma_wrap_SSYR2K

      subroutine plasma_wrap_SSYRK(UPLO,TRANS,N,K,ALPHA,A,LDA,BETA,C,LDC)
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
            real, intent(in) :: ALPHA
            real, intent(in) :: BETA
            real, intent(inout), target :: A(LDA,*)
            real, intent(inout), target :: C(LDC,*)
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
               local_TRANS = PlasmaTrans
            else
               local_TRANS = -1
            end if
            if (.not. plasma_initialized) call plasma_init(24,local_ret)
            ! write(*,*) " Calling PLASMA_SSYRK"
            call PLASMA_SSYRK(local_UPLO,local_TRANS,N,K,ALPHA,A,LDA,BETA,C,LDC,local_ret)
      end subroutine plasma_wrap_SSYRK

      subroutine plasma_wrap_STRMM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
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
            real, intent(in) :: ALPHA
            real, intent(inout), target :: A(LDA,*)
            real, intent(inout), target :: B(LDB,*)
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
               local_TRANSA = PlasmaTrans
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
            ! write(*,*) " Calling PLASMA_STRMM"
            call PLASMA_STRMM(local_SIDE,local_UPLO,local_TRANSA,local_DIAG,M,N,ALPHA,A,LDA,B,LDB,local_ret)
      end subroutine plasma_wrap_STRMM

      subroutine plasma_wrap_STRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
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
            real, intent(in) :: ALPHA
            real, intent(inout), target :: A(LDA,*)
            real, intent(inout), target :: B(LDB,*)
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
               local_TRANSA = PlasmaTrans
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
            ! write(*,*) " Calling PLASMA_STRSM"
            call PLASMA_STRSM(local_SIDE,local_UPLO,local_TRANSA,local_DIAG,M,N,ALPHA,A,LDA,B,LDB,local_ret)
      end subroutine plasma_wrap_STRSM

      subroutine plasma_wrap_SLACPY(UPLO,M,N,A,LDA,B,LDB)
            use iso_c_binding
            use plasma
            implicit none
            integer, parameter :: wp = kind(0.0d0)
            integer, intent(in) :: LDA
            integer, intent(in) :: LDB
            integer, intent(in) :: M
            integer, intent(in) :: N
            character, intent(in) :: UPLO
            real, intent(inout), target :: A(LDA,*)
            real, intent(out), target :: B(LDB,*)
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
            ! write(*,*) " Calling PLASMA_SLACPY"
            call PLASMA_SLACPY(local_UPLO,M,N,A,LDA,B,LDB,local_ret)
      end subroutine plasma_wrap_SLACPY

      subroutine plasma_wrap_SLASET(UPLO,M,N,ALPHA,BETA,A,LDA)
            use iso_c_binding
            use plasma
            implicit none
            integer, parameter :: wp = kind(0.0d0)
            integer, intent(in) :: LDA
            integer, intent(in) :: M
            integer, intent(in) :: N
            character, intent(in) :: UPLO
            real, intent(in) :: ALPHA
            real, intent(in) :: BETA
            real, intent(inout), target :: A(LDA,*)
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
            ! write(*,*) " Calling PLASMA_SLASET"
            call PLASMA_SLASET(local_UPLO,M,N,ALPHA,BETA,A,LDA,local_ret)
      end subroutine plasma_wrap_SLASET

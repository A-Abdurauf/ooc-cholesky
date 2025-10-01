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
! @file plasma_dsf90.F90
!
!  PLASMA Fortran 90 interfaces using Fortran 2003 ISO C bindings
!  PLASMA is a software package provided by Univ. of Tennessee,
!  Univ. of California Berkeley and Univ. of Colorado Denver
!
! @version 2.8.0
! @author Numerical Algorithms Group
! @author Mathieu Faverge
! @date 2011-12-15
! @generated ds Fri Apr  1 11:03:02 2016
!
module plasma_ds

      interface
         function PLASMA_dsgesv_c(N,NRHS,A,LDA,IPIV,B,LDB,X,LDX,ITER) &
          & bind(c, name='PLASMA_dsgesv')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_dsgesv_c
            integer(kind=c_int), value :: N
            integer(kind=c_int), value :: NRHS
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: IPIV
            type(c_ptr), value :: B
            integer(kind=c_int), value :: LDB
            type(c_ptr), value :: X
            integer(kind=c_int), value :: LDX
            type(c_ptr), value :: ITER
         end function PLASMA_dsgesv_c
      end interface

      interface
         function PLASMA_dsposv_c(uplo,N,NRHS,A,LDA,B,LDB,X,LDX,ITER) &
          & bind(c, name='PLASMA_dsposv')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_dsposv_c
            integer(kind=c_int), value :: uplo
            integer(kind=c_int), value :: N
            integer(kind=c_int), value :: NRHS
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: B
            integer(kind=c_int), value :: LDB
            type(c_ptr), value :: X
            integer(kind=c_int), value :: LDX
            type(c_ptr), value :: ITER
         end function PLASMA_dsposv_c
      end interface

      interface
         function PLASMA_dsungesv_c(trans,N,NRHS,A,LDA,B,LDB,X,LDX,ITER) &
          & bind(c, name='PLASMA_dsungesv')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_dsungesv_c
            integer(kind=c_int), value :: trans
            integer(kind=c_int), value :: N
            integer(kind=c_int), value :: NRHS
            type(c_ptr), value :: A
            integer(kind=c_int), value :: LDA
            type(c_ptr), value :: B
            integer(kind=c_int), value :: LDB
            type(c_ptr), value :: X
            integer(kind=c_int), value :: LDX
            type(c_ptr), value :: ITER
         end function PLASMA_dsungesv_c
      end interface

      interface
         function PLASMA_dsgesv_Tile_c(A,IPIV,B,X,ITER) &
          & bind(c, name='PLASMA_dsgesv_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_dsgesv_Tile_c
            type(c_ptr), value :: A
            type(c_ptr), value :: IPIV
            type(c_ptr), value :: B
            type(c_ptr), value :: X
            type(c_ptr), value :: ITER
         end function PLASMA_dsgesv_Tile_c
      end interface

      interface
         function PLASMA_dsposv_Tile_c(uplo,A,B,X,ITER) &
          & bind(c, name='PLASMA_dsposv_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_dsposv_Tile_c
            integer(kind=c_int), value :: uplo
            type(c_ptr), value :: A
            type(c_ptr), value :: B
            type(c_ptr), value :: X
            type(c_ptr), value :: ITER
         end function PLASMA_dsposv_Tile_c
      end interface

      interface
         function PLASMA_dsungesv_Tile_c(trans,A,T,B,X,ITER) &
          & bind(c, name='PLASMA_dsungesv_Tile')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_dsungesv_Tile_c
            integer(kind=c_int), value :: trans
            type(c_ptr), value :: A
            type(c_ptr), value :: T
            type(c_ptr), value :: B
            type(c_ptr), value :: X
            type(c_ptr), value :: ITER
         end function PLASMA_dsungesv_Tile_c
      end interface

      interface
         function PLASMA_dsgesv_Tile_Async_c(A,IPIV,B,X,ITER,sequence,request) &
          & bind(c, name='PLASMA_dsgesv_Tile_Async')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_dsgesv_Tile_Async_c
            type(c_ptr), value :: A
            type(c_ptr), value :: IPIV
            type(c_ptr), value :: B
            type(c_ptr), value :: X
            type(c_ptr), value :: ITER
            type(c_ptr), value :: sequence
            type(c_ptr), value :: request
         end function PLASMA_dsgesv_Tile_Async_c
      end interface

      interface
         function PLASMA_dsposv_Tile_Async_c(uplo,A,B,X,ITER,sequence,request) &
          & bind(c, name='PLASMA_dsposv_Tile_Async')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_dsposv_Tile_Async_c
            integer(kind=c_int), value :: uplo
            type(c_ptr), value :: A
            type(c_ptr), value :: B
            type(c_ptr), value :: X
            type(c_ptr), value :: ITER
            type(c_ptr), value :: sequence
            type(c_ptr), value :: request
         end function PLASMA_dsposv_Tile_Async_c
      end interface

      interface
         function PLASMA_dsungesv_Tile_Async_c(trans,A,T,B,X,ITER,sequence,request) &
          & bind(c, name='PLASMA_dsungesv_Tile_Async')
            use iso_c_binding
            implicit none
            integer(kind=c_int) :: PLASMA_dsungesv_Tile_Async_c
            integer(kind=c_int), value :: trans
            type(c_ptr), value :: A
            type(c_ptr), value :: T
            type(c_ptr), value :: B
            type(c_ptr), value :: X
            type(c_ptr), value :: ITER
            type(c_ptr), value :: sequence
            type(c_ptr), value :: request
         end function PLASMA_dsungesv_Tile_Async_c
      end interface

  contains

      subroutine PLASMA_dsgesv(N,NRHS,A,LDA,IPIV,B,LDB,X,LDX,ITER,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(out), target :: ITER
         integer(kind=c_int), intent(out), target :: IPIV(*)
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDB
         integer(kind=c_int), intent(in) :: LDX
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: NRHS
         real(kind=c_double), intent(in), target :: A(LDA,*)
         real(kind=c_double), intent(in), target :: B(LDB,*)
         real(kind=c_double), intent(out), target :: X(LDX,*)
         info = PLASMA_dsgesv_c(N,NRHS,c_loc(A),LDA,c_loc(IPIV),c_loc(B),LDB,c_loc(X),LDX,c_loc(ITER))
      end subroutine PLASMA_dsgesv

      subroutine PLASMA_dsposv(uplo,N,NRHS,A,LDA,B,LDB,X,LDX,ITER,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(out), target :: ITER
         integer(kind=c_int), intent(in) :: uplo
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDB
         integer(kind=c_int), intent(in) :: LDX
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: NRHS
         real(kind=c_double), intent(in), target :: A(LDA,*)
         real(kind=c_double), intent(in), target :: B(LDB,*)
         real(kind=c_double), intent(out), target :: X(LDX,*)
         info = PLASMA_dsposv_c(uplo,N,NRHS,c_loc(A),LDA,c_loc(B),LDB,c_loc(X),LDX,c_loc(ITER))
      end subroutine PLASMA_dsposv

      subroutine PLASMA_dsungesv(trans,N,NRHS,A,LDA,B,LDB,X,LDX,ITER,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(out), target :: ITER
         integer(kind=c_int), intent(in) :: trans
         integer(kind=c_int), intent(in) :: LDA
         integer(kind=c_int), intent(in) :: LDB
         integer(kind=c_int), intent(in) :: LDX
         integer(kind=c_int), intent(in) :: N
         integer(kind=c_int), intent(in) :: NRHS
         real(kind=c_double), intent(in), target :: A(LDA,*)
         real(kind=c_double), intent(in), target :: B(LDB,*)
         real(kind=c_double), intent(out), target :: X(LDX,*)
         info = PLASMA_dsungesv_c(trans,N,NRHS,c_loc(A),LDA,c_loc(B),LDB,c_loc(X),LDX,c_loc(ITER))
      end subroutine PLASMA_dsungesv

      subroutine PLASMA_dsgesv_Tile(A,IPIV,B,X,ITER,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(out), target :: IPIV(*)
          integer(kind=c_int), intent(out), target :: ITER
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: X ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_dsgesv_Tile_c(A,c_loc(IPIV),B,X,c_loc(ITER))
      end subroutine PLASMA_dsgesv_Tile

      subroutine PLASMA_dsposv_Tile(uplo,A,B,X,ITER,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(out), target :: ITER
         integer(kind=c_int), intent(in) :: uplo
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: X ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_dsposv_Tile_c(uplo,A,B,X,c_loc(ITER))
      end subroutine PLASMA_dsposv_Tile

      subroutine PLASMA_dsungesv_Tile(trans,A,T,B,X,ITER,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(out), target :: ITER
         integer(kind=c_int), intent(in) :: trans
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: X ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_dsungesv_Tile_c(trans,A,T,B,X,c_loc(ITER))
      end subroutine PLASMA_dsungesv_Tile

      subroutine PLASMA_dsgesv_Tile_Async(A,IPIV,B,X,ITER,sequence,request,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(out), target :: IPIV(*)
         integer(kind=c_int), intent(out), target :: ITER
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: X ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: request ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: sequence ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_dsgesv_Tile_Async_c(A,c_loc(IPIV),B,X,c_loc(ITER),sequence,request)
      end subroutine PLASMA_dsgesv_Tile_Async

      subroutine PLASMA_dsposv_Tile_Async(uplo,A,B,X,ITER,sequence,request,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(out), target :: ITER
         integer(kind=c_int), intent(in) :: uplo
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: X ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: request ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: sequence ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_dsposv_Tile_Async_c(uplo,A,B,X,c_loc(ITER),sequence,request)
      end subroutine PLASMA_dsposv_Tile_Async

      subroutine PLASMA_dsungesv_Tile_Async(trans,A,T,B,X,ITER,sequence,request,info)
         use iso_c_binding
         implicit none
         integer(kind=c_int), intent(out) :: info
         integer(kind=c_int), intent(out), target :: ITER
         integer(kind=c_int), intent(in) :: trans
         type(c_ptr), value :: A ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: B ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: T ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: X ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: request ! Arg managed by PLASMA: opaque to Fortran
         type(c_ptr), value :: sequence ! Arg managed by PLASMA: opaque to Fortran
         info = PLASMA_dsungesv_Tile_Async_c(trans,A,T,B,X,c_loc(ITER),sequence,request)
      end subroutine PLASMA_dsungesv_Tile_Async

end module plasma_ds

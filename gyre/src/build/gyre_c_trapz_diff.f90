!fpx3_header(0.13_3a)
!
!dependencies
!   dir: ~/Documents/MESA_SS/gyre/src/build 
!   sources: -
!   includes: ../extern/core/core.inc ../diff/gyre_trapz_diff.inc
!   uses: ISO_FORTRAN_ENV gyre_diff gyre_state gyre_linalg gyre_eqns gyre_point core_kinds core_linalg gyre_ext
!   provides: gyre_c_trapz_diff
!end dependencies
!
!end fpx3_header
! Module   : gyre_c_trapz_diff
! Purpose  : difference equations (quasi-trapezoidal, complex)
!
! Copyright 2013-2015 Rich Townsend
!
! This file is part of GYRE. GYRE is free software: you can
! redistribute it and/or modify it under the terms of the GNU General
! Public License as published by the Free Software Foundation, version 3.
!
! GYRE is distributed in the hope that it will be useful, but WITHOUT
! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
! or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
! License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.

! Incfile  : gyre_trapz_diff
! Purpose  : difference equations (quasi-trapezoidal, template)
!
! Copyright 2013-2017 Rich Townsend
!
! This file is part of GYRE. GYRE is free software: you can
! redistribute it and/or modify it under the terms of the GNU General
! Public License as published by the Free Software Foundation, version 3.
!
! GYRE is distributed in the hope that it will be useful, but WITHOUT
! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
! or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
! License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.

! Incfile  : core
! Purpose  : fpx3 macros

!****

!****

!****

!****

!****

!****

!****

!****

!****

module gyre_c_trapz_diff

  ! Uses

  use core_kinds
  use core_linalg

  use gyre_diff
  use gyre_eqns
  use gyre_ext
  use gyre_linalg
  use gyre_point
  use gyre_state

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type, extends (c_diff_t) :: c_trapz_diff_t
     private
     class(c_eqns_t), allocatable :: eq
     real(WP), allocatable           :: w(:)
     real(WP)                        :: dx
     real(WP)                        :: x_m
   contains
     private
     procedure, public :: build
  end type c_trapz_diff_t

  ! Interfaces

  interface c_trapz_diff_t
     module procedure c_trapz_diff_t_
  end interface c_trapz_diff_t

  ! Access specifiers

  private

  public :: c_trapz_diff_t

contains

  function c_trapz_diff_t_ (eq, pt_a, pt_b, w) result (df)

    class(c_eqns_t), intent(in) :: eq
    type(point_t), intent(in)      :: pt_a
    type(point_t), intent(in)      :: pt_b
    real(WP), intent(in)           :: w(:)
    type(c_trapz_diff_t)        :: df

    ! Construct the trapz_diff_t

    df%w = w

    df%dx = pt_b%x - pt_a%x
    df%x_m = 0.5_WP*(pt_a%x + pt_b%x)

    allocate(df%eq, SOURCE=eq)

    call df%eq%stencil([pt_a,pt_b])

    df%n_e = eq%n_e

    ! Finish

    return

  end function c_trapz_diff_t_

  !****

  subroutine build (this, st, E_l, E_r, scl)

    class(c_trapz_diff_t), intent(in) :: this
    class(c_state_t), intent(in)      :: st
    complex(WP), intent(out)               :: E_l(:,:)
    complex(WP), intent(out)               :: E_r(:,:)
    type(c_ext_t), intent(out)        :: scl

    complex(WP) :: xA(this%n_e,this%n_e,2)

    ! Evaluate the RHS matrix

    xA(:,:,1) = this%eq%xA(1, st)
    xA(:,:,2) = this%eq%xA(2, st)

    ! Build the difference equations

    E_l = -this%x_m*identity_matrix(this%n_e) - this%dx*MATMUL(diagonal_matrix(1._WP-this%w), xA(:,:,1))
    E_r =  this%x_m*identity_matrix(this%n_e) - this%dx*MATMUL(diagonal_matrix(      this%w), xA(:,:,2))

    scl = c_ext_t(1._WP)

    ! Finish

  end subroutine build

end module gyre_c_trapz_diff


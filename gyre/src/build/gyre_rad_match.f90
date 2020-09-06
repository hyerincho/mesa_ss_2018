!fpx3_header(0.13_3a)
!
!dependencies
!   dir: ~/Documents/MESA_SS/gyre/src/build 
!   sources: -
!   includes: ../extern/core/core.inc
!   uses: gyre_mode_par gyre_model_util gyre_state gyre_ext gyre_rad_trans gyre_model gyre_diff core_kinds gyre_osc_par ISO_FORTRAN_ENV gyre_context gyre_point
!   provides: gyre_rad_match
!end dependencies
!
!end fpx3_header
! Module   : gyre_rad_match
! Purpose  : adiabatic radial match conditions
!
! Copyright 2016-2017 Rich Townsend
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

module gyre_rad_match

  ! Uses

  use core_kinds

  use gyre_context
  use gyre_diff
  use gyre_ext
  use gyre_model
  use gyre_model_util
  use gyre_mode_par
  use gyre_osc_par
  use gyre_point
  use gyre_rad_trans
  use gyre_state

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Parameter definitions

  integer, parameter :: J_U = 1

  integer, parameter :: J_LAST = J_U

  ! Derived-type definitions

  type, extends (r_diff_t) :: rad_match_t
     private
     type(context_t), pointer :: cx => null()
     type(rad_trans_t)        :: tr
     real(WP), allocatable    :: coeff(:,:)
   contains
     private
     procedure         :: stencil_
     procedure, public :: build
  end type rad_match_t

  ! Interfaces

  interface rad_match_t
     module procedure rad_match_t_
  end interface rad_match_t

  ! Access specifiers

  private
  public :: rad_match_t

contains

  function rad_match_t_ (cx, pt_i, pt_a, pt_b, md_p, os_p) result (mt)

    type(context_t), pointer, intent(in) :: cx
    type(point_t), intent(in)            :: pt_i
    type(point_t), intent(in)            :: pt_a
    type(point_t), intent(in)            :: pt_b
    type(mode_par_t), intent(in)         :: md_p
    type(osc_par_t), intent(in)          :: os_p
    type(rad_match_t)                    :: mt

    ! Construct the rad_match_t

    mt%cx => cx

    mt%tr = rad_trans_t(cx, pt_i, md_p, os_p)

    call mt%stencil_(pt_a, pt_b)

    mt%n_e = 2

    ! Finish

    return

  end function rad_match_t_

  !****

  subroutine stencil_ (this, pt_a, pt_b)

    class(rad_match_t), intent(inout) :: this
    type(point_t), intent(in)         :: pt_a
    type(point_t), intent(in)         :: pt_b

    ! Calculate coefficients at the stencil points

    associate (ml => this%cx%ml)

      call check_model(ml, [I_U])

      allocate(this%coeff(2,J_LAST))

      this%coeff(1,J_U) = ml%coeff(I_U, pt_a)
      this%coeff(2,J_U) = ml%coeff(I_U, pt_b)

    end associate

    ! Set up stencil for the tr component

    call this%tr%stencil([pt_a,pt_b])

    ! Finish

    return

  end subroutine stencil_

  !****

  subroutine build (this, st, E_l, E_r, scl)

    class(rad_match_t), intent(in) :: this
    class(r_state_t), intent(in)   :: st
    real(WP), intent(out)          :: E_l(:,:)
    real(WP), intent(out)          :: E_r(:,:)
    type(r_ext_t), intent(out)     :: scl

    ! Build the difference equations

    ! Calculate coefficients

    associate( &
      U_l => this%coeff(1,J_U), &
      U_r => this%coeff(2,J_U))

      ! Evaluate the match conditions (y_1 continuous, y_2 not)

      E_l(1,1) = -1._WP
      E_l(1,2) = 0._WP

      E_l(2,1) = U_l
      E_l(2,2) = -U_l

      !

      E_r(1,1) = 1._WP
      E_r(1,2) = 0._WP

      E_r(2,1) = -U_r
      E_r(2,2) = U_r

      scl = r_ext_t(1._WP)

    end associate

    ! Apply the variables transformation

    call this%tr%trans_cond(E_l, 1, st)
    call this%tr%trans_cond(E_r, 2, st)

    ! Finish

    return

  end subroutine build

end module gyre_rad_match

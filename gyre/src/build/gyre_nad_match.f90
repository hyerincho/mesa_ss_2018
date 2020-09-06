!fpx3_header(0.13_3a)
!
!dependencies
!   dir: ~/Documents/MESA_SS/gyre/src/build 
!   sources: -
!   includes: ../extern/core/core.inc
!   uses: gyre_model gyre_point gyre_osc_par gyre_ext gyre_state gyre_nad_trans gyre_diff core_kinds gyre_model_util gyre_context ISO_FORTRAN_ENV gyre_mode_par
!   provides: gyre_nad_match
!end dependencies
!
!end fpx3_header
! Module   : gyre_nad_match
! Purpose  : nonadiabatic match conditions
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

module gyre_nad_match

  ! Uses

  use core_kinds

  use gyre_context
  use gyre_diff
  use gyre_ext
  use gyre_model
  use gyre_model_util
  use gyre_mode_par
  use gyre_nad_trans
  use gyre_osc_par
  use gyre_point
  use gyre_state

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Parameter definitions

  integer, parameter :: J_V = 1
  integer, parameter :: J_U = 2
  integer, parameter :: J_NABLA_AD = 3

  integer, parameter :: J_LAST = J_NABLA_AD

  ! Derived-type definitions

  type, extends (c_diff_t) :: nad_match_t
     private
     type(context_t), pointer :: cx => null()
     type(nad_trans_t)        :: tr
     real(WP), allocatable    :: coeff(:,:)
   contains
     private
     procedure         :: stencil_
     procedure, public :: build
  end type nad_match_t

  ! Interfaces

  interface nad_match_t
     module procedure nad_match_t_
  end interface nad_match_t

  ! Access specifiers

  private
  public :: nad_match_t

contains

  function nad_match_t_ (cx, pt_a, pt_b, md_p, os_p) result (mt)

    type(context_t), pointer, intent(in) :: cx
    type(point_t), intent(in)            :: pt_a
    type(point_t), intent(in)            :: pt_b
    type(mode_par_t), intent(in)         :: md_p
    type(osc_par_t), intent(in)          :: os_p
    type(nad_match_t)                    :: mt

    ! Construct the nad_match_t

    mt%cx => cx

    mt%tr = nad_trans_t(cx, md_p, os_p)

    call mt%stencil_(pt_a, pt_b)

    mt%n_e = 6

    ! Finish

    return

  end function nad_match_t_

  !****

  subroutine stencil_ (this, pt_a, pt_b)

    class(nad_match_t), intent(inout) :: this
    type(point_t), intent(in)         :: pt_a
    type(point_t), intent(in)         :: pt_b

    ! Calculate coefficients at the stencil points

    associate (ml => this%cx%ml)

      call check_model(ml, [I_V_2,I_U,I_NABLA_AD])

      allocate(this%coeff(2,J_LAST))

      this%coeff(1,J_V) = ml%coeff(I_V_2, pt_a)*pt_a%x**2
      this%coeff(2,J_V) = ml%coeff(I_V_2, pt_b)*pt_b%x**2

      this%coeff(1,J_U) = ml%coeff(I_U, pt_a)
      this%coeff(2,J_U) = ml%coeff(I_U, pt_b)

      this%coeff(1,J_NABLA_AD) = ml%coeff(I_NABLA_AD, pt_a)
      this%coeff(2,J_NABLA_AD) = ml%coeff(I_NABLA_AD, pt_b)

    end associate

    ! Set up stencil for the tr component

    call this%tr%stencil([pt_a,pt_b])

    ! Finish

    return

  end subroutine stencil_

  !****

  subroutine build (this, st, E_l, E_r, scl)

    class(nad_match_t), intent(in) :: this
    class(c_state_t), intent(in)   :: st
    complex(WP), intent(out)       :: E_l(:,:)
    complex(WP), intent(out)       :: E_r(:,:)
    type(c_ext_t), intent(out)     :: scl

    ! Build the difference equations

    associate( &
         V_l => this%coeff(1,J_V), &
         V_r => this%coeff(2,J_V), &
         U_l => this%coeff(1,J_U), &
         U_r => this%coeff(2,J_U), &
         nabla_ad_l => this%coeff(1,J_NABLA_AD), &
         nabla_ad_r => this%coeff(2,J_NABLA_AD))

      ! Evaluate the match conditions (y_1, y_3, y_6 continuous, y_2,
      ! y_4, y_5 not)

      E_l(1,1) = -1._WP
      E_l(1,2) = 0._WP
      E_l(1,3) = 0._WP
      E_l(1,4) = 0._WP
      E_l(1,5) = 0._WP
      E_l(1,6) = 0._WP

      E_l(2,1) = U_l
      E_l(2,2) = -U_l
      E_l(2,3) = 0._WP
      E_l(2,4) = 0._WP
      E_l(2,5) = 0._WP
      E_l(2,6) = 0._WP

      E_l(3,1) = 0._WP
      E_l(3,2) = 0._WP
      E_l(3,3) = -1._WP
      E_l(3,4) = 0._WP
      E_l(3,5) = 0._WP
      E_l(3,6) = 0._WP

      E_l(4,1) = -U_l
      E_l(4,2) = 0._WP
      E_l(4,3) = 0._WP
      E_l(4,4) = -1._WP
      E_l(4,5) = 0._WP
      E_l(4,6) = 0._WP

      E_l(5,1) = V_l*nabla_ad_l
      E_l(5,2) = -V_l*nabla_ad_l
      E_l(5,3) = 0._WP
      E_l(5,4) = 0._WP
      E_l(5,5) = -1._WP
      E_l(5,6) = 0._WP

      E_l(6,1) = 0._WP
      E_l(6,2) = 0._WP
      E_l(6,3) = 0._WP
      E_l(6,4) = 0._WP
      E_l(6,5) = 0._WP
      E_l(6,6) = -1._WP

      !

      E_r(1,1) = 1._WP
      E_r(1,2) = 0._WP
      E_r(1,3) = 0._WP
      E_r(1,4) = 0._WP
      E_r(1,5) = 0._WP
      E_r(1,6) = 0._WP

      E_r(2,1) = -U_r
      E_r(2,2) = U_r
      E_r(2,3) = 0._WP
      E_r(2,4) = 0._WP
      E_r(2,5) = 0._WP
      E_r(2,6) = 0._WP

      E_r(3,1) = 0._WP
      E_r(3,2) = 0._WP
      E_r(3,3) = 1._WP
      E_r(3,4) = 0._WP
      E_r(3,5) = 0._WP
      E_r(3,6) = 0._WP

      E_r(4,1) = U_r
      E_r(4,2) = 0._WP
      E_r(4,3) = 0._WP
      E_r(4,4) = 1._WP
      E_r(4,5) = 0._WP
      E_r(4,6) = 0._WP

      E_r(5,1) = -V_r*nabla_ad_r
      E_r(5,2) = V_r*nabla_ad_r
      E_r(5,3) = 0._WP
      E_r(5,4) = 0._WP
      E_r(5,5) = 1._WP
      E_r(5,6) = 0._WP

      E_r(6,1) = 0._WP
      E_r(6,2) = 0._WP
      E_r(6,3) = 0._WP
      E_r(6,4) = 0._WP
      E_r(6,5) = 0._WP
      E_r(6,6) = 1._WP

      scl = c_ext_t(1._WP)

    end associate

    ! Apply the variables transformation

    call this%tr%trans_cond(E_l, 1, st)
    call this%tr%trans_cond(E_r, 2, st)

    ! Finish

    return

  end subroutine build

end module gyre_nad_match
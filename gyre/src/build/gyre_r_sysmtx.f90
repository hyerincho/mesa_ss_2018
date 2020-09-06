!fpx3_header(0.13_3a)
!
!dependencies
!   dir: ~/Documents/MESA_SS/gyre/src/build 
!   sources: -
!   includes: ../matrix/gyre_sysmtx.inc ../extern/core/core.inc
!   uses: core_kinds gyre_ext
!   provides: gyre_r_sysmtx
!end dependencies
!
!end fpx3_header
! Module   : gyre_r_sysmtx
! Purpose  : system matrix (real)
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

! Incfile  : gyre_sysmtx
! Purpose  : system matrix (template)
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

module gyre_r_sysmtx

  ! Uses

  use core_kinds

  use gyre_ext

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type, abstract :: r_sysmtx_t
     integer :: n   ! Number of equation blocks
     integer :: n_e ! Number of equations per block
     integer :: n_i ! Number of inner boundary conditions
     integer :: n_o ! Number of outer boundary conditions
   contains
     procedure(set_B), deferred          :: set_B_i
     procedure(set_B), deferred          :: set_B_o
     procedure(set_E), deferred          :: set_E
     procedure(factor), deferred         :: factor
     procedure(det), deferred            :: det
     procedure(soln_vec_hom), deferred   :: soln_vec_hom
     procedure(soln_vec_inhom), deferred :: soln_vec_inhom
     procedure(resd_vec), deferred       :: resd_vec
  end type r_sysmtx_t

  ! Interfaces

  abstract interface

     subroutine set_B (this, B, scl)
       use core_kinds
       use gyre_ext
       import r_sysmtx_t
       class(r_sysmtx_t), intent(inout) :: this
       real(WP), intent(in)               :: B(:,:)
       real(WP), intent(in)               :: scl(:)
     end subroutine set_B

     subroutine set_E (this, k, E_l, E_r, scl)
       use core_kinds
       use gyre_ext
       import r_sysmtx_t
       class(r_sysmtx_t), intent(inout) :: this
       integer, intent(in)                 :: k
       real(WP), intent(in)               :: E_l(:,:)
       real(WP), intent(in)               :: E_r(:,:)
       type(r_ext_t), intent(in)        :: scl
     end subroutine set_E

     subroutine factor (this)
       import r_sysmtx_t
       class(r_sysmtx_t), intent(inout) :: this
     end subroutine factor

     function det (this)
       use gyre_ext
       import r_sysmtx_t
       class(r_sysmtx_t), intent(in) :: this
       type(r_ext_t)                 :: det
     end function det

     function soln_vec_hom (this) result (v)
       use core_kinds
       use gyre_ext
       import r_sysmtx_t
       class(r_sysmtx_t), intent(in) :: this
       real(WP)                        :: v(this%n_e*(this%n+1))
     end function soln_vec_hom

     function soln_vec_inhom (this, w_i, w_o) result (v)
       use core_kinds
       use gyre_ext
       import r_sysmtx_t
       class(r_sysmtx_t), intent(in) :: this
       real(WP), intent(in)            :: w_i(:)
       real(WP), intent(in)            :: w_o(:)
       real(WP)                        :: v(this%n_e*(this%n+1))
     end function soln_vec_inhom

     function resd_vec (this, v) result (w)
       use core_kinds
       use gyre_ext
       import r_sysmtx_t
       class(r_sysmtx_t), intent(in) :: this
       real(WP), intent(in)            :: v(:)
       real(WP)                        :: w(this%n_e*(this%n+1))
     end function resd_vec

  end interface

  ! Access specifiers

  private

  public :: r_sysmtx_t

end module gyre_r_sysmtx


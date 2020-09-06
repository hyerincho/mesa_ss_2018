!fpx3_header(0.13_3a)
!
!dependencies
!   dir: ~/Documents/MESA_SS/gyre/src/build 
!   sources: -
!   includes: ../diff/gyre_eqns.inc
!   uses: gyre_state core_kinds gyre_point
!   provides: gyre_c_eqns
!end dependencies
!
!end fpx3_header
! Module   : gyre_c_eqns
! Purpose  : differential equations evaluation (complex)
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

! Incfile  : gyre_eqns
! Purpose  : differential equations evaluation (template)
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

module gyre_c_eqns

  ! Uses

  use core_kinds

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type, abstract :: c_eqns_t
     integer :: n_e
   contains
     procedure(stencil), deferred :: stencil
     procedure(A), deferred       :: A
     procedure(xA), deferred      :: xA
  end type c_eqns_t

  ! Interfaces

  abstract interface

     subroutine stencil (this, pt)
       use gyre_point
       import c_eqns_t
       class(c_eqns_t), intent(inout) :: this
       type(point_t), intent(in)         :: pt(:)
     end subroutine stencil

     function A (this, i, st)
       use core_kinds
       use gyre_state
       import c_eqns_t
       class(c_eqns_t), intent(in)  :: this
       integer, intent(in)             :: i
       class(c_state_t), intent(in) :: st
       complex(WP)                       :: A(this%n_e,this%n_e)
     end function A

     function xA (this, i, st)
       use core_kinds
       use gyre_state
       import c_eqns_t
       class(c_eqns_t), intent(in)  :: this
       integer, intent(in)             :: i
       class(c_state_t), intent(in) :: st
       complex(WP)                       :: xA(this%n_e,this%n_e)
     end function xA

  end interface

  ! Access specifiers

  private

  public :: c_eqns_t

end module gyre_c_eqns


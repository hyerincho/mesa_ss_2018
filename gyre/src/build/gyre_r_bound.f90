!fpx3_header(0.13_3a)
!
!dependencies
!   dir: ~/Documents/MESA_SS/gyre/src/build 
!   sources: -
!   includes: ../extern/core/core.inc ../bvp/gyre_bound.inc
!   uses: core_kinds gyre_state
!   provides: gyre_r_bound
!end dependencies
!
!end fpx3_header
! Module   : gyre_r_bound
! Purpose  : boundary conditions (real)
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

! Incfile  : gyre_bound
! Purpose  : boundary conditions (template)
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

module gyre_r_bound

  ! Uses

  use core_kinds

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type, abstract :: r_bound_t
     integer :: n_i
     integer :: n_o
     integer :: n_e
   contains
     procedure(build), deferred :: build_i
     procedure(build), deferred :: build_o
   end type r_bound_t

  ! Interfaces

  abstract interface

     subroutine build (this, st, B, scl)
       use core_kinds
       use gyre_state
       import r_bound_t
       class(r_bound_t), intent(in) :: this
       class(r_state_t), intent(in) :: st
       real(WP), intent(out)          :: B(:,:)
       real(WP), intent(out)          :: scl(:)
     end subroutine build

  end interface

  ! Access specifiers

  private

  public :: r_bound_t

end module gyre_r_bound


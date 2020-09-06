!fpx3_header(0.13_3a)
!
!dependencies
!   dir: ~/Documents/MESA_SS/gyre/src/build 
!   sources: -
!   includes: ../extern/core/core.inc ../bvp/gyre_bound.inc
!   uses: gyre_state core_kinds
!   provides: gyre_c_bound
!end dependencies
!
!end fpx3_header
! Module   : gyre_c_bound
! Purpose  : boundary conditions (complex)
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

module gyre_c_bound

  ! Uses

  use core_kinds

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type, abstract :: c_bound_t
     integer :: n_i
     integer :: n_o
     integer :: n_e
   contains
     procedure(build), deferred :: build_i
     procedure(build), deferred :: build_o
   end type c_bound_t

  ! Interfaces

  abstract interface

     subroutine build (this, st, B, scl)
       use core_kinds
       use gyre_state
       import c_bound_t
       class(c_bound_t), intent(in) :: this
       class(c_state_t), intent(in) :: st
       complex(WP), intent(out)          :: B(:,:)
       complex(WP), intent(out)          :: scl(:)
     end subroutine build

  end interface

  ! Access specifiers

  private

  public :: c_bound_t

end module gyre_c_bound


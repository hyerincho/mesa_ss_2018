!fpx3_header(0.13_3a)
!
!dependencies
!   dir: ~/Documents/MESA_SS/gyre/src/build 
!   sources: -
!   includes: ../extern/core/core.inc ../diff/gyre_diff.inc
!   uses: core_kinds gyre_ext gyre_state
!   provides: gyre_c_diff
!end dependencies
!
!end fpx3_header
! Incfile  : gyre_c_diff
! Purpose  : difference equations (complex)
!
! Copyright 2015 Rich Townsend
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

! Incfile  : gyre_diff
! Purpose  : difference equations (template)
!
! Copyright 2015-2017 Rich Townsend
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

module gyre_c_diff

  ! Uses

  use core_kinds

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type, abstract :: c_diff_t
     integer :: n_e
   contains
     procedure(build), deferred :: build
  end type c_diff_t

  ! Interfaces

  abstract interface

     subroutine build (this, st, E_l, E_r, scl)
       use core_kinds
       use gyre_ext
       use gyre_state
       import c_diff_t
       class(c_diff_t), intent(in)  :: this
       class(c_state_t), intent(in) :: st
       complex(WP), intent(out)          :: E_l(:,:)
       complex(WP), intent(out)          :: E_r(:,:)
       type(c_ext_t), intent(out)   :: scl
     end subroutine build

  end interface

  ! Access specifiers

  private

  public :: c_diff_t

end module gyre_c_diff


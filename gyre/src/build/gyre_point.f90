!fpx3_header(0.13_3a)
!
!dependencies
!   dir: ~/Documents/MESA_SS/gyre/src/build 
!   sources: -
!   includes: ../extern/core/core.inc ../extern/core/core_parallel.inc
!   uses: core_kinds core_parallel
!   provides: gyre_point
!end dependencies
!
!end fpx3_header
! Module   : gyre_point
! Purpose  : segmented grid point
!
! Copyright 2013-2016 Rich Townsend
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

! Incfile  : core_parallel
! Purpose  : parallel support fpx3 macros

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

!****

!****

!****

!****

!****

module gyre_point

  ! Uses

  use core_kinds
  use core_parallel

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type :: point_t
     integer  :: s
     real(WP) :: x
   contains
     private
     procedure       :: op_eq_
     generic, public :: operator(==) => op_eq_
  end type point_t

  ! Interfaces

  ! Access specifiers

  private

  public :: point_t

  ! Procedures

contains

  elemental function op_eq_ (this, that) result (eq)

    class(point_t), intent(in) :: this
    class(point_t), intent(in) :: that
    logical                    :: eq

    ! Evaluate the equality operator

    eq = this%s == that%s .AND. this%x == that%x

    ! Finish

    return

  end function op_eq_

  !****

end module gyre_point

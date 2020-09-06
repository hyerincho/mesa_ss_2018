!fpx3_header(0.13_3a)
!
!dependencies
!   dir: ~/Documents/MESA_SS/gyre/src/build 
!   sources: -
!   includes: ../extern/core/core.inc
!   uses: core_kinds
!   provides: gyre_r_state
!end dependencies
!
!end fpx3_header
! Module   : gyre_r_state
! Purpose  : mutable state data for solvers (real)
!
! Copyright 2017 Rich Townsend
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

module gyre_r_state

  ! Uses

  use core_kinds

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type :: r_state_t
     real(WP) :: omega
  end type r_state_t

  ! Access specifiers

  private

  public :: r_state_t

end module gyre_r_state

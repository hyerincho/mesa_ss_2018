!fpx3_header(0.13_3a)
!
!dependencies
!   dir: ~/Documents/MESA_SS/gyre/src/build 
!   sources: -
!   includes: ../extern/core/core.inc
!   uses: gyre_r_state gyre_c_state
!   provides: gyre_state
!end dependencies
!
!end fpx3_header
! Module   : gyre_state
! Purpose  : mutable state data for solvers
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

module gyre_state

  ! Uses

  use gyre_r_state
  use gyre_c_state

  ! No implicit typing

  implicit none

  ! Access specifiers

  private

  public :: r_state_t
  public :: c_state_t

end module gyre_state

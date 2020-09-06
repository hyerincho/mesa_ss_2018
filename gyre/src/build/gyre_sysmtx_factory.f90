!fpx3_header(0.13_3a)
!
!dependencies
!   dir: ~/Documents/MESA_SS/gyre/src/build 
!   sources: -
!   includes: ../extern/core/core.inc
!   uses: gyre_num_par ISO_FORTRAN_ENV gyre_band_sysmtx core_kinds gyre_sysmtx gyre_block_sysmtx
!   provides: gyre_sysmtx_factory
!end dependencies
!
!end fpx3_header
! Module   : gyre_sysmtx_factory
! Purpose  : factory procedures for r_sysmtx_t and c_sysmtx_t types
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

module gyre_sysmtx_factory

  ! Uses

  use core_kinds

  use gyre_sysmtx
  use gyre_band_sysmtx
  use gyre_block_sysmtx
  use gyre_num_par

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Interfaces

  interface r_sysmtx_t
     module procedure r_sysmtx_t_
  end interface r_sysmtx_t

  interface c_sysmtx_t
     module procedure c_sysmtx_t_
  end interface c_sysmtx_t

  ! Access specifiers

  private

  public :: r_sysmtx_t
  public :: c_sysmtx_t

  ! Procedures

contains

  function r_sysmtx_t_ (n, n_e, n_i, n_o, np) result (sm)

    integer, intent(in)               :: n
    integer, intent(in)               :: n_e
    integer, intent(in)               :: n_i
    integer, intent(in)               :: n_o
    type(num_par_t), intent(in)       :: np
    class(r_sysmtx_t), allocatable :: sm

    ! Create a ${T}_sysmtx_t

    select case (np%matrix_type)
    case ('BAND')
       allocate(sm, SOURCE=r_band_sysmtx_t(n, n_e, n_i, n_o))
    case ('BLOCK')
       allocate(sm, SOURCE=r_block_sysmtx_t(n, n_e, n_i, n_o))
    case default

    write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 90 <gyre_sysmtx_factory:r_sysmtx_t_>:'
    write(UNIT=ERROR_UNIT, FMT=*) 'Invalid matrix_type'

  stop 'Program aborted'

    end select

    ! Finish

    return

  end function r_sysmtx_t_

  function c_sysmtx_t_ (n, n_e, n_i, n_o, np) result (sm)

    integer, intent(in)               :: n
    integer, intent(in)               :: n_e
    integer, intent(in)               :: n_i
    integer, intent(in)               :: n_o
    type(num_par_t), intent(in)       :: np
    class(c_sysmtx_t), allocatable :: sm

    ! Create a ${T}_sysmtx_t

    select case (np%matrix_type)
    case ('BAND')
       allocate(sm, SOURCE=c_band_sysmtx_t(n, n_e, n_i, n_o))
    case ('BLOCK')
       allocate(sm, SOURCE=c_block_sysmtx_t(n, n_e, n_i, n_o))
    case default

    write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 91 <gyre_sysmtx_factory:c_sysmtx_t_>:'
    write(UNIT=ERROR_UNIT, FMT=*) 'Invalid matrix_type'

  stop 'Program aborted'

    end select

    ! Finish

    return

  end function c_sysmtx_t_

end module gyre_sysmtx_factory

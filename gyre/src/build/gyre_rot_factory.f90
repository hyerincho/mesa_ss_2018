!fpx3_header(0.13_3a)
!
!dependencies
!   dir: ~/Documents/MESA_SS/gyre/src/build 
!   sources: -
!   includes: ../extern/core/core.inc
!   uses: core_kinds gyre_osc_par gyre_mode_par gyre_dopp_rot gyre_tar_rot ISO_FORTRAN_ENV gyre_rot
!   provides: gyre_rot_factory
!end dependencies
!
!end fpx3_header
! Incfile  : gyre_rot_factory
! Purpose  : factory procedures for r_rot_t and c_rot_t types
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

module gyre_rot_factory

  ! Uses

  use core_kinds

  use gyre_mode_par
  use gyre_osc_par
  use gyre_rot

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Interfaces

  interface r_rot_t
     module procedure r_rot_t_
  end interface r_rot_t

  interface c_rot_t
     module procedure c_rot_t_
  end interface c_rot_t

  ! Access specifiers

  private

  public :: r_rot_t
  public :: c_rot_t

  ! Procedures

contains

  function r_rot_t_ (md_p, os_p) result (rt)

    use gyre_dopp_rot

    use gyre_tar_rot

    type(mode_par_t), intent(in)   :: md_p
    type(osc_par_t), intent(in)    :: os_p
    class(r_rot_t), allocatable :: rt

    ! Create a ${T}_rot_t

    select case (os_p%rotation_method)
    case ('DOPPLER')
       allocate(rt, SOURCE=r_dopp_rot_t(md_p))
    case ('TAR')

       allocate(rt, SOURCE=r_tar_rot_t(md_p))

    case default

    write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 95 <gyre_rot_factory:r_rot_t_>:'
    write(UNIT=ERROR_UNIT, FMT=*) 'Invalid rotation_method'

  stop 'Program aborted'

    end select

    ! Finish

    return

  end function r_rot_t_

  function c_rot_t_ (md_p, os_p) result (rt)

    use gyre_dopp_rot

    use gyre_tar_rot

    type(mode_par_t), intent(in)   :: md_p
    type(osc_par_t), intent(in)    :: os_p
    class(c_rot_t), allocatable :: rt

    ! Create a ${T}_rot_t

    select case (os_p%rotation_method)
    case ('DOPPLER')
       allocate(rt, SOURCE=c_dopp_rot_t(md_p))
    case ('TAR')

       allocate(rt, SOURCE=c_tar_rot_t(md_p))

    case default

    write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 96 <gyre_rot_factory:c_rot_t_>:'
    write(UNIT=ERROR_UNIT, FMT=*) 'Invalid rotation_method'

  stop 'Program aborted'

    end select

    ! Finish

    return

  end function c_rot_t_

end module gyre_rot_factory

!fpx3_header(0.13_3a)
!
!dependencies
!   dir: ~/Documents/MESA_SS/gyre/src/build 
!   sources: -
!   includes: ../extern/core/core.inc
!   uses: gyre_colloc_diff gyre_mirk_diff ISO_FORTRAN_ENV gyre_diff gyre_num_par core_kinds gyre_eqns gyre_magnus_diff gyre_trapz_diff gyre_point
!   provides: gyre_diff_factory
!end dependencies
!
!end fpx3_header
! Module   : gyre_diff_factory
! Purpose  : factory procedures for r_diff_t and c_diff_t types
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

module gyre_diff_factory

  ! Uses

  use core_kinds

  use gyre_colloc_diff
  use gyre_diff
  use gyre_eqns
  use gyre_magnus_diff
  use gyre_mirk_diff
  use gyre_num_par
  use gyre_point
  use gyre_trapz_diff

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Interfaces

  interface r_diff_t
     module procedure r_diff_t_
  end interface r_diff_t

  interface c_diff_t
     module procedure c_diff_t_
  end interface c_diff_t

  ! Access specifiers

  private

  public :: r_diff_t
  public :: c_diff_t

  ! Procedures

contains

  function r_diff_t_ (eq, pt_a, pt_b, nm_p) result (df)

    class(r_eqns_t), intent(in)  :: eq
    type(point_t), intent(in)       :: pt_a
    type(point_t), intent(in)       :: pt_b
    type(num_par_t), intent(in)     :: nm_p
    class(r_diff_t), allocatable :: df

    ! Construct the ${T}_diff_t

    select case (nm_p%diff_scheme)
    case ('MAGNUS_GL2')
       allocate(df, SOURCE=r_magnus_diff_t(eq, pt_a, pt_b, 'GL2'))
    case ('MAGNUS_GL4')
       allocate(df, SOURCE=r_magnus_diff_t(eq, pt_a, pt_b, 'GL4'))
    case ('MAGNUS_GL6')
       allocate(df, SOURCE=r_magnus_diff_t(eq, pt_a, pt_b, 'GL6'))
    case ('COLLOC_GL2')
       allocate(df, SOURCE=r_colloc_diff_t(eq, pt_a, pt_b, 'GL2'))
    case ('COLLOC_GL4')
       allocate(df, SOURCE=r_colloc_diff_t(eq, pt_a, pt_b, 'GL4'))
    case ('COLLOC_GL6')
       allocate(df, SOURCE=r_colloc_diff_t(eq, pt_a, pt_b, 'GL6'))
    case ('TRAPZ')
       allocate(df, SOURCE=r_trapz_diff_t(eq, pt_a, pt_b, SPREAD(0.5_WP, DIM=1, NCOPIES=eq%n_e)))
    case ('MIRK')
       allocate(df, SOURCE=r_mirk_diff_t(eq, pt_a, pt_b))
    case default

    write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 105 <gyre_diff_factory:r_diff_t_>:'
    write(UNIT=ERROR_UNIT, FMT=*) 'Invalid diff_scheme'

  stop 'Program aborted'

    end select

    ! Finish

    return

  end function r_diff_t_

  function c_diff_t_ (eq, pt_a, pt_b, nm_p) result (df)

    class(c_eqns_t), intent(in)  :: eq
    type(point_t), intent(in)       :: pt_a
    type(point_t), intent(in)       :: pt_b
    type(num_par_t), intent(in)     :: nm_p
    class(c_diff_t), allocatable :: df

    ! Construct the ${T}_diff_t

    select case (nm_p%diff_scheme)
    case ('MAGNUS_GL2')
       allocate(df, SOURCE=c_magnus_diff_t(eq, pt_a, pt_b, 'GL2'))
    case ('MAGNUS_GL4')
       allocate(df, SOURCE=c_magnus_diff_t(eq, pt_a, pt_b, 'GL4'))
    case ('MAGNUS_GL6')
       allocate(df, SOURCE=c_magnus_diff_t(eq, pt_a, pt_b, 'GL6'))
    case ('COLLOC_GL2')
       allocate(df, SOURCE=c_colloc_diff_t(eq, pt_a, pt_b, 'GL2'))
    case ('COLLOC_GL4')
       allocate(df, SOURCE=c_colloc_diff_t(eq, pt_a, pt_b, 'GL4'))
    case ('COLLOC_GL6')
       allocate(df, SOURCE=c_colloc_diff_t(eq, pt_a, pt_b, 'GL6'))
    case ('TRAPZ')
       allocate(df, SOURCE=c_trapz_diff_t(eq, pt_a, pt_b, SPREAD(0.5_WP, DIM=1, NCOPIES=eq%n_e)))
    case ('MIRK')
       allocate(df, SOURCE=c_mirk_diff_t(eq, pt_a, pt_b))
    case default

    write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 106 <gyre_diff_factory:c_diff_t_>:'
    write(UNIT=ERROR_UNIT, FMT=*) 'Invalid diff_scheme'

  stop 'Program aborted'

    end select

    ! Finish

    return

  end function c_diff_t_

end module gyre_diff_factory

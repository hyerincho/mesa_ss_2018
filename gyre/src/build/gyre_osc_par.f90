!fpx3_header(0.13_3a)
!
!dependencies
!   dir: ~/Documents/MESA_SS/gyre/src/build 
!   sources: -
!   includes: ../extern/core/core_parallel.inc ../extern/core/core.inc
!   uses: core_parallel core_kinds core_constants
!   provides: gyre_osc_par
!end dependencies
!
!end fpx3_header
! Module   : gyre_osc_par
! Purpose  : oscillation parameters
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

module gyre_osc_par

  ! Uses

  use core_kinds
  use core_constants, only : FILENAME_LEN
  use core_parallel

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type :: osc_par_t
     real(WP)                :: x_ref
     character(64)           :: rotation_method
     character(64)           :: variables_set
     character(64)           :: inner_bound
     character(64)           :: outer_bound
     character(64)           :: inertia_norm
     character(64)           :: time_factor
     character(64)           :: conv_scheme
     character(64)           :: deps_scheme
     character(FILENAME_LEN) :: deps_file
     character(256)          :: deps_file_format
     character(2048)         :: tag_list
     logical                 :: nonadiabatic
     logical                 :: quasiad_eigfuncs
     logical                 :: cowling_approx
     logical                 :: nar_approx
     logical                 :: narf_approx
     logical                 :: eddington_approx
     logical                 :: complex_lambda
     logical                 :: reduce_order
  end type osc_par_t

  ! Interfaces

 ! Access specifiers

  private

  public :: osc_par_t
  public :: read_osc_par

  ! Procedures

contains

!****

  subroutine read_osc_par (unit, os_p)

    integer, intent(in)                       :: unit
    type(osc_par_t), allocatable, intent(out) :: os_p(:)

    integer                               :: n_os_p
    integer                               :: i
    real(WP)                              :: x_ref
    character(LEN(os_p%rotation_method))  :: rotation_method
    character(LEN(os_p%variables_set))    :: variables_set
    character(LEN(os_p%inner_bound))      :: inner_bound
    character(LEN(os_p%outer_bound))      :: outer_bound
    character(LEN(os_p%inertia_norm))     :: inertia_norm
    character(LEN(os_p%time_factor))      :: time_factor
    character(LEN(os_p%conv_scheme))      :: conv_scheme
    character(LEN(os_p%deps_scheme))      :: deps_scheme
    character(LEN(os_p%deps_file))        :: deps_file
    character(LEN(os_p%deps_file_format)) :: deps_file_format
    character(LEN(os_p%tag_list))         :: tag_list
    logical                               :: quasiad_eigfuncs
    logical                               :: nonadiabatic
    logical                               :: cowling_approx
    logical                               :: nar_approx
    logical                               :: narf_approx
    logical                               :: eddington_approx
    logical                               :: complex_lambda
    logical                               :: reduce_order

    namelist /osc/ x_ref, rotation_method, inner_bound, outer_bound, &
         variables_set, inertia_norm, time_factor, &
         conv_scheme, deps_scheme, deps_file, deps_file_format, &
         tag_list, quasiad_eigfuncs, nonadiabatic, &
         cowling_approx, nar_approx, narf_approx, eddington_approx, &
         complex_lambda, reduce_order

    ! Count the number of osc namelists

    rewind(unit)

    n_os_p = 0

    count_loop : do
       read(unit, NML=osc, END=100)
       n_os_p = n_os_p + 1
    end do count_loop

100 continue

    ! Read oscillation parameters

    rewind(unit)

    allocate(os_p(n_os_p))

    read_loop : do i = 1,n_os_p

       x_ref = HUGE(0._WP)

       rotation_method = 'DOPPLER'
       variables_set = 'GYRE'
       inner_bound = 'REGULAR'
       outer_bound = 'VACUUM'
       inertia_norm = 'BOTH'
       time_factor = 'OSC'
       conv_scheme = 'FROZEN_PESNELL_1'
       deps_scheme = 'MODEL'
       deps_file = ''
       deps_file_format = ''
       tag_list = ''

       quasiad_eigfuncs = .FALSE.
       nonadiabatic = .FALSE.
       cowling_approx = .FALSE.
       nar_approx = .FALSE.
       narf_approx = .FALSE.
       eddington_approx = .FALSE.
       complex_lambda = .FALSE.
       reduce_order = .TRUE.

       read(unit, NML=osc)

       ! Initialize the osc_par

       os_p(i) = osc_par_t(x_ref=x_ref, &
                           rotation_method=rotation_method, &
                           variables_set=variables_set, &
                           inner_bound=inner_bound, &
                           outer_bound=outer_bound, &
                           inertia_norm=inertia_norm, &
                           time_factor=time_factor, &
                           conv_scheme=conv_scheme, &
                           deps_scheme=deps_scheme, &
                           deps_file=deps_file, &
                           deps_file_format=deps_file_format, &
                           tag_list=tag_list, &
                           quasiad_eigfuncs=quasiad_eigfuncs, &
                           nonadiabatic=nonadiabatic, &
                           cowling_approx=cowling_approx, &
                           nar_approx=nar_approx, &
                           narf_approx=narf_approx, &
                           eddington_approx=eddington_approx, &
                           complex_lambda=complex_lambda, &
                           reduce_order=reduce_order)

    end do read_loop

    ! Finish

    return

  end subroutine read_osc_par

  !****

end module gyre_osc_par

!fpx3_header(0.13_3a)
!
!dependencies
!   dir: ~/Documents/MESA_SS/gyre/src/build 
!   sources: -
!   includes: ../extern/core/core.inc
!   uses: gyre_b3_file gyre_mesa_file gyre_evol_model gyre_poly_file gyre_losc_file gyre_model gyre_famdl_file gyre_model_par gyre_hom_model core_kinds gyre_poly_model gyre_fgong_file gyre_gsm_file gyre_osc_file gyre_amdl_file ISO_FORTRAN_ENV
!   provides: gyre_model_factory
!end dependencies
!
!end fpx3_header
! Module   : gyre_model_factory
! Purpose  : factory procedures for model_t
!
! Copyright 2016 Rich Townsend
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

module gyre_model_factory

  ! Uses

  use core_kinds

  use gyre_evol_model
  use gyre_hom_model
  use gyre_model
  use gyre_model_par
  use gyre_poly_model

  use gyre_amdl_file
  use gyre_famdl_file
  use gyre_fgong_file
  use gyre_losc_file
  use gyre_mesa_file
  use gyre_osc_file

  use gyre_b3_file
  use gyre_gsm_file
  use gyre_poly_file

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Interfaces

  interface model_t
     module procedure model_t_
  end interface model_t

  ! Access specifiers

  private

  public :: model_t

  ! Procedures

contains

  function model_t_ (ml_p) result (ml)

    type(model_par_t), intent(in) :: ml_p
    class(model_t), pointer       :: ml

    ! Construct the model_t

    select case (ml_p%model_type)
    case ('EVOL')

       select case (ml_p%file_format)
       case ('AMDL')
          call read_amdl_model(ml_p, ml)
       case ('B3')

          call read_b3_model(ml_p, ml)

       case ('FAMDL')
          call read_famdl_model(ml_p, ml)
       case ('FGONG')
          call read_fgong_model(ml_p, ml)
       case ('GSM')

          call read_gsm_model(ml_p, ml)

       case ('LOSC')
          call read_losc_model(ml_p, ml)
       case ('MESA')
          call read_mesa_model(ml_p, ml)
       case ('OSC')
          call read_osc_model(ml_p, ml)
       case default

    write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 103 <gyre_model_factory:model_t_>:'
    write(UNIT=ERROR_UNIT, FMT=*) 'Invalid file_format'

  stop 'Program aborted'

       end select

    case ('POLY')

       call read_poly_model(ml_p, ml)

    case ('HOM')

       allocate(ml, SOURCE=hom_model_t(ml_p))

    case default

    write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 120 <gyre_model_factory:model_t_>:'
    write(UNIT=ERROR_UNIT, FMT=*) 'Invalid model_type'

  stop 'Program aborted'

    end select

    ! Finish

    return

  end function model_t_

  !****

end module gyre_model_factory

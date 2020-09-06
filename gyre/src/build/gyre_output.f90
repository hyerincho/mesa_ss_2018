!fpx3_header(0.13_3a)
!
!dependencies
!   dir: ~/Documents/MESA_SS/gyre/src/build 
!   sources: -
!   includes: ../output/gyre_output.inc ../extern/core/core.inc
!   uses: gyre_constants ISO_FORTRAN_ENV gyre_freq gyre_evol_model gyre_model gyre_mode gyre_out_par core_string gyre_point gyre_writer gyre_txt_writer gyre_hdf_writer gyre_util gyre_poly_model core_kinds
!   provides: gyre_output
!end dependencies
!
!end fpx3_header
! Module   : gyre_output
! Purpose  : output routines
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

! Incfile  : gyre_output
! Purpose  : macros for output module
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

module gyre_output

  ! Uses

  use core_kinds
  use core_string

  use gyre_constants
  use gyre_evol_model
  use gyre_freq
  use gyre_hdf_writer
  use gyre_mode
  use gyre_model
  use gyre_out_par
  use gyre_point
  use gyre_poly_model
  use gyre_txt_writer
  use gyre_util
  use gyre_writer

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Access specifiers

  private

  public :: write_summary
  public :: write_mode

contains

  subroutine write_summary (md, ot_p)

    type(mode_t), intent(in)    :: md(:)
    type(out_par_t), intent(in) :: ot_p

    class(writer_t), allocatable                        :: wr
    character(LEN(ot_p%summary_item_list)), allocatable :: items(:)
    logical                                             :: invalid_items
    integer                                             :: n_md
    real(WP)                                            :: data_r(SIZE(md))
    complex(WP)                                         :: data_c(SIZE(md))
    integer                                             :: i
    integer                                             :: i_md

    ! Write the summary file

    if (SIZE(md) == 0 .OR. ot_p%summary_file == '') return

    ! Open the file

    select case (ot_p%summary_file_format)
    case ('HDF')
       allocate(wr, SOURCE=hdf_writer_t(ot_p%summary_file, ot_p%label))
    case ('TXT')
       allocate(wr, SOURCE=txt_writer_t(ot_p%summary_file, ot_p%label))
    case default

    write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 82 <gyre_output:write_summary>:'
    write(UNIT=ERROR_UNIT, FMT=*) 'Invalid summary_file_format'

  stop 'Program aborted'

    end select

    ! Split the item list

    items = split_list(ot_p%summary_item_list, ',')

    ! Write the items

    invalid_items = .FALSE.

    n_md = SIZE(md)

    item_loop : do i = 1, SIZE(items)

       select case (items(i))

       case ('j')

          call wr%write('j', md%j)

       case ('l')

          call wr%write('l', md%l)

       case ('l_i')

          call wr%write('l_i', md%l_i)

       case ('m')

          call wr%write('m', md%m)

       case ('n_p')

          call wr%write('n_p', md%n_p)

       case ('n_g')

          call wr%write('n_g', md%n_g)

       case ('n_pg')

          call wr%write('n_pg', md%n_pg)

       case ('omega')

          call wr%write('omega', md%omega)

       case ('omega_int')

          do i_md = 1, n_md
             data_c(i_md) = md(i_md)%omega_int()
          end do

          call wr%write('omega_int', data_c)

       case ('freq')

          do i_md = 1, n_md
             data_c(i_md) = md(i_md)%freq(ot_p%freq_units, ot_p%freq_frame)
          end do

          call wr%write('freq', data_c)

       case ('freq_units')

          call wr%write('freq_units', ot_p%freq_units)

       case ('freq_frame')

          call wr%write('freq_frame', ot_p%freq_frame)

       case ('Delta_p')

          do i_md = 1, n_md
             data_r(i_md) = md(i_md)%cx%ml%Delta_p(md(i_md)%gr%x_i(), md(i_md)%gr%x_o())
          end do

          call wr%write('Delta_p', data_r)

       case ('Delta_g')

          do i_md = 1, n_md
             data_r(i_md) = md(i_md)%cx%ml%Delta_g(md(i_md)%gr%x_i(), md(i_md)%gr%x_o(), &
                                                   md(i_md)%l*(md(i_md)%l+1._WP))
          end do

          call wr%write('Delta_g', data_r)

case ('eta')

   do i_md = 1, n_md
      data_r(i_md) = md(i_md)%eta()
   end do

   call wr%write('eta', data_r)

case ('f_T')

   do i_md = 1, n_md
      data_r(i_md) = md(i_md)%f_T()
   end do

   call wr%write('f_T', data_r)

case ('f_g')

   do i_md = 1, n_md
      data_r(i_md) = md(i_md)%f_g()
   end do

   call wr%write('f_g', data_r)

case ('psi_T')

   do i_md = 1, n_md
      data_r(i_md) = md(i_md)%psi_T()
   end do

   call wr%write('psi_T', data_r)

case ('psi_g')

   do i_md = 1, n_md
      data_r(i_md) = md(i_md)%psi_g()
   end do

   call wr%write('psi_g', data_r)

case ('E')

   do i_md = 1, n_md
      data_r(i_md) = md(i_md)%E()
   end do

   call wr%write('E', data_r)

case ('E_p')

   do i_md = 1, n_md
      data_r(i_md) = md(i_md)%E_p()
   end do

   call wr%write('E_p', data_r)

case ('E_g')

   do i_md = 1, n_md
      data_r(i_md) = md(i_md)%E_g()
   end do

   call wr%write('E_g', data_r)

case ('E_norm')

   do i_md = 1, n_md
      data_r(i_md) = md(i_md)%E_norm()
   end do

   call wr%write('E_norm', data_r)

case ('E_ratio')

   do i_md = 1, n_md
      data_r(i_md) = md(i_md)%E_ratio()
   end do

   call wr%write('E_ratio', data_r)

case ('H')

   do i_md = 1, n_md
      data_r(i_md) = md(i_md)%H()
   end do

   call wr%write('H', data_r)

case ('W')

   do i_md = 1, n_md
      data_r(i_md) = md(i_md)%W()
   end do

   call wr%write('W', data_r)

case ('W_eps')

   do i_md = 1, n_md
      data_r(i_md) = md(i_md)%W_eps()
   end do

   call wr%write('W_eps', data_r)

case ('tau_ss')

   do i_md = 1, n_md
      data_r(i_md) = md(i_md)%tau_ss()
   end do

   call wr%write('tau_ss', data_r)

case ('tau_tr')

   do i_md = 1, n_md
      data_r(i_md) = md(i_md)%tau_tr()
   end do

   call wr%write('tau_tr', data_r)

case ('beta')

   do i_md = 1, n_md
      data_r(i_md) = md(i_md)%beta()
   end do

   call wr%write('beta', data_r)

case ('x_ref')

   do i_md = 1, n_md
      data_r(i_md) = md(i_md)%gr%pt(md(i_md)%k_ref)%x
   end do

   call wr%write('x_ref', data_r)

case ('xi_r_ref')

   do i_md = 1, n_md
      data_c(i_md) = md(i_md)%xi_r(md(i_md)%k_ref)
   end do

   call wr%write('xi_r_ref', data_c)

case ('xi_h_ref')

   do i_md = 1, n_md
      data_c(i_md) = md(i_md)%xi_h(md(i_md)%k_ref)
   end do

   call wr%write('xi_h_ref', data_c)

case ('eul_phi_ref')

   do i_md = 1, n_md
      data_c(i_md) = md(i_md)%eul_phi(md(i_md)%k_ref)
   end do

   call wr%write('eul_phi_ref', data_c)

case ('deul_phi_ref')

   do i_md = 1, n_md
      data_c(i_md) = md(i_md)%deul_phi(md(i_md)%k_ref)
   end do

   call wr%write('deul_phi_ref', data_c)

case ('lag_S_ref')

   do i_md = 1, n_md
      data_c(i_md) = md(i_md)%lag_S(md(i_md)%k_ref)
   end do

   call wr%write('lag_S_ref', data_c)

case ('lag_L_ref')

   do i_md = 1, n_md
      data_c(i_md) = md(i_md)%lag_L(md(i_md)%k_ref)
   end do

   call wr%write('lag_L_ref', data_c)

       case default

          if (n_md >= 1) then
             select type (ml => md(1)%cx%ml)
             type is (evol_model_t)
                call write_summary_evol_(items(i), ml, wr)
             class default
                write(ERROR_UNIT, *) 'item:', TRIM(items(i))
                invalid_items = .TRUE.
             end select
          endif

       end select

    end do item_loop

    ! Write restart metadata to HDF5 files

    select case (ot_p%summary_file_format)
    case ('HDF')
       call wr%write('i', [(md(i_md)%md_p%i, i_md=1,n_md)])
    end select

    ! Close the file

    call wr%final()

    ! Check whether any invalid items were found

    if (invalid_items) then

    write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 226 <gyre_output:write_summary>:'
    write(UNIT=ERROR_UNIT, FMT=*) 'Invalid item(s) in summary_item_list'

  stop 'Program aborted'

    end if

    ! Finish

    return

  contains

    subroutine write_summary_evol_ (item, ml, wr)

      character(*), intent(in)       :: item
      type(evol_model_t), intent(in) :: ml
      class(writer_t), intent(inout) :: wr

      ! Write the item

      select case (items(i))
      case ('M_star')
         call wr%write('M_star', ml%M_star)
      case ('R_star')
         call wr%write('R_star', ml%R_star)
      case ('L_star')
         call wr%write('L_star', ml%L_star)
      case default
         write(ERROR_UNIT, *) 'item:', TRIM(items(i))
         invalid_items = .TRUE.
      end select

      ! Finish

      return

    end subroutine write_summary_evol_

  end subroutine write_summary

  !****

  subroutine write_mode (md, ot_p)

    type(mode_t), intent(in)    :: md
    type(out_par_t), intent(in) :: ot_p

    character(:), allocatable                        :: mode_file
    class(writer_t), allocatable                     :: wr
    character(LEN(ot_p%mode_item_list)), allocatable :: items(:)
    logical                                          :: invalid_items
    integer                                          :: i
    integer                                          :: k
    real(WP)                                         :: data_r(md%n_k)
    complex(WP)                                      :: data_c(md%n_k)

    ! Write the mode file

    if (ot_p%mode_template == '') return

    if (filter_mode_(md, ot_p)) return

    ! Set up the filename

    mode_file = ot_p%mode_template

    mode_file = subst_(mode_file, '%J', md%j, '(I5.5)')
    mode_file = subst_(mode_file, '%L', md%l, '(I3.3)')
    mode_file = subst_(mode_file, '%M', md%m, '(SP,I3.2)')
    mode_file = subst_(mode_file, '%N', md%n_pg, '(SP,I6.5)')

    mode_file = subst_(mode_file, '%j', md%j, '(I0)')
    mode_file = subst_(mode_file, '%l', md%l, '(I0)')
    mode_file = subst_(mode_file, '%m', md%m, '(SP,I0)')
    mode_file = subst_(mode_file, '%n', md%n_pg, '(SP,I0)')

    ! Open the file

    select case (ot_p%mode_file_format)
    case ('HDF')
       allocate(wr, SOURCE=hdf_writer_t(mode_file, ot_p%label))
    case ('TXT')
       allocate(wr, SOURCE=txt_writer_t(mode_file, ot_p%label))
    case default

    write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 307 <gyre_output:write_mode>:'
    write(UNIT=ERROR_UNIT, FMT=*) 'Invalid mode_file_format'

  stop 'Program aborted'

    end select

    ! Split the item list

    items = split_list(ot_p%mode_item_list, ',')

    ! Write the items

    invalid_items = .FALSE.

    associate (pt => md%gr%pt)

      item_loop : do i = 1, SIZE(items)

         select case (items(i))

         case ('n')

            call wr%write('n', md%n_k)

         case ('j')

            call wr%write('j', md%j)

         case ('l')

            call wr%write('l', md%l)

         case ('l_i')

            call wr%write('l_i', md%l_i)

         case ('m')

            call wr%write('m', md%m)

         case ('lambda')

            call wr%write('lambda', [(md%lambda(k), k=1,md%n_k)])

         case ('n_p')

            call wr%write('n_p', md%n_p)

         case ('n_g')

            call wr%write('n_g', md%n_g)

         case ('n_pg')

            call wr%write('n_pg', md%n_pg)

         case ('omega')

            call wr%write('omega', md%omega)

         case ('omega_int')

            call wr%write('omega_int', md%omega_int())

         case ('freq')

            call wr%write('freq', md%freq(ot_p%freq_units, ot_p%freq_frame))

         case ('freq_units')

            call wr%write('freq_units', ot_p%freq_units)

         case ('freq_frame')

            call wr%write('freq_frame', ot_p%freq_frame)

         case ('Delta_p')

            call wr%write('Delta_p', md%cx%ml%Delta_p(md%gr%x_i(), md%gr%x_o()))

         case ('Delta_g')

            call wr%write('Delta_g', md%cx%ml%Delta_g(md%gr%x_i(), md%gr%x_o(), md%l*(md%l+1._WP)))

         case ('eta')

            call wr%write('eta', md%eta())

         case ('f_T')

            call wr%write('f_T', md%f_T())

         case ('f_g')

            call wr%write('f_g', md%f_g())

         case ('psi_T')

            call wr%write('psi_T', md%psi_T())

         case ('psi_g')

            call wr%write('psi_g', md%psi_g())

         case ('E')

            call wr%write('E', md%E())

         case ('E_p')

            call wr%write('E_p', md%E_p())

         case ('E_g')

            call wr%write('E_g', md%E_g())

         case ('E_norm')

            call wr%write('E_norm', md%E_norm())

         case ('E_ratio')

            call wr%write('E_ratio', md%E_ratio())

         case ('H')

            call wr%write('H', md%H())

         case ('W')

            call wr%write('W', md%W())

         case ('W_eps')

            call wr%write('W_eps', md%W_eps())

         case ('tau_ss')

            call wr%write('tau_ss', md%tau_ss())

         case ('tau_tr')

            call wr%write('tau_tr', md%tau_tr())

         case ('beta')

            call wr%write('beta', md%beta())

         case ('x')

            call wr%write('x', md%gr%pt%x)

         case ('x_ref')

            call wr%write('x_ref', pt(md%k_ref)%x)

case ('V_2')

   do k = 1, SIZE(pt)
      data_r(k) = md%cx%ml%coeff(I_V_2, pt(k))
   end do

   call wr%write('V_2', data_r)

case ('As')

   do k = 1, SIZE(pt)
      data_r(k) = md%cx%ml%coeff(I_AS, pt(k))
   end do

   call wr%write('As', data_r)

case ('U')

   do k = 1, SIZE(pt)
      data_r(k) = md%cx%ml%coeff(I_U, pt(k))
   end do

   call wr%write('U', data_r)

case ('c_1')

   do k = 1, SIZE(pt)
      data_r(k) = md%cx%ml%coeff(I_C_1, pt(k))
   end do

   call wr%write('c_1', data_r)

case ('Gamma_1')

   do k = 1, SIZE(pt)
      data_r(k) = md%cx%ml%coeff(I_GAMMA_1, pt(k))
   end do

   call wr%write('Gamma_1', data_r)

case ('nabla')

   do k = 1, SIZE(pt)
      data_r(k) = md%cx%ml%coeff(I_NABLA, pt(k))
   end do

   call wr%write('nabla', data_r)

case ('nabla_ad')

   do k = 1, SIZE(pt)
      data_r(k) = md%cx%ml%coeff(I_NABLA_AD, pt(k))
   end do

   call wr%write('nabla_ad', data_r)

case ('dnabla_ad')

   do k = 1, SIZE(pt)
      data_r(k) = md%cx%ml%dcoeff(I_NABLA_AD, pt(k))
   end do

   call wr%write('dnabla_ad', data_r)

case ('delta')

   do k = 1, SIZE(pt)
      data_r(k) = md%cx%ml%coeff(I_DELTA, pt(k))
   end do

   call wr%write('delta', data_r)

case ('c_lum')

   do k = 1, SIZE(pt)
      data_r(k) = md%cx%ml%coeff(I_C_LUM, pt(k))
   end do

   call wr%write('c_lum', data_r)

case ('c_rad')

   do k = 1, SIZE(pt)
      data_r(k) = md%cx%ml%coeff(I_C_RAD, pt(k))
   end do

   call wr%write('c_rad', data_r)

case ('c_thn')

   do k = 1, SIZE(pt)
      data_r(k) = md%cx%ml%coeff(I_C_THN, pt(k))
   end do

   call wr%write('c_thn', data_r)

case ('c_thk')

   do k = 1, SIZE(pt)
      data_r(k) = md%cx%ml%coeff(I_C_THK, pt(k))
   end do

   call wr%write('c_thk', data_r)

case ('c_eps')

   do k = 1, SIZE(pt)
      data_r(k) = md%cx%ml%coeff(I_C_EPS, pt(k))
   end do

   call wr%write('c_eps', data_r)

case ('eps_rho')

   do k = 1, SIZE(pt)
      data_r(k) = md%cx%ml%coeff(I_EPS_RHO, pt(k))
   end do

   call wr%write('eps_rho', data_r)

case ('eps_T')

   do k = 1, SIZE(pt)
      data_r(k) = md%cx%ml%coeff(I_EPS_T, pt(k))
   end do

   call wr%write('eps_T', data_r)

case ('kap_rho')

   do k = 1, SIZE(pt)
      data_r(k) = md%cx%ml%coeff(I_KAP_RHO, pt(k))
   end do

   call wr%write('kap_rho', data_r)

case ('kap_T')

   do k = 1, SIZE(pt)
      data_r(k) = md%cx%ml%coeff(I_KAP_T, pt(k))
   end do

   call wr%write('kap_T', data_r)

case ('Omega_rot')

   do k = 1, SIZE(pt)
      data_r(k) = md%cx%ml%coeff(I_OMEGA_ROT, pt(k))
   end do

   call wr%write('Omega_rot', data_r)

case ('y_1')

   do k = 1, SIZE(pt)
      data_c(k) = md%y_i(1, k)
   end do

   call wr%write('y_1', data_c)

case ('y_2')

   do k = 1, SIZE(pt)
      data_c(k) = md%y_i(2, k)
   end do

   call wr%write('y_2', data_c)

case ('y_3')

   do k = 1, SIZE(pt)
      data_c(k) = md%y_i(3, k)
   end do

   call wr%write('y_3', data_c)

case ('y_4')

   do k = 1, SIZE(pt)
      data_c(k) = md%y_i(4, k)
   end do

   call wr%write('y_4', data_c)

case ('y_5')

   do k = 1, SIZE(pt)
      data_c(k) = md%y_i(5, k)
   end do

   call wr%write('y_5', data_c)

case ('y_6')

   do k = 1, SIZE(pt)
      data_c(k) = md%y_i(6, k)
   end do

   call wr%write('y_6', data_c)

case ('xi_r')

   do k = 1, SIZE(pt)
      data_c(k) = md%xi_r(k)
   end do

   call wr%write('xi_r', data_c)

case ('xi_h')

   do k = 1, SIZE(pt)
      data_c(k) = md%xi_h(k)
   end do

   call wr%write('xi_h', data_c)

case ('eul_phi')

   do k = 1, SIZE(pt)
      data_c(k) = md%eul_phi(k)
   end do

   call wr%write('eul_phi', data_c)

case ('deul_phi')

   do k = 1, SIZE(pt)
      data_c(k) = md%deul_phi(k)
   end do

   call wr%write('deul_phi', data_c)

case ('eul_P')

   do k = 1, SIZE(pt)
      data_c(k) = md%eul_P(k)
   end do

   call wr%write('eul_P', data_c)

case ('eul_rho')

   do k = 1, SIZE(pt)
      data_c(k) = md%eul_rho(k)
   end do

   call wr%write('eul_rho', data_c)

case ('eul_T')

   do k = 1, SIZE(pt)
      data_c(k) = md%eul_T(k)
   end do

   call wr%write('eul_T', data_c)

case ('lag_P')

   do k = 1, SIZE(pt)
      data_c(k) = md%lag_P(k)
   end do

   call wr%write('lag_P', data_c)

case ('lag_rho')

   do k = 1, SIZE(pt)
      data_c(k) = md%lag_rho(k)
   end do

   call wr%write('lag_rho', data_c)

case ('lag_T')

   do k = 1, SIZE(pt)
      data_c(k) = md%lag_T(k)
   end do

   call wr%write('lag_T', data_c)

case ('lag_S')

   do k = 1, SIZE(pt)
      data_c(k) = md%lag_S(k)
   end do

   call wr%write('lag_S', data_c)

case ('lag_L')

   do k = 1, SIZE(pt)
      data_c(k) = md%lag_L(k)
   end do

   call wr%write('lag_L', data_c)

case ('dE_dx')

   do k = 1, SIZE(pt)
      data_r(k) = md%dE_dx(k)
   end do

   call wr%write('dE_dx', data_r)

case ('dW_dx')

   do k = 1, SIZE(pt)
      data_r(k) = md%dW_dx(k)
   end do

   call wr%write('dW_dx', data_r)

case ('dW_eps_dx')

   do k = 1, SIZE(pt)
      data_r(k) = md%dW_eps_dx(k)
   end do

   call wr%write('dW_eps_dx', data_r)

case ('dbeta_dx')

   do k = 1, SIZE(pt)
      data_r(k) = md%dbeta_dx(k)
   end do

   call wr%write('dbeta_dx', data_r)

case ('dtau_dx_ss')

   do k = 1, SIZE(pt)
      data_r(k) = md%dtau_dx_ss(k)
   end do

   call wr%write('dtau_dx_ss', data_r)

case ('dtau_dx_tr')

   do k = 1, SIZE(pt)
      data_r(k) = md%dtau_dx_tr(k)
   end do

   call wr%write('dtau_dx_tr', data_r)

case ('Yt_1')

   do k = 1, SIZE(pt)
      data_c(k) = md%Yt_1(k)
   end do

   call wr%write('Yt_1', data_c)

case ('Yt_2')

   do k = 1, SIZE(pt)
      data_c(k) = md%Yt_2(k)
   end do

   call wr%write('Yt_2', data_c)

case ('I_0')

   do k = 1, SIZE(pt)
      data_c(k) = md%I_0(k)
   end do

   call wr%write('I_0', data_c)

case ('I_1')

   do k = 1, SIZE(pt)
      data_c(k) = md%I_1(k)
   end do

   call wr%write('I_1', data_c)

case ('alpha_0')

   do k = 1, SIZE(pt)
      data_r(k) = md%alpha_0(k)
   end do

   call wr%write('alpha_0', data_r)

case ('alpha_1')

   do k = 1, SIZE(pt)
      data_r(k) = md%alpha_1(k)
   end do

   call wr%write('alpha_1', data_r)

case ('prop_type')

   do k = 1, SIZE(pt)
      data_c(k) = md%prop_type(k)
   end do

   call wr%write('prop_type', data_c)

case ('xi_r_ref')

   call wr%write('xi_r_ref', md%xi_r(md%k_ref))

case ('xi_h_ref')

   call wr%write('xi_h_ref', md%xi_h(md%k_ref))

case ('eul_phi_ref')

   call wr%write('eul_phi_ref', md%eul_phi(md%k_ref))

case ('deul_phi_ref')

   call wr%write('deul_phi_ref', md%deul_phi(md%k_ref))

case ('lag_S_ref')

   call wr%write('lag_S_ref', md%lag_S(md%k_ref))

case ('lag_L_ref')

   call wr%write('lag_L_ref', md%lag_L(md%k_ref))

         case default

            select type (ml => md%cx%ml)
            type is (evol_model_t)
               call write_mode_evol_(ml)
            class default
               write(ERROR_UNIT, *) 'item:', TRIM(items(i))
               invalid_items = .TRUE.
            end select

         end select

      end do item_loop

    end associate

    ! Close the file

    call wr%final()

    ! Check whether any invalid items were found

    if (invalid_items) then

    write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 545 <gyre_output:write_mode>:'
    write(UNIT=ERROR_UNIT, FMT=*) 'Invalid item(s) in mode_item_list'

  stop 'Program aborted'

    end if

    ! Finish

    return

  contains

    subroutine write_mode_evol_ (ml)

      type(evol_model_t), intent(in) :: ml

      integer :: k

      ! Write the item

      associate (pt => md%gr%pt)

        select case (items(i))
        case ('M_star')
           call wr%write('M_star', ml%M_star)
        case ('R_star')
           call wr%write('R_star', ml%R_star)
        case ('L_star')
           call wr%write('L_star', ml%L_star)
        case ('M_r')
           call wr%write('M_r', [(ml%M_r(pt(k)), k=1,SIZE(pt))])
        case ('P')
           call wr%write('P', [(ml%P(pt(k)), k=1,SIZE(pt))])
        case ('rho')
           call wr%write('rho', [(ml%rho(pt(k)), k=1,SIZE(pt))])
        case ('T')
           call wr%write('T', [(ml%T(pt(k)), k=1,SIZE(pt))])
        case default
           write(ERROR_UNIT, *) 'item:', TRIM(items(i))
           invalid_items = .TRUE.
        end select

      end associate

      ! Finish

      return

    end subroutine write_mode_evol_

  end subroutine write_mode

  !****

  function filter_mode_ (md, ot_p) result (filter_mode)

    type(mode_t), intent(in)    :: md
    type(out_par_t), intent(in) :: ot_p
    logical                     :: filter_mode

    character(LEN(ot_p%mode_filter_list)), allocatable :: filters(:)
    integer                                            :: i

    ! Decide whether to filter the mode

    filters = split_list(ot_p%mode_filter_list, ',')

    filter_mode = .FALSE.

    item_loop : do i = 1, SIZE(filters)

       select case (filters(i))
       case ('stable')
          filter_mode = filter_mode .OR. AIMAG(md%omega) <= 0._WP
       case ('unstable')
          filter_mode = filter_mode .OR. AIMAG(md%omega) > 0._WP
       case default

    write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 619 <gyre_output:filter_mode_>:'
    write(UNIT=ERROR_UNIT, FMT=*) 'Unrecognized filter in mode_filter_list'

  stop 'Program aborted'

       end select

    end do item_loop

    ! Finish

    return

  end function filter_mode_

  !****

  function subst_ (string, pattern, i, format) result (new_string)

    character(*), intent(in)  :: string
    character(*), intent(in)  :: pattern
    integer, intent(in)       :: i
    character(*), intent(in)  :: format
    character(:), allocatable :: new_string

    character(64) :: substring

    ! Write i into the substring buffer

    write(substring, format) i

    ! Do the replacement

    new_string = replace(string, pattern, TRIM(substring), every=.TRUE.)

    ! Finish

    return

  end function subst_

end module gyre_output

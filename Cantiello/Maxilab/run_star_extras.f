! ***********************************************************************
!
!   Copyright (C) 2010  Bill Paxton
!
!   this file is part of mesa.
!
!   mesa is free software; you can redistribute it and/or modify
!   it under the terms of the gnu general library public license as published
!   by the free software foundation; either version 2 of the license, or
!   (at your option) any later version.
!
!   mesa is distributed in the hope that it will be useful,
!   but without any warranty; without even the im_planetlied warranty of
!   merchantability or fitness for a particular pur_planetose.  see the
!   gnu library general public license for more details.
!
!   you should have received a copy of the gnu library general public license
!   along with this software; if not, write to the free software
!   foundation, inc., 59 tem_planetle place, suite 330, boston, ma 02111-1307 usa
!
! ***********************************************************************

      module run_star_extras

      use star_lib
      use star_def
      use const_def
      use crlibm_lib

      implicit none


      ! These values can be saved in photos and restored at restarts
      real(dp) :: R_engulf, Deltar
      ! Note R_engulf is the coordinate of the planet's closest point to the stellar center. Deltar is the infall distance during a timestep


      ! the routines that take care of doing the save/restore are the following:
      ! alloc_extra_info and unpack_extra_info << called by extras_startup
      ! store_extra_info << called by extras_finish_step
      ! these routines call move_extra_info.
      ! it must know about each of your variables to be saved/restored.
      ! so edit move_extra_info when you change the set of variables.

      ! these routines are called by the standard run_star check_model
      contains

      subroutine extras_controls(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         ! this is the place to set any procedure pointers you want to change
         ! e.g., other_wind, other_mixing, other_energy  (see star_data.inc)

         ! Uncomment these lines if you wish to use the functions in this file,
         ! otherwise we use a null_ version which does nothing.
         s% other_energy => energy_routine

         s% extras_startup => extras_startup
         s% extras_start_step => extras_start_step
         s% extras_check_model => extras_check_model
         s% extras_finish_step => extras_finish_step
         s% extras_after_evolve => extras_after_evolve
         s% how_many_extra_history_columns => how_many_extra_history_columns
         s% data_for_extra_history_columns => data_for_extra_history_columns
         s% how_many_extra_profile_columns => how_many_extra_profile_columns
         s% data_for_extra_profile_columns => data_for_extra_profile_columns

         s% how_many_extra_history_header_items => how_many_extra_history_header_items
         s% data_for_extra_history_header_items => data_for_extra_history_header_items
         s% how_many_extra_profile_header_items => how_many_extra_profile_header_items
         s% data_for_extra_profile_header_items => data_for_extra_profile_header_items

         ! Once you have set the function pointers you want,
         ! then uncomment this (or set it in your star_job inlist)
         ! to disable the printed warning message,
         !  s% job% warn_run_star_extras =.false.

      end subroutine extras_controls

      subroutine energy_routine(id, ierr)
        integer, intent(in) :: id
        integer, intent(out) :: ierr
        logical :: restart, first
        type (star_info), pointer :: s
        integer :: k, nz
        integer :: krr_center, krr_bottom, krr_top
        real(dp) :: e_orbit, m_planet, r_planet, area, de
        real(dp) :: rr, v_kepler, rho_bar
        real(dp) :: dmsum, de_heat
        real(dp) :: f_disruption
        ierr = 0

        ! Reads model infos from star structure s. Initialize variables.
        call star_ptr(id, s, ierr)
        if (ierr /= 0) return

        nz = s% nz               ! Mesh size (primary)
        s% extra_heat(:) = 0     ! Initialize extra_heat vector

        ! Initialize injected energy and radial coordinate change
        de = 0d0
        deltar = 0d0

        ! Mass and radius of injested planet

        m_planet = s% x_ctrl(1) * Msun ! Read values of Planet Mass and Radius from inlist ( e.g. 0.001*Msun , 0.1*Rsun)
        r_planet = s% x_ctrl(2) * Rsun !


        ! Find gridpoint corresponding to the location of the engulfed planet bottom, center and top
        ! Not much happening here...
          krr_bottom=1
          do while (krr_bottom >= 1 .and. krr_bottom < nz .and. s% r(krr_bottom) >= R_engulf)
               krr_bottom = krr_bottom + 1
          end do
          krr_center=krr_bottom
            do while (krr_center >=2 .and. s% r(krr_center) < R_engulf+r_planet)
                krr_center = krr_center - 1
            end do
          krr_top=krr_center
          do while (krr_top >= 2 .and. s% r(krr_top) < R_engulf+2*r_planet)
              krr_top = krr_top - 1
          end do

          ! Calculate mass contained in the spherical shell occupied by the planet (shellular approximation)
          dmsum = sum(s% dm(krr_top:krr_bottom))
          rho_bar = dot_product(s% rho(krr_top:krr_bottom), s% dm(krr_top:krr_bottom))/dmsum ! Mass weighted density of the region of impact for drag calculation

          ! Calculate orbital keplerian velocity of the planet (we assume circular orbits)
          call orbital_velocity(s% m(krr_center), R_engulf+r_planet, v_kepler)

          ! Check if the planet has been destroyed by ram pressure (f>1)
          f_disruption = check_disruption(m_planet,r_planet,v_kepler,rho_bar)

          ! Calculate area for drag
          area = 0d0 ! Initialize cross section of planet for calculating aerodynamic drag
          ! Do the calculation only if this is a grazing collision and if the planet has not been destroyed yet
          if ((R_engulf <= s% r(1)) .and. (R_engulf>=(s% r(1) - 2d0*r_planet)) .and. (f_disruption <= 1d0)) then
              ! Calculate intersected area. Rstar-rr is x in sketch
              area = intercepted_area(s% r(1)-R_engulf, r_planet)
              write(*,*) 'Grazing Collision. Engulfed area fraction: ', area/(pi * pow2(r_planet))
          else
              ! Full engulfment. Cross section area = Planet area
              area = pi * pow2(r_planet)
              write(*,*) 'Full engulfment'
          end if


          ! If the planet has not been destroyed by ram pressure, deposit drag luminosity and heat the envelope
          ! Spread in the region occupied by the planet. Update radial coordinate of the engulfed planet too.
          if ( f_disruption <= 1d0 ) then
              ! Note we use s% r(krr_bottom)+r_planet instead of s% r(krr_center) because during grazing phase krr_center = krr_bottom
              call drag( s% m(krr_center), m_planet, area, rho_bar, s% dt, R_engulf+r_planet, de, Deltar)

              ! If the planet has not been destroyed by ram pressure, deposit drag energy and heat the envelope
              ! Spread in the region occupied by the planet.
              do k=krr_top, krr_bottom
                  s% extra_heat(k) = (de/dmsum/s% dt) ! Uniform heating (erg/g/sec)
              end do
            ! Calculate orbital energy
             call orbital_energy(s% m(krr_center), m_planet, R_engulf+r_planet, e_orbit)
             write(*,*) 'Injected Energy / Orbital Energy: ', abs(de/e_orbit)
          else
            write(*,*) '******************** Planet destroyed at R/Rsun = ', R_engulf/Rsun,'*********************'
            Deltar = 0d0
            s% use_other_energy = .false.
          endif

        ! Save variables for history

         s% xtra1 = v_kepler/1d5                ! Orbital velocity
         s% xtra2 = Deltar                      ! Infall distance
         s% xtra3 = de                          ! Injected energy
         s% xtra4 = f_disruption                ! Disruption factor
         s% xtra5 = area/(pi * pow2(r_planet))  ! Engulfed fraction
         s% xtra6 = dmsum/Msun                  ! Heated mass (Msun)

      end subroutine energy_routine

            ! Useful subroutines and functions
            ! Calculate Change in radial position and energy loss due to drag
            subroutine drag(m1, m2, area, rho, dt, r, de, dr)
                 real(dp) :: cdr, cde
                 real(dp), intent(in) ::  m1, m2, area, rho, dt, r
                 real(dp), intent(out) :: de, dr
                 ! use const_def, only: standard_cgrav
                 ! Calculate Deltar (infall distance due to aerodynamic drag)
                 ! We consider cross section area of planet. For compact objects (e.g. NS) one needs to use the accretion radius instead.
                 ! See e.g. equation B.3 in Tylenda & Soker 2006
                 cdr = area*sqrt(standard_cgrav*m1) / m2
                 dr = cdr * rho * sqrt(r) * dt
                 ! Calculate DeltaE (energy loss)
                 ! See e.g. equation B.2 in Tylenda & Soker 2006, where we have used v = keplerian velocity
                 cde = (0.5d0*area)*(standard_cgrav*m1)**1.5d0
                 de = cde *  rho * r**(-1.5d0) * dt
            end subroutine drag

            ! Calculate orbital velocity
            subroutine orbital_velocity(m1, r, v_kepler)
                real(dp), intent(in)  :: m1, r
                real(dp), intent(out) :: v_kepler
                !use const_def, only: standard_cgrav
                v_kepler = sqrt(standard_cgrav*m1/r)
            end subroutine orbital_velocity

            ! Calculate binary orbital energy
            subroutine orbital_energy(m1, m2, r, energy)
                real(dp), intent(in) ::  m1, m2,  r
                real(dp), intent(out) :: energy
                !use const_def, only: standard_cgrav
                energy = -standard_cgrav*m1*m2/(2d0*r)
            end subroutine orbital_energy

            ! Calculate 2D Intercepted area of planet grazing host star
            ! r_planet = planet radius, x = Rstar-R_engulf (see sketch)

            real(dp) function intercepted_area(x, r_planet) result(area)
                real(dp), intent(in) :: x, r_planet
                real(dp) :: alpha, y
                ! Case when less than half of the planet is engulfed
                if (x < r_planet) then
                    y = r_planet - x
                    alpha = acos (y/r_planet)
                    area = r_planet * (r_planet*alpha - y * sin(alpha))
               else
                ! Case when more than half of the planet is engulfed
                   y = x - r_planet
                   alpha = acos (y/r_planet)
                   area = pi*pow2(r_planet) - r_planet * (r_planet*alpha - (y * sin(alpha)))
               endif
            end function intercepted_area

            real(dp) function check_disruption(m_planet,r_planet,v_planet,rho_ambient) result(f)
                  ! f > 1 means disruption. This is expected when the ram pressure integrated
                  ! over the planet cross section approaches the planet binding energy
                 real(dp), intent(in) :: m_planet,r_planet,v_planet,rho_ambient
                 real(dp) :: v_esc_planet_square, rho_planet
                 rho_planet = 3d0*m_planet/(4d0*pi*pow3(r_planet))
                 v_esc_planet_square = standard_cgrav*m_planet/r_planet
                 ! Eq.5 in Jia & Spruit 2018  https://arxiv.org/abs/1808.00467
                 f = (rho_ambient*pow2(v_planet)) / (rho_planet*v_esc_planet_square)
            end function check_disruption


      integer function extras_startup(id, restart, ierr)
         integer, intent(in) :: id
         logical, intent(in) :: restart
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)

         if (ierr /= 0) return
         extras_startup = 0
         if (.not. restart) then
            ! If this is not a restart, set the collision radius as the stellar radius at the beginning of calculation
            ! aka a grazing collision (R_engulf is the coordinate of the planet's closest point to the stellar center)
            R_engulf = s% r(1) !
            call alloc_extra_info(s)
         else ! it is a restart -> Unpack value of R_engulf from photo
            call unpack_extra_info(s)
         end if
      end function extras_startup


      integer function extras_start_step(id, id_extra)
         integer, intent(in) :: id, id_extra
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         extras_start_step = 0
      end function extras_start_step


      ! returns either keep_going, retry, backup, or terminate.
      integer function extras_check_model(id, id_extra)
         integer, intent(in) :: id, id_extra
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_check_model = keep_going

         ! if you want to check multiple conditions, it can be useful
         ! to set a different termination code depending on which
         ! condition was triggered.  MESA provides 9 customizeable
         ! termination codes, named t_xtra1 .. t_xtra9.  You can
         ! customize the messages that will be printed upon exit by
         ! setting the corresponding termination_code_str value.
         ! termination_code_str(t_xtra1) = 'my termination condition'

         ! by default, indicate where (in the code) MESA terminated
         if (extras_check_model == terminate) s% termination_code = t_extras_check_model
      end function extras_check_model


      integer function how_many_extra_history_columns(id, id_extra)
         integer, intent(in) :: id, id_extra
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_history_columns = 9
      end function how_many_extra_history_columns


      subroutine data_for_extra_history_columns(id, id_extra, n, names, vals, ierr)
         integer, intent(in) :: id, id_extra, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         names(1) = 'R_Engulfed_Planet' ! Radial distance from stellar center of engulfed planet
         names(2) = 'Orbital_velocity' ! v_kepler
         names(3) = 'Log_Infall_distance'  ! dr
         names(4) = 'Log_Injected_energy'  ! de=dl*dt
         names(5) = 'Log_Destruction_factor'  ! Eq. 5 from Jia & Spruit 2018
         names(6) = 'Engulfed_fraction'  ! Cross section of the planet engulfed in the star
         names(7) = 'Total_mass_affected'
         names(8) = 'Planet_mass'
         names(9) = 'Planet_radius'
         vals(1) = R_engulf / Rsun
         vals(2) = s% xtra1                 ! Orbital velocity
         vals(3) = safe_log10_cr( s% xtra2) ! Infall distance
         vals(4) = safe_log10_cr( s% xtra3) ! Injected energy
         vals(5) = safe_log10_cr( s% xtra4) ! Disruption factor (ratio between ram pressure and binding energy density)
         vals(6) = s% xtra5                 ! Engulfed fraction
         vals(7) = s% xtra6
         vals(8) = s% x_ctrl(1)
         vals(9) = s% x_ctrl(2)
         ! note: do NOT add the extras names to history_columns.list
         ! the history_columns.list is only for the built-in log column options.
         ! it must not include the new column names you are adding here.
      end subroutine data_for_extra_history_columns


      integer function how_many_extra_profile_columns(id, id_extra)
         use star_def, only: star_info
         integer, intent(in) :: id, id_extra
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_profile_columns = 1
      end function how_many_extra_profile_columns


      subroutine data_for_extra_profile_columns(id, id_extra, n, nz, names, vals, ierr)
         use star_def, only: star_info, maxlen_profile_column_name
         use const_def, only: dp
         integer, intent(in) :: id, id_extra, n, nz
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(nz,n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         ! Adding extra heating so that we can plot on pgstar
         names(1) = 'engulfment_heating'
         do k = 1, nz
            vals(k,1) =  safe_log10_cr(s% extra_heat(k))
         end do

         !note: do NOT add the extra names to profile_columns.list
         ! the profile_columns.list is only for the built-in profile column options.
         ! it must not include the new column names you are adding here.
      end subroutine data_for_extra_profile_columns

      subroutine how_many_extra_history_header_items(id, id_extra, num_cols)
         integer, intent(in) :: id, id_extra
         integer, intent(out) :: num_cols
         num_cols=0
      end subroutine how_many_extra_history_header_items

      subroutine data_for_extra_history_header_items( &
                  id, id_extra, num_extra_header_items, &
                  extra_header_item_names, extra_header_item_vals, ierr)
         integer, intent(in) :: id, id_extra, num_extra_header_items
         character (len=*), pointer :: extra_header_item_names(:)
         real(dp), pointer :: extra_header_item_vals(:)
         type(star_info), pointer :: s
         integer, intent(out) :: ierr
         ierr = 0
         call star_ptr(id,s,ierr)
         if(ierr/=0) return

         !here is an example for adding an extra history header item
         !set num_cols=1 in how_many_extra_history_header_items and then unccomment these lines
         !extra_header_item_names(1) = 'mixing_length_alpha'
         !extra_header_item_vals(1) = s% mixing_length_alpha
      end subroutine data_for_extra_history_header_items


      subroutine how_many_extra_profile_header_items(id, id_extra, num_cols)
          integer, intent(in) :: id, id_extra
          integer, intent(out) :: num_cols
          num_cols = 0
      end subroutine how_many_extra_profile_header_items

      subroutine data_for_extra_profile_header_items( &
                  id, id_extra, num_extra_header_items, &
                  extra_header_item_names, extra_header_item_vals, ierr)
          integer, intent(in) :: id, id_extra, num_extra_header_items
          character (len=*), pointer :: extra_header_item_names(:)
          real(dp), pointer :: extra_header_item_vals(:)
          type(star_info), pointer :: s
          integer, intent(out) :: ierr
          ierr = 0
          call star_ptr(id,s,ierr)
          if(ierr/=0) return

          !here is an exam_planetle for adding an extra profile header item
          !set num_cols=1 in how_many_extra_profile_header_items and then unccomment these lines
          !extra_header_item_names(1) = 'mixing_length_alpha'
          !extra_header_item_vals(1) = s% mixing_length_alpha
      end subroutine data_for_extra_profile_header_items


      ! returns either keep_going or terminate.
      ! note: cannot request retry or backup; extras_check_model can do that.

      integer function extras_finish_step(id, id_extra)
         integer, intent(in) :: id, id_extra
         integer :: ierr, k, num
         logical :: grazing_phase
         real(dp) :: dr, delta_e, area, energy
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_finish_step = keep_going
         call store_extra_info(s)


         ! ##########  Determine next timestep dt_next so that the planet infall distance dr is not too large compared to the planet's radius
         delta_e = 0d0
         dr = 0d0
         energy = 0d0
         grazing_phase = .false.

         ! Only do this if the planet is still around and falling into the star
         if (s% xtra4 < 1d0 .and. Deltar > 0d0) then
           ! Calculate approx gridpoint location of planet center
           k=1
           do while (s% r(k) > (R_engulf+s% x_ctrl(2) * Rsun ))
             k=k+1
           end do
           ! Calculate predicted dr in two cases: grazing phase and full engulfment
           if (R_engulf >= (s% r(1) - 2d0*s% x_ctrl(2) * Rsun)) then
             area = intercepted_area(s% r(1)-R_engulf, s% x_ctrl(2) * Rsun)
             grazing_phase = .true.
           else
             area = pi * pow2(s% x_ctrl(2) * Rsun)
           end if
           ! Estimate dr
           call drag(s% m(k), s% x_ctrl(1)*Msun,area,s% rho(k),s% dt_next,s% r(k),delta_e,dr)

           if (grazing_phase) then                       ! Grazing Phase (requires small dr)
             do while (dr/(s% x_ctrl(2) * Rsun) > 0.05)
               s% dt_next = s% dt_next/2d0              ! There are better strategies, but this is simple enough
               call drag(s% m(k), s% x_ctrl(1)*Msun,area,s% rho(k),s% dt_next,s% r(k),delta_e,dr)
               write(*,*) 'Grazing Dr/r_p too large: ', dr/(s% x_ctrl(2) * Rsun),'Decreasing timestep to ', s% dt_next
             end do
           else
             do while (dr/(s% x_ctrl(2) * Rsun) > 1.0)    ! or Full Engulfment (allow for larger dr)
               s% dt_next = s% dt_next/2d0
               call drag(s% m(k), s% x_ctrl(1)*Msun,area,s% rho(k),s% dt_next,s% r(k),delta_e,dr)
               write(*,*) 'Engulfed Dr/r_p too large: ', dr/(s% x_ctrl(2) * Rsun),'Decreasing timestep to ', s% dt_next
             end do
           end if
         end if




         ! to save a profile,
            ! s% need_to_save_profiles_now = .true.
         ! to update the star log,
            ! s% need_to_update_history_now = .true.

         ! see extras_check_model for information about custom termination codes
         ! by default, indicate where (in the code) MESA terminated
         if (extras_finish_step == terminate) s% termination_code = t_extras_finish_step

       contains

       real(dp) function engulfment_dt(dr, rho_bar, v , r , m_star, m_planet, area) result(dt)
          real(dp), intent(in) :: dr, rho_bar, v , r , m_star, m_planet, area
          real(dp) :: cdr
          cdr = area*sqrt(standard_cgrav*m_star) / m_planet
          dt = dr / (cdr * rho_bar * sqrt(r))
       end function engulfment_dt


      end function extras_finish_step


      subroutine extras_after_evolve(id, id_extra, ierr)
         integer, intent(in) :: id, id_extra
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
      end subroutine extras_after_evolve




      ! routines for saving and restoring extra data so can do restarts

         ! put these defs at the top and delete from the following routines
         !integer, parameter :: extra_info_alloc = 1
         !integer, parameter :: extra_info_get = 2
         !integer, parameter :: extra_info_put = 3


      subroutine alloc_extra_info(s)
         integer, parameter :: extra_info_alloc = 1
         type (star_info), pointer :: s
         call move_extra_info(s,extra_info_alloc)
      end subroutine alloc_extra_info


      subroutine unpack_extra_info(s)
         integer, parameter :: extra_info_get = 2
         type (star_info), pointer :: s
         call move_extra_info(s,extra_info_get)
      end subroutine unpack_extra_info


      subroutine store_extra_info(s)
         integer, parameter :: extra_info_put = 3
         type (star_info), pointer :: s
         call move_extra_info(s,extra_info_put)
      end subroutine store_extra_info


      subroutine move_extra_info(s,op)
         integer, parameter :: extra_info_alloc = 1
         integer, parameter :: extra_info_get = 2
         integer, parameter :: extra_info_put = 3
         type (star_info), pointer :: s
         integer, intent(in) :: op

         integer :: i, j, num_ints, num_dbls, ierr

         i = 0
         ! call move_int or move_flg
         num_ints = i

         i = 0
         ! call move_dbl
         ! (TAs) Important to understand what this is and what is done here. This is essential for MESA to remember R_engulf between timesteps
         ! and to allow for restarts from photos
         call move_dbl(R_engulf)
         num_dbls = i

         if (op /= extra_info_alloc) return
         if (num_ints == 0 .and. num_dbls == 0) return

         ierr = 0
         call star_alloc_extras(s% id, num_ints, num_dbls, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in star_alloc_extras'
            write(*,*) 'alloc_extras num_ints', num_ints
            write(*,*) 'alloc_extras num_dbls', num_dbls
            stop 1
         end if

         contains

         subroutine move_dbl(dbl)
            real(dp) :: dbl
            i = i+1
            select case (op)
            case (extra_info_get)
               dbl = s% extra_work(i)
            case (extra_info_put)
               s% extra_work(i) = dbl
            end select
         end subroutine move_dbl

         subroutine move_int(int)
            integer :: int
            i = i+1
            select case (op)
            case (extra_info_get)
               int = s% extra_iwork(i)
            case (extra_info_put)
               s% extra_iwork(i) = int
            end select
         end subroutine move_int

         subroutine move_flg(flg)
            logical :: flg
            i = i+1
            select case (op)
            case (extra_info_get)
               flg = (s% extra_iwork(i) /= 0)
            case (extra_info_put)
               if (flg) then
                  s% extra_iwork(i) = 1
               else
                  s% extra_iwork(i) = 0
               end if
            end select
         end subroutine move_flg

      end subroutine move_extra_info

      end module run_star_extras

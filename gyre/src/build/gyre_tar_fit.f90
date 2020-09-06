!fpx3_header(0.13_3a)
!
!dependencies
!   dir: ~/Documents/MESA_SS/gyre/src/build 
!   sources: -
!   includes: ../extern/core/core.inc
!   uses: core_parallel core_constants ISO_FORTRAN_ENV core_hgroup gyre_cheb_fit core_kinds gyre_tar_eigen
!   provides: gyre_tar_fit
!end dependencies
!
!end fpx3_header
! Program  : gyre_tar_fit
! Purpose  : fits to traditional approximation of rotation (TAR) eigenvalues
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

module gyre_tar_fit

  ! Uses

  use core_kinds

  use core_hgroup

  use core_parallel
  use core_constants

  use gyre_tar_eigen
  use gyre_cheb_fit

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type :: tar_fit_t
     private
     type(cheb_fit_t) :: cf
     real(WP)         :: nu_0
     integer, public  :: m
     integer, public  :: k
   contains
     private
     procedure       :: lambda_r_
     procedure       :: lambda_c_
     generic, public :: lambda => lambda_r_, lambda_c_
  end type tar_fit_t

  ! Interfaces

  interface tar_fit_t
     module procedure tar_fit_t_tol_
  end interface tar_fit_t

  interface read
     module procedure read_
  end interface read
  interface write
     module procedure write_
  end interface write

  ! Access specifiers

  private

  public :: tar_fit_t

  public :: read
  public :: write

  ! Procedures

contains

  function tar_fit_t_tol_ (m, k, cheb_tol) result (tf)

    integer, intent(in)  :: m
    integer, intent(in)  :: k
    real(WP), intent(in) :: cheb_tol
    type(tar_fit_t)      :: tf

    integer :: l

    ! Construct the tar_fit_t with the specified tolerances

    tf%m = m
    tf%k = k

    if (k >= 0) then

       ! Gravito-inertial waves

       tf%nu_0 = 0._WP

       tf%cf = cheb_fit_t(-1._WP, 1._WP, cheb_tol, f_grav_)

    else

       ! Rossby waves

       l = ABS(m) + ABS(k) - 1

       tf%nu_0 = -REAL(l*(l+1), WP)/REAL(m, WP)

       if (m > 0) then

          tf%cf = cheb_fit_t(-1._WP, 0._WP, cheb_tol, f_ross_)

       elseif (m < 0) then

          tf%cf = cheb_fit_t(0._WP, 1._WP, cheb_tol, f_ross_)

       else

    write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 124 <gyre_tar_fit:tar_fit_t_tol_>:'
    write(UNIT=ERROR_UNIT, FMT=*) 'Invalid m for Rossby waves'

  stop 'Program aborted'

       endif

    endif

    ! Finish

    return

  contains

    function f_grav_ (x) result (f)

      real(WP), intent(in) :: x
      real(WP)             :: f

      ! Calculate the gravity-wave eigenvalue function

      if (m == 0 .AND. k == 0) then

         f = 1._WP

      else

         if (x == -1._WP) then

            if (m <= 0) then
               if (k == 0) then
                  f = 0._WP
               else
                  f = (2._WP*(k-1) + 1._WP)**2
               endif
            elseif (m > 0) then
               f = (2._WP*(k+1) + 1._WP)**2
            endif

         elseif (x == 1._WP) then

            if (m < 0) then
               f = (2._WP*(k+1) + 1._WP)**2
            else
               if (k == 0) then
                  f = 0._WP
               else
                  f = (2._WP*(k-1) + 1._WP)**2
               endif
            endif

         else

            associate (nu => TAN(HALFPI*x))
              f = lambda(nu, m, k)/lambda_norm_grav_(nu, m, k)
            end associate

         endif

      endif

      ! Finish

      return

    end function f_grav_

    function f_ross_ (x) result (f)

      real(WP), intent(in) :: x
      real(WP)             :: f

      ! Calculate the Rossby-wave eigenvaule function

      if (x == 0._WP) then

         f = 0._WP

      elseif (m > 0 .AND. x == -1._WP) then

         f = 1._WP

      elseif (m < 0 .AND. x == 1._WP) then

         f = 1._WP

      else

         associate (nu => TAN(HALFPI*x) + tf%nu_0)
              f = lambda(nu, m, k)/lambda_norm_ross_(nu, m, k)
         end associate

      endif

    end function f_ross_

  end function tar_fit_t_tol_

  !****

  subroutine read_ (hg, tf)

    type(hgroup_t), intent(inout) :: hg
    type(tar_fit_t), intent(out)  :: tf

    type(hgroup_t) :: hg_comp

    ! Read the tar_fit_t

    call read_attr(hg, 'm', tf%m)
    call read_attr(hg, 'k', tf%k)

    call read_attr(hg, 'nu_0', tf%nu_0)

    hg_comp = hgroup_t(hg, 'cf')
    call read(hg_comp, tf%cf)
    call hg_comp%final()

    ! Finish

    return

  end subroutine read_

  !****

  subroutine write_ (hg, tf)

    type(hgroup_t), intent(inout) :: hg
    type(tar_fit_t), intent(in)   :: tf

    type(hgroup_t) :: hg_comp

    ! Write the tar_fit_t

    call write_attr(hg, 'm', tf%m)
    call write_attr(hg, 'k', tf%k)

    call write_attr(hg, 'nu_0', tf%nu_0)

    hg_comp = hgroup_t(hg, 'cf')
    call write(hg_comp, tf%cf)
    call hg_comp%final()

    ! Finish

    return

  end subroutine write_

  !****

  function lambda_r_ (this, nu) result (lambda)

    class(tar_fit_t), intent(in), target :: this
    real(WP), intent(in)                 :: nu
    real(WP)                             :: lambda

    ! Evaluate the eigenvalue of Laplace's tidal equation (real)

    if (this%k >= 0) then

       ! Gravity waves

       associate (x => ATAN(nu)/HALFPI)
         lambda = this%cf%eval(x)*lambda_norm_grav_(nu, this%m, this%k)
       end associate

    else

       ! Rossby waves

       associate (x => ATAN(nu - this%nu_0)/HALFPI)

         if (this%m > 0) then

    if(.NOT. (nu <= this%nu_0)) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''nu <= this%nu_0'' failed at line 301 <gyre_tar_fit:lambda_r_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Invalid nu for Rossby waves'
      stop
    endif

         elseif (this%m < 0) then

    if(.NOT. (nu >= this%nu_0)) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''nu >= this%nu_0'' failed at line 303 <gyre_tar_fit:lambda_r_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Invalid nu for Rossby waves'
      stop
    endif

         else

    write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 305 <gyre_tar_fit:lambda_r_>:'
    write(UNIT=ERROR_UNIT, FMT=*) 'Invalid m for Rossby waves'

  stop 'Program aborted'

         endif

         lambda = this%cf%eval(x)*lambda_norm_ross_(nu, this%m, this%k)

       end associate

    end if

    ! Finish

    return

  end function lambda_r_

  !****

  function lambda_c_ (this, nu) result (lambda)

    class(tar_fit_t), intent(in) :: this
    complex(WP), intent(in)      :: nu
    complex(WP)                  :: lambda

    ! Evaluate the eigenvalue of Laplace's tidal equation (complex)

    if (this%k >= 0) then

       ! Gravity waves

       associate (x => ATAN(nu)/HALFPI)
         lambda = this%cf%eval(x)*lambda_norm_grav_(REAL(nu), this%m, this%k)
       end associate

    else

       ! Rossby waves

       associate (x => ATAN(nu-this%nu_0)/HALFPI)

         if (this%m > 0) then

    if(.NOT. (REAL(nu) <= this%nu_0)) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''REAL(nu) <= this%nu_0'' failed at line 345 <gyre_tar_fit:lambda_c_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Invalid nu for Rossby waves'
      stop
    endif

         elseif (this%m < 0) then

    if(.NOT. (REAL(nu) >= this%nu_0)) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''REAL(nu) >= this%nu_0'' failed at line 347 <gyre_tar_fit:lambda_c_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Invalid nu for Rossby waves'
      stop
    endif

         else

    write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 349 <gyre_tar_fit:lambda_c_>:'
    write(UNIT=ERROR_UNIT, FMT=*) 'Invalid m for Rossby waves'

  stop 'Program aborted'

         endif

         lambda = this%cf%eval(x)*lambda_norm_ross_(REAL(nu), this%m, this%k)

       end associate

    end if

    ! Finish

    return

  end function lambda_c_

  !****

  function lambda_norm_grav_ (nu, m, k) result (lambda_norm)

    real(WP), intent(in) :: nu
    integer, intent(in)  :: m
    integer, intent(in)  :: k
    real(WP)             :: lambda_norm

    integer :: l

    ! Evaluate the gravity-wave eigenvalue normalization function

    if (k >= 0) then

       l = ABS(m) + k

       lambda_norm = nu**2 + l*(l+1)

    else

    write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 385 <gyre_tar_fit:lambda_norm_grav_>:'
    write(UNIT=ERROR_UNIT, FMT=*) 'Invalid k for gravity waves'

  stop 'Program aborted'

    endif

    ! Finish

    return

  end function lambda_norm_grav_

  !****

  function lambda_norm_ross_ (nu, m, k) result (lambda_norm)

    real(WP), intent(in) :: nu
    integer, intent(in)  :: m
    integer, intent(in)  :: k
    real(WP)             :: lambda_norm

    integer :: s

    ! Evaluate the Rossby-wave eigenvalue normalization function

    if (m*nu < 0._WP) then

       if (k < -1) then
          s = -k -1
          lambda_norm = REAL(m, WP)**2/(2*s+1)**2
       elseif (k == -1) then
          lambda_norm = nu**2
       else

    write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 416 <gyre_tar_fit:lambda_norm_ross_>:'
    write(UNIT=ERROR_UNIT, FMT=*) 'Invalid k for Rossby waves'

  stop 'Program aborted'

       endif

    else

    write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 421 <gyre_tar_fit:lambda_norm_ross_>:'
    write(UNIT=ERROR_UNIT, FMT=*) 'Invalid m*nu for Rossby waves'

  stop 'Program aborted'

    endif

    ! Finish

    return

  end function lambda_norm_ross_

end module gyre_tar_fit

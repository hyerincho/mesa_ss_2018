!fpx3_header(0.13_3a)
!
!dependencies
!   dir: ~/Documents/MESA_SS/gyre/src/build 
!   sources: -
!   includes: ../extern/core/core_parallel.inc ../extern/core/core.inc ../interp/gyre_interp.inc
!   uses: core_parallel core_kinds core_linalg ISO_FORTRAN_ENV core_order gyre_calc core_hgroup
!   provides: gyre_r_interp
!end dependencies
!
!end fpx3_header
! Module   : gyre_r_interp
! Purpose  : piecewise cubic interpolators (real)
!
! Copyright 2015-2016 Rich Townsend
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

! Module   : gyre_interp
! Purpose  : piecewise cubic interpolation (template)
!
! Copyright 2015-2017 Rich Townsend
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

module gyre_r_interp

  ! Uses

  use core_kinds

  use core_hgroup

  use core_linalg
  use core_parallel
  use core_order

  use gyre_calc

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type r_interp_t
     private
     real(WP), allocatable  :: x_(:)     ! Abscissa
     real(WP), allocatable :: f_(:)     ! Ordinate
     real(WP), allocatable :: df_dx_(:) ! First derivatives
     integer                :: n         ! Number of points
   contains
     private
     procedure         :: x_n_
     generic, public   :: x => x_n_
     procedure, public :: x_min => x_min_
     procedure, public :: x_max => x_max_
     procedure         :: f_1_
     procedure         :: f_v_
     procedure         :: f_n_
     generic, public   :: f => f_1_, f_v_, f_n_
     procedure         :: df_dx_1_
     procedure         :: df_dx_v_
     procedure         :: df_dx_n_
     generic, public   :: df_dx => df_dx_1_, df_dx_v_, df_dx_n_
     procedure         :: int_f_n_
     generic, public   :: int_f => int_f_n_
  end type r_interp_t

  ! Interfaces

  interface r_interp_t
     module procedure r_interp_t_
     module procedure r_interp_t_eval_derivs_
  end interface r_interp_t

  interface read
     module procedure read_
  end interface read
  interface write
     module procedure write_
  end interface write

  ! Access specifiers

  private

  public :: r_interp_t

  public :: read
  public :: write

  ! Procedures

contains

  function r_interp_t_ (x, f, df_dx) result (in)

    real(WP), intent(in)  :: x(:)
    real(WP), intent(in) :: f(:)
    real(WP), intent(in) :: df_dx(:)
    type(r_interp_t)   :: in

    if(.NOT. (ALL(x(2:) > x(:SIZE(x)-1)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(x(2:) > x(:SIZE(x)-1))'' failed at line 20 <gyre_r_interp:r_interp_t_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Non-monotonic abscissa'
      stop
    endif

    ! Construct the interp_t

    in%x_ = x
    in%f_ = f

    in%df_dx_ = df_dx

    in%n = SIZE(x)

    ! Finish

    return

  end function r_interp_t_

  !****

  function r_interp_t_eval_derivs_ (x, f, deriv_type, df_dx_a, df_dx_b) result (in)

    real(WP), intent(in)            :: x(:)
    real(WP), intent(in)           :: f(:)
    character(*), intent(in)        :: deriv_type
    real(WP), optional, intent(in) :: df_dx_a
    real(WP), optional, intent(in) :: df_dx_b
    type(r_interp_t)             :: in

    ! Construct the interp_t, with derivatives calculated according to
    ! deriv_type

    in = r_interp_t(x, f, deriv(x, f, deriv_type, df_dx_a, df_dx_b))

    ! Finish

    return

  end function r_interp_t_eval_derivs_

  !****

  subroutine read_ (hg, in)

    type(hgroup_t), intent(inout)    :: hg
    type(r_interp_t), intent(out) :: in

    real(WP), allocatable  :: x(:)
    real(WP), allocatable :: f(:)
    real(WP), allocatable :: df_dx(:)

    ! Read the interp_t

    call read_dset_alloc(hg, 'x', x)
    call read_dset_alloc(hg, 'f', f)
    call read_dset_alloc(hg, 'df_dx', df_dx)

    in = r_interp_t(x, f, df_dx)

    ! Finish

    return

  end subroutine read_

  !****

  subroutine write_ (hg, in)

    type(hgroup_t), intent(inout)   :: hg
    type(r_interp_t), intent(in) :: in

    ! Write the interp_t

    call write_attr(hg, 'n', in%n)

    call write_dset(hg, 'x', in%x_)
    call write_dset(hg, 'f', in%f_)
    call write_dset(hg, 'df_dx', in%df_dx_)

    ! Finish

    return

  end subroutine write_

  !****

  !****

  function x_n_ (this) result (x)

    class(r_interp_t), intent(in) :: this
    real(WP)                         :: x(this%n)

    ! Return the abscissa points

    x = this%x_

    return

  end function x_n_

  !****

  function x_min_ (this) result (x_min)

    class(r_interp_t), intent(in) :: this
    real(WP)                         :: x_min

    ! Return the minimum abscissa point

    x_min = this%x_(1)

    ! FInish

    return

  end function x_min_

  !****

  function x_max_ (this) result (x_max)

    class(r_interp_t), intent(in) :: this
    real(WP)                         :: x_max

    ! Return the maximum abscissa point

    x_max = this%x_(this%n)

    ! FInish

    return

  end function x_max_

  !****

  function f_1_ (this, x) result (f)

    class(r_interp_t), intent(in) :: this
    real(WP), intent(in)             :: x
    real(WP)                        :: f

    integer  :: i
    real(WP) :: h
    real(WP) :: w

    ! Interpolate f at a single point

    ! Set up the bracketing index

    i = 1

    call locate(this%x_, x, i)

    if(.NOT. (i > 0 .AND. i < this%n)) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''i > 0 .AND. i < this%n'' failed at line 20 <gyre_r_interp:f_1_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Out-of-bounds interpolation'
      stop
    endif

    ! Set up the interpolation weights

    h = this%x_(i+1) - this%x_(i)
    w = (x - this%x_(i))/h

    ! Do the interpolation

    f =    this%f_(i  )*phi_(1._WP-w) + &
           this%f_(i+1)*phi_(w      ) - &
     h*this%df_dx_(i  )*psi_(1._WP-w) + &
     h*this%df_dx_(i+1)*psi_(w      )

    ! Finish

    return

  end function f_1_

  !****

  function f_v_ (this, x) result (f)

    class(r_interp_t), intent(in) :: this
    real(WP), intent(in)             :: x(:)
    real(WP)                        :: f(SIZE(x))

    integer  :: i
    integer  :: j
    real(WP) :: h
    real(WP) :: w

    ! Interpolate f at a vector of points

    i = 1

    x_loop : do j = 1,SIZE(x)

       ! Update the bracketing index

       call locate(this%x_, x(j), i)

    if(.NOT. (i > 0 .AND. i < this%n)) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''i > 0 .AND. i < this%n'' failed at line 20 <gyre_r_interp:f_v_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Out-of-bounds interpolation'
      stop
    endif

       ! Set up the interpolation weights

       h = this%x_(i+1) - this%x_(i)
       w = (x(j) - this%x_(i))/h

       ! Do the interpolation

       f(j) = this%f_(i  )*phi_(1._WP-w) + &
              this%f_(i+1)*phi_(w      ) - &
        h*this%df_dx_(i  )*psi_(1._WP-w) + &
        h*this%df_dx_(i+1)*psi_(w      )

    end do x_loop

    ! Finish

    return

  end function f_v_

  !****

  function f_n_ (this) result (f)

    class(r_interp_t), intent(in) :: this
    real(WP)                        :: f(this%n)

    ! Return f at abscissa points

    f = this%f_

    ! Finish

  end function f_n_

  !****

  function df_dx_1_ (this, x) result (df_dx)

    class(r_interp_t), intent(in) :: this
    real(WP), intent(in)             :: x
    real(WP)                        :: df_dx

    integer  :: i
    real(WP) :: h
    real(WP) :: w

    ! Differentiate f at a single point

    ! Set up the bracketing index

    i = 1

    call locate(this%x_, x, i)

    if(.NOT. (i > 0 .AND. i < this%n)) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''i > 0 .AND. i < this%n'' failed at line 20 <gyre_r_interp:df_dx_1_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Out-of-bounds interpolation'
      stop
    endif

    ! Set up the interpolation weights

    h = this%x_(i+1) - this%x_(i)
    w = (x - this%x_(i))/h

    ! Do the interpolation

    df_dx =    -this%f_(i  )*dphi_dt_(1._WP-w)/h + &
                this%f_(i+1)*dphi_dt_(w      )/h + &
            this%df_dx_(i  )*dpsi_dt_(1._WP-w) + &
            this%df_dx_(i+1)*dpsi_dt_(w      )

    ! Finish

    return

  end function df_dx_1_

  !****

  function df_dx_v_ (this, x) result (df_dx)

    class(r_interp_t), intent(in) :: this
    real(WP), intent(in)             :: x(:)
    real(WP)                        :: df_dx(SIZE(x))

    integer  :: i
    integer  :: j
    real(WP) :: h
    real(WP) :: w

    ! Differentiate f at a vector of points

    i = 1

    x_loop : do j = 1,SIZE(x)

       ! Update the bracketing index

       call locate(this%x_, x(j), i)

    if(.NOT. (i > 0 .AND. i < this%n)) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''i > 0 .AND. i < this%n'' failed at line 20 <gyre_r_interp:df_dx_v_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Out-of-bounds interpolation'
      stop
    endif

       ! Set up the interpolation weights

       h = this%x_(i+1) - this%x_(i)
       w = (x(j) - this%x_(i))/h

       ! Do the interpolation

       df_dx(j) = -this%f_(i  )*dphi_dt_(1._WP-w)/h + &
                   this%f_(i+1)*dphi_dt_(w      )/h + &
               this%df_dx_(i  )*dpsi_dt_(1._WP-w) + &
               this%df_dx_(i+1)*dpsi_dt_(w      )

    end do x_loop

    ! Finish

    return

  end function df_dx_v_

  !****

  function df_dx_n_ (this) result (df_dx)

    class(r_interp_t), intent(in) :: this
    real(WP)                   :: df_dx(this%n)

    ! Return df_dx at abscissa points

    df_dx = this%df_dx_

    ! Finish

    return

  end function df_dx_n_

  !****

  function int_f_n_ (this) result (int_f)

    class(r_interp_t), intent(in) :: this
    real(WP)                        :: int_f(this%n)

    integer  :: i
    real(WP) :: h

    ! Integrate f across the full domain

    int_f(1) = 0._WP

    x_loop : do i = 1,this%n-1

       h = this%x_(i+1) - this%x_(i)

       int_f(i+1) = int_f(i) + (this%f_(i) + this%f_(i+1))*h/2._WP - &
                               (this%df_dx_(i+1) - this%df_dx_(i))*h**2/12._WP

    end do x_loop

    ! Finish

    return

  end function int_f_n_

  !****

  function phi_ (t)

    real(WP), intent(in) :: t
    real(WP)             :: phi_

    ! Evaluate the phi basis function

    phi_ = 3._WP*t**2 - 2._WP*t**3

    return

  end function phi_

  !****

  function psi_ (t)

    real(WP), intent(in) :: t
    real(WP)             :: psi_

    ! Evaluate the psi basis function 

    psi_ = t**3 - t**2

    return

  end function psi_

  !****

  function dphi_dt_ (t)

    real(WP), intent(in) :: t
    real(WP)             :: dphi_dt_

    ! Evaluate the first derivative of the phi basis function

    dphi_dt_ = 6._WP*t - 6._WP*t**2

    return

  end function dphi_dt_

  !****

  function dpsi_dt_ (t)

    real(WP), intent(in) :: t
    real(WP)             :: dpsi_dt_

    ! Evaluate the first derivative of the psi basis function

    dpsi_dt_ = 3._WP*t**2 - 2._WP*t

    return

  end function dpsi_dt_

end module gyre_r_interp


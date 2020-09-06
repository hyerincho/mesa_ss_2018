!fpx3_header(0.13_3a)
!
!dependencies
!   dir: ~/Documents/MESA_SS/gyre/src/build 
!   sources: -
!   includes: ../extern/core/core_parallel.inc ../extern/core/core.inc
!   uses: core_parallel ISO_FORTRAN_ENV core_kinds gyre_r_ext
!   provides: gyre_c_ext
!end dependencies
!
!end fpx3_header
! Module   : gyre_c_ext
! Purpose  : extented-range arithmetic (complex)
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

module gyre_c_ext

  ! Uses

  use core_kinds
  use core_parallel

  use gyre_r_ext

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Parameter definitions

  real(WP), parameter :: RADIX_WP = REAL(RADIX(1._WP), WP)

  ! Derived-type definitions

  type c_ext_t
     private
     complex(WP) :: f ! Fractional part
     integer     :: e ! Exponent
   contains

    procedure             :: op_plus_cx_cx_
    procedure             :: op_plus_cx_rx_
    procedure, pass(that) :: op_plus_rx_cx_
    procedure             :: op_plus_cx_r_
    procedure, pass(that) :: op_plus_r_cx_
    procedure             :: op_plus_cx_c_
    procedure, pass(that) :: op_plus_c_cx_
    generic, public       :: operator(+) => op_plus_cx_cx_, op_plus_cx_rx_, op_plus_rx_cx_, &
                                              op_plus_cx_r_, op_plus_r_cx_, &
                                              op_plus_cx_c_, op_plus_c_cx_

     procedure       :: op_minus_cx_
     generic, public :: operator(-) => op_minus_cx_

    procedure             :: op_minus_cx_cx_
    procedure             :: op_minus_cx_rx_
    procedure, pass(that) :: op_minus_rx_cx_
    procedure             :: op_minus_cx_r_
    procedure, pass(that) :: op_minus_r_cx_
    procedure             :: op_minus_cx_c_
    procedure, pass(that) :: op_minus_c_cx_
    generic, public       :: operator(-) => op_minus_cx_cx_, op_minus_cx_rx_, op_minus_rx_cx_, &
                                              op_minus_cx_r_, op_minus_r_cx_, &
                                              op_minus_cx_c_, op_minus_c_cx_

    procedure             :: op_times_cx_cx_
    procedure             :: op_times_cx_rx_
    procedure, pass(that) :: op_times_rx_cx_
    procedure             :: op_times_cx_r_
    procedure, pass(that) :: op_times_r_cx_
    procedure             :: op_times_cx_c_
    procedure, pass(that) :: op_times_c_cx_
    generic, public       :: operator(*) => op_times_cx_cx_, op_times_cx_rx_, op_times_rx_cx_, &
                                              op_times_cx_r_, op_times_r_cx_, &
                                              op_times_cx_c_, op_times_c_cx_

    procedure             :: op_divide_cx_cx_
    procedure             :: op_divide_cx_rx_
    procedure, pass(that) :: op_divide_rx_cx_
    procedure             :: op_divide_cx_r_
    procedure, pass(that) :: op_divide_r_cx_
    procedure             :: op_divide_cx_c_
    procedure, pass(that) :: op_divide_c_cx_
    generic, public       :: operator(/) => op_divide_cx_cx_, op_divide_cx_rx_, op_divide_rx_cx_, &
                                              op_divide_cx_r_, op_divide_r_cx_, &
                                              op_divide_cx_c_, op_divide_c_cx_

    procedure             :: op_eq_cx_cx_
    procedure             :: op_eq_cx_rx_
    procedure, pass(that) :: op_eq_rx_cx_
    procedure             :: op_eq_cx_r_
    procedure, pass(that) :: op_eq_r_cx_
    procedure             :: op_eq_cx_c_
    procedure, pass(that) :: op_eq_c_cx_
    generic, public       :: operator(==) => op_eq_cx_cx_, op_eq_cx_rx_, op_eq_rx_cx_, &
                                              op_eq_cx_r_, op_eq_r_cx_, &
                                              op_eq_cx_c_, op_eq_c_cx_

    procedure             :: op_neq_cx_cx_
    procedure             :: op_neq_cx_rx_
    procedure, pass(that) :: op_neq_rx_cx_
    procedure             :: op_neq_cx_r_
    procedure, pass(that) :: op_neq_r_cx_
    procedure             :: op_neq_cx_c_
    procedure, pass(that) :: op_neq_c_cx_
    generic, public       :: operator(/=) => op_neq_cx_cx_, op_neq_cx_rx_, op_neq_rx_cx_, &
                                              op_neq_cx_r_, op_neq_r_cx_, &
                                              op_neq_cx_c_, op_neq_c_cx_

  end type c_ext_t

  ! Interface blocks

  interface c_ext_t
     module procedure c_ext_t_r_
     module procedure c_ext_t_r_i_
     module procedure c_ext_t_c_
     module procedure c_ext_t_rx_
     module procedure c_ext_t_rx_ix_
  end interface c_ext_t

  interface r_ext_t
     module procedure r_ext_t_cx_
  end interface r_ext_t

  interface real
     module procedure real_
  end interface real

  interface cmplx
     module procedure cmplx_
  end interface cmplx

  interface valid
     module procedure valid_
  end interface valid

  interface conjg
     module procedure conjg_
  end interface conjg

  interface real_part
     module procedure real_part_
  end interface real_part

  interface imag_part
     module procedure imag_part_
  end interface imag_part

  interface product
     module procedure product_
  end interface product

  interface exp
     module procedure exp_
  end interface exp

  interface sqrt
     module procedure sqrt_
  end interface sqrt

  interface abs
     module procedure abs_
  end interface abs

  interface fraction
     module procedure fraction_
  end interface fraction

  interface exponent
     module procedure exponent_
  end interface exponent

  interface scale
     module procedure scale_
  end interface scale

  ! Access specifiers

  private

  public :: c_ext_t
  public :: r_ext_t
  public :: real
  public :: cmplx
  public :: real_part
  public :: imag_part
  public :: valid
  public :: conjg
  public :: product
  public :: abs
  public :: exp
  public :: sqrt
  public :: fraction
  public :: exponent
  public :: scale

  ! Procedures

contains

  elemental function c_ext_t_r_ (r) result (cx)

    real(WP), intent(in) :: r
    type(c_ext_t)        :: cx

    ! Construct the c_ext_t from the real r

    call split_(CMPLX(r, KIND=WP), cx%f, cx%e)

    ! Finish

    return

  end function c_ext_t_r_

  !****

  elemental function c_ext_t_r_i_ (r, i) result (cx)

    real(WP), intent(in) :: r
    real(WP), intent(in) :: i
    type(c_ext_t)        :: cx

    ! Construct the c_ext_t from the real/imaginary pair r, i

    call split_(CMPLX(r, i, KIND=WP), cx%f, cx%e)

    ! Finish

    return

  end function c_ext_t_r_i_

  !****

  elemental function c_ext_t_c_ (c) result (cx)

    complex(WP), intent(in) :: c
    type(c_ext_t)           :: cx

    ! Construct the c_ext_t from the complex c

    call split_(c, cx%f, cx%e)

    ! Finish

    return

  end function c_ext_t_c_

  !****

  elemental function c_ext_t_rx_ (rx) result (cx)

    type(r_ext_t), intent(in) :: rx
    type(c_ext_t)             :: cx

    ! Construct the c_ext_t from the r_ext_t rx

    call split_(CMPLX(FRACTION(rx), KIND=WP), cx%f, cx%e)
    cx = scale(cx, EXPONENT(rx))

    ! Finish

    return

  end function c_ext_t_rx_

  !****

  elemental function c_ext_t_rx_ix_ (rx, ix) result (cx)

    type(r_ext_t), intent(in) :: rx
    type(r_ext_t), intent(in) :: ix
    type(c_ext_t)             :: cx

    ! Construct the c_ext_t from the r_ext_t real/imaginary pair rx,
    ! ix

    cx = c_ext_t(rx) + c_ext_t(0._WP, 1._WP)*c_ext_t(ix)

    ! Finish

    return

  end function c_ext_t_rx_ix_

  !****

  elemental function op_plus_cx_cx_ (this, that) result (cx)

    class(c_ext_t), intent(in) :: this
    class(c_ext_t), intent(in) :: that
    type(c_ext_t)              :: cx

    complex(WP) :: f
    integer     :: e

    ! Evaluate the plus operator

    if (this%f == 0._WP) then

       cx%f = that%f
       cx%e = that%e

    elseif (that%f == 0._WP) then

       cx%f = this%f
       cx%e = this%e

    else

       if (this%e > that%e) then
          f = this%f + cmplx(c_ext_t(that%f, that%e - this%e))
          e = this%e
       else
          f = cmplx(c_ext_t(this%f, this%e - that%e)) + that%f
          e = that%e
       endif

       call split_(f, cx%f, cx%e)
       cx%e = cx%e + e

    endif

    ! Finish

    return

  end function op_plus_cx_cx_

  !****

  elemental function op_minus_cx_ (this) result (cx)

    class(c_ext_t), intent(in) :: this
    type(c_ext_t)              :: cx

    ! Evaluate the unary minus operator

    cx%f = -this%f
    cx%e = this%e

    ! Finish

    return

  end function op_minus_cx_

  !****

  elemental function op_minus_cx_cx_ (this, that) result (cx)

    class(c_ext_t), intent(in) :: this
    class(c_ext_t), intent(in) :: that
    type(c_ext_t)              :: cx

    complex(WP) :: f
    integer     :: e

    ! Evaluate the minus operator

    if (this%f == 0._WP) then

       cx%f = -that%f
       cx%e = that%e

    elseif (that%f == 0._WP) then

       cx%f = this%f
       cx%e = this%e

    else

       if (this%e > that%e) then
          f = this%f - cmplx(c_ext_t(that%f, that%e - this%e))
          e = this%e
       else
          f = cmplx(c_ext_t(this%f, this%e - that%e)) - that%f
          e = that%e
       endif

       call split_(f, cx%f, cx%e)
       cx = scale(cx, e)

    endif

    ! Finish

    return

  end function op_minus_cx_cx_

  !****

  elemental function op_times_cx_cx_ (this, that) result (cx)

    class(c_ext_t), intent(in) :: this
    class(c_ext_t), intent(in) :: that
    type(c_ext_t)              :: cx

    complex(WP) :: f
    integer     :: e

    ! Evaluate the times operator

    if (this%f == 0._WP .OR. that%f == 0._WP) then

       cx = c_ext_t(0._WP)

    else

       f = this%f*that%f
       e = this%e + that%e

       call split_(f, cx%f, cx%e)
       cx = scale(cx, e)

    endif

    ! Finish

    return

  end function op_times_cx_cx_

  !****

  elemental function op_divide_cx_cx_ (this, that) result (cx)

    class(c_ext_t), intent(in) :: this
    class(c_ext_t), intent(in) :: that
    type(c_ext_t)              :: cx

    complex(WP) :: f
    integer     :: e

    ! Evaluate the divide operator

    if (this%f == 0._WP .AND. that%f /= 0._WP) then

       cx = c_ext_t(0._WP)

    else

       f = this%f/that%f
       e = this%e - that%e

       call split_(f, cx%f, cx%e)
       cx = scale(cx, e)

    endif

    ! Finish

    return

  end function op_divide_cx_cx_

  !****

  elemental function op_eq_cx_cx_ (this, that) result (eq)

    class(c_ext_t), intent(in) :: this
    class(c_ext_t), intent(in) :: that
    logical                    :: eq

    ! Evaluate the equality operator

    eq = this%e == that%e .AND. this%f == that%f

    ! Finish

    return

  end function op_eq_cx_cx_

  !****

  elemental function op_neq_cx_cx_ (this, that) result (neq)

    class(c_ext_t), intent(in) :: this
    class(c_ext_t), intent(in) :: that
    logical                    :: neq

    ! Evaluate the inequality operator

    neq = .NOT. this == that

    ! Finish

    return

  end function op_neq_cx_cx_

  !****

  elemental function op_plus_cx_rx_ (this, that) result (plus)

    class(c_ext_t), intent(in) :: this
    class(r_ext_t), intent(in) :: that
    type(c_ext_t)                      :: plus

    ! Evaluate the operator to mixed c_ext_t/r_ext_t types

    plus = this + c_ext_t(that)

    ! Finish

    return

  end function op_plus_cx_rx_

  elemental function op_plus_rx_cx_ (this, that) result (plus)

    class(r_ext_t), intent(in) :: this
    class(c_ext_t), intent(in) :: that
    type(c_ext_t)                      :: plus

    ! Evaluate the operator to mixed r_ext_t/c_ext_t types

    plus = c_ext_t(this) + that

    ! Finish

    return

  end function op_plus_rx_cx_

  elemental function op_plus_cx_c_ (this, that) result (plus)

    class(c_ext_t), intent(in) :: this
    complex(WP), intent(in)    :: that
    type(c_ext_t)                      :: plus

    ! Evaluate the operator to mixed c_ext_t/complex types

    plus = this + c_ext_t(that)

    ! Finish

    return

  end function op_plus_cx_c_

  elemental function op_plus_c_cx_ (this, that) result (plus)

    complex(WP), intent(in)    :: this
    class(c_ext_t), intent(in) :: that
    type(c_ext_t)                      :: plus

    ! Evaluate the operator to mixed complex/c_ext_t types

    plus = c_ext_t(this) + that

    ! Finish

    return

  end function op_plus_c_cx_

  elemental function op_plus_cx_r_ (this, that) result (plus)

    class(c_ext_t), intent(in) :: this
    real(WP), intent(in)       :: that
    type(c_ext_t)                      :: plus

    ! Evaluate the operator to mixed c_ext_t/real types

    plus = this + c_ext_t(that)

    ! Finish

    return

  end function op_plus_cx_r_

  elemental function op_plus_r_cx_ (this, that) result (plus)

    real(WP), intent(in)       :: this
    class(c_ext_t), intent(in) :: that
    type(c_ext_t)                      :: plus

    ! Evaluate the operator to mixed real/c_ext_t types

    plus = c_ext_t(this) + that

    ! Finish

    return

  end function op_plus_r_cx_

  elemental function op_minus_cx_rx_ (this, that) result (minus)

    class(c_ext_t), intent(in) :: this
    class(r_ext_t), intent(in) :: that
    type(c_ext_t)                      :: minus

    ! Evaluate the operator to mixed c_ext_t/r_ext_t types

    minus = this - c_ext_t(that)

    ! Finish

    return

  end function op_minus_cx_rx_

  elemental function op_minus_rx_cx_ (this, that) result (minus)

    class(r_ext_t), intent(in) :: this
    class(c_ext_t), intent(in) :: that
    type(c_ext_t)                      :: minus

    ! Evaluate the operator to mixed r_ext_t/c_ext_t types

    minus = c_ext_t(this) - that

    ! Finish

    return

  end function op_minus_rx_cx_

  elemental function op_minus_cx_c_ (this, that) result (minus)

    class(c_ext_t), intent(in) :: this
    complex(WP), intent(in)    :: that
    type(c_ext_t)                      :: minus

    ! Evaluate the operator to mixed c_ext_t/complex types

    minus = this - c_ext_t(that)

    ! Finish

    return

  end function op_minus_cx_c_

  elemental function op_minus_c_cx_ (this, that) result (minus)

    complex(WP), intent(in)    :: this
    class(c_ext_t), intent(in) :: that
    type(c_ext_t)                      :: minus

    ! Evaluate the operator to mixed complex/c_ext_t types

    minus = c_ext_t(this) - that

    ! Finish

    return

  end function op_minus_c_cx_

  elemental function op_minus_cx_r_ (this, that) result (minus)

    class(c_ext_t), intent(in) :: this
    real(WP), intent(in)       :: that
    type(c_ext_t)                      :: minus

    ! Evaluate the operator to mixed c_ext_t/real types

    minus = this - c_ext_t(that)

    ! Finish

    return

  end function op_minus_cx_r_

  elemental function op_minus_r_cx_ (this, that) result (minus)

    real(WP), intent(in)       :: this
    class(c_ext_t), intent(in) :: that
    type(c_ext_t)                      :: minus

    ! Evaluate the operator to mixed real/c_ext_t types

    minus = c_ext_t(this) - that

    ! Finish

    return

  end function op_minus_r_cx_

  elemental function op_times_cx_rx_ (this, that) result (times)

    class(c_ext_t), intent(in) :: this
    class(r_ext_t), intent(in) :: that
    type(c_ext_t)                      :: times

    ! Evaluate the operator to mixed c_ext_t/r_ext_t types

    times = this * c_ext_t(that)

    ! Finish

    return

  end function op_times_cx_rx_

  elemental function op_times_rx_cx_ (this, that) result (times)

    class(r_ext_t), intent(in) :: this
    class(c_ext_t), intent(in) :: that
    type(c_ext_t)                      :: times

    ! Evaluate the operator to mixed r_ext_t/c_ext_t types

    times = c_ext_t(this) * that

    ! Finish

    return

  end function op_times_rx_cx_

  elemental function op_times_cx_c_ (this, that) result (times)

    class(c_ext_t), intent(in) :: this
    complex(WP), intent(in)    :: that
    type(c_ext_t)                      :: times

    ! Evaluate the operator to mixed c_ext_t/complex types

    times = this * c_ext_t(that)

    ! Finish

    return

  end function op_times_cx_c_

  elemental function op_times_c_cx_ (this, that) result (times)

    complex(WP), intent(in)    :: this
    class(c_ext_t), intent(in) :: that
    type(c_ext_t)                      :: times

    ! Evaluate the operator to mixed complex/c_ext_t types

    times = c_ext_t(this) * that

    ! Finish

    return

  end function op_times_c_cx_

  elemental function op_times_cx_r_ (this, that) result (times)

    class(c_ext_t), intent(in) :: this
    real(WP), intent(in)       :: that
    type(c_ext_t)                      :: times

    ! Evaluate the operator to mixed c_ext_t/real types

    times = this * c_ext_t(that)

    ! Finish

    return

  end function op_times_cx_r_

  elemental function op_times_r_cx_ (this, that) result (times)

    real(WP), intent(in)       :: this
    class(c_ext_t), intent(in) :: that
    type(c_ext_t)                      :: times

    ! Evaluate the operator to mixed real/c_ext_t types

    times = c_ext_t(this) * that

    ! Finish

    return

  end function op_times_r_cx_

  elemental function op_divide_cx_rx_ (this, that) result (divide)

    class(c_ext_t), intent(in) :: this
    class(r_ext_t), intent(in) :: that
    type(c_ext_t)                      :: divide

    ! Evaluate the operator to mixed c_ext_t/r_ext_t types

    divide = this / c_ext_t(that)

    ! Finish

    return

  end function op_divide_cx_rx_

  elemental function op_divide_rx_cx_ (this, that) result (divide)

    class(r_ext_t), intent(in) :: this
    class(c_ext_t), intent(in) :: that
    type(c_ext_t)                      :: divide

    ! Evaluate the operator to mixed r_ext_t/c_ext_t types

    divide = c_ext_t(this) / that

    ! Finish

    return

  end function op_divide_rx_cx_

  elemental function op_divide_cx_c_ (this, that) result (divide)

    class(c_ext_t), intent(in) :: this
    complex(WP), intent(in)    :: that
    type(c_ext_t)                      :: divide

    ! Evaluate the operator to mixed c_ext_t/complex types

    divide = this / c_ext_t(that)

    ! Finish

    return

  end function op_divide_cx_c_

  elemental function op_divide_c_cx_ (this, that) result (divide)

    complex(WP), intent(in)    :: this
    class(c_ext_t), intent(in) :: that
    type(c_ext_t)                      :: divide

    ! Evaluate the operator to mixed complex/c_ext_t types

    divide = c_ext_t(this) / that

    ! Finish

    return

  end function op_divide_c_cx_

  elemental function op_divide_cx_r_ (this, that) result (divide)

    class(c_ext_t), intent(in) :: this
    real(WP), intent(in)       :: that
    type(c_ext_t)                      :: divide

    ! Evaluate the operator to mixed c_ext_t/real types

    divide = this / c_ext_t(that)

    ! Finish

    return

  end function op_divide_cx_r_

  elemental function op_divide_r_cx_ (this, that) result (divide)

    real(WP), intent(in)       :: this
    class(c_ext_t), intent(in) :: that
    type(c_ext_t)                      :: divide

    ! Evaluate the operator to mixed real/c_ext_t types

    divide = c_ext_t(this) / that

    ! Finish

    return

  end function op_divide_r_cx_

  elemental function op_eq_cx_rx_ (this, that) result (eq)

    class(c_ext_t), intent(in) :: this
    class(r_ext_t), intent(in) :: that
    logical                      :: eq

    ! Evaluate the operator to mixed c_ext_t/r_ext_t types

    eq = this == c_ext_t(that)

    ! Finish

    return

  end function op_eq_cx_rx_

  elemental function op_eq_rx_cx_ (this, that) result (eq)

    class(r_ext_t), intent(in) :: this
    class(c_ext_t), intent(in) :: that
    logical                      :: eq

    ! Evaluate the operator to mixed r_ext_t/c_ext_t types

    eq = c_ext_t(this) == that

    ! Finish

    return

  end function op_eq_rx_cx_

  elemental function op_eq_cx_c_ (this, that) result (eq)

    class(c_ext_t), intent(in) :: this
    complex(WP), intent(in)    :: that
    logical                      :: eq

    ! Evaluate the operator to mixed c_ext_t/complex types

    eq = this == c_ext_t(that)

    ! Finish

    return

  end function op_eq_cx_c_

  elemental function op_eq_c_cx_ (this, that) result (eq)

    complex(WP), intent(in)    :: this
    class(c_ext_t), intent(in) :: that
    logical                      :: eq

    ! Evaluate the operator to mixed complex/c_ext_t types

    eq = c_ext_t(this) == that

    ! Finish

    return

  end function op_eq_c_cx_

  elemental function op_eq_cx_r_ (this, that) result (eq)

    class(c_ext_t), intent(in) :: this
    real(WP), intent(in)       :: that
    logical                      :: eq

    ! Evaluate the operator to mixed c_ext_t/real types

    eq = this == c_ext_t(that)

    ! Finish

    return

  end function op_eq_cx_r_

  elemental function op_eq_r_cx_ (this, that) result (eq)

    real(WP), intent(in)       :: this
    class(c_ext_t), intent(in) :: that
    logical                      :: eq

    ! Evaluate the operator to mixed real/c_ext_t types

    eq = c_ext_t(this) == that

    ! Finish

    return

  end function op_eq_r_cx_

  elemental function op_neq_cx_rx_ (this, that) result (neq)

    class(c_ext_t), intent(in) :: this
    class(r_ext_t), intent(in) :: that
    logical                      :: neq

    ! Evaluate the operator to mixed c_ext_t/r_ext_t types

    neq = this /= c_ext_t(that)

    ! Finish

    return

  end function op_neq_cx_rx_

  elemental function op_neq_rx_cx_ (this, that) result (neq)

    class(r_ext_t), intent(in) :: this
    class(c_ext_t), intent(in) :: that
    logical                      :: neq

    ! Evaluate the operator to mixed r_ext_t/c_ext_t types

    neq = c_ext_t(this) /= that

    ! Finish

    return

  end function op_neq_rx_cx_

  elemental function op_neq_cx_c_ (this, that) result (neq)

    class(c_ext_t), intent(in) :: this
    complex(WP), intent(in)    :: that
    logical                      :: neq

    ! Evaluate the operator to mixed c_ext_t/complex types

    neq = this /= c_ext_t(that)

    ! Finish

    return

  end function op_neq_cx_c_

  elemental function op_neq_c_cx_ (this, that) result (neq)

    complex(WP), intent(in)    :: this
    class(c_ext_t), intent(in) :: that
    logical                      :: neq

    ! Evaluate the operator to mixed complex/c_ext_t types

    neq = c_ext_t(this) /= that

    ! Finish

    return

  end function op_neq_c_cx_

  elemental function op_neq_cx_r_ (this, that) result (neq)

    class(c_ext_t), intent(in) :: this
    real(WP), intent(in)       :: that
    logical                      :: neq

    ! Evaluate the operator to mixed c_ext_t/real types

    neq = this /= c_ext_t(that)

    ! Finish

    return

  end function op_neq_cx_r_

  elemental function op_neq_r_cx_ (this, that) result (neq)

    real(WP), intent(in)       :: this
    class(c_ext_t), intent(in) :: that
    logical                      :: neq

    ! Evaluate the operator to mixed real/c_ext_t types

    neq = c_ext_t(this) /= that

    ! Finish

    return

  end function op_neq_r_cx_

  !****

  elemental function r_ext_t_cx_ (cx) result (rx)

    type(c_ext_t), intent(in) :: cx
    type(r_ext_t)             :: rx

    ! Construct the r_ext_t from the c_ext_t

    rx = r_ext_t(REAL(cx%f))
    rx = scale(rx, cx%e)

    ! Finish

    return

  end function r_ext_t_cx_

  !****

  elemental function real_ (cx) result (r)

    type(c_ext_t), intent(in) :: cx
    real(WP)                  :: r

    ! Convert cx to real

    r = REAL(cmplx(cx))

    ! Finish

    return

  end function real_

  !****

  elemental function cmplx_ (cx) result (c)

    type(c_ext_t), intent(in) :: cx
    complex(WP)               :: c

    integer :: e_min

    ! Convert cx to complex

    if (cx%f /= 0._WP) then

       e_min = MINEXPONENT(0._WP)

       if (cx%e >= e_min) then
          c = cx%f*RADIX_WP**cx%e
       else
          c = (cx%f*RADIX_WP**MAX(cx%e-e_min, -DIGITS(0._WP)-1))*RADIX_WP**e_min
       endif

    else

       c = 0._WP

    endif

    ! Finish

    return

  end function cmplx_

  !****

  elemental function real_part_ (cx) result (rx)

    type(c_ext_t), intent(in) :: cx
    type(r_ext_t)             :: rx

    ! Extract the real part of cx

    rx = scale(r_ext_t(REAL(cx%f)), cx%e)

    ! Finish

    return

  end function real_part_

  !****

  elemental function imag_part_ (cx) result (rx)

    type(c_ext_t), intent(in) :: cx
    type(r_ext_t)             :: rx

    ! Extract the imaginary part of cx

    rx = scale(r_ext_t(AIMAG(cx%f)), cx%e)

    ! Finish

    return

  end function imag_part_

  !****

  elemental function valid_ (cx) result (valid)

    type(c_ext_t), intent(in) :: cx
    logical                   :: valid

    ! Determine if cx is valid

    valid = ABS(REAL(cx%f)) >= RADIX_WP**(-1) .AND. ABS(REAL(cx%f)) < 1._WP

    ! Finish

    return

  end function valid_

  !****

  elemental function conjg_ (cx) result (conjg_cx)

    type(c_ext_t), intent(in) :: cx
    type(c_ext_t)             :: conjg_cx

    ! Calculate the complex conjugate of cx

    conjg_cx%f = CONJG(cx%f)
    conjg_cx%e = cx%e

    ! Finish

    return

  end function conjg_

  !****

  function product_ (cx) result (prod_cx)

    type(c_ext_t), intent(in) :: cx(:)
    type(c_ext_t)             :: prod_cx

    integer     :: i
    complex(WP) :: f
    integer     :: e

    ! Calculate the product of the elements of cx

    prod_cx%f = 1._WP
    prod_cx%e = SUM(cx%e)

    do i = 1,SIZE(cx)

       call split_(prod_cx%f*cx(i)%f, f, e)

       prod_cx%f = f
       prod_cx = scale(prod_cx, e)

    end do

    ! Finish

    return

  end function product_

  !****

  elemental function abs_ (cx) result (abs_cx)

    type(c_ext_t), intent(in) :: cx
    type(r_ext_t)             :: abs_cx

    ! Calculate the absolute value of cx

    abs_cx = r_ext_t(ABS(cx%f))
    abs_cx = scale(abs_cx, cx%e)

    ! Finish

    return

  end function abs_

  !****

  elemental function exp_ (cx) result (exp_cx)

    type(c_ext_t), intent(in) :: cx
    type(c_ext_t)             :: exp_cx

    type(r_ext_t) :: exp_ex

    ! Calculate the exponential of cx

    exp_ex = exp(r_ext_t(cx))

    exp_cx%f = FRACTION(exp_ex)*EXP((0._WP,1._WP)*AIMAG(cmplx(cx)))
    exp_cx%e = EXPONENT(exp_ex)

    ! Finish

    return

  end function exp_

  !****

  elemental function sqrt_ (cx) result (sqrt_cx)

    type(c_ext_t), intent(in) :: cx
    type(c_ext_t)             :: sqrt_cx

    ! Calculate the square root of cx

    sqrt_cx = c_ext_t(SQRT(cx%f))

    if (MOD(cx%e, 2) == 0) then
       sqrt_cx = scale(sqrt_cx, cx%e/2)
    else
       sqrt_cx = scale(sqrt_cx, (cx%e-1)/2)*SQRT(2._WP)
    endif

    ! Finish

    return

  end function sqrt_

  !****

  elemental function fraction_ (cx) result (fraction_cx)

    type(c_ext_t), intent(in) :: cx
    complex(WP)               :: fraction_cx

    ! Return the fraction part of cx

    fraction_cx = cx%f

    ! Finish

    return

  end function fraction_

  !****

  elemental function exponent_ (cx) result (exponent_cx)

    type(c_ext_t), intent(in) :: cx
    integer                   :: exponent_cx

    ! Return the exponent part of cx

    exponent_cx = cx%e

    ! Finish

    return

  end function exponent_

  !****

  elemental function scale_ (cx, de) result (scale_cx)

    class(c_ext_t), intent(in) :: cx
    integer, intent(in)        :: de
    type(c_ext_t)              :: scale_cx

    ! Scale cx by RADIX_WP**de

    scale_cx%f = cx%f

    if (scale_cx%f /= 0._WP) then
       scale_cx%e = cx%e + de
    else
       scale_cx%e = 0
    endif

    ! Finish

    return

  end function scale_

  !****

  elemental subroutine split_ (c, f, e)

    complex(WP), intent(in)  :: c
    complex(WP), intent(out) :: f
    integer, intent(out)     :: e

    real(WP)    :: c_r
    real(WP)    :: c_i
    real(WP)    :: f_r
    real(WP)    :: f_i
    integer     :: e_r
    integer     :: e_i
    integer     :: e_ref

    ! Spit c into fraction and exponent parts

    c_r = REAL(c)
    c_i = AIMAG(c)

    f_r = FRACTION(c_r)
    f_i = FRACTION(c_i)

    if (f_r == 0._WP .AND. f_i == 0._WP) then

       f = 0._WP
       e = 0

    else

       e_r = EXPONENT(c_r)
       e_i = EXPONENT(c_i)

       e_ref = MAXVAL([e_r,e_i], MASK=[f_r,f_i] /= 0._WP)

       if (f_r /= 0._WP) then
          if (e_r - e_ref < MINEXPONENT(0._WP)) then
             f_r = 0
          else
             f_r = f_r*RADIX_WP**(e_r - e_ref)
          endif
       endif

       if (f_i /= 0._WP) then
          if (e_i - e_ref < MINEXPONENT(0._WP)) then
             f_i = 0
          else
             f_i = f_i*RADIX_WP**(e_i - e_ref)
          endif
       endif

       f = CMPLX(f_r, f_i, WP)
       e = e_ref

    endif

    ! Finish

    return

  end subroutine split_

  !****

end module gyre_c_ext

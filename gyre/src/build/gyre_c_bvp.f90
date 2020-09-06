!fpx3_header(0.13_3a)
!
!dependencies
!   dir: ~/Documents/MESA_SS/gyre/src/build 
!   sources: -
!   includes: ../bvp/gyre_bvp.inc ../extern/core/core.inc
!   uses: gyre_ext gyre_diff core_kinds gyre_num_par gyre_state gyre_sysmtx_factory gyre_bound gyre_status ISO_FORTRAN_ENV gyre_sysmtx
!   provides: gyre_c_bvp
!end dependencies
!
!end fpx3_header
! Module   : gyre_c_bvp
! Purpose  : parametric boundary value problems (complex)
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

! Incfile  : gyre_bvp
! Purpose  : parametric boundary value problems (template)
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

module gyre_c_bvp

  ! Uses

  use core_kinds

  use gyre_bound
  use gyre_diff
  use gyre_ext
  use gyre_num_par
  use gyre_state
  use gyre_status
  use gyre_sysmtx
  use gyre_sysmtx_factory

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type :: c_bvp_t
     private
     class(c_diff_t), allocatable   :: df(:)
     class(c_bound_t), allocatable  :: bd
     class(c_sysmtx_t), allocatable :: sm
     integer, public                   :: n_k
     integer, public                   :: n_e
     integer, public                   :: n_i
     integer, public                   :: n_o
     logical                           :: factored
   contains 
     private
     procedure, public :: build
     procedure, public :: factor
     procedure, public :: det
     procedure, public :: soln_vec_hom
     procedure, public :: soln_vec_inhom
     procedure, public :: resd_vec
  end type c_bvp_t

  ! Interfaces

  interface c_bvp_t
     module procedure c_bvp_t
  end interface c_bvp_t

  ! Access specifiers

  private

  public :: c_bvp_t
  public :: c_bvp_t_

contains

  function c_bvp_t_ (bd, df, nm_p) result (bp)

    class(c_diff_t), intent(in)   :: df(:)
    class(c_bound_t), intent(in)  :: bd
    type(num_par_t), intent(in)      :: nm_p
    type(c_bvp_t)                 :: bp

    integer :: n_k
    integer :: n_e
    integer :: n_i
    integer :: n_o

    ! Perform basic validations

    n_k = SIZE(df) + 1
    n_e = bd%n_e

    n_i = bd%n_i
    n_o = bd%n_o

    ! Construct the bvp_t

    allocate(bp%df(n_k-1), SOURCE=df)

    allocate(bp%bd, SOURCE=bd)

    allocate(bp%sm, SOURCE=c_sysmtx_t(n_k-1, n_e, n_i, n_o, nm_p))

    bp%n_k = n_k
    bp%n_e = n_e
    bp%n_i = n_i
    bp%n_o = n_o

    bp%factored = .FALSE.

    ! Finish

    return

  end function c_bvp_t_

  !****

  subroutine build (this, st)

    class(c_bvp_t), target, intent(inout) :: this
    class(c_state_t), intent(in)          :: st

    complex(WP)        :: B_i(this%n_i,this%n_e)
    complex(WP)        :: B_o(this%n_o,this%n_e)
    complex(WP)        :: E_l(this%n_e,this%n_e)
    complex(WP)        :: E_r(this%n_e,this%n_e)
    complex(WP)        :: scl_i(this%n_i)
    complex(WP)        :: scl_o(this%n_o)
    type(c_ext_t) :: scl
    integer          :: k

    ! Build the bvp for the specified state

    ! Set up boundary conditions

    call this%bd%build_i(st, B_i, scl_i)
    call this%sm%set_B_i(B_i, scl_i)

    call this%bd%build_o(st, B_o, scl_o)
    call this%sm%set_B_o(B_o, scl_o)

    ! Set up difference equations

    !$OMP PARALLEL DO PRIVATE (E_l, E_r, scl) SCHEDULE (DYNAMIC)
    sub_loop : do k = 1, this%n_k-1
       call this%df(k)%build(st, E_l, E_r, scl)
       call this%sm%set_E(k, E_l, E_r, scl)
    end do sub_loop

    ! Reset the factored flag

    this%factored = .FALSE.

    ! Finish

    return

  end subroutine build

  !****

  subroutine factor (this)

    class(c_bvp_t), intent(inout) :: this

    ! Factorize the sysmtx

    call this%sm%factor()

    this%factored = .TRUE.

    ! Finish

    return

  end subroutine factor

  !****

  function det (this)

    class(c_bvp_t), intent(inout) :: this
    type(c_ext_t)                 :: det

    if(.NOT. (this%factored)) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''this%factored'' failed at line 20 <gyre_c_bvp:det>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Matrix has not been factorized'
      stop
    endif

    ! Evaluate the determinant of the sysmtx

    det = this%sm%det()

    ! Finish

    return

  end function det

  !****

  function soln_vec_hom (this) result (y)

    class(c_bvp_t), intent(inout) :: this
    complex(WP)                        :: y(this%n_e,this%n_k)

    complex(WP) :: v(this%n_e*this%n_k)

    ! Evaluate the solution vector y of the homogeneous system

    v = this%sm%soln_vec_hom()

    y = RESHAPE(v, SHAPE(y))

    ! Finish

    return

  end function soln_vec_hom

  !****

  function soln_vec_inhom (this, w_i, w_o) result (y)

    class(c_bvp_t), intent(inout) :: this
    complex(WP), intent(in)            :: w_i(:)
    complex(WP), intent(in)            :: w_o(:)
    complex(WP)                        :: y(this%n_e,this%n_k)

    complex(WP) :: v(this%n_e*this%n_k)

    ! Evaluate the solution vector y of the inhomogeneous system

    v = this%sm%soln_vec_inhom(w_i, w_o)

    y = RESHAPE(v, SHAPE(y))

    ! Finish

    return

  end function soln_vec_inhom

  !****

  function resd_vec (this, y) result (w)

    class(c_bvp_t), intent(inout) :: this
    complex(WP), intent(in)            :: y(:,:)
    complex(WP)                        :: w(this%n_e*this%n_k)

    complex(WP) :: v(this%n_e*this%n_k)

    ! Evaluate the residuals vector

    v = RESHAPE(y, SHAPE(v))

    w = this%sm%resd_vec(v)

    ! Finish

    return

  end function resd_vec

end module gyre_c_bvp


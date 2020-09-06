!fpx3_header(0.13_3a)
!
!dependencies
!   dir: ~/Documents/MESA_SS/gyre/src/build 
!   sources: -
!   includes: ../diff/gyre_magnus_diff.inc ../extern/core/core.inc
!   uses: gyre_eqns gyre_linalg gyre_diff ISO_FORTRAN_ENV core_kinds gyre_ext gyre_point gyre_state
!   provides: gyre_r_magnus_diff
!end dependencies
!
!end fpx3_header
! Module   : gyre_r_magnus_ivp
! Purpose  : difference equations (Magnus method, real)
!
! Copyright 2013-2015 Rich Townsend
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

! Incfile  : gyre_magnus_diff
! Purpose  : difference equations (Magnus method, template)
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

module gyre_r_magnus_diff

  ! Uses

  use core_kinds

  use gyre_diff
  use gyre_eqns
  use gyre_ext
  use gyre_linalg
  use gyre_point
  use gyre_state

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Parameter definitions

  integer, parameter :: GL2_SCHEME = 1
  integer, parameter :: GL4_SCHEME = 2
  integer, parameter :: GL6_SCHEME = 3

  real(WP), parameter :: GL2_C_I(1) = [0.5_WP]

  real(WP), parameter :: GL4_C_I(2) = [ &
       0.5_WP-SQRT(3._WP)/6._WP, &
       0.5_WP+SQRT(3._WP)/6._WP]

  real(WP), parameter :: GL6_C_I(3) = [ &
       0.5_WP-SQRT(15._WP)/10._WP, &
       0.5_WP, &
       0.5_WP+SQRT(15._WP)/10._WP]

  ! Derived-type definitions

  type, extends (r_diff_t) :: r_magnus_diff_t
     private
     class(r_eqns_t), allocatable :: eq
     real(WP)                        :: dx
     integer                         :: scheme
   contains
     private
     procedure, public :: build => build
     procedure         :: dOmega_GL2_
     procedure         :: dOmega_GL4_
     procedure         :: dOmega_GL6_
  end type r_magnus_diff_t

  ! Interfaces

  interface r_magnus_diff_t
     module procedure r_magnus_diff_t_
  end interface r_magnus_diff_t

  ! Access specifiers

  private

  public :: r_magnus_diff_t

  ! Procedures

contains

  function r_magnus_diff_t_ (eq, pt_a, pt_b, scheme) result (df)

    class(r_eqns_t), intent(in) :: eq
    type(point_t), intent(in)      :: pt_a
    type(point_t), intent(in)      :: pt_b
    character(*), intent(in)       :: scheme
    type(r_magnus_diff_t)       :: df

    real(WP), allocatable      :: c_i(:)
    type(point_t), allocatable :: pt(:)

    ! Construct the magnus_diff_t

    select case (scheme)
    case ('GL2')
       df%scheme = GL2_SCHEME
       c_i = GL2_C_I
    case ('GL4')
       df%scheme = GL4_SCHEME
       c_i = GL4_C_I
    case ('GL6')
       df%scheme = GL6_SCHEME
       c_i = GL6_C_I
    case default

    write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 20 <gyre_r_magnus_diff:r_magnus_diff_t_>:'
    write(UNIT=ERROR_UNIT, FMT=*) 'Invalid scheme'

  stop 'Program aborted'

    end select

    df%dx = pt_b%x - pt_a%x

    allocate(pt(SIZE(c_i)))

    pt%s = pt_a%s
    pt%x = pt_a%x + c_i*df%dx

    allocate(df%eq, SOURCE=eq)

    call df%eq%stencil(pt)

    df%n_e = eq%n_e

    ! Finish

    return

  end function r_magnus_diff_t_

  !****

  subroutine build (this, st, E_l, E_r, scl)

    class(r_magnus_diff_t), intent(in) :: this
    class(r_state_t), intent(in)       :: st
    real(WP), intent(out)                :: E_l(:,:)
    real(WP), intent(out)                :: E_r(:,:)
    type(r_ext_t), intent(out)         :: scl

    logical, parameter :: UPWIND = .TRUE.

    real(WP)   :: dOmega(this%n_e,this%n_e)
    complex(WP) :: lambda(this%n_e)
    complex(WP) :: V_l(this%n_e,this%n_e)
    complex(WP) :: V_r(this%n_e,this%n_e)
    integer     :: i
    complex(WP) :: V_pos(this%n_e,this%n_e)
    complex(WP) :: V_neg(this%n_e,this%n_e)
    complex(WP) :: E_l_(this%n_e,this%n_e)
    complex(WP) :: E_r_(this%n_e,this%n_e)

    ! Build the difference equations

    ! Evaluate the Magnus slope matrix

    select case (this%scheme)
    case (GL2_SCHEME)
       dOmega = this%dOmega_GL2_(st)
    case (GL4_SCHEME)
       dOmega = this%dOmega_GL4_(st)
    case (GL6_SCHEME)
       dOmega = this%dOmega_GL6_(st)
    case default

    write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 20 <gyre_r_magnus_diff:build>:'
    write(UNIT=ERROR_UNIT, FMT=*) 'Invalid scheme'

  stop 'Program aborted'

    end select

    ! Decompose it

    call eigen_decompose(dOmega, lambda, V_l, V_r)

    ! Build the difference equations

    if (UPWIND) then

       ! Use 'upwinding' for stability

       do i = 1, this%n_e
          call ZCOPY(this%n_e, V_r(1,i), 1, V_pos(1,i), 1)
          if (REAL(lambda(i)) >= 0._WP) then
             call ZSCAL(this%n_e, EXP(-lambda(i)*this%dx), V_pos(1,i), 1)
          endif
       end do

       do i = 1, this%n_e
          call ZCOPY(this%n_e, V_r(1,i), 1, V_neg(1,i), 1)
          if (REAL(lambda(i)) < 0._WP) then
             call ZSCAL(this%n_e, EXP(lambda(i)*this%dx), V_neg(1,i), 1)
          endif
       end do

       call ZGEMM('N', 'N', this%n_e, this%n_e, this%n_e, CMPLX(-1._WP, KIND=WP), &
                     V_neg, this%n_e, V_l, this%n_e, CMPLX(0._WP, KIND=WP), &
                     E_l_, this%n_e)

       call ZGEMM('N', 'N', this%n_e, this%n_e, this%n_e, CMPLX(1._WP, KIND=WP), &
                     V_pos, this%n_e, V_l, this%n_e, CMPLX(0._WP, KIND=WP), &
                     E_r_, this%n_e)

       scl = exp(r_ext_t(SUM(lambda, MASK=REAL(lambda) >= 0._WP)*this%dx))

    else

       ! Use the [TowTei2013] expression (their eqns. 4, 16, 18)

       do i = 1, this%n_e
          call ZCOPY(this%n_e, V_r(1,i), 1, V_neg(1,i), 1)
          call ZSCAL(this%n_e, EXP(lambda(i)*this%dx), V_neg(1,i), 1)
       end do

       call ZGEMM('N', 'N', this%n_e, this%n_e, this%n_e, CMPLX(1._WP, KIND=WP), &
                     V_neg, this%n_e, V_l, this%n_e, CMPLX(0._WP, KIND=WP), &
                     E_l_, this%n_e)

       do i = 1, this%n_e
          E_r_(:,i) = 0._WP
          E_r_(i,i) = -1._WP
       end do

       scl = r_ext_t(1._WP)

    endif

    E_l = REAL(E_l_)
    E_r = REAL(E_r_)

    ! Finish

    return

  end subroutine build

  !****

  function dOmega_GL2_ (this, st) result (dOmega)

    class(r_magnus_diff_t), intent(in) :: this
    class(r_state_t), intent(in)       :: st
    real(WP)                             :: dOmega(this%n_e,this%n_e)

    real(WP) :: A(this%n_e,this%n_e,1)

    ! Evaluate the GL2 Magnus slope matrix

    ! Calculate the RHS matrix

    A(:,:,1) = this%eq%A(1, st)

    ! Set up the slope matrix

    dOmega = A(:,:,1)

    ! Finish

    return

  end function dOmega_GL2_

  !****

  function dOmega_GL4_ (this, st) result (dOmega)

    class(r_magnus_diff_t), intent(in) :: this
    class(r_state_t), intent(in)       :: st
    real(WP)                             :: dOmega(this%n_e,this%n_e)

    real(WP) :: A(this%n_e,this%n_e,2)
    real(WP) :: dalpha(this%n_e,this%n_e,2)

    ! Evaluate the GL4 Magnus slope matrix

    ! Calculate the RHS matrices

    A(:,:,1) = this%eq%A(1, st)
    A(:,:,2) = this%eq%A(2, st)

    ! Set up the Magnus slope matrix (Blanes et al. 2009, eqns. 243
    ! and 253; note that the 12 in the denominator of their expression
    ! for alpha_2 is erroneous)

    dalpha(:,:,1) = 0.5_WP*(A(:,:,1) + A(:,:,2))
    dalpha(:,:,2) = SQRT(3._WP)*(A(:,:,2) - A(:,:,1))

    dOmega = dalpha(:,:,1) - this%dx*commutator(dalpha(:,:,1), dalpha(:,:,2))/12._WP

    ! Finish

    return

  end function dOmega_GL4_

  !****

  function dOmega_GL6_ (this, st) result (dOmega)

    class(r_magnus_diff_t), intent(in) :: this
    class(r_state_t), intent(in)       :: st
    real(WP)                             :: dOmega(this%n_e,this%n_e)

    real(WP) :: A(this%n_e,this%n_e,3)
    real(WP) :: dalpha(this%n_e,this%n_e,3)
    real(WP) :: dC(this%n_e,this%n_e,2)

    ! Evaluate the GL6 Magnus slope matrix

    ! Calculate the RHS matrices

    A(:,:,1) = this%eq%A(1, st)
    A(:,:,2) = this%eq%A(2, st)
    A(:,:,3) = this%eq%A(3, st)

    ! Set up the Magnus slope matrix (Blanes et al. 2009, eqns. 251
    ! and 257)

    dalpha(:,:,1) = A(:,:,2)
    dalpha(:,:,2) = SQRT(15._WP)*(A(:,:,3) - A(:,:,1))/3
    dalpha(:,:,3) = 10*(A(:,:,3) - 2*A(:,:,2) + A(:,:,1))/3

    dC(:,:,1) = this%dx*commutator(dalpha(:,:,1), dalpha(:,:,2))
    dC(:,:,2) = -this%dx*commutator(dalpha(:,:,1), 2*dalpha(:,:,3)+dC(:,:,1))/60

    dOmega = dalpha(:,:,1) + dalpha(:,:,3)/12 + &
            this%dx*commutator(-20*dalpha(:,:,1)-dalpha(:,:,3)+dC(:,:,1), dalpha(:,:,2)+dC(:,:,2))/240

    ! Finish

    return

  end function dOmega_GL6_

end module gyre_r_magnus_diff


!fpx3_header(0.13_3a)
!
!dependencies
!   dir: ~/Documents/MESA_SS/gyre/src/build 
!   sources: -
!   includes: ../extern/core/core.inc
!   uses: gyre_constants gyre_out_par gyre_num_par gyre_grid_par gyre_osc_par core_parallel ISO_FORTRAN_ENV gyre_mode_par gyre_scan_par core_memory core_kinds
!   provides: gyre_util
!end dependencies
!
!end fpx3_header
! Module   : gyre_util
! Purpose  : miscellaneous utility routines
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

module gyre_util

  ! Uses

  use core_kinds
  use core_parallel
  use core_memory

  use gyre_constants
  use gyre_grid_par
  use gyre_mode_par
  use gyre_num_par
  use gyre_osc_par
  use gyre_out_par
  use gyre_scan_par

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Module variables

  character(64), save :: log_level_m

  ! Interfaces

  interface select_par
     module procedure select_par_gr_1_
     module procedure select_par_nm_1_
     module procedure select_par_os_1_
     module procedure select_par_sc_1_
     module procedure select_par_gr_v_
     module procedure select_par_nm_v_
     module procedure select_par_os_v_
     module procedure select_par_sc_v_
  end interface select_par

  interface sprint
     module procedure sprint_
  end interface sprint

  interface integrate
     module procedure integrate_r_
     module procedure Integrate_c_
  end interface integrate

  interface integral
     module procedure integral_r_
     module procedure Integral_c_
  end interface integral

  ! Access specifiers

  private

  public :: form_header
  public :: set_log_level
  public :: check_log_level
  public :: select_par
  public :: split_list
  public :: join_fmts
  public :: sprint
  public :: rjust
  public :: phase
  public :: integrate
  public :: integral

  ! Procedures

contains

  function form_header (header, underchar)

    character(*), intent(in)           :: header
    character(*), optional, intent(in) :: underchar
    character(:), allocatable          :: form_header

    ! Format the header string

    if(PRESENT(underchar)) then

       if(underchar == '') then

          form_header = TRIM(header) // NEW_LINE('') // &
                        REPEAT(' ', LEN(header)) // NEW_LINE('')

       else

          form_header = TRIM(header) // NEW_LINE('') // &
                        REPEAT(underchar, LEN(header)/LEN(underchar)) // NEW_LINE('')

       endif

    else

       form_header = TRIM(header) // NEW_LINE('')

    endif

    ! Finish

    return

  end function form_header

  !****

  subroutine set_log_level (log_level)

    character(*), intent(in) :: log_level

    ! Set the log level

    select case (log_level)
    case ('DEBUG')
    case ('INFO')
    case ('WARN')
    case default

    write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 140 <gyre_util:set_log_level>:'
    write(UNIT=ERROR_UNIT, FMT=*) 'Invalid log_level'

  stop 'Program aborted'

    end select

    log_level_m = log_level

    ! Finish

    return

  end subroutine set_log_level

  !****

  function check_log_level (log_level, rank)

    character(*), intent(in)      :: log_level
    integer, optional, intent(in) :: rank
    logical                       :: check_log_level

    integer :: rank_

    if(PRESENT(rank)) then
       rank_ = rank
    else
       rank_ = 0
    endif

    ! Check whether we should write log output

    if (MPI_RANK == rank_) then

       select case (log_level)
       case ('DEBUG')
          check_log_level = log_level_m == 'DEBUG'
       case ('INFO')
          check_log_level = log_level_m == 'INFO' .OR. &
                            log_level_m == 'DEBUG'
       case ('WARN')
          check_log_level = log_level_m == 'WARN' .OR. &
                            log_level_m == 'INFO' .OR. &
                            log_level_m == 'DEBUG'
       case default

    write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 182 <gyre_util:check_log_level>:'
    write(UNIT=ERROR_UNIT, FMT=*) 'Invalid log_level'

  stop 'Program aborted'

       end select

    else

       check_log_level = .FALSE.

    endif

    ! Finish

    return

  end function check_log_level

  !****

  subroutine select_par_gr_1_ (par, tag, par_sel)

    type(grid_par_t), intent(in)   :: par(:)
    character(*), intent(in)      :: tag
    type(grid_par_t), intent(out)  :: par_sel

    type(grid_par_t), allocatable :: par_sel_(:)

    ! Select the last parameter whose tag_list matches tag

    call select_par(par, tag, par_sel_)

    par_sel = par_sel_(SIZE(par_sel_))

    ! Finish

    return

  end subroutine select_par_gr_1_

  subroutine select_par_nm_1_ (par, tag, par_sel)

    type(num_par_t), intent(in)   :: par(:)
    character(*), intent(in)      :: tag
    type(num_par_t), intent(out)  :: par_sel

    type(num_par_t), allocatable :: par_sel_(:)

    ! Select the last parameter whose tag_list matches tag

    call select_par(par, tag, par_sel_)

    par_sel = par_sel_(SIZE(par_sel_))

    ! Finish

    return

  end subroutine select_par_nm_1_

  subroutine select_par_os_1_ (par, tag, par_sel)

    type(osc_par_t), intent(in)   :: par(:)
    character(*), intent(in)      :: tag
    type(osc_par_t), intent(out)  :: par_sel

    type(osc_par_t), allocatable :: par_sel_(:)

    ! Select the last parameter whose tag_list matches tag

    call select_par(par, tag, par_sel_)

    par_sel = par_sel_(SIZE(par_sel_))

    ! Finish

    return

  end subroutine select_par_os_1_

  subroutine select_par_sc_1_ (par, tag, par_sel)

    type(scan_par_t), intent(in)   :: par(:)
    character(*), intent(in)      :: tag
    type(scan_par_t), intent(out)  :: par_sel

    type(scan_par_t), allocatable :: par_sel_(:)

    ! Select the last parameter whose tag_list matches tag

    call select_par(par, tag, par_sel_)

    par_sel = par_sel_(SIZE(par_sel_))

    ! Finish

    return

  end subroutine select_par_sc_1_

  !****

  subroutine select_par_gr_v_ (par, tag, par_sel)

    type(grid_par_t), intent(in)               :: par(:)
    character(*), intent(in)                  :: tag
    type(grid_par_t), allocatable, intent(out) :: par_sel(:)

    integer :: i
    logical :: mask(SIZE(par))
    integer :: n_par_sel
    integer :: j

    ! Select all parameters whose $VAR_list matches $VAR

    mask_loop : do i = 1,SIZE(par)
       mask(i) = (par(i)%tag_list == '') .OR. &
                 (tag /= '' .AND. ANY(split_list(par(i)%tag_list, ',') == tag))
    end do mask_loop

    n_par_sel = COUNT(mask)

    if(.NOT. (n_par_sel >= 1)) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''n_par_sel >= 1'' failed at line 281 <gyre_util:select_par_gr_v_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'No matching &grid namelists'
      stop
    endif

    allocate(par_sel(n_par_sel))

    j = 0

    select_loop : do i = 1,SIZE(par)
       if (mask(i)) then
          j = j + 1
          par_sel(j) = par(i)
       endif
    end do select_loop

    ! Finish

    return

  end subroutine select_par_gr_v_

  subroutine select_par_nm_v_ (par, tag, par_sel)

    type(num_par_t), intent(in)               :: par(:)
    character(*), intent(in)                  :: tag
    type(num_par_t), allocatable, intent(out) :: par_sel(:)

    integer :: i
    logical :: mask(SIZE(par))
    integer :: n_par_sel
    integer :: j

    ! Select all parameters whose $VAR_list matches $VAR

    mask_loop : do i = 1,SIZE(par)
       mask(i) = (par(i)%tag_list == '') .OR. &
                 (tag /= '' .AND. ANY(split_list(par(i)%tag_list, ',') == tag))
    end do mask_loop

    n_par_sel = COUNT(mask)

    if(.NOT. (n_par_sel >= 1)) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''n_par_sel >= 1'' failed at line 282 <gyre_util:select_par_nm_v_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'No matching &num namelists'
      stop
    endif

    allocate(par_sel(n_par_sel))

    j = 0

    select_loop : do i = 1,SIZE(par)
       if (mask(i)) then
          j = j + 1
          par_sel(j) = par(i)
       endif
    end do select_loop

    ! Finish

    return

  end subroutine select_par_nm_v_

  subroutine select_par_os_v_ (par, tag, par_sel)

    type(osc_par_t), intent(in)               :: par(:)
    character(*), intent(in)                  :: tag
    type(osc_par_t), allocatable, intent(out) :: par_sel(:)

    integer :: i
    logical :: mask(SIZE(par))
    integer :: n_par_sel
    integer :: j

    ! Select all parameters whose $VAR_list matches $VAR

    mask_loop : do i = 1,SIZE(par)
       mask(i) = (par(i)%tag_list == '') .OR. &
                 (tag /= '' .AND. ANY(split_list(par(i)%tag_list, ',') == tag))
    end do mask_loop

    n_par_sel = COUNT(mask)

    if(.NOT. (n_par_sel >= 1)) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''n_par_sel >= 1'' failed at line 283 <gyre_util:select_par_os_v_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'No matching &osc namelists'
      stop
    endif

    allocate(par_sel(n_par_sel))

    j = 0

    select_loop : do i = 1,SIZE(par)
       if (mask(i)) then
          j = j + 1
          par_sel(j) = par(i)
       endif
    end do select_loop

    ! Finish

    return

  end subroutine select_par_os_v_

  subroutine select_par_sc_v_ (par, tag, par_sel)

    type(scan_par_t), intent(in)               :: par(:)
    character(*), intent(in)                  :: tag
    type(scan_par_t), allocatable, intent(out) :: par_sel(:)

    integer :: i
    logical :: mask(SIZE(par))
    integer :: n_par_sel
    integer :: j

    ! Select all parameters whose $VAR_list matches $VAR

    mask_loop : do i = 1,SIZE(par)
       mask(i) = (par(i)%tag_list == '') .OR. &
                 (tag /= '' .AND. ANY(split_list(par(i)%tag_list, ',') == tag))
    end do mask_loop

    n_par_sel = COUNT(mask)

    if(.NOT. (n_par_sel >= 1)) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''n_par_sel >= 1'' failed at line 284 <gyre_util:select_par_sc_v_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'No matching &scan namelists'
      stop
    endif

    allocate(par_sel(n_par_sel))

    j = 0

    select_loop : do i = 1,SIZE(par)
       if (mask(i)) then
          j = j + 1
          par_sel(j) = par(i)
       endif
    end do select_loop

    ! Finish

    return

  end subroutine select_par_sc_v_

  !****

  function split_list (list, delim) result (elems)

    character(*), intent(in)          :: list
    character(1), intent(in)          :: delim
    character(LEN(list)), allocatable :: elems(:)

    character(LEN(list)) :: list_
    integer              :: d
    integer              :: n
    integer              :: j

    ! Split the delimited list into an array of elements

    d = 16

    allocate(elems(d))

    n = 0

    ! Repeatedly split on delimiters

    list_ = list

    split_loop : do

       if(list_ == '') exit split_loop

       j = INDEX(list_, delim)

       if(j <= 0) then
          n = n + 1
          elems(n) = ADJUSTL(list_)
          exit split_loop
       endif

       n = n + 1

       ! Chop out the element

       elems(n) = ADJUSTL(list_(:j-1))
       list_ = list_(j+1:)

       ! If necessary, expand the array

       if(n >= d) then
          d = 2*d
          call reallocate(elems, [d])
       end if

    end do split_loop

    ! Reallocate elems to the correct length

    call reallocate(elems, [n])

    ! Finish

    return

  end function split_list

!****

  function join_fmts (fmts, n) result (fmt)

    character(*), intent(in)  :: fmts(:)
    integer, intent(in)       :: n(:)
    character(:), allocatable :: fmt

    integer :: i

    ! Join format strings with the appropriate repeat counts

    if(SUM(n) > 0) then

       do i = 1, SIZE(fmts)

          if(ALLOCATED(fmt)) then
             fmt = fmt//','//sprint(n(i))//fmts(i)
          else
             fmt = sprint(n(i))//fmts(i)
          endif

       end do

    else

       fmt = ''

    endif

    ! Add wrap-around parens

    fmt = '('//fmt//')'

    ! Finish

    return

  end function join_fmts

!****

  function sprint_ (i) result (a)

    integer, intent(in)       :: i
    character(:), allocatable :: a

    integer :: n

    ! Print an integer into a character

    ! First, determine the length

    if(i > 0) then
       n = FLOOR(LOG10(REAL(i))) + 1
    elseif(i < 0) then
       n = FLOOR(LOG10(REAL(ABS(i)))) + 2
    else
       n = 1
    endif

    allocate(character(n)::a)

    ! Do the conversion

    write(a, 100) i
100 format(I0)

    ! Finish

    return

  end function sprint_

!****

  function rjust (a, n) result (a_just)

    character(*), intent(in) :: a
    integer, intent(in)      :: n
    character(n)             :: a_just

    ! Right-justify a in a field width of n

    a_just = REPEAT(' ', MAX(n-LEN_TRIM(a), 0))//a

    ! Finish

    return

  end function rjust

!****

  function phase (z)

    complex(WP), intent(in) :: z
    real(WP)                :: phase

    ! Calculate the phase (in radians) of the complex number z

    phase = ATAN2(AIMAG(z), REAL(z))

    ! Finish

    return

  end function phase

!****

  function integrate_r_ (x, y, mask) result (int_y)

    real(WP), intent(in)          :: x(:)
    real(WP), intent(in)         :: y(:)
    logical, optional, intent(in) :: mask(:)
    real(WP)                     :: int_y

    integer :: n

    if (PRESENT(mask)) then

    endif

    ! Integrate y(x) using trapezoidal quadrature, applying the
    ! optional mask

    n = SIZE(x)

    if (PRESENT(mask)) then
       int_y = SUM(0.5_WP*(y(2:) + y(:n-1))*(x(2:) - x(:n-1)), MASK=mask(2:) .AND. mask(:n-1))
    else
       int_y = SUM(0.5_WP*(y(2:) + y(:n-1))*(x(2:) - x(:n-1)))
    endif

    ! Finish

    return

  end function integrate_r_

  function integrate_c_ (x, y, mask) result (int_y)

    real(WP), intent(in)          :: x(:)
    complex(WP), intent(in)         :: y(:)
    logical, optional, intent(in) :: mask(:)
    complex(WP)                     :: int_y

    integer :: n

    if (PRESENT(mask)) then

    endif

    ! Integrate y(x) using trapezoidal quadrature, applying the
    ! optional mask

    n = SIZE(x)

    if (PRESENT(mask)) then
       int_y = SUM(0.5_WP*(y(2:) + y(:n-1))*(x(2:) - x(:n-1)), MASK=mask(2:) .AND. mask(:n-1))
    else
       int_y = SUM(0.5_WP*(y(2:) + y(:n-1))*(x(2:) - x(:n-1)))
    endif

    ! Finish

    return

  end function integrate_c_

!****

  function integral_r_ (x, y) result (int_y)

    real(WP), intent(in)  :: x(:)
    real(WP), intent(in) :: y(:)
    real(WP)             :: int_y(SIZE(x))

    integer :: n
    integer :: i

    ! Calculate the integral of y(x) using trapezoidal quadrature

    n = SIZE(x)

    int_y(1) = 0._WP

    int_loop : do i = 2, n
       int_y(i) = int_y(i-1) + 0.5_WP*(y(i) + y(i-1))*(x(i) - x(i-1))
    end do int_loop

    ! Finish

    return

  end function integral_r_

  function integral_c_ (x, y) result (int_y)

    real(WP), intent(in)  :: x(:)
    complex(WP), intent(in) :: y(:)
    complex(WP)             :: int_y(SIZE(x))

    integer :: n
    integer :: i

    ! Calculate the integral of y(x) using trapezoidal quadrature

    n = SIZE(x)

    int_y(1) = 0._WP

    int_loop : do i = 2, n
       int_y(i) = int_y(i-1) + 0.5_WP*(y(i) + y(i-1))*(x(i) - x(i-1))
    end do int_loop

    ! Finish

    return

  end function integral_c_

end module gyre_util

!fpx3_header(0.13_3a)
!
!dependencies
!   dir: ~/Documents/MESA_SS/gyre/src/build 
!   sources: -
!   includes: ../extern/core/core.inc
!   uses: ISO_FORTRAN_ENV core_kinds
!   provides: core_system
!end dependencies
!
!end fpx3_header
! Module   : core_system
! Purpose  : operating system support

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

module core_system

  ! Uses

  use core_kinds

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Interfaces

  interface get_arg
     module procedure get_arg_i_i4_
     module procedure get_arg_i_i8_
     module procedure get_arg_r_sp_
     module procedure get_arg_r_dp_
     module procedure get_arg_c_sp_
     module procedure get_arg_c_dp_
     module procedure get_arg_a_
     module procedure get_arg_l_
  end interface get_arg

  interface get_env
     module procedure get_env_i_i4_
     module procedure get_env_i_i8_
     module procedure get_env_r_sp_
     module procedure get_env_r_dp_
     module procedure get_env_c_sp_
     module procedure get_env_c_dp_
     module procedure get_env_a_
     module procedure get_env_l_
  end interface get_env

  ! Access specifiers

  private

  public :: n_arg
  public :: get_arg
  public :: get_env

contains

  function n_arg ()

    integer :: n_arg

    ! Get the number of arguments

    n_arg = COMMAND_ARGUMENT_COUNT()

    ! Finish

    return

  end function n_arg

!****

  subroutine get_arg_i_i4_ (number, value, status)

    integer, intent(in)            :: number
    integer(I4), intent(out)       :: value
    integer, optional, intent(out) :: status

    integer                   :: length
    integer                   :: status_
    character(:), allocatable :: buffer

    ! Read the numbered command argument

    ! Determine the argument length

    call GET_COMMAND_ARGUMENT(number, LENGTH=length, STATUS=status_)

    if (PRESENT(status)) then
       status = status_
       if (status_ /= 0) return
    else

    if(.NOT. (status_ == 0)) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''status_ == 0'' failed at line 114 <core_system:get_arg_i_i4_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Error when reading command argument'
      stop
    endif

    endif

    ! Read the argument into a character buffer

    allocate(character(length)::buffer)

    call GET_COMMAND_ARGUMENT(number, buffer)

    ! Extract the value from the buffer

    read(buffer, *) value

    ! Finish

    return

  end subroutine get_arg_i_i4_

  subroutine get_arg_i_i8_ (number, value, status)

    integer, intent(in)            :: number
    integer(I8), intent(out)       :: value
    integer, optional, intent(out) :: status

    integer                   :: length
    integer                   :: status_
    character(:), allocatable :: buffer

    ! Read the numbered command argument

    ! Determine the argument length

    call GET_COMMAND_ARGUMENT(number, LENGTH=length, STATUS=status_)

    if (PRESENT(status)) then
       status = status_
       if (status_ /= 0) return
    else

    if(.NOT. (status_ == 0)) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''status_ == 0'' failed at line 115 <core_system:get_arg_i_i8_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Error when reading command argument'
      stop
    endif

    endif

    ! Read the argument into a character buffer

    allocate(character(length)::buffer)

    call GET_COMMAND_ARGUMENT(number, buffer)

    ! Extract the value from the buffer

    read(buffer, *) value

    ! Finish

    return

  end subroutine get_arg_i_i8_

  subroutine get_arg_r_sp_ (number, value, status)

    integer, intent(in)            :: number
    real(SP), intent(out)       :: value
    integer, optional, intent(out) :: status

    integer                   :: length
    integer                   :: status_
    character(:), allocatable :: buffer

    ! Read the numbered command argument

    ! Determine the argument length

    call GET_COMMAND_ARGUMENT(number, LENGTH=length, STATUS=status_)

    if (PRESENT(status)) then
       status = status_
       if (status_ /= 0) return
    else

    if(.NOT. (status_ == 0)) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''status_ == 0'' failed at line 116 <core_system:get_arg_r_sp_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Error when reading command argument'
      stop
    endif

    endif

    ! Read the argument into a character buffer

    allocate(character(length)::buffer)

    call GET_COMMAND_ARGUMENT(number, buffer)

    ! Extract the value from the buffer

    read(buffer, *) value

    ! Finish

    return

  end subroutine get_arg_r_sp_

  subroutine get_arg_r_dp_ (number, value, status)

    integer, intent(in)            :: number
    real(DP), intent(out)       :: value
    integer, optional, intent(out) :: status

    integer                   :: length
    integer                   :: status_
    character(:), allocatable :: buffer

    ! Read the numbered command argument

    ! Determine the argument length

    call GET_COMMAND_ARGUMENT(number, LENGTH=length, STATUS=status_)

    if (PRESENT(status)) then
       status = status_
       if (status_ /= 0) return
    else

    if(.NOT. (status_ == 0)) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''status_ == 0'' failed at line 117 <core_system:get_arg_r_dp_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Error when reading command argument'
      stop
    endif

    endif

    ! Read the argument into a character buffer

    allocate(character(length)::buffer)

    call GET_COMMAND_ARGUMENT(number, buffer)

    ! Extract the value from the buffer

    read(buffer, *) value

    ! Finish

    return

  end subroutine get_arg_r_dp_

  subroutine get_arg_c_sp_ (number, value, status)

    integer, intent(in)            :: number
    complex(SP), intent(out)       :: value
    integer, optional, intent(out) :: status

    integer                   :: length
    integer                   :: status_
    character(:), allocatable :: buffer

    ! Read the numbered command argument

    ! Determine the argument length

    call GET_COMMAND_ARGUMENT(number, LENGTH=length, STATUS=status_)

    if (PRESENT(status)) then
       status = status_
       if (status_ /= 0) return
    else

    if(.NOT. (status_ == 0)) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''status_ == 0'' failed at line 118 <core_system:get_arg_c_sp_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Error when reading command argument'
      stop
    endif

    endif

    ! Read the argument into a character buffer

    allocate(character(length)::buffer)

    call GET_COMMAND_ARGUMENT(number, buffer)

    ! Extract the value from the buffer

    read(buffer, *) value

    ! Finish

    return

  end subroutine get_arg_c_sp_

  subroutine get_arg_c_dp_ (number, value, status)

    integer, intent(in)            :: number
    complex(DP), intent(out)       :: value
    integer, optional, intent(out) :: status

    integer                   :: length
    integer                   :: status_
    character(:), allocatable :: buffer

    ! Read the numbered command argument

    ! Determine the argument length

    call GET_COMMAND_ARGUMENT(number, LENGTH=length, STATUS=status_)

    if (PRESENT(status)) then
       status = status_
       if (status_ /= 0) return
    else

    if(.NOT. (status_ == 0)) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''status_ == 0'' failed at line 119 <core_system:get_arg_c_dp_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Error when reading command argument'
      stop
    endif

    endif

    ! Read the argument into a character buffer

    allocate(character(length)::buffer)

    call GET_COMMAND_ARGUMENT(number, buffer)

    ! Extract the value from the buffer

    read(buffer, *) value

    ! Finish

    return

  end subroutine get_arg_c_dp_

  subroutine get_arg_l_ (number, value, status)

    integer, intent(in)            :: number
    logical, intent(out)       :: value
    integer, optional, intent(out) :: status

    integer                   :: length
    integer                   :: status_
    character(:), allocatable :: buffer

    ! Read the numbered command argument

    ! Determine the argument length

    call GET_COMMAND_ARGUMENT(number, LENGTH=length, STATUS=status_)

    if (PRESENT(status)) then
       status = status_
       if (status_ /= 0) return
    else

    if(.NOT. (status_ == 0)) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''status_ == 0'' failed at line 120 <core_system:get_arg_l_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Error when reading command argument'
      stop
    endif

    endif

    ! Read the argument into a character buffer

    allocate(character(length)::buffer)

    call GET_COMMAND_ARGUMENT(number, buffer)

    ! Extract the value from the buffer

    read(buffer, *) value

    ! Finish

    return

  end subroutine get_arg_l_

!****

  subroutine get_arg_a_ (number, value, status)

    integer, intent(in)                    :: number
    character(:), allocatable, intent(out) :: value
    integer, optional, intent(out)         :: status

    integer :: length
    integer :: status_

    ! Read the numbered command argument

    ! Determine the argument length

    call GET_COMMAND_ARGUMENT(number, LENGTH=length, STATUS=status_)

    if (PRESENT(status)) then
       status = status_
       if (status_ /= 0) return
    else

    if(.NOT. (status_ == 0)) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''status_ == 0'' failed at line 143 <core_system:get_arg_a_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Error when reading command argument'
      stop
    endif

    endif

    ! Read the argument

    if (length > 0) then
       allocate(character(length)::value)
       call GET_COMMAND_ARGUMENT(number, value)
    endif

    ! Finish

    return

  end subroutine get_arg_a_

!****

  subroutine get_env_i_i4_ (name, value, status)

    character(*), intent(in)       :: name
    integer(I4), intent(out)       :: value
    integer, optional, intent(out) :: status

    integer                   :: length
    integer                   :: status_
    character(:), allocatable :: buffer

    ! Read the named environment variable

    ! Determine the variable length

    call GET_ENVIRONMENT_VARIABLE(name, LENGTH=length, STATUS=status_)

    if (PRESENT(status)) then
       status = status_
       if (status_ /= 0) return
    else

    if(.NOT. (status_ == 0)) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''status_ == 0'' failed at line 207 <core_system:get_env_i_i4_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Error when reading environment variable'
      stop
    endif

    endif

    ! Read the variable into a character buffer

    allocate(character(length)::buffer)

    call GET_ENVIRONMENT_VARIABLE(name, buffer)

    ! Extract the value from the buffer

    read(buffer, *) value

    ! Finish

    return

  end subroutine get_env_i_i4_

  subroutine get_env_i_i8_ (name, value, status)

    character(*), intent(in)       :: name
    integer(I8), intent(out)       :: value
    integer, optional, intent(out) :: status

    integer                   :: length
    integer                   :: status_
    character(:), allocatable :: buffer

    ! Read the named environment variable

    ! Determine the variable length

    call GET_ENVIRONMENT_VARIABLE(name, LENGTH=length, STATUS=status_)

    if (PRESENT(status)) then
       status = status_
       if (status_ /= 0) return
    else

    if(.NOT. (status_ == 0)) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''status_ == 0'' failed at line 208 <core_system:get_env_i_i8_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Error when reading environment variable'
      stop
    endif

    endif

    ! Read the variable into a character buffer

    allocate(character(length)::buffer)

    call GET_ENVIRONMENT_VARIABLE(name, buffer)

    ! Extract the value from the buffer

    read(buffer, *) value

    ! Finish

    return

  end subroutine get_env_i_i8_

  subroutine get_env_r_sp_ (name, value, status)

    character(*), intent(in)       :: name
    real(SP), intent(out)       :: value
    integer, optional, intent(out) :: status

    integer                   :: length
    integer                   :: status_
    character(:), allocatable :: buffer

    ! Read the named environment variable

    ! Determine the variable length

    call GET_ENVIRONMENT_VARIABLE(name, LENGTH=length, STATUS=status_)

    if (PRESENT(status)) then
       status = status_
       if (status_ /= 0) return
    else

    if(.NOT. (status_ == 0)) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''status_ == 0'' failed at line 209 <core_system:get_env_r_sp_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Error when reading environment variable'
      stop
    endif

    endif

    ! Read the variable into a character buffer

    allocate(character(length)::buffer)

    call GET_ENVIRONMENT_VARIABLE(name, buffer)

    ! Extract the value from the buffer

    read(buffer, *) value

    ! Finish

    return

  end subroutine get_env_r_sp_

  subroutine get_env_r_dp_ (name, value, status)

    character(*), intent(in)       :: name
    real(DP), intent(out)       :: value
    integer, optional, intent(out) :: status

    integer                   :: length
    integer                   :: status_
    character(:), allocatable :: buffer

    ! Read the named environment variable

    ! Determine the variable length

    call GET_ENVIRONMENT_VARIABLE(name, LENGTH=length, STATUS=status_)

    if (PRESENT(status)) then
       status = status_
       if (status_ /= 0) return
    else

    if(.NOT. (status_ == 0)) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''status_ == 0'' failed at line 210 <core_system:get_env_r_dp_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Error when reading environment variable'
      stop
    endif

    endif

    ! Read the variable into a character buffer

    allocate(character(length)::buffer)

    call GET_ENVIRONMENT_VARIABLE(name, buffer)

    ! Extract the value from the buffer

    read(buffer, *) value

    ! Finish

    return

  end subroutine get_env_r_dp_

  subroutine get_env_c_sp_ (name, value, status)

    character(*), intent(in)       :: name
    complex(SP), intent(out)       :: value
    integer, optional, intent(out) :: status

    integer                   :: length
    integer                   :: status_
    character(:), allocatable :: buffer

    ! Read the named environment variable

    ! Determine the variable length

    call GET_ENVIRONMENT_VARIABLE(name, LENGTH=length, STATUS=status_)

    if (PRESENT(status)) then
       status = status_
       if (status_ /= 0) return
    else

    if(.NOT. (status_ == 0)) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''status_ == 0'' failed at line 211 <core_system:get_env_c_sp_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Error when reading environment variable'
      stop
    endif

    endif

    ! Read the variable into a character buffer

    allocate(character(length)::buffer)

    call GET_ENVIRONMENT_VARIABLE(name, buffer)

    ! Extract the value from the buffer

    read(buffer, *) value

    ! Finish

    return

  end subroutine get_env_c_sp_

  subroutine get_env_c_dp_ (name, value, status)

    character(*), intent(in)       :: name
    complex(DP), intent(out)       :: value
    integer, optional, intent(out) :: status

    integer                   :: length
    integer                   :: status_
    character(:), allocatable :: buffer

    ! Read the named environment variable

    ! Determine the variable length

    call GET_ENVIRONMENT_VARIABLE(name, LENGTH=length, STATUS=status_)

    if (PRESENT(status)) then
       status = status_
       if (status_ /= 0) return
    else

    if(.NOT. (status_ == 0)) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''status_ == 0'' failed at line 212 <core_system:get_env_c_dp_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Error when reading environment variable'
      stop
    endif

    endif

    ! Read the variable into a character buffer

    allocate(character(length)::buffer)

    call GET_ENVIRONMENT_VARIABLE(name, buffer)

    ! Extract the value from the buffer

    read(buffer, *) value

    ! Finish

    return

  end subroutine get_env_c_dp_

  subroutine get_env_l_ (name, value, status)

    character(*), intent(in)       :: name
    logical, intent(out)       :: value
    integer, optional, intent(out) :: status

    integer                   :: length
    integer                   :: status_
    character(:), allocatable :: buffer

    ! Read the named environment variable

    ! Determine the variable length

    call GET_ENVIRONMENT_VARIABLE(name, LENGTH=length, STATUS=status_)

    if (PRESENT(status)) then
       status = status_
       if (status_ /= 0) return
    else

    if(.NOT. (status_ == 0)) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''status_ == 0'' failed at line 213 <core_system:get_env_l_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Error when reading environment variable'
      stop
    endif

    endif

    ! Read the variable into a character buffer

    allocate(character(length)::buffer)

    call GET_ENVIRONMENT_VARIABLE(name, buffer)

    ! Extract the value from the buffer

    read(buffer, *) value

    ! Finish

    return

  end subroutine get_env_l_

!****

  subroutine get_env_a_ (name, value, status)

    character(*), intent(in)               :: name
    character(:), allocatable, intent(out) :: value
    integer, optional, intent(out)         :: status

    integer :: length
    integer :: status_

    ! Read the named environment variable

    ! Determine the variable length

    call GET_ENVIRONMENT_VARIABLE(name, LENGTH=length, STATUS=status_)

    if (PRESENT(status)) then
       status = status_
       if (status_ /= 0) return
    else

    if(.NOT. (status_ == 0)) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''status_ == 0'' failed at line 236 <core_system:get_env_a_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Error when reading environment variable'
      stop
    endif

    endif

    ! Read the variable

    if (length > 0) then
       allocate(character(length)::value)
       call GET_ENVIRONMENT_VARIABLE(name, value)
    endif

    ! Finish

    return

  end subroutine get_env_a_

end module core_system

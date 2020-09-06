!fpx3_header(0.13_3a)
!
!dependencies
!   dir: ~/Documents/MESA_SS/gyre/src/build 
!   sources: -
!   includes: ../extern/core/core.inc
!   uses: ISO_FORTRAN_ENV core_order omp_lib core_kinds
!   provides: core_parallel
!end dependencies
!
!end fpx3_header
! Module   : core_parallel
! Purpose  : parallel support

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

module core_parallel

  ! Uses

  use core_kinds
  use core_order

  use ISO_FORTRAN_ENV

  use omp_lib

  ! No implicit typing

  implicit none

  ! Module variables

  integer, save, protected :: MPI_SIZE
  integer, save, protected :: MPI_RANK
  integer, save, protected :: OMP_SIZE_MAX

  ! Interfaces

  ! Access specifiers

  private

  public :: MPI_SIZE
  public :: MPI_RANK
  public :: OMP_SIZE_MAX
  public :: init_parallel
  public :: final_parallel
  public :: omp_size
  public :: omp_rank

  public :: omp_get_thread_num

  public :: partition_tasks

contains

  subroutine init_parallel ()

    ! Initialize MPI

    MPI_SIZE = 1
    MPI_RANK = 0

    ! Initialize OpenMP

    OMP_SIZE_MAX = omp_get_max_threads()

    ! Finish

    return

  end subroutine init_parallel

!****

  subroutine final_parallel ()

    ! Finalize MPI

    MPI_SIZE = 0
    MPI_RANK = 0

    ! Finish

    return

  end subroutine final_parallel

!****

  function omp_size ()

    integer :: omp_size

    ! Get the OpenMP thread size

    omp_size = omp_get_num_threads()

    ! Finish

    return

  end function omp_size

!****

  function omp_rank ()

    integer :: omp_rank

    ! Get the OpenMP thread rank

    omp_rank = omp_get_thread_num()

    ! finish

    return

  end function omp_rank

!****

!****

  subroutine partition_tasks (n, m, k_part)

    integer, intent(in)  :: n
    integer, intent(in)  :: m
    integer, intent(out) :: k_part(:)

    integer :: n_work
    integer :: n_part(SIZE(k_part)-1)
    integer :: i

    n_work = SIZE(k_part) - 1

    ! Partition n tasks among n_work workers --- where possible, at
    ! least m tasks per worker. The resulting partitioning indices are
    ! returned in the array k_part

    n_part = m*(n/(m*n_work))

    size_loop : do i = 1, n_work
       if(SUM(n_part) >= n-m+1) then
          n_part(i) = n_part(i) + n - SUM(n_part)
          exit size_loop
       else
          n_part(i) = n_part(i) + m
       endif
    end do size_loop

    if(.NOT. (SUM(n_part) == n)) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''SUM(n_part) == n'' failed at line 1324 <core_parallel:partition_tasks>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Partitioning failed'
      stop
    endif

    k_part(1) = 1

    index_loop : do i = 1, n_work
       k_part(i+1) = k_part(i) + n_part(i)
    end do index_loop

    ! Finish

    return

  end subroutine partition_tasks

end module core_parallel

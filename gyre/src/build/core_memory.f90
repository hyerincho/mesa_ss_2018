!fpx3_header(0.13_3a)
!
!dependencies
!   dir: ~/Documents/MESA_SS/gyre/src/build 
!   sources: -
!   includes: ../extern/core/core_memory.inc ../extern/core/core.inc
!   uses: core_kinds ISO_FORTRAN_ENV
!   provides: core_memory
!end dependencies
!
!end fpx3_header
! Module   : core_memory
! Purpose  : memory management

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

! Incfile  : core_memory
! Purpose  : memory management fpx3 macros

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

module core_memory

  ! Uses

  use core_kinds

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Interfaces

  interface reallocate
     module procedure reallocate_i_i4_1_
     module procedure reallocate_i_i4_2_
     module procedure reallocate_i_i4_3_
     module procedure reallocate_i_i4_4_
     module procedure reallocate_i_i8_1_
     module procedure reallocate_i_i8_2_
     module procedure reallocate_i_i8_3_
     module procedure reallocate_i_i8_4_
     module procedure reallocate_r_sp_1_
     module procedure reallocate_r_sp_2_
     module procedure reallocate_r_sp_3_
     module procedure reallocate_r_sp_4_
     module procedure reallocate_r_dp_1_
     module procedure reallocate_r_dp_2_
     module procedure reallocate_r_dp_3_
     module procedure reallocate_r_dp_4_
     module procedure reallocate_c_sp_1_
     module procedure reallocate_c_sp_2_
     module procedure reallocate_c_sp_3_
     module procedure reallocate_c_sp_4_
     module procedure reallocate_c_dp_1_
     module procedure reallocate_c_dp_2_
     module procedure reallocate_c_dp_3_
     module procedure reallocate_c_dp_4_
     module procedure reallocate_a_1_
     module procedure reallocate_a_2_
     module procedure reallocate_a_3_
     module procedure reallocate_a_4_
     module procedure reallocate_l_1_
     module procedure reallocate_l_2_
     module procedure reallocate_l_3_
     module procedure reallocate_l_4_
  end interface reallocate

  ! Access specifiers

  private

  public :: reallocate

  ! Procedures

contains

subroutine reallocate_i_i4_1_ (array, shape_new, start)

  integer(I4), allocatable, intent(inout) :: array(:)
  integer, intent(in)                     :: shape_new(:)
  integer, intent(in), optional           :: start(:)

  integer(I4), allocatable :: array_new(:)
  integer                      :: i_a(1)
  integer                      :: i_b(1)
  integer                      :: i_c(1)

  if(PRESENT(start)) then

  end if

  ! Set up the reallocation indices

  if(PRESENT(start)) then
     i_a = start - 1
  else
     i_a = MIN(SHAPE(array), shape_new)
  endif

  i_b = MAX(i_a, i_a+SHAPE(array)-shape_new) + 1
  i_c = MAX(i_a, i_a+shape_new-SHAPE(array)) + 1

  ! Allocate the new array

  allocate(array_new(shape_new(1)))

  ! Copy the data from old to new

  array_new(:i_a(1)) = array(:i_a(1))
array_new(i_c(1):) = array(i_b(1):)

  ! Move the allocation

  call MOVE_ALLOC(array_new, array)

  ! Finish

  return

end subroutine reallocate_i_i4_1_

subroutine reallocate_i_i4_2_ (array, shape_new, start)

  integer(I4), allocatable, intent(inout) :: array(:,:)
  integer, intent(in)                     :: shape_new(:)
  integer, intent(in), optional           :: start(:)

  integer(I4), allocatable :: array_new(:,:)
  integer                      :: i_a(2)
  integer                      :: i_b(2)
  integer                      :: i_c(2)

  if(PRESENT(start)) then

  end if

  ! Set up the reallocation indices

  if(PRESENT(start)) then
     i_a = start - 1
  else
     i_a = MIN(SHAPE(array), shape_new)
  endif

  i_b = MAX(i_a, i_a+SHAPE(array)-shape_new) + 1
  i_c = MAX(i_a, i_a+shape_new-SHAPE(array)) + 1

  ! Allocate the new array

  allocate(array_new(shape_new(1),shape_new(2)))

  ! Copy the data from old to new

  array_new(:i_a(1),:i_a(2)) = array(:i_a(1),:i_a(2))
array_new(i_c(1):,:i_a(2)) = array(i_b(1):,:i_a(2))
array_new(:i_a(1),i_c(2):) = array(:i_a(1),i_b(2):)
array_new(i_c(1):,i_c(2):) = array(i_b(1):,i_b(2):)

  ! Move the allocation

  call MOVE_ALLOC(array_new, array)

  ! Finish

  return

end subroutine reallocate_i_i4_2_

subroutine reallocate_i_i4_3_ (array, shape_new, start)

  integer(I4), allocatable, intent(inout) :: array(:,:,:)
  integer, intent(in)                     :: shape_new(:)
  integer, intent(in), optional           :: start(:)

  integer(I4), allocatable :: array_new(:,:,:)
  integer                      :: i_a(3)
  integer                      :: i_b(3)
  integer                      :: i_c(3)

  if(PRESENT(start)) then

  end if

  ! Set up the reallocation indices

  if(PRESENT(start)) then
     i_a = start - 1
  else
     i_a = MIN(SHAPE(array), shape_new)
  endif

  i_b = MAX(i_a, i_a+SHAPE(array)-shape_new) + 1
  i_c = MAX(i_a, i_a+shape_new-SHAPE(array)) + 1

  ! Allocate the new array

  allocate(array_new(shape_new(1),shape_new(2),shape_new(3)))

  ! Copy the data from old to new

  array_new(:i_a(1),:i_a(2),:i_a(3)) = array(:i_a(1),:i_a(2),:i_a(3))
array_new(i_c(1):,:i_a(2),:i_a(3)) = array(i_b(1):,:i_a(2),:i_a(3))
array_new(:i_a(1),i_c(2):,:i_a(3)) = array(:i_a(1),i_b(2):,:i_a(3))
array_new(i_c(1):,i_c(2):,:i_a(3)) = array(i_b(1):,i_b(2):,:i_a(3))
array_new(:i_a(1),:i_a(2),i_c(3):) = array(:i_a(1),:i_a(2),i_b(3):)
array_new(i_c(1):,:i_a(2),i_c(3):) = array(i_b(1):,:i_a(2),i_b(3):)
array_new(:i_a(1),i_c(2):,i_c(3):) = array(:i_a(1),i_b(2):,i_b(3):)
array_new(i_c(1):,i_c(2):,i_c(3):) = array(i_b(1):,i_b(2):,i_b(3):)

  ! Move the allocation

  call MOVE_ALLOC(array_new, array)

  ! Finish

  return

end subroutine reallocate_i_i4_3_

subroutine reallocate_i_i4_4_ (array, shape_new, start)

  integer(I4), allocatable, intent(inout) :: array(:,:,:,:)
  integer, intent(in)                     :: shape_new(:)
  integer, intent(in), optional           :: start(:)

  integer(I4), allocatable :: array_new(:,:,:,:)
  integer                      :: i_a(4)
  integer                      :: i_b(4)
  integer                      :: i_c(4)

  if(PRESENT(start)) then

  end if

  ! Set up the reallocation indices

  if(PRESENT(start)) then
     i_a = start - 1
  else
     i_a = MIN(SHAPE(array), shape_new)
  endif

  i_b = MAX(i_a, i_a+SHAPE(array)-shape_new) + 1
  i_c = MAX(i_a, i_a+shape_new-SHAPE(array)) + 1

  ! Allocate the new array

  allocate(array_new(shape_new(1),shape_new(2),shape_new(3),shape_new(4)))

  ! Copy the data from old to new

  array_new(:i_a(1),:i_a(2),:i_a(3),:i_a(4)) = array(:i_a(1),:i_a(2),:i_a(3),:i_a(4))
array_new(i_c(1):,:i_a(2),:i_a(3),:i_a(4)) = array(i_b(1):,:i_a(2),:i_a(3),:i_a(4))
array_new(:i_a(1),i_c(2):,:i_a(3),:i_a(4)) = array(:i_a(1),i_b(2):,:i_a(3),:i_a(4))
array_new(i_c(1):,i_c(2):,:i_a(3),:i_a(4)) = array(i_b(1):,i_b(2):,:i_a(3),:i_a(4))
array_new(:i_a(1),:i_a(2),i_c(3):,:i_a(4)) = array(:i_a(1),:i_a(2),i_b(3):,:i_a(4))
array_new(i_c(1):,:i_a(2),i_c(3):,:i_a(4)) = array(i_b(1):,:i_a(2),i_b(3):,:i_a(4))
array_new(:i_a(1),i_c(2):,i_c(3):,:i_a(4)) = array(:i_a(1),i_b(2):,i_b(3):,:i_a(4))
array_new(i_c(1):,i_c(2):,i_c(3):,:i_a(4)) = array(i_b(1):,i_b(2):,i_b(3):,:i_a(4))
array_new(:i_a(1),:i_a(2),:i_a(3),i_c(4):) = array(:i_a(1),:i_a(2),:i_a(3),i_b(4):)
array_new(i_c(1):,:i_a(2),:i_a(3),i_c(4):) = array(i_b(1):,:i_a(2),:i_a(3),i_b(4):)
array_new(:i_a(1),i_c(2):,:i_a(3),i_c(4):) = array(:i_a(1),i_b(2):,:i_a(3),i_b(4):)
array_new(i_c(1):,i_c(2):,:i_a(3),i_c(4):) = array(i_b(1):,i_b(2):,:i_a(3),i_b(4):)
array_new(:i_a(1),:i_a(2),i_c(3):,i_c(4):) = array(:i_a(1),:i_a(2),i_b(3):,i_b(4):)
array_new(i_c(1):,:i_a(2),i_c(3):,i_c(4):) = array(i_b(1):,:i_a(2),i_b(3):,i_b(4):)
array_new(:i_a(1),i_c(2):,i_c(3):,i_c(4):) = array(:i_a(1),i_b(2):,i_b(3):,i_b(4):)
array_new(i_c(1):,i_c(2):,i_c(3):,i_c(4):) = array(i_b(1):,i_b(2):,i_b(3):,i_b(4):)

  ! Move the allocation

  call MOVE_ALLOC(array_new, array)

  ! Finish

  return

end subroutine reallocate_i_i4_4_

subroutine reallocate_i_i8_1_ (array, shape_new, start)

  integer(I8), allocatable, intent(inout) :: array(:)
  integer, intent(in)                     :: shape_new(:)
  integer, intent(in), optional           :: start(:)

  integer(I8), allocatable :: array_new(:)
  integer                      :: i_a(1)
  integer                      :: i_b(1)
  integer                      :: i_c(1)

  if(PRESENT(start)) then

  end if

  ! Set up the reallocation indices

  if(PRESENT(start)) then
     i_a = start - 1
  else
     i_a = MIN(SHAPE(array), shape_new)
  endif

  i_b = MAX(i_a, i_a+SHAPE(array)-shape_new) + 1
  i_c = MAX(i_a, i_a+shape_new-SHAPE(array)) + 1

  ! Allocate the new array

  allocate(array_new(shape_new(1)))

  ! Copy the data from old to new

  array_new(:i_a(1)) = array(:i_a(1))
array_new(i_c(1):) = array(i_b(1):)

  ! Move the allocation

  call MOVE_ALLOC(array_new, array)

  ! Finish

  return

end subroutine reallocate_i_i8_1_

subroutine reallocate_i_i8_2_ (array, shape_new, start)

  integer(I8), allocatable, intent(inout) :: array(:,:)
  integer, intent(in)                     :: shape_new(:)
  integer, intent(in), optional           :: start(:)

  integer(I8), allocatable :: array_new(:,:)
  integer                      :: i_a(2)
  integer                      :: i_b(2)
  integer                      :: i_c(2)

  if(PRESENT(start)) then

  end if

  ! Set up the reallocation indices

  if(PRESENT(start)) then
     i_a = start - 1
  else
     i_a = MIN(SHAPE(array), shape_new)
  endif

  i_b = MAX(i_a, i_a+SHAPE(array)-shape_new) + 1
  i_c = MAX(i_a, i_a+shape_new-SHAPE(array)) + 1

  ! Allocate the new array

  allocate(array_new(shape_new(1),shape_new(2)))

  ! Copy the data from old to new

  array_new(:i_a(1),:i_a(2)) = array(:i_a(1),:i_a(2))
array_new(i_c(1):,:i_a(2)) = array(i_b(1):,:i_a(2))
array_new(:i_a(1),i_c(2):) = array(:i_a(1),i_b(2):)
array_new(i_c(1):,i_c(2):) = array(i_b(1):,i_b(2):)

  ! Move the allocation

  call MOVE_ALLOC(array_new, array)

  ! Finish

  return

end subroutine reallocate_i_i8_2_

subroutine reallocate_i_i8_3_ (array, shape_new, start)

  integer(I8), allocatable, intent(inout) :: array(:,:,:)
  integer, intent(in)                     :: shape_new(:)
  integer, intent(in), optional           :: start(:)

  integer(I8), allocatable :: array_new(:,:,:)
  integer                      :: i_a(3)
  integer                      :: i_b(3)
  integer                      :: i_c(3)

  if(PRESENT(start)) then

  end if

  ! Set up the reallocation indices

  if(PRESENT(start)) then
     i_a = start - 1
  else
     i_a = MIN(SHAPE(array), shape_new)
  endif

  i_b = MAX(i_a, i_a+SHAPE(array)-shape_new) + 1
  i_c = MAX(i_a, i_a+shape_new-SHAPE(array)) + 1

  ! Allocate the new array

  allocate(array_new(shape_new(1),shape_new(2),shape_new(3)))

  ! Copy the data from old to new

  array_new(:i_a(1),:i_a(2),:i_a(3)) = array(:i_a(1),:i_a(2),:i_a(3))
array_new(i_c(1):,:i_a(2),:i_a(3)) = array(i_b(1):,:i_a(2),:i_a(3))
array_new(:i_a(1),i_c(2):,:i_a(3)) = array(:i_a(1),i_b(2):,:i_a(3))
array_new(i_c(1):,i_c(2):,:i_a(3)) = array(i_b(1):,i_b(2):,:i_a(3))
array_new(:i_a(1),:i_a(2),i_c(3):) = array(:i_a(1),:i_a(2),i_b(3):)
array_new(i_c(1):,:i_a(2),i_c(3):) = array(i_b(1):,:i_a(2),i_b(3):)
array_new(:i_a(1),i_c(2):,i_c(3):) = array(:i_a(1),i_b(2):,i_b(3):)
array_new(i_c(1):,i_c(2):,i_c(3):) = array(i_b(1):,i_b(2):,i_b(3):)

  ! Move the allocation

  call MOVE_ALLOC(array_new, array)

  ! Finish

  return

end subroutine reallocate_i_i8_3_

subroutine reallocate_i_i8_4_ (array, shape_new, start)

  integer(I8), allocatable, intent(inout) :: array(:,:,:,:)
  integer, intent(in)                     :: shape_new(:)
  integer, intent(in), optional           :: start(:)

  integer(I8), allocatable :: array_new(:,:,:,:)
  integer                      :: i_a(4)
  integer                      :: i_b(4)
  integer                      :: i_c(4)

  if(PRESENT(start)) then

  end if

  ! Set up the reallocation indices

  if(PRESENT(start)) then
     i_a = start - 1
  else
     i_a = MIN(SHAPE(array), shape_new)
  endif

  i_b = MAX(i_a, i_a+SHAPE(array)-shape_new) + 1
  i_c = MAX(i_a, i_a+shape_new-SHAPE(array)) + 1

  ! Allocate the new array

  allocate(array_new(shape_new(1),shape_new(2),shape_new(3),shape_new(4)))

  ! Copy the data from old to new

  array_new(:i_a(1),:i_a(2),:i_a(3),:i_a(4)) = array(:i_a(1),:i_a(2),:i_a(3),:i_a(4))
array_new(i_c(1):,:i_a(2),:i_a(3),:i_a(4)) = array(i_b(1):,:i_a(2),:i_a(3),:i_a(4))
array_new(:i_a(1),i_c(2):,:i_a(3),:i_a(4)) = array(:i_a(1),i_b(2):,:i_a(3),:i_a(4))
array_new(i_c(1):,i_c(2):,:i_a(3),:i_a(4)) = array(i_b(1):,i_b(2):,:i_a(3),:i_a(4))
array_new(:i_a(1),:i_a(2),i_c(3):,:i_a(4)) = array(:i_a(1),:i_a(2),i_b(3):,:i_a(4))
array_new(i_c(1):,:i_a(2),i_c(3):,:i_a(4)) = array(i_b(1):,:i_a(2),i_b(3):,:i_a(4))
array_new(:i_a(1),i_c(2):,i_c(3):,:i_a(4)) = array(:i_a(1),i_b(2):,i_b(3):,:i_a(4))
array_new(i_c(1):,i_c(2):,i_c(3):,:i_a(4)) = array(i_b(1):,i_b(2):,i_b(3):,:i_a(4))
array_new(:i_a(1),:i_a(2),:i_a(3),i_c(4):) = array(:i_a(1),:i_a(2),:i_a(3),i_b(4):)
array_new(i_c(1):,:i_a(2),:i_a(3),i_c(4):) = array(i_b(1):,:i_a(2),:i_a(3),i_b(4):)
array_new(:i_a(1),i_c(2):,:i_a(3),i_c(4):) = array(:i_a(1),i_b(2):,:i_a(3),i_b(4):)
array_new(i_c(1):,i_c(2):,:i_a(3),i_c(4):) = array(i_b(1):,i_b(2):,:i_a(3),i_b(4):)
array_new(:i_a(1),:i_a(2),i_c(3):,i_c(4):) = array(:i_a(1),:i_a(2),i_b(3):,i_b(4):)
array_new(i_c(1):,:i_a(2),i_c(3):,i_c(4):) = array(i_b(1):,:i_a(2),i_b(3):,i_b(4):)
array_new(:i_a(1),i_c(2):,i_c(3):,i_c(4):) = array(:i_a(1),i_b(2):,i_b(3):,i_b(4):)
array_new(i_c(1):,i_c(2):,i_c(3):,i_c(4):) = array(i_b(1):,i_b(2):,i_b(3):,i_b(4):)

  ! Move the allocation

  call MOVE_ALLOC(array_new, array)

  ! Finish

  return

end subroutine reallocate_i_i8_4_

subroutine reallocate_r_sp_1_ (array, shape_new, start)

  real(SP), allocatable, intent(inout) :: array(:)
  integer, intent(in)                     :: shape_new(:)
  integer, intent(in), optional           :: start(:)

  real(SP), allocatable :: array_new(:)
  integer                      :: i_a(1)
  integer                      :: i_b(1)
  integer                      :: i_c(1)

  if(PRESENT(start)) then

  end if

  ! Set up the reallocation indices

  if(PRESENT(start)) then
     i_a = start - 1
  else
     i_a = MIN(SHAPE(array), shape_new)
  endif

  i_b = MAX(i_a, i_a+SHAPE(array)-shape_new) + 1
  i_c = MAX(i_a, i_a+shape_new-SHAPE(array)) + 1

  ! Allocate the new array

  allocate(array_new(shape_new(1)))

  ! Copy the data from old to new

  array_new(:i_a(1)) = array(:i_a(1))
array_new(i_c(1):) = array(i_b(1):)

  ! Move the allocation

  call MOVE_ALLOC(array_new, array)

  ! Finish

  return

end subroutine reallocate_r_sp_1_

subroutine reallocate_r_sp_2_ (array, shape_new, start)

  real(SP), allocatable, intent(inout) :: array(:,:)
  integer, intent(in)                     :: shape_new(:)
  integer, intent(in), optional           :: start(:)

  real(SP), allocatable :: array_new(:,:)
  integer                      :: i_a(2)
  integer                      :: i_b(2)
  integer                      :: i_c(2)

  if(PRESENT(start)) then

  end if

  ! Set up the reallocation indices

  if(PRESENT(start)) then
     i_a = start - 1
  else
     i_a = MIN(SHAPE(array), shape_new)
  endif

  i_b = MAX(i_a, i_a+SHAPE(array)-shape_new) + 1
  i_c = MAX(i_a, i_a+shape_new-SHAPE(array)) + 1

  ! Allocate the new array

  allocate(array_new(shape_new(1),shape_new(2)))

  ! Copy the data from old to new

  array_new(:i_a(1),:i_a(2)) = array(:i_a(1),:i_a(2))
array_new(i_c(1):,:i_a(2)) = array(i_b(1):,:i_a(2))
array_new(:i_a(1),i_c(2):) = array(:i_a(1),i_b(2):)
array_new(i_c(1):,i_c(2):) = array(i_b(1):,i_b(2):)

  ! Move the allocation

  call MOVE_ALLOC(array_new, array)

  ! Finish

  return

end subroutine reallocate_r_sp_2_

subroutine reallocate_r_sp_3_ (array, shape_new, start)

  real(SP), allocatable, intent(inout) :: array(:,:,:)
  integer, intent(in)                     :: shape_new(:)
  integer, intent(in), optional           :: start(:)

  real(SP), allocatable :: array_new(:,:,:)
  integer                      :: i_a(3)
  integer                      :: i_b(3)
  integer                      :: i_c(3)

  if(PRESENT(start)) then

  end if

  ! Set up the reallocation indices

  if(PRESENT(start)) then
     i_a = start - 1
  else
     i_a = MIN(SHAPE(array), shape_new)
  endif

  i_b = MAX(i_a, i_a+SHAPE(array)-shape_new) + 1
  i_c = MAX(i_a, i_a+shape_new-SHAPE(array)) + 1

  ! Allocate the new array

  allocate(array_new(shape_new(1),shape_new(2),shape_new(3)))

  ! Copy the data from old to new

  array_new(:i_a(1),:i_a(2),:i_a(3)) = array(:i_a(1),:i_a(2),:i_a(3))
array_new(i_c(1):,:i_a(2),:i_a(3)) = array(i_b(1):,:i_a(2),:i_a(3))
array_new(:i_a(1),i_c(2):,:i_a(3)) = array(:i_a(1),i_b(2):,:i_a(3))
array_new(i_c(1):,i_c(2):,:i_a(3)) = array(i_b(1):,i_b(2):,:i_a(3))
array_new(:i_a(1),:i_a(2),i_c(3):) = array(:i_a(1),:i_a(2),i_b(3):)
array_new(i_c(1):,:i_a(2),i_c(3):) = array(i_b(1):,:i_a(2),i_b(3):)
array_new(:i_a(1),i_c(2):,i_c(3):) = array(:i_a(1),i_b(2):,i_b(3):)
array_new(i_c(1):,i_c(2):,i_c(3):) = array(i_b(1):,i_b(2):,i_b(3):)

  ! Move the allocation

  call MOVE_ALLOC(array_new, array)

  ! Finish

  return

end subroutine reallocate_r_sp_3_

subroutine reallocate_r_sp_4_ (array, shape_new, start)

  real(SP), allocatable, intent(inout) :: array(:,:,:,:)
  integer, intent(in)                     :: shape_new(:)
  integer, intent(in), optional           :: start(:)

  real(SP), allocatable :: array_new(:,:,:,:)
  integer                      :: i_a(4)
  integer                      :: i_b(4)
  integer                      :: i_c(4)

  if(PRESENT(start)) then

  end if

  ! Set up the reallocation indices

  if(PRESENT(start)) then
     i_a = start - 1
  else
     i_a = MIN(SHAPE(array), shape_new)
  endif

  i_b = MAX(i_a, i_a+SHAPE(array)-shape_new) + 1
  i_c = MAX(i_a, i_a+shape_new-SHAPE(array)) + 1

  ! Allocate the new array

  allocate(array_new(shape_new(1),shape_new(2),shape_new(3),shape_new(4)))

  ! Copy the data from old to new

  array_new(:i_a(1),:i_a(2),:i_a(3),:i_a(4)) = array(:i_a(1),:i_a(2),:i_a(3),:i_a(4))
array_new(i_c(1):,:i_a(2),:i_a(3),:i_a(4)) = array(i_b(1):,:i_a(2),:i_a(3),:i_a(4))
array_new(:i_a(1),i_c(2):,:i_a(3),:i_a(4)) = array(:i_a(1),i_b(2):,:i_a(3),:i_a(4))
array_new(i_c(1):,i_c(2):,:i_a(3),:i_a(4)) = array(i_b(1):,i_b(2):,:i_a(3),:i_a(4))
array_new(:i_a(1),:i_a(2),i_c(3):,:i_a(4)) = array(:i_a(1),:i_a(2),i_b(3):,:i_a(4))
array_new(i_c(1):,:i_a(2),i_c(3):,:i_a(4)) = array(i_b(1):,:i_a(2),i_b(3):,:i_a(4))
array_new(:i_a(1),i_c(2):,i_c(3):,:i_a(4)) = array(:i_a(1),i_b(2):,i_b(3):,:i_a(4))
array_new(i_c(1):,i_c(2):,i_c(3):,:i_a(4)) = array(i_b(1):,i_b(2):,i_b(3):,:i_a(4))
array_new(:i_a(1),:i_a(2),:i_a(3),i_c(4):) = array(:i_a(1),:i_a(2),:i_a(3),i_b(4):)
array_new(i_c(1):,:i_a(2),:i_a(3),i_c(4):) = array(i_b(1):,:i_a(2),:i_a(3),i_b(4):)
array_new(:i_a(1),i_c(2):,:i_a(3),i_c(4):) = array(:i_a(1),i_b(2):,:i_a(3),i_b(4):)
array_new(i_c(1):,i_c(2):,:i_a(3),i_c(4):) = array(i_b(1):,i_b(2):,:i_a(3),i_b(4):)
array_new(:i_a(1),:i_a(2),i_c(3):,i_c(4):) = array(:i_a(1),:i_a(2),i_b(3):,i_b(4):)
array_new(i_c(1):,:i_a(2),i_c(3):,i_c(4):) = array(i_b(1):,:i_a(2),i_b(3):,i_b(4):)
array_new(:i_a(1),i_c(2):,i_c(3):,i_c(4):) = array(:i_a(1),i_b(2):,i_b(3):,i_b(4):)
array_new(i_c(1):,i_c(2):,i_c(3):,i_c(4):) = array(i_b(1):,i_b(2):,i_b(3):,i_b(4):)

  ! Move the allocation

  call MOVE_ALLOC(array_new, array)

  ! Finish

  return

end subroutine reallocate_r_sp_4_

subroutine reallocate_r_dp_1_ (array, shape_new, start)

  real(DP), allocatable, intent(inout) :: array(:)
  integer, intent(in)                     :: shape_new(:)
  integer, intent(in), optional           :: start(:)

  real(DP), allocatable :: array_new(:)
  integer                      :: i_a(1)
  integer                      :: i_b(1)
  integer                      :: i_c(1)

  if(PRESENT(start)) then

  end if

  ! Set up the reallocation indices

  if(PRESENT(start)) then
     i_a = start - 1
  else
     i_a = MIN(SHAPE(array), shape_new)
  endif

  i_b = MAX(i_a, i_a+SHAPE(array)-shape_new) + 1
  i_c = MAX(i_a, i_a+shape_new-SHAPE(array)) + 1

  ! Allocate the new array

  allocate(array_new(shape_new(1)))

  ! Copy the data from old to new

  array_new(:i_a(1)) = array(:i_a(1))
array_new(i_c(1):) = array(i_b(1):)

  ! Move the allocation

  call MOVE_ALLOC(array_new, array)

  ! Finish

  return

end subroutine reallocate_r_dp_1_

subroutine reallocate_r_dp_2_ (array, shape_new, start)

  real(DP), allocatable, intent(inout) :: array(:,:)
  integer, intent(in)                     :: shape_new(:)
  integer, intent(in), optional           :: start(:)

  real(DP), allocatable :: array_new(:,:)
  integer                      :: i_a(2)
  integer                      :: i_b(2)
  integer                      :: i_c(2)

  if(PRESENT(start)) then

  end if

  ! Set up the reallocation indices

  if(PRESENT(start)) then
     i_a = start - 1
  else
     i_a = MIN(SHAPE(array), shape_new)
  endif

  i_b = MAX(i_a, i_a+SHAPE(array)-shape_new) + 1
  i_c = MAX(i_a, i_a+shape_new-SHAPE(array)) + 1

  ! Allocate the new array

  allocate(array_new(shape_new(1),shape_new(2)))

  ! Copy the data from old to new

  array_new(:i_a(1),:i_a(2)) = array(:i_a(1),:i_a(2))
array_new(i_c(1):,:i_a(2)) = array(i_b(1):,:i_a(2))
array_new(:i_a(1),i_c(2):) = array(:i_a(1),i_b(2):)
array_new(i_c(1):,i_c(2):) = array(i_b(1):,i_b(2):)

  ! Move the allocation

  call MOVE_ALLOC(array_new, array)

  ! Finish

  return

end subroutine reallocate_r_dp_2_

subroutine reallocate_r_dp_3_ (array, shape_new, start)

  real(DP), allocatable, intent(inout) :: array(:,:,:)
  integer, intent(in)                     :: shape_new(:)
  integer, intent(in), optional           :: start(:)

  real(DP), allocatable :: array_new(:,:,:)
  integer                      :: i_a(3)
  integer                      :: i_b(3)
  integer                      :: i_c(3)

  if(PRESENT(start)) then

  end if

  ! Set up the reallocation indices

  if(PRESENT(start)) then
     i_a = start - 1
  else
     i_a = MIN(SHAPE(array), shape_new)
  endif

  i_b = MAX(i_a, i_a+SHAPE(array)-shape_new) + 1
  i_c = MAX(i_a, i_a+shape_new-SHAPE(array)) + 1

  ! Allocate the new array

  allocate(array_new(shape_new(1),shape_new(2),shape_new(3)))

  ! Copy the data from old to new

  array_new(:i_a(1),:i_a(2),:i_a(3)) = array(:i_a(1),:i_a(2),:i_a(3))
array_new(i_c(1):,:i_a(2),:i_a(3)) = array(i_b(1):,:i_a(2),:i_a(3))
array_new(:i_a(1),i_c(2):,:i_a(3)) = array(:i_a(1),i_b(2):,:i_a(3))
array_new(i_c(1):,i_c(2):,:i_a(3)) = array(i_b(1):,i_b(2):,:i_a(3))
array_new(:i_a(1),:i_a(2),i_c(3):) = array(:i_a(1),:i_a(2),i_b(3):)
array_new(i_c(1):,:i_a(2),i_c(3):) = array(i_b(1):,:i_a(2),i_b(3):)
array_new(:i_a(1),i_c(2):,i_c(3):) = array(:i_a(1),i_b(2):,i_b(3):)
array_new(i_c(1):,i_c(2):,i_c(3):) = array(i_b(1):,i_b(2):,i_b(3):)

  ! Move the allocation

  call MOVE_ALLOC(array_new, array)

  ! Finish

  return

end subroutine reallocate_r_dp_3_

subroutine reallocate_r_dp_4_ (array, shape_new, start)

  real(DP), allocatable, intent(inout) :: array(:,:,:,:)
  integer, intent(in)                     :: shape_new(:)
  integer, intent(in), optional           :: start(:)

  real(DP), allocatable :: array_new(:,:,:,:)
  integer                      :: i_a(4)
  integer                      :: i_b(4)
  integer                      :: i_c(4)

  if(PRESENT(start)) then

  end if

  ! Set up the reallocation indices

  if(PRESENT(start)) then
     i_a = start - 1
  else
     i_a = MIN(SHAPE(array), shape_new)
  endif

  i_b = MAX(i_a, i_a+SHAPE(array)-shape_new) + 1
  i_c = MAX(i_a, i_a+shape_new-SHAPE(array)) + 1

  ! Allocate the new array

  allocate(array_new(shape_new(1),shape_new(2),shape_new(3),shape_new(4)))

  ! Copy the data from old to new

  array_new(:i_a(1),:i_a(2),:i_a(3),:i_a(4)) = array(:i_a(1),:i_a(2),:i_a(3),:i_a(4))
array_new(i_c(1):,:i_a(2),:i_a(3),:i_a(4)) = array(i_b(1):,:i_a(2),:i_a(3),:i_a(4))
array_new(:i_a(1),i_c(2):,:i_a(3),:i_a(4)) = array(:i_a(1),i_b(2):,:i_a(3),:i_a(4))
array_new(i_c(1):,i_c(2):,:i_a(3),:i_a(4)) = array(i_b(1):,i_b(2):,:i_a(3),:i_a(4))
array_new(:i_a(1),:i_a(2),i_c(3):,:i_a(4)) = array(:i_a(1),:i_a(2),i_b(3):,:i_a(4))
array_new(i_c(1):,:i_a(2),i_c(3):,:i_a(4)) = array(i_b(1):,:i_a(2),i_b(3):,:i_a(4))
array_new(:i_a(1),i_c(2):,i_c(3):,:i_a(4)) = array(:i_a(1),i_b(2):,i_b(3):,:i_a(4))
array_new(i_c(1):,i_c(2):,i_c(3):,:i_a(4)) = array(i_b(1):,i_b(2):,i_b(3):,:i_a(4))
array_new(:i_a(1),:i_a(2),:i_a(3),i_c(4):) = array(:i_a(1),:i_a(2),:i_a(3),i_b(4):)
array_new(i_c(1):,:i_a(2),:i_a(3),i_c(4):) = array(i_b(1):,:i_a(2),:i_a(3),i_b(4):)
array_new(:i_a(1),i_c(2):,:i_a(3),i_c(4):) = array(:i_a(1),i_b(2):,:i_a(3),i_b(4):)
array_new(i_c(1):,i_c(2):,:i_a(3),i_c(4):) = array(i_b(1):,i_b(2):,:i_a(3),i_b(4):)
array_new(:i_a(1),:i_a(2),i_c(3):,i_c(4):) = array(:i_a(1),:i_a(2),i_b(3):,i_b(4):)
array_new(i_c(1):,:i_a(2),i_c(3):,i_c(4):) = array(i_b(1):,:i_a(2),i_b(3):,i_b(4):)
array_new(:i_a(1),i_c(2):,i_c(3):,i_c(4):) = array(:i_a(1),i_b(2):,i_b(3):,i_b(4):)
array_new(i_c(1):,i_c(2):,i_c(3):,i_c(4):) = array(i_b(1):,i_b(2):,i_b(3):,i_b(4):)

  ! Move the allocation

  call MOVE_ALLOC(array_new, array)

  ! Finish

  return

end subroutine reallocate_r_dp_4_

subroutine reallocate_c_sp_1_ (array, shape_new, start)

  complex(SP), allocatable, intent(inout) :: array(:)
  integer, intent(in)                     :: shape_new(:)
  integer, intent(in), optional           :: start(:)

  complex(SP), allocatable :: array_new(:)
  integer                      :: i_a(1)
  integer                      :: i_b(1)
  integer                      :: i_c(1)

  if(PRESENT(start)) then

  end if

  ! Set up the reallocation indices

  if(PRESENT(start)) then
     i_a = start - 1
  else
     i_a = MIN(SHAPE(array), shape_new)
  endif

  i_b = MAX(i_a, i_a+SHAPE(array)-shape_new) + 1
  i_c = MAX(i_a, i_a+shape_new-SHAPE(array)) + 1

  ! Allocate the new array

  allocate(array_new(shape_new(1)))

  ! Copy the data from old to new

  array_new(:i_a(1)) = array(:i_a(1))
array_new(i_c(1):) = array(i_b(1):)

  ! Move the allocation

  call MOVE_ALLOC(array_new, array)

  ! Finish

  return

end subroutine reallocate_c_sp_1_

subroutine reallocate_c_sp_2_ (array, shape_new, start)

  complex(SP), allocatable, intent(inout) :: array(:,:)
  integer, intent(in)                     :: shape_new(:)
  integer, intent(in), optional           :: start(:)

  complex(SP), allocatable :: array_new(:,:)
  integer                      :: i_a(2)
  integer                      :: i_b(2)
  integer                      :: i_c(2)

  if(PRESENT(start)) then

  end if

  ! Set up the reallocation indices

  if(PRESENT(start)) then
     i_a = start - 1
  else
     i_a = MIN(SHAPE(array), shape_new)
  endif

  i_b = MAX(i_a, i_a+SHAPE(array)-shape_new) + 1
  i_c = MAX(i_a, i_a+shape_new-SHAPE(array)) + 1

  ! Allocate the new array

  allocate(array_new(shape_new(1),shape_new(2)))

  ! Copy the data from old to new

  array_new(:i_a(1),:i_a(2)) = array(:i_a(1),:i_a(2))
array_new(i_c(1):,:i_a(2)) = array(i_b(1):,:i_a(2))
array_new(:i_a(1),i_c(2):) = array(:i_a(1),i_b(2):)
array_new(i_c(1):,i_c(2):) = array(i_b(1):,i_b(2):)

  ! Move the allocation

  call MOVE_ALLOC(array_new, array)

  ! Finish

  return

end subroutine reallocate_c_sp_2_

subroutine reallocate_c_sp_3_ (array, shape_new, start)

  complex(SP), allocatable, intent(inout) :: array(:,:,:)
  integer, intent(in)                     :: shape_new(:)
  integer, intent(in), optional           :: start(:)

  complex(SP), allocatable :: array_new(:,:,:)
  integer                      :: i_a(3)
  integer                      :: i_b(3)
  integer                      :: i_c(3)

  if(PRESENT(start)) then

  end if

  ! Set up the reallocation indices

  if(PRESENT(start)) then
     i_a = start - 1
  else
     i_a = MIN(SHAPE(array), shape_new)
  endif

  i_b = MAX(i_a, i_a+SHAPE(array)-shape_new) + 1
  i_c = MAX(i_a, i_a+shape_new-SHAPE(array)) + 1

  ! Allocate the new array

  allocate(array_new(shape_new(1),shape_new(2),shape_new(3)))

  ! Copy the data from old to new

  array_new(:i_a(1),:i_a(2),:i_a(3)) = array(:i_a(1),:i_a(2),:i_a(3))
array_new(i_c(1):,:i_a(2),:i_a(3)) = array(i_b(1):,:i_a(2),:i_a(3))
array_new(:i_a(1),i_c(2):,:i_a(3)) = array(:i_a(1),i_b(2):,:i_a(3))
array_new(i_c(1):,i_c(2):,:i_a(3)) = array(i_b(1):,i_b(2):,:i_a(3))
array_new(:i_a(1),:i_a(2),i_c(3):) = array(:i_a(1),:i_a(2),i_b(3):)
array_new(i_c(1):,:i_a(2),i_c(3):) = array(i_b(1):,:i_a(2),i_b(3):)
array_new(:i_a(1),i_c(2):,i_c(3):) = array(:i_a(1),i_b(2):,i_b(3):)
array_new(i_c(1):,i_c(2):,i_c(3):) = array(i_b(1):,i_b(2):,i_b(3):)

  ! Move the allocation

  call MOVE_ALLOC(array_new, array)

  ! Finish

  return

end subroutine reallocate_c_sp_3_

subroutine reallocate_c_sp_4_ (array, shape_new, start)

  complex(SP), allocatable, intent(inout) :: array(:,:,:,:)
  integer, intent(in)                     :: shape_new(:)
  integer, intent(in), optional           :: start(:)

  complex(SP), allocatable :: array_new(:,:,:,:)
  integer                      :: i_a(4)
  integer                      :: i_b(4)
  integer                      :: i_c(4)

  if(PRESENT(start)) then

  end if

  ! Set up the reallocation indices

  if(PRESENT(start)) then
     i_a = start - 1
  else
     i_a = MIN(SHAPE(array), shape_new)
  endif

  i_b = MAX(i_a, i_a+SHAPE(array)-shape_new) + 1
  i_c = MAX(i_a, i_a+shape_new-SHAPE(array)) + 1

  ! Allocate the new array

  allocate(array_new(shape_new(1),shape_new(2),shape_new(3),shape_new(4)))

  ! Copy the data from old to new

  array_new(:i_a(1),:i_a(2),:i_a(3),:i_a(4)) = array(:i_a(1),:i_a(2),:i_a(3),:i_a(4))
array_new(i_c(1):,:i_a(2),:i_a(3),:i_a(4)) = array(i_b(1):,:i_a(2),:i_a(3),:i_a(4))
array_new(:i_a(1),i_c(2):,:i_a(3),:i_a(4)) = array(:i_a(1),i_b(2):,:i_a(3),:i_a(4))
array_new(i_c(1):,i_c(2):,:i_a(3),:i_a(4)) = array(i_b(1):,i_b(2):,:i_a(3),:i_a(4))
array_new(:i_a(1),:i_a(2),i_c(3):,:i_a(4)) = array(:i_a(1),:i_a(2),i_b(3):,:i_a(4))
array_new(i_c(1):,:i_a(2),i_c(3):,:i_a(4)) = array(i_b(1):,:i_a(2),i_b(3):,:i_a(4))
array_new(:i_a(1),i_c(2):,i_c(3):,:i_a(4)) = array(:i_a(1),i_b(2):,i_b(3):,:i_a(4))
array_new(i_c(1):,i_c(2):,i_c(3):,:i_a(4)) = array(i_b(1):,i_b(2):,i_b(3):,:i_a(4))
array_new(:i_a(1),:i_a(2),:i_a(3),i_c(4):) = array(:i_a(1),:i_a(2),:i_a(3),i_b(4):)
array_new(i_c(1):,:i_a(2),:i_a(3),i_c(4):) = array(i_b(1):,:i_a(2),:i_a(3),i_b(4):)
array_new(:i_a(1),i_c(2):,:i_a(3),i_c(4):) = array(:i_a(1),i_b(2):,:i_a(3),i_b(4):)
array_new(i_c(1):,i_c(2):,:i_a(3),i_c(4):) = array(i_b(1):,i_b(2):,:i_a(3),i_b(4):)
array_new(:i_a(1),:i_a(2),i_c(3):,i_c(4):) = array(:i_a(1),:i_a(2),i_b(3):,i_b(4):)
array_new(i_c(1):,:i_a(2),i_c(3):,i_c(4):) = array(i_b(1):,:i_a(2),i_b(3):,i_b(4):)
array_new(:i_a(1),i_c(2):,i_c(3):,i_c(4):) = array(:i_a(1),i_b(2):,i_b(3):,i_b(4):)
array_new(i_c(1):,i_c(2):,i_c(3):,i_c(4):) = array(i_b(1):,i_b(2):,i_b(3):,i_b(4):)

  ! Move the allocation

  call MOVE_ALLOC(array_new, array)

  ! Finish

  return

end subroutine reallocate_c_sp_4_

subroutine reallocate_c_dp_1_ (array, shape_new, start)

  complex(DP), allocatable, intent(inout) :: array(:)
  integer, intent(in)                     :: shape_new(:)
  integer, intent(in), optional           :: start(:)

  complex(DP), allocatable :: array_new(:)
  integer                      :: i_a(1)
  integer                      :: i_b(1)
  integer                      :: i_c(1)

  if(PRESENT(start)) then

  end if

  ! Set up the reallocation indices

  if(PRESENT(start)) then
     i_a = start - 1
  else
     i_a = MIN(SHAPE(array), shape_new)
  endif

  i_b = MAX(i_a, i_a+SHAPE(array)-shape_new) + 1
  i_c = MAX(i_a, i_a+shape_new-SHAPE(array)) + 1

  ! Allocate the new array

  allocate(array_new(shape_new(1)))

  ! Copy the data from old to new

  array_new(:i_a(1)) = array(:i_a(1))
array_new(i_c(1):) = array(i_b(1):)

  ! Move the allocation

  call MOVE_ALLOC(array_new, array)

  ! Finish

  return

end subroutine reallocate_c_dp_1_

subroutine reallocate_c_dp_2_ (array, shape_new, start)

  complex(DP), allocatable, intent(inout) :: array(:,:)
  integer, intent(in)                     :: shape_new(:)
  integer, intent(in), optional           :: start(:)

  complex(DP), allocatable :: array_new(:,:)
  integer                      :: i_a(2)
  integer                      :: i_b(2)
  integer                      :: i_c(2)

  if(PRESENT(start)) then

  end if

  ! Set up the reallocation indices

  if(PRESENT(start)) then
     i_a = start - 1
  else
     i_a = MIN(SHAPE(array), shape_new)
  endif

  i_b = MAX(i_a, i_a+SHAPE(array)-shape_new) + 1
  i_c = MAX(i_a, i_a+shape_new-SHAPE(array)) + 1

  ! Allocate the new array

  allocate(array_new(shape_new(1),shape_new(2)))

  ! Copy the data from old to new

  array_new(:i_a(1),:i_a(2)) = array(:i_a(1),:i_a(2))
array_new(i_c(1):,:i_a(2)) = array(i_b(1):,:i_a(2))
array_new(:i_a(1),i_c(2):) = array(:i_a(1),i_b(2):)
array_new(i_c(1):,i_c(2):) = array(i_b(1):,i_b(2):)

  ! Move the allocation

  call MOVE_ALLOC(array_new, array)

  ! Finish

  return

end subroutine reallocate_c_dp_2_

subroutine reallocate_c_dp_3_ (array, shape_new, start)

  complex(DP), allocatable, intent(inout) :: array(:,:,:)
  integer, intent(in)                     :: shape_new(:)
  integer, intent(in), optional           :: start(:)

  complex(DP), allocatable :: array_new(:,:,:)
  integer                      :: i_a(3)
  integer                      :: i_b(3)
  integer                      :: i_c(3)

  if(PRESENT(start)) then

  end if

  ! Set up the reallocation indices

  if(PRESENT(start)) then
     i_a = start - 1
  else
     i_a = MIN(SHAPE(array), shape_new)
  endif

  i_b = MAX(i_a, i_a+SHAPE(array)-shape_new) + 1
  i_c = MAX(i_a, i_a+shape_new-SHAPE(array)) + 1

  ! Allocate the new array

  allocate(array_new(shape_new(1),shape_new(2),shape_new(3)))

  ! Copy the data from old to new

  array_new(:i_a(1),:i_a(2),:i_a(3)) = array(:i_a(1),:i_a(2),:i_a(3))
array_new(i_c(1):,:i_a(2),:i_a(3)) = array(i_b(1):,:i_a(2),:i_a(3))
array_new(:i_a(1),i_c(2):,:i_a(3)) = array(:i_a(1),i_b(2):,:i_a(3))
array_new(i_c(1):,i_c(2):,:i_a(3)) = array(i_b(1):,i_b(2):,:i_a(3))
array_new(:i_a(1),:i_a(2),i_c(3):) = array(:i_a(1),:i_a(2),i_b(3):)
array_new(i_c(1):,:i_a(2),i_c(3):) = array(i_b(1):,:i_a(2),i_b(3):)
array_new(:i_a(1),i_c(2):,i_c(3):) = array(:i_a(1),i_b(2):,i_b(3):)
array_new(i_c(1):,i_c(2):,i_c(3):) = array(i_b(1):,i_b(2):,i_b(3):)

  ! Move the allocation

  call MOVE_ALLOC(array_new, array)

  ! Finish

  return

end subroutine reallocate_c_dp_3_

subroutine reallocate_c_dp_4_ (array, shape_new, start)

  complex(DP), allocatable, intent(inout) :: array(:,:,:,:)
  integer, intent(in)                     :: shape_new(:)
  integer, intent(in), optional           :: start(:)

  complex(DP), allocatable :: array_new(:,:,:,:)
  integer                      :: i_a(4)
  integer                      :: i_b(4)
  integer                      :: i_c(4)

  if(PRESENT(start)) then

  end if

  ! Set up the reallocation indices

  if(PRESENT(start)) then
     i_a = start - 1
  else
     i_a = MIN(SHAPE(array), shape_new)
  endif

  i_b = MAX(i_a, i_a+SHAPE(array)-shape_new) + 1
  i_c = MAX(i_a, i_a+shape_new-SHAPE(array)) + 1

  ! Allocate the new array

  allocate(array_new(shape_new(1),shape_new(2),shape_new(3),shape_new(4)))

  ! Copy the data from old to new

  array_new(:i_a(1),:i_a(2),:i_a(3),:i_a(4)) = array(:i_a(1),:i_a(2),:i_a(3),:i_a(4))
array_new(i_c(1):,:i_a(2),:i_a(3),:i_a(4)) = array(i_b(1):,:i_a(2),:i_a(3),:i_a(4))
array_new(:i_a(1),i_c(2):,:i_a(3),:i_a(4)) = array(:i_a(1),i_b(2):,:i_a(3),:i_a(4))
array_new(i_c(1):,i_c(2):,:i_a(3),:i_a(4)) = array(i_b(1):,i_b(2):,:i_a(3),:i_a(4))
array_new(:i_a(1),:i_a(2),i_c(3):,:i_a(4)) = array(:i_a(1),:i_a(2),i_b(3):,:i_a(4))
array_new(i_c(1):,:i_a(2),i_c(3):,:i_a(4)) = array(i_b(1):,:i_a(2),i_b(3):,:i_a(4))
array_new(:i_a(1),i_c(2):,i_c(3):,:i_a(4)) = array(:i_a(1),i_b(2):,i_b(3):,:i_a(4))
array_new(i_c(1):,i_c(2):,i_c(3):,:i_a(4)) = array(i_b(1):,i_b(2):,i_b(3):,:i_a(4))
array_new(:i_a(1),:i_a(2),:i_a(3),i_c(4):) = array(:i_a(1),:i_a(2),:i_a(3),i_b(4):)
array_new(i_c(1):,:i_a(2),:i_a(3),i_c(4):) = array(i_b(1):,:i_a(2),:i_a(3),i_b(4):)
array_new(:i_a(1),i_c(2):,:i_a(3),i_c(4):) = array(:i_a(1),i_b(2):,:i_a(3),i_b(4):)
array_new(i_c(1):,i_c(2):,:i_a(3),i_c(4):) = array(i_b(1):,i_b(2):,:i_a(3),i_b(4):)
array_new(:i_a(1),:i_a(2),i_c(3):,i_c(4):) = array(:i_a(1),:i_a(2),i_b(3):,i_b(4):)
array_new(i_c(1):,:i_a(2),i_c(3):,i_c(4):) = array(i_b(1):,:i_a(2),i_b(3):,i_b(4):)
array_new(:i_a(1),i_c(2):,i_c(3):,i_c(4):) = array(:i_a(1),i_b(2):,i_b(3):,i_b(4):)
array_new(i_c(1):,i_c(2):,i_c(3):,i_c(4):) = array(i_b(1):,i_b(2):,i_b(3):,i_b(4):)

  ! Move the allocation

  call MOVE_ALLOC(array_new, array)

  ! Finish

  return

end subroutine reallocate_c_dp_4_

subroutine reallocate_a_1_ (array, shape_new, start)

  character(*), allocatable, intent(inout) :: array(:)
  integer, intent(in)                     :: shape_new(:)
  integer, intent(in), optional           :: start(:)

  character(LEN(array)), allocatable :: array_new(:)
  integer                      :: i_a(1)
  integer                      :: i_b(1)
  integer                      :: i_c(1)

  if(PRESENT(start)) then

  end if

  ! Set up the reallocation indices

  if(PRESENT(start)) then
     i_a = start - 1
  else
     i_a = MIN(SHAPE(array), shape_new)
  endif

  i_b = MAX(i_a, i_a+SHAPE(array)-shape_new) + 1
  i_c = MAX(i_a, i_a+shape_new-SHAPE(array)) + 1

  ! Allocate the new array

  allocate(array_new(shape_new(1)))

  ! Copy the data from old to new

  array_new(:i_a(1)) = array(:i_a(1))
array_new(i_c(1):) = array(i_b(1):)

  ! Move the allocation

  call MOVE_ALLOC(array_new, array)

  ! Finish

  return

end subroutine reallocate_a_1_

subroutine reallocate_a_2_ (array, shape_new, start)

  character(*), allocatable, intent(inout) :: array(:,:)
  integer, intent(in)                     :: shape_new(:)
  integer, intent(in), optional           :: start(:)

  character(LEN(array)), allocatable :: array_new(:,:)
  integer                      :: i_a(2)
  integer                      :: i_b(2)
  integer                      :: i_c(2)

  if(PRESENT(start)) then

  end if

  ! Set up the reallocation indices

  if(PRESENT(start)) then
     i_a = start - 1
  else
     i_a = MIN(SHAPE(array), shape_new)
  endif

  i_b = MAX(i_a, i_a+SHAPE(array)-shape_new) + 1
  i_c = MAX(i_a, i_a+shape_new-SHAPE(array)) + 1

  ! Allocate the new array

  allocate(array_new(shape_new(1),shape_new(2)))

  ! Copy the data from old to new

  array_new(:i_a(1),:i_a(2)) = array(:i_a(1),:i_a(2))
array_new(i_c(1):,:i_a(2)) = array(i_b(1):,:i_a(2))
array_new(:i_a(1),i_c(2):) = array(:i_a(1),i_b(2):)
array_new(i_c(1):,i_c(2):) = array(i_b(1):,i_b(2):)

  ! Move the allocation

  call MOVE_ALLOC(array_new, array)

  ! Finish

  return

end subroutine reallocate_a_2_

subroutine reallocate_a_3_ (array, shape_new, start)

  character(*), allocatable, intent(inout) :: array(:,:,:)
  integer, intent(in)                     :: shape_new(:)
  integer, intent(in), optional           :: start(:)

  character(LEN(array)), allocatable :: array_new(:,:,:)
  integer                      :: i_a(3)
  integer                      :: i_b(3)
  integer                      :: i_c(3)

  if(PRESENT(start)) then

  end if

  ! Set up the reallocation indices

  if(PRESENT(start)) then
     i_a = start - 1
  else
     i_a = MIN(SHAPE(array), shape_new)
  endif

  i_b = MAX(i_a, i_a+SHAPE(array)-shape_new) + 1
  i_c = MAX(i_a, i_a+shape_new-SHAPE(array)) + 1

  ! Allocate the new array

  allocate(array_new(shape_new(1),shape_new(2),shape_new(3)))

  ! Copy the data from old to new

  array_new(:i_a(1),:i_a(2),:i_a(3)) = array(:i_a(1),:i_a(2),:i_a(3))
array_new(i_c(1):,:i_a(2),:i_a(3)) = array(i_b(1):,:i_a(2),:i_a(3))
array_new(:i_a(1),i_c(2):,:i_a(3)) = array(:i_a(1),i_b(2):,:i_a(3))
array_new(i_c(1):,i_c(2):,:i_a(3)) = array(i_b(1):,i_b(2):,:i_a(3))
array_new(:i_a(1),:i_a(2),i_c(3):) = array(:i_a(1),:i_a(2),i_b(3):)
array_new(i_c(1):,:i_a(2),i_c(3):) = array(i_b(1):,:i_a(2),i_b(3):)
array_new(:i_a(1),i_c(2):,i_c(3):) = array(:i_a(1),i_b(2):,i_b(3):)
array_new(i_c(1):,i_c(2):,i_c(3):) = array(i_b(1):,i_b(2):,i_b(3):)

  ! Move the allocation

  call MOVE_ALLOC(array_new, array)

  ! Finish

  return

end subroutine reallocate_a_3_

subroutine reallocate_a_4_ (array, shape_new, start)

  character(*), allocatable, intent(inout) :: array(:,:,:,:)
  integer, intent(in)                     :: shape_new(:)
  integer, intent(in), optional           :: start(:)

  character(LEN(array)), allocatable :: array_new(:,:,:,:)
  integer                      :: i_a(4)
  integer                      :: i_b(4)
  integer                      :: i_c(4)

  if(PRESENT(start)) then

  end if

  ! Set up the reallocation indices

  if(PRESENT(start)) then
     i_a = start - 1
  else
     i_a = MIN(SHAPE(array), shape_new)
  endif

  i_b = MAX(i_a, i_a+SHAPE(array)-shape_new) + 1
  i_c = MAX(i_a, i_a+shape_new-SHAPE(array)) + 1

  ! Allocate the new array

  allocate(array_new(shape_new(1),shape_new(2),shape_new(3),shape_new(4)))

  ! Copy the data from old to new

  array_new(:i_a(1),:i_a(2),:i_a(3),:i_a(4)) = array(:i_a(1),:i_a(2),:i_a(3),:i_a(4))
array_new(i_c(1):,:i_a(2),:i_a(3),:i_a(4)) = array(i_b(1):,:i_a(2),:i_a(3),:i_a(4))
array_new(:i_a(1),i_c(2):,:i_a(3),:i_a(4)) = array(:i_a(1),i_b(2):,:i_a(3),:i_a(4))
array_new(i_c(1):,i_c(2):,:i_a(3),:i_a(4)) = array(i_b(1):,i_b(2):,:i_a(3),:i_a(4))
array_new(:i_a(1),:i_a(2),i_c(3):,:i_a(4)) = array(:i_a(1),:i_a(2),i_b(3):,:i_a(4))
array_new(i_c(1):,:i_a(2),i_c(3):,:i_a(4)) = array(i_b(1):,:i_a(2),i_b(3):,:i_a(4))
array_new(:i_a(1),i_c(2):,i_c(3):,:i_a(4)) = array(:i_a(1),i_b(2):,i_b(3):,:i_a(4))
array_new(i_c(1):,i_c(2):,i_c(3):,:i_a(4)) = array(i_b(1):,i_b(2):,i_b(3):,:i_a(4))
array_new(:i_a(1),:i_a(2),:i_a(3),i_c(4):) = array(:i_a(1),:i_a(2),:i_a(3),i_b(4):)
array_new(i_c(1):,:i_a(2),:i_a(3),i_c(4):) = array(i_b(1):,:i_a(2),:i_a(3),i_b(4):)
array_new(:i_a(1),i_c(2):,:i_a(3),i_c(4):) = array(:i_a(1),i_b(2):,:i_a(3),i_b(4):)
array_new(i_c(1):,i_c(2):,:i_a(3),i_c(4):) = array(i_b(1):,i_b(2):,:i_a(3),i_b(4):)
array_new(:i_a(1),:i_a(2),i_c(3):,i_c(4):) = array(:i_a(1),:i_a(2),i_b(3):,i_b(4):)
array_new(i_c(1):,:i_a(2),i_c(3):,i_c(4):) = array(i_b(1):,:i_a(2),i_b(3):,i_b(4):)
array_new(:i_a(1),i_c(2):,i_c(3):,i_c(4):) = array(:i_a(1),i_b(2):,i_b(3):,i_b(4):)
array_new(i_c(1):,i_c(2):,i_c(3):,i_c(4):) = array(i_b(1):,i_b(2):,i_b(3):,i_b(4):)

  ! Move the allocation

  call MOVE_ALLOC(array_new, array)

  ! Finish

  return

end subroutine reallocate_a_4_

subroutine reallocate_l_1_ (array, shape_new, start)

  logical, allocatable, intent(inout) :: array(:)
  integer, intent(in)                     :: shape_new(:)
  integer, intent(in), optional           :: start(:)

  logical, allocatable :: array_new(:)
  integer                      :: i_a(1)
  integer                      :: i_b(1)
  integer                      :: i_c(1)

  if(PRESENT(start)) then

  end if

  ! Set up the reallocation indices

  if(PRESENT(start)) then
     i_a = start - 1
  else
     i_a = MIN(SHAPE(array), shape_new)
  endif

  i_b = MAX(i_a, i_a+SHAPE(array)-shape_new) + 1
  i_c = MAX(i_a, i_a+shape_new-SHAPE(array)) + 1

  ! Allocate the new array

  allocate(array_new(shape_new(1)))

  ! Copy the data from old to new

  array_new(:i_a(1)) = array(:i_a(1))
array_new(i_c(1):) = array(i_b(1):)

  ! Move the allocation

  call MOVE_ALLOC(array_new, array)

  ! Finish

  return

end subroutine reallocate_l_1_

subroutine reallocate_l_2_ (array, shape_new, start)

  logical, allocatable, intent(inout) :: array(:,:)
  integer, intent(in)                     :: shape_new(:)
  integer, intent(in), optional           :: start(:)

  logical, allocatable :: array_new(:,:)
  integer                      :: i_a(2)
  integer                      :: i_b(2)
  integer                      :: i_c(2)

  if(PRESENT(start)) then

  end if

  ! Set up the reallocation indices

  if(PRESENT(start)) then
     i_a = start - 1
  else
     i_a = MIN(SHAPE(array), shape_new)
  endif

  i_b = MAX(i_a, i_a+SHAPE(array)-shape_new) + 1
  i_c = MAX(i_a, i_a+shape_new-SHAPE(array)) + 1

  ! Allocate the new array

  allocate(array_new(shape_new(1),shape_new(2)))

  ! Copy the data from old to new

  array_new(:i_a(1),:i_a(2)) = array(:i_a(1),:i_a(2))
array_new(i_c(1):,:i_a(2)) = array(i_b(1):,:i_a(2))
array_new(:i_a(1),i_c(2):) = array(:i_a(1),i_b(2):)
array_new(i_c(1):,i_c(2):) = array(i_b(1):,i_b(2):)

  ! Move the allocation

  call MOVE_ALLOC(array_new, array)

  ! Finish

  return

end subroutine reallocate_l_2_

subroutine reallocate_l_3_ (array, shape_new, start)

  logical, allocatable, intent(inout) :: array(:,:,:)
  integer, intent(in)                     :: shape_new(:)
  integer, intent(in), optional           :: start(:)

  logical, allocatable :: array_new(:,:,:)
  integer                      :: i_a(3)
  integer                      :: i_b(3)
  integer                      :: i_c(3)

  if(PRESENT(start)) then

  end if

  ! Set up the reallocation indices

  if(PRESENT(start)) then
     i_a = start - 1
  else
     i_a = MIN(SHAPE(array), shape_new)
  endif

  i_b = MAX(i_a, i_a+SHAPE(array)-shape_new) + 1
  i_c = MAX(i_a, i_a+shape_new-SHAPE(array)) + 1

  ! Allocate the new array

  allocate(array_new(shape_new(1),shape_new(2),shape_new(3)))

  ! Copy the data from old to new

  array_new(:i_a(1),:i_a(2),:i_a(3)) = array(:i_a(1),:i_a(2),:i_a(3))
array_new(i_c(1):,:i_a(2),:i_a(3)) = array(i_b(1):,:i_a(2),:i_a(3))
array_new(:i_a(1),i_c(2):,:i_a(3)) = array(:i_a(1),i_b(2):,:i_a(3))
array_new(i_c(1):,i_c(2):,:i_a(3)) = array(i_b(1):,i_b(2):,:i_a(3))
array_new(:i_a(1),:i_a(2),i_c(3):) = array(:i_a(1),:i_a(2),i_b(3):)
array_new(i_c(1):,:i_a(2),i_c(3):) = array(i_b(1):,:i_a(2),i_b(3):)
array_new(:i_a(1),i_c(2):,i_c(3):) = array(:i_a(1),i_b(2):,i_b(3):)
array_new(i_c(1):,i_c(2):,i_c(3):) = array(i_b(1):,i_b(2):,i_b(3):)

  ! Move the allocation

  call MOVE_ALLOC(array_new, array)

  ! Finish

  return

end subroutine reallocate_l_3_

subroutine reallocate_l_4_ (array, shape_new, start)

  logical, allocatable, intent(inout) :: array(:,:,:,:)
  integer, intent(in)                     :: shape_new(:)
  integer, intent(in), optional           :: start(:)

  logical, allocatable :: array_new(:,:,:,:)
  integer                      :: i_a(4)
  integer                      :: i_b(4)
  integer                      :: i_c(4)

  if(PRESENT(start)) then

  end if

  ! Set up the reallocation indices

  if(PRESENT(start)) then
     i_a = start - 1
  else
     i_a = MIN(SHAPE(array), shape_new)
  endif

  i_b = MAX(i_a, i_a+SHAPE(array)-shape_new) + 1
  i_c = MAX(i_a, i_a+shape_new-SHAPE(array)) + 1

  ! Allocate the new array

  allocate(array_new(shape_new(1),shape_new(2),shape_new(3),shape_new(4)))

  ! Copy the data from old to new

  array_new(:i_a(1),:i_a(2),:i_a(3),:i_a(4)) = array(:i_a(1),:i_a(2),:i_a(3),:i_a(4))
array_new(i_c(1):,:i_a(2),:i_a(3),:i_a(4)) = array(i_b(1):,:i_a(2),:i_a(3),:i_a(4))
array_new(:i_a(1),i_c(2):,:i_a(3),:i_a(4)) = array(:i_a(1),i_b(2):,:i_a(3),:i_a(4))
array_new(i_c(1):,i_c(2):,:i_a(3),:i_a(4)) = array(i_b(1):,i_b(2):,:i_a(3),:i_a(4))
array_new(:i_a(1),:i_a(2),i_c(3):,:i_a(4)) = array(:i_a(1),:i_a(2),i_b(3):,:i_a(4))
array_new(i_c(1):,:i_a(2),i_c(3):,:i_a(4)) = array(i_b(1):,:i_a(2),i_b(3):,:i_a(4))
array_new(:i_a(1),i_c(2):,i_c(3):,:i_a(4)) = array(:i_a(1),i_b(2):,i_b(3):,:i_a(4))
array_new(i_c(1):,i_c(2):,i_c(3):,:i_a(4)) = array(i_b(1):,i_b(2):,i_b(3):,:i_a(4))
array_new(:i_a(1),:i_a(2),:i_a(3),i_c(4):) = array(:i_a(1),:i_a(2),:i_a(3),i_b(4):)
array_new(i_c(1):,:i_a(2),:i_a(3),i_c(4):) = array(i_b(1):,:i_a(2),:i_a(3),i_b(4):)
array_new(:i_a(1),i_c(2):,:i_a(3),i_c(4):) = array(:i_a(1),i_b(2):,:i_a(3),i_b(4):)
array_new(i_c(1):,i_c(2):,:i_a(3),i_c(4):) = array(i_b(1):,i_b(2):,:i_a(3),i_b(4):)
array_new(:i_a(1),:i_a(2),i_c(3):,i_c(4):) = array(:i_a(1),:i_a(2),i_b(3):,i_b(4):)
array_new(i_c(1):,:i_a(2),i_c(3):,i_c(4):) = array(i_b(1):,:i_a(2),i_b(3):,i_b(4):)
array_new(:i_a(1),i_c(2):,i_c(3):,i_c(4):) = array(:i_a(1),i_b(2):,i_b(3):,i_b(4):)
array_new(i_c(1):,i_c(2):,i_c(3):,i_c(4):) = array(i_b(1):,i_b(2):,i_b(3):,i_b(4):)

  ! Move the allocation

  call MOVE_ALLOC(array_new, array)

  ! Finish

  return

end subroutine reallocate_l_4_

end module core_memory

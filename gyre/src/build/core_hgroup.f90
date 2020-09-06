!fpx3_header(0.13_3a)
!
!dependencies
!   dir: ~/Documents/MESA_SS/gyre/src/build 
!   sources: -
!   includes: ../extern/core/core.inc
!   uses: hdf5 ISO_C_BINDING ISO_FORTRAN_ENV core_kinds
!   provides: core_hgroup
!end dependencies
!
!end fpx3_header
! Module   : core_hgroup
! Purpose  : HDF5 input/output

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

module core_hgroup

  ! Uses

  use core_kinds

  use hdf5

  use ISO_FORTRAN_ENV
  use ISO_C_BINDING

  ! No implicit typing

  implicit none

  ! Parameters

  integer, parameter :: CREATE_FILE = 1
  integer, parameter :: OPEN_FILE = 2

  ! Derived-type definitions

  type hgroup_t
     private
     integer(HID_T)   :: file_id = 0
     integer(HID_T)   :: group_id = 0
     integer, pointer :: ref_count => null()
   contains
     private
     procedure, public :: final => final_
     procedure, public :: exists => exists_
  end type hgroup_t

  ! Module variables

  integer, save        :: ref_count = 0
  integer(HID_T), save :: c_sp_mem_type_id
  integer(HID_T), save :: c_dp_mem_type_id
  integer(HID_T), save :: c_sp_file_type_id
  integer(HID_T), save :: c_dp_file_type_id

  ! Interfaces

  interface hgroup_t
     module procedure init_file_
     module procedure init_group_
  end interface hgroup_t

  interface read_attr
     module procedure read_attr_i_i4_0_
     module procedure read_attr_i_i4_1_
     module procedure read_attr_i_i4_2_
     module procedure read_attr_i_i4_3_
     module procedure read_attr_i_i4_4_
     module procedure read_attr_i_i8_0_
     module procedure read_attr_i_i8_1_
     module procedure read_attr_i_i8_2_
     module procedure read_attr_i_i8_3_
     module procedure read_attr_i_i8_4_
     module procedure read_attr_r_sp_0_
     module procedure read_attr_r_sp_1_
     module procedure read_attr_r_sp_2_
     module procedure read_attr_r_sp_3_
     module procedure read_attr_r_sp_4_
     module procedure read_attr_r_dp_0_
     module procedure read_attr_r_dp_1_
     module procedure read_attr_r_dp_2_
     module procedure read_attr_r_dp_3_
     module procedure read_attr_r_dp_4_
     module procedure read_attr_c_sp_0_
     module procedure read_attr_c_sp_1_
     module procedure read_attr_c_sp_2_
     module procedure read_attr_c_sp_3_
     module procedure read_attr_c_sp_4_
     module procedure read_attr_c_dp_0_
     module procedure read_attr_c_dp_1_
     module procedure read_attr_c_dp_2_
     module procedure read_attr_c_dp_3_
     module procedure read_attr_c_dp_4_
     module procedure read_attr_a_0_
     module procedure read_attr_a_1_
     module procedure read_attr_a_2_
     module procedure read_attr_a_3_
     module procedure read_attr_a_4_
     module procedure read_attr_l_0_
     module procedure read_attr_l_1_
     module procedure read_attr_l_2_
     module procedure read_attr_l_3_
     module procedure read_attr_l_4_
  end interface read_attr

  interface read_dset
     module procedure read_dset_i_i4_0_
     module procedure read_dset_i_i4_1_
     module procedure read_dset_i_i4_2_
     module procedure read_dset_i_i4_3_
     module procedure read_dset_i_i4_4_
     module procedure read_dset_i_i8_0_
     module procedure read_dset_i_i8_1_
     module procedure read_dset_i_i8_2_
     module procedure read_dset_i_i8_3_
     module procedure read_dset_i_i8_4_
     module procedure read_dset_r_sp_0_
     module procedure read_dset_r_sp_1_
     module procedure read_dset_r_sp_2_
     module procedure read_dset_r_sp_3_
     module procedure read_dset_r_sp_4_
     module procedure read_dset_r_dp_0_
     module procedure read_dset_r_dp_1_
     module procedure read_dset_r_dp_2_
     module procedure read_dset_r_dp_3_
     module procedure read_dset_r_dp_4_
     module procedure read_dset_c_sp_0_
     module procedure read_dset_c_sp_1_
     module procedure read_dset_c_sp_2_
     module procedure read_dset_c_sp_3_
     module procedure read_dset_c_sp_4_
     module procedure read_dset_c_dp_0_
     module procedure read_dset_c_dp_1_
     module procedure read_dset_c_dp_2_
     module procedure read_dset_c_dp_3_
     module procedure read_dset_c_dp_4_
     module procedure read_dset_a_0_
     module procedure read_dset_a_1_
     module procedure read_dset_a_2_
     module procedure read_dset_a_3_
     module procedure read_dset_a_4_
     module procedure read_dset_l_0_
     module procedure read_dset_l_1_
     module procedure read_dset_l_2_
     module procedure read_dset_l_3_
     module procedure read_dset_l_4_
  end interface read_dset

  interface read_attr_alloc
     module procedure read_attr_alloc_i_i4_0_
     module procedure read_attr_alloc_i_i4_1_
     module procedure read_attr_alloc_i_i4_2_
     module procedure read_attr_alloc_i_i4_3_
     module procedure read_attr_alloc_i_i4_4_
     module procedure read_attr_alloc_i_i8_0_
     module procedure read_attr_alloc_i_i8_1_
     module procedure read_attr_alloc_i_i8_2_
     module procedure read_attr_alloc_i_i8_3_
     module procedure read_attr_alloc_i_i8_4_
     module procedure read_attr_alloc_r_sp_0_
     module procedure read_attr_alloc_r_sp_1_
     module procedure read_attr_alloc_r_sp_2_
     module procedure read_attr_alloc_r_sp_3_
     module procedure read_attr_alloc_r_sp_4_
     module procedure read_attr_alloc_r_dp_0_
     module procedure read_attr_alloc_r_dp_1_
     module procedure read_attr_alloc_r_dp_2_
     module procedure read_attr_alloc_r_dp_3_
     module procedure read_attr_alloc_r_dp_4_
     module procedure read_attr_alloc_c_sp_0_
     module procedure read_attr_alloc_c_sp_1_
     module procedure read_attr_alloc_c_sp_2_
     module procedure read_attr_alloc_c_sp_3_
     module procedure read_attr_alloc_c_sp_4_
     module procedure read_attr_alloc_c_dp_0_
     module procedure read_attr_alloc_c_dp_1_
     module procedure read_attr_alloc_c_dp_2_
     module procedure read_attr_alloc_c_dp_3_
     module procedure read_attr_alloc_c_dp_4_
     module procedure read_attr_alloc_a_0_
     module procedure read_attr_alloc_a_1_
     module procedure read_attr_alloc_a_2_
     module procedure read_attr_alloc_a_3_
     module procedure read_attr_alloc_a_4_
     module procedure read_attr_alloc_l_0_
     module procedure read_attr_alloc_l_1_
     module procedure read_attr_alloc_l_2_
     module procedure read_attr_alloc_l_3_
     module procedure read_attr_alloc_l_4_
  end interface read_attr_alloc

  interface read_dset_alloc
     module procedure read_dset_alloc_i_i4_0_
     module procedure read_dset_alloc_i_i4_1_
     module procedure read_dset_alloc_i_i4_2_
     module procedure read_dset_alloc_i_i4_3_
     module procedure read_dset_alloc_i_i4_4_
     module procedure read_dset_alloc_i_i8_0_
     module procedure read_dset_alloc_i_i8_1_
     module procedure read_dset_alloc_i_i8_2_
     module procedure read_dset_alloc_i_i8_3_
     module procedure read_dset_alloc_i_i8_4_
     module procedure read_dset_alloc_r_sp_0_
     module procedure read_dset_alloc_r_sp_1_
     module procedure read_dset_alloc_r_sp_2_
     module procedure read_dset_alloc_r_sp_3_
     module procedure read_dset_alloc_r_sp_4_
     module procedure read_dset_alloc_r_dp_0_
     module procedure read_dset_alloc_r_dp_1_
     module procedure read_dset_alloc_r_dp_2_
     module procedure read_dset_alloc_r_dp_3_
     module procedure read_dset_alloc_r_dp_4_
     module procedure read_dset_alloc_c_sp_0_
     module procedure read_dset_alloc_c_sp_1_
     module procedure read_dset_alloc_c_sp_2_
     module procedure read_dset_alloc_c_sp_3_
     module procedure read_dset_alloc_c_sp_4_
     module procedure read_dset_alloc_c_dp_0_
     module procedure read_dset_alloc_c_dp_1_
     module procedure read_dset_alloc_c_dp_2_
     module procedure read_dset_alloc_c_dp_3_
     module procedure read_dset_alloc_c_dp_4_
     module procedure read_dset_alloc_a_0_
     module procedure read_dset_alloc_a_1_
     module procedure read_dset_alloc_a_2_
     module procedure read_dset_alloc_a_3_
     module procedure read_dset_alloc_a_4_
     module procedure read_dset_alloc_l_0_
     module procedure read_dset_alloc_l_1_
     module procedure read_dset_alloc_l_2_
     module procedure read_dset_alloc_l_3_
     module procedure read_dset_alloc_l_4_
  end interface read_dset_alloc

  interface write_attr
     module procedure write_attr_i_i4_0_
     module procedure write_attr_i_i4_1_
     module procedure write_attr_i_i4_2_
     module procedure write_attr_i_i4_3_
     module procedure write_attr_i_i4_4_
     module procedure write_attr_i_i8_0_
     module procedure write_attr_i_i8_1_
     module procedure write_attr_i_i8_2_
     module procedure write_attr_i_i8_3_
     module procedure write_attr_i_i8_4_
     module procedure write_attr_r_sp_0_
     module procedure write_attr_r_sp_1_
     module procedure write_attr_r_sp_2_
     module procedure write_attr_r_sp_3_
     module procedure write_attr_r_sp_4_
     module procedure write_attr_r_dp_0_
     module procedure write_attr_r_dp_1_
     module procedure write_attr_r_dp_2_
     module procedure write_attr_r_dp_3_
     module procedure write_attr_r_dp_4_
     module procedure write_attr_c_sp_0_
     module procedure write_attr_c_sp_1_
     module procedure write_attr_c_sp_2_
     module procedure write_attr_c_sp_3_
     module procedure write_attr_c_sp_4_
     module procedure write_attr_c_dp_0_
     module procedure write_attr_c_dp_1_
     module procedure write_attr_c_dp_2_
     module procedure write_attr_c_dp_3_
     module procedure write_attr_c_dp_4_
     module procedure write_attr_a_0_
     module procedure write_attr_a_1_
     module procedure write_attr_a_2_
     module procedure write_attr_a_3_
     module procedure write_attr_a_4_
     module procedure write_attr_l_0_
     module procedure write_attr_l_1_
     module procedure write_attr_l_2_
     module procedure write_attr_l_3_
     module procedure write_attr_l_4_
  end interface write_attr

  interface write_dset
     module procedure write_dset_i_i4_0_
     module procedure write_dset_i_i4_1_
     module procedure write_dset_i_i4_2_
     module procedure write_dset_i_i4_3_
     module procedure write_dset_i_i4_4_
     module procedure write_dset_i_i8_0_
     module procedure write_dset_i_i8_1_
     module procedure write_dset_i_i8_2_
     module procedure write_dset_i_i8_3_
     module procedure write_dset_i_i8_4_
     module procedure write_dset_r_sp_0_
     module procedure write_dset_r_sp_1_
     module procedure write_dset_r_sp_2_
     module procedure write_dset_r_sp_3_
     module procedure write_dset_r_sp_4_
     module procedure write_dset_r_dp_0_
     module procedure write_dset_r_dp_1_
     module procedure write_dset_r_dp_2_
     module procedure write_dset_r_dp_3_
     module procedure write_dset_r_dp_4_
     module procedure write_dset_c_sp_0_
     module procedure write_dset_c_sp_1_
     module procedure write_dset_c_sp_2_
     module procedure write_dset_c_sp_3_
     module procedure write_dset_c_sp_4_
     module procedure write_dset_c_dp_0_
     module procedure write_dset_c_dp_1_
     module procedure write_dset_c_dp_2_
     module procedure write_dset_c_dp_3_
     module procedure write_dset_c_dp_4_
     module procedure write_dset_a_0_
     module procedure write_dset_a_1_
     module procedure write_dset_a_2_
     module procedure write_dset_a_3_
     module procedure write_dset_a_4_
     module procedure write_dset_l_0_
     module procedure write_dset_l_1_
     module procedure write_dset_l_2_
     module procedure write_dset_l_3_
     module procedure write_dset_l_4_
  end interface write_dset

  interface attr_exists
     module procedure attr_exists_
  end interface attr_exists

  interface dset_exists
     module procedure dset_exists_
  end interface dset_exists

  ! Access specifiers

  private

  public :: CREATE_FILE
  public :: OPEN_FILE
  public :: hgroup_t
  public :: read_attr
  public :: read_dset
  public :: read_attr_alloc
  public :: read_dset_alloc
  public :: write_attr
  public :: write_dset
  public :: get_attr_shape
  public :: get_dset_shape
  public :: elem_group_name
  public :: attr_exists
  public :: dset_exists

  ! Procedures

contains

  function init_file_ (file_name, access_type) result (hg)

    character(*), intent(in) :: file_name
    integer, intent(in)      :: access_type
    type(hgroup_t)           :: hg

    integer        :: hdf_err
    integer(HID_T) :: file_id
    integer(HID_T) :: group_id

    ! If necessary, open the HDF5 library

    if(ref_count == 0) then
       call open_library_()
    endif

    ref_count = ref_count + 1

    ! Depending on the access_type, open or create the file

    select case(access_type)
    case(CREATE_FILE)
       call h5fcreate_f( file_name,  H5F_ACC_TRUNC_F,  file_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 386 <core_hgroup:init_file_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5fcreate_f'
   stop 'Program aborted'
endif

    case(OPEN_FILE)
       call h5fopen_f( file_name,  H5F_ACC_RDWR_F,  file_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 388 <core_hgroup:init_file_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5fopen_f'
   stop 'Program aborted'
endif

    case default

    write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 390 <core_hgroup:init_file_>:'
    write(UNIT=ERROR_UNIT, FMT=*) 'Invalid access_type'

  stop 'Program aborted'

    end select

    ! Open the root group

    call h5gopen_f( file_id,  '/',  group_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 395 <core_hgroup:init_file_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5gopen_f'
   stop 'Program aborted'
endif

    ! Construct the hgroup_t

    hg%file_id = file_id
    hg%group_id = group_id

    allocate(hg%ref_count)
    hg%ref_count = 1

    ! Finish

    return

  end function init_file_

!****

  function init_group_ (hg_parent, group_name) result (hg)

    type(hgroup_t), intent(inout) :: hg_parent
    character(*), intent(in)      :: group_name
    type(hgroup_t)                :: hg

    integer        :: hdf_err
    integer(HID_T) :: group_id

    ! Depending on whether the group already exists, open or create it

    if(hg_parent%exists(group_name)) then
       call h5gopen_f( hg_parent%group_id,  group_name,  group_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 425 <core_hgroup:init_group_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5gopen_f'
   stop 'Program aborted'
endif

    else
       call h5gcreate_f( hg_parent%group_id,  group_name,  group_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 427 <core_hgroup:init_group_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5gcreate_f'
   stop 'Program aborted'
endif

    endif

    ! Construct the hgroup_t

    hg%file_id = hg_parent%file_id
    hg%group_id = group_id

    hg%ref_count => hg_parent%ref_count
    hg%ref_count = hg%ref_count + 1

    ! Finish

    return

  end function init_group_

!****

  subroutine open_library_ ()

    integer :: hdf_err

    ! Open the HDF5 library

    call h5open_f(hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 452 <core_hgroup:open_library_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5open_f'
   stop 'Program aborted'
endif

    ! Create complex data types

    call create_complex_type_(real_mem_type_(SP), c_sp_mem_type_id)
    call create_complex_type_(real_mem_type_(DP), c_dp_mem_type_id)

    call create_complex_type_(real_mem_type_(SP), c_sp_file_type_id)
    call create_complex_type_(real_mem_type_(DP), c_dp_file_type_id)

    ! Finish

    return

  contains

    subroutine create_complex_type_ (comp_type_id, type_id)

      integer(HID_T), intent(in)  :: comp_type_id
      integer(HID_T), intent(out) :: type_id

      integer         :: hdf_err
      integer(SIZE_T) :: comp_size

      ! Create a complex data type

      call h5tget_size_f( comp_type_id,  comp_size, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 478 <core_hgroup:open_library_:create_complex_type_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5tget_size_f'
   stop 'Program aborted'
endif

      call h5tcreate_f( H5T_COMPOUND_F,  INT(2*comp_size, SIZE_T),  type_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 480 <core_hgroup:open_library_:create_complex_type_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5tcreate_f'
   stop 'Program aborted'
endif

      call h5tinsert_f( type_id,  're',  INT(0, SIZE_T),  comp_type_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 481 <core_hgroup:open_library_:create_complex_type_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5tinsert_f'
   stop 'Program aborted'
endif

      call h5tinsert_f( type_id,  'im',  INT(comp_size, SIZE_T),  comp_type_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 482 <core_hgroup:open_library_:create_complex_type_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5tinsert_f'
   stop 'Program aborted'
endif

      ! Finish

      return

    end subroutine create_complex_type_

  end subroutine open_library_

!****

  subroutine final_ (this)

    class(hgroup_t), intent(inout) :: this

    integer :: hdf_err

    ! Close the group

    call h5gclose_f( this%group_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 502 <core_hgroup:final_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5gclose_f'
   stop 'Program aborted'
endif

    this%ref_count = this%ref_count - 1

    ! If necessary, close the file also

    if(this%ref_count == 0) then
       call h5fclose_f( this%file_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 509 <core_hgroup:final_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5fclose_f'
   stop 'Program aborted'
endif

       deallocate(this%ref_count)
       ref_count = ref_count - 1
    endif

    ! If necessary, close the HDF5 library

    if(ref_count == 0) then
       call close_library_()
    endif

    ! Finish

    return

  end subroutine final_

!****

  subroutine close_library_ ()

    integer :: hdf_err

    ! Close complex data types

    call h5tclose_f( c_sp_mem_type_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 534 <core_hgroup:close_library_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5tclose_f'
   stop 'Program aborted'
endif

    call h5tclose_f( c_dp_mem_type_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 535 <core_hgroup:close_library_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5tclose_f'
   stop 'Program aborted'
endif

    call h5tclose_f( c_sp_file_type_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 537 <core_hgroup:close_library_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5tclose_f'
   stop 'Program aborted'
endif

    call h5tclose_f( c_dp_file_type_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 538 <core_hgroup:close_library_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5tclose_f'
   stop 'Program aborted'
endif

    ! Close the HDF5 library

    call h5close_f(hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 542 <core_hgroup:close_library_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5close_f'
   stop 'Program aborted'
endif

    ! Finish

    return

  end subroutine close_library_

!****

  function exists_ (this, group_name)

    class(hgroup_t), intent(inout) :: this
    character(*), intent(in)       :: group_name
    logical                        :: exists_

    integer        :: hdf_err
    integer(HID_T) :: group_id

    ! Determine whether the named group already exists

    call h5eset_auto_f(0, hdf_err)

    call h5gopen_f(this%group_id, group_name, group_id, hdf_err)

    if(hdf_err >= 0) then
       exists_ = .TRUE.
       call h5gclose_f(group_id, hdf_err)
    else
       exists_ = .FALSE.
    endif

    call h5eset_auto_f(1, hdf_err)

    ! Finish

    return

  end function exists_

!****

  subroutine read_attr_i_i4_0_ (hg, item_name, data)

    type(hgroup_t), intent(inout)               :: hg
    character(*), intent(in)                    :: item_name
    integer(I4), target, intent(out) :: data

    integer(HID_T)   :: mem_type_id
    integer          :: hdf_err
    integer(HID_T)   :: item_id
    type(C_PTR)      :: data_ptr

    ! Check the shape of the item agrees with that of data

    ! Read the item

    mem_type_id = integer_mem_type_(I4)

    data_ptr = C_LOC(data)

    call h5aopen_name_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 656 <core_hgroup:read_attr_i_i4_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aopen_name_f'
   stop 'Program aborted'
endif

    call h5aread_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 656 <core_hgroup:read_attr_i_i4_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aread_f'
   stop 'Program aborted'
endif

    call h5aclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 656 <core_hgroup:read_attr_i_i4_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aclose_f'
   stop 'Program aborted'
endif

    ! Finish

    return

  end subroutine read_attr_i_i4_0_

  subroutine read_attr_i_i4_1_ (hg, item_name, data)

    type(hgroup_t), intent(inout)               :: hg
    character(*), intent(in)                    :: item_name
    integer(I4), target, intent(out) :: data(:)

    contiguous :: data

    integer(HSIZE_T) :: item_shape(1)

    integer(HID_T)   :: mem_type_id
    integer          :: hdf_err
    integer(HID_T)   :: item_id
    type(C_PTR)      :: data_ptr

    ! Check the shape of the item agrees with that of data

    call get_attr_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 657 <core_hgroup:read_attr_i_i4_1_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

    ! Read the item

    mem_type_id = integer_mem_type_(I4)

    data_ptr = C_LOC(data)

    call h5aopen_name_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 657 <core_hgroup:read_attr_i_i4_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aopen_name_f'
   stop 'Program aborted'
endif

    call h5aread_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 657 <core_hgroup:read_attr_i_i4_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aread_f'
   stop 'Program aborted'
endif

    call h5aclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 657 <core_hgroup:read_attr_i_i4_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aclose_f'
   stop 'Program aborted'
endif

    ! Finish

    return

  end subroutine read_attr_i_i4_1_

  subroutine read_attr_i_i4_2_ (hg, item_name, data)

    type(hgroup_t), intent(inout)               :: hg
    character(*), intent(in)                    :: item_name
    integer(I4), target, intent(out) :: data(:,:)

    contiguous :: data

    integer(HSIZE_T) :: item_shape(2)

    integer(HID_T)   :: mem_type_id
    integer          :: hdf_err
    integer(HID_T)   :: item_id
    type(C_PTR)      :: data_ptr

    ! Check the shape of the item agrees with that of data

    call get_attr_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 658 <core_hgroup:read_attr_i_i4_2_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

    ! Read the item

    mem_type_id = integer_mem_type_(I4)

    data_ptr = C_LOC(data)

    call h5aopen_name_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 658 <core_hgroup:read_attr_i_i4_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aopen_name_f'
   stop 'Program aborted'
endif

    call h5aread_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 658 <core_hgroup:read_attr_i_i4_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aread_f'
   stop 'Program aborted'
endif

    call h5aclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 658 <core_hgroup:read_attr_i_i4_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aclose_f'
   stop 'Program aborted'
endif

    ! Finish

    return

  end subroutine read_attr_i_i4_2_

  subroutine read_attr_i_i4_3_ (hg, item_name, data)

    type(hgroup_t), intent(inout)               :: hg
    character(*), intent(in)                    :: item_name
    integer(I4), target, intent(out) :: data(:,:,:)

    contiguous :: data

    integer(HSIZE_T) :: item_shape(3)

    integer(HID_T)   :: mem_type_id
    integer          :: hdf_err
    integer(HID_T)   :: item_id
    type(C_PTR)      :: data_ptr

    ! Check the shape of the item agrees with that of data

    call get_attr_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 659 <core_hgroup:read_attr_i_i4_3_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

    ! Read the item

    mem_type_id = integer_mem_type_(I4)

    data_ptr = C_LOC(data)

    call h5aopen_name_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 659 <core_hgroup:read_attr_i_i4_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aopen_name_f'
   stop 'Program aborted'
endif

    call h5aread_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 659 <core_hgroup:read_attr_i_i4_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aread_f'
   stop 'Program aborted'
endif

    call h5aclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 659 <core_hgroup:read_attr_i_i4_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aclose_f'
   stop 'Program aborted'
endif

    ! Finish

    return

  end subroutine read_attr_i_i4_3_

  subroutine read_attr_i_i4_4_ (hg, item_name, data)

    type(hgroup_t), intent(inout)               :: hg
    character(*), intent(in)                    :: item_name
    integer(I4), target, intent(out) :: data(:,:,:,:)

    contiguous :: data

    integer(HSIZE_T) :: item_shape(4)

    integer(HID_T)   :: mem_type_id
    integer          :: hdf_err
    integer(HID_T)   :: item_id
    type(C_PTR)      :: data_ptr

    ! Check the shape of the item agrees with that of data

    call get_attr_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 660 <core_hgroup:read_attr_i_i4_4_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

    ! Read the item

    mem_type_id = integer_mem_type_(I4)

    data_ptr = C_LOC(data)

    call h5aopen_name_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 660 <core_hgroup:read_attr_i_i4_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aopen_name_f'
   stop 'Program aborted'
endif

    call h5aread_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 660 <core_hgroup:read_attr_i_i4_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aread_f'
   stop 'Program aborted'
endif

    call h5aclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 660 <core_hgroup:read_attr_i_i4_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aclose_f'
   stop 'Program aborted'
endif

    ! Finish

    return

  end subroutine read_attr_i_i4_4_

  subroutine read_attr_i_i8_0_ (hg, item_name, data)

    type(hgroup_t), intent(inout)               :: hg
    character(*), intent(in)                    :: item_name
    integer(I8), target, intent(out) :: data

    integer(HID_T)   :: mem_type_id
    integer          :: hdf_err
    integer(HID_T)   :: item_id
    type(C_PTR)      :: data_ptr

    ! Check the shape of the item agrees with that of data

    ! Read the item

    mem_type_id = integer_mem_type_(I8)

    data_ptr = C_LOC(data)

    call h5aopen_name_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 662 <core_hgroup:read_attr_i_i8_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aopen_name_f'
   stop 'Program aborted'
endif

    call h5aread_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 662 <core_hgroup:read_attr_i_i8_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aread_f'
   stop 'Program aborted'
endif

    call h5aclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 662 <core_hgroup:read_attr_i_i8_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aclose_f'
   stop 'Program aborted'
endif

    ! Finish

    return

  end subroutine read_attr_i_i8_0_

  subroutine read_attr_i_i8_1_ (hg, item_name, data)

    type(hgroup_t), intent(inout)               :: hg
    character(*), intent(in)                    :: item_name
    integer(I8), target, intent(out) :: data(:)

    contiguous :: data

    integer(HSIZE_T) :: item_shape(1)

    integer(HID_T)   :: mem_type_id
    integer          :: hdf_err
    integer(HID_T)   :: item_id
    type(C_PTR)      :: data_ptr

    ! Check the shape of the item agrees with that of data

    call get_attr_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 663 <core_hgroup:read_attr_i_i8_1_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

    ! Read the item

    mem_type_id = integer_mem_type_(I8)

    data_ptr = C_LOC(data)

    call h5aopen_name_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 663 <core_hgroup:read_attr_i_i8_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aopen_name_f'
   stop 'Program aborted'
endif

    call h5aread_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 663 <core_hgroup:read_attr_i_i8_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aread_f'
   stop 'Program aborted'
endif

    call h5aclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 663 <core_hgroup:read_attr_i_i8_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aclose_f'
   stop 'Program aborted'
endif

    ! Finish

    return

  end subroutine read_attr_i_i8_1_

  subroutine read_attr_i_i8_2_ (hg, item_name, data)

    type(hgroup_t), intent(inout)               :: hg
    character(*), intent(in)                    :: item_name
    integer(I8), target, intent(out) :: data(:,:)

    contiguous :: data

    integer(HSIZE_T) :: item_shape(2)

    integer(HID_T)   :: mem_type_id
    integer          :: hdf_err
    integer(HID_T)   :: item_id
    type(C_PTR)      :: data_ptr

    ! Check the shape of the item agrees with that of data

    call get_attr_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 664 <core_hgroup:read_attr_i_i8_2_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

    ! Read the item

    mem_type_id = integer_mem_type_(I8)

    data_ptr = C_LOC(data)

    call h5aopen_name_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 664 <core_hgroup:read_attr_i_i8_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aopen_name_f'
   stop 'Program aborted'
endif

    call h5aread_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 664 <core_hgroup:read_attr_i_i8_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aread_f'
   stop 'Program aborted'
endif

    call h5aclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 664 <core_hgroup:read_attr_i_i8_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aclose_f'
   stop 'Program aborted'
endif

    ! Finish

    return

  end subroutine read_attr_i_i8_2_

  subroutine read_attr_i_i8_3_ (hg, item_name, data)

    type(hgroup_t), intent(inout)               :: hg
    character(*), intent(in)                    :: item_name
    integer(I8), target, intent(out) :: data(:,:,:)

    contiguous :: data

    integer(HSIZE_T) :: item_shape(3)

    integer(HID_T)   :: mem_type_id
    integer          :: hdf_err
    integer(HID_T)   :: item_id
    type(C_PTR)      :: data_ptr

    ! Check the shape of the item agrees with that of data

    call get_attr_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 665 <core_hgroup:read_attr_i_i8_3_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

    ! Read the item

    mem_type_id = integer_mem_type_(I8)

    data_ptr = C_LOC(data)

    call h5aopen_name_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 665 <core_hgroup:read_attr_i_i8_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aopen_name_f'
   stop 'Program aborted'
endif

    call h5aread_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 665 <core_hgroup:read_attr_i_i8_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aread_f'
   stop 'Program aborted'
endif

    call h5aclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 665 <core_hgroup:read_attr_i_i8_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aclose_f'
   stop 'Program aborted'
endif

    ! Finish

    return

  end subroutine read_attr_i_i8_3_

  subroutine read_attr_i_i8_4_ (hg, item_name, data)

    type(hgroup_t), intent(inout)               :: hg
    character(*), intent(in)                    :: item_name
    integer(I8), target, intent(out) :: data(:,:,:,:)

    contiguous :: data

    integer(HSIZE_T) :: item_shape(4)

    integer(HID_T)   :: mem_type_id
    integer          :: hdf_err
    integer(HID_T)   :: item_id
    type(C_PTR)      :: data_ptr

    ! Check the shape of the item agrees with that of data

    call get_attr_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 666 <core_hgroup:read_attr_i_i8_4_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

    ! Read the item

    mem_type_id = integer_mem_type_(I8)

    data_ptr = C_LOC(data)

    call h5aopen_name_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 666 <core_hgroup:read_attr_i_i8_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aopen_name_f'
   stop 'Program aborted'
endif

    call h5aread_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 666 <core_hgroup:read_attr_i_i8_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aread_f'
   stop 'Program aborted'
endif

    call h5aclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 666 <core_hgroup:read_attr_i_i8_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aclose_f'
   stop 'Program aborted'
endif

    ! Finish

    return

  end subroutine read_attr_i_i8_4_

  subroutine read_attr_r_sp_0_ (hg, item_name, data)

    type(hgroup_t), intent(inout)               :: hg
    character(*), intent(in)                    :: item_name
    real(SP), target, intent(out) :: data

    integer(HID_T)   :: mem_type_id
    integer          :: hdf_err
    integer(HID_T)   :: item_id
    type(C_PTR)      :: data_ptr

    ! Check the shape of the item agrees with that of data

    ! Read the item

    mem_type_id = real_mem_type_(SP)

    data_ptr = C_LOC(data)

    call h5aopen_name_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 668 <core_hgroup:read_attr_r_sp_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aopen_name_f'
   stop 'Program aborted'
endif

    call h5aread_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 668 <core_hgroup:read_attr_r_sp_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aread_f'
   stop 'Program aborted'
endif

    call h5aclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 668 <core_hgroup:read_attr_r_sp_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aclose_f'
   stop 'Program aborted'
endif

    ! Finish

    return

  end subroutine read_attr_r_sp_0_

  subroutine read_attr_r_sp_1_ (hg, item_name, data)

    type(hgroup_t), intent(inout)               :: hg
    character(*), intent(in)                    :: item_name
    real(SP), target, intent(out) :: data(:)

    contiguous :: data

    integer(HSIZE_T) :: item_shape(1)

    integer(HID_T)   :: mem_type_id
    integer          :: hdf_err
    integer(HID_T)   :: item_id
    type(C_PTR)      :: data_ptr

    ! Check the shape of the item agrees with that of data

    call get_attr_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 669 <core_hgroup:read_attr_r_sp_1_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

    ! Read the item

    mem_type_id = real_mem_type_(SP)

    data_ptr = C_LOC(data)

    call h5aopen_name_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 669 <core_hgroup:read_attr_r_sp_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aopen_name_f'
   stop 'Program aborted'
endif

    call h5aread_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 669 <core_hgroup:read_attr_r_sp_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aread_f'
   stop 'Program aborted'
endif

    call h5aclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 669 <core_hgroup:read_attr_r_sp_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aclose_f'
   stop 'Program aborted'
endif

    ! Finish

    return

  end subroutine read_attr_r_sp_1_

  subroutine read_attr_r_sp_2_ (hg, item_name, data)

    type(hgroup_t), intent(inout)               :: hg
    character(*), intent(in)                    :: item_name
    real(SP), target, intent(out) :: data(:,:)

    contiguous :: data

    integer(HSIZE_T) :: item_shape(2)

    integer(HID_T)   :: mem_type_id
    integer          :: hdf_err
    integer(HID_T)   :: item_id
    type(C_PTR)      :: data_ptr

    ! Check the shape of the item agrees with that of data

    call get_attr_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 670 <core_hgroup:read_attr_r_sp_2_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

    ! Read the item

    mem_type_id = real_mem_type_(SP)

    data_ptr = C_LOC(data)

    call h5aopen_name_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 670 <core_hgroup:read_attr_r_sp_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aopen_name_f'
   stop 'Program aborted'
endif

    call h5aread_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 670 <core_hgroup:read_attr_r_sp_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aread_f'
   stop 'Program aborted'
endif

    call h5aclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 670 <core_hgroup:read_attr_r_sp_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aclose_f'
   stop 'Program aborted'
endif

    ! Finish

    return

  end subroutine read_attr_r_sp_2_

  subroutine read_attr_r_sp_3_ (hg, item_name, data)

    type(hgroup_t), intent(inout)               :: hg
    character(*), intent(in)                    :: item_name
    real(SP), target, intent(out) :: data(:,:,:)

    contiguous :: data

    integer(HSIZE_T) :: item_shape(3)

    integer(HID_T)   :: mem_type_id
    integer          :: hdf_err
    integer(HID_T)   :: item_id
    type(C_PTR)      :: data_ptr

    ! Check the shape of the item agrees with that of data

    call get_attr_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 671 <core_hgroup:read_attr_r_sp_3_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

    ! Read the item

    mem_type_id = real_mem_type_(SP)

    data_ptr = C_LOC(data)

    call h5aopen_name_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 671 <core_hgroup:read_attr_r_sp_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aopen_name_f'
   stop 'Program aborted'
endif

    call h5aread_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 671 <core_hgroup:read_attr_r_sp_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aread_f'
   stop 'Program aborted'
endif

    call h5aclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 671 <core_hgroup:read_attr_r_sp_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aclose_f'
   stop 'Program aborted'
endif

    ! Finish

    return

  end subroutine read_attr_r_sp_3_

  subroutine read_attr_r_sp_4_ (hg, item_name, data)

    type(hgroup_t), intent(inout)               :: hg
    character(*), intent(in)                    :: item_name
    real(SP), target, intent(out) :: data(:,:,:,:)

    contiguous :: data

    integer(HSIZE_T) :: item_shape(4)

    integer(HID_T)   :: mem_type_id
    integer          :: hdf_err
    integer(HID_T)   :: item_id
    type(C_PTR)      :: data_ptr

    ! Check the shape of the item agrees with that of data

    call get_attr_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 672 <core_hgroup:read_attr_r_sp_4_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

    ! Read the item

    mem_type_id = real_mem_type_(SP)

    data_ptr = C_LOC(data)

    call h5aopen_name_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 672 <core_hgroup:read_attr_r_sp_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aopen_name_f'
   stop 'Program aborted'
endif

    call h5aread_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 672 <core_hgroup:read_attr_r_sp_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aread_f'
   stop 'Program aborted'
endif

    call h5aclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 672 <core_hgroup:read_attr_r_sp_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aclose_f'
   stop 'Program aborted'
endif

    ! Finish

    return

  end subroutine read_attr_r_sp_4_

  subroutine read_attr_r_dp_0_ (hg, item_name, data)

    type(hgroup_t), intent(inout)               :: hg
    character(*), intent(in)                    :: item_name
    real(DP), target, intent(out) :: data

    integer(HID_T)   :: mem_type_id
    integer          :: hdf_err
    integer(HID_T)   :: item_id
    type(C_PTR)      :: data_ptr

    ! Check the shape of the item agrees with that of data

    ! Read the item

    mem_type_id = real_mem_type_(DP)

    data_ptr = C_LOC(data)

    call h5aopen_name_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 674 <core_hgroup:read_attr_r_dp_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aopen_name_f'
   stop 'Program aborted'
endif

    call h5aread_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 674 <core_hgroup:read_attr_r_dp_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aread_f'
   stop 'Program aborted'
endif

    call h5aclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 674 <core_hgroup:read_attr_r_dp_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aclose_f'
   stop 'Program aborted'
endif

    ! Finish

    return

  end subroutine read_attr_r_dp_0_

  subroutine read_attr_r_dp_1_ (hg, item_name, data)

    type(hgroup_t), intent(inout)               :: hg
    character(*), intent(in)                    :: item_name
    real(DP), target, intent(out) :: data(:)

    contiguous :: data

    integer(HSIZE_T) :: item_shape(1)

    integer(HID_T)   :: mem_type_id
    integer          :: hdf_err
    integer(HID_T)   :: item_id
    type(C_PTR)      :: data_ptr

    ! Check the shape of the item agrees with that of data

    call get_attr_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 675 <core_hgroup:read_attr_r_dp_1_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

    ! Read the item

    mem_type_id = real_mem_type_(DP)

    data_ptr = C_LOC(data)

    call h5aopen_name_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 675 <core_hgroup:read_attr_r_dp_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aopen_name_f'
   stop 'Program aborted'
endif

    call h5aread_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 675 <core_hgroup:read_attr_r_dp_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aread_f'
   stop 'Program aborted'
endif

    call h5aclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 675 <core_hgroup:read_attr_r_dp_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aclose_f'
   stop 'Program aborted'
endif

    ! Finish

    return

  end subroutine read_attr_r_dp_1_

  subroutine read_attr_r_dp_2_ (hg, item_name, data)

    type(hgroup_t), intent(inout)               :: hg
    character(*), intent(in)                    :: item_name
    real(DP), target, intent(out) :: data(:,:)

    contiguous :: data

    integer(HSIZE_T) :: item_shape(2)

    integer(HID_T)   :: mem_type_id
    integer          :: hdf_err
    integer(HID_T)   :: item_id
    type(C_PTR)      :: data_ptr

    ! Check the shape of the item agrees with that of data

    call get_attr_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 676 <core_hgroup:read_attr_r_dp_2_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

    ! Read the item

    mem_type_id = real_mem_type_(DP)

    data_ptr = C_LOC(data)

    call h5aopen_name_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 676 <core_hgroup:read_attr_r_dp_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aopen_name_f'
   stop 'Program aborted'
endif

    call h5aread_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 676 <core_hgroup:read_attr_r_dp_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aread_f'
   stop 'Program aborted'
endif

    call h5aclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 676 <core_hgroup:read_attr_r_dp_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aclose_f'
   stop 'Program aborted'
endif

    ! Finish

    return

  end subroutine read_attr_r_dp_2_

  subroutine read_attr_r_dp_3_ (hg, item_name, data)

    type(hgroup_t), intent(inout)               :: hg
    character(*), intent(in)                    :: item_name
    real(DP), target, intent(out) :: data(:,:,:)

    contiguous :: data

    integer(HSIZE_T) :: item_shape(3)

    integer(HID_T)   :: mem_type_id
    integer          :: hdf_err
    integer(HID_T)   :: item_id
    type(C_PTR)      :: data_ptr

    ! Check the shape of the item agrees with that of data

    call get_attr_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 677 <core_hgroup:read_attr_r_dp_3_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

    ! Read the item

    mem_type_id = real_mem_type_(DP)

    data_ptr = C_LOC(data)

    call h5aopen_name_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 677 <core_hgroup:read_attr_r_dp_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aopen_name_f'
   stop 'Program aborted'
endif

    call h5aread_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 677 <core_hgroup:read_attr_r_dp_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aread_f'
   stop 'Program aborted'
endif

    call h5aclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 677 <core_hgroup:read_attr_r_dp_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aclose_f'
   stop 'Program aborted'
endif

    ! Finish

    return

  end subroutine read_attr_r_dp_3_

  subroutine read_attr_r_dp_4_ (hg, item_name, data)

    type(hgroup_t), intent(inout)               :: hg
    character(*), intent(in)                    :: item_name
    real(DP), target, intent(out) :: data(:,:,:,:)

    contiguous :: data

    integer(HSIZE_T) :: item_shape(4)

    integer(HID_T)   :: mem_type_id
    integer          :: hdf_err
    integer(HID_T)   :: item_id
    type(C_PTR)      :: data_ptr

    ! Check the shape of the item agrees with that of data

    call get_attr_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 678 <core_hgroup:read_attr_r_dp_4_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

    ! Read the item

    mem_type_id = real_mem_type_(DP)

    data_ptr = C_LOC(data)

    call h5aopen_name_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 678 <core_hgroup:read_attr_r_dp_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aopen_name_f'
   stop 'Program aborted'
endif

    call h5aread_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 678 <core_hgroup:read_attr_r_dp_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aread_f'
   stop 'Program aborted'
endif

    call h5aclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 678 <core_hgroup:read_attr_r_dp_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aclose_f'
   stop 'Program aborted'
endif

    ! Finish

    return

  end subroutine read_attr_r_dp_4_

  subroutine read_attr_c_sp_0_ (hg, item_name, data)

    type(hgroup_t), intent(inout)               :: hg
    character(*), intent(in)                    :: item_name
    complex(SP), target, intent(out) :: data

    integer(HID_T)   :: mem_type_id
    integer          :: hdf_err
    integer(HID_T)   :: item_id
    type(C_PTR)      :: data_ptr

    ! Check the shape of the item agrees with that of data

    ! Read the item

    mem_type_id = complex_mem_type_(SP)

    data_ptr = C_LOC(data)

    call h5aopen_name_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 680 <core_hgroup:read_attr_c_sp_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aopen_name_f'
   stop 'Program aborted'
endif

    call h5aread_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 680 <core_hgroup:read_attr_c_sp_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aread_f'
   stop 'Program aborted'
endif

    call h5aclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 680 <core_hgroup:read_attr_c_sp_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aclose_f'
   stop 'Program aborted'
endif

    ! Finish

    return

  end subroutine read_attr_c_sp_0_

  subroutine read_attr_c_sp_1_ (hg, item_name, data)

    type(hgroup_t), intent(inout)               :: hg
    character(*), intent(in)                    :: item_name
    complex(SP), target, intent(out) :: data(:)

    contiguous :: data

    integer(HSIZE_T) :: item_shape(1)

    integer(HID_T)   :: mem_type_id
    integer          :: hdf_err
    integer(HID_T)   :: item_id
    type(C_PTR)      :: data_ptr

    ! Check the shape of the item agrees with that of data

    call get_attr_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 681 <core_hgroup:read_attr_c_sp_1_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

    ! Read the item

    mem_type_id = complex_mem_type_(SP)

    data_ptr = C_LOC(data)

    call h5aopen_name_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 681 <core_hgroup:read_attr_c_sp_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aopen_name_f'
   stop 'Program aborted'
endif

    call h5aread_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 681 <core_hgroup:read_attr_c_sp_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aread_f'
   stop 'Program aborted'
endif

    call h5aclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 681 <core_hgroup:read_attr_c_sp_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aclose_f'
   stop 'Program aborted'
endif

    ! Finish

    return

  end subroutine read_attr_c_sp_1_

  subroutine read_attr_c_sp_2_ (hg, item_name, data)

    type(hgroup_t), intent(inout)               :: hg
    character(*), intent(in)                    :: item_name
    complex(SP), target, intent(out) :: data(:,:)

    contiguous :: data

    integer(HSIZE_T) :: item_shape(2)

    integer(HID_T)   :: mem_type_id
    integer          :: hdf_err
    integer(HID_T)   :: item_id
    type(C_PTR)      :: data_ptr

    ! Check the shape of the item agrees with that of data

    call get_attr_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 682 <core_hgroup:read_attr_c_sp_2_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

    ! Read the item

    mem_type_id = complex_mem_type_(SP)

    data_ptr = C_LOC(data)

    call h5aopen_name_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 682 <core_hgroup:read_attr_c_sp_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aopen_name_f'
   stop 'Program aborted'
endif

    call h5aread_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 682 <core_hgroup:read_attr_c_sp_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aread_f'
   stop 'Program aborted'
endif

    call h5aclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 682 <core_hgroup:read_attr_c_sp_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aclose_f'
   stop 'Program aborted'
endif

    ! Finish

    return

  end subroutine read_attr_c_sp_2_

  subroutine read_attr_c_sp_3_ (hg, item_name, data)

    type(hgroup_t), intent(inout)               :: hg
    character(*), intent(in)                    :: item_name
    complex(SP), target, intent(out) :: data(:,:,:)

    contiguous :: data

    integer(HSIZE_T) :: item_shape(3)

    integer(HID_T)   :: mem_type_id
    integer          :: hdf_err
    integer(HID_T)   :: item_id
    type(C_PTR)      :: data_ptr

    ! Check the shape of the item agrees with that of data

    call get_attr_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 683 <core_hgroup:read_attr_c_sp_3_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

    ! Read the item

    mem_type_id = complex_mem_type_(SP)

    data_ptr = C_LOC(data)

    call h5aopen_name_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 683 <core_hgroup:read_attr_c_sp_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aopen_name_f'
   stop 'Program aborted'
endif

    call h5aread_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 683 <core_hgroup:read_attr_c_sp_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aread_f'
   stop 'Program aborted'
endif

    call h5aclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 683 <core_hgroup:read_attr_c_sp_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aclose_f'
   stop 'Program aborted'
endif

    ! Finish

    return

  end subroutine read_attr_c_sp_3_

  subroutine read_attr_c_sp_4_ (hg, item_name, data)

    type(hgroup_t), intent(inout)               :: hg
    character(*), intent(in)                    :: item_name
    complex(SP), target, intent(out) :: data(:,:,:,:)

    contiguous :: data

    integer(HSIZE_T) :: item_shape(4)

    integer(HID_T)   :: mem_type_id
    integer          :: hdf_err
    integer(HID_T)   :: item_id
    type(C_PTR)      :: data_ptr

    ! Check the shape of the item agrees with that of data

    call get_attr_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 684 <core_hgroup:read_attr_c_sp_4_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

    ! Read the item

    mem_type_id = complex_mem_type_(SP)

    data_ptr = C_LOC(data)

    call h5aopen_name_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 684 <core_hgroup:read_attr_c_sp_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aopen_name_f'
   stop 'Program aborted'
endif

    call h5aread_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 684 <core_hgroup:read_attr_c_sp_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aread_f'
   stop 'Program aborted'
endif

    call h5aclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 684 <core_hgroup:read_attr_c_sp_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aclose_f'
   stop 'Program aborted'
endif

    ! Finish

    return

  end subroutine read_attr_c_sp_4_

  subroutine read_attr_c_dp_0_ (hg, item_name, data)

    type(hgroup_t), intent(inout)               :: hg
    character(*), intent(in)                    :: item_name
    complex(DP), target, intent(out) :: data

    integer(HID_T)   :: mem_type_id
    integer          :: hdf_err
    integer(HID_T)   :: item_id
    type(C_PTR)      :: data_ptr

    ! Check the shape of the item agrees with that of data

    ! Read the item

    mem_type_id = complex_mem_type_(DP)

    data_ptr = C_LOC(data)

    call h5aopen_name_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 686 <core_hgroup:read_attr_c_dp_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aopen_name_f'
   stop 'Program aborted'
endif

    call h5aread_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 686 <core_hgroup:read_attr_c_dp_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aread_f'
   stop 'Program aborted'
endif

    call h5aclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 686 <core_hgroup:read_attr_c_dp_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aclose_f'
   stop 'Program aborted'
endif

    ! Finish

    return

  end subroutine read_attr_c_dp_0_

  subroutine read_attr_c_dp_1_ (hg, item_name, data)

    type(hgroup_t), intent(inout)               :: hg
    character(*), intent(in)                    :: item_name
    complex(DP), target, intent(out) :: data(:)

    contiguous :: data

    integer(HSIZE_T) :: item_shape(1)

    integer(HID_T)   :: mem_type_id
    integer          :: hdf_err
    integer(HID_T)   :: item_id
    type(C_PTR)      :: data_ptr

    ! Check the shape of the item agrees with that of data

    call get_attr_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 687 <core_hgroup:read_attr_c_dp_1_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

    ! Read the item

    mem_type_id = complex_mem_type_(DP)

    data_ptr = C_LOC(data)

    call h5aopen_name_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 687 <core_hgroup:read_attr_c_dp_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aopen_name_f'
   stop 'Program aborted'
endif

    call h5aread_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 687 <core_hgroup:read_attr_c_dp_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aread_f'
   stop 'Program aborted'
endif

    call h5aclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 687 <core_hgroup:read_attr_c_dp_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aclose_f'
   stop 'Program aborted'
endif

    ! Finish

    return

  end subroutine read_attr_c_dp_1_

  subroutine read_attr_c_dp_2_ (hg, item_name, data)

    type(hgroup_t), intent(inout)               :: hg
    character(*), intent(in)                    :: item_name
    complex(DP), target, intent(out) :: data(:,:)

    contiguous :: data

    integer(HSIZE_T) :: item_shape(2)

    integer(HID_T)   :: mem_type_id
    integer          :: hdf_err
    integer(HID_T)   :: item_id
    type(C_PTR)      :: data_ptr

    ! Check the shape of the item agrees with that of data

    call get_attr_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 688 <core_hgroup:read_attr_c_dp_2_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

    ! Read the item

    mem_type_id = complex_mem_type_(DP)

    data_ptr = C_LOC(data)

    call h5aopen_name_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 688 <core_hgroup:read_attr_c_dp_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aopen_name_f'
   stop 'Program aborted'
endif

    call h5aread_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 688 <core_hgroup:read_attr_c_dp_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aread_f'
   stop 'Program aborted'
endif

    call h5aclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 688 <core_hgroup:read_attr_c_dp_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aclose_f'
   stop 'Program aborted'
endif

    ! Finish

    return

  end subroutine read_attr_c_dp_2_

  subroutine read_attr_c_dp_3_ (hg, item_name, data)

    type(hgroup_t), intent(inout)               :: hg
    character(*), intent(in)                    :: item_name
    complex(DP), target, intent(out) :: data(:,:,:)

    contiguous :: data

    integer(HSIZE_T) :: item_shape(3)

    integer(HID_T)   :: mem_type_id
    integer          :: hdf_err
    integer(HID_T)   :: item_id
    type(C_PTR)      :: data_ptr

    ! Check the shape of the item agrees with that of data

    call get_attr_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 689 <core_hgroup:read_attr_c_dp_3_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

    ! Read the item

    mem_type_id = complex_mem_type_(DP)

    data_ptr = C_LOC(data)

    call h5aopen_name_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 689 <core_hgroup:read_attr_c_dp_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aopen_name_f'
   stop 'Program aborted'
endif

    call h5aread_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 689 <core_hgroup:read_attr_c_dp_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aread_f'
   stop 'Program aborted'
endif

    call h5aclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 689 <core_hgroup:read_attr_c_dp_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aclose_f'
   stop 'Program aborted'
endif

    ! Finish

    return

  end subroutine read_attr_c_dp_3_

  subroutine read_attr_c_dp_4_ (hg, item_name, data)

    type(hgroup_t), intent(inout)               :: hg
    character(*), intent(in)                    :: item_name
    complex(DP), target, intent(out) :: data(:,:,:,:)

    contiguous :: data

    integer(HSIZE_T) :: item_shape(4)

    integer(HID_T)   :: mem_type_id
    integer          :: hdf_err
    integer(HID_T)   :: item_id
    type(C_PTR)      :: data_ptr

    ! Check the shape of the item agrees with that of data

    call get_attr_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 690 <core_hgroup:read_attr_c_dp_4_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

    ! Read the item

    mem_type_id = complex_mem_type_(DP)

    data_ptr = C_LOC(data)

    call h5aopen_name_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 690 <core_hgroup:read_attr_c_dp_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aopen_name_f'
   stop 'Program aborted'
endif

    call h5aread_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 690 <core_hgroup:read_attr_c_dp_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aread_f'
   stop 'Program aborted'
endif

    call h5aclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 690 <core_hgroup:read_attr_c_dp_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aclose_f'
   stop 'Program aborted'
endif

    ! Finish

    return

  end subroutine read_attr_c_dp_4_

  subroutine read_dset_i_i4_0_ (hg, item_name, data)

    type(hgroup_t), intent(inout)               :: hg
    character(*), intent(in)                    :: item_name
    integer(I4), target, intent(out) :: data

    integer(HID_T)   :: mem_type_id
    integer          :: hdf_err
    integer(HID_T)   :: item_id
    type(C_PTR)      :: data_ptr

    ! Check the shape of the item agrees with that of data

    ! Read the item

    mem_type_id = integer_mem_type_(I4)

    data_ptr = C_LOC(data)

    call h5dopen_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 692 <core_hgroup:read_dset_i_i4_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dopen_f'
   stop 'Program aborted'
endif

    call h5dread_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 692 <core_hgroup:read_dset_i_i4_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dread_f'
   stop 'Program aborted'
endif

    call h5dclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 692 <core_hgroup:read_dset_i_i4_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dclose_f'
   stop 'Program aborted'
endif

    ! Finish

    return

  end subroutine read_dset_i_i4_0_

  subroutine read_dset_i_i4_1_ (hg, item_name, data)

    type(hgroup_t), intent(inout)               :: hg
    character(*), intent(in)                    :: item_name
    integer(I4), target, intent(out) :: data(:)

    contiguous :: data

    integer(HSIZE_T) :: item_shape(1)

    integer(HID_T)   :: mem_type_id
    integer          :: hdf_err
    integer(HID_T)   :: item_id
    type(C_PTR)      :: data_ptr

    ! Check the shape of the item agrees with that of data

    call get_dset_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 693 <core_hgroup:read_dset_i_i4_1_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

    ! Read the item

    mem_type_id = integer_mem_type_(I4)

    data_ptr = C_LOC(data)

    call h5dopen_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 693 <core_hgroup:read_dset_i_i4_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dopen_f'
   stop 'Program aborted'
endif

    call h5dread_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 693 <core_hgroup:read_dset_i_i4_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dread_f'
   stop 'Program aborted'
endif

    call h5dclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 693 <core_hgroup:read_dset_i_i4_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dclose_f'
   stop 'Program aborted'
endif

    ! Finish

    return

  end subroutine read_dset_i_i4_1_

  subroutine read_dset_i_i4_2_ (hg, item_name, data)

    type(hgroup_t), intent(inout)               :: hg
    character(*), intent(in)                    :: item_name
    integer(I4), target, intent(out) :: data(:,:)

    contiguous :: data

    integer(HSIZE_T) :: item_shape(2)

    integer(HID_T)   :: mem_type_id
    integer          :: hdf_err
    integer(HID_T)   :: item_id
    type(C_PTR)      :: data_ptr

    ! Check the shape of the item agrees with that of data

    call get_dset_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 694 <core_hgroup:read_dset_i_i4_2_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

    ! Read the item

    mem_type_id = integer_mem_type_(I4)

    data_ptr = C_LOC(data)

    call h5dopen_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 694 <core_hgroup:read_dset_i_i4_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dopen_f'
   stop 'Program aborted'
endif

    call h5dread_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 694 <core_hgroup:read_dset_i_i4_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dread_f'
   stop 'Program aborted'
endif

    call h5dclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 694 <core_hgroup:read_dset_i_i4_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dclose_f'
   stop 'Program aborted'
endif

    ! Finish

    return

  end subroutine read_dset_i_i4_2_

  subroutine read_dset_i_i4_3_ (hg, item_name, data)

    type(hgroup_t), intent(inout)               :: hg
    character(*), intent(in)                    :: item_name
    integer(I4), target, intent(out) :: data(:,:,:)

    contiguous :: data

    integer(HSIZE_T) :: item_shape(3)

    integer(HID_T)   :: mem_type_id
    integer          :: hdf_err
    integer(HID_T)   :: item_id
    type(C_PTR)      :: data_ptr

    ! Check the shape of the item agrees with that of data

    call get_dset_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 695 <core_hgroup:read_dset_i_i4_3_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

    ! Read the item

    mem_type_id = integer_mem_type_(I4)

    data_ptr = C_LOC(data)

    call h5dopen_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 695 <core_hgroup:read_dset_i_i4_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dopen_f'
   stop 'Program aborted'
endif

    call h5dread_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 695 <core_hgroup:read_dset_i_i4_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dread_f'
   stop 'Program aborted'
endif

    call h5dclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 695 <core_hgroup:read_dset_i_i4_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dclose_f'
   stop 'Program aborted'
endif

    ! Finish

    return

  end subroutine read_dset_i_i4_3_

  subroutine read_dset_i_i4_4_ (hg, item_name, data)

    type(hgroup_t), intent(inout)               :: hg
    character(*), intent(in)                    :: item_name
    integer(I4), target, intent(out) :: data(:,:,:,:)

    contiguous :: data

    integer(HSIZE_T) :: item_shape(4)

    integer(HID_T)   :: mem_type_id
    integer          :: hdf_err
    integer(HID_T)   :: item_id
    type(C_PTR)      :: data_ptr

    ! Check the shape of the item agrees with that of data

    call get_dset_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 696 <core_hgroup:read_dset_i_i4_4_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

    ! Read the item

    mem_type_id = integer_mem_type_(I4)

    data_ptr = C_LOC(data)

    call h5dopen_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 696 <core_hgroup:read_dset_i_i4_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dopen_f'
   stop 'Program aborted'
endif

    call h5dread_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 696 <core_hgroup:read_dset_i_i4_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dread_f'
   stop 'Program aborted'
endif

    call h5dclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 696 <core_hgroup:read_dset_i_i4_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dclose_f'
   stop 'Program aborted'
endif

    ! Finish

    return

  end subroutine read_dset_i_i4_4_

  subroutine read_dset_i_i8_0_ (hg, item_name, data)

    type(hgroup_t), intent(inout)               :: hg
    character(*), intent(in)                    :: item_name
    integer(I8), target, intent(out) :: data

    integer(HID_T)   :: mem_type_id
    integer          :: hdf_err
    integer(HID_T)   :: item_id
    type(C_PTR)      :: data_ptr

    ! Check the shape of the item agrees with that of data

    ! Read the item

    mem_type_id = integer_mem_type_(I8)

    data_ptr = C_LOC(data)

    call h5dopen_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 698 <core_hgroup:read_dset_i_i8_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dopen_f'
   stop 'Program aborted'
endif

    call h5dread_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 698 <core_hgroup:read_dset_i_i8_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dread_f'
   stop 'Program aborted'
endif

    call h5dclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 698 <core_hgroup:read_dset_i_i8_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dclose_f'
   stop 'Program aborted'
endif

    ! Finish

    return

  end subroutine read_dset_i_i8_0_

  subroutine read_dset_i_i8_1_ (hg, item_name, data)

    type(hgroup_t), intent(inout)               :: hg
    character(*), intent(in)                    :: item_name
    integer(I8), target, intent(out) :: data(:)

    contiguous :: data

    integer(HSIZE_T) :: item_shape(1)

    integer(HID_T)   :: mem_type_id
    integer          :: hdf_err
    integer(HID_T)   :: item_id
    type(C_PTR)      :: data_ptr

    ! Check the shape of the item agrees with that of data

    call get_dset_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 699 <core_hgroup:read_dset_i_i8_1_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

    ! Read the item

    mem_type_id = integer_mem_type_(I8)

    data_ptr = C_LOC(data)

    call h5dopen_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 699 <core_hgroup:read_dset_i_i8_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dopen_f'
   stop 'Program aborted'
endif

    call h5dread_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 699 <core_hgroup:read_dset_i_i8_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dread_f'
   stop 'Program aborted'
endif

    call h5dclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 699 <core_hgroup:read_dset_i_i8_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dclose_f'
   stop 'Program aborted'
endif

    ! Finish

    return

  end subroutine read_dset_i_i8_1_

  subroutine read_dset_i_i8_2_ (hg, item_name, data)

    type(hgroup_t), intent(inout)               :: hg
    character(*), intent(in)                    :: item_name
    integer(I8), target, intent(out) :: data(:,:)

    contiguous :: data

    integer(HSIZE_T) :: item_shape(2)

    integer(HID_T)   :: mem_type_id
    integer          :: hdf_err
    integer(HID_T)   :: item_id
    type(C_PTR)      :: data_ptr

    ! Check the shape of the item agrees with that of data

    call get_dset_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 700 <core_hgroup:read_dset_i_i8_2_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

    ! Read the item

    mem_type_id = integer_mem_type_(I8)

    data_ptr = C_LOC(data)

    call h5dopen_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 700 <core_hgroup:read_dset_i_i8_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dopen_f'
   stop 'Program aborted'
endif

    call h5dread_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 700 <core_hgroup:read_dset_i_i8_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dread_f'
   stop 'Program aborted'
endif

    call h5dclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 700 <core_hgroup:read_dset_i_i8_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dclose_f'
   stop 'Program aborted'
endif

    ! Finish

    return

  end subroutine read_dset_i_i8_2_

  subroutine read_dset_i_i8_3_ (hg, item_name, data)

    type(hgroup_t), intent(inout)               :: hg
    character(*), intent(in)                    :: item_name
    integer(I8), target, intent(out) :: data(:,:,:)

    contiguous :: data

    integer(HSIZE_T) :: item_shape(3)

    integer(HID_T)   :: mem_type_id
    integer          :: hdf_err
    integer(HID_T)   :: item_id
    type(C_PTR)      :: data_ptr

    ! Check the shape of the item agrees with that of data

    call get_dset_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 701 <core_hgroup:read_dset_i_i8_3_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

    ! Read the item

    mem_type_id = integer_mem_type_(I8)

    data_ptr = C_LOC(data)

    call h5dopen_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 701 <core_hgroup:read_dset_i_i8_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dopen_f'
   stop 'Program aborted'
endif

    call h5dread_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 701 <core_hgroup:read_dset_i_i8_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dread_f'
   stop 'Program aborted'
endif

    call h5dclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 701 <core_hgroup:read_dset_i_i8_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dclose_f'
   stop 'Program aborted'
endif

    ! Finish

    return

  end subroutine read_dset_i_i8_3_

  subroutine read_dset_i_i8_4_ (hg, item_name, data)

    type(hgroup_t), intent(inout)               :: hg
    character(*), intent(in)                    :: item_name
    integer(I8), target, intent(out) :: data(:,:,:,:)

    contiguous :: data

    integer(HSIZE_T) :: item_shape(4)

    integer(HID_T)   :: mem_type_id
    integer          :: hdf_err
    integer(HID_T)   :: item_id
    type(C_PTR)      :: data_ptr

    ! Check the shape of the item agrees with that of data

    call get_dset_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 702 <core_hgroup:read_dset_i_i8_4_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

    ! Read the item

    mem_type_id = integer_mem_type_(I8)

    data_ptr = C_LOC(data)

    call h5dopen_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 702 <core_hgroup:read_dset_i_i8_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dopen_f'
   stop 'Program aborted'
endif

    call h5dread_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 702 <core_hgroup:read_dset_i_i8_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dread_f'
   stop 'Program aborted'
endif

    call h5dclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 702 <core_hgroup:read_dset_i_i8_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dclose_f'
   stop 'Program aborted'
endif

    ! Finish

    return

  end subroutine read_dset_i_i8_4_

  subroutine read_dset_r_sp_0_ (hg, item_name, data)

    type(hgroup_t), intent(inout)               :: hg
    character(*), intent(in)                    :: item_name
    real(SP), target, intent(out) :: data

    integer(HID_T)   :: mem_type_id
    integer          :: hdf_err
    integer(HID_T)   :: item_id
    type(C_PTR)      :: data_ptr

    ! Check the shape of the item agrees with that of data

    ! Read the item

    mem_type_id = real_mem_type_(SP)

    data_ptr = C_LOC(data)

    call h5dopen_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 704 <core_hgroup:read_dset_r_sp_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dopen_f'
   stop 'Program aborted'
endif

    call h5dread_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 704 <core_hgroup:read_dset_r_sp_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dread_f'
   stop 'Program aborted'
endif

    call h5dclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 704 <core_hgroup:read_dset_r_sp_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dclose_f'
   stop 'Program aborted'
endif

    ! Finish

    return

  end subroutine read_dset_r_sp_0_

  subroutine read_dset_r_sp_1_ (hg, item_name, data)

    type(hgroup_t), intent(inout)               :: hg
    character(*), intent(in)                    :: item_name
    real(SP), target, intent(out) :: data(:)

    contiguous :: data

    integer(HSIZE_T) :: item_shape(1)

    integer(HID_T)   :: mem_type_id
    integer          :: hdf_err
    integer(HID_T)   :: item_id
    type(C_PTR)      :: data_ptr

    ! Check the shape of the item agrees with that of data

    call get_dset_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 705 <core_hgroup:read_dset_r_sp_1_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

    ! Read the item

    mem_type_id = real_mem_type_(SP)

    data_ptr = C_LOC(data)

    call h5dopen_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 705 <core_hgroup:read_dset_r_sp_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dopen_f'
   stop 'Program aborted'
endif

    call h5dread_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 705 <core_hgroup:read_dset_r_sp_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dread_f'
   stop 'Program aborted'
endif

    call h5dclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 705 <core_hgroup:read_dset_r_sp_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dclose_f'
   stop 'Program aborted'
endif

    ! Finish

    return

  end subroutine read_dset_r_sp_1_

  subroutine read_dset_r_sp_2_ (hg, item_name, data)

    type(hgroup_t), intent(inout)               :: hg
    character(*), intent(in)                    :: item_name
    real(SP), target, intent(out) :: data(:,:)

    contiguous :: data

    integer(HSIZE_T) :: item_shape(2)

    integer(HID_T)   :: mem_type_id
    integer          :: hdf_err
    integer(HID_T)   :: item_id
    type(C_PTR)      :: data_ptr

    ! Check the shape of the item agrees with that of data

    call get_dset_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 706 <core_hgroup:read_dset_r_sp_2_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

    ! Read the item

    mem_type_id = real_mem_type_(SP)

    data_ptr = C_LOC(data)

    call h5dopen_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 706 <core_hgroup:read_dset_r_sp_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dopen_f'
   stop 'Program aborted'
endif

    call h5dread_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 706 <core_hgroup:read_dset_r_sp_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dread_f'
   stop 'Program aborted'
endif

    call h5dclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 706 <core_hgroup:read_dset_r_sp_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dclose_f'
   stop 'Program aborted'
endif

    ! Finish

    return

  end subroutine read_dset_r_sp_2_

  subroutine read_dset_r_sp_3_ (hg, item_name, data)

    type(hgroup_t), intent(inout)               :: hg
    character(*), intent(in)                    :: item_name
    real(SP), target, intent(out) :: data(:,:,:)

    contiguous :: data

    integer(HSIZE_T) :: item_shape(3)

    integer(HID_T)   :: mem_type_id
    integer          :: hdf_err
    integer(HID_T)   :: item_id
    type(C_PTR)      :: data_ptr

    ! Check the shape of the item agrees with that of data

    call get_dset_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 707 <core_hgroup:read_dset_r_sp_3_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

    ! Read the item

    mem_type_id = real_mem_type_(SP)

    data_ptr = C_LOC(data)

    call h5dopen_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 707 <core_hgroup:read_dset_r_sp_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dopen_f'
   stop 'Program aborted'
endif

    call h5dread_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 707 <core_hgroup:read_dset_r_sp_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dread_f'
   stop 'Program aborted'
endif

    call h5dclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 707 <core_hgroup:read_dset_r_sp_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dclose_f'
   stop 'Program aborted'
endif

    ! Finish

    return

  end subroutine read_dset_r_sp_3_

  subroutine read_dset_r_sp_4_ (hg, item_name, data)

    type(hgroup_t), intent(inout)               :: hg
    character(*), intent(in)                    :: item_name
    real(SP), target, intent(out) :: data(:,:,:,:)

    contiguous :: data

    integer(HSIZE_T) :: item_shape(4)

    integer(HID_T)   :: mem_type_id
    integer          :: hdf_err
    integer(HID_T)   :: item_id
    type(C_PTR)      :: data_ptr

    ! Check the shape of the item agrees with that of data

    call get_dset_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 708 <core_hgroup:read_dset_r_sp_4_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

    ! Read the item

    mem_type_id = real_mem_type_(SP)

    data_ptr = C_LOC(data)

    call h5dopen_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 708 <core_hgroup:read_dset_r_sp_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dopen_f'
   stop 'Program aborted'
endif

    call h5dread_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 708 <core_hgroup:read_dset_r_sp_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dread_f'
   stop 'Program aborted'
endif

    call h5dclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 708 <core_hgroup:read_dset_r_sp_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dclose_f'
   stop 'Program aborted'
endif

    ! Finish

    return

  end subroutine read_dset_r_sp_4_

  subroutine read_dset_r_dp_0_ (hg, item_name, data)

    type(hgroup_t), intent(inout)               :: hg
    character(*), intent(in)                    :: item_name
    real(DP), target, intent(out) :: data

    integer(HID_T)   :: mem_type_id
    integer          :: hdf_err
    integer(HID_T)   :: item_id
    type(C_PTR)      :: data_ptr

    ! Check the shape of the item agrees with that of data

    ! Read the item

    mem_type_id = real_mem_type_(DP)

    data_ptr = C_LOC(data)

    call h5dopen_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 710 <core_hgroup:read_dset_r_dp_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dopen_f'
   stop 'Program aborted'
endif

    call h5dread_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 710 <core_hgroup:read_dset_r_dp_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dread_f'
   stop 'Program aborted'
endif

    call h5dclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 710 <core_hgroup:read_dset_r_dp_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dclose_f'
   stop 'Program aborted'
endif

    ! Finish

    return

  end subroutine read_dset_r_dp_0_

  subroutine read_dset_r_dp_1_ (hg, item_name, data)

    type(hgroup_t), intent(inout)               :: hg
    character(*), intent(in)                    :: item_name
    real(DP), target, intent(out) :: data(:)

    contiguous :: data

    integer(HSIZE_T) :: item_shape(1)

    integer(HID_T)   :: mem_type_id
    integer          :: hdf_err
    integer(HID_T)   :: item_id
    type(C_PTR)      :: data_ptr

    ! Check the shape of the item agrees with that of data

    call get_dset_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 711 <core_hgroup:read_dset_r_dp_1_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

    ! Read the item

    mem_type_id = real_mem_type_(DP)

    data_ptr = C_LOC(data)

    call h5dopen_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 711 <core_hgroup:read_dset_r_dp_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dopen_f'
   stop 'Program aborted'
endif

    call h5dread_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 711 <core_hgroup:read_dset_r_dp_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dread_f'
   stop 'Program aborted'
endif

    call h5dclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 711 <core_hgroup:read_dset_r_dp_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dclose_f'
   stop 'Program aborted'
endif

    ! Finish

    return

  end subroutine read_dset_r_dp_1_

  subroutine read_dset_r_dp_2_ (hg, item_name, data)

    type(hgroup_t), intent(inout)               :: hg
    character(*), intent(in)                    :: item_name
    real(DP), target, intent(out) :: data(:,:)

    contiguous :: data

    integer(HSIZE_T) :: item_shape(2)

    integer(HID_T)   :: mem_type_id
    integer          :: hdf_err
    integer(HID_T)   :: item_id
    type(C_PTR)      :: data_ptr

    ! Check the shape of the item agrees with that of data

    call get_dset_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 712 <core_hgroup:read_dset_r_dp_2_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

    ! Read the item

    mem_type_id = real_mem_type_(DP)

    data_ptr = C_LOC(data)

    call h5dopen_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 712 <core_hgroup:read_dset_r_dp_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dopen_f'
   stop 'Program aborted'
endif

    call h5dread_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 712 <core_hgroup:read_dset_r_dp_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dread_f'
   stop 'Program aborted'
endif

    call h5dclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 712 <core_hgroup:read_dset_r_dp_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dclose_f'
   stop 'Program aborted'
endif

    ! Finish

    return

  end subroutine read_dset_r_dp_2_

  subroutine read_dset_r_dp_3_ (hg, item_name, data)

    type(hgroup_t), intent(inout)               :: hg
    character(*), intent(in)                    :: item_name
    real(DP), target, intent(out) :: data(:,:,:)

    contiguous :: data

    integer(HSIZE_T) :: item_shape(3)

    integer(HID_T)   :: mem_type_id
    integer          :: hdf_err
    integer(HID_T)   :: item_id
    type(C_PTR)      :: data_ptr

    ! Check the shape of the item agrees with that of data

    call get_dset_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 713 <core_hgroup:read_dset_r_dp_3_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

    ! Read the item

    mem_type_id = real_mem_type_(DP)

    data_ptr = C_LOC(data)

    call h5dopen_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 713 <core_hgroup:read_dset_r_dp_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dopen_f'
   stop 'Program aborted'
endif

    call h5dread_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 713 <core_hgroup:read_dset_r_dp_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dread_f'
   stop 'Program aborted'
endif

    call h5dclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 713 <core_hgroup:read_dset_r_dp_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dclose_f'
   stop 'Program aborted'
endif

    ! Finish

    return

  end subroutine read_dset_r_dp_3_

  subroutine read_dset_r_dp_4_ (hg, item_name, data)

    type(hgroup_t), intent(inout)               :: hg
    character(*), intent(in)                    :: item_name
    real(DP), target, intent(out) :: data(:,:,:,:)

    contiguous :: data

    integer(HSIZE_T) :: item_shape(4)

    integer(HID_T)   :: mem_type_id
    integer          :: hdf_err
    integer(HID_T)   :: item_id
    type(C_PTR)      :: data_ptr

    ! Check the shape of the item agrees with that of data

    call get_dset_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 714 <core_hgroup:read_dset_r_dp_4_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

    ! Read the item

    mem_type_id = real_mem_type_(DP)

    data_ptr = C_LOC(data)

    call h5dopen_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 714 <core_hgroup:read_dset_r_dp_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dopen_f'
   stop 'Program aborted'
endif

    call h5dread_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 714 <core_hgroup:read_dset_r_dp_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dread_f'
   stop 'Program aborted'
endif

    call h5dclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 714 <core_hgroup:read_dset_r_dp_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dclose_f'
   stop 'Program aborted'
endif

    ! Finish

    return

  end subroutine read_dset_r_dp_4_

  subroutine read_dset_c_sp_0_ (hg, item_name, data)

    type(hgroup_t), intent(inout)               :: hg
    character(*), intent(in)                    :: item_name
    complex(SP), target, intent(out) :: data

    integer(HID_T)   :: mem_type_id
    integer          :: hdf_err
    integer(HID_T)   :: item_id
    type(C_PTR)      :: data_ptr

    ! Check the shape of the item agrees with that of data

    ! Read the item

    mem_type_id = complex_mem_type_(SP)

    data_ptr = C_LOC(data)

    call h5dopen_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 716 <core_hgroup:read_dset_c_sp_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dopen_f'
   stop 'Program aborted'
endif

    call h5dread_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 716 <core_hgroup:read_dset_c_sp_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dread_f'
   stop 'Program aborted'
endif

    call h5dclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 716 <core_hgroup:read_dset_c_sp_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dclose_f'
   stop 'Program aborted'
endif

    ! Finish

    return

  end subroutine read_dset_c_sp_0_

  subroutine read_dset_c_sp_1_ (hg, item_name, data)

    type(hgroup_t), intent(inout)               :: hg
    character(*), intent(in)                    :: item_name
    complex(SP), target, intent(out) :: data(:)

    contiguous :: data

    integer(HSIZE_T) :: item_shape(1)

    integer(HID_T)   :: mem_type_id
    integer          :: hdf_err
    integer(HID_T)   :: item_id
    type(C_PTR)      :: data_ptr

    ! Check the shape of the item agrees with that of data

    call get_dset_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 717 <core_hgroup:read_dset_c_sp_1_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

    ! Read the item

    mem_type_id = complex_mem_type_(SP)

    data_ptr = C_LOC(data)

    call h5dopen_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 717 <core_hgroup:read_dset_c_sp_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dopen_f'
   stop 'Program aborted'
endif

    call h5dread_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 717 <core_hgroup:read_dset_c_sp_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dread_f'
   stop 'Program aborted'
endif

    call h5dclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 717 <core_hgroup:read_dset_c_sp_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dclose_f'
   stop 'Program aborted'
endif

    ! Finish

    return

  end subroutine read_dset_c_sp_1_

  subroutine read_dset_c_sp_2_ (hg, item_name, data)

    type(hgroup_t), intent(inout)               :: hg
    character(*), intent(in)                    :: item_name
    complex(SP), target, intent(out) :: data(:,:)

    contiguous :: data

    integer(HSIZE_T) :: item_shape(2)

    integer(HID_T)   :: mem_type_id
    integer          :: hdf_err
    integer(HID_T)   :: item_id
    type(C_PTR)      :: data_ptr

    ! Check the shape of the item agrees with that of data

    call get_dset_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 718 <core_hgroup:read_dset_c_sp_2_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

    ! Read the item

    mem_type_id = complex_mem_type_(SP)

    data_ptr = C_LOC(data)

    call h5dopen_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 718 <core_hgroup:read_dset_c_sp_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dopen_f'
   stop 'Program aborted'
endif

    call h5dread_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 718 <core_hgroup:read_dset_c_sp_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dread_f'
   stop 'Program aborted'
endif

    call h5dclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 718 <core_hgroup:read_dset_c_sp_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dclose_f'
   stop 'Program aborted'
endif

    ! Finish

    return

  end subroutine read_dset_c_sp_2_

  subroutine read_dset_c_sp_3_ (hg, item_name, data)

    type(hgroup_t), intent(inout)               :: hg
    character(*), intent(in)                    :: item_name
    complex(SP), target, intent(out) :: data(:,:,:)

    contiguous :: data

    integer(HSIZE_T) :: item_shape(3)

    integer(HID_T)   :: mem_type_id
    integer          :: hdf_err
    integer(HID_T)   :: item_id
    type(C_PTR)      :: data_ptr

    ! Check the shape of the item agrees with that of data

    call get_dset_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 719 <core_hgroup:read_dset_c_sp_3_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

    ! Read the item

    mem_type_id = complex_mem_type_(SP)

    data_ptr = C_LOC(data)

    call h5dopen_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 719 <core_hgroup:read_dset_c_sp_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dopen_f'
   stop 'Program aborted'
endif

    call h5dread_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 719 <core_hgroup:read_dset_c_sp_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dread_f'
   stop 'Program aborted'
endif

    call h5dclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 719 <core_hgroup:read_dset_c_sp_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dclose_f'
   stop 'Program aborted'
endif

    ! Finish

    return

  end subroutine read_dset_c_sp_3_

  subroutine read_dset_c_sp_4_ (hg, item_name, data)

    type(hgroup_t), intent(inout)               :: hg
    character(*), intent(in)                    :: item_name
    complex(SP), target, intent(out) :: data(:,:,:,:)

    contiguous :: data

    integer(HSIZE_T) :: item_shape(4)

    integer(HID_T)   :: mem_type_id
    integer          :: hdf_err
    integer(HID_T)   :: item_id
    type(C_PTR)      :: data_ptr

    ! Check the shape of the item agrees with that of data

    call get_dset_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 720 <core_hgroup:read_dset_c_sp_4_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

    ! Read the item

    mem_type_id = complex_mem_type_(SP)

    data_ptr = C_LOC(data)

    call h5dopen_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 720 <core_hgroup:read_dset_c_sp_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dopen_f'
   stop 'Program aborted'
endif

    call h5dread_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 720 <core_hgroup:read_dset_c_sp_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dread_f'
   stop 'Program aborted'
endif

    call h5dclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 720 <core_hgroup:read_dset_c_sp_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dclose_f'
   stop 'Program aborted'
endif

    ! Finish

    return

  end subroutine read_dset_c_sp_4_

  subroutine read_dset_c_dp_0_ (hg, item_name, data)

    type(hgroup_t), intent(inout)               :: hg
    character(*), intent(in)                    :: item_name
    complex(DP), target, intent(out) :: data

    integer(HID_T)   :: mem_type_id
    integer          :: hdf_err
    integer(HID_T)   :: item_id
    type(C_PTR)      :: data_ptr

    ! Check the shape of the item agrees with that of data

    ! Read the item

    mem_type_id = complex_mem_type_(DP)

    data_ptr = C_LOC(data)

    call h5dopen_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 722 <core_hgroup:read_dset_c_dp_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dopen_f'
   stop 'Program aborted'
endif

    call h5dread_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 722 <core_hgroup:read_dset_c_dp_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dread_f'
   stop 'Program aborted'
endif

    call h5dclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 722 <core_hgroup:read_dset_c_dp_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dclose_f'
   stop 'Program aborted'
endif

    ! Finish

    return

  end subroutine read_dset_c_dp_0_

  subroutine read_dset_c_dp_1_ (hg, item_name, data)

    type(hgroup_t), intent(inout)               :: hg
    character(*), intent(in)                    :: item_name
    complex(DP), target, intent(out) :: data(:)

    contiguous :: data

    integer(HSIZE_T) :: item_shape(1)

    integer(HID_T)   :: mem_type_id
    integer          :: hdf_err
    integer(HID_T)   :: item_id
    type(C_PTR)      :: data_ptr

    ! Check the shape of the item agrees with that of data

    call get_dset_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 723 <core_hgroup:read_dset_c_dp_1_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

    ! Read the item

    mem_type_id = complex_mem_type_(DP)

    data_ptr = C_LOC(data)

    call h5dopen_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 723 <core_hgroup:read_dset_c_dp_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dopen_f'
   stop 'Program aborted'
endif

    call h5dread_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 723 <core_hgroup:read_dset_c_dp_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dread_f'
   stop 'Program aborted'
endif

    call h5dclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 723 <core_hgroup:read_dset_c_dp_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dclose_f'
   stop 'Program aborted'
endif

    ! Finish

    return

  end subroutine read_dset_c_dp_1_

  subroutine read_dset_c_dp_2_ (hg, item_name, data)

    type(hgroup_t), intent(inout)               :: hg
    character(*), intent(in)                    :: item_name
    complex(DP), target, intent(out) :: data(:,:)

    contiguous :: data

    integer(HSIZE_T) :: item_shape(2)

    integer(HID_T)   :: mem_type_id
    integer          :: hdf_err
    integer(HID_T)   :: item_id
    type(C_PTR)      :: data_ptr

    ! Check the shape of the item agrees with that of data

    call get_dset_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 724 <core_hgroup:read_dset_c_dp_2_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

    ! Read the item

    mem_type_id = complex_mem_type_(DP)

    data_ptr = C_LOC(data)

    call h5dopen_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 724 <core_hgroup:read_dset_c_dp_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dopen_f'
   stop 'Program aborted'
endif

    call h5dread_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 724 <core_hgroup:read_dset_c_dp_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dread_f'
   stop 'Program aborted'
endif

    call h5dclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 724 <core_hgroup:read_dset_c_dp_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dclose_f'
   stop 'Program aborted'
endif

    ! Finish

    return

  end subroutine read_dset_c_dp_2_

  subroutine read_dset_c_dp_3_ (hg, item_name, data)

    type(hgroup_t), intent(inout)               :: hg
    character(*), intent(in)                    :: item_name
    complex(DP), target, intent(out) :: data(:,:,:)

    contiguous :: data

    integer(HSIZE_T) :: item_shape(3)

    integer(HID_T)   :: mem_type_id
    integer          :: hdf_err
    integer(HID_T)   :: item_id
    type(C_PTR)      :: data_ptr

    ! Check the shape of the item agrees with that of data

    call get_dset_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 725 <core_hgroup:read_dset_c_dp_3_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

    ! Read the item

    mem_type_id = complex_mem_type_(DP)

    data_ptr = C_LOC(data)

    call h5dopen_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 725 <core_hgroup:read_dset_c_dp_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dopen_f'
   stop 'Program aborted'
endif

    call h5dread_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 725 <core_hgroup:read_dset_c_dp_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dread_f'
   stop 'Program aborted'
endif

    call h5dclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 725 <core_hgroup:read_dset_c_dp_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dclose_f'
   stop 'Program aborted'
endif

    ! Finish

    return

  end subroutine read_dset_c_dp_3_

  subroutine read_dset_c_dp_4_ (hg, item_name, data)

    type(hgroup_t), intent(inout)               :: hg
    character(*), intent(in)                    :: item_name
    complex(DP), target, intent(out) :: data(:,:,:,:)

    contiguous :: data

    integer(HSIZE_T) :: item_shape(4)

    integer(HID_T)   :: mem_type_id
    integer          :: hdf_err
    integer(HID_T)   :: item_id
    type(C_PTR)      :: data_ptr

    ! Check the shape of the item agrees with that of data

    call get_dset_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 726 <core_hgroup:read_dset_c_dp_4_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

    ! Read the item

    mem_type_id = complex_mem_type_(DP)

    data_ptr = C_LOC(data)

    call h5dopen_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 726 <core_hgroup:read_dset_c_dp_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dopen_f'
   stop 'Program aborted'
endif

    call h5dread_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 726 <core_hgroup:read_dset_c_dp_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dread_f'
   stop 'Program aborted'
endif

    call h5dclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 726 <core_hgroup:read_dset_c_dp_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dclose_f'
   stop 'Program aborted'
endif

    ! Finish

    return

  end subroutine read_dset_c_dp_4_

!****

  subroutine read_attr_a_0_ (hg, item_name, data)

    type(hgroup_t), intent(inout)     :: hg
    character(*), intent(in)          :: item_name
    character(*), target, intent(out) :: data

    integer          :: hdf_err
    integer(HID_T)   :: mem_type_id
    type(C_PTR)      :: data_ptr
    integer(HID_T)   :: item_id

    ! Check the shape of the item agrees with that of data

    ! Read the character item

    call h5tcopy_f( H5T_NATIVE_CHARACTER,  mem_type_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 850 <core_hgroup:read_attr_a_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5tcopy_f'
   stop 'Program aborted'
endif

    call h5tset_size_f( mem_type_id,  LEN(data, SIZE_T), hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 850 <core_hgroup:read_attr_a_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5tset_size_f'
   stop 'Program aborted'
endif

    data_ptr = C_LOC(data)

    call h5aopen_name_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 850 <core_hgroup:read_attr_a_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aopen_name_f'
   stop 'Program aborted'
endif

    call h5aread_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 850 <core_hgroup:read_attr_a_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aread_f'
   stop 'Program aborted'
endif

    call h5aclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 850 <core_hgroup:read_attr_a_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aclose_f'
   stop 'Program aborted'
endif

    call h5tclose_f( mem_type_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 850 <core_hgroup:read_attr_a_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5tclose_f'
   stop 'Program aborted'
endif

    ! Finish

    return

  end subroutine read_attr_a_0_

!****

  subroutine read_attr_l_0_ (hg, item_name, data)

    type(hgroup_t), intent(inout) :: hg
    character(*), intent(in)      :: item_name
    logical, intent(out)          :: data

    integer          :: data_i

    ! Check the shape of the item agrees with that of data

    ! Read the logical item

    call read_attr(hg, item_name, data_i)

    data = data_i /= 0

    ! Finish

    return

  end subroutine read_attr_l_0_

  subroutine read_attr_a_1_ (hg, item_name, data)

    type(hgroup_t), intent(inout)     :: hg
    character(*), intent(in)          :: item_name
    character(*), target, intent(out) :: data(:)

    contiguous :: data

    integer(HSIZE_T) :: item_shape(1)

    integer          :: hdf_err
    integer(HID_T)   :: mem_type_id
    type(C_PTR)      :: data_ptr
    integer(HID_T)   :: item_id

    ! Check the shape of the item agrees with that of data

    call get_attr_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 851 <core_hgroup:read_attr_a_1_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

    ! Read the character item

    call h5tcopy_f( H5T_NATIVE_CHARACTER,  mem_type_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 851 <core_hgroup:read_attr_a_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5tcopy_f'
   stop 'Program aborted'
endif

    call h5tset_size_f( mem_type_id,  LEN(data, SIZE_T), hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 851 <core_hgroup:read_attr_a_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5tset_size_f'
   stop 'Program aborted'
endif

    data_ptr = C_LOC(data)

    call h5aopen_name_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 851 <core_hgroup:read_attr_a_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aopen_name_f'
   stop 'Program aborted'
endif

    call h5aread_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 851 <core_hgroup:read_attr_a_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aread_f'
   stop 'Program aborted'
endif

    call h5aclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 851 <core_hgroup:read_attr_a_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aclose_f'
   stop 'Program aborted'
endif

    call h5tclose_f( mem_type_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 851 <core_hgroup:read_attr_a_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5tclose_f'
   stop 'Program aborted'
endif

    ! Finish

    return

  end subroutine read_attr_a_1_

!****

  subroutine read_attr_l_1_ (hg, item_name, data)

    type(hgroup_t), intent(inout) :: hg
    character(*), intent(in)      :: item_name
    logical, intent(out)          :: data(:)

    integer(HSIZE_T) :: item_shape(1)

    integer          :: data_i(SIZE(data,1))

    ! Check the shape of the item agrees with that of data

    call get_attr_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 851 <core_hgroup:read_attr_l_1_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

    ! Read the logical item

    call read_attr(hg, item_name, data_i)

    data = data_i /= 0

    ! Finish

    return

  end subroutine read_attr_l_1_

  subroutine read_attr_a_2_ (hg, item_name, data)

    type(hgroup_t), intent(inout)     :: hg
    character(*), intent(in)          :: item_name
    character(*), target, intent(out) :: data(:,:)

    contiguous :: data

    integer(HSIZE_T) :: item_shape(2)

    integer          :: hdf_err
    integer(HID_T)   :: mem_type_id
    type(C_PTR)      :: data_ptr
    integer(HID_T)   :: item_id

    ! Check the shape of the item agrees with that of data

    call get_attr_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 852 <core_hgroup:read_attr_a_2_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

    ! Read the character item

    call h5tcopy_f( H5T_NATIVE_CHARACTER,  mem_type_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 852 <core_hgroup:read_attr_a_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5tcopy_f'
   stop 'Program aborted'
endif

    call h5tset_size_f( mem_type_id,  LEN(data, SIZE_T), hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 852 <core_hgroup:read_attr_a_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5tset_size_f'
   stop 'Program aborted'
endif

    data_ptr = C_LOC(data)

    call h5aopen_name_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 852 <core_hgroup:read_attr_a_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aopen_name_f'
   stop 'Program aborted'
endif

    call h5aread_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 852 <core_hgroup:read_attr_a_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aread_f'
   stop 'Program aborted'
endif

    call h5aclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 852 <core_hgroup:read_attr_a_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aclose_f'
   stop 'Program aborted'
endif

    call h5tclose_f( mem_type_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 852 <core_hgroup:read_attr_a_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5tclose_f'
   stop 'Program aborted'
endif

    ! Finish

    return

  end subroutine read_attr_a_2_

!****

  subroutine read_attr_l_2_ (hg, item_name, data)

    type(hgroup_t), intent(inout) :: hg
    character(*), intent(in)      :: item_name
    logical, intent(out)          :: data(:,:)

    integer(HSIZE_T) :: item_shape(2)

    integer          :: data_i(SIZE(data,1),SIZE(data,2))

    ! Check the shape of the item agrees with that of data

    call get_attr_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 852 <core_hgroup:read_attr_l_2_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

    ! Read the logical item

    call read_attr(hg, item_name, data_i)

    data = data_i /= 0

    ! Finish

    return

  end subroutine read_attr_l_2_

  subroutine read_attr_a_3_ (hg, item_name, data)

    type(hgroup_t), intent(inout)     :: hg
    character(*), intent(in)          :: item_name
    character(*), target, intent(out) :: data(:,:,:)

    contiguous :: data

    integer(HSIZE_T) :: item_shape(3)

    integer          :: hdf_err
    integer(HID_T)   :: mem_type_id
    type(C_PTR)      :: data_ptr
    integer(HID_T)   :: item_id

    ! Check the shape of the item agrees with that of data

    call get_attr_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 853 <core_hgroup:read_attr_a_3_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

    ! Read the character item

    call h5tcopy_f( H5T_NATIVE_CHARACTER,  mem_type_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 853 <core_hgroup:read_attr_a_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5tcopy_f'
   stop 'Program aborted'
endif

    call h5tset_size_f( mem_type_id,  LEN(data, SIZE_T), hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 853 <core_hgroup:read_attr_a_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5tset_size_f'
   stop 'Program aborted'
endif

    data_ptr = C_LOC(data)

    call h5aopen_name_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 853 <core_hgroup:read_attr_a_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aopen_name_f'
   stop 'Program aborted'
endif

    call h5aread_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 853 <core_hgroup:read_attr_a_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aread_f'
   stop 'Program aborted'
endif

    call h5aclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 853 <core_hgroup:read_attr_a_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aclose_f'
   stop 'Program aborted'
endif

    call h5tclose_f( mem_type_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 853 <core_hgroup:read_attr_a_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5tclose_f'
   stop 'Program aborted'
endif

    ! Finish

    return

  end subroutine read_attr_a_3_

!****

  subroutine read_attr_l_3_ (hg, item_name, data)

    type(hgroup_t), intent(inout) :: hg
    character(*), intent(in)      :: item_name
    logical, intent(out)          :: data(:,:,:)

    integer(HSIZE_T) :: item_shape(3)

    integer          :: data_i(SIZE(data,1),SIZE(data,2),SIZE(data,3))

    ! Check the shape of the item agrees with that of data

    call get_attr_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 853 <core_hgroup:read_attr_l_3_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

    ! Read the logical item

    call read_attr(hg, item_name, data_i)

    data = data_i /= 0

    ! Finish

    return

  end subroutine read_attr_l_3_

  subroutine read_attr_a_4_ (hg, item_name, data)

    type(hgroup_t), intent(inout)     :: hg
    character(*), intent(in)          :: item_name
    character(*), target, intent(out) :: data(:,:,:,:)

    contiguous :: data

    integer(HSIZE_T) :: item_shape(4)

    integer          :: hdf_err
    integer(HID_T)   :: mem_type_id
    type(C_PTR)      :: data_ptr
    integer(HID_T)   :: item_id

    ! Check the shape of the item agrees with that of data

    call get_attr_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 854 <core_hgroup:read_attr_a_4_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

    ! Read the character item

    call h5tcopy_f( H5T_NATIVE_CHARACTER,  mem_type_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 854 <core_hgroup:read_attr_a_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5tcopy_f'
   stop 'Program aborted'
endif

    call h5tset_size_f( mem_type_id,  LEN(data, SIZE_T), hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 854 <core_hgroup:read_attr_a_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5tset_size_f'
   stop 'Program aborted'
endif

    data_ptr = C_LOC(data)

    call h5aopen_name_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 854 <core_hgroup:read_attr_a_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aopen_name_f'
   stop 'Program aborted'
endif

    call h5aread_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 854 <core_hgroup:read_attr_a_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aread_f'
   stop 'Program aborted'
endif

    call h5aclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 854 <core_hgroup:read_attr_a_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aclose_f'
   stop 'Program aborted'
endif

    call h5tclose_f( mem_type_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 854 <core_hgroup:read_attr_a_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5tclose_f'
   stop 'Program aborted'
endif

    ! Finish

    return

  end subroutine read_attr_a_4_

!****

  subroutine read_attr_l_4_ (hg, item_name, data)

    type(hgroup_t), intent(inout) :: hg
    character(*), intent(in)      :: item_name
    logical, intent(out)          :: data(:,:,:,:)

    integer(HSIZE_T) :: item_shape(4)

    integer          :: data_i(SIZE(data,1),SIZE(data,2),SIZE(data,3),SIZE(data,4))

    ! Check the shape of the item agrees with that of data

    call get_attr_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 854 <core_hgroup:read_attr_l_4_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

    ! Read the logical item

    call read_attr(hg, item_name, data_i)

    data = data_i /= 0

    ! Finish

    return

  end subroutine read_attr_l_4_

  subroutine read_dset_a_0_ (hg, item_name, data)

    type(hgroup_t), intent(inout)     :: hg
    character(*), intent(in)          :: item_name
    character(*), target, intent(out) :: data

    integer          :: hdf_err
    integer(HID_T)   :: mem_type_id
    type(C_PTR)      :: data_ptr
    integer(HID_T)   :: item_id

    ! Check the shape of the item agrees with that of data

    ! Read the character item

    call h5tcopy_f( H5T_NATIVE_CHARACTER,  mem_type_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 856 <core_hgroup:read_dset_a_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5tcopy_f'
   stop 'Program aborted'
endif

    call h5tset_size_f( mem_type_id,  LEN(data, SIZE_T), hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 856 <core_hgroup:read_dset_a_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5tset_size_f'
   stop 'Program aborted'
endif

    data_ptr = C_LOC(data)

    call h5dopen_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 856 <core_hgroup:read_dset_a_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dopen_f'
   stop 'Program aborted'
endif

    call h5dread_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 856 <core_hgroup:read_dset_a_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dread_f'
   stop 'Program aborted'
endif

    call h5dclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 856 <core_hgroup:read_dset_a_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dclose_f'
   stop 'Program aborted'
endif

    call h5tclose_f( mem_type_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 856 <core_hgroup:read_dset_a_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5tclose_f'
   stop 'Program aborted'
endif

    ! Finish

    return

  end subroutine read_dset_a_0_

!****

  subroutine read_dset_l_0_ (hg, item_name, data)

    type(hgroup_t), intent(inout) :: hg
    character(*), intent(in)      :: item_name
    logical, intent(out)          :: data

    integer          :: data_i

    ! Check the shape of the item agrees with that of data

    ! Read the logical item

    call read_dset(hg, item_name, data_i)

    data = data_i /= 0

    ! Finish

    return

  end subroutine read_dset_l_0_

  subroutine read_dset_a_1_ (hg, item_name, data)

    type(hgroup_t), intent(inout)     :: hg
    character(*), intent(in)          :: item_name
    character(*), target, intent(out) :: data(:)

    contiguous :: data

    integer(HSIZE_T) :: item_shape(1)

    integer          :: hdf_err
    integer(HID_T)   :: mem_type_id
    type(C_PTR)      :: data_ptr
    integer(HID_T)   :: item_id

    ! Check the shape of the item agrees with that of data

    call get_dset_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 857 <core_hgroup:read_dset_a_1_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

    ! Read the character item

    call h5tcopy_f( H5T_NATIVE_CHARACTER,  mem_type_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 857 <core_hgroup:read_dset_a_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5tcopy_f'
   stop 'Program aborted'
endif

    call h5tset_size_f( mem_type_id,  LEN(data, SIZE_T), hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 857 <core_hgroup:read_dset_a_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5tset_size_f'
   stop 'Program aborted'
endif

    data_ptr = C_LOC(data)

    call h5dopen_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 857 <core_hgroup:read_dset_a_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dopen_f'
   stop 'Program aborted'
endif

    call h5dread_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 857 <core_hgroup:read_dset_a_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dread_f'
   stop 'Program aborted'
endif

    call h5dclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 857 <core_hgroup:read_dset_a_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dclose_f'
   stop 'Program aborted'
endif

    call h5tclose_f( mem_type_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 857 <core_hgroup:read_dset_a_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5tclose_f'
   stop 'Program aborted'
endif

    ! Finish

    return

  end subroutine read_dset_a_1_

!****

  subroutine read_dset_l_1_ (hg, item_name, data)

    type(hgroup_t), intent(inout) :: hg
    character(*), intent(in)      :: item_name
    logical, intent(out)          :: data(:)

    integer(HSIZE_T) :: item_shape(1)

    integer          :: data_i(SIZE(data,1))

    ! Check the shape of the item agrees with that of data

    call get_dset_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 857 <core_hgroup:read_dset_l_1_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

    ! Read the logical item

    call read_dset(hg, item_name, data_i)

    data = data_i /= 0

    ! Finish

    return

  end subroutine read_dset_l_1_

  subroutine read_dset_a_2_ (hg, item_name, data)

    type(hgroup_t), intent(inout)     :: hg
    character(*), intent(in)          :: item_name
    character(*), target, intent(out) :: data(:,:)

    contiguous :: data

    integer(HSIZE_T) :: item_shape(2)

    integer          :: hdf_err
    integer(HID_T)   :: mem_type_id
    type(C_PTR)      :: data_ptr
    integer(HID_T)   :: item_id

    ! Check the shape of the item agrees with that of data

    call get_dset_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 858 <core_hgroup:read_dset_a_2_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

    ! Read the character item

    call h5tcopy_f( H5T_NATIVE_CHARACTER,  mem_type_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 858 <core_hgroup:read_dset_a_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5tcopy_f'
   stop 'Program aborted'
endif

    call h5tset_size_f( mem_type_id,  LEN(data, SIZE_T), hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 858 <core_hgroup:read_dset_a_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5tset_size_f'
   stop 'Program aborted'
endif

    data_ptr = C_LOC(data)

    call h5dopen_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 858 <core_hgroup:read_dset_a_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dopen_f'
   stop 'Program aborted'
endif

    call h5dread_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 858 <core_hgroup:read_dset_a_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dread_f'
   stop 'Program aborted'
endif

    call h5dclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 858 <core_hgroup:read_dset_a_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dclose_f'
   stop 'Program aborted'
endif

    call h5tclose_f( mem_type_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 858 <core_hgroup:read_dset_a_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5tclose_f'
   stop 'Program aborted'
endif

    ! Finish

    return

  end subroutine read_dset_a_2_

!****

  subroutine read_dset_l_2_ (hg, item_name, data)

    type(hgroup_t), intent(inout) :: hg
    character(*), intent(in)      :: item_name
    logical, intent(out)          :: data(:,:)

    integer(HSIZE_T) :: item_shape(2)

    integer          :: data_i(SIZE(data,1),SIZE(data,2))

    ! Check the shape of the item agrees with that of data

    call get_dset_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 858 <core_hgroup:read_dset_l_2_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

    ! Read the logical item

    call read_dset(hg, item_name, data_i)

    data = data_i /= 0

    ! Finish

    return

  end subroutine read_dset_l_2_

  subroutine read_dset_a_3_ (hg, item_name, data)

    type(hgroup_t), intent(inout)     :: hg
    character(*), intent(in)          :: item_name
    character(*), target, intent(out) :: data(:,:,:)

    contiguous :: data

    integer(HSIZE_T) :: item_shape(3)

    integer          :: hdf_err
    integer(HID_T)   :: mem_type_id
    type(C_PTR)      :: data_ptr
    integer(HID_T)   :: item_id

    ! Check the shape of the item agrees with that of data

    call get_dset_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 859 <core_hgroup:read_dset_a_3_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

    ! Read the character item

    call h5tcopy_f( H5T_NATIVE_CHARACTER,  mem_type_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 859 <core_hgroup:read_dset_a_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5tcopy_f'
   stop 'Program aborted'
endif

    call h5tset_size_f( mem_type_id,  LEN(data, SIZE_T), hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 859 <core_hgroup:read_dset_a_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5tset_size_f'
   stop 'Program aborted'
endif

    data_ptr = C_LOC(data)

    call h5dopen_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 859 <core_hgroup:read_dset_a_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dopen_f'
   stop 'Program aborted'
endif

    call h5dread_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 859 <core_hgroup:read_dset_a_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dread_f'
   stop 'Program aborted'
endif

    call h5dclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 859 <core_hgroup:read_dset_a_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dclose_f'
   stop 'Program aborted'
endif

    call h5tclose_f( mem_type_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 859 <core_hgroup:read_dset_a_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5tclose_f'
   stop 'Program aborted'
endif

    ! Finish

    return

  end subroutine read_dset_a_3_

!****

  subroutine read_dset_l_3_ (hg, item_name, data)

    type(hgroup_t), intent(inout) :: hg
    character(*), intent(in)      :: item_name
    logical, intent(out)          :: data(:,:,:)

    integer(HSIZE_T) :: item_shape(3)

    integer          :: data_i(SIZE(data,1),SIZE(data,2),SIZE(data,3))

    ! Check the shape of the item agrees with that of data

    call get_dset_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 859 <core_hgroup:read_dset_l_3_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

    ! Read the logical item

    call read_dset(hg, item_name, data_i)

    data = data_i /= 0

    ! Finish

    return

  end subroutine read_dset_l_3_

  subroutine read_dset_a_4_ (hg, item_name, data)

    type(hgroup_t), intent(inout)     :: hg
    character(*), intent(in)          :: item_name
    character(*), target, intent(out) :: data(:,:,:,:)

    contiguous :: data

    integer(HSIZE_T) :: item_shape(4)

    integer          :: hdf_err
    integer(HID_T)   :: mem_type_id
    type(C_PTR)      :: data_ptr
    integer(HID_T)   :: item_id

    ! Check the shape of the item agrees with that of data

    call get_dset_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 860 <core_hgroup:read_dset_a_4_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

    ! Read the character item

    call h5tcopy_f( H5T_NATIVE_CHARACTER,  mem_type_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 860 <core_hgroup:read_dset_a_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5tcopy_f'
   stop 'Program aborted'
endif

    call h5tset_size_f( mem_type_id,  LEN(data, SIZE_T), hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 860 <core_hgroup:read_dset_a_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5tset_size_f'
   stop 'Program aborted'
endif

    data_ptr = C_LOC(data)

    call h5dopen_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 860 <core_hgroup:read_dset_a_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dopen_f'
   stop 'Program aborted'
endif

    call h5dread_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 860 <core_hgroup:read_dset_a_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dread_f'
   stop 'Program aborted'
endif

    call h5dclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 860 <core_hgroup:read_dset_a_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dclose_f'
   stop 'Program aborted'
endif

    call h5tclose_f( mem_type_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 860 <core_hgroup:read_dset_a_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5tclose_f'
   stop 'Program aborted'
endif

    ! Finish

    return

  end subroutine read_dset_a_4_

!****

  subroutine read_dset_l_4_ (hg, item_name, data)

    type(hgroup_t), intent(inout) :: hg
    character(*), intent(in)      :: item_name
    logical, intent(out)          :: data(:,:,:,:)

    integer(HSIZE_T) :: item_shape(4)

    integer          :: data_i(SIZE(data,1),SIZE(data,2),SIZE(data,3),SIZE(data,4))

    ! Check the shape of the item agrees with that of data

    call get_dset_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 860 <core_hgroup:read_dset_l_4_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

    ! Read the logical item

    call read_dset(hg, item_name, data_i)

    data = data_i /= 0

    ! Finish

    return

  end subroutine read_dset_l_4_

!****

  subroutine read_attr_alloc_i_i4_0_ (hg, item_name, data)

    type(hgroup_t), intent(inout)          :: hg
    character(*), intent(in)               :: item_name
    integer(I4), intent(inout), allocatable :: data

    ! If necessary, allocate the item

    if(ALLOCATED(data)) deallocate(data)

    allocate(data)

    ! Read the item

    call read_attr(hg, item_name, data)

    ! Finish

    return

  end subroutine read_attr_alloc_i_i4_0_

  subroutine read_attr_alloc_i_i4_1_ (hg, item_name, data)

    type(hgroup_t), intent(inout)          :: hg
    character(*), intent(in)               :: item_name
    integer(I4), intent(inout), allocatable :: data(:)

    integer(HSIZE_T) :: item_shape(1)

    ! If necessary, allocate the item

    if(ALLOCATED(data)) deallocate(data)

    call get_attr_shape(hg, item_name, item_shape)
    allocate(data(item_shape(1)))

    ! Read the item

    call read_attr(hg, item_name, data)

    ! Finish

    return

  end subroutine read_attr_alloc_i_i4_1_

  subroutine read_attr_alloc_i_i4_2_ (hg, item_name, data)

    type(hgroup_t), intent(inout)          :: hg
    character(*), intent(in)               :: item_name
    integer(I4), intent(inout), allocatable :: data(:,:)

    integer(HSIZE_T) :: item_shape(2)

    ! If necessary, allocate the item

    if(ALLOCATED(data)) deallocate(data)

    call get_attr_shape(hg, item_name, item_shape)
    allocate(data(item_shape(1),item_shape(2)))

    ! Read the item

    call read_attr(hg, item_name, data)

    ! Finish

    return

  end subroutine read_attr_alloc_i_i4_2_

  subroutine read_attr_alloc_i_i4_3_ (hg, item_name, data)

    type(hgroup_t), intent(inout)          :: hg
    character(*), intent(in)               :: item_name
    integer(I4), intent(inout), allocatable :: data(:,:,:)

    integer(HSIZE_T) :: item_shape(3)

    ! If necessary, allocate the item

    if(ALLOCATED(data)) deallocate(data)

    call get_attr_shape(hg, item_name, item_shape)
    allocate(data(item_shape(1),item_shape(2),item_shape(3)))

    ! Read the item

    call read_attr(hg, item_name, data)

    ! Finish

    return

  end subroutine read_attr_alloc_i_i4_3_

  subroutine read_attr_alloc_i_i4_4_ (hg, item_name, data)

    type(hgroup_t), intent(inout)          :: hg
    character(*), intent(in)               :: item_name
    integer(I4), intent(inout), allocatable :: data(:,:,:,:)

    integer(HSIZE_T) :: item_shape(4)

    ! If necessary, allocate the item

    if(ALLOCATED(data)) deallocate(data)

    call get_attr_shape(hg, item_name, item_shape)
    allocate(data(item_shape(1),item_shape(2),item_shape(3),item_shape(4)))

    ! Read the item

    call read_attr(hg, item_name, data)

    ! Finish

    return

  end subroutine read_attr_alloc_i_i4_4_

  subroutine read_attr_alloc_i_i8_0_ (hg, item_name, data)

    type(hgroup_t), intent(inout)          :: hg
    character(*), intent(in)               :: item_name
    integer(I8), intent(inout), allocatable :: data

    ! If necessary, allocate the item

    if(ALLOCATED(data)) deallocate(data)

    allocate(data)

    ! Read the item

    call read_attr(hg, item_name, data)

    ! Finish

    return

  end subroutine read_attr_alloc_i_i8_0_

  subroutine read_attr_alloc_i_i8_1_ (hg, item_name, data)

    type(hgroup_t), intent(inout)          :: hg
    character(*), intent(in)               :: item_name
    integer(I8), intent(inout), allocatable :: data(:)

    integer(HSIZE_T) :: item_shape(1)

    ! If necessary, allocate the item

    if(ALLOCATED(data)) deallocate(data)

    call get_attr_shape(hg, item_name, item_shape)
    allocate(data(item_shape(1)))

    ! Read the item

    call read_attr(hg, item_name, data)

    ! Finish

    return

  end subroutine read_attr_alloc_i_i8_1_

  subroutine read_attr_alloc_i_i8_2_ (hg, item_name, data)

    type(hgroup_t), intent(inout)          :: hg
    character(*), intent(in)               :: item_name
    integer(I8), intent(inout), allocatable :: data(:,:)

    integer(HSIZE_T) :: item_shape(2)

    ! If necessary, allocate the item

    if(ALLOCATED(data)) deallocate(data)

    call get_attr_shape(hg, item_name, item_shape)
    allocate(data(item_shape(1),item_shape(2)))

    ! Read the item

    call read_attr(hg, item_name, data)

    ! Finish

    return

  end subroutine read_attr_alloc_i_i8_2_

  subroutine read_attr_alloc_i_i8_3_ (hg, item_name, data)

    type(hgroup_t), intent(inout)          :: hg
    character(*), intent(in)               :: item_name
    integer(I8), intent(inout), allocatable :: data(:,:,:)

    integer(HSIZE_T) :: item_shape(3)

    ! If necessary, allocate the item

    if(ALLOCATED(data)) deallocate(data)

    call get_attr_shape(hg, item_name, item_shape)
    allocate(data(item_shape(1),item_shape(2),item_shape(3)))

    ! Read the item

    call read_attr(hg, item_name, data)

    ! Finish

    return

  end subroutine read_attr_alloc_i_i8_3_

  subroutine read_attr_alloc_i_i8_4_ (hg, item_name, data)

    type(hgroup_t), intent(inout)          :: hg
    character(*), intent(in)               :: item_name
    integer(I8), intent(inout), allocatable :: data(:,:,:,:)

    integer(HSIZE_T) :: item_shape(4)

    ! If necessary, allocate the item

    if(ALLOCATED(data)) deallocate(data)

    call get_attr_shape(hg, item_name, item_shape)
    allocate(data(item_shape(1),item_shape(2),item_shape(3),item_shape(4)))

    ! Read the item

    call read_attr(hg, item_name, data)

    ! Finish

    return

  end subroutine read_attr_alloc_i_i8_4_

  subroutine read_attr_alloc_r_sp_0_ (hg, item_name, data)

    type(hgroup_t), intent(inout)          :: hg
    character(*), intent(in)               :: item_name
    real(SP), intent(inout), allocatable :: data

    ! If necessary, allocate the item

    if(ALLOCATED(data)) deallocate(data)

    allocate(data)

    ! Read the item

    call read_attr(hg, item_name, data)

    ! Finish

    return

  end subroutine read_attr_alloc_r_sp_0_

  subroutine read_attr_alloc_r_sp_1_ (hg, item_name, data)

    type(hgroup_t), intent(inout)          :: hg
    character(*), intent(in)               :: item_name
    real(SP), intent(inout), allocatable :: data(:)

    integer(HSIZE_T) :: item_shape(1)

    ! If necessary, allocate the item

    if(ALLOCATED(data)) deallocate(data)

    call get_attr_shape(hg, item_name, item_shape)
    allocate(data(item_shape(1)))

    ! Read the item

    call read_attr(hg, item_name, data)

    ! Finish

    return

  end subroutine read_attr_alloc_r_sp_1_

  subroutine read_attr_alloc_r_sp_2_ (hg, item_name, data)

    type(hgroup_t), intent(inout)          :: hg
    character(*), intent(in)               :: item_name
    real(SP), intent(inout), allocatable :: data(:,:)

    integer(HSIZE_T) :: item_shape(2)

    ! If necessary, allocate the item

    if(ALLOCATED(data)) deallocate(data)

    call get_attr_shape(hg, item_name, item_shape)
    allocate(data(item_shape(1),item_shape(2)))

    ! Read the item

    call read_attr(hg, item_name, data)

    ! Finish

    return

  end subroutine read_attr_alloc_r_sp_2_

  subroutine read_attr_alloc_r_sp_3_ (hg, item_name, data)

    type(hgroup_t), intent(inout)          :: hg
    character(*), intent(in)               :: item_name
    real(SP), intent(inout), allocatable :: data(:,:,:)

    integer(HSIZE_T) :: item_shape(3)

    ! If necessary, allocate the item

    if(ALLOCATED(data)) deallocate(data)

    call get_attr_shape(hg, item_name, item_shape)
    allocate(data(item_shape(1),item_shape(2),item_shape(3)))

    ! Read the item

    call read_attr(hg, item_name, data)

    ! Finish

    return

  end subroutine read_attr_alloc_r_sp_3_

  subroutine read_attr_alloc_r_sp_4_ (hg, item_name, data)

    type(hgroup_t), intent(inout)          :: hg
    character(*), intent(in)               :: item_name
    real(SP), intent(inout), allocatable :: data(:,:,:,:)

    integer(HSIZE_T) :: item_shape(4)

    ! If necessary, allocate the item

    if(ALLOCATED(data)) deallocate(data)

    call get_attr_shape(hg, item_name, item_shape)
    allocate(data(item_shape(1),item_shape(2),item_shape(3),item_shape(4)))

    ! Read the item

    call read_attr(hg, item_name, data)

    ! Finish

    return

  end subroutine read_attr_alloc_r_sp_4_

  subroutine read_attr_alloc_r_dp_0_ (hg, item_name, data)

    type(hgroup_t), intent(inout)          :: hg
    character(*), intent(in)               :: item_name
    real(DP), intent(inout), allocatable :: data

    ! If necessary, allocate the item

    if(ALLOCATED(data)) deallocate(data)

    allocate(data)

    ! Read the item

    call read_attr(hg, item_name, data)

    ! Finish

    return

  end subroutine read_attr_alloc_r_dp_0_

  subroutine read_attr_alloc_r_dp_1_ (hg, item_name, data)

    type(hgroup_t), intent(inout)          :: hg
    character(*), intent(in)               :: item_name
    real(DP), intent(inout), allocatable :: data(:)

    integer(HSIZE_T) :: item_shape(1)

    ! If necessary, allocate the item

    if(ALLOCATED(data)) deallocate(data)

    call get_attr_shape(hg, item_name, item_shape)
    allocate(data(item_shape(1)))

    ! Read the item

    call read_attr(hg, item_name, data)

    ! Finish

    return

  end subroutine read_attr_alloc_r_dp_1_

  subroutine read_attr_alloc_r_dp_2_ (hg, item_name, data)

    type(hgroup_t), intent(inout)          :: hg
    character(*), intent(in)               :: item_name
    real(DP), intent(inout), allocatable :: data(:,:)

    integer(HSIZE_T) :: item_shape(2)

    ! If necessary, allocate the item

    if(ALLOCATED(data)) deallocate(data)

    call get_attr_shape(hg, item_name, item_shape)
    allocate(data(item_shape(1),item_shape(2)))

    ! Read the item

    call read_attr(hg, item_name, data)

    ! Finish

    return

  end subroutine read_attr_alloc_r_dp_2_

  subroutine read_attr_alloc_r_dp_3_ (hg, item_name, data)

    type(hgroup_t), intent(inout)          :: hg
    character(*), intent(in)               :: item_name
    real(DP), intent(inout), allocatable :: data(:,:,:)

    integer(HSIZE_T) :: item_shape(3)

    ! If necessary, allocate the item

    if(ALLOCATED(data)) deallocate(data)

    call get_attr_shape(hg, item_name, item_shape)
    allocate(data(item_shape(1),item_shape(2),item_shape(3)))

    ! Read the item

    call read_attr(hg, item_name, data)

    ! Finish

    return

  end subroutine read_attr_alloc_r_dp_3_

  subroutine read_attr_alloc_r_dp_4_ (hg, item_name, data)

    type(hgroup_t), intent(inout)          :: hg
    character(*), intent(in)               :: item_name
    real(DP), intent(inout), allocatable :: data(:,:,:,:)

    integer(HSIZE_T) :: item_shape(4)

    ! If necessary, allocate the item

    if(ALLOCATED(data)) deallocate(data)

    call get_attr_shape(hg, item_name, item_shape)
    allocate(data(item_shape(1),item_shape(2),item_shape(3),item_shape(4)))

    ! Read the item

    call read_attr(hg, item_name, data)

    ! Finish

    return

  end subroutine read_attr_alloc_r_dp_4_

  subroutine read_attr_alloc_c_sp_0_ (hg, item_name, data)

    type(hgroup_t), intent(inout)          :: hg
    character(*), intent(in)               :: item_name
    complex(SP), intent(inout), allocatable :: data

    ! If necessary, allocate the item

    if(ALLOCATED(data)) deallocate(data)

    allocate(data)

    ! Read the item

    call read_attr(hg, item_name, data)

    ! Finish

    return

  end subroutine read_attr_alloc_c_sp_0_

  subroutine read_attr_alloc_c_sp_1_ (hg, item_name, data)

    type(hgroup_t), intent(inout)          :: hg
    character(*), intent(in)               :: item_name
    complex(SP), intent(inout), allocatable :: data(:)

    integer(HSIZE_T) :: item_shape(1)

    ! If necessary, allocate the item

    if(ALLOCATED(data)) deallocate(data)

    call get_attr_shape(hg, item_name, item_shape)
    allocate(data(item_shape(1)))

    ! Read the item

    call read_attr(hg, item_name, data)

    ! Finish

    return

  end subroutine read_attr_alloc_c_sp_1_

  subroutine read_attr_alloc_c_sp_2_ (hg, item_name, data)

    type(hgroup_t), intent(inout)          :: hg
    character(*), intent(in)               :: item_name
    complex(SP), intent(inout), allocatable :: data(:,:)

    integer(HSIZE_T) :: item_shape(2)

    ! If necessary, allocate the item

    if(ALLOCATED(data)) deallocate(data)

    call get_attr_shape(hg, item_name, item_shape)
    allocate(data(item_shape(1),item_shape(2)))

    ! Read the item

    call read_attr(hg, item_name, data)

    ! Finish

    return

  end subroutine read_attr_alloc_c_sp_2_

  subroutine read_attr_alloc_c_sp_3_ (hg, item_name, data)

    type(hgroup_t), intent(inout)          :: hg
    character(*), intent(in)               :: item_name
    complex(SP), intent(inout), allocatable :: data(:,:,:)

    integer(HSIZE_T) :: item_shape(3)

    ! If necessary, allocate the item

    if(ALLOCATED(data)) deallocate(data)

    call get_attr_shape(hg, item_name, item_shape)
    allocate(data(item_shape(1),item_shape(2),item_shape(3)))

    ! Read the item

    call read_attr(hg, item_name, data)

    ! Finish

    return

  end subroutine read_attr_alloc_c_sp_3_

  subroutine read_attr_alloc_c_sp_4_ (hg, item_name, data)

    type(hgroup_t), intent(inout)          :: hg
    character(*), intent(in)               :: item_name
    complex(SP), intent(inout), allocatable :: data(:,:,:,:)

    integer(HSIZE_T) :: item_shape(4)

    ! If necessary, allocate the item

    if(ALLOCATED(data)) deallocate(data)

    call get_attr_shape(hg, item_name, item_shape)
    allocate(data(item_shape(1),item_shape(2),item_shape(3),item_shape(4)))

    ! Read the item

    call read_attr(hg, item_name, data)

    ! Finish

    return

  end subroutine read_attr_alloc_c_sp_4_

  subroutine read_attr_alloc_c_dp_0_ (hg, item_name, data)

    type(hgroup_t), intent(inout)          :: hg
    character(*), intent(in)               :: item_name
    complex(DP), intent(inout), allocatable :: data

    ! If necessary, allocate the item

    if(ALLOCATED(data)) deallocate(data)

    allocate(data)

    ! Read the item

    call read_attr(hg, item_name, data)

    ! Finish

    return

  end subroutine read_attr_alloc_c_dp_0_

  subroutine read_attr_alloc_c_dp_1_ (hg, item_name, data)

    type(hgroup_t), intent(inout)          :: hg
    character(*), intent(in)               :: item_name
    complex(DP), intent(inout), allocatable :: data(:)

    integer(HSIZE_T) :: item_shape(1)

    ! If necessary, allocate the item

    if(ALLOCATED(data)) deallocate(data)

    call get_attr_shape(hg, item_name, item_shape)
    allocate(data(item_shape(1)))

    ! Read the item

    call read_attr(hg, item_name, data)

    ! Finish

    return

  end subroutine read_attr_alloc_c_dp_1_

  subroutine read_attr_alloc_c_dp_2_ (hg, item_name, data)

    type(hgroup_t), intent(inout)          :: hg
    character(*), intent(in)               :: item_name
    complex(DP), intent(inout), allocatable :: data(:,:)

    integer(HSIZE_T) :: item_shape(2)

    ! If necessary, allocate the item

    if(ALLOCATED(data)) deallocate(data)

    call get_attr_shape(hg, item_name, item_shape)
    allocate(data(item_shape(1),item_shape(2)))

    ! Read the item

    call read_attr(hg, item_name, data)

    ! Finish

    return

  end subroutine read_attr_alloc_c_dp_2_

  subroutine read_attr_alloc_c_dp_3_ (hg, item_name, data)

    type(hgroup_t), intent(inout)          :: hg
    character(*), intent(in)               :: item_name
    complex(DP), intent(inout), allocatable :: data(:,:,:)

    integer(HSIZE_T) :: item_shape(3)

    ! If necessary, allocate the item

    if(ALLOCATED(data)) deallocate(data)

    call get_attr_shape(hg, item_name, item_shape)
    allocate(data(item_shape(1),item_shape(2),item_shape(3)))

    ! Read the item

    call read_attr(hg, item_name, data)

    ! Finish

    return

  end subroutine read_attr_alloc_c_dp_3_

  subroutine read_attr_alloc_c_dp_4_ (hg, item_name, data)

    type(hgroup_t), intent(inout)          :: hg
    character(*), intent(in)               :: item_name
    complex(DP), intent(inout), allocatable :: data(:,:,:,:)

    integer(HSIZE_T) :: item_shape(4)

    ! If necessary, allocate the item

    if(ALLOCATED(data)) deallocate(data)

    call get_attr_shape(hg, item_name, item_shape)
    allocate(data(item_shape(1),item_shape(2),item_shape(3),item_shape(4)))

    ! Read the item

    call read_attr(hg, item_name, data)

    ! Finish

    return

  end subroutine read_attr_alloc_c_dp_4_

  subroutine read_attr_alloc_a_0_ (hg, item_name, data)

    type(hgroup_t), intent(inout)          :: hg
    character(*), intent(in)               :: item_name
    character(*), intent(inout), allocatable :: data

    ! If necessary, allocate the item

    if(ALLOCATED(data)) deallocate(data)

    allocate(data)

    ! Read the item

    call read_attr(hg, item_name, data)

    ! Finish

    return

  end subroutine read_attr_alloc_a_0_

  subroutine read_attr_alloc_a_1_ (hg, item_name, data)

    type(hgroup_t), intent(inout)          :: hg
    character(*), intent(in)               :: item_name
    character(*), intent(inout), allocatable :: data(:)

    integer(HSIZE_T) :: item_shape(1)

    ! If necessary, allocate the item

    if(ALLOCATED(data)) deallocate(data)

    call get_attr_shape(hg, item_name, item_shape)
    allocate(data(item_shape(1)))

    ! Read the item

    call read_attr(hg, item_name, data)

    ! Finish

    return

  end subroutine read_attr_alloc_a_1_

  subroutine read_attr_alloc_a_2_ (hg, item_name, data)

    type(hgroup_t), intent(inout)          :: hg
    character(*), intent(in)               :: item_name
    character(*), intent(inout), allocatable :: data(:,:)

    integer(HSIZE_T) :: item_shape(2)

    ! If necessary, allocate the item

    if(ALLOCATED(data)) deallocate(data)

    call get_attr_shape(hg, item_name, item_shape)
    allocate(data(item_shape(1),item_shape(2)))

    ! Read the item

    call read_attr(hg, item_name, data)

    ! Finish

    return

  end subroutine read_attr_alloc_a_2_

  subroutine read_attr_alloc_a_3_ (hg, item_name, data)

    type(hgroup_t), intent(inout)          :: hg
    character(*), intent(in)               :: item_name
    character(*), intent(inout), allocatable :: data(:,:,:)

    integer(HSIZE_T) :: item_shape(3)

    ! If necessary, allocate the item

    if(ALLOCATED(data)) deallocate(data)

    call get_attr_shape(hg, item_name, item_shape)
    allocate(data(item_shape(1),item_shape(2),item_shape(3)))

    ! Read the item

    call read_attr(hg, item_name, data)

    ! Finish

    return

  end subroutine read_attr_alloc_a_3_

  subroutine read_attr_alloc_a_4_ (hg, item_name, data)

    type(hgroup_t), intent(inout)          :: hg
    character(*), intent(in)               :: item_name
    character(*), intent(inout), allocatable :: data(:,:,:,:)

    integer(HSIZE_T) :: item_shape(4)

    ! If necessary, allocate the item

    if(ALLOCATED(data)) deallocate(data)

    call get_attr_shape(hg, item_name, item_shape)
    allocate(data(item_shape(1),item_shape(2),item_shape(3),item_shape(4)))

    ! Read the item

    call read_attr(hg, item_name, data)

    ! Finish

    return

  end subroutine read_attr_alloc_a_4_

  subroutine read_attr_alloc_l_0_ (hg, item_name, data)

    type(hgroup_t), intent(inout)          :: hg
    character(*), intent(in)               :: item_name
    logical, intent(inout), allocatable :: data

    ! If necessary, allocate the item

    if(ALLOCATED(data)) deallocate(data)

    allocate(data)

    ! Read the item

    call read_attr(hg, item_name, data)

    ! Finish

    return

  end subroutine read_attr_alloc_l_0_

  subroutine read_attr_alloc_l_1_ (hg, item_name, data)

    type(hgroup_t), intent(inout)          :: hg
    character(*), intent(in)               :: item_name
    logical, intent(inout), allocatable :: data(:)

    integer(HSIZE_T) :: item_shape(1)

    ! If necessary, allocate the item

    if(ALLOCATED(data)) deallocate(data)

    call get_attr_shape(hg, item_name, item_shape)
    allocate(data(item_shape(1)))

    ! Read the item

    call read_attr(hg, item_name, data)

    ! Finish

    return

  end subroutine read_attr_alloc_l_1_

  subroutine read_attr_alloc_l_2_ (hg, item_name, data)

    type(hgroup_t), intent(inout)          :: hg
    character(*), intent(in)               :: item_name
    logical, intent(inout), allocatable :: data(:,:)

    integer(HSIZE_T) :: item_shape(2)

    ! If necessary, allocate the item

    if(ALLOCATED(data)) deallocate(data)

    call get_attr_shape(hg, item_name, item_shape)
    allocate(data(item_shape(1),item_shape(2)))

    ! Read the item

    call read_attr(hg, item_name, data)

    ! Finish

    return

  end subroutine read_attr_alloc_l_2_

  subroutine read_attr_alloc_l_3_ (hg, item_name, data)

    type(hgroup_t), intent(inout)          :: hg
    character(*), intent(in)               :: item_name
    logical, intent(inout), allocatable :: data(:,:,:)

    integer(HSIZE_T) :: item_shape(3)

    ! If necessary, allocate the item

    if(ALLOCATED(data)) deallocate(data)

    call get_attr_shape(hg, item_name, item_shape)
    allocate(data(item_shape(1),item_shape(2),item_shape(3)))

    ! Read the item

    call read_attr(hg, item_name, data)

    ! Finish

    return

  end subroutine read_attr_alloc_l_3_

  subroutine read_attr_alloc_l_4_ (hg, item_name, data)

    type(hgroup_t), intent(inout)          :: hg
    character(*), intent(in)               :: item_name
    logical, intent(inout), allocatable :: data(:,:,:,:)

    integer(HSIZE_T) :: item_shape(4)

    ! If necessary, allocate the item

    if(ALLOCATED(data)) deallocate(data)

    call get_attr_shape(hg, item_name, item_shape)
    allocate(data(item_shape(1),item_shape(2),item_shape(3),item_shape(4)))

    ! Read the item

    call read_attr(hg, item_name, data)

    ! Finish

    return

  end subroutine read_attr_alloc_l_4_

  subroutine read_dset_alloc_i_i4_0_ (hg, item_name, data)

    type(hgroup_t), intent(inout)          :: hg
    character(*), intent(in)               :: item_name
    integer(I4), intent(inout), allocatable :: data

    ! If necessary, allocate the item

    if(ALLOCATED(data)) deallocate(data)

    allocate(data)

    ! Read the item

    call read_dset(hg, item_name, data)

    ! Finish

    return

  end subroutine read_dset_alloc_i_i4_0_

  subroutine read_dset_alloc_i_i4_1_ (hg, item_name, data)

    type(hgroup_t), intent(inout)          :: hg
    character(*), intent(in)               :: item_name
    integer(I4), intent(inout), allocatable :: data(:)

    integer(HSIZE_T) :: item_shape(1)

    ! If necessary, allocate the item

    if(ALLOCATED(data)) deallocate(data)

    call get_dset_shape(hg, item_name, item_shape)
    allocate(data(item_shape(1)))

    ! Read the item

    call read_dset(hg, item_name, data)

    ! Finish

    return

  end subroutine read_dset_alloc_i_i4_1_

  subroutine read_dset_alloc_i_i4_2_ (hg, item_name, data)

    type(hgroup_t), intent(inout)          :: hg
    character(*), intent(in)               :: item_name
    integer(I4), intent(inout), allocatable :: data(:,:)

    integer(HSIZE_T) :: item_shape(2)

    ! If necessary, allocate the item

    if(ALLOCATED(data)) deallocate(data)

    call get_dset_shape(hg, item_name, item_shape)
    allocate(data(item_shape(1),item_shape(2)))

    ! Read the item

    call read_dset(hg, item_name, data)

    ! Finish

    return

  end subroutine read_dset_alloc_i_i4_2_

  subroutine read_dset_alloc_i_i4_3_ (hg, item_name, data)

    type(hgroup_t), intent(inout)          :: hg
    character(*), intent(in)               :: item_name
    integer(I4), intent(inout), allocatable :: data(:,:,:)

    integer(HSIZE_T) :: item_shape(3)

    ! If necessary, allocate the item

    if(ALLOCATED(data)) deallocate(data)

    call get_dset_shape(hg, item_name, item_shape)
    allocate(data(item_shape(1),item_shape(2),item_shape(3)))

    ! Read the item

    call read_dset(hg, item_name, data)

    ! Finish

    return

  end subroutine read_dset_alloc_i_i4_3_

  subroutine read_dset_alloc_i_i4_4_ (hg, item_name, data)

    type(hgroup_t), intent(inout)          :: hg
    character(*), intent(in)               :: item_name
    integer(I4), intent(inout), allocatable :: data(:,:,:,:)

    integer(HSIZE_T) :: item_shape(4)

    ! If necessary, allocate the item

    if(ALLOCATED(data)) deallocate(data)

    call get_dset_shape(hg, item_name, item_shape)
    allocate(data(item_shape(1),item_shape(2),item_shape(3),item_shape(4)))

    ! Read the item

    call read_dset(hg, item_name, data)

    ! Finish

    return

  end subroutine read_dset_alloc_i_i4_4_

  subroutine read_dset_alloc_i_i8_0_ (hg, item_name, data)

    type(hgroup_t), intent(inout)          :: hg
    character(*), intent(in)               :: item_name
    integer(I8), intent(inout), allocatable :: data

    ! If necessary, allocate the item

    if(ALLOCATED(data)) deallocate(data)

    allocate(data)

    ! Read the item

    call read_dset(hg, item_name, data)

    ! Finish

    return

  end subroutine read_dset_alloc_i_i8_0_

  subroutine read_dset_alloc_i_i8_1_ (hg, item_name, data)

    type(hgroup_t), intent(inout)          :: hg
    character(*), intent(in)               :: item_name
    integer(I8), intent(inout), allocatable :: data(:)

    integer(HSIZE_T) :: item_shape(1)

    ! If necessary, allocate the item

    if(ALLOCATED(data)) deallocate(data)

    call get_dset_shape(hg, item_name, item_shape)
    allocate(data(item_shape(1)))

    ! Read the item

    call read_dset(hg, item_name, data)

    ! Finish

    return

  end subroutine read_dset_alloc_i_i8_1_

  subroutine read_dset_alloc_i_i8_2_ (hg, item_name, data)

    type(hgroup_t), intent(inout)          :: hg
    character(*), intent(in)               :: item_name
    integer(I8), intent(inout), allocatable :: data(:,:)

    integer(HSIZE_T) :: item_shape(2)

    ! If necessary, allocate the item

    if(ALLOCATED(data)) deallocate(data)

    call get_dset_shape(hg, item_name, item_shape)
    allocate(data(item_shape(1),item_shape(2)))

    ! Read the item

    call read_dset(hg, item_name, data)

    ! Finish

    return

  end subroutine read_dset_alloc_i_i8_2_

  subroutine read_dset_alloc_i_i8_3_ (hg, item_name, data)

    type(hgroup_t), intent(inout)          :: hg
    character(*), intent(in)               :: item_name
    integer(I8), intent(inout), allocatable :: data(:,:,:)

    integer(HSIZE_T) :: item_shape(3)

    ! If necessary, allocate the item

    if(ALLOCATED(data)) deallocate(data)

    call get_dset_shape(hg, item_name, item_shape)
    allocate(data(item_shape(1),item_shape(2),item_shape(3)))

    ! Read the item

    call read_dset(hg, item_name, data)

    ! Finish

    return

  end subroutine read_dset_alloc_i_i8_3_

  subroutine read_dset_alloc_i_i8_4_ (hg, item_name, data)

    type(hgroup_t), intent(inout)          :: hg
    character(*), intent(in)               :: item_name
    integer(I8), intent(inout), allocatable :: data(:,:,:,:)

    integer(HSIZE_T) :: item_shape(4)

    ! If necessary, allocate the item

    if(ALLOCATED(data)) deallocate(data)

    call get_dset_shape(hg, item_name, item_shape)
    allocate(data(item_shape(1),item_shape(2),item_shape(3),item_shape(4)))

    ! Read the item

    call read_dset(hg, item_name, data)

    ! Finish

    return

  end subroutine read_dset_alloc_i_i8_4_

  subroutine read_dset_alloc_r_sp_0_ (hg, item_name, data)

    type(hgroup_t), intent(inout)          :: hg
    character(*), intent(in)               :: item_name
    real(SP), intent(inout), allocatable :: data

    ! If necessary, allocate the item

    if(ALLOCATED(data)) deallocate(data)

    allocate(data)

    ! Read the item

    call read_dset(hg, item_name, data)

    ! Finish

    return

  end subroutine read_dset_alloc_r_sp_0_

  subroutine read_dset_alloc_r_sp_1_ (hg, item_name, data)

    type(hgroup_t), intent(inout)          :: hg
    character(*), intent(in)               :: item_name
    real(SP), intent(inout), allocatable :: data(:)

    integer(HSIZE_T) :: item_shape(1)

    ! If necessary, allocate the item

    if(ALLOCATED(data)) deallocate(data)

    call get_dset_shape(hg, item_name, item_shape)
    allocate(data(item_shape(1)))

    ! Read the item

    call read_dset(hg, item_name, data)

    ! Finish

    return

  end subroutine read_dset_alloc_r_sp_1_

  subroutine read_dset_alloc_r_sp_2_ (hg, item_name, data)

    type(hgroup_t), intent(inout)          :: hg
    character(*), intent(in)               :: item_name
    real(SP), intent(inout), allocatable :: data(:,:)

    integer(HSIZE_T) :: item_shape(2)

    ! If necessary, allocate the item

    if(ALLOCATED(data)) deallocate(data)

    call get_dset_shape(hg, item_name, item_shape)
    allocate(data(item_shape(1),item_shape(2)))

    ! Read the item

    call read_dset(hg, item_name, data)

    ! Finish

    return

  end subroutine read_dset_alloc_r_sp_2_

  subroutine read_dset_alloc_r_sp_3_ (hg, item_name, data)

    type(hgroup_t), intent(inout)          :: hg
    character(*), intent(in)               :: item_name
    real(SP), intent(inout), allocatable :: data(:,:,:)

    integer(HSIZE_T) :: item_shape(3)

    ! If necessary, allocate the item

    if(ALLOCATED(data)) deallocate(data)

    call get_dset_shape(hg, item_name, item_shape)
    allocate(data(item_shape(1),item_shape(2),item_shape(3)))

    ! Read the item

    call read_dset(hg, item_name, data)

    ! Finish

    return

  end subroutine read_dset_alloc_r_sp_3_

  subroutine read_dset_alloc_r_sp_4_ (hg, item_name, data)

    type(hgroup_t), intent(inout)          :: hg
    character(*), intent(in)               :: item_name
    real(SP), intent(inout), allocatable :: data(:,:,:,:)

    integer(HSIZE_T) :: item_shape(4)

    ! If necessary, allocate the item

    if(ALLOCATED(data)) deallocate(data)

    call get_dset_shape(hg, item_name, item_shape)
    allocate(data(item_shape(1),item_shape(2),item_shape(3),item_shape(4)))

    ! Read the item

    call read_dset(hg, item_name, data)

    ! Finish

    return

  end subroutine read_dset_alloc_r_sp_4_

  subroutine read_dset_alloc_r_dp_0_ (hg, item_name, data)

    type(hgroup_t), intent(inout)          :: hg
    character(*), intent(in)               :: item_name
    real(DP), intent(inout), allocatable :: data

    ! If necessary, allocate the item

    if(ALLOCATED(data)) deallocate(data)

    allocate(data)

    ! Read the item

    call read_dset(hg, item_name, data)

    ! Finish

    return

  end subroutine read_dset_alloc_r_dp_0_

  subroutine read_dset_alloc_r_dp_1_ (hg, item_name, data)

    type(hgroup_t), intent(inout)          :: hg
    character(*), intent(in)               :: item_name
    real(DP), intent(inout), allocatable :: data(:)

    integer(HSIZE_T) :: item_shape(1)

    ! If necessary, allocate the item

    if(ALLOCATED(data)) deallocate(data)

    call get_dset_shape(hg, item_name, item_shape)
    allocate(data(item_shape(1)))

    ! Read the item

    call read_dset(hg, item_name, data)

    ! Finish

    return

  end subroutine read_dset_alloc_r_dp_1_

  subroutine read_dset_alloc_r_dp_2_ (hg, item_name, data)

    type(hgroup_t), intent(inout)          :: hg
    character(*), intent(in)               :: item_name
    real(DP), intent(inout), allocatable :: data(:,:)

    integer(HSIZE_T) :: item_shape(2)

    ! If necessary, allocate the item

    if(ALLOCATED(data)) deallocate(data)

    call get_dset_shape(hg, item_name, item_shape)
    allocate(data(item_shape(1),item_shape(2)))

    ! Read the item

    call read_dset(hg, item_name, data)

    ! Finish

    return

  end subroutine read_dset_alloc_r_dp_2_

  subroutine read_dset_alloc_r_dp_3_ (hg, item_name, data)

    type(hgroup_t), intent(inout)          :: hg
    character(*), intent(in)               :: item_name
    real(DP), intent(inout), allocatable :: data(:,:,:)

    integer(HSIZE_T) :: item_shape(3)

    ! If necessary, allocate the item

    if(ALLOCATED(data)) deallocate(data)

    call get_dset_shape(hg, item_name, item_shape)
    allocate(data(item_shape(1),item_shape(2),item_shape(3)))

    ! Read the item

    call read_dset(hg, item_name, data)

    ! Finish

    return

  end subroutine read_dset_alloc_r_dp_3_

  subroutine read_dset_alloc_r_dp_4_ (hg, item_name, data)

    type(hgroup_t), intent(inout)          :: hg
    character(*), intent(in)               :: item_name
    real(DP), intent(inout), allocatable :: data(:,:,:,:)

    integer(HSIZE_T) :: item_shape(4)

    ! If necessary, allocate the item

    if(ALLOCATED(data)) deallocate(data)

    call get_dset_shape(hg, item_name, item_shape)
    allocate(data(item_shape(1),item_shape(2),item_shape(3),item_shape(4)))

    ! Read the item

    call read_dset(hg, item_name, data)

    ! Finish

    return

  end subroutine read_dset_alloc_r_dp_4_

  subroutine read_dset_alloc_c_sp_0_ (hg, item_name, data)

    type(hgroup_t), intent(inout)          :: hg
    character(*), intent(in)               :: item_name
    complex(SP), intent(inout), allocatable :: data

    ! If necessary, allocate the item

    if(ALLOCATED(data)) deallocate(data)

    allocate(data)

    ! Read the item

    call read_dset(hg, item_name, data)

    ! Finish

    return

  end subroutine read_dset_alloc_c_sp_0_

  subroutine read_dset_alloc_c_sp_1_ (hg, item_name, data)

    type(hgroup_t), intent(inout)          :: hg
    character(*), intent(in)               :: item_name
    complex(SP), intent(inout), allocatable :: data(:)

    integer(HSIZE_T) :: item_shape(1)

    ! If necessary, allocate the item

    if(ALLOCATED(data)) deallocate(data)

    call get_dset_shape(hg, item_name, item_shape)
    allocate(data(item_shape(1)))

    ! Read the item

    call read_dset(hg, item_name, data)

    ! Finish

    return

  end subroutine read_dset_alloc_c_sp_1_

  subroutine read_dset_alloc_c_sp_2_ (hg, item_name, data)

    type(hgroup_t), intent(inout)          :: hg
    character(*), intent(in)               :: item_name
    complex(SP), intent(inout), allocatable :: data(:,:)

    integer(HSIZE_T) :: item_shape(2)

    ! If necessary, allocate the item

    if(ALLOCATED(data)) deallocate(data)

    call get_dset_shape(hg, item_name, item_shape)
    allocate(data(item_shape(1),item_shape(2)))

    ! Read the item

    call read_dset(hg, item_name, data)

    ! Finish

    return

  end subroutine read_dset_alloc_c_sp_2_

  subroutine read_dset_alloc_c_sp_3_ (hg, item_name, data)

    type(hgroup_t), intent(inout)          :: hg
    character(*), intent(in)               :: item_name
    complex(SP), intent(inout), allocatable :: data(:,:,:)

    integer(HSIZE_T) :: item_shape(3)

    ! If necessary, allocate the item

    if(ALLOCATED(data)) deallocate(data)

    call get_dset_shape(hg, item_name, item_shape)
    allocate(data(item_shape(1),item_shape(2),item_shape(3)))

    ! Read the item

    call read_dset(hg, item_name, data)

    ! Finish

    return

  end subroutine read_dset_alloc_c_sp_3_

  subroutine read_dset_alloc_c_sp_4_ (hg, item_name, data)

    type(hgroup_t), intent(inout)          :: hg
    character(*), intent(in)               :: item_name
    complex(SP), intent(inout), allocatable :: data(:,:,:,:)

    integer(HSIZE_T) :: item_shape(4)

    ! If necessary, allocate the item

    if(ALLOCATED(data)) deallocate(data)

    call get_dset_shape(hg, item_name, item_shape)
    allocate(data(item_shape(1),item_shape(2),item_shape(3),item_shape(4)))

    ! Read the item

    call read_dset(hg, item_name, data)

    ! Finish

    return

  end subroutine read_dset_alloc_c_sp_4_

  subroutine read_dset_alloc_c_dp_0_ (hg, item_name, data)

    type(hgroup_t), intent(inout)          :: hg
    character(*), intent(in)               :: item_name
    complex(DP), intent(inout), allocatable :: data

    ! If necessary, allocate the item

    if(ALLOCATED(data)) deallocate(data)

    allocate(data)

    ! Read the item

    call read_dset(hg, item_name, data)

    ! Finish

    return

  end subroutine read_dset_alloc_c_dp_0_

  subroutine read_dset_alloc_c_dp_1_ (hg, item_name, data)

    type(hgroup_t), intent(inout)          :: hg
    character(*), intent(in)               :: item_name
    complex(DP), intent(inout), allocatable :: data(:)

    integer(HSIZE_T) :: item_shape(1)

    ! If necessary, allocate the item

    if(ALLOCATED(data)) deallocate(data)

    call get_dset_shape(hg, item_name, item_shape)
    allocate(data(item_shape(1)))

    ! Read the item

    call read_dset(hg, item_name, data)

    ! Finish

    return

  end subroutine read_dset_alloc_c_dp_1_

  subroutine read_dset_alloc_c_dp_2_ (hg, item_name, data)

    type(hgroup_t), intent(inout)          :: hg
    character(*), intent(in)               :: item_name
    complex(DP), intent(inout), allocatable :: data(:,:)

    integer(HSIZE_T) :: item_shape(2)

    ! If necessary, allocate the item

    if(ALLOCATED(data)) deallocate(data)

    call get_dset_shape(hg, item_name, item_shape)
    allocate(data(item_shape(1),item_shape(2)))

    ! Read the item

    call read_dset(hg, item_name, data)

    ! Finish

    return

  end subroutine read_dset_alloc_c_dp_2_

  subroutine read_dset_alloc_c_dp_3_ (hg, item_name, data)

    type(hgroup_t), intent(inout)          :: hg
    character(*), intent(in)               :: item_name
    complex(DP), intent(inout), allocatable :: data(:,:,:)

    integer(HSIZE_T) :: item_shape(3)

    ! If necessary, allocate the item

    if(ALLOCATED(data)) deallocate(data)

    call get_dset_shape(hg, item_name, item_shape)
    allocate(data(item_shape(1),item_shape(2),item_shape(3)))

    ! Read the item

    call read_dset(hg, item_name, data)

    ! Finish

    return

  end subroutine read_dset_alloc_c_dp_3_

  subroutine read_dset_alloc_c_dp_4_ (hg, item_name, data)

    type(hgroup_t), intent(inout)          :: hg
    character(*), intent(in)               :: item_name
    complex(DP), intent(inout), allocatable :: data(:,:,:,:)

    integer(HSIZE_T) :: item_shape(4)

    ! If necessary, allocate the item

    if(ALLOCATED(data)) deallocate(data)

    call get_dset_shape(hg, item_name, item_shape)
    allocate(data(item_shape(1),item_shape(2),item_shape(3),item_shape(4)))

    ! Read the item

    call read_dset(hg, item_name, data)

    ! Finish

    return

  end subroutine read_dset_alloc_c_dp_4_

  subroutine read_dset_alloc_a_0_ (hg, item_name, data)

    type(hgroup_t), intent(inout)          :: hg
    character(*), intent(in)               :: item_name
    character(*), intent(inout), allocatable :: data

    ! If necessary, allocate the item

    if(ALLOCATED(data)) deallocate(data)

    allocate(data)

    ! Read the item

    call read_dset(hg, item_name, data)

    ! Finish

    return

  end subroutine read_dset_alloc_a_0_

  subroutine read_dset_alloc_a_1_ (hg, item_name, data)

    type(hgroup_t), intent(inout)          :: hg
    character(*), intent(in)               :: item_name
    character(*), intent(inout), allocatable :: data(:)

    integer(HSIZE_T) :: item_shape(1)

    ! If necessary, allocate the item

    if(ALLOCATED(data)) deallocate(data)

    call get_dset_shape(hg, item_name, item_shape)
    allocate(data(item_shape(1)))

    ! Read the item

    call read_dset(hg, item_name, data)

    ! Finish

    return

  end subroutine read_dset_alloc_a_1_

  subroutine read_dset_alloc_a_2_ (hg, item_name, data)

    type(hgroup_t), intent(inout)          :: hg
    character(*), intent(in)               :: item_name
    character(*), intent(inout), allocatable :: data(:,:)

    integer(HSIZE_T) :: item_shape(2)

    ! If necessary, allocate the item

    if(ALLOCATED(data)) deallocate(data)

    call get_dset_shape(hg, item_name, item_shape)
    allocate(data(item_shape(1),item_shape(2)))

    ! Read the item

    call read_dset(hg, item_name, data)

    ! Finish

    return

  end subroutine read_dset_alloc_a_2_

  subroutine read_dset_alloc_a_3_ (hg, item_name, data)

    type(hgroup_t), intent(inout)          :: hg
    character(*), intent(in)               :: item_name
    character(*), intent(inout), allocatable :: data(:,:,:)

    integer(HSIZE_T) :: item_shape(3)

    ! If necessary, allocate the item

    if(ALLOCATED(data)) deallocate(data)

    call get_dset_shape(hg, item_name, item_shape)
    allocate(data(item_shape(1),item_shape(2),item_shape(3)))

    ! Read the item

    call read_dset(hg, item_name, data)

    ! Finish

    return

  end subroutine read_dset_alloc_a_3_

  subroutine read_dset_alloc_a_4_ (hg, item_name, data)

    type(hgroup_t), intent(inout)          :: hg
    character(*), intent(in)               :: item_name
    character(*), intent(inout), allocatable :: data(:,:,:,:)

    integer(HSIZE_T) :: item_shape(4)

    ! If necessary, allocate the item

    if(ALLOCATED(data)) deallocate(data)

    call get_dset_shape(hg, item_name, item_shape)
    allocate(data(item_shape(1),item_shape(2),item_shape(3),item_shape(4)))

    ! Read the item

    call read_dset(hg, item_name, data)

    ! Finish

    return

  end subroutine read_dset_alloc_a_4_

  subroutine read_dset_alloc_l_0_ (hg, item_name, data)

    type(hgroup_t), intent(inout)          :: hg
    character(*), intent(in)               :: item_name
    logical, intent(inout), allocatable :: data

    ! If necessary, allocate the item

    if(ALLOCATED(data)) deallocate(data)

    allocate(data)

    ! Read the item

    call read_dset(hg, item_name, data)

    ! Finish

    return

  end subroutine read_dset_alloc_l_0_

  subroutine read_dset_alloc_l_1_ (hg, item_name, data)

    type(hgroup_t), intent(inout)          :: hg
    character(*), intent(in)               :: item_name
    logical, intent(inout), allocatable :: data(:)

    integer(HSIZE_T) :: item_shape(1)

    ! If necessary, allocate the item

    if(ALLOCATED(data)) deallocate(data)

    call get_dset_shape(hg, item_name, item_shape)
    allocate(data(item_shape(1)))

    ! Read the item

    call read_dset(hg, item_name, data)

    ! Finish

    return

  end subroutine read_dset_alloc_l_1_

  subroutine read_dset_alloc_l_2_ (hg, item_name, data)

    type(hgroup_t), intent(inout)          :: hg
    character(*), intent(in)               :: item_name
    logical, intent(inout), allocatable :: data(:,:)

    integer(HSIZE_T) :: item_shape(2)

    ! If necessary, allocate the item

    if(ALLOCATED(data)) deallocate(data)

    call get_dset_shape(hg, item_name, item_shape)
    allocate(data(item_shape(1),item_shape(2)))

    ! Read the item

    call read_dset(hg, item_name, data)

    ! Finish

    return

  end subroutine read_dset_alloc_l_2_

  subroutine read_dset_alloc_l_3_ (hg, item_name, data)

    type(hgroup_t), intent(inout)          :: hg
    character(*), intent(in)               :: item_name
    logical, intent(inout), allocatable :: data(:,:,:)

    integer(HSIZE_T) :: item_shape(3)

    ! If necessary, allocate the item

    if(ALLOCATED(data)) deallocate(data)

    call get_dset_shape(hg, item_name, item_shape)
    allocate(data(item_shape(1),item_shape(2),item_shape(3)))

    ! Read the item

    call read_dset(hg, item_name, data)

    ! Finish

    return

  end subroutine read_dset_alloc_l_3_

  subroutine read_dset_alloc_l_4_ (hg, item_name, data)

    type(hgroup_t), intent(inout)          :: hg
    character(*), intent(in)               :: item_name
    logical, intent(inout), allocatable :: data(:,:,:,:)

    integer(HSIZE_T) :: item_shape(4)

    ! If necessary, allocate the item

    if(ALLOCATED(data)) deallocate(data)

    call get_dset_shape(hg, item_name, item_shape)
    allocate(data(item_shape(1),item_shape(2),item_shape(3),item_shape(4)))

    ! Read the item

    call read_dset(hg, item_name, data)

    ! Finish

    return

  end subroutine read_dset_alloc_l_4_

!****

  subroutine write_attr_i_i4_0_ (hg, item_name, data, overwrite, comp_level)

    type(hgroup_t), intent(inout)              :: hg
    character(*), intent(in)                   :: item_name
    integer(I4), target, intent(in) :: data

    logical, intent(in), optional              :: overwrite
    integer, intent(in), optional              :: comp_level

    logical          :: overwrite_
    integer          :: hdf_err
    logical          :: item_exists
    integer(HSIZE_T) :: item_shape(0)
    integer(HID_T)   :: dspace_id
    integer(HID_T)   :: plist_id
    integer(HID_T)   :: mem_type_id
    integer(HID_T)   :: file_type_id
    type(C_PTR)      :: data_ptr
    integer(HID_T)   :: item_id

    if(PRESENT(overwrite)) then
       overwrite_ = overwrite
    else
       overwrite_ = .FALSE.
    endif

    ! Write the item

    item_exists = attr_exists_(hg, item_name)
    overwrite_ = overwrite_ .AND. item_exists

    if(overwrite_) then

       call get_attr_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 1143 <core_hgroup:write_attr_i_i4' &
 & // '_0_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

       mem_type_id = integer_mem_type_(I4)

       data_ptr = C_LOC(data)

       call h5aopen_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1143 <core_hgroup:write_attr_i_i4_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aopen_f'
   stop 'Program aborted'
endif

       call h5awrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1143 <core_hgroup:write_attr_i_i4_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5awrite_f'
   stop 'Program aborted'
endif

       call h5aclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1143 <core_hgroup:write_attr_i_i4_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aclose_f'
   stop 'Program aborted'
endif

    else

       call h5screate_f( H5S_SCALAR_F,  dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1143 <core_hgroup:write_attr_i_i4_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5screate_f'
   stop 'Program aborted'
endif

       call h5pcreate_f( H5P_ATTRIBUTE_CREATE_F,  plist_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1143 <core_hgroup:write_attr_i_i4_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pcreate_f'
   stop 'Program aborted'
endif

    if(.NOT. (.NOT. PRESENT(comp_level))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''.NOT. PRESENT(comp_level)'' failed at line 1143 <core_hgroup:write_attr_i_i4_0_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Attributes cannot be compressed'
      stop
    endif

       mem_type_id = integer_mem_type_(I4)
       file_type_id = integer_file_type_(I4)

       data_ptr = C_LOC(data)

       call h5acreate_f( hg%group_id,  item_name,  file_type_id,  dspace_id,  item_id, hdf_err,  acpl_id=plist_id)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1143 <core_hgroup:write_attr_i_i4_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5acreate_f'
   stop 'Program aborted'
endif

       call h5awrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1143 <core_hgroup:write_attr_i_i4_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5awrite_f'
   stop 'Program aborted'
endif

       call h5pclose_f( plist_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1143 <core_hgroup:write_attr_i_i4_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pclose_f'
   stop 'Program aborted'
endif

       call h5aclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1143 <core_hgroup:write_attr_i_i4_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aclose_f'
   stop 'Program aborted'
endif

       call h5sclose_f( dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1143 <core_hgroup:write_attr_i_i4_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5sclose_f'
   stop 'Program aborted'
endif

    endif

    ! Finish

    return

  end subroutine write_attr_i_i4_0_

  subroutine write_attr_i_i4_1_ (hg, item_name, data, overwrite, comp_level)

    type(hgroup_t), intent(inout)              :: hg
    character(*), intent(in)                   :: item_name
    integer(I4), target, intent(in) :: data(:)

    contiguous :: data

    logical, intent(in), optional              :: overwrite
    integer, intent(in), optional              :: comp_level

    logical          :: overwrite_
    integer          :: hdf_err
    logical          :: item_exists
    integer(HSIZE_T) :: item_shape(1)
    integer(HID_T)   :: dspace_id
    integer(HID_T)   :: plist_id
    integer(HID_T)   :: mem_type_id
    integer(HID_T)   :: file_type_id
    type(C_PTR)      :: data_ptr
    integer(HID_T)   :: item_id

    if(PRESENT(overwrite)) then
       overwrite_ = overwrite
    else
       overwrite_ = .FALSE.
    endif

    ! Write the item

    item_exists = attr_exists_(hg, item_name)
    overwrite_ = overwrite_ .AND. item_exists

    if(overwrite_) then

       call get_attr_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 1144 <core_hgroup:write_attr_i_i4' &
 & // '_1_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

       mem_type_id = integer_mem_type_(I4)

       data_ptr = C_LOC(data)

       call h5aopen_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1144 <core_hgroup:write_attr_i_i4_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aopen_f'
   stop 'Program aborted'
endif

       call h5awrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1144 <core_hgroup:write_attr_i_i4_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5awrite_f'
   stop 'Program aborted'
endif

       call h5aclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1144 <core_hgroup:write_attr_i_i4_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aclose_f'
   stop 'Program aborted'
endif

    else

       call h5screate_simple_f( 1,  INT(SHAPE(data), HSIZE_T),  dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1144 <core_hgroup:write_attr_i_i4_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5screate_simple_f'
   stop 'Program aborted'
endif

       call h5pcreate_f( H5P_ATTRIBUTE_CREATE_F,  plist_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1144 <core_hgroup:write_attr_i_i4_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pcreate_f'
   stop 'Program aborted'
endif

    if(.NOT. (.NOT. PRESENT(comp_level))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''.NOT. PRESENT(comp_level)'' failed at line 1144 <core_hgroup:write_attr_i_i4_1_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Attributes cannot be compressed'
      stop
    endif

       mem_type_id = integer_mem_type_(I4)
       file_type_id = integer_file_type_(I4)

       data_ptr = C_LOC(data)

       call h5acreate_f( hg%group_id,  item_name,  file_type_id,  dspace_id,  item_id, hdf_err,  acpl_id=plist_id)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1144 <core_hgroup:write_attr_i_i4_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5acreate_f'
   stop 'Program aborted'
endif

       call h5awrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1144 <core_hgroup:write_attr_i_i4_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5awrite_f'
   stop 'Program aborted'
endif

       call h5pclose_f( plist_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1144 <core_hgroup:write_attr_i_i4_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pclose_f'
   stop 'Program aborted'
endif

       call h5aclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1144 <core_hgroup:write_attr_i_i4_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aclose_f'
   stop 'Program aborted'
endif

       call h5sclose_f( dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1144 <core_hgroup:write_attr_i_i4_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5sclose_f'
   stop 'Program aborted'
endif

    endif

    ! Finish

    return

  end subroutine write_attr_i_i4_1_

  subroutine write_attr_i_i4_2_ (hg, item_name, data, overwrite, comp_level)

    type(hgroup_t), intent(inout)              :: hg
    character(*), intent(in)                   :: item_name
    integer(I4), target, intent(in) :: data(:,:)

    contiguous :: data

    logical, intent(in), optional              :: overwrite
    integer, intent(in), optional              :: comp_level

    logical          :: overwrite_
    integer          :: hdf_err
    logical          :: item_exists
    integer(HSIZE_T) :: item_shape(2)
    integer(HID_T)   :: dspace_id
    integer(HID_T)   :: plist_id
    integer(HID_T)   :: mem_type_id
    integer(HID_T)   :: file_type_id
    type(C_PTR)      :: data_ptr
    integer(HID_T)   :: item_id

    if(PRESENT(overwrite)) then
       overwrite_ = overwrite
    else
       overwrite_ = .FALSE.
    endif

    ! Write the item

    item_exists = attr_exists_(hg, item_name)
    overwrite_ = overwrite_ .AND. item_exists

    if(overwrite_) then

       call get_attr_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 1145 <core_hgroup:write_attr_i_i4' &
 & // '_2_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

       mem_type_id = integer_mem_type_(I4)

       data_ptr = C_LOC(data)

       call h5aopen_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1145 <core_hgroup:write_attr_i_i4_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aopen_f'
   stop 'Program aborted'
endif

       call h5awrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1145 <core_hgroup:write_attr_i_i4_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5awrite_f'
   stop 'Program aborted'
endif

       call h5aclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1145 <core_hgroup:write_attr_i_i4_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aclose_f'
   stop 'Program aborted'
endif

    else

       call h5screate_simple_f( 2,  INT(SHAPE(data), HSIZE_T),  dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1145 <core_hgroup:write_attr_i_i4_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5screate_simple_f'
   stop 'Program aborted'
endif

       call h5pcreate_f( H5P_ATTRIBUTE_CREATE_F,  plist_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1145 <core_hgroup:write_attr_i_i4_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pcreate_f'
   stop 'Program aborted'
endif

    if(.NOT. (.NOT. PRESENT(comp_level))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''.NOT. PRESENT(comp_level)'' failed at line 1145 <core_hgroup:write_attr_i_i4_2_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Attributes cannot be compressed'
      stop
    endif

       mem_type_id = integer_mem_type_(I4)
       file_type_id = integer_file_type_(I4)

       data_ptr = C_LOC(data)

       call h5acreate_f( hg%group_id,  item_name,  file_type_id,  dspace_id,  item_id, hdf_err,  acpl_id=plist_id)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1145 <core_hgroup:write_attr_i_i4_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5acreate_f'
   stop 'Program aborted'
endif

       call h5awrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1145 <core_hgroup:write_attr_i_i4_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5awrite_f'
   stop 'Program aborted'
endif

       call h5pclose_f( plist_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1145 <core_hgroup:write_attr_i_i4_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pclose_f'
   stop 'Program aborted'
endif

       call h5aclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1145 <core_hgroup:write_attr_i_i4_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aclose_f'
   stop 'Program aborted'
endif

       call h5sclose_f( dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1145 <core_hgroup:write_attr_i_i4_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5sclose_f'
   stop 'Program aborted'
endif

    endif

    ! Finish

    return

  end subroutine write_attr_i_i4_2_

  subroutine write_attr_i_i4_3_ (hg, item_name, data, overwrite, comp_level)

    type(hgroup_t), intent(inout)              :: hg
    character(*), intent(in)                   :: item_name
    integer(I4), target, intent(in) :: data(:,:,:)

    contiguous :: data

    logical, intent(in), optional              :: overwrite
    integer, intent(in), optional              :: comp_level

    logical          :: overwrite_
    integer          :: hdf_err
    logical          :: item_exists
    integer(HSIZE_T) :: item_shape(3)
    integer(HID_T)   :: dspace_id
    integer(HID_T)   :: plist_id
    integer(HID_T)   :: mem_type_id
    integer(HID_T)   :: file_type_id
    type(C_PTR)      :: data_ptr
    integer(HID_T)   :: item_id

    if(PRESENT(overwrite)) then
       overwrite_ = overwrite
    else
       overwrite_ = .FALSE.
    endif

    ! Write the item

    item_exists = attr_exists_(hg, item_name)
    overwrite_ = overwrite_ .AND. item_exists

    if(overwrite_) then

       call get_attr_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 1146 <core_hgroup:write_attr_i_i4' &
 & // '_3_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

       mem_type_id = integer_mem_type_(I4)

       data_ptr = C_LOC(data)

       call h5aopen_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1146 <core_hgroup:write_attr_i_i4_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aopen_f'
   stop 'Program aborted'
endif

       call h5awrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1146 <core_hgroup:write_attr_i_i4_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5awrite_f'
   stop 'Program aborted'
endif

       call h5aclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1146 <core_hgroup:write_attr_i_i4_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aclose_f'
   stop 'Program aborted'
endif

    else

       call h5screate_simple_f( 3,  INT(SHAPE(data), HSIZE_T),  dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1146 <core_hgroup:write_attr_i_i4_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5screate_simple_f'
   stop 'Program aborted'
endif

       call h5pcreate_f( H5P_ATTRIBUTE_CREATE_F,  plist_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1146 <core_hgroup:write_attr_i_i4_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pcreate_f'
   stop 'Program aborted'
endif

    if(.NOT. (.NOT. PRESENT(comp_level))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''.NOT. PRESENT(comp_level)'' failed at line 1146 <core_hgroup:write_attr_i_i4_3_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Attributes cannot be compressed'
      stop
    endif

       mem_type_id = integer_mem_type_(I4)
       file_type_id = integer_file_type_(I4)

       data_ptr = C_LOC(data)

       call h5acreate_f( hg%group_id,  item_name,  file_type_id,  dspace_id,  item_id, hdf_err,  acpl_id=plist_id)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1146 <core_hgroup:write_attr_i_i4_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5acreate_f'
   stop 'Program aborted'
endif

       call h5awrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1146 <core_hgroup:write_attr_i_i4_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5awrite_f'
   stop 'Program aborted'
endif

       call h5pclose_f( plist_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1146 <core_hgroup:write_attr_i_i4_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pclose_f'
   stop 'Program aborted'
endif

       call h5aclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1146 <core_hgroup:write_attr_i_i4_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aclose_f'
   stop 'Program aborted'
endif

       call h5sclose_f( dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1146 <core_hgroup:write_attr_i_i4_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5sclose_f'
   stop 'Program aborted'
endif

    endif

    ! Finish

    return

  end subroutine write_attr_i_i4_3_

  subroutine write_attr_i_i4_4_ (hg, item_name, data, overwrite, comp_level)

    type(hgroup_t), intent(inout)              :: hg
    character(*), intent(in)                   :: item_name
    integer(I4), target, intent(in) :: data(:,:,:,:)

    contiguous :: data

    logical, intent(in), optional              :: overwrite
    integer, intent(in), optional              :: comp_level

    logical          :: overwrite_
    integer          :: hdf_err
    logical          :: item_exists
    integer(HSIZE_T) :: item_shape(4)
    integer(HID_T)   :: dspace_id
    integer(HID_T)   :: plist_id
    integer(HID_T)   :: mem_type_id
    integer(HID_T)   :: file_type_id
    type(C_PTR)      :: data_ptr
    integer(HID_T)   :: item_id

    if(PRESENT(overwrite)) then
       overwrite_ = overwrite
    else
       overwrite_ = .FALSE.
    endif

    ! Write the item

    item_exists = attr_exists_(hg, item_name)
    overwrite_ = overwrite_ .AND. item_exists

    if(overwrite_) then

       call get_attr_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 1147 <core_hgroup:write_attr_i_i4' &
 & // '_4_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

       mem_type_id = integer_mem_type_(I4)

       data_ptr = C_LOC(data)

       call h5aopen_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1147 <core_hgroup:write_attr_i_i4_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aopen_f'
   stop 'Program aborted'
endif

       call h5awrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1147 <core_hgroup:write_attr_i_i4_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5awrite_f'
   stop 'Program aborted'
endif

       call h5aclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1147 <core_hgroup:write_attr_i_i4_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aclose_f'
   stop 'Program aborted'
endif

    else

       call h5screate_simple_f( 4,  INT(SHAPE(data), HSIZE_T),  dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1147 <core_hgroup:write_attr_i_i4_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5screate_simple_f'
   stop 'Program aborted'
endif

       call h5pcreate_f( H5P_ATTRIBUTE_CREATE_F,  plist_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1147 <core_hgroup:write_attr_i_i4_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pcreate_f'
   stop 'Program aborted'
endif

    if(.NOT. (.NOT. PRESENT(comp_level))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''.NOT. PRESENT(comp_level)'' failed at line 1147 <core_hgroup:write_attr_i_i4_4_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Attributes cannot be compressed'
      stop
    endif

       mem_type_id = integer_mem_type_(I4)
       file_type_id = integer_file_type_(I4)

       data_ptr = C_LOC(data)

       call h5acreate_f( hg%group_id,  item_name,  file_type_id,  dspace_id,  item_id, hdf_err,  acpl_id=plist_id)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1147 <core_hgroup:write_attr_i_i4_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5acreate_f'
   stop 'Program aborted'
endif

       call h5awrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1147 <core_hgroup:write_attr_i_i4_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5awrite_f'
   stop 'Program aborted'
endif

       call h5pclose_f( plist_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1147 <core_hgroup:write_attr_i_i4_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pclose_f'
   stop 'Program aborted'
endif

       call h5aclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1147 <core_hgroup:write_attr_i_i4_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aclose_f'
   stop 'Program aborted'
endif

       call h5sclose_f( dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1147 <core_hgroup:write_attr_i_i4_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5sclose_f'
   stop 'Program aborted'
endif

    endif

    ! Finish

    return

  end subroutine write_attr_i_i4_4_

  subroutine write_attr_i_i8_0_ (hg, item_name, data, overwrite, comp_level)

    type(hgroup_t), intent(inout)              :: hg
    character(*), intent(in)                   :: item_name
    integer(I8), target, intent(in) :: data

    logical, intent(in), optional              :: overwrite
    integer, intent(in), optional              :: comp_level

    logical          :: overwrite_
    integer          :: hdf_err
    logical          :: item_exists
    integer(HSIZE_T) :: item_shape(0)
    integer(HID_T)   :: dspace_id
    integer(HID_T)   :: plist_id
    integer(HID_T)   :: mem_type_id
    integer(HID_T)   :: file_type_id
    type(C_PTR)      :: data_ptr
    integer(HID_T)   :: item_id

    if(PRESENT(overwrite)) then
       overwrite_ = overwrite
    else
       overwrite_ = .FALSE.
    endif

    ! Write the item

    item_exists = attr_exists_(hg, item_name)
    overwrite_ = overwrite_ .AND. item_exists

    if(overwrite_) then

       call get_attr_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 1149 <core_hgroup:write_attr_i_i8' &
 & // '_0_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

       mem_type_id = integer_mem_type_(I8)

       data_ptr = C_LOC(data)

       call h5aopen_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1149 <core_hgroup:write_attr_i_i8_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aopen_f'
   stop 'Program aborted'
endif

       call h5awrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1149 <core_hgroup:write_attr_i_i8_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5awrite_f'
   stop 'Program aborted'
endif

       call h5aclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1149 <core_hgroup:write_attr_i_i8_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aclose_f'
   stop 'Program aborted'
endif

    else

       call h5screate_f( H5S_SCALAR_F,  dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1149 <core_hgroup:write_attr_i_i8_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5screate_f'
   stop 'Program aborted'
endif

       call h5pcreate_f( H5P_ATTRIBUTE_CREATE_F,  plist_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1149 <core_hgroup:write_attr_i_i8_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pcreate_f'
   stop 'Program aborted'
endif

    if(.NOT. (.NOT. PRESENT(comp_level))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''.NOT. PRESENT(comp_level)'' failed at line 1149 <core_hgroup:write_attr_i_i8_0_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Attributes cannot be compressed'
      stop
    endif

       mem_type_id = integer_mem_type_(I8)
       file_type_id = integer_file_type_(I8)

       data_ptr = C_LOC(data)

       call h5acreate_f( hg%group_id,  item_name,  file_type_id,  dspace_id,  item_id, hdf_err,  acpl_id=plist_id)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1149 <core_hgroup:write_attr_i_i8_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5acreate_f'
   stop 'Program aborted'
endif

       call h5awrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1149 <core_hgroup:write_attr_i_i8_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5awrite_f'
   stop 'Program aborted'
endif

       call h5pclose_f( plist_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1149 <core_hgroup:write_attr_i_i8_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pclose_f'
   stop 'Program aborted'
endif

       call h5aclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1149 <core_hgroup:write_attr_i_i8_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aclose_f'
   stop 'Program aborted'
endif

       call h5sclose_f( dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1149 <core_hgroup:write_attr_i_i8_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5sclose_f'
   stop 'Program aborted'
endif

    endif

    ! Finish

    return

  end subroutine write_attr_i_i8_0_

  subroutine write_attr_i_i8_1_ (hg, item_name, data, overwrite, comp_level)

    type(hgroup_t), intent(inout)              :: hg
    character(*), intent(in)                   :: item_name
    integer(I8), target, intent(in) :: data(:)

    contiguous :: data

    logical, intent(in), optional              :: overwrite
    integer, intent(in), optional              :: comp_level

    logical          :: overwrite_
    integer          :: hdf_err
    logical          :: item_exists
    integer(HSIZE_T) :: item_shape(1)
    integer(HID_T)   :: dspace_id
    integer(HID_T)   :: plist_id
    integer(HID_T)   :: mem_type_id
    integer(HID_T)   :: file_type_id
    type(C_PTR)      :: data_ptr
    integer(HID_T)   :: item_id

    if(PRESENT(overwrite)) then
       overwrite_ = overwrite
    else
       overwrite_ = .FALSE.
    endif

    ! Write the item

    item_exists = attr_exists_(hg, item_name)
    overwrite_ = overwrite_ .AND. item_exists

    if(overwrite_) then

       call get_attr_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 1150 <core_hgroup:write_attr_i_i8' &
 & // '_1_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

       mem_type_id = integer_mem_type_(I8)

       data_ptr = C_LOC(data)

       call h5aopen_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1150 <core_hgroup:write_attr_i_i8_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aopen_f'
   stop 'Program aborted'
endif

       call h5awrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1150 <core_hgroup:write_attr_i_i8_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5awrite_f'
   stop 'Program aborted'
endif

       call h5aclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1150 <core_hgroup:write_attr_i_i8_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aclose_f'
   stop 'Program aborted'
endif

    else

       call h5screate_simple_f( 1,  INT(SHAPE(data), HSIZE_T),  dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1150 <core_hgroup:write_attr_i_i8_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5screate_simple_f'
   stop 'Program aborted'
endif

       call h5pcreate_f( H5P_ATTRIBUTE_CREATE_F,  plist_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1150 <core_hgroup:write_attr_i_i8_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pcreate_f'
   stop 'Program aborted'
endif

    if(.NOT. (.NOT. PRESENT(comp_level))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''.NOT. PRESENT(comp_level)'' failed at line 1150 <core_hgroup:write_attr_i_i8_1_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Attributes cannot be compressed'
      stop
    endif

       mem_type_id = integer_mem_type_(I8)
       file_type_id = integer_file_type_(I8)

       data_ptr = C_LOC(data)

       call h5acreate_f( hg%group_id,  item_name,  file_type_id,  dspace_id,  item_id, hdf_err,  acpl_id=plist_id)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1150 <core_hgroup:write_attr_i_i8_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5acreate_f'
   stop 'Program aborted'
endif

       call h5awrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1150 <core_hgroup:write_attr_i_i8_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5awrite_f'
   stop 'Program aborted'
endif

       call h5pclose_f( plist_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1150 <core_hgroup:write_attr_i_i8_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pclose_f'
   stop 'Program aborted'
endif

       call h5aclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1150 <core_hgroup:write_attr_i_i8_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aclose_f'
   stop 'Program aborted'
endif

       call h5sclose_f( dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1150 <core_hgroup:write_attr_i_i8_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5sclose_f'
   stop 'Program aborted'
endif

    endif

    ! Finish

    return

  end subroutine write_attr_i_i8_1_

  subroutine write_attr_i_i8_2_ (hg, item_name, data, overwrite, comp_level)

    type(hgroup_t), intent(inout)              :: hg
    character(*), intent(in)                   :: item_name
    integer(I8), target, intent(in) :: data(:,:)

    contiguous :: data

    logical, intent(in), optional              :: overwrite
    integer, intent(in), optional              :: comp_level

    logical          :: overwrite_
    integer          :: hdf_err
    logical          :: item_exists
    integer(HSIZE_T) :: item_shape(2)
    integer(HID_T)   :: dspace_id
    integer(HID_T)   :: plist_id
    integer(HID_T)   :: mem_type_id
    integer(HID_T)   :: file_type_id
    type(C_PTR)      :: data_ptr
    integer(HID_T)   :: item_id

    if(PRESENT(overwrite)) then
       overwrite_ = overwrite
    else
       overwrite_ = .FALSE.
    endif

    ! Write the item

    item_exists = attr_exists_(hg, item_name)
    overwrite_ = overwrite_ .AND. item_exists

    if(overwrite_) then

       call get_attr_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 1151 <core_hgroup:write_attr_i_i8' &
 & // '_2_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

       mem_type_id = integer_mem_type_(I8)

       data_ptr = C_LOC(data)

       call h5aopen_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1151 <core_hgroup:write_attr_i_i8_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aopen_f'
   stop 'Program aborted'
endif

       call h5awrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1151 <core_hgroup:write_attr_i_i8_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5awrite_f'
   stop 'Program aborted'
endif

       call h5aclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1151 <core_hgroup:write_attr_i_i8_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aclose_f'
   stop 'Program aborted'
endif

    else

       call h5screate_simple_f( 2,  INT(SHAPE(data), HSIZE_T),  dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1151 <core_hgroup:write_attr_i_i8_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5screate_simple_f'
   stop 'Program aborted'
endif

       call h5pcreate_f( H5P_ATTRIBUTE_CREATE_F,  plist_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1151 <core_hgroup:write_attr_i_i8_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pcreate_f'
   stop 'Program aborted'
endif

    if(.NOT. (.NOT. PRESENT(comp_level))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''.NOT. PRESENT(comp_level)'' failed at line 1151 <core_hgroup:write_attr_i_i8_2_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Attributes cannot be compressed'
      stop
    endif

       mem_type_id = integer_mem_type_(I8)
       file_type_id = integer_file_type_(I8)

       data_ptr = C_LOC(data)

       call h5acreate_f( hg%group_id,  item_name,  file_type_id,  dspace_id,  item_id, hdf_err,  acpl_id=plist_id)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1151 <core_hgroup:write_attr_i_i8_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5acreate_f'
   stop 'Program aborted'
endif

       call h5awrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1151 <core_hgroup:write_attr_i_i8_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5awrite_f'
   stop 'Program aborted'
endif

       call h5pclose_f( plist_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1151 <core_hgroup:write_attr_i_i8_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pclose_f'
   stop 'Program aborted'
endif

       call h5aclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1151 <core_hgroup:write_attr_i_i8_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aclose_f'
   stop 'Program aborted'
endif

       call h5sclose_f( dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1151 <core_hgroup:write_attr_i_i8_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5sclose_f'
   stop 'Program aborted'
endif

    endif

    ! Finish

    return

  end subroutine write_attr_i_i8_2_

  subroutine write_attr_i_i8_3_ (hg, item_name, data, overwrite, comp_level)

    type(hgroup_t), intent(inout)              :: hg
    character(*), intent(in)                   :: item_name
    integer(I8), target, intent(in) :: data(:,:,:)

    contiguous :: data

    logical, intent(in), optional              :: overwrite
    integer, intent(in), optional              :: comp_level

    logical          :: overwrite_
    integer          :: hdf_err
    logical          :: item_exists
    integer(HSIZE_T) :: item_shape(3)
    integer(HID_T)   :: dspace_id
    integer(HID_T)   :: plist_id
    integer(HID_T)   :: mem_type_id
    integer(HID_T)   :: file_type_id
    type(C_PTR)      :: data_ptr
    integer(HID_T)   :: item_id

    if(PRESENT(overwrite)) then
       overwrite_ = overwrite
    else
       overwrite_ = .FALSE.
    endif

    ! Write the item

    item_exists = attr_exists_(hg, item_name)
    overwrite_ = overwrite_ .AND. item_exists

    if(overwrite_) then

       call get_attr_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 1152 <core_hgroup:write_attr_i_i8' &
 & // '_3_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

       mem_type_id = integer_mem_type_(I8)

       data_ptr = C_LOC(data)

       call h5aopen_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1152 <core_hgroup:write_attr_i_i8_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aopen_f'
   stop 'Program aborted'
endif

       call h5awrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1152 <core_hgroup:write_attr_i_i8_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5awrite_f'
   stop 'Program aborted'
endif

       call h5aclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1152 <core_hgroup:write_attr_i_i8_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aclose_f'
   stop 'Program aborted'
endif

    else

       call h5screate_simple_f( 3,  INT(SHAPE(data), HSIZE_T),  dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1152 <core_hgroup:write_attr_i_i8_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5screate_simple_f'
   stop 'Program aborted'
endif

       call h5pcreate_f( H5P_ATTRIBUTE_CREATE_F,  plist_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1152 <core_hgroup:write_attr_i_i8_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pcreate_f'
   stop 'Program aborted'
endif

    if(.NOT. (.NOT. PRESENT(comp_level))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''.NOT. PRESENT(comp_level)'' failed at line 1152 <core_hgroup:write_attr_i_i8_3_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Attributes cannot be compressed'
      stop
    endif

       mem_type_id = integer_mem_type_(I8)
       file_type_id = integer_file_type_(I8)

       data_ptr = C_LOC(data)

       call h5acreate_f( hg%group_id,  item_name,  file_type_id,  dspace_id,  item_id, hdf_err,  acpl_id=plist_id)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1152 <core_hgroup:write_attr_i_i8_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5acreate_f'
   stop 'Program aborted'
endif

       call h5awrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1152 <core_hgroup:write_attr_i_i8_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5awrite_f'
   stop 'Program aborted'
endif

       call h5pclose_f( plist_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1152 <core_hgroup:write_attr_i_i8_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pclose_f'
   stop 'Program aborted'
endif

       call h5aclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1152 <core_hgroup:write_attr_i_i8_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aclose_f'
   stop 'Program aborted'
endif

       call h5sclose_f( dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1152 <core_hgroup:write_attr_i_i8_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5sclose_f'
   stop 'Program aborted'
endif

    endif

    ! Finish

    return

  end subroutine write_attr_i_i8_3_

  subroutine write_attr_i_i8_4_ (hg, item_name, data, overwrite, comp_level)

    type(hgroup_t), intent(inout)              :: hg
    character(*), intent(in)                   :: item_name
    integer(I8), target, intent(in) :: data(:,:,:,:)

    contiguous :: data

    logical, intent(in), optional              :: overwrite
    integer, intent(in), optional              :: comp_level

    logical          :: overwrite_
    integer          :: hdf_err
    logical          :: item_exists
    integer(HSIZE_T) :: item_shape(4)
    integer(HID_T)   :: dspace_id
    integer(HID_T)   :: plist_id
    integer(HID_T)   :: mem_type_id
    integer(HID_T)   :: file_type_id
    type(C_PTR)      :: data_ptr
    integer(HID_T)   :: item_id

    if(PRESENT(overwrite)) then
       overwrite_ = overwrite
    else
       overwrite_ = .FALSE.
    endif

    ! Write the item

    item_exists = attr_exists_(hg, item_name)
    overwrite_ = overwrite_ .AND. item_exists

    if(overwrite_) then

       call get_attr_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 1153 <core_hgroup:write_attr_i_i8' &
 & // '_4_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

       mem_type_id = integer_mem_type_(I8)

       data_ptr = C_LOC(data)

       call h5aopen_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1153 <core_hgroup:write_attr_i_i8_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aopen_f'
   stop 'Program aborted'
endif

       call h5awrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1153 <core_hgroup:write_attr_i_i8_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5awrite_f'
   stop 'Program aborted'
endif

       call h5aclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1153 <core_hgroup:write_attr_i_i8_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aclose_f'
   stop 'Program aborted'
endif

    else

       call h5screate_simple_f( 4,  INT(SHAPE(data), HSIZE_T),  dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1153 <core_hgroup:write_attr_i_i8_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5screate_simple_f'
   stop 'Program aborted'
endif

       call h5pcreate_f( H5P_ATTRIBUTE_CREATE_F,  plist_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1153 <core_hgroup:write_attr_i_i8_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pcreate_f'
   stop 'Program aborted'
endif

    if(.NOT. (.NOT. PRESENT(comp_level))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''.NOT. PRESENT(comp_level)'' failed at line 1153 <core_hgroup:write_attr_i_i8_4_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Attributes cannot be compressed'
      stop
    endif

       mem_type_id = integer_mem_type_(I8)
       file_type_id = integer_file_type_(I8)

       data_ptr = C_LOC(data)

       call h5acreate_f( hg%group_id,  item_name,  file_type_id,  dspace_id,  item_id, hdf_err,  acpl_id=plist_id)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1153 <core_hgroup:write_attr_i_i8_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5acreate_f'
   stop 'Program aborted'
endif

       call h5awrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1153 <core_hgroup:write_attr_i_i8_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5awrite_f'
   stop 'Program aborted'
endif

       call h5pclose_f( plist_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1153 <core_hgroup:write_attr_i_i8_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pclose_f'
   stop 'Program aborted'
endif

       call h5aclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1153 <core_hgroup:write_attr_i_i8_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aclose_f'
   stop 'Program aborted'
endif

       call h5sclose_f( dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1153 <core_hgroup:write_attr_i_i8_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5sclose_f'
   stop 'Program aborted'
endif

    endif

    ! Finish

    return

  end subroutine write_attr_i_i8_4_

  subroutine write_attr_r_sp_0_ (hg, item_name, data, overwrite, comp_level)

    type(hgroup_t), intent(inout)              :: hg
    character(*), intent(in)                   :: item_name
    real(SP), target, intent(in) :: data

    logical, intent(in), optional              :: overwrite
    integer, intent(in), optional              :: comp_level

    logical          :: overwrite_
    integer          :: hdf_err
    logical          :: item_exists
    integer(HSIZE_T) :: item_shape(0)
    integer(HID_T)   :: dspace_id
    integer(HID_T)   :: plist_id
    integer(HID_T)   :: mem_type_id
    integer(HID_T)   :: file_type_id
    type(C_PTR)      :: data_ptr
    integer(HID_T)   :: item_id

    if(PRESENT(overwrite)) then
       overwrite_ = overwrite
    else
       overwrite_ = .FALSE.
    endif

    ! Write the item

    item_exists = attr_exists_(hg, item_name)
    overwrite_ = overwrite_ .AND. item_exists

    if(overwrite_) then

       call get_attr_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 1155 <core_hgroup:write_attr_r_sp' &
 & // '_0_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

       mem_type_id = real_mem_type_(SP)

       data_ptr = C_LOC(data)

       call h5aopen_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1155 <core_hgroup:write_attr_r_sp_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aopen_f'
   stop 'Program aborted'
endif

       call h5awrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1155 <core_hgroup:write_attr_r_sp_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5awrite_f'
   stop 'Program aborted'
endif

       call h5aclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1155 <core_hgroup:write_attr_r_sp_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aclose_f'
   stop 'Program aborted'
endif

    else

       call h5screate_f( H5S_SCALAR_F,  dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1155 <core_hgroup:write_attr_r_sp_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5screate_f'
   stop 'Program aborted'
endif

       call h5pcreate_f( H5P_ATTRIBUTE_CREATE_F,  plist_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1155 <core_hgroup:write_attr_r_sp_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pcreate_f'
   stop 'Program aborted'
endif

    if(.NOT. (.NOT. PRESENT(comp_level))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''.NOT. PRESENT(comp_level)'' failed at line 1155 <core_hgroup:write_attr_r_sp_0_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Attributes cannot be compressed'
      stop
    endif

       mem_type_id = real_mem_type_(SP)
       file_type_id = real_file_type_(SP)

       data_ptr = C_LOC(data)

       call h5acreate_f( hg%group_id,  item_name,  file_type_id,  dspace_id,  item_id, hdf_err,  acpl_id=plist_id)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1155 <core_hgroup:write_attr_r_sp_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5acreate_f'
   stop 'Program aborted'
endif

       call h5awrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1155 <core_hgroup:write_attr_r_sp_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5awrite_f'
   stop 'Program aborted'
endif

       call h5pclose_f( plist_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1155 <core_hgroup:write_attr_r_sp_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pclose_f'
   stop 'Program aborted'
endif

       call h5aclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1155 <core_hgroup:write_attr_r_sp_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aclose_f'
   stop 'Program aborted'
endif

       call h5sclose_f( dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1155 <core_hgroup:write_attr_r_sp_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5sclose_f'
   stop 'Program aborted'
endif

    endif

    ! Finish

    return

  end subroutine write_attr_r_sp_0_

  subroutine write_attr_r_sp_1_ (hg, item_name, data, overwrite, comp_level)

    type(hgroup_t), intent(inout)              :: hg
    character(*), intent(in)                   :: item_name
    real(SP), target, intent(in) :: data(:)

    contiguous :: data

    logical, intent(in), optional              :: overwrite
    integer, intent(in), optional              :: comp_level

    logical          :: overwrite_
    integer          :: hdf_err
    logical          :: item_exists
    integer(HSIZE_T) :: item_shape(1)
    integer(HID_T)   :: dspace_id
    integer(HID_T)   :: plist_id
    integer(HID_T)   :: mem_type_id
    integer(HID_T)   :: file_type_id
    type(C_PTR)      :: data_ptr
    integer(HID_T)   :: item_id

    if(PRESENT(overwrite)) then
       overwrite_ = overwrite
    else
       overwrite_ = .FALSE.
    endif

    ! Write the item

    item_exists = attr_exists_(hg, item_name)
    overwrite_ = overwrite_ .AND. item_exists

    if(overwrite_) then

       call get_attr_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 1156 <core_hgroup:write_attr_r_sp' &
 & // '_1_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

       mem_type_id = real_mem_type_(SP)

       data_ptr = C_LOC(data)

       call h5aopen_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1156 <core_hgroup:write_attr_r_sp_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aopen_f'
   stop 'Program aborted'
endif

       call h5awrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1156 <core_hgroup:write_attr_r_sp_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5awrite_f'
   stop 'Program aborted'
endif

       call h5aclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1156 <core_hgroup:write_attr_r_sp_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aclose_f'
   stop 'Program aborted'
endif

    else

       call h5screate_simple_f( 1,  INT(SHAPE(data), HSIZE_T),  dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1156 <core_hgroup:write_attr_r_sp_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5screate_simple_f'
   stop 'Program aborted'
endif

       call h5pcreate_f( H5P_ATTRIBUTE_CREATE_F,  plist_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1156 <core_hgroup:write_attr_r_sp_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pcreate_f'
   stop 'Program aborted'
endif

    if(.NOT. (.NOT. PRESENT(comp_level))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''.NOT. PRESENT(comp_level)'' failed at line 1156 <core_hgroup:write_attr_r_sp_1_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Attributes cannot be compressed'
      stop
    endif

       mem_type_id = real_mem_type_(SP)
       file_type_id = real_file_type_(SP)

       data_ptr = C_LOC(data)

       call h5acreate_f( hg%group_id,  item_name,  file_type_id,  dspace_id,  item_id, hdf_err,  acpl_id=plist_id)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1156 <core_hgroup:write_attr_r_sp_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5acreate_f'
   stop 'Program aborted'
endif

       call h5awrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1156 <core_hgroup:write_attr_r_sp_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5awrite_f'
   stop 'Program aborted'
endif

       call h5pclose_f( plist_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1156 <core_hgroup:write_attr_r_sp_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pclose_f'
   stop 'Program aborted'
endif

       call h5aclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1156 <core_hgroup:write_attr_r_sp_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aclose_f'
   stop 'Program aborted'
endif

       call h5sclose_f( dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1156 <core_hgroup:write_attr_r_sp_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5sclose_f'
   stop 'Program aborted'
endif

    endif

    ! Finish

    return

  end subroutine write_attr_r_sp_1_

  subroutine write_attr_r_sp_2_ (hg, item_name, data, overwrite, comp_level)

    type(hgroup_t), intent(inout)              :: hg
    character(*), intent(in)                   :: item_name
    real(SP), target, intent(in) :: data(:,:)

    contiguous :: data

    logical, intent(in), optional              :: overwrite
    integer, intent(in), optional              :: comp_level

    logical          :: overwrite_
    integer          :: hdf_err
    logical          :: item_exists
    integer(HSIZE_T) :: item_shape(2)
    integer(HID_T)   :: dspace_id
    integer(HID_T)   :: plist_id
    integer(HID_T)   :: mem_type_id
    integer(HID_T)   :: file_type_id
    type(C_PTR)      :: data_ptr
    integer(HID_T)   :: item_id

    if(PRESENT(overwrite)) then
       overwrite_ = overwrite
    else
       overwrite_ = .FALSE.
    endif

    ! Write the item

    item_exists = attr_exists_(hg, item_name)
    overwrite_ = overwrite_ .AND. item_exists

    if(overwrite_) then

       call get_attr_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 1157 <core_hgroup:write_attr_r_sp' &
 & // '_2_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

       mem_type_id = real_mem_type_(SP)

       data_ptr = C_LOC(data)

       call h5aopen_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1157 <core_hgroup:write_attr_r_sp_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aopen_f'
   stop 'Program aborted'
endif

       call h5awrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1157 <core_hgroup:write_attr_r_sp_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5awrite_f'
   stop 'Program aborted'
endif

       call h5aclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1157 <core_hgroup:write_attr_r_sp_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aclose_f'
   stop 'Program aborted'
endif

    else

       call h5screate_simple_f( 2,  INT(SHAPE(data), HSIZE_T),  dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1157 <core_hgroup:write_attr_r_sp_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5screate_simple_f'
   stop 'Program aborted'
endif

       call h5pcreate_f( H5P_ATTRIBUTE_CREATE_F,  plist_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1157 <core_hgroup:write_attr_r_sp_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pcreate_f'
   stop 'Program aborted'
endif

    if(.NOT. (.NOT. PRESENT(comp_level))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''.NOT. PRESENT(comp_level)'' failed at line 1157 <core_hgroup:write_attr_r_sp_2_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Attributes cannot be compressed'
      stop
    endif

       mem_type_id = real_mem_type_(SP)
       file_type_id = real_file_type_(SP)

       data_ptr = C_LOC(data)

       call h5acreate_f( hg%group_id,  item_name,  file_type_id,  dspace_id,  item_id, hdf_err,  acpl_id=plist_id)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1157 <core_hgroup:write_attr_r_sp_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5acreate_f'
   stop 'Program aborted'
endif

       call h5awrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1157 <core_hgroup:write_attr_r_sp_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5awrite_f'
   stop 'Program aborted'
endif

       call h5pclose_f( plist_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1157 <core_hgroup:write_attr_r_sp_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pclose_f'
   stop 'Program aborted'
endif

       call h5aclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1157 <core_hgroup:write_attr_r_sp_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aclose_f'
   stop 'Program aborted'
endif

       call h5sclose_f( dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1157 <core_hgroup:write_attr_r_sp_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5sclose_f'
   stop 'Program aborted'
endif

    endif

    ! Finish

    return

  end subroutine write_attr_r_sp_2_

  subroutine write_attr_r_sp_3_ (hg, item_name, data, overwrite, comp_level)

    type(hgroup_t), intent(inout)              :: hg
    character(*), intent(in)                   :: item_name
    real(SP), target, intent(in) :: data(:,:,:)

    contiguous :: data

    logical, intent(in), optional              :: overwrite
    integer, intent(in), optional              :: comp_level

    logical          :: overwrite_
    integer          :: hdf_err
    logical          :: item_exists
    integer(HSIZE_T) :: item_shape(3)
    integer(HID_T)   :: dspace_id
    integer(HID_T)   :: plist_id
    integer(HID_T)   :: mem_type_id
    integer(HID_T)   :: file_type_id
    type(C_PTR)      :: data_ptr
    integer(HID_T)   :: item_id

    if(PRESENT(overwrite)) then
       overwrite_ = overwrite
    else
       overwrite_ = .FALSE.
    endif

    ! Write the item

    item_exists = attr_exists_(hg, item_name)
    overwrite_ = overwrite_ .AND. item_exists

    if(overwrite_) then

       call get_attr_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 1158 <core_hgroup:write_attr_r_sp' &
 & // '_3_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

       mem_type_id = real_mem_type_(SP)

       data_ptr = C_LOC(data)

       call h5aopen_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1158 <core_hgroup:write_attr_r_sp_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aopen_f'
   stop 'Program aborted'
endif

       call h5awrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1158 <core_hgroup:write_attr_r_sp_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5awrite_f'
   stop 'Program aborted'
endif

       call h5aclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1158 <core_hgroup:write_attr_r_sp_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aclose_f'
   stop 'Program aborted'
endif

    else

       call h5screate_simple_f( 3,  INT(SHAPE(data), HSIZE_T),  dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1158 <core_hgroup:write_attr_r_sp_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5screate_simple_f'
   stop 'Program aborted'
endif

       call h5pcreate_f( H5P_ATTRIBUTE_CREATE_F,  plist_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1158 <core_hgroup:write_attr_r_sp_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pcreate_f'
   stop 'Program aborted'
endif

    if(.NOT. (.NOT. PRESENT(comp_level))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''.NOT. PRESENT(comp_level)'' failed at line 1158 <core_hgroup:write_attr_r_sp_3_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Attributes cannot be compressed'
      stop
    endif

       mem_type_id = real_mem_type_(SP)
       file_type_id = real_file_type_(SP)

       data_ptr = C_LOC(data)

       call h5acreate_f( hg%group_id,  item_name,  file_type_id,  dspace_id,  item_id, hdf_err,  acpl_id=plist_id)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1158 <core_hgroup:write_attr_r_sp_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5acreate_f'
   stop 'Program aborted'
endif

       call h5awrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1158 <core_hgroup:write_attr_r_sp_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5awrite_f'
   stop 'Program aborted'
endif

       call h5pclose_f( plist_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1158 <core_hgroup:write_attr_r_sp_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pclose_f'
   stop 'Program aborted'
endif

       call h5aclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1158 <core_hgroup:write_attr_r_sp_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aclose_f'
   stop 'Program aborted'
endif

       call h5sclose_f( dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1158 <core_hgroup:write_attr_r_sp_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5sclose_f'
   stop 'Program aborted'
endif

    endif

    ! Finish

    return

  end subroutine write_attr_r_sp_3_

  subroutine write_attr_r_sp_4_ (hg, item_name, data, overwrite, comp_level)

    type(hgroup_t), intent(inout)              :: hg
    character(*), intent(in)                   :: item_name
    real(SP), target, intent(in) :: data(:,:,:,:)

    contiguous :: data

    logical, intent(in), optional              :: overwrite
    integer, intent(in), optional              :: comp_level

    logical          :: overwrite_
    integer          :: hdf_err
    logical          :: item_exists
    integer(HSIZE_T) :: item_shape(4)
    integer(HID_T)   :: dspace_id
    integer(HID_T)   :: plist_id
    integer(HID_T)   :: mem_type_id
    integer(HID_T)   :: file_type_id
    type(C_PTR)      :: data_ptr
    integer(HID_T)   :: item_id

    if(PRESENT(overwrite)) then
       overwrite_ = overwrite
    else
       overwrite_ = .FALSE.
    endif

    ! Write the item

    item_exists = attr_exists_(hg, item_name)
    overwrite_ = overwrite_ .AND. item_exists

    if(overwrite_) then

       call get_attr_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 1159 <core_hgroup:write_attr_r_sp' &
 & // '_4_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

       mem_type_id = real_mem_type_(SP)

       data_ptr = C_LOC(data)

       call h5aopen_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1159 <core_hgroup:write_attr_r_sp_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aopen_f'
   stop 'Program aborted'
endif

       call h5awrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1159 <core_hgroup:write_attr_r_sp_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5awrite_f'
   stop 'Program aborted'
endif

       call h5aclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1159 <core_hgroup:write_attr_r_sp_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aclose_f'
   stop 'Program aborted'
endif

    else

       call h5screate_simple_f( 4,  INT(SHAPE(data), HSIZE_T),  dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1159 <core_hgroup:write_attr_r_sp_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5screate_simple_f'
   stop 'Program aborted'
endif

       call h5pcreate_f( H5P_ATTRIBUTE_CREATE_F,  plist_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1159 <core_hgroup:write_attr_r_sp_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pcreate_f'
   stop 'Program aborted'
endif

    if(.NOT. (.NOT. PRESENT(comp_level))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''.NOT. PRESENT(comp_level)'' failed at line 1159 <core_hgroup:write_attr_r_sp_4_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Attributes cannot be compressed'
      stop
    endif

       mem_type_id = real_mem_type_(SP)
       file_type_id = real_file_type_(SP)

       data_ptr = C_LOC(data)

       call h5acreate_f( hg%group_id,  item_name,  file_type_id,  dspace_id,  item_id, hdf_err,  acpl_id=plist_id)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1159 <core_hgroup:write_attr_r_sp_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5acreate_f'
   stop 'Program aborted'
endif

       call h5awrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1159 <core_hgroup:write_attr_r_sp_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5awrite_f'
   stop 'Program aborted'
endif

       call h5pclose_f( plist_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1159 <core_hgroup:write_attr_r_sp_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pclose_f'
   stop 'Program aborted'
endif

       call h5aclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1159 <core_hgroup:write_attr_r_sp_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aclose_f'
   stop 'Program aborted'
endif

       call h5sclose_f( dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1159 <core_hgroup:write_attr_r_sp_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5sclose_f'
   stop 'Program aborted'
endif

    endif

    ! Finish

    return

  end subroutine write_attr_r_sp_4_

  subroutine write_attr_r_dp_0_ (hg, item_name, data, overwrite, comp_level)

    type(hgroup_t), intent(inout)              :: hg
    character(*), intent(in)                   :: item_name
    real(DP), target, intent(in) :: data

    logical, intent(in), optional              :: overwrite
    integer, intent(in), optional              :: comp_level

    logical          :: overwrite_
    integer          :: hdf_err
    logical          :: item_exists
    integer(HSIZE_T) :: item_shape(0)
    integer(HID_T)   :: dspace_id
    integer(HID_T)   :: plist_id
    integer(HID_T)   :: mem_type_id
    integer(HID_T)   :: file_type_id
    type(C_PTR)      :: data_ptr
    integer(HID_T)   :: item_id

    if(PRESENT(overwrite)) then
       overwrite_ = overwrite
    else
       overwrite_ = .FALSE.
    endif

    ! Write the item

    item_exists = attr_exists_(hg, item_name)
    overwrite_ = overwrite_ .AND. item_exists

    if(overwrite_) then

       call get_attr_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 1161 <core_hgroup:write_attr_r_dp' &
 & // '_0_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

       mem_type_id = real_mem_type_(DP)

       data_ptr = C_LOC(data)

       call h5aopen_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1161 <core_hgroup:write_attr_r_dp_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aopen_f'
   stop 'Program aborted'
endif

       call h5awrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1161 <core_hgroup:write_attr_r_dp_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5awrite_f'
   stop 'Program aborted'
endif

       call h5aclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1161 <core_hgroup:write_attr_r_dp_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aclose_f'
   stop 'Program aborted'
endif

    else

       call h5screate_f( H5S_SCALAR_F,  dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1161 <core_hgroup:write_attr_r_dp_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5screate_f'
   stop 'Program aborted'
endif

       call h5pcreate_f( H5P_ATTRIBUTE_CREATE_F,  plist_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1161 <core_hgroup:write_attr_r_dp_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pcreate_f'
   stop 'Program aborted'
endif

    if(.NOT. (.NOT. PRESENT(comp_level))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''.NOT. PRESENT(comp_level)'' failed at line 1161 <core_hgroup:write_attr_r_dp_0_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Attributes cannot be compressed'
      stop
    endif

       mem_type_id = real_mem_type_(DP)
       file_type_id = real_file_type_(DP)

       data_ptr = C_LOC(data)

       call h5acreate_f( hg%group_id,  item_name,  file_type_id,  dspace_id,  item_id, hdf_err,  acpl_id=plist_id)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1161 <core_hgroup:write_attr_r_dp_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5acreate_f'
   stop 'Program aborted'
endif

       call h5awrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1161 <core_hgroup:write_attr_r_dp_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5awrite_f'
   stop 'Program aborted'
endif

       call h5pclose_f( plist_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1161 <core_hgroup:write_attr_r_dp_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pclose_f'
   stop 'Program aborted'
endif

       call h5aclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1161 <core_hgroup:write_attr_r_dp_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aclose_f'
   stop 'Program aborted'
endif

       call h5sclose_f( dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1161 <core_hgroup:write_attr_r_dp_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5sclose_f'
   stop 'Program aborted'
endif

    endif

    ! Finish

    return

  end subroutine write_attr_r_dp_0_

  subroutine write_attr_r_dp_1_ (hg, item_name, data, overwrite, comp_level)

    type(hgroup_t), intent(inout)              :: hg
    character(*), intent(in)                   :: item_name
    real(DP), target, intent(in) :: data(:)

    contiguous :: data

    logical, intent(in), optional              :: overwrite
    integer, intent(in), optional              :: comp_level

    logical          :: overwrite_
    integer          :: hdf_err
    logical          :: item_exists
    integer(HSIZE_T) :: item_shape(1)
    integer(HID_T)   :: dspace_id
    integer(HID_T)   :: plist_id
    integer(HID_T)   :: mem_type_id
    integer(HID_T)   :: file_type_id
    type(C_PTR)      :: data_ptr
    integer(HID_T)   :: item_id

    if(PRESENT(overwrite)) then
       overwrite_ = overwrite
    else
       overwrite_ = .FALSE.
    endif

    ! Write the item

    item_exists = attr_exists_(hg, item_name)
    overwrite_ = overwrite_ .AND. item_exists

    if(overwrite_) then

       call get_attr_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 1162 <core_hgroup:write_attr_r_dp' &
 & // '_1_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

       mem_type_id = real_mem_type_(DP)

       data_ptr = C_LOC(data)

       call h5aopen_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1162 <core_hgroup:write_attr_r_dp_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aopen_f'
   stop 'Program aborted'
endif

       call h5awrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1162 <core_hgroup:write_attr_r_dp_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5awrite_f'
   stop 'Program aborted'
endif

       call h5aclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1162 <core_hgroup:write_attr_r_dp_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aclose_f'
   stop 'Program aborted'
endif

    else

       call h5screate_simple_f( 1,  INT(SHAPE(data), HSIZE_T),  dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1162 <core_hgroup:write_attr_r_dp_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5screate_simple_f'
   stop 'Program aborted'
endif

       call h5pcreate_f( H5P_ATTRIBUTE_CREATE_F,  plist_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1162 <core_hgroup:write_attr_r_dp_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pcreate_f'
   stop 'Program aborted'
endif

    if(.NOT. (.NOT. PRESENT(comp_level))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''.NOT. PRESENT(comp_level)'' failed at line 1162 <core_hgroup:write_attr_r_dp_1_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Attributes cannot be compressed'
      stop
    endif

       mem_type_id = real_mem_type_(DP)
       file_type_id = real_file_type_(DP)

       data_ptr = C_LOC(data)

       call h5acreate_f( hg%group_id,  item_name,  file_type_id,  dspace_id,  item_id, hdf_err,  acpl_id=plist_id)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1162 <core_hgroup:write_attr_r_dp_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5acreate_f'
   stop 'Program aborted'
endif

       call h5awrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1162 <core_hgroup:write_attr_r_dp_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5awrite_f'
   stop 'Program aborted'
endif

       call h5pclose_f( plist_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1162 <core_hgroup:write_attr_r_dp_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pclose_f'
   stop 'Program aborted'
endif

       call h5aclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1162 <core_hgroup:write_attr_r_dp_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aclose_f'
   stop 'Program aborted'
endif

       call h5sclose_f( dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1162 <core_hgroup:write_attr_r_dp_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5sclose_f'
   stop 'Program aborted'
endif

    endif

    ! Finish

    return

  end subroutine write_attr_r_dp_1_

  subroutine write_attr_r_dp_2_ (hg, item_name, data, overwrite, comp_level)

    type(hgroup_t), intent(inout)              :: hg
    character(*), intent(in)                   :: item_name
    real(DP), target, intent(in) :: data(:,:)

    contiguous :: data

    logical, intent(in), optional              :: overwrite
    integer, intent(in), optional              :: comp_level

    logical          :: overwrite_
    integer          :: hdf_err
    logical          :: item_exists
    integer(HSIZE_T) :: item_shape(2)
    integer(HID_T)   :: dspace_id
    integer(HID_T)   :: plist_id
    integer(HID_T)   :: mem_type_id
    integer(HID_T)   :: file_type_id
    type(C_PTR)      :: data_ptr
    integer(HID_T)   :: item_id

    if(PRESENT(overwrite)) then
       overwrite_ = overwrite
    else
       overwrite_ = .FALSE.
    endif

    ! Write the item

    item_exists = attr_exists_(hg, item_name)
    overwrite_ = overwrite_ .AND. item_exists

    if(overwrite_) then

       call get_attr_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 1163 <core_hgroup:write_attr_r_dp' &
 & // '_2_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

       mem_type_id = real_mem_type_(DP)

       data_ptr = C_LOC(data)

       call h5aopen_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1163 <core_hgroup:write_attr_r_dp_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aopen_f'
   stop 'Program aborted'
endif

       call h5awrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1163 <core_hgroup:write_attr_r_dp_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5awrite_f'
   stop 'Program aborted'
endif

       call h5aclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1163 <core_hgroup:write_attr_r_dp_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aclose_f'
   stop 'Program aborted'
endif

    else

       call h5screate_simple_f( 2,  INT(SHAPE(data), HSIZE_T),  dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1163 <core_hgroup:write_attr_r_dp_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5screate_simple_f'
   stop 'Program aborted'
endif

       call h5pcreate_f( H5P_ATTRIBUTE_CREATE_F,  plist_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1163 <core_hgroup:write_attr_r_dp_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pcreate_f'
   stop 'Program aborted'
endif

    if(.NOT. (.NOT. PRESENT(comp_level))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''.NOT. PRESENT(comp_level)'' failed at line 1163 <core_hgroup:write_attr_r_dp_2_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Attributes cannot be compressed'
      stop
    endif

       mem_type_id = real_mem_type_(DP)
       file_type_id = real_file_type_(DP)

       data_ptr = C_LOC(data)

       call h5acreate_f( hg%group_id,  item_name,  file_type_id,  dspace_id,  item_id, hdf_err,  acpl_id=plist_id)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1163 <core_hgroup:write_attr_r_dp_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5acreate_f'
   stop 'Program aborted'
endif

       call h5awrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1163 <core_hgroup:write_attr_r_dp_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5awrite_f'
   stop 'Program aborted'
endif

       call h5pclose_f( plist_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1163 <core_hgroup:write_attr_r_dp_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pclose_f'
   stop 'Program aborted'
endif

       call h5aclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1163 <core_hgroup:write_attr_r_dp_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aclose_f'
   stop 'Program aborted'
endif

       call h5sclose_f( dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1163 <core_hgroup:write_attr_r_dp_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5sclose_f'
   stop 'Program aborted'
endif

    endif

    ! Finish

    return

  end subroutine write_attr_r_dp_2_

  subroutine write_attr_r_dp_3_ (hg, item_name, data, overwrite, comp_level)

    type(hgroup_t), intent(inout)              :: hg
    character(*), intent(in)                   :: item_name
    real(DP), target, intent(in) :: data(:,:,:)

    contiguous :: data

    logical, intent(in), optional              :: overwrite
    integer, intent(in), optional              :: comp_level

    logical          :: overwrite_
    integer          :: hdf_err
    logical          :: item_exists
    integer(HSIZE_T) :: item_shape(3)
    integer(HID_T)   :: dspace_id
    integer(HID_T)   :: plist_id
    integer(HID_T)   :: mem_type_id
    integer(HID_T)   :: file_type_id
    type(C_PTR)      :: data_ptr
    integer(HID_T)   :: item_id

    if(PRESENT(overwrite)) then
       overwrite_ = overwrite
    else
       overwrite_ = .FALSE.
    endif

    ! Write the item

    item_exists = attr_exists_(hg, item_name)
    overwrite_ = overwrite_ .AND. item_exists

    if(overwrite_) then

       call get_attr_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 1164 <core_hgroup:write_attr_r_dp' &
 & // '_3_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

       mem_type_id = real_mem_type_(DP)

       data_ptr = C_LOC(data)

       call h5aopen_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1164 <core_hgroup:write_attr_r_dp_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aopen_f'
   stop 'Program aborted'
endif

       call h5awrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1164 <core_hgroup:write_attr_r_dp_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5awrite_f'
   stop 'Program aborted'
endif

       call h5aclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1164 <core_hgroup:write_attr_r_dp_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aclose_f'
   stop 'Program aborted'
endif

    else

       call h5screate_simple_f( 3,  INT(SHAPE(data), HSIZE_T),  dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1164 <core_hgroup:write_attr_r_dp_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5screate_simple_f'
   stop 'Program aborted'
endif

       call h5pcreate_f( H5P_ATTRIBUTE_CREATE_F,  plist_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1164 <core_hgroup:write_attr_r_dp_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pcreate_f'
   stop 'Program aborted'
endif

    if(.NOT. (.NOT. PRESENT(comp_level))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''.NOT. PRESENT(comp_level)'' failed at line 1164 <core_hgroup:write_attr_r_dp_3_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Attributes cannot be compressed'
      stop
    endif

       mem_type_id = real_mem_type_(DP)
       file_type_id = real_file_type_(DP)

       data_ptr = C_LOC(data)

       call h5acreate_f( hg%group_id,  item_name,  file_type_id,  dspace_id,  item_id, hdf_err,  acpl_id=plist_id)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1164 <core_hgroup:write_attr_r_dp_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5acreate_f'
   stop 'Program aborted'
endif

       call h5awrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1164 <core_hgroup:write_attr_r_dp_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5awrite_f'
   stop 'Program aborted'
endif

       call h5pclose_f( plist_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1164 <core_hgroup:write_attr_r_dp_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pclose_f'
   stop 'Program aborted'
endif

       call h5aclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1164 <core_hgroup:write_attr_r_dp_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aclose_f'
   stop 'Program aborted'
endif

       call h5sclose_f( dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1164 <core_hgroup:write_attr_r_dp_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5sclose_f'
   stop 'Program aborted'
endif

    endif

    ! Finish

    return

  end subroutine write_attr_r_dp_3_

  subroutine write_attr_r_dp_4_ (hg, item_name, data, overwrite, comp_level)

    type(hgroup_t), intent(inout)              :: hg
    character(*), intent(in)                   :: item_name
    real(DP), target, intent(in) :: data(:,:,:,:)

    contiguous :: data

    logical, intent(in), optional              :: overwrite
    integer, intent(in), optional              :: comp_level

    logical          :: overwrite_
    integer          :: hdf_err
    logical          :: item_exists
    integer(HSIZE_T) :: item_shape(4)
    integer(HID_T)   :: dspace_id
    integer(HID_T)   :: plist_id
    integer(HID_T)   :: mem_type_id
    integer(HID_T)   :: file_type_id
    type(C_PTR)      :: data_ptr
    integer(HID_T)   :: item_id

    if(PRESENT(overwrite)) then
       overwrite_ = overwrite
    else
       overwrite_ = .FALSE.
    endif

    ! Write the item

    item_exists = attr_exists_(hg, item_name)
    overwrite_ = overwrite_ .AND. item_exists

    if(overwrite_) then

       call get_attr_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 1165 <core_hgroup:write_attr_r_dp' &
 & // '_4_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

       mem_type_id = real_mem_type_(DP)

       data_ptr = C_LOC(data)

       call h5aopen_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1165 <core_hgroup:write_attr_r_dp_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aopen_f'
   stop 'Program aborted'
endif

       call h5awrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1165 <core_hgroup:write_attr_r_dp_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5awrite_f'
   stop 'Program aborted'
endif

       call h5aclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1165 <core_hgroup:write_attr_r_dp_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aclose_f'
   stop 'Program aborted'
endif

    else

       call h5screate_simple_f( 4,  INT(SHAPE(data), HSIZE_T),  dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1165 <core_hgroup:write_attr_r_dp_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5screate_simple_f'
   stop 'Program aborted'
endif

       call h5pcreate_f( H5P_ATTRIBUTE_CREATE_F,  plist_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1165 <core_hgroup:write_attr_r_dp_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pcreate_f'
   stop 'Program aborted'
endif

    if(.NOT. (.NOT. PRESENT(comp_level))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''.NOT. PRESENT(comp_level)'' failed at line 1165 <core_hgroup:write_attr_r_dp_4_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Attributes cannot be compressed'
      stop
    endif

       mem_type_id = real_mem_type_(DP)
       file_type_id = real_file_type_(DP)

       data_ptr = C_LOC(data)

       call h5acreate_f( hg%group_id,  item_name,  file_type_id,  dspace_id,  item_id, hdf_err,  acpl_id=plist_id)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1165 <core_hgroup:write_attr_r_dp_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5acreate_f'
   stop 'Program aborted'
endif

       call h5awrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1165 <core_hgroup:write_attr_r_dp_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5awrite_f'
   stop 'Program aborted'
endif

       call h5pclose_f( plist_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1165 <core_hgroup:write_attr_r_dp_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pclose_f'
   stop 'Program aborted'
endif

       call h5aclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1165 <core_hgroup:write_attr_r_dp_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aclose_f'
   stop 'Program aborted'
endif

       call h5sclose_f( dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1165 <core_hgroup:write_attr_r_dp_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5sclose_f'
   stop 'Program aborted'
endif

    endif

    ! Finish

    return

  end subroutine write_attr_r_dp_4_

  subroutine write_attr_c_sp_0_ (hg, item_name, data, overwrite, comp_level)

    type(hgroup_t), intent(inout)              :: hg
    character(*), intent(in)                   :: item_name
    complex(SP), target, intent(in) :: data

    logical, intent(in), optional              :: overwrite
    integer, intent(in), optional              :: comp_level

    logical          :: overwrite_
    integer          :: hdf_err
    logical          :: item_exists
    integer(HSIZE_T) :: item_shape(0)
    integer(HID_T)   :: dspace_id
    integer(HID_T)   :: plist_id
    integer(HID_T)   :: mem_type_id
    integer(HID_T)   :: file_type_id
    type(C_PTR)      :: data_ptr
    integer(HID_T)   :: item_id

    if(PRESENT(overwrite)) then
       overwrite_ = overwrite
    else
       overwrite_ = .FALSE.
    endif

    ! Write the item

    item_exists = attr_exists_(hg, item_name)
    overwrite_ = overwrite_ .AND. item_exists

    if(overwrite_) then

       call get_attr_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 1167 <core_hgroup:write_attr_c_sp' &
 & // '_0_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

       mem_type_id = complex_mem_type_(SP)

       data_ptr = C_LOC(data)

       call h5aopen_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1167 <core_hgroup:write_attr_c_sp_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aopen_f'
   stop 'Program aborted'
endif

       call h5awrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1167 <core_hgroup:write_attr_c_sp_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5awrite_f'
   stop 'Program aborted'
endif

       call h5aclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1167 <core_hgroup:write_attr_c_sp_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aclose_f'
   stop 'Program aborted'
endif

    else

       call h5screate_f( H5S_SCALAR_F,  dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1167 <core_hgroup:write_attr_c_sp_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5screate_f'
   stop 'Program aborted'
endif

       call h5pcreate_f( H5P_ATTRIBUTE_CREATE_F,  plist_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1167 <core_hgroup:write_attr_c_sp_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pcreate_f'
   stop 'Program aborted'
endif

    if(.NOT. (.NOT. PRESENT(comp_level))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''.NOT. PRESENT(comp_level)'' failed at line 1167 <core_hgroup:write_attr_c_sp_0_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Attributes cannot be compressed'
      stop
    endif

       mem_type_id = complex_mem_type_(SP)
       file_type_id = complex_file_type_(SP)

       data_ptr = C_LOC(data)

       call h5acreate_f( hg%group_id,  item_name,  file_type_id,  dspace_id,  item_id, hdf_err,  acpl_id=plist_id)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1167 <core_hgroup:write_attr_c_sp_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5acreate_f'
   stop 'Program aborted'
endif

       call h5awrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1167 <core_hgroup:write_attr_c_sp_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5awrite_f'
   stop 'Program aborted'
endif

       call h5pclose_f( plist_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1167 <core_hgroup:write_attr_c_sp_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pclose_f'
   stop 'Program aborted'
endif

       call h5aclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1167 <core_hgroup:write_attr_c_sp_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aclose_f'
   stop 'Program aborted'
endif

       call h5sclose_f( dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1167 <core_hgroup:write_attr_c_sp_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5sclose_f'
   stop 'Program aborted'
endif

    endif

    ! Finish

    return

  end subroutine write_attr_c_sp_0_

  subroutine write_attr_c_sp_1_ (hg, item_name, data, overwrite, comp_level)

    type(hgroup_t), intent(inout)              :: hg
    character(*), intent(in)                   :: item_name
    complex(SP), target, intent(in) :: data(:)

    contiguous :: data

    logical, intent(in), optional              :: overwrite
    integer, intent(in), optional              :: comp_level

    logical          :: overwrite_
    integer          :: hdf_err
    logical          :: item_exists
    integer(HSIZE_T) :: item_shape(1)
    integer(HID_T)   :: dspace_id
    integer(HID_T)   :: plist_id
    integer(HID_T)   :: mem_type_id
    integer(HID_T)   :: file_type_id
    type(C_PTR)      :: data_ptr
    integer(HID_T)   :: item_id

    if(PRESENT(overwrite)) then
       overwrite_ = overwrite
    else
       overwrite_ = .FALSE.
    endif

    ! Write the item

    item_exists = attr_exists_(hg, item_name)
    overwrite_ = overwrite_ .AND. item_exists

    if(overwrite_) then

       call get_attr_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 1168 <core_hgroup:write_attr_c_sp' &
 & // '_1_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

       mem_type_id = complex_mem_type_(SP)

       data_ptr = C_LOC(data)

       call h5aopen_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1168 <core_hgroup:write_attr_c_sp_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aopen_f'
   stop 'Program aborted'
endif

       call h5awrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1168 <core_hgroup:write_attr_c_sp_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5awrite_f'
   stop 'Program aborted'
endif

       call h5aclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1168 <core_hgroup:write_attr_c_sp_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aclose_f'
   stop 'Program aborted'
endif

    else

       call h5screate_simple_f( 1,  INT(SHAPE(data), HSIZE_T),  dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1168 <core_hgroup:write_attr_c_sp_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5screate_simple_f'
   stop 'Program aborted'
endif

       call h5pcreate_f( H5P_ATTRIBUTE_CREATE_F,  plist_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1168 <core_hgroup:write_attr_c_sp_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pcreate_f'
   stop 'Program aborted'
endif

    if(.NOT. (.NOT. PRESENT(comp_level))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''.NOT. PRESENT(comp_level)'' failed at line 1168 <core_hgroup:write_attr_c_sp_1_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Attributes cannot be compressed'
      stop
    endif

       mem_type_id = complex_mem_type_(SP)
       file_type_id = complex_file_type_(SP)

       data_ptr = C_LOC(data)

       call h5acreate_f( hg%group_id,  item_name,  file_type_id,  dspace_id,  item_id, hdf_err,  acpl_id=plist_id)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1168 <core_hgroup:write_attr_c_sp_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5acreate_f'
   stop 'Program aborted'
endif

       call h5awrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1168 <core_hgroup:write_attr_c_sp_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5awrite_f'
   stop 'Program aborted'
endif

       call h5pclose_f( plist_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1168 <core_hgroup:write_attr_c_sp_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pclose_f'
   stop 'Program aborted'
endif

       call h5aclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1168 <core_hgroup:write_attr_c_sp_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aclose_f'
   stop 'Program aborted'
endif

       call h5sclose_f( dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1168 <core_hgroup:write_attr_c_sp_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5sclose_f'
   stop 'Program aborted'
endif

    endif

    ! Finish

    return

  end subroutine write_attr_c_sp_1_

  subroutine write_attr_c_sp_2_ (hg, item_name, data, overwrite, comp_level)

    type(hgroup_t), intent(inout)              :: hg
    character(*), intent(in)                   :: item_name
    complex(SP), target, intent(in) :: data(:,:)

    contiguous :: data

    logical, intent(in), optional              :: overwrite
    integer, intent(in), optional              :: comp_level

    logical          :: overwrite_
    integer          :: hdf_err
    logical          :: item_exists
    integer(HSIZE_T) :: item_shape(2)
    integer(HID_T)   :: dspace_id
    integer(HID_T)   :: plist_id
    integer(HID_T)   :: mem_type_id
    integer(HID_T)   :: file_type_id
    type(C_PTR)      :: data_ptr
    integer(HID_T)   :: item_id

    if(PRESENT(overwrite)) then
       overwrite_ = overwrite
    else
       overwrite_ = .FALSE.
    endif

    ! Write the item

    item_exists = attr_exists_(hg, item_name)
    overwrite_ = overwrite_ .AND. item_exists

    if(overwrite_) then

       call get_attr_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 1169 <core_hgroup:write_attr_c_sp' &
 & // '_2_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

       mem_type_id = complex_mem_type_(SP)

       data_ptr = C_LOC(data)

       call h5aopen_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1169 <core_hgroup:write_attr_c_sp_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aopen_f'
   stop 'Program aborted'
endif

       call h5awrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1169 <core_hgroup:write_attr_c_sp_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5awrite_f'
   stop 'Program aborted'
endif

       call h5aclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1169 <core_hgroup:write_attr_c_sp_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aclose_f'
   stop 'Program aborted'
endif

    else

       call h5screate_simple_f( 2,  INT(SHAPE(data), HSIZE_T),  dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1169 <core_hgroup:write_attr_c_sp_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5screate_simple_f'
   stop 'Program aborted'
endif

       call h5pcreate_f( H5P_ATTRIBUTE_CREATE_F,  plist_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1169 <core_hgroup:write_attr_c_sp_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pcreate_f'
   stop 'Program aborted'
endif

    if(.NOT. (.NOT. PRESENT(comp_level))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''.NOT. PRESENT(comp_level)'' failed at line 1169 <core_hgroup:write_attr_c_sp_2_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Attributes cannot be compressed'
      stop
    endif

       mem_type_id = complex_mem_type_(SP)
       file_type_id = complex_file_type_(SP)

       data_ptr = C_LOC(data)

       call h5acreate_f( hg%group_id,  item_name,  file_type_id,  dspace_id,  item_id, hdf_err,  acpl_id=plist_id)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1169 <core_hgroup:write_attr_c_sp_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5acreate_f'
   stop 'Program aborted'
endif

       call h5awrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1169 <core_hgroup:write_attr_c_sp_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5awrite_f'
   stop 'Program aborted'
endif

       call h5pclose_f( plist_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1169 <core_hgroup:write_attr_c_sp_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pclose_f'
   stop 'Program aborted'
endif

       call h5aclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1169 <core_hgroup:write_attr_c_sp_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aclose_f'
   stop 'Program aborted'
endif

       call h5sclose_f( dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1169 <core_hgroup:write_attr_c_sp_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5sclose_f'
   stop 'Program aborted'
endif

    endif

    ! Finish

    return

  end subroutine write_attr_c_sp_2_

  subroutine write_attr_c_sp_3_ (hg, item_name, data, overwrite, comp_level)

    type(hgroup_t), intent(inout)              :: hg
    character(*), intent(in)                   :: item_name
    complex(SP), target, intent(in) :: data(:,:,:)

    contiguous :: data

    logical, intent(in), optional              :: overwrite
    integer, intent(in), optional              :: comp_level

    logical          :: overwrite_
    integer          :: hdf_err
    logical          :: item_exists
    integer(HSIZE_T) :: item_shape(3)
    integer(HID_T)   :: dspace_id
    integer(HID_T)   :: plist_id
    integer(HID_T)   :: mem_type_id
    integer(HID_T)   :: file_type_id
    type(C_PTR)      :: data_ptr
    integer(HID_T)   :: item_id

    if(PRESENT(overwrite)) then
       overwrite_ = overwrite
    else
       overwrite_ = .FALSE.
    endif

    ! Write the item

    item_exists = attr_exists_(hg, item_name)
    overwrite_ = overwrite_ .AND. item_exists

    if(overwrite_) then

       call get_attr_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 1170 <core_hgroup:write_attr_c_sp' &
 & // '_3_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

       mem_type_id = complex_mem_type_(SP)

       data_ptr = C_LOC(data)

       call h5aopen_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1170 <core_hgroup:write_attr_c_sp_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aopen_f'
   stop 'Program aborted'
endif

       call h5awrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1170 <core_hgroup:write_attr_c_sp_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5awrite_f'
   stop 'Program aborted'
endif

       call h5aclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1170 <core_hgroup:write_attr_c_sp_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aclose_f'
   stop 'Program aborted'
endif

    else

       call h5screate_simple_f( 3,  INT(SHAPE(data), HSIZE_T),  dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1170 <core_hgroup:write_attr_c_sp_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5screate_simple_f'
   stop 'Program aborted'
endif

       call h5pcreate_f( H5P_ATTRIBUTE_CREATE_F,  plist_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1170 <core_hgroup:write_attr_c_sp_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pcreate_f'
   stop 'Program aborted'
endif

    if(.NOT. (.NOT. PRESENT(comp_level))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''.NOT. PRESENT(comp_level)'' failed at line 1170 <core_hgroup:write_attr_c_sp_3_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Attributes cannot be compressed'
      stop
    endif

       mem_type_id = complex_mem_type_(SP)
       file_type_id = complex_file_type_(SP)

       data_ptr = C_LOC(data)

       call h5acreate_f( hg%group_id,  item_name,  file_type_id,  dspace_id,  item_id, hdf_err,  acpl_id=plist_id)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1170 <core_hgroup:write_attr_c_sp_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5acreate_f'
   stop 'Program aborted'
endif

       call h5awrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1170 <core_hgroup:write_attr_c_sp_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5awrite_f'
   stop 'Program aborted'
endif

       call h5pclose_f( plist_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1170 <core_hgroup:write_attr_c_sp_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pclose_f'
   stop 'Program aborted'
endif

       call h5aclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1170 <core_hgroup:write_attr_c_sp_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aclose_f'
   stop 'Program aborted'
endif

       call h5sclose_f( dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1170 <core_hgroup:write_attr_c_sp_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5sclose_f'
   stop 'Program aborted'
endif

    endif

    ! Finish

    return

  end subroutine write_attr_c_sp_3_

  subroutine write_attr_c_sp_4_ (hg, item_name, data, overwrite, comp_level)

    type(hgroup_t), intent(inout)              :: hg
    character(*), intent(in)                   :: item_name
    complex(SP), target, intent(in) :: data(:,:,:,:)

    contiguous :: data

    logical, intent(in), optional              :: overwrite
    integer, intent(in), optional              :: comp_level

    logical          :: overwrite_
    integer          :: hdf_err
    logical          :: item_exists
    integer(HSIZE_T) :: item_shape(4)
    integer(HID_T)   :: dspace_id
    integer(HID_T)   :: plist_id
    integer(HID_T)   :: mem_type_id
    integer(HID_T)   :: file_type_id
    type(C_PTR)      :: data_ptr
    integer(HID_T)   :: item_id

    if(PRESENT(overwrite)) then
       overwrite_ = overwrite
    else
       overwrite_ = .FALSE.
    endif

    ! Write the item

    item_exists = attr_exists_(hg, item_name)
    overwrite_ = overwrite_ .AND. item_exists

    if(overwrite_) then

       call get_attr_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 1171 <core_hgroup:write_attr_c_sp' &
 & // '_4_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

       mem_type_id = complex_mem_type_(SP)

       data_ptr = C_LOC(data)

       call h5aopen_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1171 <core_hgroup:write_attr_c_sp_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aopen_f'
   stop 'Program aborted'
endif

       call h5awrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1171 <core_hgroup:write_attr_c_sp_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5awrite_f'
   stop 'Program aborted'
endif

       call h5aclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1171 <core_hgroup:write_attr_c_sp_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aclose_f'
   stop 'Program aborted'
endif

    else

       call h5screate_simple_f( 4,  INT(SHAPE(data), HSIZE_T),  dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1171 <core_hgroup:write_attr_c_sp_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5screate_simple_f'
   stop 'Program aborted'
endif

       call h5pcreate_f( H5P_ATTRIBUTE_CREATE_F,  plist_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1171 <core_hgroup:write_attr_c_sp_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pcreate_f'
   stop 'Program aborted'
endif

    if(.NOT. (.NOT. PRESENT(comp_level))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''.NOT. PRESENT(comp_level)'' failed at line 1171 <core_hgroup:write_attr_c_sp_4_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Attributes cannot be compressed'
      stop
    endif

       mem_type_id = complex_mem_type_(SP)
       file_type_id = complex_file_type_(SP)

       data_ptr = C_LOC(data)

       call h5acreate_f( hg%group_id,  item_name,  file_type_id,  dspace_id,  item_id, hdf_err,  acpl_id=plist_id)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1171 <core_hgroup:write_attr_c_sp_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5acreate_f'
   stop 'Program aborted'
endif

       call h5awrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1171 <core_hgroup:write_attr_c_sp_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5awrite_f'
   stop 'Program aborted'
endif

       call h5pclose_f( plist_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1171 <core_hgroup:write_attr_c_sp_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pclose_f'
   stop 'Program aborted'
endif

       call h5aclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1171 <core_hgroup:write_attr_c_sp_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aclose_f'
   stop 'Program aborted'
endif

       call h5sclose_f( dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1171 <core_hgroup:write_attr_c_sp_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5sclose_f'
   stop 'Program aborted'
endif

    endif

    ! Finish

    return

  end subroutine write_attr_c_sp_4_

  subroutine write_attr_c_dp_0_ (hg, item_name, data, overwrite, comp_level)

    type(hgroup_t), intent(inout)              :: hg
    character(*), intent(in)                   :: item_name
    complex(DP), target, intent(in) :: data

    logical, intent(in), optional              :: overwrite
    integer, intent(in), optional              :: comp_level

    logical          :: overwrite_
    integer          :: hdf_err
    logical          :: item_exists
    integer(HSIZE_T) :: item_shape(0)
    integer(HID_T)   :: dspace_id
    integer(HID_T)   :: plist_id
    integer(HID_T)   :: mem_type_id
    integer(HID_T)   :: file_type_id
    type(C_PTR)      :: data_ptr
    integer(HID_T)   :: item_id

    if(PRESENT(overwrite)) then
       overwrite_ = overwrite
    else
       overwrite_ = .FALSE.
    endif

    ! Write the item

    item_exists = attr_exists_(hg, item_name)
    overwrite_ = overwrite_ .AND. item_exists

    if(overwrite_) then

       call get_attr_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 1173 <core_hgroup:write_attr_c_dp' &
 & // '_0_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

       mem_type_id = complex_mem_type_(DP)

       data_ptr = C_LOC(data)

       call h5aopen_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1173 <core_hgroup:write_attr_c_dp_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aopen_f'
   stop 'Program aborted'
endif

       call h5awrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1173 <core_hgroup:write_attr_c_dp_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5awrite_f'
   stop 'Program aborted'
endif

       call h5aclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1173 <core_hgroup:write_attr_c_dp_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aclose_f'
   stop 'Program aborted'
endif

    else

       call h5screate_f( H5S_SCALAR_F,  dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1173 <core_hgroup:write_attr_c_dp_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5screate_f'
   stop 'Program aborted'
endif

       call h5pcreate_f( H5P_ATTRIBUTE_CREATE_F,  plist_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1173 <core_hgroup:write_attr_c_dp_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pcreate_f'
   stop 'Program aborted'
endif

    if(.NOT. (.NOT. PRESENT(comp_level))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''.NOT. PRESENT(comp_level)'' failed at line 1173 <core_hgroup:write_attr_c_dp_0_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Attributes cannot be compressed'
      stop
    endif

       mem_type_id = complex_mem_type_(DP)
       file_type_id = complex_file_type_(DP)

       data_ptr = C_LOC(data)

       call h5acreate_f( hg%group_id,  item_name,  file_type_id,  dspace_id,  item_id, hdf_err,  acpl_id=plist_id)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1173 <core_hgroup:write_attr_c_dp_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5acreate_f'
   stop 'Program aborted'
endif

       call h5awrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1173 <core_hgroup:write_attr_c_dp_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5awrite_f'
   stop 'Program aborted'
endif

       call h5pclose_f( plist_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1173 <core_hgroup:write_attr_c_dp_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pclose_f'
   stop 'Program aborted'
endif

       call h5aclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1173 <core_hgroup:write_attr_c_dp_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aclose_f'
   stop 'Program aborted'
endif

       call h5sclose_f( dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1173 <core_hgroup:write_attr_c_dp_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5sclose_f'
   stop 'Program aborted'
endif

    endif

    ! Finish

    return

  end subroutine write_attr_c_dp_0_

  subroutine write_attr_c_dp_1_ (hg, item_name, data, overwrite, comp_level)

    type(hgroup_t), intent(inout)              :: hg
    character(*), intent(in)                   :: item_name
    complex(DP), target, intent(in) :: data(:)

    contiguous :: data

    logical, intent(in), optional              :: overwrite
    integer, intent(in), optional              :: comp_level

    logical          :: overwrite_
    integer          :: hdf_err
    logical          :: item_exists
    integer(HSIZE_T) :: item_shape(1)
    integer(HID_T)   :: dspace_id
    integer(HID_T)   :: plist_id
    integer(HID_T)   :: mem_type_id
    integer(HID_T)   :: file_type_id
    type(C_PTR)      :: data_ptr
    integer(HID_T)   :: item_id

    if(PRESENT(overwrite)) then
       overwrite_ = overwrite
    else
       overwrite_ = .FALSE.
    endif

    ! Write the item

    item_exists = attr_exists_(hg, item_name)
    overwrite_ = overwrite_ .AND. item_exists

    if(overwrite_) then

       call get_attr_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 1174 <core_hgroup:write_attr_c_dp' &
 & // '_1_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

       mem_type_id = complex_mem_type_(DP)

       data_ptr = C_LOC(data)

       call h5aopen_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1174 <core_hgroup:write_attr_c_dp_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aopen_f'
   stop 'Program aborted'
endif

       call h5awrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1174 <core_hgroup:write_attr_c_dp_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5awrite_f'
   stop 'Program aborted'
endif

       call h5aclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1174 <core_hgroup:write_attr_c_dp_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aclose_f'
   stop 'Program aborted'
endif

    else

       call h5screate_simple_f( 1,  INT(SHAPE(data), HSIZE_T),  dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1174 <core_hgroup:write_attr_c_dp_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5screate_simple_f'
   stop 'Program aborted'
endif

       call h5pcreate_f( H5P_ATTRIBUTE_CREATE_F,  plist_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1174 <core_hgroup:write_attr_c_dp_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pcreate_f'
   stop 'Program aborted'
endif

    if(.NOT. (.NOT. PRESENT(comp_level))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''.NOT. PRESENT(comp_level)'' failed at line 1174 <core_hgroup:write_attr_c_dp_1_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Attributes cannot be compressed'
      stop
    endif

       mem_type_id = complex_mem_type_(DP)
       file_type_id = complex_file_type_(DP)

       data_ptr = C_LOC(data)

       call h5acreate_f( hg%group_id,  item_name,  file_type_id,  dspace_id,  item_id, hdf_err,  acpl_id=plist_id)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1174 <core_hgroup:write_attr_c_dp_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5acreate_f'
   stop 'Program aborted'
endif

       call h5awrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1174 <core_hgroup:write_attr_c_dp_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5awrite_f'
   stop 'Program aborted'
endif

       call h5pclose_f( plist_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1174 <core_hgroup:write_attr_c_dp_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pclose_f'
   stop 'Program aborted'
endif

       call h5aclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1174 <core_hgroup:write_attr_c_dp_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aclose_f'
   stop 'Program aborted'
endif

       call h5sclose_f( dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1174 <core_hgroup:write_attr_c_dp_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5sclose_f'
   stop 'Program aborted'
endif

    endif

    ! Finish

    return

  end subroutine write_attr_c_dp_1_

  subroutine write_attr_c_dp_2_ (hg, item_name, data, overwrite, comp_level)

    type(hgroup_t), intent(inout)              :: hg
    character(*), intent(in)                   :: item_name
    complex(DP), target, intent(in) :: data(:,:)

    contiguous :: data

    logical, intent(in), optional              :: overwrite
    integer, intent(in), optional              :: comp_level

    logical          :: overwrite_
    integer          :: hdf_err
    logical          :: item_exists
    integer(HSIZE_T) :: item_shape(2)
    integer(HID_T)   :: dspace_id
    integer(HID_T)   :: plist_id
    integer(HID_T)   :: mem_type_id
    integer(HID_T)   :: file_type_id
    type(C_PTR)      :: data_ptr
    integer(HID_T)   :: item_id

    if(PRESENT(overwrite)) then
       overwrite_ = overwrite
    else
       overwrite_ = .FALSE.
    endif

    ! Write the item

    item_exists = attr_exists_(hg, item_name)
    overwrite_ = overwrite_ .AND. item_exists

    if(overwrite_) then

       call get_attr_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 1175 <core_hgroup:write_attr_c_dp' &
 & // '_2_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

       mem_type_id = complex_mem_type_(DP)

       data_ptr = C_LOC(data)

       call h5aopen_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1175 <core_hgroup:write_attr_c_dp_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aopen_f'
   stop 'Program aborted'
endif

       call h5awrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1175 <core_hgroup:write_attr_c_dp_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5awrite_f'
   stop 'Program aborted'
endif

       call h5aclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1175 <core_hgroup:write_attr_c_dp_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aclose_f'
   stop 'Program aborted'
endif

    else

       call h5screate_simple_f( 2,  INT(SHAPE(data), HSIZE_T),  dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1175 <core_hgroup:write_attr_c_dp_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5screate_simple_f'
   stop 'Program aborted'
endif

       call h5pcreate_f( H5P_ATTRIBUTE_CREATE_F,  plist_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1175 <core_hgroup:write_attr_c_dp_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pcreate_f'
   stop 'Program aborted'
endif

    if(.NOT. (.NOT. PRESENT(comp_level))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''.NOT. PRESENT(comp_level)'' failed at line 1175 <core_hgroup:write_attr_c_dp_2_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Attributes cannot be compressed'
      stop
    endif

       mem_type_id = complex_mem_type_(DP)
       file_type_id = complex_file_type_(DP)

       data_ptr = C_LOC(data)

       call h5acreate_f( hg%group_id,  item_name,  file_type_id,  dspace_id,  item_id, hdf_err,  acpl_id=plist_id)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1175 <core_hgroup:write_attr_c_dp_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5acreate_f'
   stop 'Program aborted'
endif

       call h5awrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1175 <core_hgroup:write_attr_c_dp_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5awrite_f'
   stop 'Program aborted'
endif

       call h5pclose_f( plist_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1175 <core_hgroup:write_attr_c_dp_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pclose_f'
   stop 'Program aborted'
endif

       call h5aclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1175 <core_hgroup:write_attr_c_dp_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aclose_f'
   stop 'Program aborted'
endif

       call h5sclose_f( dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1175 <core_hgroup:write_attr_c_dp_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5sclose_f'
   stop 'Program aborted'
endif

    endif

    ! Finish

    return

  end subroutine write_attr_c_dp_2_

  subroutine write_attr_c_dp_3_ (hg, item_name, data, overwrite, comp_level)

    type(hgroup_t), intent(inout)              :: hg
    character(*), intent(in)                   :: item_name
    complex(DP), target, intent(in) :: data(:,:,:)

    contiguous :: data

    logical, intent(in), optional              :: overwrite
    integer, intent(in), optional              :: comp_level

    logical          :: overwrite_
    integer          :: hdf_err
    logical          :: item_exists
    integer(HSIZE_T) :: item_shape(3)
    integer(HID_T)   :: dspace_id
    integer(HID_T)   :: plist_id
    integer(HID_T)   :: mem_type_id
    integer(HID_T)   :: file_type_id
    type(C_PTR)      :: data_ptr
    integer(HID_T)   :: item_id

    if(PRESENT(overwrite)) then
       overwrite_ = overwrite
    else
       overwrite_ = .FALSE.
    endif

    ! Write the item

    item_exists = attr_exists_(hg, item_name)
    overwrite_ = overwrite_ .AND. item_exists

    if(overwrite_) then

       call get_attr_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 1176 <core_hgroup:write_attr_c_dp' &
 & // '_3_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

       mem_type_id = complex_mem_type_(DP)

       data_ptr = C_LOC(data)

       call h5aopen_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1176 <core_hgroup:write_attr_c_dp_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aopen_f'
   stop 'Program aborted'
endif

       call h5awrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1176 <core_hgroup:write_attr_c_dp_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5awrite_f'
   stop 'Program aborted'
endif

       call h5aclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1176 <core_hgroup:write_attr_c_dp_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aclose_f'
   stop 'Program aborted'
endif

    else

       call h5screate_simple_f( 3,  INT(SHAPE(data), HSIZE_T),  dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1176 <core_hgroup:write_attr_c_dp_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5screate_simple_f'
   stop 'Program aborted'
endif

       call h5pcreate_f( H5P_ATTRIBUTE_CREATE_F,  plist_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1176 <core_hgroup:write_attr_c_dp_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pcreate_f'
   stop 'Program aborted'
endif

    if(.NOT. (.NOT. PRESENT(comp_level))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''.NOT. PRESENT(comp_level)'' failed at line 1176 <core_hgroup:write_attr_c_dp_3_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Attributes cannot be compressed'
      stop
    endif

       mem_type_id = complex_mem_type_(DP)
       file_type_id = complex_file_type_(DP)

       data_ptr = C_LOC(data)

       call h5acreate_f( hg%group_id,  item_name,  file_type_id,  dspace_id,  item_id, hdf_err,  acpl_id=plist_id)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1176 <core_hgroup:write_attr_c_dp_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5acreate_f'
   stop 'Program aborted'
endif

       call h5awrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1176 <core_hgroup:write_attr_c_dp_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5awrite_f'
   stop 'Program aborted'
endif

       call h5pclose_f( plist_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1176 <core_hgroup:write_attr_c_dp_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pclose_f'
   stop 'Program aborted'
endif

       call h5aclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1176 <core_hgroup:write_attr_c_dp_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aclose_f'
   stop 'Program aborted'
endif

       call h5sclose_f( dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1176 <core_hgroup:write_attr_c_dp_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5sclose_f'
   stop 'Program aborted'
endif

    endif

    ! Finish

    return

  end subroutine write_attr_c_dp_3_

  subroutine write_attr_c_dp_4_ (hg, item_name, data, overwrite, comp_level)

    type(hgroup_t), intent(inout)              :: hg
    character(*), intent(in)                   :: item_name
    complex(DP), target, intent(in) :: data(:,:,:,:)

    contiguous :: data

    logical, intent(in), optional              :: overwrite
    integer, intent(in), optional              :: comp_level

    logical          :: overwrite_
    integer          :: hdf_err
    logical          :: item_exists
    integer(HSIZE_T) :: item_shape(4)
    integer(HID_T)   :: dspace_id
    integer(HID_T)   :: plist_id
    integer(HID_T)   :: mem_type_id
    integer(HID_T)   :: file_type_id
    type(C_PTR)      :: data_ptr
    integer(HID_T)   :: item_id

    if(PRESENT(overwrite)) then
       overwrite_ = overwrite
    else
       overwrite_ = .FALSE.
    endif

    ! Write the item

    item_exists = attr_exists_(hg, item_name)
    overwrite_ = overwrite_ .AND. item_exists

    if(overwrite_) then

       call get_attr_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 1177 <core_hgroup:write_attr_c_dp' &
 & // '_4_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

       mem_type_id = complex_mem_type_(DP)

       data_ptr = C_LOC(data)

       call h5aopen_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1177 <core_hgroup:write_attr_c_dp_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aopen_f'
   stop 'Program aborted'
endif

       call h5awrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1177 <core_hgroup:write_attr_c_dp_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5awrite_f'
   stop 'Program aborted'
endif

       call h5aclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1177 <core_hgroup:write_attr_c_dp_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aclose_f'
   stop 'Program aborted'
endif

    else

       call h5screate_simple_f( 4,  INT(SHAPE(data), HSIZE_T),  dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1177 <core_hgroup:write_attr_c_dp_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5screate_simple_f'
   stop 'Program aborted'
endif

       call h5pcreate_f( H5P_ATTRIBUTE_CREATE_F,  plist_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1177 <core_hgroup:write_attr_c_dp_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pcreate_f'
   stop 'Program aborted'
endif

    if(.NOT. (.NOT. PRESENT(comp_level))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''.NOT. PRESENT(comp_level)'' failed at line 1177 <core_hgroup:write_attr_c_dp_4_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Attributes cannot be compressed'
      stop
    endif

       mem_type_id = complex_mem_type_(DP)
       file_type_id = complex_file_type_(DP)

       data_ptr = C_LOC(data)

       call h5acreate_f( hg%group_id,  item_name,  file_type_id,  dspace_id,  item_id, hdf_err,  acpl_id=plist_id)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1177 <core_hgroup:write_attr_c_dp_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5acreate_f'
   stop 'Program aborted'
endif

       call h5awrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1177 <core_hgroup:write_attr_c_dp_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5awrite_f'
   stop 'Program aborted'
endif

       call h5pclose_f( plist_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1177 <core_hgroup:write_attr_c_dp_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pclose_f'
   stop 'Program aborted'
endif

       call h5aclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1177 <core_hgroup:write_attr_c_dp_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aclose_f'
   stop 'Program aborted'
endif

       call h5sclose_f( dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1177 <core_hgroup:write_attr_c_dp_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5sclose_f'
   stop 'Program aborted'
endif

    endif

    ! Finish

    return

  end subroutine write_attr_c_dp_4_

  subroutine write_dset_i_i4_0_ (hg, item_name, data, overwrite, comp_level)

    type(hgroup_t), intent(inout)              :: hg
    character(*), intent(in)                   :: item_name
    integer(I4), target, intent(in) :: data

    logical, intent(in), optional              :: overwrite
    integer, intent(in), optional              :: comp_level

    logical          :: overwrite_
    integer          :: hdf_err
    logical          :: item_exists
    integer(HSIZE_T) :: item_shape(0)
    integer(HID_T)   :: dspace_id
    integer(HID_T)   :: plist_id
    integer(HID_T)   :: mem_type_id
    integer(HID_T)   :: file_type_id
    type(C_PTR)      :: data_ptr
    integer(HID_T)   :: item_id

    if(PRESENT(overwrite)) then
       overwrite_ = overwrite
    else
       overwrite_ = .FALSE.
    endif

    ! Write the item

    item_exists = dset_exists_(hg, item_name)
    overwrite_ = overwrite_ .AND. item_exists

    if(overwrite_) then

       call get_dset_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 1179 <core_hgroup:write_dset_i_i4' &
 & // '_0_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

       mem_type_id = integer_mem_type_(I4)

       data_ptr = C_LOC(data)

       call h5dopen_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1179 <core_hgroup:write_dset_i_i4_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dopen_f'
   stop 'Program aborted'
endif

       call h5dwrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1179 <core_hgroup:write_dset_i_i4_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dwrite_f'
   stop 'Program aborted'
endif

       call h5dclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1179 <core_hgroup:write_dset_i_i4_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dclose_f'
   stop 'Program aborted'
endif

    else

       call h5screate_f( H5S_SCALAR_F,  dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1179 <core_hgroup:write_dset_i_i4_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5screate_f'
   stop 'Program aborted'
endif

       call h5pcreate_f( H5P_DATASET_CREATE_F,  plist_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1179 <core_hgroup:write_dset_i_i4_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pcreate_f'
   stop 'Program aborted'
endif

       if (PRESENT(comp_level)) then
          call h5pset_chunk_f( plist_id,  0,  INT(SHAPE(data), HSIZE_T), hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1179 <core_hgroup:write_dset_i_i4_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pset_chunk_f'
   stop 'Program aborted'
endif

          call h5pset_deflate_f( plist_id,  comp_level, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1179 <core_hgroup:write_dset_i_i4_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pset_deflate_f'
   stop 'Program aborted'
endif

       endif

       mem_type_id = integer_mem_type_(I4)
       file_type_id = integer_file_type_(I4)

       data_ptr = C_LOC(data)

       call h5dcreate_f( hg%group_id,  item_name,  file_type_id,  dspace_id,  item_id, hdf_err,  dcpl_id=plist_id)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1179 <core_hgroup:write_dset_i_i4_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dcreate_f'
   stop 'Program aborted'
endif

       call h5dwrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1179 <core_hgroup:write_dset_i_i4_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dwrite_f'
   stop 'Program aborted'
endif

       call h5pclose_f( plist_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1179 <core_hgroup:write_dset_i_i4_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pclose_f'
   stop 'Program aborted'
endif

       call h5dclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1179 <core_hgroup:write_dset_i_i4_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dclose_f'
   stop 'Program aborted'
endif

       call h5sclose_f( dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1179 <core_hgroup:write_dset_i_i4_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5sclose_f'
   stop 'Program aborted'
endif

    endif

    ! Finish

    return

  end subroutine write_dset_i_i4_0_

  subroutine write_dset_i_i4_1_ (hg, item_name, data, overwrite, comp_level)

    type(hgroup_t), intent(inout)              :: hg
    character(*), intent(in)                   :: item_name
    integer(I4), target, intent(in) :: data(:)

    contiguous :: data

    logical, intent(in), optional              :: overwrite
    integer, intent(in), optional              :: comp_level

    logical          :: overwrite_
    integer          :: hdf_err
    logical          :: item_exists
    integer(HSIZE_T) :: item_shape(1)
    integer(HID_T)   :: dspace_id
    integer(HID_T)   :: plist_id
    integer(HID_T)   :: mem_type_id
    integer(HID_T)   :: file_type_id
    type(C_PTR)      :: data_ptr
    integer(HID_T)   :: item_id

    if(PRESENT(overwrite)) then
       overwrite_ = overwrite
    else
       overwrite_ = .FALSE.
    endif

    ! Write the item

    item_exists = dset_exists_(hg, item_name)
    overwrite_ = overwrite_ .AND. item_exists

    if(overwrite_) then

       call get_dset_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 1180 <core_hgroup:write_dset_i_i4' &
 & // '_1_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

       mem_type_id = integer_mem_type_(I4)

       data_ptr = C_LOC(data)

       call h5dopen_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1180 <core_hgroup:write_dset_i_i4_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dopen_f'
   stop 'Program aborted'
endif

       call h5dwrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1180 <core_hgroup:write_dset_i_i4_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dwrite_f'
   stop 'Program aborted'
endif

       call h5dclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1180 <core_hgroup:write_dset_i_i4_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dclose_f'
   stop 'Program aborted'
endif

    else

       call h5screate_simple_f( 1,  INT(SHAPE(data), HSIZE_T),  dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1180 <core_hgroup:write_dset_i_i4_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5screate_simple_f'
   stop 'Program aborted'
endif

       call h5pcreate_f( H5P_DATASET_CREATE_F,  plist_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1180 <core_hgroup:write_dset_i_i4_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pcreate_f'
   stop 'Program aborted'
endif

       if (PRESENT(comp_level)) then
          call h5pset_chunk_f( plist_id,  1,  INT(SHAPE(data), HSIZE_T), hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1180 <core_hgroup:write_dset_i_i4_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pset_chunk_f'
   stop 'Program aborted'
endif

          call h5pset_deflate_f( plist_id,  comp_level, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1180 <core_hgroup:write_dset_i_i4_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pset_deflate_f'
   stop 'Program aborted'
endif

       endif

       mem_type_id = integer_mem_type_(I4)
       file_type_id = integer_file_type_(I4)

       data_ptr = C_LOC(data)

       call h5dcreate_f( hg%group_id,  item_name,  file_type_id,  dspace_id,  item_id, hdf_err,  dcpl_id=plist_id)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1180 <core_hgroup:write_dset_i_i4_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dcreate_f'
   stop 'Program aborted'
endif

       call h5dwrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1180 <core_hgroup:write_dset_i_i4_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dwrite_f'
   stop 'Program aborted'
endif

       call h5pclose_f( plist_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1180 <core_hgroup:write_dset_i_i4_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pclose_f'
   stop 'Program aborted'
endif

       call h5dclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1180 <core_hgroup:write_dset_i_i4_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dclose_f'
   stop 'Program aborted'
endif

       call h5sclose_f( dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1180 <core_hgroup:write_dset_i_i4_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5sclose_f'
   stop 'Program aborted'
endif

    endif

    ! Finish

    return

  end subroutine write_dset_i_i4_1_

  subroutine write_dset_i_i4_2_ (hg, item_name, data, overwrite, comp_level)

    type(hgroup_t), intent(inout)              :: hg
    character(*), intent(in)                   :: item_name
    integer(I4), target, intent(in) :: data(:,:)

    contiguous :: data

    logical, intent(in), optional              :: overwrite
    integer, intent(in), optional              :: comp_level

    logical          :: overwrite_
    integer          :: hdf_err
    logical          :: item_exists
    integer(HSIZE_T) :: item_shape(2)
    integer(HID_T)   :: dspace_id
    integer(HID_T)   :: plist_id
    integer(HID_T)   :: mem_type_id
    integer(HID_T)   :: file_type_id
    type(C_PTR)      :: data_ptr
    integer(HID_T)   :: item_id

    if(PRESENT(overwrite)) then
       overwrite_ = overwrite
    else
       overwrite_ = .FALSE.
    endif

    ! Write the item

    item_exists = dset_exists_(hg, item_name)
    overwrite_ = overwrite_ .AND. item_exists

    if(overwrite_) then

       call get_dset_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 1181 <core_hgroup:write_dset_i_i4' &
 & // '_2_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

       mem_type_id = integer_mem_type_(I4)

       data_ptr = C_LOC(data)

       call h5dopen_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1181 <core_hgroup:write_dset_i_i4_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dopen_f'
   stop 'Program aborted'
endif

       call h5dwrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1181 <core_hgroup:write_dset_i_i4_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dwrite_f'
   stop 'Program aborted'
endif

       call h5dclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1181 <core_hgroup:write_dset_i_i4_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dclose_f'
   stop 'Program aborted'
endif

    else

       call h5screate_simple_f( 2,  INT(SHAPE(data), HSIZE_T),  dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1181 <core_hgroup:write_dset_i_i4_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5screate_simple_f'
   stop 'Program aborted'
endif

       call h5pcreate_f( H5P_DATASET_CREATE_F,  plist_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1181 <core_hgroup:write_dset_i_i4_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pcreate_f'
   stop 'Program aborted'
endif

       if (PRESENT(comp_level)) then
          call h5pset_chunk_f( plist_id,  2,  INT(SHAPE(data), HSIZE_T), hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1181 <core_hgroup:write_dset_i_i4_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pset_chunk_f'
   stop 'Program aborted'
endif

          call h5pset_deflate_f( plist_id,  comp_level, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1181 <core_hgroup:write_dset_i_i4_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pset_deflate_f'
   stop 'Program aborted'
endif

       endif

       mem_type_id = integer_mem_type_(I4)
       file_type_id = integer_file_type_(I4)

       data_ptr = C_LOC(data)

       call h5dcreate_f( hg%group_id,  item_name,  file_type_id,  dspace_id,  item_id, hdf_err,  dcpl_id=plist_id)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1181 <core_hgroup:write_dset_i_i4_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dcreate_f'
   stop 'Program aborted'
endif

       call h5dwrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1181 <core_hgroup:write_dset_i_i4_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dwrite_f'
   stop 'Program aborted'
endif

       call h5pclose_f( plist_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1181 <core_hgroup:write_dset_i_i4_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pclose_f'
   stop 'Program aborted'
endif

       call h5dclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1181 <core_hgroup:write_dset_i_i4_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dclose_f'
   stop 'Program aborted'
endif

       call h5sclose_f( dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1181 <core_hgroup:write_dset_i_i4_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5sclose_f'
   stop 'Program aborted'
endif

    endif

    ! Finish

    return

  end subroutine write_dset_i_i4_2_

  subroutine write_dset_i_i4_3_ (hg, item_name, data, overwrite, comp_level)

    type(hgroup_t), intent(inout)              :: hg
    character(*), intent(in)                   :: item_name
    integer(I4), target, intent(in) :: data(:,:,:)

    contiguous :: data

    logical, intent(in), optional              :: overwrite
    integer, intent(in), optional              :: comp_level

    logical          :: overwrite_
    integer          :: hdf_err
    logical          :: item_exists
    integer(HSIZE_T) :: item_shape(3)
    integer(HID_T)   :: dspace_id
    integer(HID_T)   :: plist_id
    integer(HID_T)   :: mem_type_id
    integer(HID_T)   :: file_type_id
    type(C_PTR)      :: data_ptr
    integer(HID_T)   :: item_id

    if(PRESENT(overwrite)) then
       overwrite_ = overwrite
    else
       overwrite_ = .FALSE.
    endif

    ! Write the item

    item_exists = dset_exists_(hg, item_name)
    overwrite_ = overwrite_ .AND. item_exists

    if(overwrite_) then

       call get_dset_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 1182 <core_hgroup:write_dset_i_i4' &
 & // '_3_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

       mem_type_id = integer_mem_type_(I4)

       data_ptr = C_LOC(data)

       call h5dopen_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1182 <core_hgroup:write_dset_i_i4_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dopen_f'
   stop 'Program aborted'
endif

       call h5dwrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1182 <core_hgroup:write_dset_i_i4_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dwrite_f'
   stop 'Program aborted'
endif

       call h5dclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1182 <core_hgroup:write_dset_i_i4_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dclose_f'
   stop 'Program aborted'
endif

    else

       call h5screate_simple_f( 3,  INT(SHAPE(data), HSIZE_T),  dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1182 <core_hgroup:write_dset_i_i4_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5screate_simple_f'
   stop 'Program aborted'
endif

       call h5pcreate_f( H5P_DATASET_CREATE_F,  plist_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1182 <core_hgroup:write_dset_i_i4_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pcreate_f'
   stop 'Program aborted'
endif

       if (PRESENT(comp_level)) then
          call h5pset_chunk_f( plist_id,  3,  INT(SHAPE(data), HSIZE_T), hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1182 <core_hgroup:write_dset_i_i4_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pset_chunk_f'
   stop 'Program aborted'
endif

          call h5pset_deflate_f( plist_id,  comp_level, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1182 <core_hgroup:write_dset_i_i4_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pset_deflate_f'
   stop 'Program aborted'
endif

       endif

       mem_type_id = integer_mem_type_(I4)
       file_type_id = integer_file_type_(I4)

       data_ptr = C_LOC(data)

       call h5dcreate_f( hg%group_id,  item_name,  file_type_id,  dspace_id,  item_id, hdf_err,  dcpl_id=plist_id)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1182 <core_hgroup:write_dset_i_i4_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dcreate_f'
   stop 'Program aborted'
endif

       call h5dwrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1182 <core_hgroup:write_dset_i_i4_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dwrite_f'
   stop 'Program aborted'
endif

       call h5pclose_f( plist_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1182 <core_hgroup:write_dset_i_i4_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pclose_f'
   stop 'Program aborted'
endif

       call h5dclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1182 <core_hgroup:write_dset_i_i4_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dclose_f'
   stop 'Program aborted'
endif

       call h5sclose_f( dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1182 <core_hgroup:write_dset_i_i4_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5sclose_f'
   stop 'Program aborted'
endif

    endif

    ! Finish

    return

  end subroutine write_dset_i_i4_3_

  subroutine write_dset_i_i4_4_ (hg, item_name, data, overwrite, comp_level)

    type(hgroup_t), intent(inout)              :: hg
    character(*), intent(in)                   :: item_name
    integer(I4), target, intent(in) :: data(:,:,:,:)

    contiguous :: data

    logical, intent(in), optional              :: overwrite
    integer, intent(in), optional              :: comp_level

    logical          :: overwrite_
    integer          :: hdf_err
    logical          :: item_exists
    integer(HSIZE_T) :: item_shape(4)
    integer(HID_T)   :: dspace_id
    integer(HID_T)   :: plist_id
    integer(HID_T)   :: mem_type_id
    integer(HID_T)   :: file_type_id
    type(C_PTR)      :: data_ptr
    integer(HID_T)   :: item_id

    if(PRESENT(overwrite)) then
       overwrite_ = overwrite
    else
       overwrite_ = .FALSE.
    endif

    ! Write the item

    item_exists = dset_exists_(hg, item_name)
    overwrite_ = overwrite_ .AND. item_exists

    if(overwrite_) then

       call get_dset_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 1183 <core_hgroup:write_dset_i_i4' &
 & // '_4_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

       mem_type_id = integer_mem_type_(I4)

       data_ptr = C_LOC(data)

       call h5dopen_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1183 <core_hgroup:write_dset_i_i4_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dopen_f'
   stop 'Program aborted'
endif

       call h5dwrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1183 <core_hgroup:write_dset_i_i4_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dwrite_f'
   stop 'Program aborted'
endif

       call h5dclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1183 <core_hgroup:write_dset_i_i4_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dclose_f'
   stop 'Program aborted'
endif

    else

       call h5screate_simple_f( 4,  INT(SHAPE(data), HSIZE_T),  dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1183 <core_hgroup:write_dset_i_i4_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5screate_simple_f'
   stop 'Program aborted'
endif

       call h5pcreate_f( H5P_DATASET_CREATE_F,  plist_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1183 <core_hgroup:write_dset_i_i4_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pcreate_f'
   stop 'Program aborted'
endif

       if (PRESENT(comp_level)) then
          call h5pset_chunk_f( plist_id,  4,  INT(SHAPE(data), HSIZE_T), hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1183 <core_hgroup:write_dset_i_i4_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pset_chunk_f'
   stop 'Program aborted'
endif

          call h5pset_deflate_f( plist_id,  comp_level, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1183 <core_hgroup:write_dset_i_i4_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pset_deflate_f'
   stop 'Program aborted'
endif

       endif

       mem_type_id = integer_mem_type_(I4)
       file_type_id = integer_file_type_(I4)

       data_ptr = C_LOC(data)

       call h5dcreate_f( hg%group_id,  item_name,  file_type_id,  dspace_id,  item_id, hdf_err,  dcpl_id=plist_id)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1183 <core_hgroup:write_dset_i_i4_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dcreate_f'
   stop 'Program aborted'
endif

       call h5dwrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1183 <core_hgroup:write_dset_i_i4_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dwrite_f'
   stop 'Program aborted'
endif

       call h5pclose_f( plist_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1183 <core_hgroup:write_dset_i_i4_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pclose_f'
   stop 'Program aborted'
endif

       call h5dclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1183 <core_hgroup:write_dset_i_i4_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dclose_f'
   stop 'Program aborted'
endif

       call h5sclose_f( dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1183 <core_hgroup:write_dset_i_i4_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5sclose_f'
   stop 'Program aborted'
endif

    endif

    ! Finish

    return

  end subroutine write_dset_i_i4_4_

  subroutine write_dset_i_i8_0_ (hg, item_name, data, overwrite, comp_level)

    type(hgroup_t), intent(inout)              :: hg
    character(*), intent(in)                   :: item_name
    integer(I8), target, intent(in) :: data

    logical, intent(in), optional              :: overwrite
    integer, intent(in), optional              :: comp_level

    logical          :: overwrite_
    integer          :: hdf_err
    logical          :: item_exists
    integer(HSIZE_T) :: item_shape(0)
    integer(HID_T)   :: dspace_id
    integer(HID_T)   :: plist_id
    integer(HID_T)   :: mem_type_id
    integer(HID_T)   :: file_type_id
    type(C_PTR)      :: data_ptr
    integer(HID_T)   :: item_id

    if(PRESENT(overwrite)) then
       overwrite_ = overwrite
    else
       overwrite_ = .FALSE.
    endif

    ! Write the item

    item_exists = dset_exists_(hg, item_name)
    overwrite_ = overwrite_ .AND. item_exists

    if(overwrite_) then

       call get_dset_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 1185 <core_hgroup:write_dset_i_i8' &
 & // '_0_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

       mem_type_id = integer_mem_type_(I8)

       data_ptr = C_LOC(data)

       call h5dopen_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1185 <core_hgroup:write_dset_i_i8_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dopen_f'
   stop 'Program aborted'
endif

       call h5dwrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1185 <core_hgroup:write_dset_i_i8_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dwrite_f'
   stop 'Program aborted'
endif

       call h5dclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1185 <core_hgroup:write_dset_i_i8_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dclose_f'
   stop 'Program aborted'
endif

    else

       call h5screate_f( H5S_SCALAR_F,  dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1185 <core_hgroup:write_dset_i_i8_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5screate_f'
   stop 'Program aborted'
endif

       call h5pcreate_f( H5P_DATASET_CREATE_F,  plist_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1185 <core_hgroup:write_dset_i_i8_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pcreate_f'
   stop 'Program aborted'
endif

       if (PRESENT(comp_level)) then
          call h5pset_chunk_f( plist_id,  0,  INT(SHAPE(data), HSIZE_T), hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1185 <core_hgroup:write_dset_i_i8_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pset_chunk_f'
   stop 'Program aborted'
endif

          call h5pset_deflate_f( plist_id,  comp_level, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1185 <core_hgroup:write_dset_i_i8_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pset_deflate_f'
   stop 'Program aborted'
endif

       endif

       mem_type_id = integer_mem_type_(I8)
       file_type_id = integer_file_type_(I8)

       data_ptr = C_LOC(data)

       call h5dcreate_f( hg%group_id,  item_name,  file_type_id,  dspace_id,  item_id, hdf_err,  dcpl_id=plist_id)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1185 <core_hgroup:write_dset_i_i8_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dcreate_f'
   stop 'Program aborted'
endif

       call h5dwrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1185 <core_hgroup:write_dset_i_i8_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dwrite_f'
   stop 'Program aborted'
endif

       call h5pclose_f( plist_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1185 <core_hgroup:write_dset_i_i8_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pclose_f'
   stop 'Program aborted'
endif

       call h5dclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1185 <core_hgroup:write_dset_i_i8_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dclose_f'
   stop 'Program aborted'
endif

       call h5sclose_f( dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1185 <core_hgroup:write_dset_i_i8_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5sclose_f'
   stop 'Program aborted'
endif

    endif

    ! Finish

    return

  end subroutine write_dset_i_i8_0_

  subroutine write_dset_i_i8_1_ (hg, item_name, data, overwrite, comp_level)

    type(hgroup_t), intent(inout)              :: hg
    character(*), intent(in)                   :: item_name
    integer(I8), target, intent(in) :: data(:)

    contiguous :: data

    logical, intent(in), optional              :: overwrite
    integer, intent(in), optional              :: comp_level

    logical          :: overwrite_
    integer          :: hdf_err
    logical          :: item_exists
    integer(HSIZE_T) :: item_shape(1)
    integer(HID_T)   :: dspace_id
    integer(HID_T)   :: plist_id
    integer(HID_T)   :: mem_type_id
    integer(HID_T)   :: file_type_id
    type(C_PTR)      :: data_ptr
    integer(HID_T)   :: item_id

    if(PRESENT(overwrite)) then
       overwrite_ = overwrite
    else
       overwrite_ = .FALSE.
    endif

    ! Write the item

    item_exists = dset_exists_(hg, item_name)
    overwrite_ = overwrite_ .AND. item_exists

    if(overwrite_) then

       call get_dset_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 1186 <core_hgroup:write_dset_i_i8' &
 & // '_1_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

       mem_type_id = integer_mem_type_(I8)

       data_ptr = C_LOC(data)

       call h5dopen_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1186 <core_hgroup:write_dset_i_i8_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dopen_f'
   stop 'Program aborted'
endif

       call h5dwrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1186 <core_hgroup:write_dset_i_i8_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dwrite_f'
   stop 'Program aborted'
endif

       call h5dclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1186 <core_hgroup:write_dset_i_i8_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dclose_f'
   stop 'Program aborted'
endif

    else

       call h5screate_simple_f( 1,  INT(SHAPE(data), HSIZE_T),  dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1186 <core_hgroup:write_dset_i_i8_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5screate_simple_f'
   stop 'Program aborted'
endif

       call h5pcreate_f( H5P_DATASET_CREATE_F,  plist_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1186 <core_hgroup:write_dset_i_i8_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pcreate_f'
   stop 'Program aborted'
endif

       if (PRESENT(comp_level)) then
          call h5pset_chunk_f( plist_id,  1,  INT(SHAPE(data), HSIZE_T), hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1186 <core_hgroup:write_dset_i_i8_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pset_chunk_f'
   stop 'Program aborted'
endif

          call h5pset_deflate_f( plist_id,  comp_level, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1186 <core_hgroup:write_dset_i_i8_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pset_deflate_f'
   stop 'Program aborted'
endif

       endif

       mem_type_id = integer_mem_type_(I8)
       file_type_id = integer_file_type_(I8)

       data_ptr = C_LOC(data)

       call h5dcreate_f( hg%group_id,  item_name,  file_type_id,  dspace_id,  item_id, hdf_err,  dcpl_id=plist_id)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1186 <core_hgroup:write_dset_i_i8_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dcreate_f'
   stop 'Program aborted'
endif

       call h5dwrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1186 <core_hgroup:write_dset_i_i8_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dwrite_f'
   stop 'Program aborted'
endif

       call h5pclose_f( plist_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1186 <core_hgroup:write_dset_i_i8_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pclose_f'
   stop 'Program aborted'
endif

       call h5dclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1186 <core_hgroup:write_dset_i_i8_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dclose_f'
   stop 'Program aborted'
endif

       call h5sclose_f( dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1186 <core_hgroup:write_dset_i_i8_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5sclose_f'
   stop 'Program aborted'
endif

    endif

    ! Finish

    return

  end subroutine write_dset_i_i8_1_

  subroutine write_dset_i_i8_2_ (hg, item_name, data, overwrite, comp_level)

    type(hgroup_t), intent(inout)              :: hg
    character(*), intent(in)                   :: item_name
    integer(I8), target, intent(in) :: data(:,:)

    contiguous :: data

    logical, intent(in), optional              :: overwrite
    integer, intent(in), optional              :: comp_level

    logical          :: overwrite_
    integer          :: hdf_err
    logical          :: item_exists
    integer(HSIZE_T) :: item_shape(2)
    integer(HID_T)   :: dspace_id
    integer(HID_T)   :: plist_id
    integer(HID_T)   :: mem_type_id
    integer(HID_T)   :: file_type_id
    type(C_PTR)      :: data_ptr
    integer(HID_T)   :: item_id

    if(PRESENT(overwrite)) then
       overwrite_ = overwrite
    else
       overwrite_ = .FALSE.
    endif

    ! Write the item

    item_exists = dset_exists_(hg, item_name)
    overwrite_ = overwrite_ .AND. item_exists

    if(overwrite_) then

       call get_dset_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 1187 <core_hgroup:write_dset_i_i8' &
 & // '_2_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

       mem_type_id = integer_mem_type_(I8)

       data_ptr = C_LOC(data)

       call h5dopen_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1187 <core_hgroup:write_dset_i_i8_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dopen_f'
   stop 'Program aborted'
endif

       call h5dwrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1187 <core_hgroup:write_dset_i_i8_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dwrite_f'
   stop 'Program aborted'
endif

       call h5dclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1187 <core_hgroup:write_dset_i_i8_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dclose_f'
   stop 'Program aborted'
endif

    else

       call h5screate_simple_f( 2,  INT(SHAPE(data), HSIZE_T),  dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1187 <core_hgroup:write_dset_i_i8_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5screate_simple_f'
   stop 'Program aborted'
endif

       call h5pcreate_f( H5P_DATASET_CREATE_F,  plist_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1187 <core_hgroup:write_dset_i_i8_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pcreate_f'
   stop 'Program aborted'
endif

       if (PRESENT(comp_level)) then
          call h5pset_chunk_f( plist_id,  2,  INT(SHAPE(data), HSIZE_T), hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1187 <core_hgroup:write_dset_i_i8_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pset_chunk_f'
   stop 'Program aborted'
endif

          call h5pset_deflate_f( plist_id,  comp_level, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1187 <core_hgroup:write_dset_i_i8_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pset_deflate_f'
   stop 'Program aborted'
endif

       endif

       mem_type_id = integer_mem_type_(I8)
       file_type_id = integer_file_type_(I8)

       data_ptr = C_LOC(data)

       call h5dcreate_f( hg%group_id,  item_name,  file_type_id,  dspace_id,  item_id, hdf_err,  dcpl_id=plist_id)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1187 <core_hgroup:write_dset_i_i8_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dcreate_f'
   stop 'Program aborted'
endif

       call h5dwrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1187 <core_hgroup:write_dset_i_i8_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dwrite_f'
   stop 'Program aborted'
endif

       call h5pclose_f( plist_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1187 <core_hgroup:write_dset_i_i8_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pclose_f'
   stop 'Program aborted'
endif

       call h5dclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1187 <core_hgroup:write_dset_i_i8_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dclose_f'
   stop 'Program aborted'
endif

       call h5sclose_f( dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1187 <core_hgroup:write_dset_i_i8_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5sclose_f'
   stop 'Program aborted'
endif

    endif

    ! Finish

    return

  end subroutine write_dset_i_i8_2_

  subroutine write_dset_i_i8_3_ (hg, item_name, data, overwrite, comp_level)

    type(hgroup_t), intent(inout)              :: hg
    character(*), intent(in)                   :: item_name
    integer(I8), target, intent(in) :: data(:,:,:)

    contiguous :: data

    logical, intent(in), optional              :: overwrite
    integer, intent(in), optional              :: comp_level

    logical          :: overwrite_
    integer          :: hdf_err
    logical          :: item_exists
    integer(HSIZE_T) :: item_shape(3)
    integer(HID_T)   :: dspace_id
    integer(HID_T)   :: plist_id
    integer(HID_T)   :: mem_type_id
    integer(HID_T)   :: file_type_id
    type(C_PTR)      :: data_ptr
    integer(HID_T)   :: item_id

    if(PRESENT(overwrite)) then
       overwrite_ = overwrite
    else
       overwrite_ = .FALSE.
    endif

    ! Write the item

    item_exists = dset_exists_(hg, item_name)
    overwrite_ = overwrite_ .AND. item_exists

    if(overwrite_) then

       call get_dset_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 1188 <core_hgroup:write_dset_i_i8' &
 & // '_3_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

       mem_type_id = integer_mem_type_(I8)

       data_ptr = C_LOC(data)

       call h5dopen_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1188 <core_hgroup:write_dset_i_i8_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dopen_f'
   stop 'Program aborted'
endif

       call h5dwrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1188 <core_hgroup:write_dset_i_i8_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dwrite_f'
   stop 'Program aborted'
endif

       call h5dclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1188 <core_hgroup:write_dset_i_i8_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dclose_f'
   stop 'Program aborted'
endif

    else

       call h5screate_simple_f( 3,  INT(SHAPE(data), HSIZE_T),  dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1188 <core_hgroup:write_dset_i_i8_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5screate_simple_f'
   stop 'Program aborted'
endif

       call h5pcreate_f( H5P_DATASET_CREATE_F,  plist_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1188 <core_hgroup:write_dset_i_i8_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pcreate_f'
   stop 'Program aborted'
endif

       if (PRESENT(comp_level)) then
          call h5pset_chunk_f( plist_id,  3,  INT(SHAPE(data), HSIZE_T), hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1188 <core_hgroup:write_dset_i_i8_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pset_chunk_f'
   stop 'Program aborted'
endif

          call h5pset_deflate_f( plist_id,  comp_level, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1188 <core_hgroup:write_dset_i_i8_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pset_deflate_f'
   stop 'Program aborted'
endif

       endif

       mem_type_id = integer_mem_type_(I8)
       file_type_id = integer_file_type_(I8)

       data_ptr = C_LOC(data)

       call h5dcreate_f( hg%group_id,  item_name,  file_type_id,  dspace_id,  item_id, hdf_err,  dcpl_id=plist_id)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1188 <core_hgroup:write_dset_i_i8_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dcreate_f'
   stop 'Program aborted'
endif

       call h5dwrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1188 <core_hgroup:write_dset_i_i8_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dwrite_f'
   stop 'Program aborted'
endif

       call h5pclose_f( plist_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1188 <core_hgroup:write_dset_i_i8_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pclose_f'
   stop 'Program aborted'
endif

       call h5dclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1188 <core_hgroup:write_dset_i_i8_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dclose_f'
   stop 'Program aborted'
endif

       call h5sclose_f( dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1188 <core_hgroup:write_dset_i_i8_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5sclose_f'
   stop 'Program aborted'
endif

    endif

    ! Finish

    return

  end subroutine write_dset_i_i8_3_

  subroutine write_dset_i_i8_4_ (hg, item_name, data, overwrite, comp_level)

    type(hgroup_t), intent(inout)              :: hg
    character(*), intent(in)                   :: item_name
    integer(I8), target, intent(in) :: data(:,:,:,:)

    contiguous :: data

    logical, intent(in), optional              :: overwrite
    integer, intent(in), optional              :: comp_level

    logical          :: overwrite_
    integer          :: hdf_err
    logical          :: item_exists
    integer(HSIZE_T) :: item_shape(4)
    integer(HID_T)   :: dspace_id
    integer(HID_T)   :: plist_id
    integer(HID_T)   :: mem_type_id
    integer(HID_T)   :: file_type_id
    type(C_PTR)      :: data_ptr
    integer(HID_T)   :: item_id

    if(PRESENT(overwrite)) then
       overwrite_ = overwrite
    else
       overwrite_ = .FALSE.
    endif

    ! Write the item

    item_exists = dset_exists_(hg, item_name)
    overwrite_ = overwrite_ .AND. item_exists

    if(overwrite_) then

       call get_dset_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 1189 <core_hgroup:write_dset_i_i8' &
 & // '_4_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

       mem_type_id = integer_mem_type_(I8)

       data_ptr = C_LOC(data)

       call h5dopen_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1189 <core_hgroup:write_dset_i_i8_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dopen_f'
   stop 'Program aborted'
endif

       call h5dwrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1189 <core_hgroup:write_dset_i_i8_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dwrite_f'
   stop 'Program aborted'
endif

       call h5dclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1189 <core_hgroup:write_dset_i_i8_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dclose_f'
   stop 'Program aborted'
endif

    else

       call h5screate_simple_f( 4,  INT(SHAPE(data), HSIZE_T),  dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1189 <core_hgroup:write_dset_i_i8_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5screate_simple_f'
   stop 'Program aborted'
endif

       call h5pcreate_f( H5P_DATASET_CREATE_F,  plist_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1189 <core_hgroup:write_dset_i_i8_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pcreate_f'
   stop 'Program aborted'
endif

       if (PRESENT(comp_level)) then
          call h5pset_chunk_f( plist_id,  4,  INT(SHAPE(data), HSIZE_T), hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1189 <core_hgroup:write_dset_i_i8_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pset_chunk_f'
   stop 'Program aborted'
endif

          call h5pset_deflate_f( plist_id,  comp_level, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1189 <core_hgroup:write_dset_i_i8_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pset_deflate_f'
   stop 'Program aborted'
endif

       endif

       mem_type_id = integer_mem_type_(I8)
       file_type_id = integer_file_type_(I8)

       data_ptr = C_LOC(data)

       call h5dcreate_f( hg%group_id,  item_name,  file_type_id,  dspace_id,  item_id, hdf_err,  dcpl_id=plist_id)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1189 <core_hgroup:write_dset_i_i8_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dcreate_f'
   stop 'Program aborted'
endif

       call h5dwrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1189 <core_hgroup:write_dset_i_i8_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dwrite_f'
   stop 'Program aborted'
endif

       call h5pclose_f( plist_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1189 <core_hgroup:write_dset_i_i8_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pclose_f'
   stop 'Program aborted'
endif

       call h5dclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1189 <core_hgroup:write_dset_i_i8_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dclose_f'
   stop 'Program aborted'
endif

       call h5sclose_f( dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1189 <core_hgroup:write_dset_i_i8_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5sclose_f'
   stop 'Program aborted'
endif

    endif

    ! Finish

    return

  end subroutine write_dset_i_i8_4_

  subroutine write_dset_r_sp_0_ (hg, item_name, data, overwrite, comp_level)

    type(hgroup_t), intent(inout)              :: hg
    character(*), intent(in)                   :: item_name
    real(SP), target, intent(in) :: data

    logical, intent(in), optional              :: overwrite
    integer, intent(in), optional              :: comp_level

    logical          :: overwrite_
    integer          :: hdf_err
    logical          :: item_exists
    integer(HSIZE_T) :: item_shape(0)
    integer(HID_T)   :: dspace_id
    integer(HID_T)   :: plist_id
    integer(HID_T)   :: mem_type_id
    integer(HID_T)   :: file_type_id
    type(C_PTR)      :: data_ptr
    integer(HID_T)   :: item_id

    if(PRESENT(overwrite)) then
       overwrite_ = overwrite
    else
       overwrite_ = .FALSE.
    endif

    ! Write the item

    item_exists = dset_exists_(hg, item_name)
    overwrite_ = overwrite_ .AND. item_exists

    if(overwrite_) then

       call get_dset_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 1191 <core_hgroup:write_dset_r_sp' &
 & // '_0_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

       mem_type_id = real_mem_type_(SP)

       data_ptr = C_LOC(data)

       call h5dopen_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1191 <core_hgroup:write_dset_r_sp_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dopen_f'
   stop 'Program aborted'
endif

       call h5dwrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1191 <core_hgroup:write_dset_r_sp_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dwrite_f'
   stop 'Program aborted'
endif

       call h5dclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1191 <core_hgroup:write_dset_r_sp_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dclose_f'
   stop 'Program aborted'
endif

    else

       call h5screate_f( H5S_SCALAR_F,  dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1191 <core_hgroup:write_dset_r_sp_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5screate_f'
   stop 'Program aborted'
endif

       call h5pcreate_f( H5P_DATASET_CREATE_F,  plist_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1191 <core_hgroup:write_dset_r_sp_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pcreate_f'
   stop 'Program aborted'
endif

       if (PRESENT(comp_level)) then
          call h5pset_chunk_f( plist_id,  0,  INT(SHAPE(data), HSIZE_T), hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1191 <core_hgroup:write_dset_r_sp_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pset_chunk_f'
   stop 'Program aborted'
endif

          call h5pset_deflate_f( plist_id,  comp_level, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1191 <core_hgroup:write_dset_r_sp_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pset_deflate_f'
   stop 'Program aborted'
endif

       endif

       mem_type_id = real_mem_type_(SP)
       file_type_id = real_file_type_(SP)

       data_ptr = C_LOC(data)

       call h5dcreate_f( hg%group_id,  item_name,  file_type_id,  dspace_id,  item_id, hdf_err,  dcpl_id=plist_id)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1191 <core_hgroup:write_dset_r_sp_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dcreate_f'
   stop 'Program aborted'
endif

       call h5dwrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1191 <core_hgroup:write_dset_r_sp_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dwrite_f'
   stop 'Program aborted'
endif

       call h5pclose_f( plist_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1191 <core_hgroup:write_dset_r_sp_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pclose_f'
   stop 'Program aborted'
endif

       call h5dclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1191 <core_hgroup:write_dset_r_sp_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dclose_f'
   stop 'Program aborted'
endif

       call h5sclose_f( dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1191 <core_hgroup:write_dset_r_sp_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5sclose_f'
   stop 'Program aborted'
endif

    endif

    ! Finish

    return

  end subroutine write_dset_r_sp_0_

  subroutine write_dset_r_sp_1_ (hg, item_name, data, overwrite, comp_level)

    type(hgroup_t), intent(inout)              :: hg
    character(*), intent(in)                   :: item_name
    real(SP), target, intent(in) :: data(:)

    contiguous :: data

    logical, intent(in), optional              :: overwrite
    integer, intent(in), optional              :: comp_level

    logical          :: overwrite_
    integer          :: hdf_err
    logical          :: item_exists
    integer(HSIZE_T) :: item_shape(1)
    integer(HID_T)   :: dspace_id
    integer(HID_T)   :: plist_id
    integer(HID_T)   :: mem_type_id
    integer(HID_T)   :: file_type_id
    type(C_PTR)      :: data_ptr
    integer(HID_T)   :: item_id

    if(PRESENT(overwrite)) then
       overwrite_ = overwrite
    else
       overwrite_ = .FALSE.
    endif

    ! Write the item

    item_exists = dset_exists_(hg, item_name)
    overwrite_ = overwrite_ .AND. item_exists

    if(overwrite_) then

       call get_dset_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 1192 <core_hgroup:write_dset_r_sp' &
 & // '_1_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

       mem_type_id = real_mem_type_(SP)

       data_ptr = C_LOC(data)

       call h5dopen_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1192 <core_hgroup:write_dset_r_sp_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dopen_f'
   stop 'Program aborted'
endif

       call h5dwrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1192 <core_hgroup:write_dset_r_sp_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dwrite_f'
   stop 'Program aborted'
endif

       call h5dclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1192 <core_hgroup:write_dset_r_sp_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dclose_f'
   stop 'Program aborted'
endif

    else

       call h5screate_simple_f( 1,  INT(SHAPE(data), HSIZE_T),  dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1192 <core_hgroup:write_dset_r_sp_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5screate_simple_f'
   stop 'Program aborted'
endif

       call h5pcreate_f( H5P_DATASET_CREATE_F,  plist_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1192 <core_hgroup:write_dset_r_sp_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pcreate_f'
   stop 'Program aborted'
endif

       if (PRESENT(comp_level)) then
          call h5pset_chunk_f( plist_id,  1,  INT(SHAPE(data), HSIZE_T), hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1192 <core_hgroup:write_dset_r_sp_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pset_chunk_f'
   stop 'Program aborted'
endif

          call h5pset_deflate_f( plist_id,  comp_level, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1192 <core_hgroup:write_dset_r_sp_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pset_deflate_f'
   stop 'Program aborted'
endif

       endif

       mem_type_id = real_mem_type_(SP)
       file_type_id = real_file_type_(SP)

       data_ptr = C_LOC(data)

       call h5dcreate_f( hg%group_id,  item_name,  file_type_id,  dspace_id,  item_id, hdf_err,  dcpl_id=plist_id)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1192 <core_hgroup:write_dset_r_sp_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dcreate_f'
   stop 'Program aborted'
endif

       call h5dwrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1192 <core_hgroup:write_dset_r_sp_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dwrite_f'
   stop 'Program aborted'
endif

       call h5pclose_f( plist_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1192 <core_hgroup:write_dset_r_sp_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pclose_f'
   stop 'Program aborted'
endif

       call h5dclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1192 <core_hgroup:write_dset_r_sp_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dclose_f'
   stop 'Program aborted'
endif

       call h5sclose_f( dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1192 <core_hgroup:write_dset_r_sp_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5sclose_f'
   stop 'Program aborted'
endif

    endif

    ! Finish

    return

  end subroutine write_dset_r_sp_1_

  subroutine write_dset_r_sp_2_ (hg, item_name, data, overwrite, comp_level)

    type(hgroup_t), intent(inout)              :: hg
    character(*), intent(in)                   :: item_name
    real(SP), target, intent(in) :: data(:,:)

    contiguous :: data

    logical, intent(in), optional              :: overwrite
    integer, intent(in), optional              :: comp_level

    logical          :: overwrite_
    integer          :: hdf_err
    logical          :: item_exists
    integer(HSIZE_T) :: item_shape(2)
    integer(HID_T)   :: dspace_id
    integer(HID_T)   :: plist_id
    integer(HID_T)   :: mem_type_id
    integer(HID_T)   :: file_type_id
    type(C_PTR)      :: data_ptr
    integer(HID_T)   :: item_id

    if(PRESENT(overwrite)) then
       overwrite_ = overwrite
    else
       overwrite_ = .FALSE.
    endif

    ! Write the item

    item_exists = dset_exists_(hg, item_name)
    overwrite_ = overwrite_ .AND. item_exists

    if(overwrite_) then

       call get_dset_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 1193 <core_hgroup:write_dset_r_sp' &
 & // '_2_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

       mem_type_id = real_mem_type_(SP)

       data_ptr = C_LOC(data)

       call h5dopen_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1193 <core_hgroup:write_dset_r_sp_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dopen_f'
   stop 'Program aborted'
endif

       call h5dwrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1193 <core_hgroup:write_dset_r_sp_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dwrite_f'
   stop 'Program aborted'
endif

       call h5dclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1193 <core_hgroup:write_dset_r_sp_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dclose_f'
   stop 'Program aborted'
endif

    else

       call h5screate_simple_f( 2,  INT(SHAPE(data), HSIZE_T),  dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1193 <core_hgroup:write_dset_r_sp_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5screate_simple_f'
   stop 'Program aborted'
endif

       call h5pcreate_f( H5P_DATASET_CREATE_F,  plist_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1193 <core_hgroup:write_dset_r_sp_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pcreate_f'
   stop 'Program aborted'
endif

       if (PRESENT(comp_level)) then
          call h5pset_chunk_f( plist_id,  2,  INT(SHAPE(data), HSIZE_T), hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1193 <core_hgroup:write_dset_r_sp_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pset_chunk_f'
   stop 'Program aborted'
endif

          call h5pset_deflate_f( plist_id,  comp_level, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1193 <core_hgroup:write_dset_r_sp_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pset_deflate_f'
   stop 'Program aborted'
endif

       endif

       mem_type_id = real_mem_type_(SP)
       file_type_id = real_file_type_(SP)

       data_ptr = C_LOC(data)

       call h5dcreate_f( hg%group_id,  item_name,  file_type_id,  dspace_id,  item_id, hdf_err,  dcpl_id=plist_id)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1193 <core_hgroup:write_dset_r_sp_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dcreate_f'
   stop 'Program aborted'
endif

       call h5dwrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1193 <core_hgroup:write_dset_r_sp_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dwrite_f'
   stop 'Program aborted'
endif

       call h5pclose_f( plist_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1193 <core_hgroup:write_dset_r_sp_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pclose_f'
   stop 'Program aborted'
endif

       call h5dclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1193 <core_hgroup:write_dset_r_sp_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dclose_f'
   stop 'Program aborted'
endif

       call h5sclose_f( dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1193 <core_hgroup:write_dset_r_sp_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5sclose_f'
   stop 'Program aborted'
endif

    endif

    ! Finish

    return

  end subroutine write_dset_r_sp_2_

  subroutine write_dset_r_sp_3_ (hg, item_name, data, overwrite, comp_level)

    type(hgroup_t), intent(inout)              :: hg
    character(*), intent(in)                   :: item_name
    real(SP), target, intent(in) :: data(:,:,:)

    contiguous :: data

    logical, intent(in), optional              :: overwrite
    integer, intent(in), optional              :: comp_level

    logical          :: overwrite_
    integer          :: hdf_err
    logical          :: item_exists
    integer(HSIZE_T) :: item_shape(3)
    integer(HID_T)   :: dspace_id
    integer(HID_T)   :: plist_id
    integer(HID_T)   :: mem_type_id
    integer(HID_T)   :: file_type_id
    type(C_PTR)      :: data_ptr
    integer(HID_T)   :: item_id

    if(PRESENT(overwrite)) then
       overwrite_ = overwrite
    else
       overwrite_ = .FALSE.
    endif

    ! Write the item

    item_exists = dset_exists_(hg, item_name)
    overwrite_ = overwrite_ .AND. item_exists

    if(overwrite_) then

       call get_dset_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 1194 <core_hgroup:write_dset_r_sp' &
 & // '_3_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

       mem_type_id = real_mem_type_(SP)

       data_ptr = C_LOC(data)

       call h5dopen_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1194 <core_hgroup:write_dset_r_sp_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dopen_f'
   stop 'Program aborted'
endif

       call h5dwrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1194 <core_hgroup:write_dset_r_sp_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dwrite_f'
   stop 'Program aborted'
endif

       call h5dclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1194 <core_hgroup:write_dset_r_sp_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dclose_f'
   stop 'Program aborted'
endif

    else

       call h5screate_simple_f( 3,  INT(SHAPE(data), HSIZE_T),  dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1194 <core_hgroup:write_dset_r_sp_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5screate_simple_f'
   stop 'Program aborted'
endif

       call h5pcreate_f( H5P_DATASET_CREATE_F,  plist_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1194 <core_hgroup:write_dset_r_sp_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pcreate_f'
   stop 'Program aborted'
endif

       if (PRESENT(comp_level)) then
          call h5pset_chunk_f( plist_id,  3,  INT(SHAPE(data), HSIZE_T), hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1194 <core_hgroup:write_dset_r_sp_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pset_chunk_f'
   stop 'Program aborted'
endif

          call h5pset_deflate_f( plist_id,  comp_level, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1194 <core_hgroup:write_dset_r_sp_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pset_deflate_f'
   stop 'Program aborted'
endif

       endif

       mem_type_id = real_mem_type_(SP)
       file_type_id = real_file_type_(SP)

       data_ptr = C_LOC(data)

       call h5dcreate_f( hg%group_id,  item_name,  file_type_id,  dspace_id,  item_id, hdf_err,  dcpl_id=plist_id)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1194 <core_hgroup:write_dset_r_sp_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dcreate_f'
   stop 'Program aborted'
endif

       call h5dwrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1194 <core_hgroup:write_dset_r_sp_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dwrite_f'
   stop 'Program aborted'
endif

       call h5pclose_f( plist_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1194 <core_hgroup:write_dset_r_sp_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pclose_f'
   stop 'Program aborted'
endif

       call h5dclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1194 <core_hgroup:write_dset_r_sp_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dclose_f'
   stop 'Program aborted'
endif

       call h5sclose_f( dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1194 <core_hgroup:write_dset_r_sp_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5sclose_f'
   stop 'Program aborted'
endif

    endif

    ! Finish

    return

  end subroutine write_dset_r_sp_3_

  subroutine write_dset_r_sp_4_ (hg, item_name, data, overwrite, comp_level)

    type(hgroup_t), intent(inout)              :: hg
    character(*), intent(in)                   :: item_name
    real(SP), target, intent(in) :: data(:,:,:,:)

    contiguous :: data

    logical, intent(in), optional              :: overwrite
    integer, intent(in), optional              :: comp_level

    logical          :: overwrite_
    integer          :: hdf_err
    logical          :: item_exists
    integer(HSIZE_T) :: item_shape(4)
    integer(HID_T)   :: dspace_id
    integer(HID_T)   :: plist_id
    integer(HID_T)   :: mem_type_id
    integer(HID_T)   :: file_type_id
    type(C_PTR)      :: data_ptr
    integer(HID_T)   :: item_id

    if(PRESENT(overwrite)) then
       overwrite_ = overwrite
    else
       overwrite_ = .FALSE.
    endif

    ! Write the item

    item_exists = dset_exists_(hg, item_name)
    overwrite_ = overwrite_ .AND. item_exists

    if(overwrite_) then

       call get_dset_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 1195 <core_hgroup:write_dset_r_sp' &
 & // '_4_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

       mem_type_id = real_mem_type_(SP)

       data_ptr = C_LOC(data)

       call h5dopen_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1195 <core_hgroup:write_dset_r_sp_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dopen_f'
   stop 'Program aborted'
endif

       call h5dwrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1195 <core_hgroup:write_dset_r_sp_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dwrite_f'
   stop 'Program aborted'
endif

       call h5dclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1195 <core_hgroup:write_dset_r_sp_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dclose_f'
   stop 'Program aborted'
endif

    else

       call h5screate_simple_f( 4,  INT(SHAPE(data), HSIZE_T),  dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1195 <core_hgroup:write_dset_r_sp_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5screate_simple_f'
   stop 'Program aborted'
endif

       call h5pcreate_f( H5P_DATASET_CREATE_F,  plist_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1195 <core_hgroup:write_dset_r_sp_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pcreate_f'
   stop 'Program aborted'
endif

       if (PRESENT(comp_level)) then
          call h5pset_chunk_f( plist_id,  4,  INT(SHAPE(data), HSIZE_T), hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1195 <core_hgroup:write_dset_r_sp_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pset_chunk_f'
   stop 'Program aborted'
endif

          call h5pset_deflate_f( plist_id,  comp_level, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1195 <core_hgroup:write_dset_r_sp_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pset_deflate_f'
   stop 'Program aborted'
endif

       endif

       mem_type_id = real_mem_type_(SP)
       file_type_id = real_file_type_(SP)

       data_ptr = C_LOC(data)

       call h5dcreate_f( hg%group_id,  item_name,  file_type_id,  dspace_id,  item_id, hdf_err,  dcpl_id=plist_id)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1195 <core_hgroup:write_dset_r_sp_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dcreate_f'
   stop 'Program aborted'
endif

       call h5dwrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1195 <core_hgroup:write_dset_r_sp_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dwrite_f'
   stop 'Program aborted'
endif

       call h5pclose_f( plist_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1195 <core_hgroup:write_dset_r_sp_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pclose_f'
   stop 'Program aborted'
endif

       call h5dclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1195 <core_hgroup:write_dset_r_sp_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dclose_f'
   stop 'Program aborted'
endif

       call h5sclose_f( dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1195 <core_hgroup:write_dset_r_sp_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5sclose_f'
   stop 'Program aborted'
endif

    endif

    ! Finish

    return

  end subroutine write_dset_r_sp_4_

  subroutine write_dset_r_dp_0_ (hg, item_name, data, overwrite, comp_level)

    type(hgroup_t), intent(inout)              :: hg
    character(*), intent(in)                   :: item_name
    real(DP), target, intent(in) :: data

    logical, intent(in), optional              :: overwrite
    integer, intent(in), optional              :: comp_level

    logical          :: overwrite_
    integer          :: hdf_err
    logical          :: item_exists
    integer(HSIZE_T) :: item_shape(0)
    integer(HID_T)   :: dspace_id
    integer(HID_T)   :: plist_id
    integer(HID_T)   :: mem_type_id
    integer(HID_T)   :: file_type_id
    type(C_PTR)      :: data_ptr
    integer(HID_T)   :: item_id

    if(PRESENT(overwrite)) then
       overwrite_ = overwrite
    else
       overwrite_ = .FALSE.
    endif

    ! Write the item

    item_exists = dset_exists_(hg, item_name)
    overwrite_ = overwrite_ .AND. item_exists

    if(overwrite_) then

       call get_dset_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 1197 <core_hgroup:write_dset_r_dp' &
 & // '_0_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

       mem_type_id = real_mem_type_(DP)

       data_ptr = C_LOC(data)

       call h5dopen_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1197 <core_hgroup:write_dset_r_dp_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dopen_f'
   stop 'Program aborted'
endif

       call h5dwrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1197 <core_hgroup:write_dset_r_dp_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dwrite_f'
   stop 'Program aborted'
endif

       call h5dclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1197 <core_hgroup:write_dset_r_dp_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dclose_f'
   stop 'Program aborted'
endif

    else

       call h5screate_f( H5S_SCALAR_F,  dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1197 <core_hgroup:write_dset_r_dp_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5screate_f'
   stop 'Program aborted'
endif

       call h5pcreate_f( H5P_DATASET_CREATE_F,  plist_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1197 <core_hgroup:write_dset_r_dp_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pcreate_f'
   stop 'Program aborted'
endif

       if (PRESENT(comp_level)) then
          call h5pset_chunk_f( plist_id,  0,  INT(SHAPE(data), HSIZE_T), hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1197 <core_hgroup:write_dset_r_dp_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pset_chunk_f'
   stop 'Program aborted'
endif

          call h5pset_deflate_f( plist_id,  comp_level, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1197 <core_hgroup:write_dset_r_dp_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pset_deflate_f'
   stop 'Program aborted'
endif

       endif

       mem_type_id = real_mem_type_(DP)
       file_type_id = real_file_type_(DP)

       data_ptr = C_LOC(data)

       call h5dcreate_f( hg%group_id,  item_name,  file_type_id,  dspace_id,  item_id, hdf_err,  dcpl_id=plist_id)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1197 <core_hgroup:write_dset_r_dp_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dcreate_f'
   stop 'Program aborted'
endif

       call h5dwrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1197 <core_hgroup:write_dset_r_dp_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dwrite_f'
   stop 'Program aborted'
endif

       call h5pclose_f( plist_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1197 <core_hgroup:write_dset_r_dp_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pclose_f'
   stop 'Program aborted'
endif

       call h5dclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1197 <core_hgroup:write_dset_r_dp_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dclose_f'
   stop 'Program aborted'
endif

       call h5sclose_f( dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1197 <core_hgroup:write_dset_r_dp_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5sclose_f'
   stop 'Program aborted'
endif

    endif

    ! Finish

    return

  end subroutine write_dset_r_dp_0_

  subroutine write_dset_r_dp_1_ (hg, item_name, data, overwrite, comp_level)

    type(hgroup_t), intent(inout)              :: hg
    character(*), intent(in)                   :: item_name
    real(DP), target, intent(in) :: data(:)

    contiguous :: data

    logical, intent(in), optional              :: overwrite
    integer, intent(in), optional              :: comp_level

    logical          :: overwrite_
    integer          :: hdf_err
    logical          :: item_exists
    integer(HSIZE_T) :: item_shape(1)
    integer(HID_T)   :: dspace_id
    integer(HID_T)   :: plist_id
    integer(HID_T)   :: mem_type_id
    integer(HID_T)   :: file_type_id
    type(C_PTR)      :: data_ptr
    integer(HID_T)   :: item_id

    if(PRESENT(overwrite)) then
       overwrite_ = overwrite
    else
       overwrite_ = .FALSE.
    endif

    ! Write the item

    item_exists = dset_exists_(hg, item_name)
    overwrite_ = overwrite_ .AND. item_exists

    if(overwrite_) then

       call get_dset_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 1198 <core_hgroup:write_dset_r_dp' &
 & // '_1_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

       mem_type_id = real_mem_type_(DP)

       data_ptr = C_LOC(data)

       call h5dopen_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1198 <core_hgroup:write_dset_r_dp_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dopen_f'
   stop 'Program aborted'
endif

       call h5dwrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1198 <core_hgroup:write_dset_r_dp_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dwrite_f'
   stop 'Program aborted'
endif

       call h5dclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1198 <core_hgroup:write_dset_r_dp_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dclose_f'
   stop 'Program aborted'
endif

    else

       call h5screate_simple_f( 1,  INT(SHAPE(data), HSIZE_T),  dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1198 <core_hgroup:write_dset_r_dp_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5screate_simple_f'
   stop 'Program aborted'
endif

       call h5pcreate_f( H5P_DATASET_CREATE_F,  plist_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1198 <core_hgroup:write_dset_r_dp_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pcreate_f'
   stop 'Program aborted'
endif

       if (PRESENT(comp_level)) then
          call h5pset_chunk_f( plist_id,  1,  INT(SHAPE(data), HSIZE_T), hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1198 <core_hgroup:write_dset_r_dp_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pset_chunk_f'
   stop 'Program aborted'
endif

          call h5pset_deflate_f( plist_id,  comp_level, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1198 <core_hgroup:write_dset_r_dp_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pset_deflate_f'
   stop 'Program aborted'
endif

       endif

       mem_type_id = real_mem_type_(DP)
       file_type_id = real_file_type_(DP)

       data_ptr = C_LOC(data)

       call h5dcreate_f( hg%group_id,  item_name,  file_type_id,  dspace_id,  item_id, hdf_err,  dcpl_id=plist_id)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1198 <core_hgroup:write_dset_r_dp_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dcreate_f'
   stop 'Program aborted'
endif

       call h5dwrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1198 <core_hgroup:write_dset_r_dp_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dwrite_f'
   stop 'Program aborted'
endif

       call h5pclose_f( plist_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1198 <core_hgroup:write_dset_r_dp_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pclose_f'
   stop 'Program aborted'
endif

       call h5dclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1198 <core_hgroup:write_dset_r_dp_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dclose_f'
   stop 'Program aborted'
endif

       call h5sclose_f( dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1198 <core_hgroup:write_dset_r_dp_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5sclose_f'
   stop 'Program aborted'
endif

    endif

    ! Finish

    return

  end subroutine write_dset_r_dp_1_

  subroutine write_dset_r_dp_2_ (hg, item_name, data, overwrite, comp_level)

    type(hgroup_t), intent(inout)              :: hg
    character(*), intent(in)                   :: item_name
    real(DP), target, intent(in) :: data(:,:)

    contiguous :: data

    logical, intent(in), optional              :: overwrite
    integer, intent(in), optional              :: comp_level

    logical          :: overwrite_
    integer          :: hdf_err
    logical          :: item_exists
    integer(HSIZE_T) :: item_shape(2)
    integer(HID_T)   :: dspace_id
    integer(HID_T)   :: plist_id
    integer(HID_T)   :: mem_type_id
    integer(HID_T)   :: file_type_id
    type(C_PTR)      :: data_ptr
    integer(HID_T)   :: item_id

    if(PRESENT(overwrite)) then
       overwrite_ = overwrite
    else
       overwrite_ = .FALSE.
    endif

    ! Write the item

    item_exists = dset_exists_(hg, item_name)
    overwrite_ = overwrite_ .AND. item_exists

    if(overwrite_) then

       call get_dset_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 1199 <core_hgroup:write_dset_r_dp' &
 & // '_2_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

       mem_type_id = real_mem_type_(DP)

       data_ptr = C_LOC(data)

       call h5dopen_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1199 <core_hgroup:write_dset_r_dp_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dopen_f'
   stop 'Program aborted'
endif

       call h5dwrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1199 <core_hgroup:write_dset_r_dp_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dwrite_f'
   stop 'Program aborted'
endif

       call h5dclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1199 <core_hgroup:write_dset_r_dp_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dclose_f'
   stop 'Program aborted'
endif

    else

       call h5screate_simple_f( 2,  INT(SHAPE(data), HSIZE_T),  dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1199 <core_hgroup:write_dset_r_dp_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5screate_simple_f'
   stop 'Program aborted'
endif

       call h5pcreate_f( H5P_DATASET_CREATE_F,  plist_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1199 <core_hgroup:write_dset_r_dp_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pcreate_f'
   stop 'Program aborted'
endif

       if (PRESENT(comp_level)) then
          call h5pset_chunk_f( plist_id,  2,  INT(SHAPE(data), HSIZE_T), hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1199 <core_hgroup:write_dset_r_dp_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pset_chunk_f'
   stop 'Program aborted'
endif

          call h5pset_deflate_f( plist_id,  comp_level, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1199 <core_hgroup:write_dset_r_dp_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pset_deflate_f'
   stop 'Program aborted'
endif

       endif

       mem_type_id = real_mem_type_(DP)
       file_type_id = real_file_type_(DP)

       data_ptr = C_LOC(data)

       call h5dcreate_f( hg%group_id,  item_name,  file_type_id,  dspace_id,  item_id, hdf_err,  dcpl_id=plist_id)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1199 <core_hgroup:write_dset_r_dp_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dcreate_f'
   stop 'Program aborted'
endif

       call h5dwrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1199 <core_hgroup:write_dset_r_dp_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dwrite_f'
   stop 'Program aborted'
endif

       call h5pclose_f( plist_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1199 <core_hgroup:write_dset_r_dp_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pclose_f'
   stop 'Program aborted'
endif

       call h5dclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1199 <core_hgroup:write_dset_r_dp_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dclose_f'
   stop 'Program aborted'
endif

       call h5sclose_f( dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1199 <core_hgroup:write_dset_r_dp_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5sclose_f'
   stop 'Program aborted'
endif

    endif

    ! Finish

    return

  end subroutine write_dset_r_dp_2_

  subroutine write_dset_r_dp_3_ (hg, item_name, data, overwrite, comp_level)

    type(hgroup_t), intent(inout)              :: hg
    character(*), intent(in)                   :: item_name
    real(DP), target, intent(in) :: data(:,:,:)

    contiguous :: data

    logical, intent(in), optional              :: overwrite
    integer, intent(in), optional              :: comp_level

    logical          :: overwrite_
    integer          :: hdf_err
    logical          :: item_exists
    integer(HSIZE_T) :: item_shape(3)
    integer(HID_T)   :: dspace_id
    integer(HID_T)   :: plist_id
    integer(HID_T)   :: mem_type_id
    integer(HID_T)   :: file_type_id
    type(C_PTR)      :: data_ptr
    integer(HID_T)   :: item_id

    if(PRESENT(overwrite)) then
       overwrite_ = overwrite
    else
       overwrite_ = .FALSE.
    endif

    ! Write the item

    item_exists = dset_exists_(hg, item_name)
    overwrite_ = overwrite_ .AND. item_exists

    if(overwrite_) then

       call get_dset_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 1200 <core_hgroup:write_dset_r_dp' &
 & // '_3_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

       mem_type_id = real_mem_type_(DP)

       data_ptr = C_LOC(data)

       call h5dopen_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1200 <core_hgroup:write_dset_r_dp_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dopen_f'
   stop 'Program aborted'
endif

       call h5dwrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1200 <core_hgroup:write_dset_r_dp_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dwrite_f'
   stop 'Program aborted'
endif

       call h5dclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1200 <core_hgroup:write_dset_r_dp_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dclose_f'
   stop 'Program aborted'
endif

    else

       call h5screate_simple_f( 3,  INT(SHAPE(data), HSIZE_T),  dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1200 <core_hgroup:write_dset_r_dp_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5screate_simple_f'
   stop 'Program aborted'
endif

       call h5pcreate_f( H5P_DATASET_CREATE_F,  plist_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1200 <core_hgroup:write_dset_r_dp_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pcreate_f'
   stop 'Program aborted'
endif

       if (PRESENT(comp_level)) then
          call h5pset_chunk_f( plist_id,  3,  INT(SHAPE(data), HSIZE_T), hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1200 <core_hgroup:write_dset_r_dp_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pset_chunk_f'
   stop 'Program aborted'
endif

          call h5pset_deflate_f( plist_id,  comp_level, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1200 <core_hgroup:write_dset_r_dp_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pset_deflate_f'
   stop 'Program aborted'
endif

       endif

       mem_type_id = real_mem_type_(DP)
       file_type_id = real_file_type_(DP)

       data_ptr = C_LOC(data)

       call h5dcreate_f( hg%group_id,  item_name,  file_type_id,  dspace_id,  item_id, hdf_err,  dcpl_id=plist_id)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1200 <core_hgroup:write_dset_r_dp_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dcreate_f'
   stop 'Program aborted'
endif

       call h5dwrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1200 <core_hgroup:write_dset_r_dp_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dwrite_f'
   stop 'Program aborted'
endif

       call h5pclose_f( plist_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1200 <core_hgroup:write_dset_r_dp_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pclose_f'
   stop 'Program aborted'
endif

       call h5dclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1200 <core_hgroup:write_dset_r_dp_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dclose_f'
   stop 'Program aborted'
endif

       call h5sclose_f( dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1200 <core_hgroup:write_dset_r_dp_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5sclose_f'
   stop 'Program aborted'
endif

    endif

    ! Finish

    return

  end subroutine write_dset_r_dp_3_

  subroutine write_dset_r_dp_4_ (hg, item_name, data, overwrite, comp_level)

    type(hgroup_t), intent(inout)              :: hg
    character(*), intent(in)                   :: item_name
    real(DP), target, intent(in) :: data(:,:,:,:)

    contiguous :: data

    logical, intent(in), optional              :: overwrite
    integer, intent(in), optional              :: comp_level

    logical          :: overwrite_
    integer          :: hdf_err
    logical          :: item_exists
    integer(HSIZE_T) :: item_shape(4)
    integer(HID_T)   :: dspace_id
    integer(HID_T)   :: plist_id
    integer(HID_T)   :: mem_type_id
    integer(HID_T)   :: file_type_id
    type(C_PTR)      :: data_ptr
    integer(HID_T)   :: item_id

    if(PRESENT(overwrite)) then
       overwrite_ = overwrite
    else
       overwrite_ = .FALSE.
    endif

    ! Write the item

    item_exists = dset_exists_(hg, item_name)
    overwrite_ = overwrite_ .AND. item_exists

    if(overwrite_) then

       call get_dset_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 1201 <core_hgroup:write_dset_r_dp' &
 & // '_4_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

       mem_type_id = real_mem_type_(DP)

       data_ptr = C_LOC(data)

       call h5dopen_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1201 <core_hgroup:write_dset_r_dp_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dopen_f'
   stop 'Program aborted'
endif

       call h5dwrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1201 <core_hgroup:write_dset_r_dp_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dwrite_f'
   stop 'Program aborted'
endif

       call h5dclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1201 <core_hgroup:write_dset_r_dp_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dclose_f'
   stop 'Program aborted'
endif

    else

       call h5screate_simple_f( 4,  INT(SHAPE(data), HSIZE_T),  dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1201 <core_hgroup:write_dset_r_dp_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5screate_simple_f'
   stop 'Program aborted'
endif

       call h5pcreate_f( H5P_DATASET_CREATE_F,  plist_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1201 <core_hgroup:write_dset_r_dp_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pcreate_f'
   stop 'Program aborted'
endif

       if (PRESENT(comp_level)) then
          call h5pset_chunk_f( plist_id,  4,  INT(SHAPE(data), HSIZE_T), hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1201 <core_hgroup:write_dset_r_dp_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pset_chunk_f'
   stop 'Program aborted'
endif

          call h5pset_deflate_f( plist_id,  comp_level, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1201 <core_hgroup:write_dset_r_dp_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pset_deflate_f'
   stop 'Program aborted'
endif

       endif

       mem_type_id = real_mem_type_(DP)
       file_type_id = real_file_type_(DP)

       data_ptr = C_LOC(data)

       call h5dcreate_f( hg%group_id,  item_name,  file_type_id,  dspace_id,  item_id, hdf_err,  dcpl_id=plist_id)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1201 <core_hgroup:write_dset_r_dp_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dcreate_f'
   stop 'Program aborted'
endif

       call h5dwrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1201 <core_hgroup:write_dset_r_dp_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dwrite_f'
   stop 'Program aborted'
endif

       call h5pclose_f( plist_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1201 <core_hgroup:write_dset_r_dp_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pclose_f'
   stop 'Program aborted'
endif

       call h5dclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1201 <core_hgroup:write_dset_r_dp_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dclose_f'
   stop 'Program aborted'
endif

       call h5sclose_f( dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1201 <core_hgroup:write_dset_r_dp_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5sclose_f'
   stop 'Program aborted'
endif

    endif

    ! Finish

    return

  end subroutine write_dset_r_dp_4_

  subroutine write_dset_c_sp_0_ (hg, item_name, data, overwrite, comp_level)

    type(hgroup_t), intent(inout)              :: hg
    character(*), intent(in)                   :: item_name
    complex(SP), target, intent(in) :: data

    logical, intent(in), optional              :: overwrite
    integer, intent(in), optional              :: comp_level

    logical          :: overwrite_
    integer          :: hdf_err
    logical          :: item_exists
    integer(HSIZE_T) :: item_shape(0)
    integer(HID_T)   :: dspace_id
    integer(HID_T)   :: plist_id
    integer(HID_T)   :: mem_type_id
    integer(HID_T)   :: file_type_id
    type(C_PTR)      :: data_ptr
    integer(HID_T)   :: item_id

    if(PRESENT(overwrite)) then
       overwrite_ = overwrite
    else
       overwrite_ = .FALSE.
    endif

    ! Write the item

    item_exists = dset_exists_(hg, item_name)
    overwrite_ = overwrite_ .AND. item_exists

    if(overwrite_) then

       call get_dset_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 1203 <core_hgroup:write_dset_c_sp' &
 & // '_0_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

       mem_type_id = complex_mem_type_(SP)

       data_ptr = C_LOC(data)

       call h5dopen_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1203 <core_hgroup:write_dset_c_sp_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dopen_f'
   stop 'Program aborted'
endif

       call h5dwrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1203 <core_hgroup:write_dset_c_sp_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dwrite_f'
   stop 'Program aborted'
endif

       call h5dclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1203 <core_hgroup:write_dset_c_sp_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dclose_f'
   stop 'Program aborted'
endif

    else

       call h5screate_f( H5S_SCALAR_F,  dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1203 <core_hgroup:write_dset_c_sp_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5screate_f'
   stop 'Program aborted'
endif

       call h5pcreate_f( H5P_DATASET_CREATE_F,  plist_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1203 <core_hgroup:write_dset_c_sp_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pcreate_f'
   stop 'Program aborted'
endif

       if (PRESENT(comp_level)) then
          call h5pset_chunk_f( plist_id,  0,  INT(SHAPE(data), HSIZE_T), hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1203 <core_hgroup:write_dset_c_sp_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pset_chunk_f'
   stop 'Program aborted'
endif

          call h5pset_deflate_f( plist_id,  comp_level, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1203 <core_hgroup:write_dset_c_sp_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pset_deflate_f'
   stop 'Program aborted'
endif

       endif

       mem_type_id = complex_mem_type_(SP)
       file_type_id = complex_file_type_(SP)

       data_ptr = C_LOC(data)

       call h5dcreate_f( hg%group_id,  item_name,  file_type_id,  dspace_id,  item_id, hdf_err,  dcpl_id=plist_id)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1203 <core_hgroup:write_dset_c_sp_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dcreate_f'
   stop 'Program aborted'
endif

       call h5dwrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1203 <core_hgroup:write_dset_c_sp_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dwrite_f'
   stop 'Program aborted'
endif

       call h5pclose_f( plist_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1203 <core_hgroup:write_dset_c_sp_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pclose_f'
   stop 'Program aborted'
endif

       call h5dclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1203 <core_hgroup:write_dset_c_sp_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dclose_f'
   stop 'Program aborted'
endif

       call h5sclose_f( dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1203 <core_hgroup:write_dset_c_sp_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5sclose_f'
   stop 'Program aborted'
endif

    endif

    ! Finish

    return

  end subroutine write_dset_c_sp_0_

  subroutine write_dset_c_sp_1_ (hg, item_name, data, overwrite, comp_level)

    type(hgroup_t), intent(inout)              :: hg
    character(*), intent(in)                   :: item_name
    complex(SP), target, intent(in) :: data(:)

    contiguous :: data

    logical, intent(in), optional              :: overwrite
    integer, intent(in), optional              :: comp_level

    logical          :: overwrite_
    integer          :: hdf_err
    logical          :: item_exists
    integer(HSIZE_T) :: item_shape(1)
    integer(HID_T)   :: dspace_id
    integer(HID_T)   :: plist_id
    integer(HID_T)   :: mem_type_id
    integer(HID_T)   :: file_type_id
    type(C_PTR)      :: data_ptr
    integer(HID_T)   :: item_id

    if(PRESENT(overwrite)) then
       overwrite_ = overwrite
    else
       overwrite_ = .FALSE.
    endif

    ! Write the item

    item_exists = dset_exists_(hg, item_name)
    overwrite_ = overwrite_ .AND. item_exists

    if(overwrite_) then

       call get_dset_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 1204 <core_hgroup:write_dset_c_sp' &
 & // '_1_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

       mem_type_id = complex_mem_type_(SP)

       data_ptr = C_LOC(data)

       call h5dopen_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1204 <core_hgroup:write_dset_c_sp_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dopen_f'
   stop 'Program aborted'
endif

       call h5dwrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1204 <core_hgroup:write_dset_c_sp_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dwrite_f'
   stop 'Program aborted'
endif

       call h5dclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1204 <core_hgroup:write_dset_c_sp_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dclose_f'
   stop 'Program aborted'
endif

    else

       call h5screate_simple_f( 1,  INT(SHAPE(data), HSIZE_T),  dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1204 <core_hgroup:write_dset_c_sp_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5screate_simple_f'
   stop 'Program aborted'
endif

       call h5pcreate_f( H5P_DATASET_CREATE_F,  plist_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1204 <core_hgroup:write_dset_c_sp_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pcreate_f'
   stop 'Program aborted'
endif

       if (PRESENT(comp_level)) then
          call h5pset_chunk_f( plist_id,  1,  INT(SHAPE(data), HSIZE_T), hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1204 <core_hgroup:write_dset_c_sp_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pset_chunk_f'
   stop 'Program aborted'
endif

          call h5pset_deflate_f( plist_id,  comp_level, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1204 <core_hgroup:write_dset_c_sp_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pset_deflate_f'
   stop 'Program aborted'
endif

       endif

       mem_type_id = complex_mem_type_(SP)
       file_type_id = complex_file_type_(SP)

       data_ptr = C_LOC(data)

       call h5dcreate_f( hg%group_id,  item_name,  file_type_id,  dspace_id,  item_id, hdf_err,  dcpl_id=plist_id)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1204 <core_hgroup:write_dset_c_sp_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dcreate_f'
   stop 'Program aborted'
endif

       call h5dwrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1204 <core_hgroup:write_dset_c_sp_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dwrite_f'
   stop 'Program aborted'
endif

       call h5pclose_f( plist_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1204 <core_hgroup:write_dset_c_sp_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pclose_f'
   stop 'Program aborted'
endif

       call h5dclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1204 <core_hgroup:write_dset_c_sp_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dclose_f'
   stop 'Program aborted'
endif

       call h5sclose_f( dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1204 <core_hgroup:write_dset_c_sp_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5sclose_f'
   stop 'Program aborted'
endif

    endif

    ! Finish

    return

  end subroutine write_dset_c_sp_1_

  subroutine write_dset_c_sp_2_ (hg, item_name, data, overwrite, comp_level)

    type(hgroup_t), intent(inout)              :: hg
    character(*), intent(in)                   :: item_name
    complex(SP), target, intent(in) :: data(:,:)

    contiguous :: data

    logical, intent(in), optional              :: overwrite
    integer, intent(in), optional              :: comp_level

    logical          :: overwrite_
    integer          :: hdf_err
    logical          :: item_exists
    integer(HSIZE_T) :: item_shape(2)
    integer(HID_T)   :: dspace_id
    integer(HID_T)   :: plist_id
    integer(HID_T)   :: mem_type_id
    integer(HID_T)   :: file_type_id
    type(C_PTR)      :: data_ptr
    integer(HID_T)   :: item_id

    if(PRESENT(overwrite)) then
       overwrite_ = overwrite
    else
       overwrite_ = .FALSE.
    endif

    ! Write the item

    item_exists = dset_exists_(hg, item_name)
    overwrite_ = overwrite_ .AND. item_exists

    if(overwrite_) then

       call get_dset_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 1205 <core_hgroup:write_dset_c_sp' &
 & // '_2_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

       mem_type_id = complex_mem_type_(SP)

       data_ptr = C_LOC(data)

       call h5dopen_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1205 <core_hgroup:write_dset_c_sp_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dopen_f'
   stop 'Program aborted'
endif

       call h5dwrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1205 <core_hgroup:write_dset_c_sp_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dwrite_f'
   stop 'Program aborted'
endif

       call h5dclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1205 <core_hgroup:write_dset_c_sp_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dclose_f'
   stop 'Program aborted'
endif

    else

       call h5screate_simple_f( 2,  INT(SHAPE(data), HSIZE_T),  dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1205 <core_hgroup:write_dset_c_sp_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5screate_simple_f'
   stop 'Program aborted'
endif

       call h5pcreate_f( H5P_DATASET_CREATE_F,  plist_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1205 <core_hgroup:write_dset_c_sp_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pcreate_f'
   stop 'Program aborted'
endif

       if (PRESENT(comp_level)) then
          call h5pset_chunk_f( plist_id,  2,  INT(SHAPE(data), HSIZE_T), hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1205 <core_hgroup:write_dset_c_sp_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pset_chunk_f'
   stop 'Program aborted'
endif

          call h5pset_deflate_f( plist_id,  comp_level, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1205 <core_hgroup:write_dset_c_sp_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pset_deflate_f'
   stop 'Program aborted'
endif

       endif

       mem_type_id = complex_mem_type_(SP)
       file_type_id = complex_file_type_(SP)

       data_ptr = C_LOC(data)

       call h5dcreate_f( hg%group_id,  item_name,  file_type_id,  dspace_id,  item_id, hdf_err,  dcpl_id=plist_id)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1205 <core_hgroup:write_dset_c_sp_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dcreate_f'
   stop 'Program aborted'
endif

       call h5dwrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1205 <core_hgroup:write_dset_c_sp_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dwrite_f'
   stop 'Program aborted'
endif

       call h5pclose_f( plist_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1205 <core_hgroup:write_dset_c_sp_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pclose_f'
   stop 'Program aborted'
endif

       call h5dclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1205 <core_hgroup:write_dset_c_sp_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dclose_f'
   stop 'Program aborted'
endif

       call h5sclose_f( dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1205 <core_hgroup:write_dset_c_sp_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5sclose_f'
   stop 'Program aborted'
endif

    endif

    ! Finish

    return

  end subroutine write_dset_c_sp_2_

  subroutine write_dset_c_sp_3_ (hg, item_name, data, overwrite, comp_level)

    type(hgroup_t), intent(inout)              :: hg
    character(*), intent(in)                   :: item_name
    complex(SP), target, intent(in) :: data(:,:,:)

    contiguous :: data

    logical, intent(in), optional              :: overwrite
    integer, intent(in), optional              :: comp_level

    logical          :: overwrite_
    integer          :: hdf_err
    logical          :: item_exists
    integer(HSIZE_T) :: item_shape(3)
    integer(HID_T)   :: dspace_id
    integer(HID_T)   :: plist_id
    integer(HID_T)   :: mem_type_id
    integer(HID_T)   :: file_type_id
    type(C_PTR)      :: data_ptr
    integer(HID_T)   :: item_id

    if(PRESENT(overwrite)) then
       overwrite_ = overwrite
    else
       overwrite_ = .FALSE.
    endif

    ! Write the item

    item_exists = dset_exists_(hg, item_name)
    overwrite_ = overwrite_ .AND. item_exists

    if(overwrite_) then

       call get_dset_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 1206 <core_hgroup:write_dset_c_sp' &
 & // '_3_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

       mem_type_id = complex_mem_type_(SP)

       data_ptr = C_LOC(data)

       call h5dopen_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1206 <core_hgroup:write_dset_c_sp_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dopen_f'
   stop 'Program aborted'
endif

       call h5dwrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1206 <core_hgroup:write_dset_c_sp_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dwrite_f'
   stop 'Program aborted'
endif

       call h5dclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1206 <core_hgroup:write_dset_c_sp_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dclose_f'
   stop 'Program aborted'
endif

    else

       call h5screate_simple_f( 3,  INT(SHAPE(data), HSIZE_T),  dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1206 <core_hgroup:write_dset_c_sp_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5screate_simple_f'
   stop 'Program aborted'
endif

       call h5pcreate_f( H5P_DATASET_CREATE_F,  plist_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1206 <core_hgroup:write_dset_c_sp_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pcreate_f'
   stop 'Program aborted'
endif

       if (PRESENT(comp_level)) then
          call h5pset_chunk_f( plist_id,  3,  INT(SHAPE(data), HSIZE_T), hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1206 <core_hgroup:write_dset_c_sp_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pset_chunk_f'
   stop 'Program aborted'
endif

          call h5pset_deflate_f( plist_id,  comp_level, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1206 <core_hgroup:write_dset_c_sp_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pset_deflate_f'
   stop 'Program aborted'
endif

       endif

       mem_type_id = complex_mem_type_(SP)
       file_type_id = complex_file_type_(SP)

       data_ptr = C_LOC(data)

       call h5dcreate_f( hg%group_id,  item_name,  file_type_id,  dspace_id,  item_id, hdf_err,  dcpl_id=plist_id)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1206 <core_hgroup:write_dset_c_sp_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dcreate_f'
   stop 'Program aborted'
endif

       call h5dwrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1206 <core_hgroup:write_dset_c_sp_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dwrite_f'
   stop 'Program aborted'
endif

       call h5pclose_f( plist_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1206 <core_hgroup:write_dset_c_sp_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pclose_f'
   stop 'Program aborted'
endif

       call h5dclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1206 <core_hgroup:write_dset_c_sp_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dclose_f'
   stop 'Program aborted'
endif

       call h5sclose_f( dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1206 <core_hgroup:write_dset_c_sp_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5sclose_f'
   stop 'Program aborted'
endif

    endif

    ! Finish

    return

  end subroutine write_dset_c_sp_3_

  subroutine write_dset_c_sp_4_ (hg, item_name, data, overwrite, comp_level)

    type(hgroup_t), intent(inout)              :: hg
    character(*), intent(in)                   :: item_name
    complex(SP), target, intent(in) :: data(:,:,:,:)

    contiguous :: data

    logical, intent(in), optional              :: overwrite
    integer, intent(in), optional              :: comp_level

    logical          :: overwrite_
    integer          :: hdf_err
    logical          :: item_exists
    integer(HSIZE_T) :: item_shape(4)
    integer(HID_T)   :: dspace_id
    integer(HID_T)   :: plist_id
    integer(HID_T)   :: mem_type_id
    integer(HID_T)   :: file_type_id
    type(C_PTR)      :: data_ptr
    integer(HID_T)   :: item_id

    if(PRESENT(overwrite)) then
       overwrite_ = overwrite
    else
       overwrite_ = .FALSE.
    endif

    ! Write the item

    item_exists = dset_exists_(hg, item_name)
    overwrite_ = overwrite_ .AND. item_exists

    if(overwrite_) then

       call get_dset_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 1207 <core_hgroup:write_dset_c_sp' &
 & // '_4_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

       mem_type_id = complex_mem_type_(SP)

       data_ptr = C_LOC(data)

       call h5dopen_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1207 <core_hgroup:write_dset_c_sp_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dopen_f'
   stop 'Program aborted'
endif

       call h5dwrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1207 <core_hgroup:write_dset_c_sp_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dwrite_f'
   stop 'Program aborted'
endif

       call h5dclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1207 <core_hgroup:write_dset_c_sp_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dclose_f'
   stop 'Program aborted'
endif

    else

       call h5screate_simple_f( 4,  INT(SHAPE(data), HSIZE_T),  dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1207 <core_hgroup:write_dset_c_sp_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5screate_simple_f'
   stop 'Program aborted'
endif

       call h5pcreate_f( H5P_DATASET_CREATE_F,  plist_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1207 <core_hgroup:write_dset_c_sp_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pcreate_f'
   stop 'Program aborted'
endif

       if (PRESENT(comp_level)) then
          call h5pset_chunk_f( plist_id,  4,  INT(SHAPE(data), HSIZE_T), hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1207 <core_hgroup:write_dset_c_sp_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pset_chunk_f'
   stop 'Program aborted'
endif

          call h5pset_deflate_f( plist_id,  comp_level, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1207 <core_hgroup:write_dset_c_sp_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pset_deflate_f'
   stop 'Program aborted'
endif

       endif

       mem_type_id = complex_mem_type_(SP)
       file_type_id = complex_file_type_(SP)

       data_ptr = C_LOC(data)

       call h5dcreate_f( hg%group_id,  item_name,  file_type_id,  dspace_id,  item_id, hdf_err,  dcpl_id=plist_id)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1207 <core_hgroup:write_dset_c_sp_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dcreate_f'
   stop 'Program aborted'
endif

       call h5dwrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1207 <core_hgroup:write_dset_c_sp_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dwrite_f'
   stop 'Program aborted'
endif

       call h5pclose_f( plist_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1207 <core_hgroup:write_dset_c_sp_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pclose_f'
   stop 'Program aborted'
endif

       call h5dclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1207 <core_hgroup:write_dset_c_sp_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dclose_f'
   stop 'Program aborted'
endif

       call h5sclose_f( dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1207 <core_hgroup:write_dset_c_sp_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5sclose_f'
   stop 'Program aborted'
endif

    endif

    ! Finish

    return

  end subroutine write_dset_c_sp_4_

  subroutine write_dset_c_dp_0_ (hg, item_name, data, overwrite, comp_level)

    type(hgroup_t), intent(inout)              :: hg
    character(*), intent(in)                   :: item_name
    complex(DP), target, intent(in) :: data

    logical, intent(in), optional              :: overwrite
    integer, intent(in), optional              :: comp_level

    logical          :: overwrite_
    integer          :: hdf_err
    logical          :: item_exists
    integer(HSIZE_T) :: item_shape(0)
    integer(HID_T)   :: dspace_id
    integer(HID_T)   :: plist_id
    integer(HID_T)   :: mem_type_id
    integer(HID_T)   :: file_type_id
    type(C_PTR)      :: data_ptr
    integer(HID_T)   :: item_id

    if(PRESENT(overwrite)) then
       overwrite_ = overwrite
    else
       overwrite_ = .FALSE.
    endif

    ! Write the item

    item_exists = dset_exists_(hg, item_name)
    overwrite_ = overwrite_ .AND. item_exists

    if(overwrite_) then

       call get_dset_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 1209 <core_hgroup:write_dset_c_dp' &
 & // '_0_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

       mem_type_id = complex_mem_type_(DP)

       data_ptr = C_LOC(data)

       call h5dopen_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1209 <core_hgroup:write_dset_c_dp_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dopen_f'
   stop 'Program aborted'
endif

       call h5dwrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1209 <core_hgroup:write_dset_c_dp_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dwrite_f'
   stop 'Program aborted'
endif

       call h5dclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1209 <core_hgroup:write_dset_c_dp_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dclose_f'
   stop 'Program aborted'
endif

    else

       call h5screate_f( H5S_SCALAR_F,  dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1209 <core_hgroup:write_dset_c_dp_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5screate_f'
   stop 'Program aborted'
endif

       call h5pcreate_f( H5P_DATASET_CREATE_F,  plist_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1209 <core_hgroup:write_dset_c_dp_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pcreate_f'
   stop 'Program aborted'
endif

       if (PRESENT(comp_level)) then
          call h5pset_chunk_f( plist_id,  0,  INT(SHAPE(data), HSIZE_T), hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1209 <core_hgroup:write_dset_c_dp_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pset_chunk_f'
   stop 'Program aborted'
endif

          call h5pset_deflate_f( plist_id,  comp_level, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1209 <core_hgroup:write_dset_c_dp_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pset_deflate_f'
   stop 'Program aborted'
endif

       endif

       mem_type_id = complex_mem_type_(DP)
       file_type_id = complex_file_type_(DP)

       data_ptr = C_LOC(data)

       call h5dcreate_f( hg%group_id,  item_name,  file_type_id,  dspace_id,  item_id, hdf_err,  dcpl_id=plist_id)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1209 <core_hgroup:write_dset_c_dp_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dcreate_f'
   stop 'Program aborted'
endif

       call h5dwrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1209 <core_hgroup:write_dset_c_dp_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dwrite_f'
   stop 'Program aborted'
endif

       call h5pclose_f( plist_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1209 <core_hgroup:write_dset_c_dp_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pclose_f'
   stop 'Program aborted'
endif

       call h5dclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1209 <core_hgroup:write_dset_c_dp_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dclose_f'
   stop 'Program aborted'
endif

       call h5sclose_f( dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1209 <core_hgroup:write_dset_c_dp_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5sclose_f'
   stop 'Program aborted'
endif

    endif

    ! Finish

    return

  end subroutine write_dset_c_dp_0_

  subroutine write_dset_c_dp_1_ (hg, item_name, data, overwrite, comp_level)

    type(hgroup_t), intent(inout)              :: hg
    character(*), intent(in)                   :: item_name
    complex(DP), target, intent(in) :: data(:)

    contiguous :: data

    logical, intent(in), optional              :: overwrite
    integer, intent(in), optional              :: comp_level

    logical          :: overwrite_
    integer          :: hdf_err
    logical          :: item_exists
    integer(HSIZE_T) :: item_shape(1)
    integer(HID_T)   :: dspace_id
    integer(HID_T)   :: plist_id
    integer(HID_T)   :: mem_type_id
    integer(HID_T)   :: file_type_id
    type(C_PTR)      :: data_ptr
    integer(HID_T)   :: item_id

    if(PRESENT(overwrite)) then
       overwrite_ = overwrite
    else
       overwrite_ = .FALSE.
    endif

    ! Write the item

    item_exists = dset_exists_(hg, item_name)
    overwrite_ = overwrite_ .AND. item_exists

    if(overwrite_) then

       call get_dset_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 1210 <core_hgroup:write_dset_c_dp' &
 & // '_1_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

       mem_type_id = complex_mem_type_(DP)

       data_ptr = C_LOC(data)

       call h5dopen_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1210 <core_hgroup:write_dset_c_dp_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dopen_f'
   stop 'Program aborted'
endif

       call h5dwrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1210 <core_hgroup:write_dset_c_dp_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dwrite_f'
   stop 'Program aborted'
endif

       call h5dclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1210 <core_hgroup:write_dset_c_dp_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dclose_f'
   stop 'Program aborted'
endif

    else

       call h5screate_simple_f( 1,  INT(SHAPE(data), HSIZE_T),  dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1210 <core_hgroup:write_dset_c_dp_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5screate_simple_f'
   stop 'Program aborted'
endif

       call h5pcreate_f( H5P_DATASET_CREATE_F,  plist_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1210 <core_hgroup:write_dset_c_dp_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pcreate_f'
   stop 'Program aborted'
endif

       if (PRESENT(comp_level)) then
          call h5pset_chunk_f( plist_id,  1,  INT(SHAPE(data), HSIZE_T), hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1210 <core_hgroup:write_dset_c_dp_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pset_chunk_f'
   stop 'Program aborted'
endif

          call h5pset_deflate_f( plist_id,  comp_level, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1210 <core_hgroup:write_dset_c_dp_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pset_deflate_f'
   stop 'Program aborted'
endif

       endif

       mem_type_id = complex_mem_type_(DP)
       file_type_id = complex_file_type_(DP)

       data_ptr = C_LOC(data)

       call h5dcreate_f( hg%group_id,  item_name,  file_type_id,  dspace_id,  item_id, hdf_err,  dcpl_id=plist_id)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1210 <core_hgroup:write_dset_c_dp_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dcreate_f'
   stop 'Program aborted'
endif

       call h5dwrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1210 <core_hgroup:write_dset_c_dp_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dwrite_f'
   stop 'Program aborted'
endif

       call h5pclose_f( plist_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1210 <core_hgroup:write_dset_c_dp_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pclose_f'
   stop 'Program aborted'
endif

       call h5dclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1210 <core_hgroup:write_dset_c_dp_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dclose_f'
   stop 'Program aborted'
endif

       call h5sclose_f( dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1210 <core_hgroup:write_dset_c_dp_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5sclose_f'
   stop 'Program aborted'
endif

    endif

    ! Finish

    return

  end subroutine write_dset_c_dp_1_

  subroutine write_dset_c_dp_2_ (hg, item_name, data, overwrite, comp_level)

    type(hgroup_t), intent(inout)              :: hg
    character(*), intent(in)                   :: item_name
    complex(DP), target, intent(in) :: data(:,:)

    contiguous :: data

    logical, intent(in), optional              :: overwrite
    integer, intent(in), optional              :: comp_level

    logical          :: overwrite_
    integer          :: hdf_err
    logical          :: item_exists
    integer(HSIZE_T) :: item_shape(2)
    integer(HID_T)   :: dspace_id
    integer(HID_T)   :: plist_id
    integer(HID_T)   :: mem_type_id
    integer(HID_T)   :: file_type_id
    type(C_PTR)      :: data_ptr
    integer(HID_T)   :: item_id

    if(PRESENT(overwrite)) then
       overwrite_ = overwrite
    else
       overwrite_ = .FALSE.
    endif

    ! Write the item

    item_exists = dset_exists_(hg, item_name)
    overwrite_ = overwrite_ .AND. item_exists

    if(overwrite_) then

       call get_dset_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 1211 <core_hgroup:write_dset_c_dp' &
 & // '_2_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

       mem_type_id = complex_mem_type_(DP)

       data_ptr = C_LOC(data)

       call h5dopen_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1211 <core_hgroup:write_dset_c_dp_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dopen_f'
   stop 'Program aborted'
endif

       call h5dwrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1211 <core_hgroup:write_dset_c_dp_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dwrite_f'
   stop 'Program aborted'
endif

       call h5dclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1211 <core_hgroup:write_dset_c_dp_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dclose_f'
   stop 'Program aborted'
endif

    else

       call h5screate_simple_f( 2,  INT(SHAPE(data), HSIZE_T),  dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1211 <core_hgroup:write_dset_c_dp_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5screate_simple_f'
   stop 'Program aborted'
endif

       call h5pcreate_f( H5P_DATASET_CREATE_F,  plist_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1211 <core_hgroup:write_dset_c_dp_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pcreate_f'
   stop 'Program aborted'
endif

       if (PRESENT(comp_level)) then
          call h5pset_chunk_f( plist_id,  2,  INT(SHAPE(data), HSIZE_T), hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1211 <core_hgroup:write_dset_c_dp_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pset_chunk_f'
   stop 'Program aborted'
endif

          call h5pset_deflate_f( plist_id,  comp_level, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1211 <core_hgroup:write_dset_c_dp_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pset_deflate_f'
   stop 'Program aborted'
endif

       endif

       mem_type_id = complex_mem_type_(DP)
       file_type_id = complex_file_type_(DP)

       data_ptr = C_LOC(data)

       call h5dcreate_f( hg%group_id,  item_name,  file_type_id,  dspace_id,  item_id, hdf_err,  dcpl_id=plist_id)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1211 <core_hgroup:write_dset_c_dp_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dcreate_f'
   stop 'Program aborted'
endif

       call h5dwrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1211 <core_hgroup:write_dset_c_dp_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dwrite_f'
   stop 'Program aborted'
endif

       call h5pclose_f( plist_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1211 <core_hgroup:write_dset_c_dp_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pclose_f'
   stop 'Program aborted'
endif

       call h5dclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1211 <core_hgroup:write_dset_c_dp_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dclose_f'
   stop 'Program aborted'
endif

       call h5sclose_f( dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1211 <core_hgroup:write_dset_c_dp_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5sclose_f'
   stop 'Program aborted'
endif

    endif

    ! Finish

    return

  end subroutine write_dset_c_dp_2_

  subroutine write_dset_c_dp_3_ (hg, item_name, data, overwrite, comp_level)

    type(hgroup_t), intent(inout)              :: hg
    character(*), intent(in)                   :: item_name
    complex(DP), target, intent(in) :: data(:,:,:)

    contiguous :: data

    logical, intent(in), optional              :: overwrite
    integer, intent(in), optional              :: comp_level

    logical          :: overwrite_
    integer          :: hdf_err
    logical          :: item_exists
    integer(HSIZE_T) :: item_shape(3)
    integer(HID_T)   :: dspace_id
    integer(HID_T)   :: plist_id
    integer(HID_T)   :: mem_type_id
    integer(HID_T)   :: file_type_id
    type(C_PTR)      :: data_ptr
    integer(HID_T)   :: item_id

    if(PRESENT(overwrite)) then
       overwrite_ = overwrite
    else
       overwrite_ = .FALSE.
    endif

    ! Write the item

    item_exists = dset_exists_(hg, item_name)
    overwrite_ = overwrite_ .AND. item_exists

    if(overwrite_) then

       call get_dset_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 1212 <core_hgroup:write_dset_c_dp' &
 & // '_3_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

       mem_type_id = complex_mem_type_(DP)

       data_ptr = C_LOC(data)

       call h5dopen_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1212 <core_hgroup:write_dset_c_dp_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dopen_f'
   stop 'Program aborted'
endif

       call h5dwrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1212 <core_hgroup:write_dset_c_dp_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dwrite_f'
   stop 'Program aborted'
endif

       call h5dclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1212 <core_hgroup:write_dset_c_dp_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dclose_f'
   stop 'Program aborted'
endif

    else

       call h5screate_simple_f( 3,  INT(SHAPE(data), HSIZE_T),  dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1212 <core_hgroup:write_dset_c_dp_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5screate_simple_f'
   stop 'Program aborted'
endif

       call h5pcreate_f( H5P_DATASET_CREATE_F,  plist_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1212 <core_hgroup:write_dset_c_dp_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pcreate_f'
   stop 'Program aborted'
endif

       if (PRESENT(comp_level)) then
          call h5pset_chunk_f( plist_id,  3,  INT(SHAPE(data), HSIZE_T), hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1212 <core_hgroup:write_dset_c_dp_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pset_chunk_f'
   stop 'Program aborted'
endif

          call h5pset_deflate_f( plist_id,  comp_level, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1212 <core_hgroup:write_dset_c_dp_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pset_deflate_f'
   stop 'Program aborted'
endif

       endif

       mem_type_id = complex_mem_type_(DP)
       file_type_id = complex_file_type_(DP)

       data_ptr = C_LOC(data)

       call h5dcreate_f( hg%group_id,  item_name,  file_type_id,  dspace_id,  item_id, hdf_err,  dcpl_id=plist_id)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1212 <core_hgroup:write_dset_c_dp_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dcreate_f'
   stop 'Program aborted'
endif

       call h5dwrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1212 <core_hgroup:write_dset_c_dp_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dwrite_f'
   stop 'Program aborted'
endif

       call h5pclose_f( plist_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1212 <core_hgroup:write_dset_c_dp_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pclose_f'
   stop 'Program aborted'
endif

       call h5dclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1212 <core_hgroup:write_dset_c_dp_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dclose_f'
   stop 'Program aborted'
endif

       call h5sclose_f( dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1212 <core_hgroup:write_dset_c_dp_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5sclose_f'
   stop 'Program aborted'
endif

    endif

    ! Finish

    return

  end subroutine write_dset_c_dp_3_

  subroutine write_dset_c_dp_4_ (hg, item_name, data, overwrite, comp_level)

    type(hgroup_t), intent(inout)              :: hg
    character(*), intent(in)                   :: item_name
    complex(DP), target, intent(in) :: data(:,:,:,:)

    contiguous :: data

    logical, intent(in), optional              :: overwrite
    integer, intent(in), optional              :: comp_level

    logical          :: overwrite_
    integer          :: hdf_err
    logical          :: item_exists
    integer(HSIZE_T) :: item_shape(4)
    integer(HID_T)   :: dspace_id
    integer(HID_T)   :: plist_id
    integer(HID_T)   :: mem_type_id
    integer(HID_T)   :: file_type_id
    type(C_PTR)      :: data_ptr
    integer(HID_T)   :: item_id

    if(PRESENT(overwrite)) then
       overwrite_ = overwrite
    else
       overwrite_ = .FALSE.
    endif

    ! Write the item

    item_exists = dset_exists_(hg, item_name)
    overwrite_ = overwrite_ .AND. item_exists

    if(overwrite_) then

       call get_dset_shape(hg, item_name, item_shape)

    if(.NOT. (ALL(item_shape == SHAPE(data)))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''ALL(item_shape == SHAPE(data))'' failed at line 1213 <core_hgroup:write_dset_c_dp' &
 & // '_4_>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Array shape mismatch'
      stop
    endif

       mem_type_id = complex_mem_type_(DP)

       data_ptr = C_LOC(data)

       call h5dopen_f( hg%group_id,  item_name,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1213 <core_hgroup:write_dset_c_dp_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dopen_f'
   stop 'Program aborted'
endif

       call h5dwrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1213 <core_hgroup:write_dset_c_dp_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dwrite_f'
   stop 'Program aborted'
endif

       call h5dclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1213 <core_hgroup:write_dset_c_dp_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dclose_f'
   stop 'Program aborted'
endif

    else

       call h5screate_simple_f( 4,  INT(SHAPE(data), HSIZE_T),  dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1213 <core_hgroup:write_dset_c_dp_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5screate_simple_f'
   stop 'Program aborted'
endif

       call h5pcreate_f( H5P_DATASET_CREATE_F,  plist_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1213 <core_hgroup:write_dset_c_dp_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pcreate_f'
   stop 'Program aborted'
endif

       if (PRESENT(comp_level)) then
          call h5pset_chunk_f( plist_id,  4,  INT(SHAPE(data), HSIZE_T), hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1213 <core_hgroup:write_dset_c_dp_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pset_chunk_f'
   stop 'Program aborted'
endif

          call h5pset_deflate_f( plist_id,  comp_level, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1213 <core_hgroup:write_dset_c_dp_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pset_deflate_f'
   stop 'Program aborted'
endif

       endif

       mem_type_id = complex_mem_type_(DP)
       file_type_id = complex_file_type_(DP)

       data_ptr = C_LOC(data)

       call h5dcreate_f( hg%group_id,  item_name,  file_type_id,  dspace_id,  item_id, hdf_err,  dcpl_id=plist_id)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1213 <core_hgroup:write_dset_c_dp_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dcreate_f'
   stop 'Program aborted'
endif

       call h5dwrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1213 <core_hgroup:write_dset_c_dp_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dwrite_f'
   stop 'Program aborted'
endif

       call h5pclose_f( plist_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1213 <core_hgroup:write_dset_c_dp_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5pclose_f'
   stop 'Program aborted'
endif

       call h5dclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1213 <core_hgroup:write_dset_c_dp_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dclose_f'
   stop 'Program aborted'
endif

       call h5sclose_f( dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1213 <core_hgroup:write_dset_c_dp_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5sclose_f'
   stop 'Program aborted'
endif

    endif

    ! Finish

    return

  end subroutine write_dset_c_dp_4_

!****

  subroutine write_attr_a_0_ (hg, item_name, data)

    type(hgroup_t), intent(inout)    :: hg
    character(*), intent(in)         :: item_name
    character(*), target, intent(in) :: data

    integer        :: hdf_err
    integer(HID_T) :: mem_type_id
    integer(HID_T) :: file_type_id
    integer(HID_T) :: dspace_id
    type(C_PTR)    :: data_ptr
    integer(HID_T) :: item_id

    ! Write the character item

    call h5tcopy_f( H5T_NATIVE_CHARACTER,  mem_type_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1307 <core_hgroup:write_attr_a_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5tcopy_f'
   stop 'Program aborted'
endif

    call h5tset_size_f( mem_type_id,  LEN(data, SIZE_T), hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1307 <core_hgroup:write_attr_a_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5tset_size_f'
   stop 'Program aborted'
endif

    call h5tcopy_f( H5T_NATIVE_CHARACTER,  file_type_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1307 <core_hgroup:write_attr_a_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5tcopy_f'
   stop 'Program aborted'
endif

    call h5tset_size_f( file_type_id,  LEN(data, SIZE_T), hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1307 <core_hgroup:write_attr_a_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5tset_size_f'
   stop 'Program aborted'
endif

    call h5screate_f( H5S_SCALAR_F,  dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1307 <core_hgroup:write_attr_a_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5screate_f'
   stop 'Program aborted'
endif

    data_ptr = C_LOC(data)

    call h5acreate_f( hg%group_id,  item_name,  file_type_id,  dspace_id,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1307 <core_hgroup:write_attr_a_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5acreate_f'
   stop 'Program aborted'
endif

    call h5awrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1307 <core_hgroup:write_attr_a_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5awrite_f'
   stop 'Program aborted'
endif

    call h5aclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1307 <core_hgroup:write_attr_a_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aclose_f'
   stop 'Program aborted'
endif

    call h5sclose_f( dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1307 <core_hgroup:write_attr_a_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5sclose_f'
   stop 'Program aborted'
endif

    call h5tclose_f( mem_type_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1307 <core_hgroup:write_attr_a_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5tclose_f'
   stop 'Program aborted'
endif

    call h5tclose_f( file_type_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1307 <core_hgroup:write_attr_a_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5tclose_f'
   stop 'Program aborted'
endif

    ! Finish

    return

  end subroutine write_attr_a_0_

!****

  subroutine write_attr_l_0_ (hg, item_name, data)

    type(hgroup_t), intent(inout) :: hg
    character(*), intent(in)      :: item_name
    logical, intent(in)           :: data

    ! Write the logical item

    call write_attr(hg, item_name, MERGE(1, 0, MASK=data))

    ! Finish

    return

  end subroutine write_attr_l_0_

  subroutine write_attr_a_1_ (hg, item_name, data)

    type(hgroup_t), intent(inout)    :: hg
    character(*), intent(in)         :: item_name
    character(*), target, intent(in) :: data(:)

    contiguous :: data

    integer        :: hdf_err
    integer(HID_T) :: mem_type_id
    integer(HID_T) :: file_type_id
    integer(HID_T) :: dspace_id
    type(C_PTR)    :: data_ptr
    integer(HID_T) :: item_id

    ! Write the character item

    call h5tcopy_f( H5T_NATIVE_CHARACTER,  mem_type_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1308 <core_hgroup:write_attr_a_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5tcopy_f'
   stop 'Program aborted'
endif

    call h5tset_size_f( mem_type_id,  LEN(data, SIZE_T), hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1308 <core_hgroup:write_attr_a_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5tset_size_f'
   stop 'Program aborted'
endif

    call h5tcopy_f( H5T_NATIVE_CHARACTER,  file_type_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1308 <core_hgroup:write_attr_a_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5tcopy_f'
   stop 'Program aborted'
endif

    call h5tset_size_f( file_type_id,  LEN(data, SIZE_T), hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1308 <core_hgroup:write_attr_a_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5tset_size_f'
   stop 'Program aborted'
endif

    call h5screate_simple_f( 1,  INT(SHAPE(data), HSIZE_T),  dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1308 <core_hgroup:write_attr_a_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5screate_simple_f'
   stop 'Program aborted'
endif

    data_ptr = C_LOC(data)

    call h5acreate_f( hg%group_id,  item_name,  file_type_id,  dspace_id,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1308 <core_hgroup:write_attr_a_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5acreate_f'
   stop 'Program aborted'
endif

    call h5awrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1308 <core_hgroup:write_attr_a_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5awrite_f'
   stop 'Program aborted'
endif

    call h5aclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1308 <core_hgroup:write_attr_a_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aclose_f'
   stop 'Program aborted'
endif

    call h5sclose_f( dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1308 <core_hgroup:write_attr_a_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5sclose_f'
   stop 'Program aborted'
endif

    call h5tclose_f( mem_type_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1308 <core_hgroup:write_attr_a_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5tclose_f'
   stop 'Program aborted'
endif

    call h5tclose_f( file_type_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1308 <core_hgroup:write_attr_a_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5tclose_f'
   stop 'Program aborted'
endif

    ! Finish

    return

  end subroutine write_attr_a_1_

!****

  subroutine write_attr_l_1_ (hg, item_name, data)

    type(hgroup_t), intent(inout) :: hg
    character(*), intent(in)      :: item_name
    logical, intent(in)           :: data(:)

    ! Write the logical item

    call write_attr(hg, item_name, MERGE(1, 0, MASK=data))

    ! Finish

    return

  end subroutine write_attr_l_1_

  subroutine write_attr_a_2_ (hg, item_name, data)

    type(hgroup_t), intent(inout)    :: hg
    character(*), intent(in)         :: item_name
    character(*), target, intent(in) :: data(:,:)

    contiguous :: data

    integer        :: hdf_err
    integer(HID_T) :: mem_type_id
    integer(HID_T) :: file_type_id
    integer(HID_T) :: dspace_id
    type(C_PTR)    :: data_ptr
    integer(HID_T) :: item_id

    ! Write the character item

    call h5tcopy_f( H5T_NATIVE_CHARACTER,  mem_type_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1309 <core_hgroup:write_attr_a_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5tcopy_f'
   stop 'Program aborted'
endif

    call h5tset_size_f( mem_type_id,  LEN(data, SIZE_T), hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1309 <core_hgroup:write_attr_a_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5tset_size_f'
   stop 'Program aborted'
endif

    call h5tcopy_f( H5T_NATIVE_CHARACTER,  file_type_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1309 <core_hgroup:write_attr_a_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5tcopy_f'
   stop 'Program aborted'
endif

    call h5tset_size_f( file_type_id,  LEN(data, SIZE_T), hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1309 <core_hgroup:write_attr_a_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5tset_size_f'
   stop 'Program aborted'
endif

    call h5screate_simple_f( 2,  INT(SHAPE(data), HSIZE_T),  dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1309 <core_hgroup:write_attr_a_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5screate_simple_f'
   stop 'Program aborted'
endif

    data_ptr = C_LOC(data)

    call h5acreate_f( hg%group_id,  item_name,  file_type_id,  dspace_id,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1309 <core_hgroup:write_attr_a_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5acreate_f'
   stop 'Program aborted'
endif

    call h5awrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1309 <core_hgroup:write_attr_a_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5awrite_f'
   stop 'Program aborted'
endif

    call h5aclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1309 <core_hgroup:write_attr_a_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aclose_f'
   stop 'Program aborted'
endif

    call h5sclose_f( dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1309 <core_hgroup:write_attr_a_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5sclose_f'
   stop 'Program aborted'
endif

    call h5tclose_f( mem_type_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1309 <core_hgroup:write_attr_a_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5tclose_f'
   stop 'Program aborted'
endif

    call h5tclose_f( file_type_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1309 <core_hgroup:write_attr_a_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5tclose_f'
   stop 'Program aborted'
endif

    ! Finish

    return

  end subroutine write_attr_a_2_

!****

  subroutine write_attr_l_2_ (hg, item_name, data)

    type(hgroup_t), intent(inout) :: hg
    character(*), intent(in)      :: item_name
    logical, intent(in)           :: data(:,:)

    ! Write the logical item

    call write_attr(hg, item_name, MERGE(1, 0, MASK=data))

    ! Finish

    return

  end subroutine write_attr_l_2_

  subroutine write_attr_a_3_ (hg, item_name, data)

    type(hgroup_t), intent(inout)    :: hg
    character(*), intent(in)         :: item_name
    character(*), target, intent(in) :: data(:,:,:)

    contiguous :: data

    integer        :: hdf_err
    integer(HID_T) :: mem_type_id
    integer(HID_T) :: file_type_id
    integer(HID_T) :: dspace_id
    type(C_PTR)    :: data_ptr
    integer(HID_T) :: item_id

    ! Write the character item

    call h5tcopy_f( H5T_NATIVE_CHARACTER,  mem_type_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1310 <core_hgroup:write_attr_a_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5tcopy_f'
   stop 'Program aborted'
endif

    call h5tset_size_f( mem_type_id,  LEN(data, SIZE_T), hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1310 <core_hgroup:write_attr_a_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5tset_size_f'
   stop 'Program aborted'
endif

    call h5tcopy_f( H5T_NATIVE_CHARACTER,  file_type_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1310 <core_hgroup:write_attr_a_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5tcopy_f'
   stop 'Program aborted'
endif

    call h5tset_size_f( file_type_id,  LEN(data, SIZE_T), hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1310 <core_hgroup:write_attr_a_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5tset_size_f'
   stop 'Program aborted'
endif

    call h5screate_simple_f( 3,  INT(SHAPE(data), HSIZE_T),  dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1310 <core_hgroup:write_attr_a_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5screate_simple_f'
   stop 'Program aborted'
endif

    data_ptr = C_LOC(data)

    call h5acreate_f( hg%group_id,  item_name,  file_type_id,  dspace_id,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1310 <core_hgroup:write_attr_a_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5acreate_f'
   stop 'Program aborted'
endif

    call h5awrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1310 <core_hgroup:write_attr_a_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5awrite_f'
   stop 'Program aborted'
endif

    call h5aclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1310 <core_hgroup:write_attr_a_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aclose_f'
   stop 'Program aborted'
endif

    call h5sclose_f( dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1310 <core_hgroup:write_attr_a_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5sclose_f'
   stop 'Program aborted'
endif

    call h5tclose_f( mem_type_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1310 <core_hgroup:write_attr_a_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5tclose_f'
   stop 'Program aborted'
endif

    call h5tclose_f( file_type_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1310 <core_hgroup:write_attr_a_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5tclose_f'
   stop 'Program aborted'
endif

    ! Finish

    return

  end subroutine write_attr_a_3_

!****

  subroutine write_attr_l_3_ (hg, item_name, data)

    type(hgroup_t), intent(inout) :: hg
    character(*), intent(in)      :: item_name
    logical, intent(in)           :: data(:,:,:)

    ! Write the logical item

    call write_attr(hg, item_name, MERGE(1, 0, MASK=data))

    ! Finish

    return

  end subroutine write_attr_l_3_

  subroutine write_attr_a_4_ (hg, item_name, data)

    type(hgroup_t), intent(inout)    :: hg
    character(*), intent(in)         :: item_name
    character(*), target, intent(in) :: data(:,:,:,:)

    contiguous :: data

    integer        :: hdf_err
    integer(HID_T) :: mem_type_id
    integer(HID_T) :: file_type_id
    integer(HID_T) :: dspace_id
    type(C_PTR)    :: data_ptr
    integer(HID_T) :: item_id

    ! Write the character item

    call h5tcopy_f( H5T_NATIVE_CHARACTER,  mem_type_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1311 <core_hgroup:write_attr_a_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5tcopy_f'
   stop 'Program aborted'
endif

    call h5tset_size_f( mem_type_id,  LEN(data, SIZE_T), hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1311 <core_hgroup:write_attr_a_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5tset_size_f'
   stop 'Program aborted'
endif

    call h5tcopy_f( H5T_NATIVE_CHARACTER,  file_type_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1311 <core_hgroup:write_attr_a_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5tcopy_f'
   stop 'Program aborted'
endif

    call h5tset_size_f( file_type_id,  LEN(data, SIZE_T), hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1311 <core_hgroup:write_attr_a_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5tset_size_f'
   stop 'Program aborted'
endif

    call h5screate_simple_f( 4,  INT(SHAPE(data), HSIZE_T),  dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1311 <core_hgroup:write_attr_a_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5screate_simple_f'
   stop 'Program aborted'
endif

    data_ptr = C_LOC(data)

    call h5acreate_f( hg%group_id,  item_name,  file_type_id,  dspace_id,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1311 <core_hgroup:write_attr_a_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5acreate_f'
   stop 'Program aborted'
endif

    call h5awrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1311 <core_hgroup:write_attr_a_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5awrite_f'
   stop 'Program aborted'
endif

    call h5aclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1311 <core_hgroup:write_attr_a_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aclose_f'
   stop 'Program aborted'
endif

    call h5sclose_f( dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1311 <core_hgroup:write_attr_a_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5sclose_f'
   stop 'Program aborted'
endif

    call h5tclose_f( mem_type_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1311 <core_hgroup:write_attr_a_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5tclose_f'
   stop 'Program aborted'
endif

    call h5tclose_f( file_type_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1311 <core_hgroup:write_attr_a_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5tclose_f'
   stop 'Program aborted'
endif

    ! Finish

    return

  end subroutine write_attr_a_4_

!****

  subroutine write_attr_l_4_ (hg, item_name, data)

    type(hgroup_t), intent(inout) :: hg
    character(*), intent(in)      :: item_name
    logical, intent(in)           :: data(:,:,:,:)

    ! Write the logical item

    call write_attr(hg, item_name, MERGE(1, 0, MASK=data))

    ! Finish

    return

  end subroutine write_attr_l_4_

  subroutine write_dset_a_0_ (hg, item_name, data)

    type(hgroup_t), intent(inout)    :: hg
    character(*), intent(in)         :: item_name
    character(*), target, intent(in) :: data

    integer        :: hdf_err
    integer(HID_T) :: mem_type_id
    integer(HID_T) :: file_type_id
    integer(HID_T) :: dspace_id
    type(C_PTR)    :: data_ptr
    integer(HID_T) :: item_id

    ! Write the character item

    call h5tcopy_f( H5T_NATIVE_CHARACTER,  mem_type_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1313 <core_hgroup:write_dset_a_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5tcopy_f'
   stop 'Program aborted'
endif

    call h5tset_size_f( mem_type_id,  LEN(data, SIZE_T), hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1313 <core_hgroup:write_dset_a_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5tset_size_f'
   stop 'Program aborted'
endif

    call h5tcopy_f( H5T_NATIVE_CHARACTER,  file_type_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1313 <core_hgroup:write_dset_a_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5tcopy_f'
   stop 'Program aborted'
endif

    call h5tset_size_f( file_type_id,  LEN(data, SIZE_T), hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1313 <core_hgroup:write_dset_a_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5tset_size_f'
   stop 'Program aborted'
endif

    call h5screate_f( H5S_SCALAR_F,  dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1313 <core_hgroup:write_dset_a_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5screate_f'
   stop 'Program aborted'
endif

    data_ptr = C_LOC(data)

    call h5dcreate_f( hg%group_id,  item_name,  file_type_id,  dspace_id,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1313 <core_hgroup:write_dset_a_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dcreate_f'
   stop 'Program aborted'
endif

    call h5dwrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1313 <core_hgroup:write_dset_a_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dwrite_f'
   stop 'Program aborted'
endif

    call h5dclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1313 <core_hgroup:write_dset_a_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dclose_f'
   stop 'Program aborted'
endif

    call h5sclose_f( dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1313 <core_hgroup:write_dset_a_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5sclose_f'
   stop 'Program aborted'
endif

    call h5tclose_f( mem_type_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1313 <core_hgroup:write_dset_a_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5tclose_f'
   stop 'Program aborted'
endif

    call h5tclose_f( file_type_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1313 <core_hgroup:write_dset_a_0_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5tclose_f'
   stop 'Program aborted'
endif

    ! Finish

    return

  end subroutine write_dset_a_0_

!****

  subroutine write_dset_l_0_ (hg, item_name, data)

    type(hgroup_t), intent(inout) :: hg
    character(*), intent(in)      :: item_name
    logical, intent(in)           :: data

    ! Write the logical item

    call write_dset(hg, item_name, MERGE(1, 0, MASK=data))

    ! Finish

    return

  end subroutine write_dset_l_0_

  subroutine write_dset_a_1_ (hg, item_name, data)

    type(hgroup_t), intent(inout)    :: hg
    character(*), intent(in)         :: item_name
    character(*), target, intent(in) :: data(:)

    contiguous :: data

    integer        :: hdf_err
    integer(HID_T) :: mem_type_id
    integer(HID_T) :: file_type_id
    integer(HID_T) :: dspace_id
    type(C_PTR)    :: data_ptr
    integer(HID_T) :: item_id

    ! Write the character item

    call h5tcopy_f( H5T_NATIVE_CHARACTER,  mem_type_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1314 <core_hgroup:write_dset_a_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5tcopy_f'
   stop 'Program aborted'
endif

    call h5tset_size_f( mem_type_id,  LEN(data, SIZE_T), hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1314 <core_hgroup:write_dset_a_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5tset_size_f'
   stop 'Program aborted'
endif

    call h5tcopy_f( H5T_NATIVE_CHARACTER,  file_type_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1314 <core_hgroup:write_dset_a_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5tcopy_f'
   stop 'Program aborted'
endif

    call h5tset_size_f( file_type_id,  LEN(data, SIZE_T), hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1314 <core_hgroup:write_dset_a_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5tset_size_f'
   stop 'Program aborted'
endif

    call h5screate_simple_f( 1,  INT(SHAPE(data), HSIZE_T),  dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1314 <core_hgroup:write_dset_a_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5screate_simple_f'
   stop 'Program aborted'
endif

    data_ptr = C_LOC(data)

    call h5dcreate_f( hg%group_id,  item_name,  file_type_id,  dspace_id,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1314 <core_hgroup:write_dset_a_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dcreate_f'
   stop 'Program aborted'
endif

    call h5dwrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1314 <core_hgroup:write_dset_a_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dwrite_f'
   stop 'Program aborted'
endif

    call h5dclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1314 <core_hgroup:write_dset_a_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dclose_f'
   stop 'Program aborted'
endif

    call h5sclose_f( dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1314 <core_hgroup:write_dset_a_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5sclose_f'
   stop 'Program aborted'
endif

    call h5tclose_f( mem_type_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1314 <core_hgroup:write_dset_a_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5tclose_f'
   stop 'Program aborted'
endif

    call h5tclose_f( file_type_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1314 <core_hgroup:write_dset_a_1_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5tclose_f'
   stop 'Program aborted'
endif

    ! Finish

    return

  end subroutine write_dset_a_1_

!****

  subroutine write_dset_l_1_ (hg, item_name, data)

    type(hgroup_t), intent(inout) :: hg
    character(*), intent(in)      :: item_name
    logical, intent(in)           :: data(:)

    ! Write the logical item

    call write_dset(hg, item_name, MERGE(1, 0, MASK=data))

    ! Finish

    return

  end subroutine write_dset_l_1_

  subroutine write_dset_a_2_ (hg, item_name, data)

    type(hgroup_t), intent(inout)    :: hg
    character(*), intent(in)         :: item_name
    character(*), target, intent(in) :: data(:,:)

    contiguous :: data

    integer        :: hdf_err
    integer(HID_T) :: mem_type_id
    integer(HID_T) :: file_type_id
    integer(HID_T) :: dspace_id
    type(C_PTR)    :: data_ptr
    integer(HID_T) :: item_id

    ! Write the character item

    call h5tcopy_f( H5T_NATIVE_CHARACTER,  mem_type_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1315 <core_hgroup:write_dset_a_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5tcopy_f'
   stop 'Program aborted'
endif

    call h5tset_size_f( mem_type_id,  LEN(data, SIZE_T), hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1315 <core_hgroup:write_dset_a_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5tset_size_f'
   stop 'Program aborted'
endif

    call h5tcopy_f( H5T_NATIVE_CHARACTER,  file_type_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1315 <core_hgroup:write_dset_a_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5tcopy_f'
   stop 'Program aborted'
endif

    call h5tset_size_f( file_type_id,  LEN(data, SIZE_T), hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1315 <core_hgroup:write_dset_a_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5tset_size_f'
   stop 'Program aborted'
endif

    call h5screate_simple_f( 2,  INT(SHAPE(data), HSIZE_T),  dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1315 <core_hgroup:write_dset_a_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5screate_simple_f'
   stop 'Program aborted'
endif

    data_ptr = C_LOC(data)

    call h5dcreate_f( hg%group_id,  item_name,  file_type_id,  dspace_id,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1315 <core_hgroup:write_dset_a_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dcreate_f'
   stop 'Program aborted'
endif

    call h5dwrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1315 <core_hgroup:write_dset_a_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dwrite_f'
   stop 'Program aborted'
endif

    call h5dclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1315 <core_hgroup:write_dset_a_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dclose_f'
   stop 'Program aborted'
endif

    call h5sclose_f( dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1315 <core_hgroup:write_dset_a_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5sclose_f'
   stop 'Program aborted'
endif

    call h5tclose_f( mem_type_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1315 <core_hgroup:write_dset_a_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5tclose_f'
   stop 'Program aborted'
endif

    call h5tclose_f( file_type_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1315 <core_hgroup:write_dset_a_2_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5tclose_f'
   stop 'Program aborted'
endif

    ! Finish

    return

  end subroutine write_dset_a_2_

!****

  subroutine write_dset_l_2_ (hg, item_name, data)

    type(hgroup_t), intent(inout) :: hg
    character(*), intent(in)      :: item_name
    logical, intent(in)           :: data(:,:)

    ! Write the logical item

    call write_dset(hg, item_name, MERGE(1, 0, MASK=data))

    ! Finish

    return

  end subroutine write_dset_l_2_

  subroutine write_dset_a_3_ (hg, item_name, data)

    type(hgroup_t), intent(inout)    :: hg
    character(*), intent(in)         :: item_name
    character(*), target, intent(in) :: data(:,:,:)

    contiguous :: data

    integer        :: hdf_err
    integer(HID_T) :: mem_type_id
    integer(HID_T) :: file_type_id
    integer(HID_T) :: dspace_id
    type(C_PTR)    :: data_ptr
    integer(HID_T) :: item_id

    ! Write the character item

    call h5tcopy_f( H5T_NATIVE_CHARACTER,  mem_type_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1316 <core_hgroup:write_dset_a_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5tcopy_f'
   stop 'Program aborted'
endif

    call h5tset_size_f( mem_type_id,  LEN(data, SIZE_T), hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1316 <core_hgroup:write_dset_a_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5tset_size_f'
   stop 'Program aborted'
endif

    call h5tcopy_f( H5T_NATIVE_CHARACTER,  file_type_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1316 <core_hgroup:write_dset_a_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5tcopy_f'
   stop 'Program aborted'
endif

    call h5tset_size_f( file_type_id,  LEN(data, SIZE_T), hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1316 <core_hgroup:write_dset_a_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5tset_size_f'
   stop 'Program aborted'
endif

    call h5screate_simple_f( 3,  INT(SHAPE(data), HSIZE_T),  dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1316 <core_hgroup:write_dset_a_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5screate_simple_f'
   stop 'Program aborted'
endif

    data_ptr = C_LOC(data)

    call h5dcreate_f( hg%group_id,  item_name,  file_type_id,  dspace_id,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1316 <core_hgroup:write_dset_a_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dcreate_f'
   stop 'Program aborted'
endif

    call h5dwrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1316 <core_hgroup:write_dset_a_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dwrite_f'
   stop 'Program aborted'
endif

    call h5dclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1316 <core_hgroup:write_dset_a_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dclose_f'
   stop 'Program aborted'
endif

    call h5sclose_f( dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1316 <core_hgroup:write_dset_a_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5sclose_f'
   stop 'Program aborted'
endif

    call h5tclose_f( mem_type_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1316 <core_hgroup:write_dset_a_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5tclose_f'
   stop 'Program aborted'
endif

    call h5tclose_f( file_type_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1316 <core_hgroup:write_dset_a_3_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5tclose_f'
   stop 'Program aborted'
endif

    ! Finish

    return

  end subroutine write_dset_a_3_

!****

  subroutine write_dset_l_3_ (hg, item_name, data)

    type(hgroup_t), intent(inout) :: hg
    character(*), intent(in)      :: item_name
    logical, intent(in)           :: data(:,:,:)

    ! Write the logical item

    call write_dset(hg, item_name, MERGE(1, 0, MASK=data))

    ! Finish

    return

  end subroutine write_dset_l_3_

  subroutine write_dset_a_4_ (hg, item_name, data)

    type(hgroup_t), intent(inout)    :: hg
    character(*), intent(in)         :: item_name
    character(*), target, intent(in) :: data(:,:,:,:)

    contiguous :: data

    integer        :: hdf_err
    integer(HID_T) :: mem_type_id
    integer(HID_T) :: file_type_id
    integer(HID_T) :: dspace_id
    type(C_PTR)    :: data_ptr
    integer(HID_T) :: item_id

    ! Write the character item

    call h5tcopy_f( H5T_NATIVE_CHARACTER,  mem_type_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1317 <core_hgroup:write_dset_a_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5tcopy_f'
   stop 'Program aborted'
endif

    call h5tset_size_f( mem_type_id,  LEN(data, SIZE_T), hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1317 <core_hgroup:write_dset_a_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5tset_size_f'
   stop 'Program aborted'
endif

    call h5tcopy_f( H5T_NATIVE_CHARACTER,  file_type_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1317 <core_hgroup:write_dset_a_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5tcopy_f'
   stop 'Program aborted'
endif

    call h5tset_size_f( file_type_id,  LEN(data, SIZE_T), hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1317 <core_hgroup:write_dset_a_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5tset_size_f'
   stop 'Program aborted'
endif

    call h5screate_simple_f( 4,  INT(SHAPE(data), HSIZE_T),  dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1317 <core_hgroup:write_dset_a_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5screate_simple_f'
   stop 'Program aborted'
endif

    data_ptr = C_LOC(data)

    call h5dcreate_f( hg%group_id,  item_name,  file_type_id,  dspace_id,  item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1317 <core_hgroup:write_dset_a_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dcreate_f'
   stop 'Program aborted'
endif

    call h5dwrite_f( item_id,  mem_type_id,  data_ptr, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1317 <core_hgroup:write_dset_a_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dwrite_f'
   stop 'Program aborted'
endif

    call h5dclose_f( item_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1317 <core_hgroup:write_dset_a_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dclose_f'
   stop 'Program aborted'
endif

    call h5sclose_f( dspace_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1317 <core_hgroup:write_dset_a_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5sclose_f'
   stop 'Program aborted'
endif

    call h5tclose_f( mem_type_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1317 <core_hgroup:write_dset_a_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5tclose_f'
   stop 'Program aborted'
endif

    call h5tclose_f( file_type_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1317 <core_hgroup:write_dset_a_4_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5tclose_f'
   stop 'Program aborted'
endif

    ! Finish

    return

  end subroutine write_dset_a_4_

!****

  subroutine write_dset_l_4_ (hg, item_name, data)

    type(hgroup_t), intent(inout) :: hg
    character(*), intent(in)      :: item_name
    logical, intent(in)           :: data(:,:,:,:)

    ! Write the logical item

    call write_dset(hg, item_name, MERGE(1, 0, MASK=data))

    ! Finish

    return

  end subroutine write_dset_l_4_

!****

  subroutine get_attr_shape (hg, attr_name, shape)

    type(hgroup_t), intent(inout) :: hg
    character(*), intent(in)      :: attr_name
    integer(HSIZE_T), intent(out) :: shape(:)

    integer          :: hdf_err
    integer(HID_T)   :: attr_id
    integer(HID_T)   :: space_id
    integer          :: rank
    integer(HSIZE_T) :: max_shape(SIZE(shape))

    ! Get the shape of the attribute

    call h5aopen_name_f( hg%group_id,  attr_name,  attr_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1335 <core_hgroup:get_attr_shape>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aopen_name_f'
   stop 'Program aborted'
endif

    call h5aget_space_f( attr_id,  space_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1337 <core_hgroup:get_attr_shape>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aget_space_f'
   stop 'Program aborted'
endif

    call h5sget_simple_extent_ndims_f( space_id,  rank, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1339 <core_hgroup:get_attr_shape>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5sget_simple_extent_ndims_f'
   stop 'Program aborted'
endif

    if(.NOT. (rank == SIZE(shape))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''rank == SIZE(shape)'' failed at line 1340 <core_hgroup:get_attr_shape>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Rank mismatch'
      stop
    endif

    call h5sget_simple_extent_dims_f( space_id,  shape,  max_shape, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1342 <core_hgroup:get_attr_shape>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5sget_simple_extent_dims_f'
   stop 'Program aborted'
endif

    call h5sclose_f( space_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1344 <core_hgroup:get_attr_shape>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5sclose_f'
   stop 'Program aborted'
endif

    call h5aclose_f( attr_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1345 <core_hgroup:get_attr_shape>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aclose_f'
   stop 'Program aborted'
endif

    ! Finish

    return

  end subroutine get_attr_shape

!****

  subroutine get_dset_shape (hg, dset_name, shape)

    type(hgroup_t), intent(inout) :: hg
    character(*), intent(in)      :: dset_name
    integer(HSIZE_T), intent(out) :: shape(:)

    integer          :: hdf_err
    integer(HID_T)   :: dset_id
    integer(HID_T)   :: space_id
    integer          :: rank
    integer(HSIZE_T) :: max_shape(SIZE(shape))

    ! Get the shape of the dataset

    call h5dopen_f( hg%group_id,  dset_name,  dset_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1369 <core_hgroup:get_dset_shape>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dopen_f'
   stop 'Program aborted'
endif

    call h5dget_space_f( dset_id,  space_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1371 <core_hgroup:get_dset_shape>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dget_space_f'
   stop 'Program aborted'
endif

    call h5sget_simple_extent_ndims_f( space_id,  rank, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1373 <core_hgroup:get_dset_shape>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5sget_simple_extent_ndims_f'
   stop 'Program aborted'
endif

    if(.NOT. (rank == SIZE(shape))) then
      write(UNIT=ERROR_UNIT, FMT=*) 'ASSERT ''rank == SIZE(shape)'' failed at line 1374 <core_hgroup:get_dset_shape>:'
      write(UNIT=ERROR_UNIT, FMT=*) 'Rank mismatch'
      stop
    endif

    call h5sget_simple_extent_dims_f( space_id,  shape,  max_shape, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1376 <core_hgroup:get_dset_shape>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5sget_simple_extent_dims_f'
   stop 'Program aborted'
endif

    call h5sclose_f( space_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1378 <core_hgroup:get_dset_shape>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5sclose_f'
   stop 'Program aborted'
endif

    call h5dclose_f( dset_id, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1379 <core_hgroup:get_dset_shape>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5dclose_f'
   stop 'Program aborted'
endif

    ! Finish

    return

  end subroutine get_dset_shape

!****

  function attr_exists_ (hg, attr_name) result (attr_exists)

    type(hgroup_t), intent(inout) :: hg
    character(*), intent(in)      :: attr_name
    logical                       :: attr_exists

    integer :: hdf_err

    ! Check if the attribute exists

    call h5aexists_f( hg%group_id,  attr_name,  attr_exists, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1399 <core_hgroup:attr_exists_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5aexists_f'
   stop 'Program aborted'
endif

    ! Finish

    return

  end function attr_exists_

!****

  function dset_exists_ (hg, dset_name) result (dset_exists)

    type(hgroup_t), intent(inout) :: hg
    character(*), intent(in)      :: dset_name
    logical                       :: dset_exists

    integer          :: hdf_err
    logical          :: link_exists
    type(h5o_info_t) :: obj_info

    ! Check if the dataset exists

    call h5lexists_f( hg%group_id,  dset_name,  link_exists, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1421 <core_hgroup:dset_exists_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5lexists_f'
   stop 'Program aborted'
endif

    if(link_exists) then

       call h5oget_info_by_name_f( hg%group_id,  dset_name,  obj_info, hdf_err)
if(hdf_err == -1) then
   call h5eprint_f (hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1425 <core_hgroup:dset_exists_>:'
   write(UNIT=ERROR_UNIT, FMT=*) 'Error in call to h5oget_info_by_name_f'
   stop 'Program aborted'
endif

       dset_exists = obj_info%type == H5O_TYPE_DATASET_F

    else

       dset_exists = .FALSE.

    endif

    ! Finish

    return

  end function dset_exists_

!****

  function integer_mem_type_ (kind) result (mem_type_id)

    integer, intent(in) :: kind
    integer(HID_T)      :: mem_type_id

    ! Determine the memory type for integers

    mem_type_id = h5kind_to_type(kind, H5_INTEGER_KIND)

    ! Finish

    return

  end function integer_mem_type_

!****

  function real_mem_type_ (kind) result (mem_type_id)

    integer, intent(in) :: kind
    integer(HID_T)      :: mem_type_id

    ! Determine the memory type for reals

    mem_type_id = h5kind_to_type(kind, H5_REAL_KIND)

    ! Finish

    return

  end function real_mem_type_

!****

  function complex_mem_type_ (kind) result (mem_type_id)

    integer, intent(in) :: kind
    integer(HID_T)      :: mem_type_id

    ! Determine the memory type for complexes

    select case (kind)
    case (SP)
       mem_type_id = c_sp_mem_type_id
    case (DP)
       mem_type_id = c_dp_mem_type_id
    case default

    write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1490 <core_hgroup:complex_mem_type_>:'
    write(UNIT=ERROR_UNIT, FMT=*) 'Unsupported kind'

  stop 'Program aborted'

    end select

    ! Finish

    return

  end function complex_mem_type_

!****

  function integer_file_type_ (kind) result (file_type_id)

    integer, intent(in) :: kind
    integer(HID_T)      :: file_type_id

    ! Determine the file type for integers

    select case (kind)
    case (I4)
       file_type_id = H5T_STD_I32LE
    case (I8)
       file_type_id = H5T_STD_I64LE
    case default

    write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1514 <core_hgroup:integer_file_type_>:'
    write(UNIT=ERROR_UNIT, FMT=*) 'Unsupported kind'

  stop 'Program aborted'

    end select

    ! Finish

    return

  end function integer_file_type_

!****

  function real_file_type_ (kind) result (file_type_id)

    integer, intent(in) :: kind
    integer(HID_T)      :: file_type_id

    ! Determine the file type for reals

    select case (kind)
    case (SP)
       file_type_id = H5T_IEEE_F32LE
    case (DP)
       file_type_id = H5T_IEEE_F64LE
    case default

    write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1538 <core_hgroup:real_file_type_>:'
    write(UNIT=ERROR_UNIT, FMT=*) 'Unsupported kind'

  stop 'Program aborted'

    end select

    ! Finish

    return

  end function real_file_type_

!****

  function complex_file_type_ (kind) result (file_type_id)

    integer, intent(in) :: kind
    integer(HID_T)      :: file_type_id

    ! Determine the file type for complexes

    select case (kind)
    case (SP)
       file_type_id = c_sp_file_type_id
    case (DP)
       file_type_id = c_dp_file_type_id
    case default

    write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 1562 <core_hgroup:complex_file_type_>:'
    write(UNIT=ERROR_UNIT, FMT=*) 'Unsupported kind'

  stop 'Program aborted'

    end select

    ! Finish

    return

  end function complex_file_type_

!****

  function elem_group_name (prefix, indices) result (group_name)

    character(*), intent(in)  :: prefix
    integer, intent(in)       :: indices(:)
    character(:), allocatable :: group_name

    integer                   :: n_indices
    character(:), allocatable :: name_format
    integer                   :: name_len
    integer                   :: i

    ! Set up an array-element group name

    n_indices = SIZE(indices)

    select case(n_indices)
    case(0)
       name_format = '(A,''()'')'
    case(1)
       name_format = '(A,''('',I0,'')'')'
    case default
       name_format = '(A,''('''//REPEAT('I0,'',''', n_indices-1)//'I0,'')'')'
    end select

    name_len = LEN_TRIM(prefix) + n_indices + 1

    do i = 1,SIZE(indices)
       if(indices(i) < 0) then
          name_len = name_len + FLOOR(LOG10(REAL(ABS(indices(i))))) + 2
       elseif(indices(i) > 0) then
          name_len = name_len + FLOOR(LOG10(REAL(indices(i)))) + 1
       else
          name_len = name_len + 1
       endif
    end do

    allocate(character(name_len) :: group_name)

    write(group_name, name_format) TRIM(prefix), indices

    ! Finish

    return

  end function elem_group_name

end module core_hgroup

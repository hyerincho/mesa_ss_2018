!fpx3_header(0.13_3a)
!
!dependencies
!   dir: ~/Documents/MESA_SS/gyre/src/build 
!   sources: -
!   includes: ../extern/core/core.inc
!   uses: core_kinds
!   provides: gyre_writer
!end dependencies
!
!end fpx3_header
! Module   : gyre_writer
! Purpose  : write data (interface)
!
! Copyright 2013 Rich Townsend
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

module gyre_writer

  ! Uses

  use core_kinds

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type, abstract :: writer_t
   contains
     procedure(final), deferred     :: final
     procedure(write_i_0), deferred :: write_i_0
     procedure(write_i_1), deferred :: write_i_1
     procedure(write_r_0), deferred :: write_r_0
     procedure(write_r_1), deferred :: write_r_1
     procedure(write_c_0), deferred :: write_c_0
     procedure(write_c_1), deferred :: write_c_1
     procedure(write_a_0), deferred :: write_a_0
     procedure(write_a_1), deferred :: write_a_1
     generic, public                :: write => write_i_0, write_i_1, write_r_0, write_r_1, &
                                                write_c_0, write_c_1, write_a_0, write_a_1
  end type writer_t

  ! Interfaces

  abstract interface

     subroutine final (this)
       import writer_t
       class(writer_t), intent(inout) :: this
     end subroutine final

     subroutine write_i_0 (this, name, data)
       use core_kinds
       import writer_t
       class(writer_t), intent(inout) :: this
       character(*), intent(in)       :: name
       integer, intent(in)         :: data
     end subroutine write_i_0

     subroutine write_i_1 (this, name, data)
       use core_kinds
       import writer_t
       class(writer_t), intent(inout) :: this
       character(*), intent(in)       :: name
       integer, intent(in)         :: data(:)
     end subroutine write_i_1

     subroutine write_r_0 (this, name, data)
       use core_kinds
       import writer_t
       class(writer_t), intent(inout) :: this
       character(*), intent(in)       :: name
       real(WP), intent(in)         :: data
     end subroutine write_r_0

     subroutine write_r_1 (this, name, data)
       use core_kinds
       import writer_t
       class(writer_t), intent(inout) :: this
       character(*), intent(in)       :: name
       real(WP), intent(in)         :: data(:)
     end subroutine write_r_1

     subroutine write_c_0 (this, name, data)
       use core_kinds
       import writer_t
       class(writer_t), intent(inout) :: this
       character(*), intent(in)       :: name
       complex(WP), intent(in)         :: data
     end subroutine write_c_0

     subroutine write_c_1 (this, name, data)
       use core_kinds
       import writer_t
       class(writer_t), intent(inout) :: this
       character(*), intent(in)       :: name
       complex(WP), intent(in)         :: data(:)
     end subroutine write_c_1

     subroutine write_a_0 (this, name, data)
       use core_kinds
       import writer_t
       class(writer_t), intent(inout) :: this
       character(*), intent(in)       :: name
       character(*), intent(in)         :: data
     end subroutine write_a_0

     subroutine write_a_1 (this, name, data)
       use core_kinds
       import writer_t
       class(writer_t), intent(inout) :: this
       character(*), intent(in)       :: name
       character(*), intent(in)         :: data(:)
     end subroutine write_a_1

  end interface

  ! Access specifiers

  private

  public :: writer_t

end module gyre_writer

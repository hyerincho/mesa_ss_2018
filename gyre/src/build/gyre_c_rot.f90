!fpx3_header(0.13_3a)
!
!dependencies
!   dir: ~/Documents/MESA_SS/gyre/src/build 
!   sources: -
!   includes: ../rot/gyre_rot.inc ../extern/core/core.inc
!   uses: core_kinds
!   provides: gyre_c_rot
!end dependencies
!
!end fpx3_header
! Incfile  : gyre_r_rot
! Purpose  : rotational effects (complex)
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

! Incfile  : gyre_rot
! Purpose  : rotational effects (template)
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

module gyre_c_rot

  ! Uses

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type, abstract :: c_rot_t
   contains
     procedure(l_e), deferred    :: l_e
     procedure(lambda), deferred :: lambda
  end type c_rot_t

  ! Interfaces

  abstract interface

     function l_e (this, Omega_rot, omega)
       use core_kinds
       import c_rot_t
       class(c_rot_t), intent(in) :: this
       real(WP), intent(in)          :: Omega_rot
       complex(WP), intent(in)         :: omega
       complex(WP)                     :: l_e
     end function l_e

     function lambda (this, Omega_rot, omega)
       use core_kinds
       import c_rot_t
       class(c_rot_t), intent(in) :: this
       real(WP), intent(in)          :: Omega_rot
       complex(WP), intent(in)         :: omega
       complex(WP)                     :: lambda
     end function lambda

  end interface

  ! Access specifiers

  private

  public :: c_rot_t

end module gyre_c_rot


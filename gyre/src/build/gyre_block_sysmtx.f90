!fpx3_header(0.13_3a)
!
!dependencies
!   dir: ~/Documents/MESA_SS/gyre/src/build 
!   sources: -
!   includes: 
!   uses: gyre_c_block_sysmtx gyre_r_block_sysmtx
!   provides: gyre_block_sysmtx
!end dependencies
!
!end fpx3_header
! Module   : gyre_block_sysmtx
! Purpose  : system matrix (block storage)
!
! Copyright 2013-215 Rich Townsend
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

module gyre_block_sysmtx

  ! Uses

  use gyre_r_block_sysmtx
  use gyre_c_block_sysmtx

end module gyre_block_sysmtx

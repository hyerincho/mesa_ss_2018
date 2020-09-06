!fpx3_header(0.13_3a)
!
!dependencies
!   dir: ~/Documents/MESA_SS/gyre/src/build 
!   sources: -
!   includes: 
!   uses: 
!   provides: core_kinds
!end dependencies
!
!end fpx3_header
! Module   : core_kinds
! Purpose  : kind types

module core_kinds

  ! No implicit typing

  implicit none

  ! Parameter definitions

  integer, parameter :: I4 = SELECTED_INT_KIND(9)
  integer, parameter :: I8 = SELECTED_INT_KIND(14)

  integer, parameter :: SP = KIND(0.)
  integer, parameter :: DP = KIND(0.D0)
  integer, parameter :: QP = SELECTED_REAL_KIND(32)

  integer, parameter :: WP = DP

  ! Acess specifiers

  private

  public :: I4
  public :: I8
  public :: SP
  public :: DP
  public :: WP
  public :: QP

end module core_kinds

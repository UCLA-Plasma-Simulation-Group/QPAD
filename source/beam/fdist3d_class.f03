module fdist3d_class

use parallel_module
use options_class
use input_class
use iso_c_binding
use param

implicit none

private

public :: fdist3d
public :: p_pf3d_std

integer, parameter :: p_pf3d_std = 1, p_pf3d_rnd = 2, p_pf3d_file = 3

type, abstract :: fdist3d

  ! Maximum number of particles in this partition
  integer(kind=LG) :: npmax

  ! If update beam momentum
  logical :: evol = .true.

  ! If enable spin dynamics
  logical :: has_spin = .false.

  ! If enable quiet start initialization
  Logical :: quiet = .true.

  ! Charge & charge/mass ratio
  real :: qm, qbm

  ! Anomalous magnet moment
  real :: amm

  ! Cell sizes
  real :: dr, dz

  ! Grid offset and number of cells
  integer :: noff_r, noff_z, nr, nz, nrp, nzp

  ! Position of the lower boundary in the longitudinal direction
  real :: z0

  contains

  procedure(init_fdist3d_intf), deferred :: new
  procedure(end_fdist3d_intf), deferred :: del
  procedure(inject_fdist3d_intf), deferred :: inject

end type fdist3d

abstract interface

subroutine init_fdist3d_intf( this, input, opts, sect_id )
  import fdist3d, input_json, options
  implicit none
  class( fdist3d ), intent(inout) :: this
  type( input_json ), intent(inout) :: input
  type( options ), intent(in) :: opts
  integer, intent(in) :: sect_id
end subroutine init_fdist3d_intf

subroutine end_fdist3d_intf( this )
  import fdist3d
  implicit none
  class( fdist3d ), intent(inout) :: this
end subroutine end_fdist3d_intf

subroutine inject_fdist3d_intf( this, x, p, s, q, npp )
  import fdist3d, LG
  implicit none
  class( fdist3d ), intent(inout) :: this
  real, intent(inout), dimension(:,:) :: x, p, s
  real, intent(inout), dimension(:) :: q
  integer(kind=LG), intent(inout) :: npp
end subroutine inject_fdist3d_intf

end interface

end module fdist3d_class
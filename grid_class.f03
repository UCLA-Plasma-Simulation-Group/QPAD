module grid_class

use parallel_class
use parallel_pipe_class
use input_class
use system

implicit none

private

character(len=18), save :: cls_name = 'grid'
integer, parameter :: cls_level = 0

public :: grid

type :: grid

private

integer, dimension(2) :: nd ! number of global grid points
integer, dimension(2) :: ndp ! number of local grid points
integer, dimension(2) :: nvp ! number of processors
integer, dimension(2) :: noff ! grid index offset

contains

generic :: new => init_grid
generic :: del => end_grid
generic :: noff => get_noff, get_noff_dim
generic :: nd => get_nd, get_nd_dim
generic :: ndp => get_ndp, get_ndp_dim

end type grid

contains

subroutine init_grid( this, pp, input )

  implicit none

  class( grid ), intent(inout) :: this
  type( input_json ), pointer, intent(inout) :: input

  integer :: lidproc, stageid, local_size, extra
  character(len=18), save :: sname = 'init_grid'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  call input%get( 'simulation.nr', this%nd(1) )
  call input%get( 'simulation.nz', this%nd(2) )

  this%nvp(1) = pp%getlnvp()
  this%nvp(2) = pp%getnstage()

  lidproc = pp%getlidproc()
  stageid = pp%getstageid()

  ! the un-evenly distributed grid points among processors are accounted for
  local_size = this%nd(1) / this%nvp(1)
  extra = this%nd(1) - local_size * this%nvp(1)
  this%noff(1) = local_size * lidproc + min( lidproc, extra )
  this%ndp(1) = local_size + merge( 1, 0, lidproc<extra )

  local_size = this%nd(2) / this%nvp(2)
  extra = this%nd(2) - local_size * this%nvp(2)
  this%noff(2) = local_size * stageid + min( stageid, extra )
  this%ndp(2) = local_size + merge( 1, 0, stageid<extra )

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine init_grid

subroutine end_grid( this )

  implicit none

  class( grid ), intent(inout) :: this
  
  character(len=18), save :: sname = 'init_grid'

  call write_dbg( cls_name, sname, cls_level, 'starts' )
  ! this routine is only a place-holder now, do nothing.
  call write_dbg( cls_name, sname, cls_level, 'ends' )
  
end subroutine end_grid

function get_nd( this )

  implicit none

  class( grid ), intent(in) :: this
  integer, dimension(2) :: get_nd

  get_nd = this%nd

end function get_nd

function get_nd_dim( this, dim )

  implicit none

  class( grid ), intent(in) :: this
  integer, intent(in) :: dim
  integer :: get_nd_dim

  get_nd_dim = this%nd(dim)
  
end function get_nd_dim

function get_ndp( this )

  implicit none

  class( grid ), intent(in) :: this
  integer, dimension(2) :: get_ndp

  get_ndp = this%ndp

end function get_ndp

function get_ndp_dim( this, dim )

  implicit none

  class( grid ), intent(in) :: this
  integer, intent(in) :: dim
  integer :: get_ndp_dim

  get_ndp_dim = this%ndp(dim)
  
end function get_ndp_dim

function get_noff( this )

  implicit none

  class( grid ), intent(in) :: this
  integer, dimension(2) :: get_noff

  get_noff = this%noff

end function get_noff

function get_noff_dim( this, dim )

  implicit none

  class( grid ), intent(in) :: this
  integer, intent(in) :: dim
  integer :: get_noff_dim

  get_noff_dim = this%noff(dim)
  
end function get_noff_dim

end module grid_class
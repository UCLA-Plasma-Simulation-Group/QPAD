module options_class

use parallel_module
use sysutil_module
use param
use input_class

implicit none

private

character(len=32), save :: cls_name = 'options'
integer, parameter :: cls_level = 0

public :: options

type :: options

  private

  integer, dimension(2) :: nd ! number of global grid points
  integer, dimension(2) :: ndp ! number of local grid points
  integer, dimension(2) :: noff ! grid index offset
  real :: dr, dxi
  integer, public :: algorithm

contains

  generic :: get_noff => get_noff_all, get_noff_dim
  generic :: get_nd   => get_nd_all, get_nd_dim
  generic :: get_ndp  => get_ndp_all, get_ndp_dim

  procedure :: new => init_options
  procedure :: del => end_options
  procedure :: get_dr, get_dxi
  procedure, private :: get_noff_all, get_noff_dim
  procedure, private :: get_nd_all, get_nd_dim
  procedure, private :: get_ndp_all, get_ndp_dim

end type options

contains

subroutine init_options( this, input_file )

  implicit none

  class( options ), intent(inout) :: this
  type( input_json ), intent(inout) :: input_file

  real :: min_val, max_val
  integer :: nr, nz, local_size, extra
  character(len=:), allocatable :: read_str

  ! read simulation algorithm
  call input_file%get( 'simulation.algorithm', read_str )

  select case ( trim(read_str) )
  case ("standard")
    this%algorithm = p_sim_standard

  case ("popas")
    this%algorithm = p_sim_popas

  case ("template")
    this%algorithm = p_sim_tmplt

  ! add new entries here if there are other algorithm modules
  ! case ("some_other_algorithm")
  !   this%algorithm = p_sim_some_other

  case default
    this%algorithm = p_sim_standard
  end select

  deallocate(read_str)

  ! read grid settings
  call input_file%get( 'simulation.grid(1)', nr )
  call input_file%get( 'simulation.grid(2)', nz )
  this%nd(1) = nr
  this%nd(2) = nz

  call input_file%get( 'simulation.box.r(1)', min_val )
  call input_file%get( 'simulation.box.r(2)', max_val )
  this%dr = ( max_val - min_val ) / nr

  call input_file%get( 'simulation.box.z(1)', min_val )
  call input_file%get( 'simulation.box.z(2)', max_val )
  this%dxi = ( max_val - min_val ) / nz

  ! the un-evenly distributed grid points among processors are accounted for
  local_size   = nr / num_procs_loc()
  extra        = nr - local_size * num_procs_loc()
  this%noff(1) = local_size * id_proc_loc() + min( id_proc_loc(), extra )
  this%ndp(1)  = local_size + merge( 1, 0, id_proc_loc() < extra )

  local_size   = nz / num_stages()
  extra        = nz - local_size * num_stages()
  this%noff(2) = local_size * id_stage() + min( id_stage(), extra )
  this%ndp(2)  = local_size + merge( 1, 0, id_stage() < extra )

  ! read particle interpolation order
  ! call in%get('simulation.interpolation', read_str)
  ! select case (trim(read_str))
  ! case ("linear")
  !   this%interpolation = p_interp_linear
  ! case default
  !   this%interpolation = p_interp_linear
  ! end select

end subroutine init_options

subroutine end_options( this )

  implicit none

  class( options ), intent(inout) :: this

  character(len=32), save :: sname = 'end_options'

  call write_dbg( cls_name, sname, cls_level, 'starts' )
  ! this routine is only a place-holder now, do nothing.
  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine end_options

function get_nd_all( this )

  implicit none

  class( options ), intent(in) :: this
  integer, dimension(2) :: get_nd_all

  get_nd_all = this%nd

end function get_nd_all

function get_nd_dim( this, dim )

  implicit none

  class( options ), intent(in) :: this
  integer, intent(in) :: dim
  integer :: get_nd_dim

  get_nd_dim = this%nd(dim)

end function get_nd_dim

function get_ndp_all( this )

  implicit none

  class( options ), intent(in) :: this
  integer, dimension(2) :: get_ndp_all

  get_ndp_all = this%ndp

end function get_ndp_all

function get_ndp_dim( this, dim )

  implicit none

  class( options ), intent(in) :: this
  integer, intent(in) :: dim
  integer :: get_ndp_dim

  get_ndp_dim = this%ndp(dim)

end function get_ndp_dim

function get_noff_all( this )

  implicit none

  class( options ), intent(in) :: this
  integer, dimension(2) :: get_noff_all

  get_noff_all = this%noff

end function get_noff_all

function get_noff_dim( this, dim )

  implicit none

  class( options ), intent(in) :: this
  integer, intent(in) :: dim
  integer :: get_noff_dim

  get_noff_dim = this%noff(dim)

end function get_noff_dim

function get_dr( this )

  implicit none

  class( options ), intent(in) :: this
  real :: get_dr

  get_dr = this%dr

end function get_dr

function get_dxi( this )

  implicit none

  class( options ), intent(in) :: this
  real :: get_dxi

  get_dxi = this%dxi

end function get_dxi

end module options_class
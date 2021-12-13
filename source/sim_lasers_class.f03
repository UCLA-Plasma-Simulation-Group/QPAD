module sim_lasers_class

use parallel_module
use options_class
use sysutil_module
use param
use input_class
use field_class
use field_laser_class
use species2d_class
use kwargs_class
use ppmsg_class

implicit none

private

public :: sim_lasers

type, public :: sim_lasers
  
  class( field_laser ), dimension(:), pointer :: laser => null()
  class( field_laser ), pointer :: laser_all => null()
  class( field ), pointer :: chi => null()
  integer :: num_lasers
  type( ppmsg ), dimension(:), allocatable :: pp_msg

  contains

  procedure :: alloc       => alloc_sim_lasers
  procedure :: new         => init_sim_lasers
  procedure :: del         => end_sim_lasers
  procedure :: deposit_chi => deposit_chi_sim_lasers
  procedure :: advance     => advance_sim_lasers

end type sim_lasers

character(len=32), save :: cls_name = 'sim_lasers'
integer, save :: cls_level = 2

contains

subroutine alloc_sim_lasers( this, input, opts )

  implicit none
  class( sim_lasers ), intent(inout) :: this
  type( input_json ), intent(inout) :: input
  type( options ), intent(in) :: opts
  ! local data
  integer :: max_mode, i
  character(len=32), save :: sname = 'alloc_sim_lasers'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  call input%get( 'simulation.nlasers', this%num_lasers )
  call input%get( 'simulation.max_mode', max_mode )

  if ( .not. associated( this%laser ) ) then
    allocate( field_laser :: this%laser( this%num_lasers ) )
  endif

  if ( .not. associated( this%laser_all ) ) then
    allocate( field_laser :: this%laser_all )
  endif

  if ( .not. associated( this%chi ) ) then
    allocate( field :: this%chi )
  endif

  do i = 1, this%num_lasers
    call this%laser(i)%alloc( input, opts, i )
  enddo

  if ( .not. allocated( this%pp_msg ) ) then
    allocate( this%pp_msg( this%num_lasers ) )
  endif

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine alloc_sim_lasers

subroutine init_sim_lasers( this, input, opts )

  implicit none

  class( sim_lasers ), intent(inout) :: this
  type( input_json ), intent(inout) :: input
  type( options ), intent(in) :: opts

  ! local data
  character(len=32), save :: sname = 'init_sim_lasers'
  character(len=:), allocatable :: str
  integer :: max_mode, iter, i
  integer, dimension(2,2) :: gc_num
  real :: k0 
  type( kw_list ) :: kwargs
  logical :: read_rst

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  call input%get( 'simulation.max_mode', max_mode )
  call input%get( 'simulation.read_restart', read_rst )
  call kwargs%append( 'iter', 0 )
  call kwargs%append( 'k0', 10.0 )

  gc_num(:,1) = (/1, 1/)
  gc_num(:,2) = (/2, 1/)

  do i = 1, this%num_lasers
    call input%get( 'laser(' // num2str(i) // ').iteration', iter )
    call input%get( 'laser(' // num2str(i) // ').k0', k0 )
    call kwargs%set( 'iter', iter )
    call kwargs%set( 'k0', k0 )
    call this%laser(i)%new( opts, 1, max_mode, gc_num, only_f1=.false., kwargs=kwargs )

    if ( read_rst ) then
      call write_err( 'Restarting is currently not available for laser field.' )
    endif

  enddo

  call this%laser_all%new_aux( opts, 1, max_mode, gc_num, only_f1=.true. )
  call this%chi%new( opts, 1, max_mode, gc_num )

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine init_sim_lasers

subroutine end_sim_lasers( this )

  implicit none
  class(sim_lasers), intent(inout) :: this
  ! local data
  character(len=32), save :: sname = 'end_sim_lasers'
  integer :: i

  call write_dbg( cls_name, sname, cls_level, 'starts' )
  do i = 1, this%num_lasers
    call this%laser(i)%del()
  enddo
  call this%laser_all%del_aux()
  call this%chi%del()
  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine end_sim_lasers

subroutine deposit_chi_sim_lasers( this, species, slice_idx )

  implicit none
  class( sim_lasers ), intent(inout) :: this
  class( species2d ), intent(inout), dimension(:), pointer :: species
  integer, intent(in) :: slice_idx

  character(len=32), save :: sname = 'deposit_chi_sim_lasers'
  integer :: k, num_species

  call write_dbg( cls_name, sname, cls_level, 'starts' )
  
  num_species = size( species )
  this%chi = 0.0
  do k = 1, num_species
    call species(k)%deposit_chi( this%chi )
  enddo
  call this%chi%acopy_gc_f1( dir=p_mpi_forward )
  call this%chi%smooth_f1()
  call this%chi%copy_gc_f1()
  call this%chi%copy_slice( slice_idx, dir=p_copy_1to2 )

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine deposit_chi_sim_lasers

subroutine advance_sim_lasers( this )

  implicit none
  class( sim_lasers ), intent(inout) :: this

  integer :: k, gc
  character(len=32), save :: sname = 'advance_sim_lasers'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  do k = 1, this%num_lasers
    gc = this%laser(k)%gc_num(1,2)
    call this%laser(k)%set_rhs( this%chi )
    call this%pp_msg(k)%get_tag()
    call this%laser(k)%pipe_recv( this%pp_msg(k), 'forward', 'guard', 'replace', gc )
    call this%laser(k)%solve( this%chi )
    call this%laser(k)%pipe_send( this%pp_msg(k), 'forward', 'inner', gc )
  enddo

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine advance_sim_lasers

end module sim_lasers_class
module field_class

use parallel_pipe_class
use grid_class
use ufield_class
use ufield_smooth_class
use hdf5io_class
use param
use sys

implicit none

private

interface add_f1
  module procedure add_f1_binary
  module procedure add_f1_binary_dim
  module procedure add_f1_unitary
  module procedure add_f1_unitary_dim
end interface

interface add_f2
  module procedure add_f2_binary
  module procedure add_f2_unitary
end interface

interface sub_f1
  module procedure sub_f1_binary
  module procedure sub_f1_binary_dim
  module procedure sub_f1_unitary
  module procedure sub_f1_unitary_dim
end interface

interface sub_f2
  module procedure sub_f2_binary
  module procedure sub_f2_unitary
end interface

interface dot_f1
  module procedure dot_f1_unitary
  module procedure dot_f1_unitary_dim
end interface

public :: field
public :: add_f1, add_f2, sub_f1, sub_f2, dot_f1

character(len=20), parameter :: cls_name = "field"
integer, parameter :: cls_level = 3

type :: field

  ! private

  class( ufield ), dimension(:), pointer :: rf_re => null()
  class( ufield ), dimension(:), pointer :: rf_im => null()
  type( ufield_smooth ) :: smooth

  real :: dr, dxi
  integer :: num_modes
  integer :: entity

  contains

  generic :: new => init_field, init_field_cp
  procedure :: del => end_field
  generic :: get_rf_re => get_rf_re_all, get_rf_re_mode
  generic :: get_rf_im => get_rf_im_all, get_rf_im_mode
  generic :: write_hdf5 => write_hdf5_single, write_hdf5_pipe
  procedure :: smooth_f1
  procedure :: copy_slice
  procedure :: get_dr, get_dxi, get_num_modes
  procedure :: copy_gc_f1, copy_gc_f2, copy_gc_stage
  procedure :: acopy_gc_f1, acopy_gc_f2, acopy_gc_stage

  procedure, private :: init_field, init_field_cp, end_field
  procedure, private :: get_rf_re_all, get_rf_re_mode, get_rf_im_all, get_rf_im_mode
  procedure, private :: write_hdf5_single, write_hdf5_pipe

  generic :: assignment(=)   => assign_f1
  generic :: as              => assign_f2
  ! generic :: operator(+)     => add_f1_v1, add_f1_v2
  ! generic :: operator(-)     => sub_f1_v1, sub_f1_v2
  ! generic :: operator(*)     => dot_f1_v1, dot_f1_v2
  ! generic :: operator(.add.) => add_f2_v1, add_f2_v2
  ! generic :: operator(.sub.) => sub_f2_v1, sub_f2_v2
  ! generic :: operator(.dot.) => dot_f2_v1, dot_f2_v2

  ! procedure, private, pass(a1) :: add_f1_v1, add_f1_v2
  ! procedure, private, pass(a1) :: dot_f1_v1, dot_f1_v2
  ! procedure, private, pass(a1) :: sub_f1_v1, sub_f1_v2
  ! procedure, private, pass(a1) :: add_f2_v1, add_f2_v2
  ! procedure, private, pass(a1) :: dot_f2_v1, dot_f2_v2
  ! procedure, private, pass(a1) :: sub_f2_v1, sub_f2_v2
  procedure, private :: assign_f1, assign_f2

end type field


contains

! =====================================================================
! Class field implementation
! =====================================================================
subroutine init_field( this, pp, gp, dim, num_modes, gc_num, &
  entity, smooth_type, smooth_order )

  implicit none

  class( field ), intent(inout) :: this
  class( parallel_pipe ), intent(in), pointer :: pp
  class( grid ), intent(in), pointer :: gp
  integer, intent(in) :: num_modes, dim
  integer, intent(in), dimension(2,2) :: gc_num
  integer, intent(in), optional :: entity, smooth_type, smooth_order

  integer :: i
  integer, dimension(2,2) :: gc_num_new
  character(len=20), save :: sname = "init_field"

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  this%num_modes = num_modes
  this%dr  = gp%get_dr()
  this%dxi = gp%get_dxi()

  if ( present(entity) ) then
    this%entity = entity
  else
    this%entity = p_entity_none
  endif

  if ( present(smooth_type) .and. present(smooth_order) ) then
    call this%smooth%new( smooth_type, smooth_order )
    gc_num_new(1,1) = max( gc_num(1,1), smooth_order )
    gc_num_new(2,1) = max( gc_num(2,1), smooth_order )
    gc_num_new(:,2) = gc_num(:,2)
  else
    call this%smooth%new( p_smooth_none, 0 )
  endif

  allocate( this%rf_re(0:num_modes) )
  allocate( this%rf_im(num_modes) )
  do i = 0, this%num_modes
    call this%rf_re(i)%new( pp, gp, dim, i, gc_num, has_2d=.true. )
    if (i==0) cycle
    call this%rf_im(i)%new( pp, gp, dim, i, gc_num, has_2d=.true. )
  enddo

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine init_field

subroutine init_field_cp( this, that )

  implicit none

  class( field ), intent(inout) :: this
  class( field ), intent(in) :: that

  integer :: i

  this%num_modes = that%get_num_modes()
  this%dr        = that%get_dr()
  this%dxi       = that%get_dxi()

  allocate( this%rf_re(0:this%num_modes) )
  allocate( this%rf_im(this%num_modes) )
  do i = 0, this%num_modes
    call this%rf_re(i)%new( that%get_rf_re(i) )
    if (i==0) cycle
    call this%rf_im(i)%new( that%get_rf_im(i) )
  enddo

end subroutine init_field_cp

subroutine end_field( this )

  implicit none

  class( field ), intent(inout) :: this

  integer :: i
  character(len=20), save :: sname = "end_field"

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  call this%rf_re(0)%del()
  do i = 1, this%num_modes
    call this%rf_re(i)%del()
    call this%rf_im(i)%del()
  enddo
  deallocate( this%rf_re, this%rf_im )

  call this%smooth%del()

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine end_field

subroutine copy_slice( this, idx, dir )

  implicit none

  class( field ), intent(inout) :: this
  integer, intent(in) :: idx, dir

  integer :: i
  character(len=20), save :: sname = "copy_slice"

  call write_dbg( cls_name, sname, cls_level, 'starts' )
  call start_tprof( 'copy slices' )

  do i = 0, this%num_modes
    call this%rf_re(i)%copy_slice( idx, dir )
    if ( i == 0 ) cycle
    call this%rf_im(i)%copy_slice( idx, dir )
  enddo

  call stop_tprof( 'copy slices' )
  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine copy_slice

subroutine copy_gc_f1( this, bnd_ax )

  implicit none

  class( field ), intent(inout) :: this
  logical, intent(in) :: bnd_ax

  integer :: i
  character(len=20), save :: sname = "copy_gc_f1"

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  do i = 0, this%num_modes
    call this%rf_re(i)%copy_gc_f1( bnd_ax )
    if ( i == 0 ) cycle
    call this%rf_im(i)%copy_gc_f1( bnd_ax )
  enddo

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine copy_gc_f1

subroutine copy_gc_f2( this, bnd_ax )

  implicit none

  class( field ), intent(inout) :: this
  logical, intent(in) :: bnd_ax

  integer :: i
  character(len=20), save :: sname = "copy_gc_f2"

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  do i = 0, this%num_modes
    call this%rf_re(i)%copy_gc_f2( bnd_ax )
    if ( i == 0 ) cycle
    call this%rf_im(i)%copy_gc_f2( bnd_ax )
  enddo

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine copy_gc_f2

subroutine copy_gc_stage( this, dir )

  implicit none

  class( field ), intent(inout) :: this
  integer, intent(in) :: dir

  integer :: i
  character(len=20), save :: sname = "copy_gc_stage"

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  do i = 0, this%num_modes
    call this%rf_re(i)%copy_gc_stage(dir)
    if ( i == 0 ) cycle
    call this%rf_im(i)%copy_gc_stage(dir)
  enddo

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine copy_gc_stage

subroutine acopy_gc_f1( this )

  implicit none

  class( field ), intent(inout) :: this

  integer :: i
  character(len=20), save :: sname = "acopy_gc_f1"

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  do i = 0, this%num_modes
    call this%rf_re(i)%acopy_gc_f1()
    if ( i == 0 ) cycle
    call this%rf_im(i)%acopy_gc_f1()
  enddo

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine acopy_gc_f1

subroutine acopy_gc_f2( this )

  implicit none

  class( field ), intent(inout) :: this

  integer :: i
  character(len=20), save :: sname = "acopy_gc_f2"

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  do i = 0, this%num_modes
    call this%rf_re(i)%acopy_gc_f2()
    if ( i == 0 ) cycle
    call this%rf_im(i)%acopy_gc_f2()
  enddo

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine acopy_gc_f2

subroutine acopy_gc_stage( this )

  implicit none

  class( field ), intent(inout) :: this

  integer :: i
  character(len=20), save :: sname = "acopy_gc_stage"

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  do i = 0, this%num_modes
    call this%rf_re(i)%acopy_gc_stage()
    if ( i == 0 ) cycle
    call this%rf_im(i)%acopy_gc_stage()
  enddo

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine acopy_gc_stage

subroutine write_hdf5_single( this, files, dim )

  implicit none

  class( field ), intent(inout) :: this
  class( hdf5file ), intent(in), dimension(:) :: files
  integer, intent(in) :: dim

  integer :: i
  character(len=32), save :: sname = 'write_hdf5_single'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  do i = 0, this%num_modes

    if ( i == 0 ) then
      call this%rf_re(i)%write_hdf5( files(1), dim )
      cycle
    endif

    call this%rf_re(i)%write_hdf5( files(2*i), dim )
    call this%rf_im(i)%write_hdf5( files(2*i+1), dim )

  enddo

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine write_hdf5_single

subroutine write_hdf5_pipe( this, files, dim, rtag, stag, id )

  implicit none

  class( field ), intent(inout) :: this
  class( hdf5file ), intent(in), dimension(:) :: files
  integer, intent(in) :: dim, rtag, stag
  integer, intent(inout) :: id

  integer :: i
  character(len=32), save :: sname = 'write_hdf5_pipe'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  do i = 0, this%num_modes

    if ( i == 0 ) then
      call this%rf_re(i)%write_hdf5( files(1), dim, rtag, stag, id )
      cycle
    endif

    call this%rf_re(i)%write_hdf5( files(2*i), dim, rtag, stag, id )
    call this%rf_im(i)%write_hdf5( files(2*i+1), dim, rtag, stag, id )

  enddo

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine write_hdf5_pipe

subroutine smooth_f1( this )

  implicit none

  class( field ), intent(inout) :: this

  integer :: i
  character(len=32), save :: sname = 'smooth_f1'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  if ( .not. this%smooth%if_smooth() ) return
  
  do i = 0, this%num_modes

    if ( i == 0 ) then
      call this%smooth%smooth_f1( this%rf_re(i), this%dr )
      cycle
    endif

    call this%smooth%smooth_f1( this%rf_re(i), this%dr )
    call this%smooth%smooth_f1( this%rf_im(i), this%dr )

  enddo

  call this%copy_gc_f1( bnd_ax = .false. )

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine smooth_f1

function get_dr( this )

  implicit none

  class( field ), intent(in) :: this
  real :: get_dr

  get_dr = this%dr

end function get_dr

function get_dxi( this )

  implicit none

  class( field ), intent(in) :: this
  real :: get_dxi

  get_dxi = this%dxi
  
end function get_dxi

function get_num_modes( this )

  implicit none

  class( field ), intent(in) :: this
  integer :: get_num_modes

  get_num_modes = this%num_modes
  
end function get_num_modes

function get_rf_re_all( this )

  implicit none

  class( field ), intent(in) :: this
  type( ufield ), dimension(:), pointer :: get_rf_re_all

  get_rf_re_all => this%rf_re

end function get_rf_re_all

function get_rf_re_mode( this, mode )

  implicit none

  class( field ), intent(in) :: this
  integer, intent(in) :: mode
  type( ufield ), pointer :: get_rf_re_mode

  get_rf_re_mode => this%rf_re(mode)

end function get_rf_re_mode

function get_rf_im_all( this )

  implicit none

  class( field ), intent(in) :: this
  type( ufield ), dimension(:), pointer :: get_rf_im_all

  get_rf_im_all => this%rf_im

end function get_rf_im_all

function get_rf_im_mode( this, mode )

  implicit none

  class( field ), intent(in) :: this
  integer, intent(in) :: mode
  type( ufield ), pointer :: get_rf_im_mode

  get_rf_im_mode => this%rf_im(mode)

end function get_rf_im_mode

subroutine assign_f1( this, that )

  implicit none

  class( field ), intent(inout) :: this
  class(*), intent(in) :: that

  integer :: i

  select type (that)

  type is (real)

    do i = 0, this%num_modes
      this%rf_re(i) = that
      if (i == 0) cycle
      this%rf_im(i) = that
    enddo

  class is (field)

    do i = 0, this%num_modes
      this%rf_re(i) = that%rf_re(i)
      if (i == 0) cycle
      this%rf_im(i) = that%rf_im(i)
    enddo

  class default

    call write_err( "invalid assignment type!" )

  end select

end subroutine assign_f1

subroutine assign_f2( this, that )

  implicit none

  class( field ), intent(inout) :: this
  class(*), intent(in) :: that

  integer :: i

  select type (that)

  type is (real)

    do i = 0, this%num_modes
      call this%rf_re(i)%as( that )
      if (i == 0) cycle
      call this%rf_im(i)%as( that )
    enddo

  class is (field)

    do i = 0, this%num_modes
      call this%rf_re(i)%as( that%get_rf_re(i) )
      if (i == 0) cycle
      call this%rf_im(i)%as( that%get_rf_im(i) )
    enddo

  class default

    call write_err( "invalid assignment type!" )

  end select

end subroutine assign_f2

subroutine add_f1_binary_dim( a1, a2, a3, dim1, dim2, dim3 )

  implicit none

  class( field ), intent(in) :: a1, a2
  class( field ), intent(inout) :: a3
  integer, intent(in), dimension(:) :: dim1, dim2, dim3

  class( ufield ), dimension(:), pointer :: ua3_re => null(), ua3_im => null()
  integer :: i

  ua3_re => a3%get_rf_re()
  ua3_im => a3%get_rf_im()

  do i = 0, a1%num_modes
    call add_f1( a1%rf_re(i), a2%rf_re(i), ua3_re(i), dim1, dim2, dim3 )
    if (i==0) cycle
    call add_f1( a1%rf_im(i), a2%rf_im(i), ua3_im(i), dim1, dim2, dim3 )
  enddo

end subroutine add_f1_binary_dim

subroutine add_f1_unitary_dim( a1, a2, dim1, dim2 )

  implicit none

  class( field ), intent(in) :: a1
  class( field ), intent(inout) :: a2
  integer, intent(in), dimension(:) :: dim1, dim2

  class( ufield ), dimension(:), pointer :: ua2_re => null(), ua2_im => null()
  integer :: i

  ua2_re => a2%get_rf_re()
  ua2_im => a2%get_rf_im()

  do i = 0, a1%num_modes
    call add_f1( a1%rf_re(i), ua2_re(i), dim1, dim2 )
    if (i==0) cycle
    call add_f1( a1%rf_im(i), ua2_im(i), dim1, dim2 )
  enddo

end subroutine add_f1_unitary_dim

subroutine add_f1_binary( a1, a2, a3 )

  implicit none

  class( field ), intent(in) :: a1, a2
  class( field ), intent(inout) :: a3

  class( ufield ), dimension(:), pointer :: ua3_re => null(), ua3_im => null()
  integer :: i

  ua3_re => a3%get_rf_re()
  ua3_im => a3%get_rf_im()

  do i = 0, a1%num_modes
    call add_f1( a1%rf_re(i), a2%rf_re(i), ua3_re(i) )
    if (i==0) cycle
    call add_f1( a1%rf_im(i), a2%rf_im(i), ua3_im(i) )
  enddo

end subroutine add_f1_binary

subroutine add_f1_unitary( a1, a2 )

  implicit none

  class( field ), intent(in) :: a1
  class( field ), intent(inout) :: a2

  class( ufield ), dimension(:), pointer :: ua2_re => null(), ua2_im => null()
  integer :: i

  ua2_re => a2%get_rf_re()
  ua2_im => a2%get_rf_im()

  do i = 0, a1%num_modes
    call add_f1( a1%rf_re(i), ua2_re(i) )
    if (i==0) cycle
    call add_f1( a1%rf_im(i), ua2_im(i) )
  enddo

end subroutine add_f1_unitary

subroutine sub_f1_binary_dim( a1, a2, a3, dim1, dim2, dim3 )

  implicit none

  class( field ), intent(in) :: a1, a2
  class( field ), intent(inout) :: a3
  integer, intent(in), dimension(:) :: dim1, dim2, dim3

  class( ufield ), dimension(:), pointer :: ua3_re => null(), ua3_im => null()
  integer :: i

  ua3_re => a3%get_rf_re()
  ua3_im => a3%get_rf_im()

  do i = 0, a1%num_modes
    call sub_f1( a1%rf_re(i), a2%rf_re(i), ua3_re(i), dim1, dim2, dim3 )
    if (i==0) cycle
    call sub_f1( a1%rf_im(i), a2%rf_im(i), ua3_im(i), dim1, dim2, dim3 )
  enddo

end subroutine sub_f1_binary_dim

subroutine sub_f1_unitary_dim( a1, a2, dim1, dim2 )

  implicit none

  class( field ), intent(in) :: a1
  class( field ), intent(inout) :: a2
  integer, intent(in), dimension(:) :: dim1, dim2

  class( ufield ), dimension(:), pointer :: ua2_re => null(), ua2_im => null()
  integer :: i

  ua2_re => a2%get_rf_re()
  ua2_im => a2%get_rf_im()

  do i = 0, a1%num_modes
    call sub_f1( a1%rf_re(i), ua2_re(i), dim1, dim2 )
    if (i==0) cycle
    call sub_f1( a1%rf_im(i), ua2_im(i), dim1, dim2 )
  enddo

end subroutine sub_f1_unitary_dim

subroutine sub_f1_binary( a1, a2, a3 )

  implicit none

  class( field ), intent(in) :: a1, a2
  class( field ), intent(inout) :: a3

  class( ufield ), dimension(:), pointer :: ua3_re => null(), ua3_im => null()
  integer :: i

  ua3_re => a3%get_rf_re()
  ua3_im => a3%get_rf_im()

  do i = 0, a1%num_modes
    call sub_f1( a1%rf_re(i), a2%rf_re(i), ua3_re(i) )
    if (i==0) cycle
    call sub_f1( a1%rf_im(i), a2%rf_im(i), ua3_im(i) )
  enddo

end subroutine sub_f1_binary

subroutine sub_f1_unitary( a1, a2 )

  implicit none

  class( field ), intent(in) :: a1
  class( field ), intent(inout) :: a2

  class( ufield ), dimension(:), pointer :: ua2_re => null(), ua2_im => null()
  integer :: i

  ua2_re => a2%get_rf_re()
  ua2_im => a2%get_rf_im()

  do i = 0, a1%num_modes
    call sub_f1( a1%rf_re(i), ua2_re(i) )
    if (i==0) cycle
    call sub_f1( a1%rf_im(i), ua2_im(i) )
  enddo

end subroutine sub_f1_unitary

subroutine dot_f1_unitary( a1, a2 )

  implicit none

  real, intent(in) :: a1
  class( field ), intent(inout) :: a2

  class( ufield ), dimension(:), pointer :: ua2_re => null(), ua2_im => null()
  integer :: i

  ua2_re => a2%get_rf_re()
  ua2_im => a2%get_rf_im()

  do i = 0, a2%num_modes
    call dot_f1( a1, ua2_re(i) )
    if (i==0) cycle
    call dot_f1( a1, ua2_im(i) )
  enddo

end subroutine dot_f1_unitary

subroutine dot_f1_unitary_dim( a1, a2, dim )

  implicit none

  real, intent(in) :: a1
  class( field ), intent(inout) :: a2
  integer, intent(in), dimension(:) :: dim

  class( ufield ), dimension(:), pointer :: ua2_re => null(), ua2_im => null()
  integer :: i

  ua2_re => a2%get_rf_re()
  ua2_im => a2%get_rf_im()

  do i = 0, a2%num_modes
    call dot_f1( a1, ua2_re(i), dim )
    if (i==0) cycle
    call dot_f1( a1, ua2_im(i), dim )
  enddo

end subroutine dot_f1_unitary_dim

subroutine add_f2_binary( a1, a2, a3 )

  implicit none

  class( field ), intent(in) :: a1, a2
  class( field ), intent(inout) :: a3

  class( ufield ), dimension(:), pointer :: ua3_re => null(), ua3_im => null()
  integer :: i

  ua3_re => a3%get_rf_re()
  ua3_im => a3%get_rf_im()

  do i = 0, a1%num_modes
    call add_f2( a1%rf_re(i), a2%rf_re(i), ua3_re(i) )
    if (i==0) cycle
    call add_f2( a1%rf_im(i), a2%rf_im(i), ua3_im(i) )
  enddo

end subroutine add_f2_binary

subroutine add_f2_unitary( a1, a2 )

  implicit none

  class( field ), intent(in) :: a1
  class( field ), intent(inout) ::a2

  class( ufield ), dimension(:), pointer :: ua2_re => null(), ua2_im => null()
  integer :: i

  ua2_re => a2%get_rf_re()
  ua2_im => a2%get_rf_im()

  do i = 0, a1%num_modes
    call add_f2( a1%rf_re(i), ua2_re(i) )
    if (i==0) cycle
    call add_f2( a1%rf_im(i), ua2_im(i) )
  enddo

end subroutine add_f2_unitary

subroutine sub_f2_binary( a1, a2, a3 )

  implicit none

  class( field ), intent(in) :: a1, a2
  class( field ), intent(inout) :: a3

  class( ufield ), dimension(:), pointer :: ua3_re => null(), ua3_im => null()
  integer :: i

  ua3_re => a3%get_rf_re()
  ua3_im => a3%get_rf_im()

  do i = 0, a1%num_modes
    call sub_f2( a1%rf_re(i), a2%rf_re(i), ua3_re(i) )
    if (i==0) cycle
    call sub_f2( a1%rf_im(i), a2%rf_im(i), ua3_im(i) )
  enddo

end subroutine sub_f2_binary

subroutine sub_f2_unitary( a1, a2 )

  implicit none

  class( field ), intent(in) :: a1
  class( field ), intent(inout) ::a2

  class( ufield ), dimension(:), pointer :: ua2_re => null(), ua2_im => null()
  integer :: i

  ua2_re => a2%get_rf_re()
  ua2_im => a2%get_rf_im()

  do i = 0, a1%num_modes
    call sub_f2( a1%rf_re(i), ua2_re(i) )
    if (i==0) cycle
    call sub_f2( a1%rf_im(i), ua2_im(i) )
  enddo

end subroutine sub_f2_unitary

! function add_f1_v1( a1, a2 ) result( a3 )

!   implicit none

!   class( field ), intent(in) :: a1
!   class(*), intent(in) :: a2
!   class( field ), allocatable :: a3

!   class( ufield ), dimension(:), pointer :: ua3_re => null(), ua3_im => null()
!   integer :: i

!   allocate(a3)
!   call a3%new(a1)
!   ua3_re => a3%get_rf_re()
!   ua3_im => a3%get_rf_im()

!   select type (a2)
!   type is (real)
!     do i = 0, a1%num_modes
!       ua3_re(i) = a1%rf_re(i) + a2
!       if (i==0) cycle
!       ua3_im(i) = a1%rf_im(i) + a2
!     enddo
!   class is (field)
!     do i = 0, a1%num_modes
!       ua3_re(i) = a1%rf_re(i) + a2%get_rf_re(i)
!       if (i==0) cycle
!       ua3_im(i) = a1%rf_im(i) + a2%get_rf_im(i)
!     enddo
!   class default
!     call write_err( "invalid assignment type!" )
!   end select

! end function add_f1_v1

! function add_f1_v2( a2, a1 ) result( a3 )

!   implicit none

!   class( field ), intent(in) :: a1
!   real, intent(in) :: a2
!   class( field ), allocatable :: a3

!   class( ufield ), dimension(:), pointer :: ua3_re => null(), ua3_im => null()
!   integer :: i

!   allocate(a3)
!   call a3%new(a1)
!   ua3_re => a3%get_rf_re()
!   ua3_im => a3%get_rf_im()

!   do i = 0, a1%num_modes
!     ua3_re(i) = a1%rf_re(i) + a2
!     if (i==0) cycle
!     ua3_im(i) = a1%rf_im(i) + a2
!   enddo

! end function add_f1_v2

! function dot_f1_v1( a1, a2 ) result( a3 )

!   implicit none

!   class( field ), intent(in) :: a1
!   class(*), intent(in) :: a2
!   class( field ), allocatable :: a3

!   class( ufield ), dimension(:), pointer :: ua3_re => null(), ua3_im => null()
!   integer :: i

!   allocate(a3)
!   call a3%new(a1)
!   ua3_re => a3%get_rf_re()
!   ua3_im => a3%get_rf_im()

!   select type (a2)
!   type is (real)
!     do i = 0, a1%num_modes
!       ua3_re(i) = a1%rf_re(i) * a2
!       if (i==0) cycle
!       ua3_im(i) = a1%rf_im(i) * a2
!     enddo
!   class is (field)
!     do i = 0, a1%num_modes
!       ua3_re(i) = a1%rf_re(i) * a2%get_rf_re(i)
!       if (i==0) cycle
!       ua3_im(i) = a1%rf_im(i) * a2%get_rf_im(i)
!     enddo
!   class default
!     call write_err( "invalid assignment type!" )
!   end select

! end function dot_f1_v1

! function dot_f1_v2( a2, a1 ) result( a3 )

!   implicit none

!   class( field ), intent(in) :: a1
!   real, intent(in) :: a2
!   class( field ), allocatable :: a3

!   class( ufield ), dimension(:), pointer :: ua3_re => null(), ua3_im => null()
!   integer :: i

!   allocate(a3)
!   call a3%new(a1)
!   ua3_re => a3%get_rf_re()
!   ua3_im => a3%get_rf_im()

!   do i = 0, a1%num_modes
!     ua3_re(i) = a1%rf_re(i) * a2
!     if (i==0) cycle
!     ua3_im(i) = a1%rf_im(i) * a2
!   enddo

! end function dot_f1_v2

! function sub_f1_v1( a1, a2 ) result( a3 )

!   implicit none

!   class( field ), intent(in) :: a1
!   class(*), intent(in) :: a2
!   class( field ), allocatable :: a3

!   class( ufield ), dimension(:), pointer :: ua3_re => null(), ua3_im => null()
!   integer :: i

!   allocate(a3)
!   call a3%new(a1)
!   ua3_re => a3%get_rf_re()
!   ua3_im => a3%get_rf_im()

!   select type (a2)
!   type is (real)
!     do i = 0, a1%num_modes
!       ua3_re(i) = a1%rf_re(i) - a2
!       if (i==0) cycle
!       ua3_im(i) = a1%rf_im(i) - a2
!     enddo
!   class is (field)
!     do i = 0, a1%num_modes
!       ua3_re(i) = a1%rf_re(i) - a2%get_rf_re(i)
!       if (i==0) cycle
!       ua3_im(i) = a1%rf_im(i) - a2%get_rf_im(i)
!     enddo
!   class default
!     call write_err( "invalid assignment type!" )
!   end select

! end function sub_f1_v1

! function sub_f1_v2( a2, a1 ) result( a3 )

!   implicit none

!   class( field ), intent(in) :: a1
!   real, intent(in) :: a2
!   class( field ), allocatable :: a3

!   class( ufield ), dimension(:), pointer :: ua3_re => null(), ua3_im => null()
!   integer :: i

!   allocate(a3)
!   call a3%new(a1)
!   ua3_re => a3%get_rf_re()
!   ua3_im => a3%get_rf_im()

!   do i = 0, a1%num_modes
!     ua3_re(i) = a1%rf_re(i) - a2
!     if (i==0) cycle
!     ua3_im(i) = a1%rf_im(i) - a2
!   enddo

! end function sub_f1_v2

! subroutine add_f2_dim( a1, a2, a3, dim1, dim2, dim3 )

!   implicit none

!   class( field ), intent(in) :: a1, a2
!   class( field ), intent(inout) :: a3
!   integer, intent(in), dimension(:) :: dim1, dim2, dim3

!   class( ufield ), dimension(:), pointer :: ua3_re => null(), ua3_im => null()
!   integer :: i

!   ua3_re => a3%get_rf_re()
!   ua3_im => a3%get_rf_im()

!   do i = 0, a1%num_modes
!     call add_f2( a1%rf_re(i), a2%rf_re(i), ua3_re(i), dim1, dim2, dim3 )
!     if (i==0) cycle
!     call add_f2( a1%rf_im(i), a2%rf_im(i), ua3_im(i), dim1, dim2, dim3 )
!   enddo

! end subroutine add_f2_dim



! subroutine sub_f2_dim( a1, a2, a3, dim1, dim2, dim3 )

!   implicit none

!   class( field ), intent(in) :: a1, a2
!   class( field ), intent(inout) :: a3
!   integer, intent(in), dimension(:) :: dim1, dim2, dim3

!   class( ufield ), dimension(:), pointer :: ua3_re => null(), ua3_im => null()
!   integer :: i

!   ua3_re => a3%get_rf_re()
!   ua3_im => a3%get_rf_im()

!   do i = 0, a1%num_modes
!     call sub_f2( a1%rf_re(i), a2%rf_re(i), ua3_re(i), dim1, dim2, dim3 )
!     if (i==0) cycle
!     call sub_f2( a1%rf_im(i), a2%rf_im(i), ua3_im(i), dim1, dim2, dim3 )
!   enddo

! end subroutine sub_f2_dim

! subroutine sub_f2_all( a1, a2, a3 )

!   implicit none

!   class( field ), intent(in) :: a1, a2
!   class( field ), intent(inout) :: a3

!   class( ufield ), dimension(:), pointer :: ua3_re => null(), ua3_im => null()
!   integer :: i

!   ua3_re => a3%get_rf_re()
!   ua3_im => a3%get_rf_im()

!   do i = 0, a1%num_modes
!     call sub_f2( a1%rf_re(i), a2%rf_re(i), ua3_re(i) )
!     if (i==0) cycle
!     call sub_f2( a1%rf_im(i), a2%rf_im(i), ua3_im(i) )
!   enddo

! end subroutine sub_f2_all

! subroutine dot_f2_dim( a1, a2, a3, dim1, dim2, dim3 )

!   implicit none

!   class( field ), intent(in) :: a1, a2
!   class( field ), intent(inout) :: a3
!   integer, intent(in), dimension(:) :: dim1, dim2, dim3

!   class( ufield ), dimension(:), pointer :: ua3_re => null(), ua3_im => null()
!   integer :: i

!   ua3_re => a3%get_rf_re()
!   ua3_im => a3%get_rf_im()

!   do i = 0, a1%num_modes
!     call dot_f2( a1%rf_re(i), a2%rf_re(i), ua3_re(i), dim1, dim2, dim3 )
!     if (i==0) cycle
!     call dot_f2( a1%rf_im(i), a2%rf_im(i), ua3_im(i), dim1, dim2, dim3 )
!   enddo

! end subroutine dot_f2_dim

! subroutine dot_f2_all( a1, a2, a3 )

!   implicit none

!   class( * ), intent(in) :: a1, a2
!   class( field ), intent(inout) :: a3

!   class( ufield ), dimension(:), pointer :: ua3_re => null(), ua3_im => null()
!   integer :: i

!   ua3_re => a3%get_rf_re()
!   ua3_im => a3%get_rf_im()

!   select type (a1)

!   type is (real)

!     select type (a2)

!       type is  (real)

!         do i = 0, a3%num_modes
!           call dot_f2( a1, a2, ua3_re(i) )
!           if (i==0) cycle
!           call dot_f2( a1, a2, ua3_im(i) )
!         enddo

!       class is (field)

!         do i = 0, a3%num_modes
!           call dot_f2( a1, a2%rf_re(i), ua3_re(i) )
!           if (i==0) cycle
!           call dot_f2( a1, a2%rf_im(i), ua3_im(i) )
!         enddo

!       class default

!         call write_err( "invalid assignment type!" )

!     end select

!   class is (field)

!     select type (a2)

!       type is (real)

!         do i = 0, a3%num_modes
!           call dot_f2( a1%rf_re(i), a2, ua3_re(i) )
!           if (i==0) cycle
!           call dot_f2( a1%rf_im(i), a2, ua3_im(i) )
!         enddo

!       class is (field)

!         do i = 0, a3%num_modes
!           call dot_f2( a1%rf_re(i), a2%rf_re(i), ua3_re(i) )
!           if (i==0) cycle
!           call dot_f2( a1%rf_im(i), a2%rf_im(i), ua3_im(i) )
!         enddo

!       class default

!         call write_err( "invalid assignment type!" )

!     end select

!   class default
!     call write_err( "invalid assignment type!" )
!   end select

! end subroutine dot_f2_all

! function add_f2_v1( a1, a2 ) result( a3 )

!   implicit none

!   class( field ), intent(in) :: a1
!   class(*), intent(in) :: a2
!   class( field ), allocatable :: a3

!   class( ufield ), dimension(:), pointer :: ua3_re => null(), ua3_im => null()
!   integer :: i

!   allocate(a3)
!   call a3%new(a1)
!   ua3_re => a3%get_rf_re()
!   ua3_im => a3%get_rf_im()

!   select type (a2)
!   type is (real)
!     do i = 0, a1%num_modes
!       call ua3_re(i)%as( a1%rf_re(i) .add. a2 )
!       if (i==0) cycle
!       call ua3_im(i)%as( a1%rf_im(i) .add. a2 )
!     enddo
!   class is (field)
!     do i = 0, a1%num_modes
!       call ua3_re(i)%as( a1%rf_re(i) .add. a2%get_rf_re(i) )
!       if (i==0) cycle
!       call ua3_im(i)%as( a1%rf_im(i) .add. a2%get_rf_im(i) )
!     enddo
!   class default
!     call write_err( "invalid assignment type!" )
!   end select

! end function add_f2_v1

! function add_f2_v2( a2, a1 ) result( a3 )

!   implicit none

!   class( field ), intent(in) :: a1
!   real, intent(in) :: a2
!   class( field ), allocatable :: a3

!   class( ufield ), dimension(:), pointer :: ua3_re => null(), ua3_im => null()
!   integer :: i

!   allocate(a3)
!   call a3%new(a1)
!   ua3_re => a3%get_rf_re()
!   ua3_im => a3%get_rf_im()

!   do i = 0, a1%num_modes
!     call ua3_re(i)%as( a1%rf_re(i) .add. a2 )
!     if (i==0) cycle
!     call ua3_im(i)%as( a1%rf_im(i) .add. a2 )
!   enddo

! end function add_f2_v2

! function dot_f2_v1( a1, a2 ) result( a3 )

!   implicit none

!   class( field ), intent(in) :: a1
!   class(*), intent(in) :: a2
!   class( field ), allocatable :: a3

!   class( ufield ), dimension(:), pointer :: ua3_re => null(), ua3_im => null()
!   integer :: i

!   allocate(a3)
!   call a3%new(a1)
!   ua3_re => a3%get_rf_re()
!   ua3_im => a3%get_rf_im()

!   select type (a2)
!   type is (real)
!     do i = 0, a1%num_modes
!       call ua3_re(i)%as( a1%rf_re(i) .dot. a2 )
!       if (i==0) cycle
!       call ua3_im(i)%as( a1%rf_im(i) .dot. a2 )
!     enddo
!   class is (field)
!     do i = 0, a1%num_modes
!       call ua3_re(i)%as( a1%rf_re(i) .dot. a2%get_rf_re(i) )
!       if (i==0) cycle
!       call ua3_im(i)%as( a1%rf_im(i) .dot. a2%get_rf_im(i) )
!     enddo
!   class default
!     call write_err( "invalid assignment type!" )
!   end select

! end function dot_f2_v1

! function dot_f2_v2( a2, a1 ) result( a3 )

!   implicit none

!   class( field ), intent(in) :: a1
!   real, intent(in) :: a2
!   class( field ), allocatable :: a3

!   class( ufield ), dimension(:), pointer :: ua3_re => null(), ua3_im => null()
!   integer :: i

!   allocate(a3)
!   call a3%new(a1)
!   ua3_re => a3%get_rf_re()
!   ua3_im => a3%get_rf_im()

!   do i = 0, a1%num_modes
!     call ua3_re(i)%as( a1%rf_re(i) .dot. a2 )
!     if (i==0) cycle
!     call ua3_im(i)%as( a1%rf_im(i) .dot. a2 )
!   enddo

! end function dot_f2_v2

! function sub_f2_v1( a1, a2 ) result( a3 )

!   implicit none

!   class( field ), intent(in) :: a1
!   class(*), intent(in) :: a2
!   class( field ), allocatable :: a3

!   class( ufield ), dimension(:), pointer :: ua3_re => null(), ua3_im => null()
!   integer :: i

!   allocate(a3)
!   call a3%new(a1)
!   ua3_re => a3%get_rf_re()
!   ua3_im => a3%get_rf_im()

!   select type (a2)
!   type is (real)
!     do i = 0, a1%num_modes
!       call ua3_re(i)%as( a1%rf_re(i) .sub. a2 )
!       if (i==0) cycle
!       call ua3_im(i)%as( a1%rf_im(i) .sub. a2 )
!     enddo
!   class is (field)
!     do i = 0, a1%num_modes
!       call ua3_re(i)%as( a1%rf_re(i) .sub. a2%get_rf_re(i) )
!       if (i==0) cycle
!       call ua3_im(i)%as( a1%rf_im(i) .sub. a2%get_rf_im(i) )
!     enddo
!   class default
!     call write_err( "invalid assignment type!" )
!   end select

! end function sub_f2_v1

! function sub_f2_v2( a2, a1 ) result( a3 )

!   implicit none

!   class( field ), intent(in) :: a1
!   real, intent(in) :: a2
!   class( field ), allocatable :: a3

!   class( ufield ), dimension(:), pointer :: ua3_re => null(), ua3_im => null()
!   integer :: i

!   allocate(a3)
!   call a3%new(a1)
!   ua3_re => a3%get_rf_re()
!   ua3_im => a3%get_rf_im()

!   do i = 0, a1%num_modes
!     call ua3_re(i)%as( a2 .sub. a1%rf_re(i) )
!     if (i==0) cycle
!     call ua3_im(i)%as( a2 .sub. a1%rf_im(i) )
!   enddo

! end function sub_f2_v2

end module field_class
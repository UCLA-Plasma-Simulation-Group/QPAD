module fdist3d_file_class

use parallel_module
use options_class
use input_class
use iso_c_binding
use fdist3d_class
use part3d_class
use sysutil_module
use param
use hdf5
use hdf5io_class

implicit none

private

type, extends(fdist3d) :: fdist3d_file

  ! private

  ! File name
  character(len=:), allocatable :: filename

  ! Beam center
  real, dimension(p_x_dim) :: beam_ctr

  ! Beam center in the file
  real, dimension(p_x_dim) :: file_ctr

  ! Length unit converting factor
  real :: xconv_fac

  ! Charge unit converting factor
  real :: qconv_fac

  contains

  procedure :: new => init_fdist3d_file
  procedure :: del => end_fdist3d_file
  procedure :: inject => inject_fdist3d_file

end type fdist3d_file

character(len=32), save :: cls_name = 'fdist3d_file'
integer, save :: cls_level = 2

public :: fdist3d_file

contains

subroutine init_fdist3d_file( this, input, opts, sect_id )

  implicit none

  class( fdist3d_file ), intent(inout) :: this
  type( input_json ), intent(inout) :: input
  type( options ), intent(in) :: opts
  integer, intent(in) :: sect_id

  real :: xtra
  integer :: npmax_tmp
  integer(kind=LG) :: npmax_guess
  character(len=20) :: sect_name
  character(len=:), allocatable :: read_str
  character(len=18), save :: sname = 'init_fdist3d_file'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  this%noff_r = opts%get_noff(1)
  this%noff_z = opts%get_noff(2)
  this%nr     = opts%get_nd(1)
  this%nz     = opts%get_nd(2)
  this%nrp    = opts%get_ndp(1)
  this%nzp    = opts%get_ndp(2)
  this%dr     = opts%get_dr()
  this%dz     = opts%get_dxi()

  sect_name = 'beam(' // num2str(sect_id) // ')'

  call input%get( trim(sect_name) // '.filename', this%filename )
  call input%get( 'simulation.box.xi(1)', this%z0 )
  call input%get( trim(sect_name) // '.beam_center(1)', this%beam_ctr(1) )
  call input%get( trim(sect_name) // '.beam_center(2)', this%beam_ctr(2) )
  call input%get( trim(sect_name) // '.beam_center(3)', this%beam_ctr(3) )
  this%beam_ctr(3) = this%beam_ctr(3) - this%z0
  call input%get( trim(sect_name) // '.file_center(1)', this%file_ctr(1) )
  call input%get( trim(sect_name) // '.file_center(2)', this%file_ctr(2) )
  call input%get( trim(sect_name) // '.file_center(3)', this%file_ctr(3) )

  this%xconv_fac = 1.0
  if ( input%found( trim(sect_name) // '.length_conv_fac' ) ) then
    call input%get( trim(sect_name) // '.length_conv_fac', this%xconv_fac )
  endif

  this%qconv_fac = 1.0
  if ( input%found( trim(sect_name) // '.charge_conv_fac' ) ) then
    call input%get( trim(sect_name) // '.charge_conv_fac', this%qconv_fac )
  endif

  call input%get( trim(sect_name) // '.q', this%qm )
  call input%get( trim(sect_name) // '.m', this%qbm )
  this%qbm = this%qm / this%qbm

  this%evol = .true.
  if ( input%found( trim(sect_name) // '.evolution' ) ) then
    call input%get( trim(sect_name) // '.evolution', this%evol )
  endif

  this%quiet = .false.
  this%has_spin = .false.
  if ( input%found( trim(sect_name) // '.has_spin' ) ) then
    call input%get( trim(sect_name) // '.has_spin', this%has_spin )
    if ( this%has_spin ) then
      call input%get( trim(sect_name) // '.anom_mag_moment', this%amm )
    endif
  endif

  ! if npmax is not given, guess a value for npmax
  xtra = 1.5
  ! this guess value must be smaller than the practical value, so the buffer reallocation
  ! will be invoked. Need a more smarter estimation algorithm.
  npmax_guess = p_npmax_min
  this%npmax = int( npmax_guess * xtra, kind=LG )

  ! if npmax is given, set it as the maximum of the given value and p_npmax_min
  if ( input%found( trim(sect_name) // '.npmax' ) ) then
    call input%get( trim(sect_name) // '.npmax', npmax_tmp )
    this%npmax = int( npmax_tmp, kind=LG )
    if ( this%npmax < p_npmax_min ) then
      call write_wrn( 'The npmax may be too small for the initialization. It has &
        &been automatically changed to ' // num2str(p_npmax_min) // '.' )
      this%npmax = int( p_npmax_min, kind=LG )
    endif
  endif

  call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine init_fdist3d_file

subroutine end_fdist3d_file( this )

  implicit none

  class( fdist3d_file ), intent(inout) :: this
  character(len=18), save :: sname = 'end_fdist3d_file'

  call write_dbg( cls_name, sname, cls_level, 'starts' )
  ! placeholder
  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine end_fdist3d_file

subroutine inject_fdist3d_file( this, part )

  use iso_c_binding
  implicit none

  class( fdist3d_file ), intent(inout) :: this
  class( part3d ), intent(inout) :: part
  ! local
  real :: rbuf, ratio
  real, dimension(:,:), allocatable :: xbuf, pbuf, sbuf
  real, dimension(:), allocatable :: qbuf
  real, dimension(2) :: edge_r, edge_z
  integer :: i, ierr
  integer(kind=LG) :: ip, ptrcur
  integer, parameter :: BLOCK_SIZE = 65536
  integer, parameter :: real_kind_8 = kind(1.0d0)
  integer(hsize_t), dimension(2) :: dims, maxdims, chunk_size, offset
  integer(hid_t) :: file_id, grp_id, treal, dtype_id
  integer(hid_t), dimension(p_x_dim) :: xdset_id
  integer(hid_t), dimension(p_p_dim) :: pdset_id
  integer(hid_t), dimension(p_s_dim) :: sdset_id
  integer(hid_t) :: qdset_id, fspace_id, mspace_id
  character(len=18), save :: sname = 'inject_fdist3d_file'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  allocate( xbuf( BLOCK_SIZE, p_x_dim ) )
  allocate( pbuf( BLOCK_SIZE, p_p_dim ) )
  allocate( qbuf( BLOCK_SIZE) )
  if ( this%has_spin ) allocate( sbuf( BLOCK_SIZE, p_s_dim ) )

  treal = detect_precision()

  edge_r(p_lower) = this%noff_r * this%dr
  edge_z(p_lower) = this%noff_z * this%dz
  edge_r(p_upper) = edge_r(p_lower) + this%nrp * this%dr
  edge_z(p_upper) = edge_z(p_lower) + this%nzp * this%dz

  call h5open_f( ierr )
  call h5fopen_f( this%filename, H5F_ACC_RDONLY_F, file_id, ierr )
  call h5gopen_f( file_id, '/', grp_id, ierr )

  ! get metadata of dataset to be read
  call h5dopen_f( grp_id, 'x1', xdset_id(1), ierr )
  call h5dget_space_f( xdset_id(1), fspace_id, ierr )
  call h5dget_type_f( xdset_id(1), dtype_id, ierr )
  call h5sget_simple_extent_dims_f( fspace_id, dims, maxdims, ierr )
  call h5sclose_f( fspace_id, ierr )
  call h5dclose_f( xdset_id(1), ierr )

  ! create memory dataspace
  call h5screate_simple_f( 2, (/int(BLOCK_SIZE, hsize_t), 1_hsize_t/), &
    mspace_id, ierr )

  ! open all the datasets to be read
  call h5dopen_f( grp_id, 'x1', xdset_id(1), ierr )
  call h5dopen_f( grp_id, 'x2', xdset_id(2), ierr )
  call h5dopen_f( grp_id, 'x3', xdset_id(3), ierr )
  call h5dopen_f( grp_id, 'p1', pdset_id(1), ierr )
  call h5dopen_f( grp_id, 'p2', pdset_id(2), ierr )
  call h5dopen_f( grp_id, 'p3', pdset_id(3), ierr )
  call h5dopen_f( grp_id, 'q', qdset_id, ierr )
  if ( this%has_spin ) then
    call h5dopen_f( grp_id, 's1', sdset_id(1), ierr )
    call h5dopen_f( grp_id, 's2', sdset_id(2), ierr )
    call h5dopen_f( grp_id, 's3', sdset_id(3), ierr )
  endif

  offset = 0
  chunk_size = (/0, 1/)
  ip = 0; part%npp = 0
  do ptrcur = 1, dims(1), BLOCK_SIZE

    ! check if last copy of table and set np
    if( ptrcur + BLOCK_SIZE > dims(1) ) then
      chunk_size(1) = dims(1) - ptrcur + 1
      call h5sset_extent_simple_f( mspace_id, 2, (/int(chunk_size(1), hsize_t), 1_hsize_t/), &
       (/int(chunk_size(1), hsize_t), 1_hsize_t/), ierr )
    else
      chunk_size(1) = BLOCK_SIZE
    endif

    if ( part%npp + chunk_size(1) > part%npmax ) then
      ratio = real(part%npp + chunk_size(1)) / part%npmax
      call part%realloc( ratio = p_buf_incr * ratio, buf_type = 'particle' )
    endif

    ! read chunk from datasets x1, x2, x3
    do i = 1, p_x_dim
      call h5dget_space_f( xdset_id(i), fspace_id, ierr )
      call h5sselect_hyperslab_f( fspace_id, H5S_SELECT_SET_F, offset, &
         chunk_size, ierr) 
      call h5dread_f( xdset_id(i), h5kind_to_type(real_kind_8, H5_REAL_KIND), &
         xbuf(1,i), chunk_size, ierr, mspace_id, fspace_id )
      call h5sclose_f( fspace_id, ierr )
    enddo

    ! read chunk from datasets p1, p2, p3
    do i = 1, p_p_dim
       call h5dget_space_f( pdset_id(i), fspace_id, ierr )
       call h5sselect_hyperslab_f( fspace_id, H5S_SELECT_SET_F, offset, &
          chunk_size, ierr) 
       call h5dread_f( pdset_id(i), h5kind_to_type(real_kind_8, H5_REAL_KIND), &
          pbuf(1,i), chunk_size, ierr, mspace_id, fspace_id )
       call h5sclose_f( fspace_id, ierr )
    enddo

    ! read chunk from dataset q
    call h5dget_space_f( qdset_id, fspace_id, ierr )
    call h5sselect_hyperslab_f( fspace_id, H5S_SELECT_SET_F, offset, &
       chunk_size, ierr) 
    call h5dread_f( qdset_id, h5kind_to_type(real_kind_8, H5_REAL_KIND), &
       qbuf, chunk_size, ierr, mspace_id, fspace_id )
    call h5sclose_f( fspace_id, ierr )

    ! read chunk from datasets s1, s2, s3
    if ( this%has_spin ) then
    do i = 1, p_s_dim
       call h5dget_space_f( sdset_id(i), fspace_id, ierr )
       call h5sselect_hyperslab_f( fspace_id, H5S_SELECT_SET_F, offset, &
          chunk_size, ierr) 
       call h5dread_f( sdset_id(i), h5kind_to_type(real_kind_8, H5_REAL_KIND), &
          sbuf(1,i), chunk_size, ierr, mspace_id, fspace_id )
       call h5sclose_f( fspace_id, ierr )
    enddo
    endif

    ! convert and store particle data
    do i = 1, chunk_size(1)

       ! centralize particle position
       xbuf(i,1) = ( xbuf(i,1) - this%file_ctr(1) ) * this%xconv_fac
       xbuf(i,2) = ( xbuf(i,2) - this%file_ctr(2) ) * this%xconv_fac
       xbuf(i,3) = ( xbuf(i,3) - this%file_ctr(3) ) * this%xconv_fac

       ! shift in QPAD space
       xbuf(i,1) = this%beam_ctr(1) + xbuf(i,1)
       xbuf(i,2) = this%beam_ctr(2) + xbuf(i,2)
       xbuf(i,3) = this%beam_ctr(3) - xbuf(i,3)

       rbuf = sqrt( xbuf(i,1)**2 + xbuf(i,2)**2 )
       if ( rbuf >= edge_r(p_lower) .and. rbuf < edge_r(p_upper) .and. &
            xbuf(i,3) >= edge_z(p_lower) .and. xbuf(i,3) < edge_z(p_upper) ) then
          ! pp = pp + 1
          part%npp = part%npp + 1
          ip = part%npp
          part%x(1,ip) = xbuf(i,1)
          part%x(2,ip) = xbuf(i,2)
          part%x(3,ip) = xbuf(i,3)
          part%p(1,ip) = pbuf(i,1)
          part%p(2,ip) = pbuf(i,2)
          part%p(3,ip) = pbuf(i,3)
          part%q(ip) = qbuf(i) * this%qconv_fac
          if ( this%has_spin ) then
             part%s(1,ip) = sbuf(i,1)
             part%s(2,ip) = sbuf(i,2)
             part%s(3,ip) = sbuf(i,3)
          endif
       endif

    enddo

    offset(1) = offset(1) + chunk_size(1)

  enddo

  call h5sclose_f( mspace_id, ierr )
  call h5dclose_f( xdset_id(1), ierr )
  call h5dclose_f( xdset_id(2), ierr )
  call h5dclose_f( xdset_id(3), ierr )
  call h5dclose_f( pdset_id(1), ierr )
  call h5dclose_f( pdset_id(2), ierr )
  call h5dclose_f( pdset_id(3), ierr )
  call h5dclose_f( qdset_id, ierr )
  if ( this%has_spin ) then
    call h5dclose_f( sdset_id(1), ierr )
    call h5dclose_f( sdset_id(2), ierr )
    call h5dclose_f( sdset_id(3), ierr )
  endif

  call h5gclose_f( grp_id, ierr )
  call h5fclose_f( file_id, ierr )
  call h5close_f( ierr )

  deallocate( xbuf, pbuf, qbuf )
  if ( this%has_spin ) deallocate( sbuf )

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine inject_fdist3d_file

end module fdist3d_file_class
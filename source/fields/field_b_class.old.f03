module field_b_class

use field_class
use field_src_class
use field_solver_class
use ufield_class
use param
use system
use parallel_pipe_class
use grid_class
use debug_tool
use mpi

implicit none

private

character(len=20), parameter :: cls_name = "field_b"
integer, parameter :: cls_level = 3

public :: field_b

type, extends( field ) :: field_b

  ! private

  class( field_solver ), dimension(:), pointer :: solver_bz => null()
  class( field_solver ), dimension(:), pointer :: solver_bperp => null()
  class( field_solver ), dimension(:), pointer :: solver_bperp_iter => null()

  real, dimension(:), pointer :: buf_re => null(), buf_im => null()
  real, dimension(:), pointer :: buf => null()

  contains

  generic :: new => init_field_b
  procedure :: del => end_field_b
  ! generic :: read_input => read_input_field_b
  generic :: solve => solve_field_bz, solve_field_bperp, solve_field_bperp_iter
  generic :: solve_old => solve_field_bperp_old

  procedure, private :: init_field_b
  procedure, private :: end_field_b
  procedure, private :: set_source_bz
  procedure, private :: set_source_bperp
  procedure, private :: set_source_bperp_old
  procedure, private :: set_source_bperp_iter
  procedure, private :: get_solution_bz
  procedure, private :: get_solution_bperp
  procedure, private :: get_solution_bperp_old
  procedure, private :: get_solution_bperp_iter
  procedure, private :: solve_field_bz
  procedure, private :: solve_field_bperp
  procedure, private :: solve_field_bperp_old
  procedure, private :: solve_field_bperp_iter

end type field_b

contains

subroutine init_field_b( this, pp, gp, dr, dxi, num_modes, part_shape, entity, iter_tol )

  implicit none

  class( field_b ), intent(inout) :: this
  class( parallel_pipe ), intent(in), pointer :: pp
  class( grid ), intent(in), pointer :: gp
  integer, intent(in) :: num_modes, part_shape, entity
  real, intent(in) :: dr, dxi, iter_tol

  integer, dimension(2,2) :: gc_num
  integer :: dim, i, nrp
  character(len=20), save :: sname = "init_field_b"

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  nrp = gp%get_ndp(1)

  select case ( part_shape )
  
  case ( p_ps_linear )
  
    gc_num(:,1) = (/1, 1/)
    gc_num(:,2) = (/0, 1/)
  
  case ( p_ps_quadratic )

    call write_err( "Quadratic particle shape not implemented." )

  case default

    call write_err( "Invalid particle shape." )

  end select

  dim = 3
  ! call initialization routine of the parent class
  call this%field%new( pp, gp, dim, dr, dxi, num_modes, gc_num, entity )

  ! initialize solver
  select case ( entity )

  case ( p_entity_plasma )

    allocate( this%solver_bz( 0:num_modes ) )
    allocate( this%solver_bperp_iter( 0:num_modes ) )
    do i = 0, num_modes
      call this%solver_bz(i)%new( pp, gp, i, dr, kind=p_fk_bz, &
        stype=p_hypre_cycred, tol=iter_tol )
      call this%solver_bperp_iter(i)%new( pp, gp, i, dr, kind=p_fk_bperp_iter, &
        stype=p_hypre_amg, tol=iter_tol )
    enddo

    allocate( this%buf(4*nrp) )
    allocate( this%buf_re(nrp), this%buf_im(nrp) )

  case ( p_entity_beam_old )

    allocate( this%solver_bperp( 0:num_modes ) )
    do i = 0, num_modes
      call this%solver_bperp(i)%new( pp, gp, i, dr, kind=p_fk_bperp_old, &
        stype=p_hypre_amg, tol=iter_tol )
    enddo

    allocate( this%buf( nrp*4 ) )

  case ( p_entity_beam )

    allocate( this%solver_bperp( 0:num_modes ) )
    do i = 0, num_modes
      call this%solver_bperp(i)%new( pp, gp, i, dr, kind=p_fk_bperp, &
        stype=p_hypre_cycred, tol=iter_tol )
    enddo

    allocate( this%buf_re(nrp), this%buf_im(nrp) )

  case default

    call write_err( 'Invalid field entity type.' )

  end select

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine init_field_b

subroutine end_field_b( this )

  implicit none

  class( field_b ), intent(inout) :: this

  integer :: i
  character(len=20), save :: sname = 'end_field_b'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  select case ( this%entity )
  case ( p_entity_plasma )
    do i = 0, this%num_modes
      call this%solver_bz(i)%del()
      call this%solver_bperp_iter(i)%del()
    enddo
    deallocate( this%solver_bz )
    deallocate( this%solver_bperp_iter )
  case ( p_entity_beam )
    do i = 0, this%num_modes
      call this%solver_bperp(i)%del()
    enddo
    deallocate( this%solver_bperp )
  end select

  if ( associated( this%buf_re ) ) deallocate( this%buf_re )
  if ( associated( this%buf_im ) ) deallocate( this%buf_im )
  if ( associated( this%buf ) ) deallocate( this%buf )

  call this%field%del()

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine end_field_b

subroutine set_source_bz( this, mode, jay_re, jay_im )

  implicit none

  class( field_b ), intent(inout) :: this
  class( ufield ), intent(in) :: jay_re
  class( ufield ), intent(in), optional :: jay_im
  integer, intent(in) :: mode

  integer :: i, nrp, noff, idproc, nvp, ierr
  real, dimension(:,:), pointer :: f1_re => null(), f1_im => null()
  real :: idrh, idr, k0, a1, a2, a3, b
  character(len=32) :: filename
  character(len=20), save :: sname = 'set_source_bz'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  nrp    = jay_re%get_ndp(1)
  idr    = 1.0 / this%dr
  idrh   = 0.5 * idr
  noff   = jay_re%get_noff(1)
  nvp    = jay_re%pp%getlnvp()
  idproc = jay_re%pp%getlidproc()
  
  f1_re => jay_re%get_f1()
  if ( present(jay_im) ) then
    f1_im => jay_im%get_f1()
  endif

  this%buf_re = 0.0
  if ( present(jay_im) ) this%buf_im = 0.0
  if ( mode == 0 ) then
    
    local_sum = 0.0
    do i = 1, nrp

      k0 = real(i+noff) - 0.5

      a1 =  idrh * (k0-0.5) / k0
      a2 = -idrh / k0
      a3 = -idrh * (k0+0.5) / k0

      this%buf_re(i) = a1 * f1_re(2,i-1) + a2 * f1_re(2,i) + a3 * f1_re(2,i+1)

    enddo

    ! calculate the derivatives at the boundary and axis
    if ( idproc == 0 ) then
      this%buf_re(1) = -idr * ( f1_re(2,1) + f1_re(2,2) )
    endif
    if ( idproc == nvp-1 ) then
      a2 = -idr * (nrp+noff+0.5) / (nrp+noff-0.5)
      a3 = idr**2 * (nrp+noff) / (nrp+noff-0.5)
      this%buf_re(nrp) = idr * f1_re(2,nrp-1) + a2 * f1_re(2,nrp)
    endif

  elseif ( mode > 0 .and. present( jay_im ) ) then
    
    do i = 1, nrp

      k0 = real(i+noff) - 0.5

      a1 =  idrh * (k0-0.5) / k0
      a2 = -idrh / k0
      a3 = -idrh * (k0+0.5) / k0
      b  =  idr * real(mode) / k0

      this%buf_re(i) = a1 * f1_re(2,i-1) + a2 * f1_re(2,i) + a3 * f1_re(2,i+1) - &
                        b * f1_im(1,i)
      this%buf_im(i) = a1 * f1_im(2,i-1) + a2 * f1_im(2,i) + a3 * f1_im(2,i+1) + &
                        b * f1_re(1,i)

    enddo

    ! calculate the derivatives at the boundary and axis
    if ( idproc == 0 ) then
      this%buf_re(1) = -idr * ( f1_re(2,1) + f1_re(2,2) + 2.0 * real(mode) * f1_im(1,1) )
      this%buf_im(1) = -idr * ( f1_im(2,1) + f1_im(2,2) - 2.0 * real(mode) * f1_re(1,1) )
    endif
    if ( idproc == nvp-1 ) then
      k0 = real(nrp+noff) - 0.5
      a2 = -idr * (k0+1.0) / k0
      b  =  idr * real(mode) / k0
      this%buf_re(nrp) = idr * f1_re(2,nrp-1) + a2 * f1_re(2,nrp) - b * f1_im(1,nrp)
      this%buf_im(nrp) = idr * f1_im(2,nrp-1) + a2 * f1_im(2,nrp) + b * f1_re(1,nrp)
    endif

  else

    call write_err( 'Invalid input arguments!' )

  endif

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine set_source_bz

subroutine set_source_bperp_old( this, mode, q_re, q_im )

  implicit none

  class( field_b ), intent(inout) :: this
  class( ufield ), intent(in) :: q_re
  class( ufield ), intent(in), optional :: q_im
  integer, intent(in) :: mode

  integer :: i, nrp, idproc, nvp, noff
  real, dimension(:,:), pointer :: f1_re => null(), f1_im => null()
  real :: idrh, idr, a1, a2, a3, b, ir
  character(len=20), save :: sname = 'set_source_bperp_old'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  idproc = q_re%pp%getlidproc()
  nvp = q_re%pp%getlnvp()
  nrp = q_re%get_ndp(1)
  noff = q_re%get_noff(1)
  idr = 1.0 / this%dr
  idrh = 0.5 * idr
  
  f1_re => q_re%get_f1()
  if ( present(q_im) ) then
    f1_im => q_im%get_f1()
  endif

  this%buf = 0.0
  if ( mode == 0 ) then
    
    do i = 1, nrp

      ! Re(Br)
      this%buf(4*i-3) = 0.0
      ! Im(Br)
      this%buf(4*i-2) = 0.0
      ! Re(Bphi)
      this%buf(4*i-1) = idrh * ( f1_re(1,i+1)-f1_re(1,i-1) )
      ! Im(Bphi)
      this%buf(4*i) = 0.0

    enddo

    ! calculate the derivatives at the boundary and axis
    if ( idproc == 0 ) then
      ! this%buf(3) = idr * ( f1_re(1,2)-f1_re(1,1) )
      this%buf(3) = idrh * ( -3.0 * f1_re(1,1) + 4.0 * f1_re(1,2) - f1_re(1,3) )
    endif
    if ( idproc == nvp-1 ) then
      this%buf(4*nrp-1) = idrh * ( 3.0 * f1_re(1,nrp) - 4.0 * f1_re(1,nrp-1) + f1_re(1,nrp-2) )
    endif

  elseif ( mode > 0 .and. present( q_im ) ) then
    
    do i = 1, nrp

      ir = idr / (real(i+noff)-0.5)
      this%buf(4*i-3) =  mode * f1_im(1,i) * ir
      this%buf(4*i-2) = -mode * f1_re(1,i) * ir
      this%buf(4*i-1) = idrh * ( f1_re(1,i+1)-f1_re(1,i-1) )
      this%buf(4*i)   = idrh * ( f1_im(1,i+1)-f1_im(1,i-1) )
      
    enddo

    ! calculate the derivatives at the boundary and axis
    if ( idproc == 0 ) then
      this%buf(3) = idrh * ( -3.0 * f1_re(1,1) + 4.0 * f1_re(1,2) - f1_re(1,3) )
      this%buf(4) = idrh * ( -3.0 * f1_im(1,1) + 4.0 * f1_im(1,2) - f1_im(1,3) )
    endif
    if ( idproc == nvp-1 ) then
      this%buf(4*nrp-1) = idrh * ( 3.0 * f1_re(1,nrp) - 4.0 * f1_re(1,nrp-1) + f1_re(1,nrp-2) )
      this%buf(4*nrp)   = idrh * ( 3.0 * f1_im(1,nrp) - 4.0 * f1_im(1,nrp-1) + f1_im(1,nrp-2) )
    endif

  else

    call write_err( 'Invalid input arguments!' )

  endif

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine set_source_bperp_old

subroutine set_source_bperp( this, mode, q_re, q_im )

  implicit none

  class( field_b ), intent(inout) :: this
  class( ufield ), intent(in) :: q_re
  class( ufield ), intent(in), optional :: q_im
  integer, intent(in) :: mode

  integer :: i, nrp, noff, dtype, ierr, comm
  real, dimension(:,:), pointer :: f1_re => null(), f1_im => null()
  character(len=20), save :: sname = 'set_source_bperp'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  nrp   = q_re%get_ndp(1)
  noff  = q_re%get_noff(1)
  
  f1_re => q_re%get_f1()
  if ( present(q_im) ) then
    f1_im => q_im%get_f1()
  endif

  this%buf_re = 0.0
  if ( present(q_im) ) this%buf_im = 0.0
  if ( mode == 0 ) then
    local_sum = 0.0
    do i = 1, nrp
      this%buf_re(i) = -1.0 * f1_re(1,i)
    enddo
  elseif ( mode > 0 .and. present(q_im) ) then
    do i = 1, nrp
      this%buf_re(i) = -1.0 * f1_re(1,i)
      this%buf_im(i) = -1.0 * f1_im(1,i)
    enddo
  else
    call write_err( 'Invalid input arguments!' )
  endif

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine set_source_bperp

subroutine set_source_bperp_iter( this, mode, djdxi_re, jay_re, djdxi_im, jay_im )

  implicit none

  class( field_b ), intent(inout) :: this
  class( ufield ), intent(in) :: djdxi_re, jay_re
  class( ufield ), intent(in), optional :: djdxi_im, jay_im
  integer, intent(in) :: mode

  integer :: i, nrp, nvp, idproc, noff, ierr
  real, dimension(:,:), pointer :: f1_re => null(), f1_im => null()
  real, dimension(:,:), pointer :: f2_re => null(), f2_im => null()
  real, dimension(:,:), pointer :: f3_re => null(), f3_im => null()
  real :: idrh, idr, a1, a2, a3, b, ir
  character(len=20), save :: sname = 'set_source_bperp_iter'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  nvp = jay_re%pp%getlnvp()
  idproc = jay_re%pp%getlidproc()
  nrp = jay_re%get_ndp(1)
  noff = jay_re%get_noff(1)
  idr = 1.0 / this%dr
  idrh = 0.5 * idr
  
  f1_re => djdxi_re%get_f1()
  f2_re => jay_re%get_f1()
  f3_re => this%rf_re(mode)%get_f1()

  if ( present(djdxi_im) .and. present(jay_im) ) then
    f1_im => djdxi_im%get_f1()
    f2_im => jay_im%get_f1()
    f3_im => this%rf_im(mode)%get_f1()
  endif

  this%buf = 0.0
  if ( mode == 0 ) then
    
    local_sum = 0.0
    do i = 1, nrp

      ! Re(Br)
      this%buf(4*i-3) = -f1_re(2,i) - f3_re(1,i)
      ! Im(Br)
      this%buf(4*i-2) = 0.0
      ! Re(Bphi)
      this%buf(4*i-1) = f1_re(1,i) + idrh * ( f2_re(3,i+1)-f2_re(3,i-1) ) - f3_re(2,i)
      ! Im(Bphi)
      this%buf(4*i) = 0.0

    enddo

    ! calculate the derivatives at the boundary and axis
    if ( idproc == 0 ) then
      this%buf(1) = -f1_re(2,1) - f3_re(1,1)
      this%buf(2) = 0.0
      this%buf(3) = f1_re(1,1) + idrh * ( -3.0 * f2_re(3,1) + 4.0 * f2_re(3,2) - f2_re(3,3) ) - f3_re(2,1)
      this%buf(4) = 0.0
    endif
    if ( idproc == nvp-1 ) then
      this%buf(4*nrp-3) = -f1_re(2,nrp) - f3_re(1,nrp)
      this%buf(4*nrp-2) = 0.0
      this%buf(4*nrp-1) = f1_re(1,nrp) + idrh * ( 3.0 * f2_re(3,nrp) - 4.0 * f2_re(3,nrp-1) + f2_re(3,nrp-2) ) - f3_re(2,nrp)
      this%buf(4*nrp)   = 0.0
    endif

  elseif ( mode > 0 .and. present( jay_im ) .and. present( djdxi_im ) ) then
    
    do i = 1, nrp

      ir = idr / (real(i+noff)-0.5)
      this%buf(4*i-3) = -f1_re(2,i) + mode * f2_im(3,i) * ir - f3_re(1,i)
      this%buf(4*i-2) = -f1_im(2,i) - mode * f2_re(3,i) * ir - f3_im(1,i)
      this%buf(4*i-1) = f1_re(1,i) + idrh * ( f2_re(3,i+1)-f2_re(3,i-1) ) - f3_re(2,i)
      this%buf(4*i)   = f1_im(1,i) + idrh * ( f2_im(3,i+1)-f2_im(3,i-1) ) - f3_im(2,i)
      
    enddo

    ! calculate the derivatives at the boundary and axis
    if ( idproc == 0 ) then
      ir = 2.0 * idr
      this%buf(1) = -f1_re(2,1) + mode * f2_im(3,1) * ir - f3_re(1,1)
      this%buf(2) = -f1_im(2,1) - mode * f2_re(3,1) * ir - f3_im(1,1)
      this%buf(3) = f1_re(1,1) + idr * ( -3.0 * f2_re(3,1) + 4.0 * f2_re(3,2) - f2_re(3,3) ) - f3_re(2,1)
      this%buf(4) = f1_im(1,1) + idr * ( 3.0 * f2_im(3,1) - 4.0 * f2_im(3,2) + f2_im(3,3) ) - f3_im(2,1)
    endif
    if ( idproc == nvp-1 ) then
      ir = idr / (real(nrp+noff)-0.5)
      this%buf(4*nrp-3) = -f1_re(2,nrp) + mode * f2_im(3,nrp) * ir - f3_re(1,nrp)
      this%buf(4*nrp-2) = -f1_im(2,nrp) - mode * f2_re(3,nrp) * ir - f3_im(1,nrp)
      this%buf(4*nrp-1) = f1_re(1,nrp) + idr * ( -3.0 * f2_re(3,nrp) + 4.0 * f2_re(3,nrp-1) - f2_re(3,nrp-2) ) - f3_re(2,nrp)
      this%buf(4*nrp)   = f1_im(1,nrp) + idr * ( 3.0 * f2_im(3,nrp) - 4.0 * f2_im(3,nrp-1) + f2_im(3,nrp-2) ) - f3_im(2,nrp)
    endif

  else

    call write_err( 'Invalid input arguments!' )

  endif

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine set_source_bperp_iter

subroutine get_solution_bz( this, mode )

  implicit none

  class( field_b ), intent(inout) :: this
  integer, intent(in) :: mode

  integer :: i, nd1p
  real, dimension(:,:), pointer :: f1_re => null(), f1_im => null()
  character(len=20), save :: sname = 'get_solution_bz'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  nd1p = this%rf_re(mode)%get_ndp(1)

  f1_re => this%rf_re(mode)%get_f1()
  do i = 1, nd1p
    f1_re(3,i) = this%buf_re(i)
  enddo

  if ( mode > 0 ) then
    f1_im => this%rf_im(mode)%get_f1()
    do i = 1, nd1p
      f1_im(3,i) = this%buf_im(i)
    enddo
  endif

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine get_solution_bz

subroutine get_solution_bperp_old( this, mode )

  implicit none

  class( field_b ), intent(inout) :: this
  integer, intent(in) :: mode

  integer :: i, nd1p
  real, dimension(:,:), pointer :: f1_re => null(), f1_im => null()
  character(len=20), save :: sname = 'get_solution_bperp_old'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  nd1p = this%rf_re(mode)%get_ndp(1)

  f1_re => this%rf_re(mode)%get_f1()
  do i = 1, nd1p
    f1_re(1,i) = this%buf(4*i-3)
    f1_re(2,i) = this%buf(4*i-1)
  enddo

  if ( mode > 0 ) then
    f1_im => this%rf_im(mode)%get_f1()
    do i = 1, nd1p
      f1_im(1,i) = this%buf(4*i-2)
      f1_im(2,i) = this%buf(4*i)
    enddo
  endif

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine get_solution_bperp_old

subroutine get_solution_bperp( this, mode )

  implicit none

  class( field_b ), intent(inout) :: this
  integer, intent(in) :: mode

  integer :: i, nrp, idproc, nvp, noff, ierr, msgid1, msgid2
  real :: idr, idrh, ir, r2
  real, dimension(2), save :: lbuf, ubuf
  real, dimension(:,:), pointer :: f1_re => null(), f1_im => null()
  integer, dimension(MPI_STATUS_SIZE) :: stat
  character(len=20), save :: sname = 'get_solution_bperp'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  nrp    = this%rf_re(mode)%get_ndp(1)
  idproc = this%rf_re(mode)%pp%getlidproc()
  nvp    = this%rf_re(mode)%pp%getlnvp()
  noff   = this%rf_re(mode)%get_noff(1)
  idr    = 1.0 / this%dr
  idrh   = 0.5 * idr

  lbuf = 0.0; ubuf = 0.0

  ! copy the guard cell of buffer

  ! forward message passing
  ! receiver
  if ( idproc > 0 ) then
    call MPI_IRECV( lbuf(1), 1, dtype, idproc-1, 1, comm, msgid1, ierr )
    call MPI_IRECV( lbuf(2), 1, dtype, idproc-1, 2, comm, msgid2, ierr )
  endif
  ! sender
  if ( idproc < nvp-1 ) then
    call MPI_SEND( this%buf_re(nrp), 1, dtype, idproc+1, 1, comm, ierr )
    call MPI_SEND( this%buf_im(nrp), 1, dtype, idproc+1, 2, comm, ierr )
  endif
  ! wait receiving finish
  if ( idproc > 0 ) then
    call MPI_WAIT( msgid1, stat, ierr )
    call MPI_WAIT( msgid2, stat, ierr )
  endif


  ! backward message passing
  ! receiver
  if ( idproc < nvp-1 ) then
    call MPI_IRECV( ubuf(1), 1, dtype, idproc+1, 1, comm, msgid1, ierr )
    call MPI_IRECV( ubuf(2), 1, dtype, idproc+1, 2, comm, msgid2, ierr )
  endif
  ! sender
  if ( idproc > 0 ) then
    call MPI_SEND( this%buf_re(1), 1, dtype, idproc-1, 1, comm, ierr )
    call MPI_SEND( this%buf_im(1), 1, dtype, idproc-1, 2, comm, ierr )
  endif
  ! wait receiving finish
  if ( idproc < nvp-1 ) then
    call MPI_WAIT( msgid1, stat, ierr )
    call MPI_WAIT( msgid2, stat, ierr )
  endif


  if ( mode == 0 ) then

    f1_re => this%rf_re(mode)%get_f1()
    do i = 2, nrp-1
      ir = idr / ( real(i+noff) - 0.5 )
      f1_re(1,i) = 0.0
      f1_re(2,i) = -idrh * ( this%buf_re(i+1) - this%buf_re(i-1) )
    enddo
    
    if ( idproc == 0 ) then
      f1_re(1,1) = 0.0
      f1_re(2,1) = -idrh * ( -3.0 * this%buf_re(1) + 4.0 * this%buf_re(2) - this%buf_re(3) )
    else
      f1_re(1,1) = 0.0
      f1_re(2,1) = -idrh * ( this%buf_re(2) - lbuf(1) )
    endif

    if ( idproc == nvp-1 ) then
      f1_re(1,nrp) = 0.0
      f1_re(2,nrp) = -idrh * ( 3.0 * this%buf_re(nrp) - 4.0 * this%buf_re(nrp-1) + this%buf_re(nrp-2) )
    else
      f1_re(1,nrp) = 0.0
      f1_re(2,nrp) = -idrh * ( ubuf(1) - this%buf_re(nrp-1) )
    endif

  else

    f1_re => this%rf_re(mode)%get_f1()
    f1_im => this%rf_im(mode)%get_f1()
    do i = 2, nrp-1
      ir = idr / ( real(i+noff) - 0.5 )
      f1_re(1,i) = -ir * mode * this%buf_im(i)
      f1_re(2,i) = -idrh * ( this%buf_re(i+1) - this%buf_re(i-1) )
      f1_im(1,i) = ir * mode * this%buf_re(i)
      f1_im(2,i) = -idrh * ( this%buf_im(i+1) - this%buf_im(i-1) )
    enddo

    if ( idproc == 0 ) then
      ir = 2.0 * idr
      f1_re(1,1) = -ir * mode * this%buf_im(1)
      f1_re(2,1) = -idrh * ( -3.0 * this%buf_re(1) + 4.0 * this%buf_re(2) - this%buf_re(3) )
      f1_im(1,1) = ir * mode * this%buf_re(1)
      f1_im(2,1) = -idrh * ( -3.0 * this%buf_im(1) + 4.0 * this%buf_im(2) - this%buf_im(3) )
    else
      ir = idr / ( real(1+noff) - 0.5 )
      f1_re(1,1) = -ir * mode * this%buf_im(1)
      f1_re(2,1) = -idrh * ( this%buf_re(2) - lbuf(1) )
      f1_im(1,1) = ir * mode * this%buf_re(1)
      f1_im(2,1) = -idrh * ( this%buf_im(2) - lbuf(2) )
    endif

    if ( idproc == nvp-1 ) then
      ir = idr / ( real(nrp+noff) - 0.5 )
      f1_re(1,nrp) = -ir * mode * this%buf_im(nrp)
      f1_re(2,nrp) = -idrh * ( 3.0 * this%buf_re(nrp) - 4.0 * this%buf_re(nrp-1) + this%buf_re(nrp-2) )
      f1_im(1,nrp) = ir * mode * this%buf_re(nrp)
      f1_im(2,nrp) = -idrh * ( 3.0 * this%buf_im(nrp) - 4.0 * this%buf_im(nrp-1) + this%buf_im(nrp-2) )
    else
      ir = idr / ( real(nrp+noff) - 0.5 )
      f1_re(1,nrp) = -ir * mode * this%buf_im(nrp)
      f1_re(2,nrp) = -idrh * ( ubuf(1) - this%buf_re(nrp-1) )
      f1_im(1,nrp) = ir * mode * this%buf_re(nrp)
      f1_im(2,nrp) = -idrh * ( ubuf(2) - this%buf_im(nrp-1) )
    endif

  endif

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine get_solution_bperp

subroutine get_solution_bperp_iter( this, mode )
! this is totally the same as get_solution_bperp()

  implicit none

  class( field_b ), intent(inout) :: this
  integer, intent(in) :: mode

  integer :: i, nd1p
  real, dimension(:,:), pointer :: f1_re => null(), f1_im => null()
  character(len=20), save :: sname = 'get_solution_bperp_iter'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  nd1p = this%rf_re(mode)%get_ndp(1)

  f1_re => this%rf_re(mode)%get_f1()
  do i = 1, nd1p
    f1_re(1,i) = this%buf(4*i-3)
    f1_re(2,i) = this%buf(4*i-1)
  enddo

  if ( mode > 0 ) then
    f1_im => this%rf_im(mode)%get_f1()
    do i = 1, nd1p
      f1_im(1,i) = this%buf(4*i-2)
      f1_im(2,i) = this%buf(4*i)
    enddo
  endif

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine get_solution_bperp_iter

subroutine solve_field_bz( this, jay )

  implicit none

  class( field_b ), intent(inout) :: this
  class( field_jay ), intent(inout) :: jay

  type( ufield ), dimension(:), pointer :: jay_re => null(), jay_im => null()
  integer :: i
  character(len=20), save :: sname = 'solve_field_bz'

  call write_dbg( cls_name, sname, cls_level, 'starts' )
  call start_tprof( 'solve bz' )

  ! call jay%copy_gc_f1()

  jay_re => jay%get_rf_re()
  jay_im => jay%get_rf_im()

  do i = 0, this%num_modes

    if ( i == 0 ) then
      call this%set_source_bz( i, jay_re(i) )
      call this%solver_bz(i)%solve( this%buf_re )
      call this%get_solution_bz(i)
      cycle
    endif

    call this%set_source_bz( i, jay_re(i), jay_im(i) )
    call this%solver_bz(i)%solve( this%buf_re )
    call this%solver_bz(i)%solve( this%buf_im )
    call this%get_solution_bz(i)

  enddo

  call this%copy_gc_f1()

  call stop_tprof( 'solve bz' )
  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine solve_field_bz

subroutine solve_field_bperp_old( this, rho )

  implicit none

  class( field_b ), intent(inout) :: this
  class( field_rho ), intent(inout) :: rho

  type( ufield ), dimension(:), pointer :: rho_re => null(), rho_im => null()
  integer :: i
  character(len=20), save :: sname = 'solve_field_bperp_old'

  call write_dbg( cls_name, sname, cls_level, 'starts' )
  call start_tprof( 'solve beam bperp (old)' )

  ! call rho%copy_gc_f1()

  rho_re => rho%get_rf_re()
  rho_im => rho%get_rf_im()

  do i = 0, this%num_modes

    if ( i == 0 ) then
      call this%set_source_bperp_old( i, rho_re(i) )
      call this%solver_bperp(i)%solve( this%buf )
      call this%get_solution_bperp_old(i)
      cycle
    endif

    call this%set_source_bperp_old( i, rho_re(i), rho_im(i) )
    call this%solver_bperp(i)%solve( this%buf )
    call this%get_solution_bperp_old(i)

  enddo

  call this%copy_gc_f1()

  call stop_tprof( 'solve beam bperp (old)' )
  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine solve_field_bperp_old

subroutine solve_field_bperp( this, rho )

  implicit none

  class( field_b ), intent(inout) :: this
  class( field_rho ), intent(inout) :: rho

  type( ufield ), dimension(:), pointer :: rho_re => null(), rho_im => null()
  integer :: i
  character(len=20), save :: sname = 'solve_field_bperp'

  call write_dbg( cls_name, sname, cls_level, 'starts' )
  call start_tprof( 'solve beam bperp' )

  ! call rho%copy_gc_f1()

  rho_re => rho%get_rf_re()
  rho_im => rho%get_rf_im()

  do i = 0, this%num_modes

    if ( i == 0 ) then
      call this%set_source_bperp( i, rho_re(i) )
      call this%solver_bperp(i)%solve( this%buf_re )
      call this%get_solution_bperp(i)
      cycle
    endif

    call this%set_source_bperp( i, rho_re(i), rho_im(i) )
    call this%solver_bperp(i)%solve( this%buf_re )
    call this%solver_bperp(i)%solve( this%buf_im )
    call this%get_solution_bperp(i)

  enddo

  call this%copy_gc_f1()

  call stop_tprof( 'solve beam bperp' )
  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine solve_field_bperp

subroutine solve_field_bperp_iter( this, djdxi, jay )

  implicit none

  class( field_b ), intent(inout) :: this
  class( field_djdxi ), intent(in) :: djdxi
  class( field_jay ), intent(inout) :: jay

  type( ufield ), dimension(:), pointer :: jay_re => null(), jay_im => null()
  type( ufield ), dimension(:), pointer :: djdxi_re => null(), djdxi_im => null()
  integer :: i
  character(len=20), save :: sname = 'solve_field_bperp_iter'

  call write_dbg( cls_name, sname, cls_level, 'starts' )
  call start_tprof( 'solve plasma bperp' )

  ! call jay%copy_gc_f1()
  ! call djdxi%copy_gc_f1() ! no need for djdxi to copy guard cells

  djdxi_re => djdxi%get_rf_re()
  djdxi_im => djdxi%get_rf_im()
  jay_re => jay%get_rf_re()
  jay_im => jay%get_rf_im()

  do i = 0, this%num_modes

    if ( i == 0 ) then
      call this%set_source_bperp_iter( i, djdxi_re(i), jay_re(i) )
      call this%solver_bperp_iter(i)%solve( this%buf )
      call this%get_solution_bperp_iter(i)
      cycle
    endif

    call this%set_source_bperp_iter( i, djdxi_re(i), jay_re(i), djdxi_im(i), jay_im(i) )
    call this%solver_bperp_iter(i)%solve( this%buf )
    call this%get_solution_bperp_iter(i)

  enddo

  call this%copy_gc_f1()

  call stop_tprof( 'solve plasma bperp' )
  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine solve_field_bperp_iter

end module field_b_class
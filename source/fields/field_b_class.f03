module field_b_class

use field_class
use field_src_class
use field_solver_class
use ufield_class
use param
use sys
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

  class( field_solver ), dimension(:), pointer :: solver_bz      => null()
  class( field_solver ), dimension(:), pointer :: solver_bt      => null()
  class( field_solver ), dimension(:), pointer :: solver_bt_iter => null()
  class( field_solver ), dimension(:), pointer :: solver_bplus   => null()
  class( field_solver ), dimension(:), pointer :: solver_bminus  => null()

  real, dimension(:), pointer :: buf1_re => null(), buf1_im => null()
  real, dimension(:), pointer :: buf2_re => null(), buf2_im => null()
  real, dimension(:), pointer :: buf    => null()

  contains

  generic :: new        => init_field_b
  procedure :: del      => end_field_b
  generic :: solve      => solve_field_bz, solve_field_bt, solve_field_bt_iter
  ! generic :: solve_old  => solve_field_bt_old, solve_field_bt_iter_old

  procedure, private :: init_field_b
  procedure, private :: end_field_b
  procedure, private :: set_source_bz
  procedure, private :: set_source_bt
  ! procedure, private :: set_source_bt_old
  procedure, private :: set_source_bt_iter
  ! procedure, private :: set_source_bt_iter_old
  procedure, private :: get_solution_bz
  procedure, private :: get_solution_bt
  ! procedure, private :: get_solution_bt_old
  procedure, private :: get_solution_bt_iter
  ! procedure, private :: get_solution_bt_iter_old
  procedure, private :: solve_field_bz
  procedure, private :: solve_field_bt
  ! procedure, private :: solve_field_bt_old
  procedure, private :: solve_field_bt_iter
  ! procedure, private :: solve_field_bt_iter_old

end type field_b

contains

subroutine init_field_b( this, pp, gp, num_modes, part_shape, boundary, entity, iter_tol )

  implicit none

  class( field_b ), intent(inout) :: this
  class( parallel_pipe ), intent(in), pointer :: pp
  class( grid ), intent(in), pointer :: gp
  integer, intent(in) :: num_modes, part_shape, entity, boundary
  real, intent(in) :: iter_tol

  integer, dimension(2,2) :: gc_num
  integer :: dim, i, nrp
  real :: dr
  character(len=20), save :: sname = "init_field_b"

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  nrp = gp%get_ndp(1)
  dr  = gp%get_dr()

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
  call this%field%new( pp, gp, dim, num_modes, gc_num, entity )

  ! initialize solver
  select case ( entity )

  ! case ( p_entity_plasma_old )

  !   allocate( this%solver_bz( 0:num_modes ) )
  !   allocate( this%solver_bt_iter( 0:num_modes ) )
  !   do i = 0, num_modes
  !     call this%solver_bz(i)%new( pp, gp, i, dr, kind=p_fk_bz, &
  !       bnd=boundary, stype=p_hypre_cycred, tol=iter_tol )
  !     call this%solver_bt_iter(i)%new( pp, gp, i, dr, kind=p_fk_bt_iter, &
  !       bnd=boundary, stype=p_hypre_amg, tol=iter_tol )
  !   enddo

  !   allocate( this%buf(4*nrp) )
  !   allocate( this%buf1_re(nrp), this%buf1_im(nrp) )

  case ( p_entity_plasma )

    allocate( this%solver_bz( 0:num_modes ) )
    allocate( this%solver_bplus( 0:num_modes ) )
    allocate( this%solver_bminus( 0:num_modes ) )
    do i = 0, num_modes
      call this%solver_bz(i)%new( pp, gp, i, dr, kind=p_fk_bz, &
        bnd=boundary, stype=p_hypre_cycred, tol=iter_tol )
      call this%solver_bplus(i)%new( pp, gp, i, dr, kind=p_fk_bplus, &
        bnd=boundary, stype=p_hypre_cycred, tol=iter_tol )
      call this%solver_bminus(i)%new( pp, gp, i, dr, kind=p_fk_bminus, &
        bnd=boundary, stype=p_hypre_cycred, tol=iter_tol )
    enddo

    allocate( this%buf1_re(nrp), this%buf1_im(nrp) )
    allocate( this%buf2_re(nrp), this%buf2_im(nrp) )

  case ( p_entity_beam_old )

    allocate( this%solver_bt( 0:num_modes ) )
    do i = 0, num_modes
      call this%solver_bt(i)%new( pp, gp, i, dr, kind=p_fk_bt_old, &
        bnd=boundary, stype=p_hypre_amg, tol=iter_tol )
    enddo

    allocate( this%buf( nrp*4 ) )

  case ( p_entity_beam )

    allocate( this%solver_bt( 0:num_modes ) )
    do i = 0, num_modes
      call this%solver_bt(i)%new( pp, gp, i, dr, kind=p_fk_bt, &
        bnd=boundary, stype=p_hypre_cycred, tol=iter_tol )
    enddo

    allocate( this%buf1_re(nrp), this%buf1_im(nrp) )

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
  ! case ( p_entity_plasma_old )
  !   do i = 0, this%num_modes
  !     call this%solver_bz(i)%del()
  !     call this%solver_bt_iter(i)%del()
  !   enddo
  !   deallocate( this%solver_bz )
  !   deallocate( this%solver_bt_iter )
  case ( p_entity_plasma )
    do i = 0, this%num_modes
      call this%solver_bz(i)%del()
      call this%solver_bplus(i)%del()
      call this%solver_bminus(i)%del()
    enddo
    deallocate( this%solver_bz )
    deallocate( this%solver_bplus )
    deallocate( this%solver_bminus )
  case ( p_entity_beam, p_entity_beam_old )
    do i = 0, this%num_modes
      call this%solver_bt(i)%del()
    enddo
    deallocate( this%solver_bt )
  end select

  if ( associated( this%buf1_re ) ) deallocate( this%buf1_re )
  if ( associated( this%buf1_im ) ) deallocate( this%buf1_im )
  if ( associated( this%buf2_re ) ) deallocate( this%buf2_re )
  if ( associated( this%buf2_im ) ) deallocate( this%buf2_im )
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
  real :: idrh, idr, dr, dr2, ir
  character(len=32) :: filename
  character(len=20), save :: sname = 'set_source_bz'

  call write_dbg( cls_name, sname, cls_level, 'starts' )
  call start_tprof( 'solve bz' )

  nrp    = jay_re%get_ndp(1)
  idr    = 1.0 / this%dr
  idrh   = 0.5 * idr
  noff   = jay_re%get_noff(1)
  nvp    = jay_re%pp%getlnvp()
  idproc = jay_re%pp%getlidproc()
  dr     = this%dr
  dr2    = dr*dr

  f1_re => jay_re%get_f1()
  if ( present(jay_im) ) then
    f1_im => jay_im%get_f1()
  endif

  this%buf1_re = 0.0
  if ( present(jay_im) ) this%buf1_im = 0.0

  if ( mode == 0 ) then

    do i = 2, nrp
      ir = idr / real(i+noff-1)
      this%buf1_re(i) = -idrh * ( f1_re(2,i+1) - f1_re(2,i-1) ) - ir * f1_re(2,i)
    enddo

    ! calculate the derivatives at the boundary and axis
    if ( idproc == 0 ) then
      this%buf1_re(1) = -2.0 * idr * f1_re(2,2)
    else
      ir = idr / real(noff)
      this%buf1_re(1) = -idrh * ( f1_re(2,2) - f1_re(2,0) ) - ir * f1_re(2,1)
    endif

    if ( idproc == nvp-1 ) then
      ir = idr / real(nrp+noff-1)
      this%buf1_re(nrp) = -idrh * ( 3.0 * f1_re(2,nrp) - 4.0 * f1_re(2,nrp-1) + f1_re(2,nrp-2) ) - ir * f1_re(2,nrp)
    endif

  elseif ( mode > 0 .and. present( jay_im ) ) then

    do i = 2, nrp
      ir = idr / real(i+noff-1)
      this%buf1_re(i) = -idrh * ( f1_re(2,i+1) - f1_re(2,i-1) ) - ir * f1_re(2,i) - mode * ir * f1_im(1,i)
      this%buf1_im(i) = -idrh * ( f1_im(2,i+1) - f1_im(2,i-1) ) - ir * f1_im(2,i) + mode * ir * f1_re(1,i)
    enddo

    ! calculate the derivatives at the boundary and axis
    if ( idproc == 0 ) then
      this%buf1_re(1) = -2.0 * idr * f1_re(2,2) - mode * idr * f1_im(1,2)
      this%buf1_im(1) = -2.0 * idr * f1_im(2,2) + mode * idr * f1_re(1,2)
    else
      ir = idr / real(noff)
      this%buf1_re(1) = -idrh * ( f1_re(2,2) - f1_re(2,0) ) - ir * f1_re(2,1) - mode * ir * f1_im(1,1)
      this%buf1_im(1) = -idrh * ( f1_im(2,2) - f1_im(2,0) ) - ir * f1_im(2,1) + mode * ir * f1_re(1,1)
    endif

    if ( idproc == nvp-1 ) then
      ir = idr / real(nrp+noff-1)
      this%buf1_re(nrp) = -idrh * ( 3.0 * f1_re(2,nrp) - 4.0 * f1_re(2,nrp-1) + f1_re(2,nrp-2) ) - ir * f1_re(2,nrp) &
                         -mode * ir * f1_im(1,nrp)
      this%buf1_im(nrp) = -idrh * ( 3.0 * f1_im(2,nrp) - 4.0 * f1_im(2,nrp-1) + f1_im(2,nrp-2) ) - ir * f1_im(2,nrp) &
                         +mode * ir * f1_re(1,nrp)
    endif

  else

    call write_err( 'Invalid input arguments!' )

  endif

  call stop_tprof( 'solve bz' )
  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine set_source_bz

! subroutine set_source_bt_old( this, mode, q_re, q_im )

!   implicit none

!   class( field_b ), intent(inout) :: this
!   class( ufield ), intent(in) :: q_re
!   class( ufield ), intent(in), optional :: q_im
!   integer, intent(in) :: mode

!   integer :: i, nrp, idproc, nvp, noff, comm, dtype, ierr
!   real, dimension(:,:), pointer :: f1_re => null(), f1_im => null()
!   real :: idrh, idr, a1, a2, a3, b, ir, dr2, rmax
!   character(len=20), save :: sname = 'set_source_bt_old'

!   call write_dbg( cls_name, sname, cls_level, 'starts' )
!   call start_tprof( 'solve beam bt' )

!   idproc = q_re%pp%getlidproc()
!   nvp = q_re%pp%getlnvp()
!   nrp = q_re%get_ndp(1)
!   noff = q_re%get_noff(1)
!   idr = 1.0 / this%dr
!   idrh = 0.5 * idr
!   dtype = q_re%pp%getmreal()
!   comm = q_re%pp%getlgrp()
!   dr2 = (this%dr)**2
!   rmax = (q_re%get_nd(1)-0.5) * this%dr

!   f1_re => q_re%get_f1()
!   if ( present(q_im) ) then
!     f1_im => q_im%get_f1()
!   endif

!   this%buf = 0.0
!   if ( mode == 0 ) then

!     select case ( this%solver_bt(0)%bnd )

!     case ( p_bnd_zero, p_bnd_open )

!       do i = 1, nrp

!         ! Re(Br)
!         this%buf(4*i-3) = 0.0
!         ! Im(Br)
!         this%buf(4*i-2) = 0.0
!         ! Re(Bphi)
!         this%buf(4*i-1) = idrh * ( f1_re(1,i+1)-f1_re(1,i-1) )
!         ! Im(Bphi)
!         this%buf(4*i) = 0.0

!       enddo

!     end select

!     ! calculate the derivatives at the boundary and axis
!     if ( idproc == 0 ) then
!       this%buf(3) = idrh * ( -3.0 * f1_re(1,1) + 4.0 * f1_re(1,2) - f1_re(1,3) )
!     endif
!     if ( idproc == nvp-1 ) then
!       this%buf(4*nrp-1) = idrh * ( 3.0 * f1_re(1,nrp) - 4.0 * f1_re(1,nrp-1) + f1_re(1,nrp-2) )
!       select case ( this%solver_bt(0)%bnd )
!       case ( p_bnd_zero, p_bnd_open )
!         ! do nothing
!       end select
!     endif

!   elseif ( mode > 0 .and. present( q_im ) ) then

!     do i = 1, nrp

!       ir = idr / (real(i+noff)-0.5)
!       this%buf(4*i-3) =  mode * f1_im(1,i) * ir
!       this%buf(4*i-2) = -mode * f1_re(1,i) * ir
!       this%buf(4*i-1) = idrh * ( f1_re(1,i+1)-f1_re(1,i-1) )
!       this%buf(4*i)   = idrh * ( f1_im(1,i+1)-f1_im(1,i-1) )

!     enddo

!     ! calculate the derivatives at the boundary and axis
!     if ( idproc == 0 ) then
!       this%buf(3) = idrh * ( -3.0 * f1_re(1,1) + 4.0 * f1_re(1,2) - f1_re(1,3) )
!       this%buf(4) = idrh * ( -3.0 * f1_im(1,1) + 4.0 * f1_im(1,2) - f1_im(1,3) )
!     endif
!     if ( idproc == nvp-1 ) then
!       this%buf(4*nrp-1) = idrh * ( 3.0 * f1_re(1,nrp) - 4.0 * f1_re(1,nrp-1) + f1_re(1,nrp-2) )
!       this%buf(4*nrp)   = idrh * ( 3.0 * f1_im(1,nrp) - 4.0 * f1_im(1,nrp-1) + f1_im(1,nrp-2) )
!     endif

!   else

!     call write_err( 'Invalid input arguments!' )

!   endif

!   call stop_tprof( 'solve beam bt' )
!   call write_dbg( cls_name, sname, cls_level, 'ends' )

! end subroutine set_source_bt_old

subroutine set_source_bt( this, mode, q_re, q_im )

  implicit none

  class( field_b ), intent(inout) :: this
  class( ufield ), intent(in) :: q_re
  class( ufield ), intent(in), optional :: q_im
  integer, intent(in) :: mode

  integer :: i, nrp, noff, dtype, ierr, comm
  real, dimension(:,:), pointer :: f1_re => null(), f1_im => null()
  real ::dr, dr2, rmax
  character(len=20), save :: sname = 'set_source_bt'

  call write_dbg( cls_name, sname, cls_level, 'starts' )
  call start_tprof( 'solve beam bt' )

  nrp   = q_re%get_ndp(1)
  noff  = q_re%get_noff(1)
  dr    = this%dr
  dr2   = dr**2
  dtype = q_re%pp%getmreal()
  comm  = q_re%pp%getlgrp()
  rmax  = (q_re%get_nd(1)-0.5) * dr

  f1_re => q_re%get_f1()
  this%buf1_re = 0.0
  if ( present(q_im) ) then
    f1_im => q_im%get_f1()
    this%buf1_im = 0.0
  endif

  if ( mode == 0 ) then

    do i = 1, nrp
      this%buf1_re(i) = -1.0 * f1_re(1,i)
    enddo

  elseif ( mode > 0 .and. present(q_im) ) then
    do i = 1, nrp
      this%buf1_re(i) = -1.0 * f1_re(1,i)
      this%buf1_im(i) = -1.0 * f1_im(1,i)
    enddo
  else
    call write_err( 'Invalid input arguments!' )
  endif

  call stop_tprof( 'solve beam bt' )
  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine set_source_bt

! subroutine set_source_bt_iter_old( this, mode, djdxi_re, jay_re, djdxi_im, jay_im )

!   implicit none

!   class( field_b ), intent(inout) :: this
!   class( ufield ), intent(in) :: djdxi_re, jay_re
!   class( ufield ), intent(in), optional :: djdxi_im, jay_im
!   integer, intent(in) :: mode

!   integer :: i, nrp, nvp, idproc, noff, dtype, ierr, comm
!   real, dimension(:,:), pointer :: f1_re => null(), f1_im => null()
!   real, dimension(:,:), pointer :: f2_re => null(), f2_im => null()
!   real, dimension(:,:), pointer :: f3_re => null(), f3_im => null()
!   real :: idrh, idr, a1, a2, a3, b, ir, dr2, dr, rmax
!   character(len=20), save :: sname = 'set_source_bt_iter_old'

!   call write_dbg( cls_name, sname, cls_level, 'starts' )
!   call start_tprof( 'solve plasma bt' )

!   dtype  = jay_re%pp%getmreal()
!   comm   = jay_re%pp%getlgrp()
!   nvp    = jay_re%pp%getlnvp()
!   idproc = jay_re%pp%getlidproc()
!   nrp    = jay_re%get_ndp(1)
!   noff   = jay_re%get_noff(1)
!   idr    = 1.0 / this%dr
!   idrh   = 0.5 * idr
!   dr     = this%dr
!   dr2    = dr**2
!   rmax   = (jay_re%get_nd(1)-0.5) * dr

!   f1_re => djdxi_re%get_f1()
!   f2_re => jay_re%get_f1()
!   f3_re => this%rf_re(mode)%get_f1()

!   if ( present(djdxi_im) .and. present(jay_im) ) then
!     f1_im => djdxi_im%get_f1()
!     f2_im => jay_im%get_f1()
!     f3_im => this%rf_im(mode)%get_f1()
!   endif

!   this%buf = 0.0
!   if ( mode == 0 ) then

!     select case ( this%solver_bt_iter(0)%bnd )

!     case ( p_bnd_zero, p_bnd_open )

!       do i = 1, nrp

!         ! Re(Br)
!         this%buf(4*i-3) = -f1_re(2,i) - f3_re(1,i)
!         ! Im(Br)
!         this%buf(4*i-2) = 0.0
!         ! Re(Bphi)
!         this%buf(4*i-1) = f1_re(1,i) + idrh * ( f2_re(3,i+1)-f2_re(3,i-1) ) - f3_re(2,i)
!         ! Im(Bphi)
!         this%buf(4*i) = 0.0

!       enddo

!     end select

!     ! calculate the derivatives at the boundary and axis
!     if ( idproc == 0 ) then
!       this%buf(1) = -f1_re(2,1) - f3_re(1,1)
!       this%buf(2) = 0.0
!       this%buf(3) = f1_re(1,1) + idrh * ( -3.0 * f2_re(3,1) + 4.0 * f2_re(3,2) - f2_re(3,3) ) - f3_re(2,1)
!       this%buf(4) = 0.0
!     endif
!     if ( idproc == nvp-1 ) then
!       this%buf(4*nrp-3) = -f1_re(2,nrp) - f3_re(1,nrp)
!       this%buf(4*nrp-2) = 0.0
!       this%buf(4*nrp-1) = f1_re(1,nrp) + idrh * ( 3.0 * f2_re(3,nrp) - 4.0 * f2_re(3,nrp-1) + f2_re(3,nrp-2) ) - f3_re(2,nrp)
!       this%buf(4*nrp)   = 0.0

!       select case ( this%solver_bt_iter(0)%bnd )
!       case ( p_bnd_zero, p_bnd_open )
!         ! do nothing
!       end select
!     endif

!   elseif ( mode > 0 .and. present( jay_im ) .and. present( djdxi_im ) ) then

!     do i = 1, nrp

!       ir = idr / (real(i+noff)-0.5)
!       this%buf(4*i-3) = -f1_re(2,i) + mode * f2_im(3,i) * ir - f3_re(1,i)
!       this%buf(4*i-2) = -f1_im(2,i) - mode * f2_re(3,i) * ir - f3_im(1,i)
!       this%buf(4*i-1) = f1_re(1,i) + idrh * ( f2_re(3,i+1)-f2_re(3,i-1) ) - f3_re(2,i)
!       this%buf(4*i)   = f1_im(1,i) + idrh * ( f2_im(3,i+1)-f2_im(3,i-1) ) - f3_im(2,i)

!     enddo

!     ! calculate the derivatives at the boundary and axis
!     if ( idproc == 0 ) then
!       ir = 2.0 * idr
!       this%buf(1) = -f1_re(2,1) + mode * f2_im(3,1) * ir - f3_re(1,1)
!       this%buf(2) = -f1_im(2,1) - mode * f2_re(3,1) * ir - f3_im(1,1)
!       this%buf(3) = f1_re(1,1) - idrh * ( 3.0 * f2_re(3,1) - 4.0 * f2_re(3,2) + f2_re(3,3) ) - f3_re(2,1)
!       this%buf(4) = f1_im(1,1) - idrh * ( 3.0 * f2_im(3,1) - 4.0 * f2_im(3,2) + f2_im(3,3) ) - f3_im(2,1)
!     endif
!     if ( idproc == nvp-1 ) then
!       ir = idr / (real(nrp+noff)-0.5)
!       this%buf(4*nrp-3) = -f1_re(2,nrp) + mode * f2_im(3,nrp) * ir - f3_re(1,nrp)
!       this%buf(4*nrp-2) = -f1_im(2,nrp) - mode * f2_re(3,nrp) * ir - f3_im(1,nrp)
!       this%buf(4*nrp-1) = f1_re(1,nrp) + idrh * ( 3.0 * f2_re(3,nrp) - 4.0 * f2_re(3,nrp-1) + f2_re(3,nrp-2) ) - f3_re(2,nrp)
!       this%buf(4*nrp)   = f1_im(1,nrp) + idrh * ( 3.0 * f2_im(3,nrp) - 4.0 * f2_im(3,nrp-1) + f2_im(3,nrp-2) ) - f3_im(2,nrp)
!     endif

!   else

!     call write_err( 'Invalid input arguments!' )

!   endif

!   call stop_tprof( 'solve plasma bt' )
!   call write_dbg( cls_name, sname, cls_level, 'ends' )

! end subroutine set_source_bt_iter_old

subroutine set_source_bt_iter( this, mode, djdxi_re, jay_re, djdxi_im, jay_im )

  implicit none

  class( field_b ), intent(inout) :: this
  class( ufield ), intent(in) :: djdxi_re, jay_re
  class( ufield ), intent(in), optional :: djdxi_im, jay_im
  integer, intent(in) :: mode

  integer :: i, nrp, nvp, idproc, noff, ierr
  real, dimension(:,:), pointer :: f1_re => null(), f1_im => null()
  real, dimension(:,:), pointer :: f2_re => null(), f2_im => null()
  real, dimension(:,:), pointer :: f3_re => null(), f3_im => null()
  real :: idrh, idr, ir
  real :: s1_re, s1_im, s2_re, s2_im
  character(len=20), save :: sname = 'set_source_bt_iter'

  call write_dbg( cls_name, sname, cls_level, 'starts' )
  call start_tprof( 'solve plasma bt' )

  nvp    = jay_re%pp%getlnvp()
  idproc = jay_re%pp%getlidproc()
  nrp    = jay_re%get_ndp(1)
  noff   = jay_re%get_noff(1)
  idr    = 1.0 / this%dr
  idrh   = 0.5 * idr

  f1_re => djdxi_re%get_f1()
  f2_re => jay_re%get_f1()
  f3_re => this%rf_re(mode)%get_f1()
  this%buf1_re = 0.0
  this%buf2_re = 0.0

  if ( present(djdxi_im) .and. present(jay_im) ) then
    f1_im => djdxi_im%get_f1()
    f2_im => jay_im%get_f1()
    f3_im => this%rf_im(mode)%get_f1()
    this%buf1_im = 0.0
    this%buf2_im = 0.0
  endif

  if ( mode == 0 ) then

    do i = 2, nrp
      this%buf1_re(i) = -f1_re(2,i) - f3_re(1,i) ! Re(Br)
      this%buf2_re(i) = f1_re(1,i) + idrh * ( f2_re(3,i+1) - f2_re(3,i-1) ) - f3_re(2,i) ! Re(Bphi)
    enddo

    ! calculate the derivatives at the boundary and axis
    if ( idproc == 0 ) then
      this%buf1_re(1) = -f1_re(2,1) - f3_re(1,1)
      this%buf2_re(1) =  f1_re(1,1) - f3_re(2,1)
    else
      this%buf1_re(1) = -f1_re(2,1) - f3_re(1,1)
      this%buf2_re(1) =  f1_re(1,1) + idrh * ( f2_re(3,2) - f2_re(3,0) ) - f3_re(2,1)
    endif

    if ( idproc == nvp-1 ) then
      this%buf1_re(nrp) = -f1_re(2,nrp) - f3_re(1,nrp)
      this%buf2_re(nrp) =  f1_re(1,nrp) + idrh * ( 3.0 * f2_re(3,nrp) - 4.0 * f2_re(3,nrp-1) + f2_re(3,nrp-2) ) - f3_re(2,nrp)
    endif

  elseif ( mode > 0 .and. present( jay_im ) .and. present( djdxi_im ) ) then

    do i = 2, nrp
      ir = idr / real(i+noff-1)
      s1_re = -f1_re(2,i) + mode * f2_im(3,i) * ir
      s1_im = -f1_im(2,i) - mode * f2_re(3,i) * ir
      s2_re =  f1_re(1,i) + idrh * ( f2_re(3,i+1) - f2_re(3,i-1) )
      s2_im =  f1_im(1,i) + idrh * ( f2_im(3,i+1) - f2_im(3,i-1) )
      this%buf1_re(i) = s1_re - s2_im - f3_re(1,i) + f3_im(2,i) ! Re(B_plus)
      this%buf1_im(i) = s1_im + s2_re - f3_im(1,i) - f3_re(2,i) ! Im(B_plus)
      this%buf2_re(i) = s1_re + s2_im - f3_re(1,i) - f3_im(2,i) ! Re(B_minus)
      this%buf2_im(i) = s1_im - s2_re - f3_im(1,i) + f3_re(2,i) ! Im(B_minus)
    enddo

    ! calculate the derivatives at the boundary and axis
    if ( idproc == 0 ) then

      if ( mod(mode,2) == 0 ) then
        s1_re = -f1_re(2,1)
        s1_im = -f1_im(2,1)
        s2_re =  f1_re(1,1)
        s2_im =  f1_im(1,1)
      else
        s1_re = -f1_re(2,1) + idr * mode * f2_im(3,2)
        s1_im = -f1_im(2,1) - idr * mode * f2_re(3,2)
        s2_re =  f1_re(1,1) + idr * f2_re(3,2)
        s2_im =  f1_im(1,1) + idr * f2_im(3,2)
      endif

      this%buf1_re(1) = s1_re - s2_im - f3_re(1,1) + f3_im(2,1) ! Re(B_plus)
      this%buf1_im(1) = s1_im + s2_re - f3_im(1,1) - f3_re(2,1) ! Im(B_plus)
      this%buf2_re(1) = s1_re + s2_im - f3_re(1,1) - f3_im(2,1) ! Re(B_minus)
      this%buf2_im(1) = s1_im - s2_re - f3_im(1,1) + f3_re(2,1) ! Im(B_minus)

    else

      ir = idr / real(noff)
      s1_re = -f1_re(2,1) + mode * f2_im(3,1) * ir
      s1_im = -f1_im(2,1) - mode * f2_re(3,1) * ir
      s2_re =  f1_re(1,1) + idrh * ( f2_re(3,2) - f2_re(3,0) )
      s2_im =  f1_im(1,1) + idrh * ( f2_im(3,2) - f2_im(3,0) )

      this%buf1_re(1) = s1_re - s2_im - f3_re(1,1) + f3_im(2,1) ! Re(B_plus)
      this%buf1_im(1) = s1_im + s2_re - f3_im(1,1) - f3_re(2,1) ! Im(B_plus)
      this%buf2_re(1) = s1_re + s2_im - f3_re(1,1) - f3_im(2,1) ! Re(B_minus)
      this%buf2_im(1) = s1_im - s2_re - f3_im(1,1) + f3_re(2,1) ! Im(B_minus)

    endif

    if ( idproc == nvp-1 ) then

      ir = idr / real(nrp+noff-1)
      s1_re = -f1_re(2,nrp) + mode * f2_im(3,nrp) * ir
      s1_im = -f1_im(2,nrp) - mode * f2_re(3,nrp) * ir
      s2_re =  f1_re(1,nrp) + idrh * ( 3.0 * f2_re(3,nrp) - 4.0 * f2_re(3,nrp-1) + f2_re(3,nrp-2) )
      s2_im =  f1_im(1,nrp) + idrh * ( 3.0 * f2_im(3,nrp) - 4.0 * f2_im(3,nrp-1) + f2_im(3,nrp-2) )
      this%buf1_re(nrp) = s1_re - s2_im - f3_re(1,nrp) + f3_im(2,nrp) ! Re(B_plus)
      this%buf1_im(nrp) = s1_im + s2_re - f3_im(1,nrp) - f3_re(2,nrp) ! Im(B_plus)
      this%buf2_re(nrp) = s1_re + s2_im - f3_re(1,nrp) - f3_im(2,nrp) ! Re(B_minus)
      this%buf2_im(nrp) = s1_im - s2_re - f3_im(1,nrp) + f3_re(2,nrp) ! Im(B_minus)

    endif

  else

    call write_err( 'Invalid input arguments!' )

  endif

  call stop_tprof( 'solve plasma bt' )
  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine set_source_bt_iter

subroutine get_solution_bz( this, mode )

  implicit none

  class( field_b ), intent(inout) :: this
  integer, intent(in) :: mode

  integer :: i, nd1p
  real, dimension(:,:), pointer :: f1_re => null(), f1_im => null()
  character(len=20), save :: sname = 'get_solution_bz'

  call write_dbg( cls_name, sname, cls_level, 'starts' )
  call start_tprof( 'solve bz' )

  nd1p = this%rf_re(mode)%get_ndp(1)

  f1_re => this%rf_re(mode)%get_f1()
  do i = 1, nd1p
    f1_re(3,i) = this%buf1_re(i)
  enddo

  if ( mode > 0 ) then
    f1_im => this%rf_im(mode)%get_f1()
    do i = 1, nd1p
      f1_im(3,i) = this%buf1_im(i)
    enddo

    ! Bz(m>0) vanishes on axis
    f1_re(1,1) = 0.0
    f1_im(1,1) = 0.0
  endif

  call stop_tprof( 'solve bz' )
  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine get_solution_bz

! subroutine get_solution_bt_old( this, mode )

!   implicit none

!   class( field_b ), intent(inout) :: this
!   integer, intent(in) :: mode

!   integer :: i, nd1p
!   real, dimension(:,:), pointer :: f1_re => null(), f1_im => null()
!   character(len=20), save :: sname = 'get_solution_bt_old'

!   call write_dbg( cls_name, sname, cls_level, 'starts' )
!   call start_tprof( 'solve beam bt' )

!   nd1p = this%rf_re(mode)%get_ndp(1)

!   f1_re => this%rf_re(mode)%get_f1()
!   do i = 1, nd1p
!     f1_re(1,i) = this%buf(4*i-3)
!     f1_re(2,i) = this%buf(4*i-1)
!   enddo

!   if ( mode > 0 ) then
!     f1_im => this%rf_im(mode)%get_f1()
!     do i = 1, nd1p
!       f1_im(1,i) = this%buf(4*i-2)
!       f1_im(2,i) = this%buf(4*i)
!     enddo
!   endif

!   call stop_tprof( 'solve beam bt' )
!   call write_dbg( cls_name, sname, cls_level, 'ends' )

! end subroutine get_solution_bt_old

subroutine get_solution_bt( this, mode )

  implicit none

  class( field_b ), intent(inout) :: this
  integer, intent(in) :: mode

  integer :: i, nrp, idproc, nvp, noff, dtype, ierr, comm, msgid1, msgid2
  real :: idr, idrh, ir
  real, dimension(2), save :: lbuf, ubuf
  real, dimension(:,:), pointer :: f1_re => null(), f1_im => null()
  integer, dimension(MPI_STATUS_SIZE) :: stat
  character(len=20), save :: sname = 'get_solution_bt'

  call write_dbg( cls_name, sname, cls_level, 'starts' )
  call start_tprof( 'solve beam bt' )

  nrp    = this%rf_re(mode)%get_ndp(1)
  idproc = this%rf_re(mode)%pp%getlidproc()
  nvp    = this%rf_re(mode)%pp%getlnvp()
  noff   = this%rf_re(mode)%get_noff(1)
  comm   = this%rf_re(mode)%pp%getlgrp()
  dtype  = this%rf_re(mode)%pp%getmreal()
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
    call MPI_SEND( this%buf1_re(nrp), 1, dtype, idproc+1, 1, comm, ierr )
    call MPI_SEND( this%buf1_im(nrp), 1, dtype, idproc+1, 2, comm, ierr )
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
    call MPI_SEND( this%buf1_re(1), 1, dtype, idproc-1, 1, comm, ierr )
    call MPI_SEND( this%buf1_im(1), 1, dtype, idproc-1, 2, comm, ierr )
  endif
  ! wait receiving finish
  if ( idproc < nvp-1 ) then
    call MPI_WAIT( msgid1, stat, ierr )
    call MPI_WAIT( msgid2, stat, ierr )
  endif


  if ( mode == 0 ) then

    f1_re => this%rf_re(mode)%get_f1()
    do i = 2, nrp-1
      ir = idr / real(i+noff-1)
      f1_re(1,i) = 0.0
      f1_re(2,i) = -idrh * ( this%buf1_re(i+1) - this%buf1_re(i-1) )
    enddo

    if ( idproc == 0 ) then
      f1_re(1,1) = 0.0
      f1_re(2,1) = 0.0
    else
      f1_re(1,1) = 0.0
      f1_re(2,1) = -idrh * ( this%buf1_re(2) - lbuf(1) )
    endif

    if ( idproc == nvp-1 ) then
      f1_re(1,nrp) = 0.0
      f1_re(2,nrp) = -idrh * ( 3.0 * this%buf1_re(nrp) - 4.0 * this%buf1_re(nrp-1) + this%buf1_re(nrp-2) )
    else
      f1_re(1,nrp) = 0.0
      f1_re(2,nrp) = -idrh * ( ubuf(1) - this%buf1_re(nrp-1) )
    endif

  else

    f1_re => this%rf_re(mode)%get_f1()
    f1_im => this%rf_im(mode)%get_f1()
    do i = 2, nrp-1
      ir = idr / real(i+noff-1)
      f1_re(1,i) = -ir * mode * this%buf1_im(i)
      f1_im(1,i) =  ir * mode * this%buf1_re(i)
      f1_re(2,i) = -idrh * ( this%buf1_re(i+1) - this%buf1_re(i-1) )
      f1_im(2,i) = -idrh * ( this%buf1_im(i+1) - this%buf1_im(i-1) )
    enddo

    if ( idproc == 0 ) then

      if ( mod(mode,2) == 0 ) then
        f1_re(1,1) = 0.0
        f1_im(1,1) = 0.0
        f1_re(2,1) = 0.0
        f1_im(2,1) = 0.0
      else
        f1_re(1,1) = -idr * mode * this%buf1_im(2)
        f1_im(1,1) =  idr * mode * this%buf1_re(2)
        f1_re(2,1) = -idr * this%buf1_re(2)
        f1_im(2,1) = -idr * this%buf1_im(2)
      endif

    else

      ir = idr / real(noff)
      f1_re(1,1) = -ir * mode * this%buf1_im(1)
      f1_im(1,1) =  ir * mode * this%buf1_re(1)
      f1_re(2,1) = -idrh * ( this%buf1_re(2) - lbuf(1) )
      f1_im(2,1) = -idrh * ( this%buf1_im(2) - lbuf(2) )

    endif

    if ( idproc == nvp-1 ) then
      ir = idr / real(nrp+noff-1)
      f1_re(1,nrp) = -ir * mode * this%buf1_im(nrp)
      f1_im(1,nrp) =  ir * mode * this%buf1_re(nrp)
      f1_re(2,nrp) = -idrh * ( 3.0 * this%buf1_re(nrp) - 4.0 * this%buf1_re(nrp-1) + this%buf1_re(nrp-2) )
      f1_im(2,nrp) = -idrh * ( 3.0 * this%buf1_im(nrp) - 4.0 * this%buf1_im(nrp-1) + this%buf1_im(nrp-2) )
    else
      ir = idr / real(nrp+noff-1)
      f1_re(1,nrp) = -ir * mode * this%buf1_im(nrp)
      f1_im(1,nrp) =  ir * mode * this%buf1_re(nrp)
      f1_re(2,nrp) = -idrh * ( ubuf(1) - this%buf1_re(nrp-1) )
      f1_im(2,nrp) = -idrh * ( ubuf(2) - this%buf1_im(nrp-1) )
    endif

  endif

  call stop_tprof( 'solve beam bt' )
  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine get_solution_bt

! subroutine get_solution_bt_iter_old( this, mode )
! ! this is totally the same as get_solution_bt()

!   implicit none

!   class( field_b ), intent(inout) :: this
!   integer, intent(in) :: mode

!   integer :: i, nd1p
!   real, dimension(:,:), pointer :: f1_re => null(), f1_im => null()
!   character(len=20), save :: sname = 'get_solution_bt_iter_old'

!   call write_dbg( cls_name, sname, cls_level, 'starts' )
!   call start_tprof( 'solve plasma bt' )

!   nd1p = this%rf_re(mode)%get_ndp(1)

!   f1_re => this%rf_re(mode)%get_f1()
!   do i = 1, nd1p
!     f1_re(1,i) = this%buf(4*i-3)
!     f1_re(2,i) = this%buf(4*i-1)
!   enddo

!   if ( mode > 0 ) then
!     f1_im => this%rf_im(mode)%get_f1()
!     do i = 1, nd1p
!       f1_im(1,i) = this%buf(4*i-2)
!       f1_im(2,i) = this%buf(4*i)
!     enddo
!   endif

!   call stop_tprof( 'solve plasma bt' )
!   call write_dbg( cls_name, sname, cls_level, 'ends' )

! end subroutine get_solution_bt_iter_old

subroutine get_solution_bt_iter( this, mode )

  implicit none

  class( field_b ), intent(inout) :: this
  integer, intent(in) :: mode

  integer :: i, nrp, idproc
  real, dimension(:,:), pointer :: f1_re => null(), f1_im => null()
  character(len=20), save :: sname = 'get_solution_bt_iter'

  call write_dbg( cls_name, sname, cls_level, 'starts' )
  call start_tprof( 'solve plasma bt' )

  nrp    = this%rf_re(mode)%get_ndp(1)
  idproc = this%rf_re(mode)%pp%getlidproc()

  if ( mode == 0 ) then

    f1_re => this%rf_re(mode)%get_f1()
    do i = 1, nrp
      f1_re(1,i) = this%buf1_re(i) ! Re(Br)
      f1_re(2,i) = this%buf2_re(i) ! Re(Bphi)
    enddo

    if ( idproc == 0 ) then
      f1_re(1,1) = 0.0
      f1_re(2,1) = 0.0
    endif

  else

    f1_re => this%rf_re(mode)%get_f1()
    f1_im => this%rf_im(mode)%get_f1()
    do i = 1, nrp
      f1_re(1,i) = 0.5 * ( this%buf1_re(i) + this%buf2_re(i) ) ! Re(Br)
      f1_im(1,i) = 0.5 * ( this%buf1_im(i) + this%buf2_im(i) ) ! Im(Br)
      f1_re(2,i) = 0.5 * ( this%buf1_im(i) - this%buf2_im(i) ) ! Re(Bphi)
      f1_im(2,i) = 0.5 * (-this%buf1_re(i) + this%buf2_re(i) ) ! Im(Bphi)
    enddo

    if ( idproc == 0 ) then
      if ( mode /= 1 ) then
        f1_re(1,1) = 0.0
        f1_im(1,1) = 0.0
        f1_re(2,1) = 0.0
        f1_im(2,1) = 0.0
      endif
    endif

  endif

  call stop_tprof( 'solve plasma bt' )
  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine get_solution_bt_iter

subroutine solve_field_bz( this, jay )

  implicit none

  class( field_b ), intent(inout) :: this
  class( field_jay ), intent(inout) :: jay

  type( ufield ), dimension(:), pointer :: jay_re => null(), jay_im => null()
  integer :: i
  character(len=20), save :: sname = 'solve_field_bz'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  jay_re => jay%get_rf_re()
  jay_im => jay%get_rf_im()

  do i = 0, this%num_modes

    if ( i == 0 ) then
      call this%set_source_bz( i, jay_re(i) )
      call this%solver_bz(i)%solve( this%buf1_re )
      call this%get_solution_bz(i)
      cycle
    endif

    call this%set_source_bz( i, jay_re(i), jay_im(i) )
    call this%solver_bz(i)%solve( this%buf1_re )
    call this%solver_bz(i)%solve( this%buf1_im )
    call this%get_solution_bz(i)

  enddo

  call this%copy_gc_f1( bnd_ax = .true. ) ! to be finished

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine solve_field_bz

! subroutine solve_field_bt_old( this, rho )

!   implicit none

!   class( field_b ), intent(inout) :: this
!   class( field_rho ), intent(inout) :: rho

!   type( ufield ), dimension(:), pointer :: rho_re => null(), rho_im => null()
!   integer :: i
!   character(len=20), save :: sname = 'solve_field_bt_old'

!   call write_dbg( cls_name, sname, cls_level, 'starts' )

!   rho_re => rho%get_rf_re()
!   rho_im => rho%get_rf_im()

!   do i = 0, this%num_modes

!     if ( i == 0 ) then
!       call this%set_source_bt_old( i, rho_re(i) )
!       call this%solver_bt(i)%solve( this%buf )
!       call this%get_solution_bt_old(i)
!       cycle
!     endif

!     call this%set_source_bt_old( i, rho_re(i), rho_im(i) )
!     call this%solver_bt(i)%solve( this%buf )
!     call this%get_solution_bt_old(i)

!   enddo

!   call this%copy_gc_f1( bnd_ax = .true. )

!   call write_dbg( cls_name, sname, cls_level, 'ends' )

! end subroutine solve_field_bt_old

subroutine solve_field_bt( this, rho )

  implicit none

  class( field_b ), intent(inout) :: this
  class( field_rho ), intent(inout) :: rho

  type( ufield ), dimension(:), pointer :: rho_re => null(), rho_im => null()
  integer :: i
  character(len=20), save :: sname = 'solve_field_bt'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  rho_re => rho%get_rf_re()
  rho_im => rho%get_rf_im()

  do i = 0, this%num_modes

    if ( i == 0 ) then
      call this%set_source_bt( i, rho_re(i) )
      call this%solver_bt(i)%solve( this%buf1_re )
      call this%get_solution_bt(i)
      cycle
    endif

    call this%set_source_bt( i, rho_re(i), rho_im(i) )
    call this%solver_bt(i)%solve( this%buf1_re )
    call this%solver_bt(i)%solve( this%buf1_im )
    call this%get_solution_bt(i)

  enddo

  call this%copy_gc_f1( bnd_ax = .true. ) ! to be finished

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine solve_field_bt

! subroutine solve_field_bt_iter_old( this, djdxi, jay )

!   implicit none

!   class( field_b ), intent(inout) :: this
!   class( field_djdxi ), intent(in) :: djdxi
!   class( field_jay ), intent(inout) :: jay

!   type( ufield ), dimension(:), pointer :: jay_re => null(), jay_im => null()
!   type( ufield ), dimension(:), pointer :: djdxi_re => null(), djdxi_im => null()
!   integer :: i
!   character(len=20), save :: sname = 'solve_field_bt_iter_old'

!   call write_dbg( cls_name, sname, cls_level, 'starts' )

!   djdxi_re => djdxi%get_rf_re()
!   djdxi_im => djdxi%get_rf_im()
!   jay_re => jay%get_rf_re()
!   jay_im => jay%get_rf_im()

!   do i = 0, this%num_modes

!     if ( i == 0 ) then
!       call this%set_source_bt_iter_old( i, djdxi_re(i), jay_re(i) )
!       call this%solver_bt_iter(i)%solve( this%buf )
!       call this%get_solution_bt_iter_old(i)
!       cycle
!     endif

!     call this%set_source_bt_iter_old( i, djdxi_re(i), jay_re(i), djdxi_im(i), jay_im(i) )
!     call this%solver_bt_iter(i)%solve( this%buf )
!     call this%get_solution_bt_iter_old(i)

!   enddo

!   call this%copy_gc_f1( bnd_ax = .true. )

!   call write_dbg( cls_name, sname, cls_level, 'ends' )

! end subroutine solve_field_bt_iter_old

subroutine solve_field_bt_iter( this, djdxi, jay )

  implicit none

  class( field_b ), intent(inout) :: this
  class( field_djdxi ), intent(in) :: djdxi
  class( field_jay ), intent(inout) :: jay

  type( ufield ), dimension(:), pointer :: jay_re => null(), jay_im => null()
  type( ufield ), dimension(:), pointer :: djdxi_re => null(), djdxi_im => null()
  integer :: i
  character(len=20), save :: sname = 'solve_field_bt_iter'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  djdxi_re => djdxi%get_rf_re()
  djdxi_im => djdxi%get_rf_im()
  jay_re => jay%get_rf_re()
  jay_im => jay%get_rf_im()

  do i = 0, this%num_modes

    if ( i == 0 ) then
      call this%set_source_bt_iter( i, djdxi_re(i), jay_re(i) )
      call this%solver_bplus(i)%solve( this%buf1_re )
      call this%solver_bminus(i)%solve( this%buf2_re )
      call this%get_solution_bt_iter(i)
      cycle
    endif

    call this%set_source_bt_iter( i, djdxi_re(i), jay_re(i), djdxi_im(i), jay_im(i) )
    call this%solver_bplus(i)%solve( this%buf1_re )
    call this%solver_bplus(i)%solve( this%buf1_im )
    call this%solver_bminus(i)%solve( this%buf2_re )
    call this%solver_bminus(i)%solve( this%buf2_im )
    call this%get_solution_bt_iter(i)

  enddo

  call this%copy_gc_f1( bnd_ax = .true. ) ! to be finished

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine solve_field_bt_iter

end module field_b_class
module field_psi_class

use parallel_pipe_class
use grid_class
use field_class
use field_solver_class
use field_src_class
use ufield_class
use param
use sysutil
use mpi

implicit none

private

character(len=20), parameter :: cls_name = "field_psi"
integer, parameter :: cls_level = 3

public :: field_psi

type, extends( field ) :: field_psi

  ! private

  class( field_solver ), dimension(:), pointer :: solver => null()
  real, dimension(:), pointer :: buf_re => null(), buf_im => null() ! buffer for source term

  contains

  generic :: new => init_field_psi
  procedure :: del => end_field_psi
  ! generic :: read_input => read_input_field_psi
  procedure :: solve => solve_field_psi

  procedure, private :: init_field_psi
  procedure, private :: end_field_psi
  procedure, private :: set_source
  procedure, private :: get_solution

end type field_psi

contains

subroutine init_field_psi( this, pp, gp, num_modes, part_shape, boundary, iter_tol )

  implicit none

  class( field_psi ), intent(inout) :: this
  class( parallel_pipe ), intent(in), pointer :: pp
  class( grid ), intent(in), pointer :: gp
  integer, intent(in) :: num_modes, part_shape, boundary
  real, intent(in) :: iter_tol

  integer, dimension(2,2) :: gc_num
  integer :: dim, i, nrp
  real :: dr
  character(len=20), save :: sname = "init_field_psi"

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  nrp = gp%get_ndp(1)
  dr = gp%get_dr()

  select case ( part_shape )

  case ( p_ps_linear )

    gc_num(:,1) = (/1, 1/)
    gc_num(:,2) = (/2, 1/)

  case ( p_ps_quadratic )

    call write_err( "Quadratic particle shape not implemented." )

  case default

    call write_err( "Invalid particle shape." )

  end select

  dim = 1
  ! call initialization routine of the parent class
  call this%field%new( pp, gp, dim, num_modes, gc_num, entity=p_entity_plasma )

  ! initialize solver
  allocate( this%solver( 0:num_modes ) )
  do i = 0, num_modes
    call this%solver(i)%new( pp, gp, i, dr, &
      kind=p_fk_psi, bnd=boundary, stype=p_hypre_cycred, tol=iter_tol )
  enddo

  allocate( this%buf_re(nrp), this%buf_im(nrp) )

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine init_field_psi

subroutine end_field_psi( this )

  implicit none

  class( field_psi ), intent(inout) :: this

  integer :: i
  character(len=20), save :: sname = 'end_field_psi'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  do i = 0, this%num_modes
    call this%solver(i)%del()
  enddo
  deallocate( this%solver )

  if ( associated( this%buf_re ) ) deallocate( this%buf_re )
  if ( associated( this%buf_im ) ) deallocate( this%buf_im )

  call this%field%del()

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine end_field_psi

subroutine set_source( this, mode, q_re, q_im )

  implicit none

  class( field_psi ), intent(inout) :: this
  class( ufield ), intent(in) :: q_re
  class( ufield ), intent(in), optional :: q_im
  integer, intent(in) :: mode

  integer :: i, nrp
  real, dimension(:,:), pointer :: f1_re => null(), f1_im => null()
  real :: dr2
  character(len=20), save :: sname = 'set_source'

  call write_dbg( cls_name, sname, cls_level, 'starts' )
  call start_tprof( 'solve psi' )

  nrp = q_re%get_ndp(1)
  dr2 = this%dr**2

  f1_re => q_re%get_f1()
  this%buf_re = 0.0
  if ( present(q_im) ) then
    f1_im => q_im%get_f1()
    this%buf_im = 0.0
  endif

  if ( mode == 0 ) then

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

  call stop_tprof( 'solve psi' )
  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine set_source

subroutine get_solution( this, mode )

  implicit none

  class( field_psi ), intent(inout) :: this
  integer, intent(in) :: mode

  integer :: i, nrp, idproc
  real, dimension(:,:), pointer :: f1_re => null(), f1_im => null()
  character(len=20), save :: sname = 'get_solution'

  call write_dbg( cls_name, sname, cls_level, 'starts' )
  call start_tprof( 'solve psi' )

  nrp    = this%rf_re(mode)%get_ndp(1)
  idproc = this%rf_re(mode)%pp%getlidproc()

  f1_re => this%rf_re(mode)%get_f1()
  do i = 1, nrp
    f1_re(1,i) = this%buf_re(i)
  enddo

  if ( mode > 0 ) then
    f1_im => this%rf_im(mode)%get_f1()
    do i = 1, nrp
      f1_im(1,i) = this%buf_im(i)
    enddo

    ! psi(m>0) vanishes on axis
    if ( idproc == 0 ) then
      f1_re(1,1) = 0.0
      f1_im(1,1) = 0.0
    endif
  endif

  call stop_tprof( 'solve psi' )
  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine get_solution

! subroutine get_solution( this, mode, idx, q_ax )

!   implicit none

!   class( field_psi ), intent(inout) :: this
!   integer, intent(in) :: mode, idx
!   real, intent(in), optional :: q_ax

!   integer :: i, j, nrp, noff
!   real, dimension(:,:), pointer :: f1_re => null(), f1_im => null()
!   real, dimension(:,:,:), pointer :: f2_re => null(), f2_im => null()
!   real :: r, r2, rmax, rmax2
!   ! longitudinal smooth
!   ! real, dimension(0:4), save :: s = (/ 0.429447852760736, &
!   !                                    0.343558282208589, &
!   !                                    0.171779141104294, &
!   !                                    0.049079754601227, &
!   !                                    0.006134969325153 /)
!   real, dimension(0:4), save :: s = (/ 0.38, 0.32, 0.2, 0.08, 0.02 /)
!   character(len=20), save :: sname = 'get_solution'

!   call write_dbg( cls_name, sname, cls_level, 'starts' )
!   call start_tprof( 'solve psi' )

!   nrp   = this%rf_re(mode)%get_ndp(1)
!   noff  = this%rf_re(mode)%get_noff(1)
!   rmax  = (this%rf_re(mode)%get_nd(1)+0.5) * this%dr
!   rmax2 = rmax**2

!   if ( mode == 0 .and. present(q_ax) ) then
!     select case ( this%solver(0)%bnd )

!     case ( p_bnd_zero, p_bnd_open )

!       ! add contribution of the source terms on axis
!       do i = 1, nrp
!         r = (i+noff-0.5) * this%dr
!         ! this%buf_re(i) = this%buf_re(i) + q_ax * log(r/rmax)
!         this%buf_re(i) = this%buf_re(i)
!       enddo

!     end select
!   endif

!   f1_re => this%rf_re(mode)%get_f1()
!   f2_re => this%rf_re(mode)%get_f2()
!   do i = 1, nrp
!     ! f1_re(1,i) = this%buf_re(i)
!     f1_re(1,i) = s(0) * this%buf_re(i)
!     do j = 1, 4
!       f1_re(1,i) = f1_re(1,i) + s(j) * f2_re(1,i,idx-j)
!     enddo
!   enddo

!   if ( mode > 0 ) then
!     f1_im => this%rf_im(mode)%get_f1()
!     f2_im => this%rf_im(mode)%get_f2()
!     do i = 1, nrp
!       ! f1_im(1,i) = this%buf_im(i)
!       f1_im(1,i) = s(0) * this%buf_im(i)
!       f1_im(1,i) = f1_im(1,i) + s(j) * f2_im(1,i,idx-j)
!     enddo
!   endif

!   call stop_tprof( 'solve psi' )
!   call write_dbg( cls_name, sname, cls_level, 'ends' )

! end subroutine get_solution

subroutine solve_field_psi( this, q )

  implicit none

  class( field_psi ), intent(inout) :: this
  class( field_rho ), intent(inout) :: q

  type( ufield ), dimension(:), pointer :: q_re => null(), q_im => null()
  integer :: mode
  character(len=20), save :: sname = 'solve_field_psi'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  q_re => q%get_rf_re()
  q_im => q%get_rf_im()

  do mode = 0, this%num_modes

    if ( mode == 0 ) then
      call this%set_source( mode, q_re(mode) )
      call this%solver(mode)%solve( this%buf_re )
      call this%get_solution(mode)
      cycle
    endif

    call this%set_source( mode, q_re(mode), q_im(mode) )
    call this%solver(mode)%solve( this%buf_re )
    call this%solver(mode)%solve( this%buf_im )
    call this%get_solution(mode)

  enddo

  call this%copy_gc_f1()

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine solve_field_psi

! subroutine solve_field_psi( this, q, idx )

!   implicit none

!   class( field_psi ), intent(inout) :: this
!   class( field_rho ), intent(inout) :: q
!   integer, intent(in) :: idx

!   type( ufield ), dimension(:), pointer :: q_re => null(), q_im => null()
!   integer :: i
!   character(len=20), save :: sname = 'solve_field_psi'

!   call write_dbg( cls_name, sname, cls_level, 'starts' )

!   call q%get_q_ax()

!   q_re => q%get_rf_re()
!   q_im => q%get_rf_im()

!   do i = 0, this%num_modes

!     if ( i == 0 ) then
!       call this%set_source( i, q%q_ax, q_re(i) )
!       call this%solver(i)%solve( this%buf_re )
!       call this%get_solution( i, idx, q%q_ax )
!       cycle
!     endif

!     call this%set_source( i, q%q_ax, q_re(i), q_im(i) )
!     call this%solver(i)%solve( this%buf_re )
!     call this%solver(i)%solve( this%buf_im )
!     call this%get_solution(i, idx)

!   enddo

!   call this%copy_gc_f1()

!   call write_dbg( cls_name, sname, cls_level, 'ends' )

! end subroutine solve_field_psi

end module field_psi_class
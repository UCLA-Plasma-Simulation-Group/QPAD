module field_b_class

use field_class
use field_psi_class
use field_src_class
use field_solver_class
use field_solver_all_modes_class
use ufield_class
use param
use sysutil_module
use parallel_module
use options_class
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
  class( field_solver ), dimension(:), pointer :: solver_bplus   => null()
  class( field_solver ), dimension(:), pointer :: solver_bminus  => null()
  class( field_solver_all_modes ), dimension(:), pointer :: solver_coef  => null()
  class( field_solver_all_modes ), pointer :: solver_bp => null()
  class( field_solver_all_modes ), pointer :: solver_bm => null()

  real, dimension(:), pointer :: buf1_re => null(), buf1_im => null()
  real, dimension(:), pointer :: buf2_re => null(), buf2_im => null()
  real, dimension(:,:,:), pointer :: buf3 => null()
  real, dimension(:,:,:), pointer :: buf4 => null()
  real, dimension(:,:,:), pointer :: buf5 => null()
  real, dimension(:,:,:), pointer :: buf6 => null()
  real, dimension(:), pointer :: buf    => null()

  contains

  generic :: solve => solve_field_bz, solve_field_bt, solve_field_bt_iter
  generic :: new   => init_field_b

  procedure :: init_field_b
  procedure :: del => end_field_b
  procedure :: alloc => alloc_field_b
  procedure, private :: set_source_bz
  procedure, private :: set_source_bt
  procedure, private :: set_source_bt_iter
  procedure, private :: get_solution_bz
  procedure, private :: get_solution_bt
  procedure, private :: get_solution_bt_iter
  procedure, private :: solve_field_bz
  procedure, private :: solve_field_bt
  procedure, private :: solve_field_bt_iter

end type field_b

contains

subroutine alloc_field_b( this, max_mode )

  implicit none

  class( field_b ), intent(inout) :: this
  integer, intent(in) :: max_mode

  if ( .not. associated( this%solver_bz ) ) then
    allocate( field_solver :: this%solver_bz(0:max_mode) )
  endif

  if ( .not. associated( this%solver_bt ) ) then
    allocate( field_solver :: this%solver_bt(0:max_mode) )
  endif

  if ( .not. associated( this%solver_bplus ) ) then
    allocate( field_solver :: this%solver_bplus(0:max_mode) )
  endif

  if ( .not. associated( this%solver_bminus ) ) then
    allocate( field_solver :: this%solver_bminus(0:max_mode) )
  endif

  if ( .not. associated( this%solver_coef ) ) then
    allocate( field_solver_all_modes :: this%solver_coef(0:1) )
  endif

  if ( .not. associated( this%solver_bp ) ) then
    allocate( field_solver_all_modes :: this%solver_bp)
  endif

  if (.not. associated(this%solver_bm)) then
    allocate(field_solver_all_modes :: this%solver_bm)
  endif

end subroutine alloc_field_b

subroutine init_field_b( this, opts, max_mode, part_shape, boundary, entity )

  implicit none

  class( field_b ), intent(inout) :: this
  type( options ), intent(in) :: opts
  integer, intent(in) :: max_mode, part_shape, entity, boundary

  integer, dimension(2,2) :: gc_num
  integer :: dim, i, nrp
  real :: dr
  character(len=20), save :: sname = "init_field_b"

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  nrp = opts%get_ndp(1)
  dr  = opts%get_dr()

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
  call this%field%new( opts, dim, max_mode, gc_num, entity )

  ! initialize solver
  select case ( entity )

  case ( p_entity_plasma )

    do i = 0, max_mode
      call this%solver_bz(i)%new( opts, i, dr, kind=p_fk_bz, &
        bnd=boundary, stype=p_hypre_cycred )
      call this%solver_bplus(i)%new( opts, i, dr, kind=p_fk_bplus, &
        bnd=boundary, stype=p_hypre_cycred )
      call this%solver_bminus(i)%new( opts, i, dr, kind=p_fk_bminus, &
        bnd=boundary, stype=p_hypre_cycred )
    enddo 

    if (max_mode > 0) then
      write(2,*) 'B Initializing 2'
      call this%solver_coef(0)%new( opts, max_mode, dr, kind=p_fk_coef, &
          bnd=boundary, stype=p_hypre_cycred )
      call this%solver_coef(1)%new( opts, max_mode, dr, kind=p_fk_coef, &
          bnd=boundary, stype=p_hypre_cycred )
    endif

    write(2,*) 'B Initializing 3'
    call this%solver_bm%new( opts, max_mode, dr, kind=p_fk_all_B_minus, &
        bnd=boundary, stype=p_hypre_cycred )

    write(2,*) 'B Initializing 4'
    call this%solver_bp%new( opts, max_mode, dr, kind=p_fk_all_B_plus, &
        bnd=boundary, stype=p_hypre_cycred )
!     endif

    allocate( this%buf1_re(nrp), this%buf1_im(nrp) )
    allocate( this%buf2_re(nrp), this%buf2_im(nrp) )
    allocate( this%buf3(1 + 2*max_mode,nrp,2) )
    allocate( this%buf4(1 + 2*max_mode,nrp,2) )
    allocate( this%buf5(1 + 2*max_mode,nrp,2) )
    allocate( this%buf6(1 + 2*max_mode,nrp,2) )

  case ( p_entity_beam )

    do i = 0, max_mode
      call this%solver_bt(i)%new( opts, i, dr, kind=p_fk_bt, &
        bnd=boundary, stype=p_hypre_cycred )
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

  case ( p_entity_plasma )
    do i = 0, this%max_mode
      call this%solver_bz(i)%del()
      call this%solver_bplus(i)%del()
      call this%solver_bminus(i)%del()
    enddo
    deallocate( this%solver_bz )
    deallocate( this%solver_bplus )
    deallocate( this%solver_bminus )
    if (this%max_mode > 0 ) then
      call this%solver_coef(0)%del()
      call this%solver_coef(1)%del()
      deallocate( this%solver_coef )
    endif
    call this%solver_bp%del()
    call this%solver_bm%del()
!     deallocate( this%solver_coef )
    deallocate( this%solver_bp )
    deallocate( this%solver_bm )
!     endif

  case ( p_entity_beam )
    do i = 0, this%max_mode
      call this%solver_bt(i)%del()
    enddo
    deallocate( this%solver_bt )
  end select

  if ( associated( this%buf1_re ) ) deallocate( this%buf1_re )
  if ( associated( this%buf1_im ) ) deallocate( this%buf1_im )
  if ( associated( this%buf2_re ) ) deallocate( this%buf2_re )
  if ( associated( this%buf2_im ) ) deallocate( this%buf2_im )
  if ( associated( this%buf ) ) deallocate( this%buf )
  if ( associated( this%buf3 ) ) deallocate( this%buf3 )
  if ( associated( this%buf4 ) ) deallocate( this%buf4 )
  if ( associated( this%buf5 ) ) deallocate( this%buf5 )
  if ( associated( this%buf6 ) ) deallocate( this%buf6 )

  call this%field%del()

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine end_field_b

subroutine set_source_bz( this, mode, jay_re, jay_im )

  implicit none

  class( field_b ), intent(inout) :: this
  class( ufield ), intent(in) :: jay_re
  class( ufield ), intent(in), optional :: jay_im
  integer, intent(in) :: mode

  integer :: i, nrp, noff, idproc, nvp
  real, dimension(:,:), pointer :: f1_re => null(), f1_im => null()
  real :: idrh, idr, dr, dr2, ir
  character(len=20), save :: sname = 'set_source_bz'

  call write_dbg( cls_name, sname, cls_level, 'starts' )
  call start_tprof( 'solve bz' )

  nrp    = jay_re%get_ndp(1)
  idr    = 1.0 / this%dr
  idrh   = 0.5 * idr
  noff   = jay_re%get_noff(1)
  nvp    = num_procs_loc()
  idproc = id_proc_loc()
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
      if ( mod(mode,2) == 0 ) then
        this%buf1_re(1) = -2.0 * idr * f1_re(2,2) - mode * idr * f1_im(1,2)
        this%buf1_im(1) = -2.0 * idr * f1_im(2,2) + mode * idr * f1_re(1,2)
      else
        this%buf1_re(1) = 0.0
        this%buf1_im(1) = 0.0

        ! since j_phi(m=1) is multiplied by factor 8 on axis, the derivative on index=2 is
        ! calculated using forward difference
        if ( mode == 1 ) then
          ir = idr
          this%buf1_re(2) = -idr * ( f1_re(2,3) - f1_re(2,2) ) - ir * f1_re(2,2) - mode * ir * f1_im(1,2)
          this%buf1_im(2) = -idr * ( f1_im(2,3) - f1_im(2,2) ) - ir * f1_im(2,2) + mode * ir * f1_re(1,2)
        endif
      endif
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

subroutine set_source_bt( this, mode, q_re, q_im )

  implicit none

  class( field_b ), intent(inout) :: this
  class( ufield ), intent(in) :: q_re
  class( ufield ), intent(in), optional :: q_im
  integer, intent(in) :: mode

  integer :: i, nrp, noff, comm
  real, dimension(:,:), pointer :: f1_re => null(), f1_im => null()
  real ::dr, dr2, rmax
  character(len=20), save :: sname = 'set_source_bt'

  call write_dbg( cls_name, sname, cls_level, 'starts' )
  call start_tprof( 'solve beam bt' )

  nrp   = q_re%get_ndp(1)
  noff  = q_re%get_noff(1)
  dr    = this%dr
  dr2   = dr**2
  comm  = comm_loc()
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

subroutine set_source_bt_iter( this, mode, djdxi_re, jay_re, djdxi_im, jay_im )

  implicit none

  class( field_b ), intent(inout) :: this
  class( ufield ), intent(in) :: djdxi_re, jay_re
  class( ufield ), intent(in), optional :: djdxi_im, jay_im
  integer, intent(in) :: mode

  integer :: i, nrp, nvp, idproc, noff
  real, dimension(:,:), pointer :: f1_re => null(), f1_im => null()
  real, dimension(:,:), pointer :: f2_re => null(), f2_im => null()
  real, dimension(:,:), pointer :: f3_re => null(), f3_im => null()
  real :: idrh, idr, ir
  real :: s1_re, s1_im, s2_re, s2_im
  character(len=20), save :: sname = 'set_source_bt_iter'

  nvp    = num_procs_loc()
  idproc = id_proc_loc()
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
      this%buf1_re(i) = -f1_re(2,i)  ! Re(Br)
!       this%buf1_re(i) = -f1_re(2,i) - f3_re(1,i) ! Re(Br)
      this%buf2_re(i) = f1_re(1,i) + idrh * ( f2_re(3,i+1) - f2_re(3,i-1) )! Re(Bphi)
!       this%buf2_re(i) = f1_re(1,i) + idrh * ( f2_re(3,i+1) - f2_re(3,i-1) ) - f3_re(2,i) ! Re(Bphi)
    enddo

    ! calculate the derivatives at the boundary and axis
    if ( idproc == 0 ) then
      this%buf1_re(1) = 0.0
      this%buf2_re(1) = 0.0
      ! since Jz(m=0) is multiplied by factor 8 on axis, the derivative on index=2 is
      ! calculated using forward difference
      this%buf2_re(2) =  f1_re(1,2) + idr * ( f2_re(3,3) - f2_re(3,2) ) 
!       this%buf2_re(2) =  f1_re(1,2) + idr * ( f2_re(3,3) - f2_re(3,2) ) - f3_re(2,2)
    else
      this%buf1_re(1) = -f1_re(2,1) 
!       this%buf1_re(1) = -f1_re(2,1) - f3_re(1,1)
      this%buf2_re(1) =  f1_re(1,1) + idrh * ( f2_re(3,2) - f2_re(3,0) ) 
!       this%buf2_re(1) =  f1_re(1,1) + idrh * ( f2_re(3,2) - f2_re(3,0) ) - f3_re(2,1)
    endif

    if ( idproc == nvp-1 ) then
      this%buf1_re(nrp) = -f1_re(2,nrp) 
!       this%buf1_re(nrp) = -f1_re(2,nrp) - f3_re(1,nrp)
      this%buf2_re(nrp) =  f1_re(1,nrp) + idrh * ( 3.0 * f2_re(3,nrp) - 4.0 * f2_re(3,nrp-1) + f2_re(3,nrp-2) )
!       this%buf2_re(nrp) =  f1_re(1,nrp) + idrh * ( 3.0 * f2_re(3,nrp) - 4.0 * f2_re(3,nrp-1) + f2_re(3,nrp-2) ) - f3_re(2,nrp)
    endif

  elseif ( mode > 0 .and. present( jay_im ) .and. present( djdxi_im ) ) then

    do i = 2, nrp
      ir = idr / real(i+noff-1)
      s1_re = -f1_re(2,i) + mode * f2_im(3,i) * ir
      s1_im = -f1_im(2,i) - mode * f2_re(3,i) * ir
      s2_re =  f1_re(1,i) + idrh * ( f2_re(3,i+1) - f2_re(3,i-1) )
      s2_im =  f1_im(1,i) + idrh * ( f2_im(3,i+1) - f2_im(3,i-1) )
!       this%buf1_re(i) = s1_re - s2_im - f3_re(1,i) + f3_im(2,i) ! Re(B_plus)
!       this%buf1_im(i) = s1_im + s2_re - f3_im(1,i) - f3_re(2,i) ! Im(B_plus)
!       this%buf2_re(i) = s1_re + s2_im - f3_re(1,i) - f3_im(2,i) ! Re(B_minus)
!       this%buf2_im(i) = s1_im - s2_re - f3_im(1,i) + f3_re(2,i) ! Im(B_minus)
      this%buf1_re(i) = s1_re - s2_im ! Re(B_plus)
      this%buf1_im(i) = s1_im + s2_re ! Im(B_plus)
      this%buf2_re(i) = s1_re + s2_im ! Re(B_minus)
      this%buf2_im(i) = s1_im - s2_re ! Im(B_minus)
    enddo

    ! calculate the derivatives at the boundary and axis
    if ( idproc == 0 ) then

      if ( mode == 1 ) then
        s1_re = -f1_re(2,1) + idr * mode * f2_im(3,2)
        s1_im = -f1_im(2,1) - idr * mode * f2_re(3,2)
        s2_re =  f1_re(1,1) + idr * f2_re(3,2)
        s2_im =  f1_im(1,1) + idr * f2_im(3,2)
      else
        if ( mod(mode,2) == 0 ) then
          s1_re = 0.0
          s1_im = 0.0
          s2_re = 0.0
          s2_im = 0.0
        else
          s1_re = idr * mode * f2_im(3,2)
          s1_im = idr * mode * f2_re(3,2)
          s2_re = idr * f2_re(3,2)
          s2_im = idr * f2_im(3,2)
        endif
      endif

!       this%buf1_re(1) = s1_re - s2_im - f3_re(1,1) + f3_im(2,1) ! Re(B_plus)
!       this%buf1_im(1) = s1_im + s2_re - f3_im(1,1) - f3_re(2,1) ! Im(B_plus)
!       this%buf2_re(1) = s1_re + s2_im - f3_re(1,1) - f3_im(2,1) ! Re(B_minus)
!       this%buf2_im(1) = s1_im - s2_re - f3_im(1,1) + f3_re(2,1) ! Im(B_minus)
      this%buf1_re(1) = s1_re - s2_im  ! Re(B_plus)
      this%buf1_im(1) = s1_im + s2_re  ! Im(B_plus)
      this%buf2_re(1) = s1_re + s2_im  ! Re(B_minus)
      this%buf2_im(1) = s1_im - s2_re  ! Im(B_minus)

    else

      ir = idr / real(noff)
      s1_re = -f1_re(2,1) + mode * f2_im(3,1) * ir
      s1_im = -f1_im(2,1) - mode * f2_re(3,1) * ir
      s2_re =  f1_re(1,1) + idrh * ( f2_re(3,2) - f2_re(3,0) )
      s2_im =  f1_im(1,1) + idrh * ( f2_im(3,2) - f2_im(3,0) )

!       this%buf1_re(1) = s1_re - s2_im - f3_re(1,1) + f3_im(2,1) ! Re(B_plus)
!       this%buf1_im(1) = s1_im + s2_re - f3_im(1,1) - f3_re(2,1) ! Im(B_plus)
!       this%buf2_re(1) = s1_re + s2_im - f3_re(1,1) - f3_im(2,1) ! Re(B_minus)
!       this%buf2_im(1) = s1_im - s2_re - f3_im(1,1) + f3_re(2,1) ! Im(B_minus)
      this%buf1_re(1) = s1_re - s2_im  ! Re(B_plus)
      this%buf1_im(1) = s1_im + s2_re  ! Im(B_plus)
      this%buf2_re(1) = s1_re + s2_im  ! Re(B_minus)
      this%buf2_im(1) = s1_im - s2_re  ! Im(B_minus)

    endif

    if ( idproc == nvp-1 ) then

      ir = idr / real(nrp+noff-1)
      s1_re = -f1_re(2,nrp) + mode * f2_im(3,nrp) * ir
      s1_im = -f1_im(2,nrp) - mode * f2_re(3,nrp) * ir
      s2_re =  f1_re(1,nrp) + idrh * ( 3.0 * f2_re(3,nrp) - 4.0 * f2_re(3,nrp-1) + f2_re(3,nrp-2) )
      s2_im =  f1_im(1,nrp) + idrh * ( 3.0 * f2_im(3,nrp) - 4.0 * f2_im(3,nrp-1) + f2_im(3,nrp-2) )
!       this%buf1_re(nrp) = s1_re - s2_im - f3_re(1,nrp) + f3_im(2,nrp) ! Re(B_plus)
!       this%buf1_im(nrp) = s1_im + s2_re - f3_im(1,nrp) - f3_re(2,nrp) ! Im(B_plus)
!       this%buf2_re(nrp) = s1_re + s2_im - f3_re(1,nrp) - f3_im(2,nrp) ! Re(B_minus)
!       this%buf2_im(nrp) = s1_im - s2_re - f3_im(1,nrp) + f3_re(2,nrp) ! Im(B_minus)
      this%buf1_re(nrp) = s1_re - s2_im  ! Re(B_plus)
      this%buf1_im(nrp) = s1_im + s2_re  ! Im(B_plus)
      this%buf2_re(nrp) = s1_re + s2_im  ! Re(B_minus)
      this%buf2_im(nrp) = s1_im - s2_re  ! Im(B_minus)

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
    f1_re(3,1) = 0.0
    f1_im(3,1) = 0.0
  endif

  call stop_tprof( 'solve bz' )
  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine get_solution_bz

subroutine get_solution_bt( this, mode )

  implicit none

  class( field_b ), intent(inout) :: this
  integer, intent(in) :: mode

  integer :: i, nrp, idproc, nvp, noff, ierr, comm, msgid1, msgid2
  real :: idr, idrh, ir
  real, dimension(2), save :: lbuf, ubuf
  real, dimension(:,:), pointer :: f1_re => null(), f1_im => null()
  integer, dimension(MPI_STATUS_SIZE) :: stat
  character(len=20), save :: sname = 'get_solution_bt'

  call write_dbg( cls_name, sname, cls_level, 'starts' )
  call start_tprof( 'solve beam bt' )

  nrp    = this%rf_re(mode)%get_ndp(1)
  idproc = id_proc_loc()
  nvp    = num_procs_loc()
  noff   = this%rf_re(mode)%get_noff(1)
  comm   = comm_loc()
  idr    = 1.0 / this%dr
  idrh   = 0.5 * idr

  lbuf = 0.0; ubuf = 0.0

  ! copy the guard cell of buffer

  ! forward message passing
  ! receiver
  if ( idproc > 0 ) then
    call MPI_IRECV( lbuf(1), 1, p_dtype_real, idproc-1, 1, comm, msgid1, ierr )
    call MPI_IRECV( lbuf(2), 1, p_dtype_real, idproc-1, 2, comm, msgid2, ierr )
  endif
  ! sender
  if ( idproc < nvp-1 ) then
    call MPI_SEND( this%buf1_re(nrp), 1, p_dtype_real, idproc+1, 1, comm, ierr )
    call MPI_SEND( this%buf1_im(nrp), 1, p_dtype_real, idproc+1, 2, comm, ierr )
  endif
  ! wait receiving finish
  if ( idproc > 0 ) then
    call MPI_WAIT( msgid1, stat, ierr )
    call MPI_WAIT( msgid2, stat, ierr )
  endif

  ! backward message passing
  ! receiver
  if ( idproc < nvp-1 ) then
    call MPI_IRECV( ubuf(1), 1, p_dtype_real, idproc+1, 1, comm, msgid1, ierr )
    call MPI_IRECV( ubuf(2), 1, p_dtype_real, idproc+1, 2, comm, msgid2, ierr )
  endif
  ! sender
  if ( idproc > 0 ) then
    call MPI_SEND( this%buf1_re(1), 1, p_dtype_real, idproc-1, 1, comm, ierr )
    call MPI_SEND( this%buf1_im(1), 1, p_dtype_real, idproc-1, 2, comm, ierr )
  endif
  ! wait receiving finish
  if ( idproc < nvp-1 ) then
    call MPI_WAIT( msgid1, stat, ierr )
    call MPI_WAIT( msgid2, stat, ierr )
  endif


  if ( mode == 0 ) then

    f1_re => this%rf_re(mode)%get_f1()
    do i = 2, nrp-1
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

      if ( mode == 1 ) then
        f1_re(1,1) = -idr * mode * this%buf1_im(2)
        f1_im(1,1) =  idr * mode * this%buf1_re(2)
        f1_re(2,1) = -idr * this%buf1_re(2)
        f1_im(2,1) = -idr * this%buf1_im(2)
      else
        f1_re(1,1) = 0.0
        f1_im(1,1) = 0.0
        f1_re(2,1) = 0.0
        f1_im(2,1) = 0.0
        ! if ( mod(mode,2) == 0 ) then
        !   f1_re(1,1) = 0.0
        !   f1_im(1,1) = 0.0
        !   f1_re(2,1) = 0.0
        !   f1_im(2,1) = 0.0
        ! else
        !   f1_re(1,1) = -idr * mode * this%buf1_im(2)
        !   f1_im(1,1) =  idr * mode * this%buf1_re(2)
        !   f1_re(2,1) = -idr * this%buf1_re(2)
        !   f1_im(2,1) = -idr * this%buf1_im(2)
        ! endif
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
  idproc = id_proc_loc()

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

  do i = 0, this%max_mode

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

  call this%copy_gc_f1()

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine solve_field_bz

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

  do i = 0, this%max_mode

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

  call this%copy_gc_f1()

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine solve_field_bt

subroutine solve_field_bt_iter( this, djdxi, jay, psi, q, qn) 

  implicit none

  class( field_b ), intent(inout) :: this
  class( field_djdxi ), intent(in) :: djdxi
  class( field_jay ), intent(inout) :: jay
  class( field_psi ), intent(in) :: psi
  class( field_rho ), intent(in) :: q
  class( field_rho ), intent(in) :: qn

  type( ufield ), dimension(:), pointer :: jay_re => null(), jay_im => null()
  type( ufield ), dimension(:), pointer :: djdxi_re => null(), djdxi_im => null()
  type( ufield ), dimension(:), pointer :: psi_re => null(), psi_im => null()
  type( ufield ), dimension(:), pointer :: q_re => null(), q_im => null()
  type( ufield ), dimension(:), pointer :: qn_re => null(), qn_im => null()
  type( ufield ), pointer :: psi_re1 => null(), psi_im1 => null()
  type( ufield ), pointer :: q_re1 => null(), q_im1 => null()
  type( ufield ), pointer :: qn_re1 => null(), qn_im1 => null()
  real, dimension(:,:), pointer :: f1_re => null()
  real, dimension(:,:), pointer :: f2_re => null()
  real, dimension(:,:), pointer :: f3_re => null()
  real, dimension(3,250,2) :: a1
  real, dimension(3,250,2) :: a2
  real, dimension(3,250,2) :: a3
  real :: qbm
  integer :: i
  character(len=20), save :: sname = 'solve_field_bt_iter'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  djdxi_re => djdxi%get_rf_re()
  djdxi_im => djdxi%get_rf_im()
  jay_re => jay%get_rf_re()
  jay_im => jay%get_rf_im()
  psi_re => psi%get_rf_re()
  psi_im => psi%get_rf_im()
  q_re => q%get_rf_re()
  q_im => q%get_rf_im()
  qn_re => qn%get_rf_re()
  qn_im => qn%get_rf_im()

  this%buf3 = 0.0
  this%buf4 = 0.0
  this%buf5 = 0.0
  this%buf6 = 0.0
  a1 = 0.0
  a2 = 0.0
  a3 = 0.0
!   if ( this%max_mode == 0 ) then
  f1_re => psi_re(0)%get_f1()
  f2_re => q_re(0)%get_f1()
  f3_re => qn_re(0)%get_f1()
  write(2,*) 'this%solver_coef%solve max_mode=0'
!   this%buf3(1+this%max_mode,:,1) = (f2_re(1,:)-f3_re(1,:))/(1+f1_re(1,:)) - 1/1836.5 * f3_re(1,:)/(1-1/1836.5*f1_re(1,:))
  a1(2,:,1) = (f2_re(1,:)-f3_re(1,:))/(1+f1_re(1,:)) - 1/1836.5 * f3_re(1,:)/(1-1/1836.5*f1_re(1,:))
!   else
!   solve ele coef
  qbm = -1.0
  write(2,*) 'this%solver_coef(0)%solve'
!     call this%solver_coef(0)%solve(src = this%buf3, psi_re = psi_re, psi_im = psi_im, q_re = q_re, q_im = q_im, qbm = qbm )
  call this%solver_coef(0)%solve(src = this%buf3, psi_re = psi_re, psi_im = psi_im, q_re = q_re, q_im = q_im,&
    qn_re = qn_re, qn_im = qn_im, qbm = qbm )
!     !   solve ion coef
  qbm = 1/1836.5
  write(2,*) 'this%solver_coef(1)%solve'
  call this%solver_coef(1)%solve(src = this%buf4, psi_re = psi_re, psi_im = psi_im, q_re = qn_re, q_im = qn_im, qbm = qbm )
  ! add coef
  write(2,*) sum(this%buf3),'this%buf3'
  write(2,*) sum(this%buf4),'this%buf4'
  this%buf3(:,:,:) = this%buf3(:,:,:) + this%buf4(:,:,:)
!   endif 
  a2(2,:,1) = this%buf3(2,:,1)
  a3(2,:,1) = a2(2,:,1) - a1(2,:,1)
  write(2,*) sum(this%buf3),'this%buf3 sum'
  write(2,*) sum(this%buf3(2,:,:)),'this%buf3 sumv1'
  write(2,*) a3(2,:,1),'a3'
    !   set source
  do i = 0, this%max_mode
    if ( i == 0 ) then
      call this%set_source_bt_iter( i, djdxi_re(i), jay_re(i) )
      this%buf5(1+this%max_mode,:,1) = this%buf1_re(:)
      this%buf6(1+this%max_mode,:,1) = this%buf2_re(:)
    else
      call this%set_source_bt_iter( i, djdxi_re(i), jay_re(i), djdxi_im(i), jay_im(i) )
      this%buf5(1+this%max_mode + i, :, 1) = this%buf1_re
      this%buf5(1+this%max_mode + i, :, 2) = this%buf1_im
      this%buf5(1+this%max_mode - i, :, 1) = this%buf1_re
      this%buf5(1+this%max_mode - i, :, 2) = -this%buf1_im
      this%buf6(1+this%max_mode + i, :, 1) = this%buf2_re
      this%buf6(1+this%max_mode + i, :, 2) = this%buf2_im
      this%buf6(1+this%max_mode - i, :, 1) = this%buf2_re
      this%buf6(1+this%max_mode - i, :, 2) = -this%buf2_im
    endif
  enddo
  call this%solver_bm%solve(src = this%buf5, u = this%buf3)
  call this%solver_bp%solve(src = this%buf6, u = this%buf3)
  do i = 0, this%max_mode
    if ( i == 0 ) then
      this%buf1_re(:) = this%buf5(1+this%max_mode,:,1)
      this%buf2_re(:) = this%buf6(1+this%max_mode,:,1)
      call this%get_solution_bt_iter(i)
    else
      this%buf1_re(:) = this%buf5(1+this%max_mode + i,:,1)
      this%buf1_im(:) = this%buf5(1+this%max_mode + i,:,2)
      this%buf2_re(:) = this%buf6(1+this%max_mode + i,:,1)
      this%buf2_im(:) = this%buf6(1+this%max_mode + i,:,2)
      call this%get_solution_bt_iter(i)
    endif
  enddo
!   endif

  call this%copy_gc_f1()

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine solve_field_bt_iter

end module field_b_class
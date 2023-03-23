module field_em_class

use field_class
use field_psi_class
use field_src_class
use field_solver_class
use ufield_class
use param
use sysutil_module
use parallel_module
use options_class
use debug_tool
use mpi

implicit none

private

character(len=20), parameter :: cls_name = "field_em"
integer, parameter :: cls_level = 3

public :: field_b, field_e

type, extends( field ) :: field_b

  ! private

  class( field_solver ), dimension(:), pointer :: solver_bz      => null()
  class( field_solver ), dimension(:), pointer :: solver_bt      => null()
  class( field_solver ), dimension(:), pointer :: solver_bphi    => null()
  class( field_solver ), dimension(:), pointer :: solver_bplus   => null()
  class( field_solver ), dimension(:), pointer :: solver_bminus  => null()
  type( options ) :: opts

  real, dimension(:), pointer :: buf1_re => null(), buf1_im => null()
  real, dimension(:), pointer :: buf2_re => null(), buf2_im => null()
  real, dimension(:), pointer :: buf    => null()
  integer :: part_shape, boundary

  contains

  generic :: solve => solve_field_bz, solve_field_bt, solve_field_bt_iter,solve_field_bphi
  generic :: new   => init_field_b

  procedure :: init_field_b
  procedure :: del => end_field_b
  procedure :: alloc => alloc_field_b
  procedure, private :: set_source_bz
  procedure, private :: set_source_bt
  procedure, private :: set_source_bphi
  procedure, private :: set_source_bt_iter
  procedure, private :: get_solution_bz
  procedure, private :: get_solution_bt
  procedure, private :: get_solution_bphi
  procedure, private :: get_solution_bt_iter
  procedure, private :: solve_field_bz
  procedure, private :: solve_field_bt
  procedure, private :: solve_field_bphi
  procedure, private :: solve_field_bt_iter

end type field_b

type, extends( field ) :: field_e

  ! private

  class( field_solver ), dimension(:), pointer :: solver_ez => null()
  real, dimension(:), pointer :: buf_re => null(), buf_im => null()

  contains

  generic :: solve => solve_field_ez_fast, solve_field_ez, &
                      solve_field_et, solve_field_et_beam
  generic :: new => init_field_e

  procedure :: init_field_e
  procedure :: alloc => alloc_field_e
  procedure :: del => end_field_e
  procedure, private :: set_source_ez
  procedure, private :: get_solution_ez
  procedure, private :: solve_field_ez
  procedure, private :: solve_field_ez_fast ! to be deleted in future
  procedure, private :: solve_field_et
  procedure, private :: solve_field_et_beam

end type field_e

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

  if ( .not. associated( this%solver_bphi ) ) then
    allocate( field_solver :: this%solver_bphi(0:max_mode) )
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
  this%opts = opts
  this%max_mode = max_mode
  this%part_shape = part_shape
  this%boundary = boundary

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
!       call this%solver_bphi(i)%new( opts, i, dr, kind=p_fk_bphi, &
!         bnd=boundary, stype=p_hypre_cycred )
    enddo

    allocate( this%buf1_re(nrp), this%buf1_im(nrp) )
    allocate( this%buf2_re(nrp), this%buf2_im(nrp) )

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
      call this%solver_bphi(i)%del()
    enddo
    deallocate( this%solver_bz )
    deallocate( this%solver_bplus )
    deallocate( this%solver_bminus )
    deallocate( this%solver_bphi )

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

subroutine set_source_bphi( this, mode, ef, psi, rho_re, rho_im, &
  cu_re, amu_re, gam_re, ve_re, cu_im, amu_im, gam_im, ve_im )

  implicit none

  class( field_b ), intent(inout) :: this
  class( field_e ), intent(in) :: ef
  class( field_psi ), intent(in) :: psi
  class( ufield ), intent(in) :: cu_re, amu_re, gam_re, rho_re,ve_re
  class( ufield ), intent(in), optional :: cu_im, amu_im, gam_im, rho_im, ve_im
  integer, intent(in) :: mode

  integer :: i, nrp, nvp, idproc, noff
  real, dimension(:,:), pointer :: f1_re => null(), f1_im => null()
  real, dimension(:,:), pointer :: f2_re => null(), f2_im => null()
  real, dimension(:,:), pointer :: f3_re => null(), f3_im => null()
  real, dimension(:,:), pointer :: f4_re => null(), f4_im => null()
  real, dimension(:,:), pointer :: f5_re => null(), f5_im => null()
  real, dimension(:,:), pointer :: f6_re => null(), f6_im => null()
  real, dimension(:,:), pointer :: f7_re => null(), f7_im => null()
  real :: idrh, idr, ir
  real :: n, ipsi
  real :: s1_re, s1_im, s2_re, s2_im, s3_re, s3_im, s4_re, s4_im, s5_re, s5_im, s6_re, s6_im
  character(len=20), save :: sname = 'set_source_bphi'

  call write_dbg( cls_name, sname, cls_level, 'starts' )
  call start_tprof( 'solve plasma bphi' )

  nvp    = num_procs_loc()
  idproc = id_proc_loc()
  nrp    = cu_re%get_ndp(1)
  noff   = cu_re%get_noff(1)
  idr    = 1.0 / this%dr
  idrh   = 0.5 * idr

  f1_re => cu_re%get_f1()
  f2_re => amu_re%get_f1()
  f3_re => gam_re%get_f1()
  f4_re => ef%rf_re(mode)%get_f1()
  f5_re => psi%rf_re(mode)%get_f1()
  f6_re => rho_re%get_f1()
  f7_re => ve_re%get_f1()
  this%buf1_re = 0.0
  this%buf2_re = 0.0


  if ( present(cu_im) .and. present(amu_im) .and. present(gam_im) ) then
    f1_im => cu_im%get_f1()
    f2_im => amu_im%get_f1()
    f3_im => gam_im%get_f1()
    f4_im => ef%rf_im(mode)%get_f1()
    f5_im => psi%rf_im(mode)%get_f1()
    f6_im => rho_im%get_f1()
    f7_im => ve_im%get_f1()
    this%buf1_im = 0.0
    this%buf2_im = 0.0
  endif

  if ( mode == 0 ) then

    do i = 2, nrp
      ir = idr / real(i+noff-1)
      n = 1 - f6_re(1,i) 
      ipsi = 1 / (1 + f5_re(1,i))
      this%buf1_re(i) = idrh * ( f1_re(3,i+1) - f1_re(3,i-1) ) + n * idrh * ( f2_re(1,i+1) - f2_re(1,i-1) ) &
                        + n * ir * f2_re(1,i) + n * ipsi * f7_re(1,i) * f4_re(3,i) + n * ipsi * ( f2_re(1,i) - &
                        f3_re(1,i)*ipsi) * idrh * ( f5_re(1,i+1) - f5_re(1,i - 1) ) 
      this%buf2_re(i) = n*ipsi
!       this%buf1_re(i) = 0.0
!       this%buf2_re(i) = 0.0
    enddo

    ! calculate the derivatives at the boundary and axis
    if ( idproc == 0 ) then
      this%buf1_re(1) = 0.0
      this%buf2_re(1) = 0.0
      ! since Jz(m=0) is multiplied by factor 8 on axis, the derivative on index=2 is
      ! calculated using forward difference
      ir = idr / real(2+noff-1)
      n = 1 - f6_re(1,2) 
      ipsi = 1 / (1 + f5_re(1,2))
      this%buf1_re(2) =   idr * ( f1_re(3,3) - f1_re(3,2) ) + n * idr * ( f2_re(1,3) - f2_re(1,2) ) &
                         + n * ir * f2_re(1,2) + n * ipsi * f7_re(1,2) * f4_re(3,2) + n * ipsi * ( f2_re(1,2) - &
                         f3_re(1,2) * ipsi ) * idr * ( f5_re(1,3) - f5_re(1,2) ) 
      this%buf2_re(2) = n*ipsi
!       this%buf1_re(2) = 0.0
!       this%buf2_re(2) = 0.0
    else
      ir = idr / real(1+noff-1)
      n = 1 - f6_re(1,1) 
      ipsi = 1 / (1 + f5_re(1,1))
      this%buf1_re(1) =  idrh * ( f1_re(3,2)-f1_re(3,0) ) + n * idrh * ( f2_re(1,2) - f2_re(1,0) ) &
                        + n * ir * f2_re(1,1) + n * ipsi * f7_re(1,1) * f4_re(3,1) + n * ipsi * ( f2_re(1,1) - &
                        f3_re(1,1)*ipsi) * idrh * ( f5_re(1,2) - f5_re(1,0) ) 
      this%buf2_re(1) = n*ipsi
!       this%buf1_re(1) = 0.0
!       this%buf2_re(1) = 0.0
    endif

    if ( idproc == nvp-1 ) then
      ir = idr / real(nrp+noff-1)
      n = 1 - f6_re(1,nrp) 
      ipsi = 1 / (1 + f5_re(1,nrp))
      this%buf1_re(nrp) = idrh * ( 3.0 * f1_re(3,nrp) - 4.0 * f1_re(3,nrp-1) + f1_re(3,nrp-2) ) + &
                          n * idrh * ( 3.0 * f2_re(1,nrp) - 4.0 * f2_re(1,nrp-1) + f2_re(1,nrp-2) ) &
                          + n * ir * f2_re(1,nrp) + n * ipsi * f7_re(1,nrp) * f4_re(3,nrp) + n * ipsi * ( f2_re(1,nrp) &
                          + f3_re(1,nrp) * ipsi) * idrh * ( 3.0 * f5_re(1,nrp) - 4.0 * f5_re(1,nrp - 1) + &
                          f5_re(1,nrp-2)) 
      this%buf2_re(nrp) = n*ipsi
!       this%buf1_re(nrp) = 0.0
!       this%buf2_re(nrp) = 0.0
    endif

  else

    call write_err( 'Invalid input arguments!' )

  endif

  call stop_tprof( 'solve plasma bphi' )
  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine set_source_bphi

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

  call write_dbg( cls_name, sname, cls_level, 'starts' )
  call start_tprof( 'solve plasma bt' )

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
      this%buf1_re(i) = -f1_re(2,i) - f3_re(1,i) ! Re(Br)
      this%buf2_re(i) = f1_re(1,i) + idrh * ( f2_re(3,i+1) - f2_re(3,i-1) ) - f3_re(2,i) ! Re(Bphi)
    enddo

    ! calculate the derivatives at the boundary and axis
    if ( idproc == 0 ) then
      this%buf1_re(1) = 0.0
      this%buf2_re(1) = 0.0
      ! since Jz(m=0) is multiplied by factor 8 on axis, the derivative on index=2 is
      ! calculated using forward difference
      this%buf2_re(2) =  f1_re(1,2) + idr * ( f2_re(3,3) - f2_re(3,2) ) - f3_re(2,2)
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

subroutine get_solution_bphi( this, mode )

  implicit none

  class( field_b ), intent(inout) :: this
  integer, intent(in) :: mode

  integer :: i, nrp, idproc
  real, dimension(:,:), pointer :: f1_re => null(), f1_im => null()
  character(len=20), save :: sname = 'get_solution_bphi'

  call write_dbg( cls_name, sname, cls_level, 'starts' )
  call start_tprof( 'solve plasma bphi' )

  nrp    = this%rf_re(mode)%get_ndp(1)
  idproc = id_proc_loc()

  if ( mode == 0 ) then

    f1_re => this%rf_re(mode)%get_f1()
    do i = 1, nrp
      f1_re(2,i) = this%buf1_re(i) ! Re(Bphi)
    enddo

    if ( idproc == 0 ) then
      f1_re(2,1) = 0.0
    endif

  endif

  call stop_tprof( 'solve plasma bphi' )
  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine get_solution_bphi


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

subroutine solve_field_bphi( this, ef, psi, rho, cu, amu, gam, ve )

  implicit none

  class( field_b ), intent(inout) :: this
  class( field_e ), intent(inout) :: ef
  class( field_psi ), intent(inout) :: psi
  class( field_rho ), intent(inout) :: rho
  class( field_jay ), intent(inout) :: cu
  class( field_jay ), intent(inout) :: amu
  class( field_jay ), intent(inout) :: ve
  class( field_gam ), intent(inout) :: gam

  type( ufield ), dimension(:), pointer :: cu_re => null(), cu_im => null()
  type( ufield ), dimension(:), pointer :: ve_re => null(), ve_im => null()
  type( ufield ), dimension(:), pointer :: amu_re => null(), amu_im => null()
  type( ufield ), dimension(:), pointer :: gam_re => null(), gam_im => null()
  type( ufield ), dimension(:), pointer :: rho_re => null(), rho_im => null()
  type( ufield ), dimension(:), pointer :: ef_re => null(), ef_im => null()
  type( ufield ), dimension(:), pointer :: psi_re => null(), psi_im => null()

  integer :: i
  real :: dr
  character(len=20), save :: sname = 'solve_field_bphi'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  cu_re => cu%get_rf_re()
  cu_im => cu%get_rf_im()
  ve_re => ve%get_rf_re()
  ve_im => ve%get_rf_im()
  amu_re => amu%get_rf_re()
  amu_im => amu%get_rf_im()
  gam_re => gam%get_rf_re()
  gam_im => gam%get_rf_im()
  rho_re => rho%get_rf_re()
  rho_im => rho%get_rf_im()

  dr = this%opts%get_dr()

  do i = 0, this%max_mode

    if ( i == 0 ) then 
      call this%set_source_bphi( i, ef, psi, rho_re(i), rho_im(i), cu_re(i), amu_re(i), gam_re(i), ve_re(i))
      call this%solver_bphi(i)%init( this%opts, i, dr, kind=p_fk_bphi, bnd=this%boundary, stype=p_hypre_cycred, coef=this%buf2_re ) 
      call this%solver_bphi(i)%solve( this%buf1_re ) 
      call this%get_solution_bphi(i) 
      cycle
    endif

!     call this%set_source_bphi( i, cu_re(i), amu_re(i), gam_re(i), cu_im(i), amu_im(i), gam_im(i) )
!     call this%solver_bphi(i)%solve( this%buf1_re )
!     call this%solver_bphi(i)%solve( this%buf1_im )
!     call this%get_solution_bphi(i)

  enddo

  call this%copy_gc_f1()

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine solve_field_bphi

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

  do i = 0, this%max_mode

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

  call this%copy_gc_f1()

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine solve_field_bt_iter

subroutine alloc_field_e( this, max_mode )

  implicit none

  class( field_e ), intent(inout) :: this
  integer, intent(in) :: max_mode

  if ( .not. associated( this%solver_ez ) ) then
    allocate( field_solver :: this%solver_ez(0:max_mode) )
  endif

end subroutine alloc_field_e

subroutine init_field_e( this, opts, max_mode, part_shape, boundary, entity )

  implicit none

  class( field_e ), intent(inout) :: this
  type( options ), intent(in) :: opts
  integer, intent(in) :: max_mode, part_shape, entity, boundary

  integer, dimension(2,2) :: gc_num
  integer :: dim, i, nrp
  real :: dr
  character(len=20), save :: sname = "init_field_e"

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
      call this%solver_ez(i)%new( opts, i, dr, kind=p_fk_ez, &
        bnd=boundary, stype=p_hypre_cycred )
    enddo
    allocate( this%buf_re(nrp), this%buf_im(nrp) )
  case ( p_entity_beam )
    ! do nothing
  case default
    call write_err( 'Invalid field entity type.' )
  end select

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine init_field_e

subroutine end_field_e( this )

  implicit none

  class( field_e ), intent(inout) :: this

  integer :: i
  character(len=20), save :: sname = 'end_field_e'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  if ( this%entity == p_entity_plasma ) then
    do i = 0, this%max_mode
      call this%solver_ez(i)%del()
    enddo
    deallocate( this%solver_ez )
  endif

  if ( associated( this%buf_re ) ) deallocate( this%buf_re )
  if ( associated( this%buf_im ) ) deallocate( this%buf_im )

  call this%field%del()

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine end_field_e

subroutine set_source_ez( this, mode, jay_re, jay_im )

  implicit none

  class( field_e ), intent(inout) :: this
  class( ufield ), intent(in) :: jay_re
  class( ufield ), intent(in), optional :: jay_im
  integer, intent(in) :: mode

  integer :: i, nrp, noff, idproc, nvp, ierr, i1, i2
  real, dimension(:,:), pointer :: f1_re => null(), f1_im => null()
  real :: idrh, idr, dr, dr2, ir, div
  character(len=20), save :: sname = 'set_source_ez'

  call write_dbg( cls_name, sname, cls_level, 'starts' )
  call start_tprof( 'solve ez' )

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

  this%buf_re = 0.0
  if ( present(jay_im) ) this%buf_im = 0.0

  if ( mode == 0 ) then

    i1 = 1
    i2 = nrp
    if ( idproc == 0 ) i1 = 2
    if ( idproc == nvp-1 ) i2 = nrp-2

    div = 0.0
    do i = i1, i2
      ir = idr / real(i+noff-1)
      this%buf_re(i) = idrh * ( f1_re(1,i+1) - f1_re(1,i-1) ) + ir * f1_re(1,i)
      div = div + this%buf_re(i) * real(i+noff-1)
    enddo

    if ( idproc == nvp-1 ) then
      ir = idr / real(nrp+noff-2)
      this%buf_re(nrp-1) = idrh * ( f1_re(1,nrp) - f1_re(1,nrp-2) ) + ir * f1_re(1,nrp-1)
      ir = idr / real(nrp+noff-1)
      this%buf_re(nrp) = idr * ( f1_re(1,nrp) - f1_re(1,nrp-1) ) + ir * f1_re(1,nrp)
      div = div - idrh * (f1_re(1,nrp-2) + f1_re(1,nrp-1)) * (real(nrp+noff) - 2.5)
    ! else
    !   ir = idr / real(nrp+noff-1)
    !   this%buf_re(nrp) = idrh * ( f1_re(1,nrp+1) - f1_re(1,nrp-1) ) + ir * f1_re(1,nrp)
    !   div = div + this%buf_re(nrp) * real(nrp+noff-1)
    endif

    call MPI_REDUCE( div, this%buf_re(1), 1, p_dtype_real, MPI_SUM, 0, comm_loc(), ierr )
    if ( idproc == 0 ) this%buf_re(1) = -8.0 * this%buf_re(1)

  elseif ( mode > 0 .and. present( jay_im ) ) then

    do i = 2, nrp
      ir = idr / real(i+noff-1)
      this%buf_re(i) = idrh * ( f1_re(1,i+1) - f1_re(1,i-1) ) + ir * f1_re(1,i) - mode * ir * f1_im(2,i)
      this%buf_im(i) = idrh * ( f1_im(1,i+1) - f1_im(1,i-1) ) + ir * f1_im(1,i) + mode * ir * f1_re(2,i)
    enddo

    ! calculate the derivatives at the boundary and axis
    if ( idproc == 0 ) then
      if ( mod(mode,2) == 0 ) then
        this%buf_re(1) = 2.0 * idr * f1_re(1,2) - mode * idr * f1_im(2,2)
        this%buf_im(1) = 2.0 * idr * f1_im(1,2) + mode * idr * f1_re(2,2)
      else
        this%buf_re(1) = 0.0
        this%buf_im(1) = 0.0

        ! since j_r(m=1) is multiplied by factor 8 on axis, the derivative on index=2 is
        ! calculated using forward difference
        if ( mode == 1 ) then
          ir = idr
          this%buf_re(2) = idr * ( f1_re(1,3) - f1_re(1,2) ) + ir * f1_re(1,2) - mode * ir * f1_im(2,2)
          this%buf_im(2) = idr * ( f1_im(1,3) - f1_im(1,2) ) + ir * f1_im(1,2) + mode * ir * f1_re(2,2)
        endif
      endif
    else
      ir = idr / real(noff)
      this%buf_re(1) = idrh * ( f1_re(1,2) - f1_re(1,0) ) + ir * f1_re(1,1) - mode * ir * f1_im(2,1)
      this%buf_im(1) = idrh * ( f1_im(1,2) - f1_im(1,0) ) + ir * f1_im(1,1) + mode * ir * f1_re(2,1)
    endif

    if ( idproc == nvp-1 ) then
      ir = idr / real(nrp+noff-1)
      this%buf_re(nrp) = idrh * ( 3.0 * f1_re(1,nrp) - 4.0 * f1_re(1,nrp-1) + f1_re(1,nrp-2) ) + ir * f1_re(1,nrp) &
                         -mode * ir * f1_im(2,nrp)
      this%buf_im(nrp) = idrh * ( 3.0 * f1_im(1,nrp) - 4.0 * f1_im(1,nrp-1) + f1_im(1,nrp-2) ) + ir * f1_im(1,nrp) &
                         +mode * ir * f1_re(2,nrp)
    endif

  else

    call write_err( 'Invalid input arguments!' )

  endif

  call stop_tprof( 'solve ez' )
  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine set_source_ez

subroutine get_solution_ez( this, mode )

  implicit none

  class( field_e ), intent(inout) :: this
  integer, intent(in) :: mode

  integer :: i, nrp
  real, dimension(:,:), pointer :: f1_re => null(), f1_im => null()
  character(len=20), save :: sname = 'get_solution_ez'

  call write_dbg( cls_name, sname, cls_level, 'starts' )
  call start_tprof( 'solve ez' )

  nrp = this%rf_re(mode)%get_ndp(1)

  f1_re => this%rf_re(mode)%get_f1()
  do i = 1, nrp
    f1_re(3,i) = this%buf_re(i)
  enddo

  if ( mode > 0 ) then
    f1_im => this%rf_im(mode)%get_f1()
    do i = 1, nrp
      f1_im(3,i) = this%buf_im(i)
    enddo

    ! Ez(m>0) vanishes on axis
    f1_re(3,1) = 0.0
    f1_im(3,1) = 0.0
  endif

  call stop_tprof( 'solve ez' )
  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine get_solution_ez

subroutine solve_field_ez( this, jay )

  implicit none

  class( field_e ), intent(inout) :: this
  class( field_jay ), intent(inout) :: jay

  type( ufield ), dimension(:), pointer :: jay_re => null(), jay_im => null()
  integer :: i
  character(len=20), save :: sname = 'solve_field_ez'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  jay_re => jay%get_rf_re()
  jay_im => jay%get_rf_im()

  do i = 0, this%max_mode

    if ( i == 0 ) then
      call this%set_source_ez( i, jay_re(i) )
      call this%solver_ez(i)%solve( this%buf_re )
      call this%get_solution_ez(i) 
    endif 
!       cycle
    

!     call this%set_source_ez( i, jay_re(i), jay_im(i) )
!     call this%solver_ez(i)%solve( this%buf_re )
!     call this%solver_ez(i)%solve( this%buf_im )
!     call this%get_solution_ez(i)

  enddo

  call this%copy_gc_f1()

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine solve_field_ez

subroutine solve_field_ez_fast( this, psi, idx )

  implicit none

  class( field_e ), intent(inout) :: this
  class( field_psi ), intent(in) :: psi
  integer, intent(in) :: idx

  type( ufield ), dimension(:), pointer :: psi_re => null(), psi_im => null()
  type( ufield ), dimension(:), pointer :: e_re => null(), e_im => null()
  real, dimension(:,:), pointer :: e_f1_re => null(), e_f1_im => null()
  real, dimension(:,:,:), pointer :: e_f2_re => null(), e_f2_im => null()
  real, dimension(:,:), pointer :: psi_f1_re => null(), psi_f1_im => null()
  real, dimension(:,:,:), pointer :: psi_f2_re => null(), psi_f2_im => null()
  integer :: i, j, nrp
  real :: idxih, idxi
  character(len=20), save :: sname = 'solve_field_ez_fast'

  call write_dbg( cls_name, sname, cls_level, 'starts' )
  call start_tprof( 'solve ez' )

  psi_re => psi%get_rf_re()
  psi_im => psi%get_rf_im()
  e_re => this%get_rf_re()
  e_im => this%get_rf_im()

  nrp = e_re(0)%get_ndp(1)
  idxih = 0.5 / this%dxi
  idxi = 1.0 / this%dxi

  do i = 0, this%max_mode

    e_f1_re => e_re(i)%get_f1()
    e_f2_re => e_re(i)%get_f2()
    psi_f1_re => psi_re(i)%get_f1()
    psi_f2_re => psi_re(i)%get_f2()

    do j = 1, nrp
      ! quadratic extrapolation
      e_f1_re(3,j) = idxih * ( 3.0 * psi_f1_re(1,j) - 4.0 * psi_f2_re(1,j,idx-1) + psi_f2_re(1,j,idx-2) )

      ! linear extrapolation
      ! e_f1_re(3,j) = idxi * ( psi_f1_re(1,j) - psi_f2_re(1,j,idx-1) )

      ! smooth
      ! e_f1_re(3,j) = ( e_f1_re(3,j) + e_f2_re(3,j,idx-1) + e_f2_re(3,j,idx-2) ) / 3
    enddo

    if ( i == 0 ) cycle

    e_f1_im => e_im(i)%get_f1()
    e_f2_im => e_im(i)%get_f2()
    psi_f1_im => psi_im(i)%get_f1()
    psi_f2_im => psi_im(i)%get_f2()

    do j = 1, nrp
      ! quadratic extrapolation
      e_f1_im(3,j) = idxih * ( 3.0 * psi_f1_im(1,j) - 4.0 * psi_f2_im(1,j,idx-1) + psi_f2_im(1,j,idx-2) )

      ! linear extrapolation
      ! e_f1_im(3,j) = idxi * ( psi_f1_im(1,j) - psi_f2_im(1,j,idx-1) )

      ! smooth
      ! e_f1_im(3,j) = ( e_f1_im(3,j) + e_f2_im(3,j,idx-1) + e_f2_im(3,j,idx-2) ) / 3
    enddo

  enddo

  call stop_tprof( 'solve ez' )

  call this%copy_gc_f1()

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine solve_field_ez_fast

subroutine solve_field_et( this, b, psi )

  implicit none

  class( field_e ), intent(inout) :: this
  class( field_b ), intent(in) :: b
  class( field_psi ), intent(inout) :: psi

  type( ufield ), dimension(:), pointer :: b_re => null(), b_im => null()
  type( ufield ), dimension(:), pointer :: psi_re => null(), psi_im => null()
  real, dimension(:,:), pointer :: ub_re => null(), ub_im => null()
  real, dimension(:,:), pointer :: upsi_re => null(), upsi_im => null()
  real, dimension(:,:), pointer :: ue_re => null(), ue_im => null()
  integer :: mode, i, nrp, noff, idproc, nvp
  real :: idr, idrh, ir
  character(len=20), save :: sname = 'solve_field_et'

  call write_dbg( cls_name, sname, cls_level, 'starts' )
  call start_tprof( 'solve plasma et' )

  idr = 1.0 / this%dr
  idrh = idr * 0.5
  nrp = this%rf_re(0)%get_ndp(1)

  b_re   => b%get_rf_re()
  b_im   => b%get_rf_im()
  psi_re => psi%get_rf_re()
  psi_im => psi%get_rf_im()

  noff   = this%rf_re(0)%get_noff(1)
  nvp    = num_procs_loc()
  idproc = id_proc_loc()

  ! m=0 mode
  ub_re   => b_re(0)%get_f1()
  upsi_re => psi_re(0)%get_f1()
  ue_re   => this%rf_re(0)%get_f1()

  do i = 1, nrp
    ue_re(1,i) = ub_re(2,i) - idrh * ( upsi_re(1,i+1) - upsi_re(1,i-1) )
    ue_re(2,i) = -ub_re(1,i)
  enddo
  if ( idproc == 0 ) then
    ue_re(1,1) = 0.0
    ue_re(2,1) = 0.0
  endif
  if ( idproc == nvp-1 ) then
    ue_re(1,nrp) = ub_re(2,nrp) + idrh * ( 4.0 * upsi_re(1,nrp-1) - upsi_re(1,nrp-2) - 3.0 * upsi_re(1,nrp) )
  endif

  ! m>0 mode
  do mode = 1, this%max_mode

    ub_re   => b_re(mode)%get_f1()
    ub_im   => b_im(mode)%get_f1()
    upsi_re => psi_re(mode)%get_f1()
    upsi_im => psi_im(mode)%get_f1()
    ue_re   => this%rf_re(mode)%get_f1()
    ue_im   => this%rf_im(mode)%get_f1()

    do i = 2, nrp
      ir = idr / real(i+noff-1)
      ue_re(1,i) = ub_re(2,i) - idrh * ( upsi_re(1,i+1) - upsi_re(1,i-1) )
      ue_im(1,i) = ub_im(2,i) - idrh * ( upsi_im(1,i+1) - upsi_im(1,i-1) )
      ue_re(2,i) = -ub_re(1,i) + ir * mode * upsi_im(1,i)
      ue_im(2,i) = -ub_im(1,i) - ir * mode * upsi_re(1,i)
    enddo

    if ( idproc == 0 ) then
      if ( mode == 1 ) then
        ue_re(1,1) = ub_re(2,1) - idr * upsi_re(1,2)
        ue_im(1,1) = ub_im(2,1) - idr * upsi_im(1,2)
        ue_re(2,1) = -ub_re(1,1) + idr * upsi_im(1,2)
        ue_im(2,1) = -ub_im(1,1) - idr * upsi_re(1,2)
      else
        ue_re(1,1) = 0.0
        ue_im(1,1) = 0.0
        ue_re(2,1) = 0.0
        ue_im(2,1) = 0.0
      endif
    else
      ir = idr / real(noff)
      ue_re(1,1) = ub_re(2,1) - idrh * ( upsi_re(1,2) - upsi_re(1,0) )
      ue_im(1,1) = ub_im(2,1) - idrh * ( upsi_im(1,2) - upsi_im(1,0) )
      ue_re(2,1) = -ub_re(1,1) + ir * mode * upsi_im(1,1)
      ue_im(2,1) = -ub_im(1,1) - ir * mode * upsi_re(1,1)
    endif

    if ( idproc == nvp-1 ) then
      ir = idr / real(nrp+noff-1)
      ue_re(1,nrp) = ub_re(2,nrp) + idrh * ( 4.0 * upsi_re(1,nrp-1) - upsi_re(1,nrp-2) - 3.0 * upsi_re(1,nrp) )
      ue_im(1,nrp) = ub_im(2,nrp) + idrh * ( 4.0 * upsi_im(1,nrp-1) - upsi_im(1,nrp-2) - 3.0 * upsi_im(1,nrp) )
    endif

  enddo

  call stop_tprof( 'solve plasma et' )

  call this%copy_gc_f1()

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine solve_field_et

subroutine solve_field_et_beam( this, b )

  implicit none

  class( field_e ), intent(inout) :: this
  class( field_b ), intent(in) :: b

  type( ufield ), dimension(:), pointer :: b_re => null(), b_im => null()
  real, dimension(:,:), pointer :: ub_re => null(), ub_im => null()
  real, dimension(:,:), pointer :: ue_re => null(), ue_im => null()
  integer :: mode, i, nrp
  character(len=20), save :: sname = 'solve_field_et_beam'

  call write_dbg( cls_name, sname, cls_level, 'starts' )
  call start_tprof( 'solve beam et' )

  b_re => b%get_rf_re()
  b_im => b%get_rf_im()
  nrp = this%rf_re(0)%get_ndp(1)

  do mode = 0, this%max_mode

    ub_re => b_re(mode)%get_f1()
    ue_re => this%rf_re(mode)%get_f1()

    do i = 1, nrp
      ue_re(1,i) = ub_re(2,i)
      ue_re(2,i) = -ub_re(1,i)
    enddo

    if ( mode == 0 ) cycle

    ub_im => b_im(mode)%get_f1()
    ue_im => this%rf_im(mode)%get_f1()

    do i = 1, nrp
      ue_im(1,i) = ub_im(2,i)
      ue_im(2,i) = -ub_im(1,i)
    enddo

  enddo

  call stop_tprof( 'solve beam et' )

  call this%copy_gc_f1()

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine solve_field_et_beam

end module field_em_class
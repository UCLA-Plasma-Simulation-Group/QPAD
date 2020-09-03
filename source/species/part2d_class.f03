! part2d class for QPAD

module part2d_class

use param
use sysutil
use parallel_pipe_class
use grid_class
use field_class
use ufield_class
use fdist2d_class
use hdf5io_class
use part2d_lib
use mpi
use interpolation

implicit none

private

public :: part2d

type part2d

   class(parallel_pipe), pointer :: pp => null()

   ! qbm = particle charge/mass ratio
   ! dt = time interval between successive calculations
   ! dr = radial cell size
   real :: qbm, dt, dr

   ! nbmax = size of buffer for passing particles between processors
   ! npp = number of particles in current partition
   ! npmax = maximum number of particles in each partition
   integer(kind=LG) :: npmax, nbmax, npp

   ! dimension of particle coordinates
   integer :: part_dim

   ! array for particle position
   real, dimension(:,:), pointer :: x => null()
   ! array for particle momenta
   real, dimension(:,:), pointer :: p => null()
   ! array for time-centered gamma
   real, dimension(:), pointer :: gamma => null()
   ! array for particle charge
   real, dimension(:), pointer :: q => null()
   ! array for psi
   real, dimension(:), pointer :: psi => null()

   ! particle upper boundaries
   real :: edge

   contains

   procedure :: new         => init_part2d
   procedure :: renew       => renew_part2d
   procedure :: del         => end_part2d
   procedure :: qdeposit    => qdeposit_part2d
   procedure :: amjdeposit  => amjdeposit_part2d
   procedure :: push        => push_part2d
   procedure :: extract_psi => extract_psi_part2d
   procedure :: pipesend    => pipesend_part2d
   procedure :: piperecv    => piperecv_part2d

   ! TODO: particle manager to be rewritten
   procedure :: pmv => pmove

   procedure :: wr => writehdf5_part2d

end type

save

character(len=20), parameter :: cls_name = "part2d"
integer, parameter :: cls_level = 2

! TODO: data communication, to be deleted
real, dimension(:,:), allocatable :: sbufl, sbufr, rbufl, rbufr
integer(kind=LG), dimension(:), allocatable :: ihole

contains

subroutine init_part2d( this, pp, gp, pf, qbm, dt, s )

   implicit none

   class(part2d), intent(inout) :: this
   class(parallel_pipe), intent(in), pointer :: pp
   class(grid), intent(in), pointer :: gp
   class(fdist2d), intent(inout) :: pf
   real, intent(in) :: qbm, dt, s

   ! local data
   character(len=18), save :: sname = 'init_part2d'
   integer :: npmax, nbmax

   call write_dbg( cls_name, sname, cls_level, 'starts' )

   this%pp  => pp
   this%qbm = qbm
   this%dt  = dt
   this%part_dim = 2 + p_p_dim + 3

   npmax      = pf%getnpmax()
   nbmax      = max(int(0.01*npmax),100)
   this%npmax = npmax
   this%nbmax = nbmax
   this%npp   = 0

   this%dr   = pf%getdex()
   this%edge = gp%get_nd(1) * this%dr
   
   allocate( this%x( 2, npmax ) )
   allocate( this%p( p_p_dim, npmax ) )
   allocate( this%gamma( npmax ), this%q( npmax ), this%psi( npmax ) )

   ! initialize particle coordinates according to specified profile
   call pf%dist( this%x, this%p, this%gamma, this%q, this%psi, this%npp, s )
   if (.not. allocated(sbufl)) then
      allocate( sbufl( this%part_dim, nbmax ), sbufr( this%part_dim, nbmax ) )
      allocate( rbufl( this%part_dim, nbmax ), rbufr( this%part_dim, nbmax ) )
      allocate( ihole( nbmax*2 ) )
   end if

   call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine init_part2d

subroutine end_part2d(this)

   implicit none

   class(part2d), intent(inout) :: this
   character(len=18), save :: sname = 'end_part2d'

   call write_dbg(cls_name, sname, cls_level, 'starts')

   deallocate( this%x, this%p, this%gamma, this%q, this%psi )

   call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine end_part2d

subroutine renew_part2d( this, pf, s )

   implicit none

   class(part2d), intent(inout) :: this
   class(fdist2d), intent(inout) :: pf
   real, intent(in) :: s
   ! local data
   character(len=18), save :: sname = 'renew_part2d'

   call write_dbg(cls_name, sname, cls_level, 'starts')

   call pf%dist( this%x, this%p, this%gamma, this%q, this%psi, this%npp, s )

   call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine renew_part2d

subroutine qdeposit_part2d( this, q )
! deposit the charge density (rho - Jz)

  implicit none

  class(part2d), intent(in) :: this
  class(field), intent(in) :: q
  ! local data
  class(ufield), dimension(:), pointer :: q_re => null(), q_im => null()
  complex(kind=DB), dimension(p_cache_size) :: phase0
  complex(kind=DB) :: phase
  real, dimension(:,:), pointer :: q0 => null(), qr => null(), qi => null()
  real, dimension(0:1) :: wt ! interpolation weight
  real, dimension(p_cache_size) :: pos ! normalized position
  real :: idr, r
  integer(kind=LG) :: ptrcur, pp
  integer :: i, j, k, nn, noff, nrp, np, mode, max_mode

  character(len=18), save :: sname = 'qdeposit_part2d'

  call write_dbg(cls_name, sname, cls_level, 'starts')
  call start_tprof( 'deposit 2D particles' )

  idr = 1.0 / this%dr

  q_re => q%get_rf_re()
  q_im => q%get_rf_im()

  max_mode = q%get_num_modes()

  noff = q_re(0)%get_noff(1)
  nrp  = q_re(0)%get_ndp(1)
  q0   => q_re(0)%get_f1()

  do ptrcur = 1, this%npp, p_cache_size

    ! check if last copy of table and set np
    if( ptrcur + p_cache_size > this%npp ) then
      np = this%npp - ptrcur + 1
    else
      np = p_cache_size
    endif

    pp = ptrcur
    do i = 1, np
       ! get normalized position
       pos(i) = sqrt( this%x(1, pp)**2 + this%x(2, pp)**2 ) * idr
       phase0(i) = cmplx( this%x(1, pp), -this%x(2, pp) ) / pos(i) * idr
       pp = pp + 1
    enddo

    pp = ptrcur
    do i = 1, np
      nn = floor( pos(i) )
      ! in-cell position
      pos(i) = pos(i) - real(nn)
      nn = nn - noff + 1

      call spline_linear( pos(i), wt )

      phase = cmplx( 1.0, 0.0 ) * this%q(pp)
      ! deposit m=0 mode
      do j = 0, 1
        q0( 1, nn+j ) = q0( 1, nn+j ) + wt(j) * real(phase)
      enddo

      ! deposit m>0 mode
      do mode = 1, max_mode
        qr => q_re(mode)%get_f1()
        qi => q_im(mode)%get_f1()
        phase = phase * phase0(i)

        do j = 0, 1
          qr( 1, nn+j ) = qr( 1, nn+j ) + wt(j) * real(phase)
          qi( 1, nn+j ) = qi( 1, nn+j ) + wt(j) * aimag(phase)
        enddo
      enddo
       pp = pp + 1
    enddo

  enddo

  if (noff == 0) then

    q0(1,0) = 0.0 ! guard cell is useless on axis
    q0(1,1) = 8.0 * q0(1,1)
    do j = 2, nrp + 1
      r = j + noff - 1
      q0(1,j) = q0(1,j) / r
    enddo

    do i = 1, max_mode
      qr => q_re(i)%get_f1()
      qi => q_im(i)%get_f1()
      qr(1,0) = 0.0
      qi(1,0) = 0.0
      qr(1,1) = 0.0
      qi(1,1) = 0.0
      do j = 2, nrp + 1
        r = j + noff - 1
        qr(1,j) = qr(1,j) / r
        qi(1,j) = qi(1,j) / r
      enddo
    enddo

  else

    do j = 0, nrp + 1
      r = j + noff -1
      q0(1,j) = q0(1,j) / r
    enddo

    do i = 1, max_mode
      qr => q_re(i)%get_f1()
      qi => q_im(i)%get_f1()
      do j = 0, nrp + 1
        r = j + noff -1
        qr(1,j) = qr(1,j) / r
        qi(1,j) = qi(1,j) / r
      enddo
    enddo

  endif

  call stop_tprof( 'deposit 2D particles' )
  call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine qdeposit_part2d

subroutine amjdeposit_part2d(this,ef,bf,cu,amu,dcu)
! deposit the current, acceleration and momentum flux

   implicit none

   class(part2d), intent(inout) :: this
   class(field), intent(in) :: cu, amu, dcu
   class(field), intent(in) :: ef, bf
! local data
   character(len=18), save :: sname = 'amjdeposit'
   class(ufield), dimension(:), pointer :: ef_re => null(), ef_im => null()
   class(ufield), dimension(:), pointer :: bf_re => null(), bf_im => null()
   class(ufield), dimension(:), pointer :: cu_re => null(), cu_im => null()
   class(ufield), dimension(:), pointer :: dcu_re => null(), dcu_im => null()
   class(ufield), dimension(:), pointer :: amu_re => null(), amu_im => null()

   call write_dbg(cls_name, sname, cls_level, 'starts')

   ef_re => ef%get_rf_re()
   ef_im => ef%get_rf_im()
   bf_re => bf%get_rf_re()
   bf_im => bf%get_rf_im()
   cu_re => cu%get_rf_re()
   cu_im => cu%get_rf_im()
   dcu_re => dcu%get_rf_re()
   dcu_im => dcu%get_rf_im()
   amu_re => amu%get_rf_re()
   amu_im => amu%get_rf_im()

   ! call part2d_amjdeposit(this%part,this%npp,this%dr,this%dt,this%qbm,&
   ! &ef_re,ef_im,bf_re,bf_im,cu_re,cu_im,dcu_re,dcu_im,amu_re,amu_im,&
   ! &ef%get_num_modes())

   call part2d_amjdeposit(this%x, this%p, this%gamma, this%q, this%psi,&
   &this%npp,this%dr,this%dt,this%qbm,&
   &ef_re,ef_im,bf_re,bf_im,cu_re,cu_im,dcu_re,dcu_im,amu_re,amu_im,&
   &ef%get_num_modes())

   call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine amjdeposit_part2d
!
subroutine push_part2d(this,ef,bf)

   implicit none

   class(part2d), intent(inout) :: this
   class(field), intent(in) :: ef, bf
! local data
   character(len=18), save :: sname = 'partpush'
   class(ufield), dimension(:), pointer :: ef_re => null(), ef_im => null()
   class(ufield), dimension(:), pointer :: bf_re => null(), bf_im => null()

   call write_dbg(cls_name, sname, cls_level, 'starts')

   ef_re => ef%get_rf_re()
   ef_im => ef%get_rf_im()
   bf_re => bf%get_rf_re()
   bf_im => bf%get_rf_im()

   ! call part2d_push(this%part,this%npp,this%dr,this%xdim,this%dt,this%qbm,&
   ! &ef_re,ef_im,bf_re,bf_im,ef%get_num_modes())

   call part2d_push(this%x, this%p, this%gamma, this%q, this%psi,&
   &this%npp,this%dr,this%part_dim,this%dt,this%qbm,&
   &ef_re,ef_im,bf_re,bf_im,ef%get_num_modes())

   call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine push_part2d
!
subroutine pmove(this,fd)

   implicit none

   class(part2d), intent(inout) :: this
   class(field), intent(in) :: fd
! local data
   character(len=18), save :: sname = 'pmove'
   class(ufield), pointer :: ud

   call write_dbg(cls_name, sname, cls_level, 'starts')

   ud => fd%get_rf_re(0)
   ! call part2d_pmove(this%part,this%pp,this%npp,this%dr,this%xdim,this%npmax,&
   ! &this%nbmax,ud,sbufl,sbufr,rbufl,rbufr,ihole)

   call part2d_pmove(this%x, this%p, this%gamma, this%q, this%psi,&
   &this%pp,this%npp,this%dr,this%part_dim,this%npmax,&
   &this%nbmax,ud,sbufl,sbufr,rbufl,rbufr,ihole)

   call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine pmove
!
subroutine extract_psi_part2d(this,psi)

   implicit none

   class(part2d), intent(inout) :: this
   class(field), intent(in) :: psi
! local data
   character(len=18), save :: sname = 'extractpsi'
   class(ufield), dimension(:), pointer :: psi_re => null(), psi_im => null()

   call write_dbg(cls_name, sname, cls_level, 'starts')

   psi_re => psi%get_rf_re()
   psi_im => psi%get_rf_im()
   ! call part2d_extractpsi(this%part,this%npp,this%dr,this%qbm,psi_re,psi_im,psi%get_num_modes())
   call part2d_extractpsi(this%x, this%p, this%gamma, this%q, this%psi,&
    this%npp,this%dr,this%qbm,psi_re,psi_im,psi%get_num_modes())

   call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine extract_psi_part2d
!
subroutine pipesend_part2d(this, tag, id)

  implicit none

  class(part2d), intent(inout) :: this
  integer, intent(in) :: tag
  integer, intent(inout) :: id

  ! local data
  character(len=18), save :: sname = 'pipesend_part2d'
  integer :: des, ierr, i
  real, dimension(:,:), allocatable, save :: sbuf

  call write_dbg(cls_name, sname, cls_level, 'starts')

  if ( .not. allocated(sbuf) ) then
    allocate( sbuf( this%part_dim, this%npmax ) )
  endif

  des = this%pp%getidproc() + this%pp%getlnvp()

  if (des >= this%pp%getnvp()) then
    id = MPI_REQUEST_NULL
    call write_dbg(cls_name, sname, cls_level, 'ends')
    return
  endif

  ! to be implemented if using tile
  ! call this%pcb()

  do i = 1, this%npp
    sbuf(1:2,i) = this%x(:,i)
    sbuf(3:5,i) = this%p(:,i)
    sbuf(6,i)   = this%gamma(i)
    sbuf(7,i)   = this%psi(i)
    sbuf(8,i)   = this%q(i)
  enddo

  ! NOTE: npp*xdim might be larger than MAX_INT32
  call MPI_ISEND(sbuf, int(this%npp*this%part_dim), this%pp%getmreal(), &
    des, tag, this%pp%getlworld(), id, ierr)

  ! check for errors
  if (ierr /= 0) then
    call write_err('MPI_ISEND failed')
  endif

  call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine pipesend_part2d

subroutine piperecv_part2d(this, tag)

  implicit none

  class(part2d), intent(inout) :: this
  integer, intent(in) :: tag
  ! local data
  character(len=18), save :: sname = 'piperecv_part2d'
  integer, dimension(MPI_STATUS_SIZE) :: istat
  integer :: nps, des, ierr, i
  real, dimension(:,:), allocatable, save :: rbuf

  call write_dbg(cls_name, sname, cls_level, 'starts')

  if ( .not. allocated(rbuf) ) then
    allocate( rbuf( this%part_dim, this%npmax ) )
  endif

  des = this%pp%getidproc() - this%pp%getlnvp()

  if (des < 0) then
    call write_dbg(cls_name, sname, cls_level, 'ends')
    return
  endif

  ! NOTE: npp*xdim might be larger than MAX_INT32
  call MPI_RECV(rbuf, int(this%npmax*this%part_dim), this%pp%getmreal(), &
    des, tag, this%pp%getlworld(), istat, ierr)

  call MPI_GET_COUNT(istat, this%pp%getmreal(), nps, ierr)

  this%npp = nps/this%part_dim

  do i = 1, this%npp
    this%x(:,i)   = rbuf(1:2,i)
    this%p(:,i)   = rbuf(3:5,i)
    this%gamma(i) = rbuf(6,i)
    this%psi(i)   = rbuf(7,i)
    this%q(i)     = rbuf(8,i)
  enddo

  ! to be implemented if using tile
  ! call this%pcp(fd)

  ! check for errors
  if (ierr /= 0) then
    call write_err('MPI failed')
  endif

   call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine piperecv_part2d
!
subroutine writehdf5_part2d(this,file)

   implicit none

   class(part2d), intent(inout) :: this
   class(hdf5file), intent(in) :: file
! local data
   character(len=18), save :: sname = 'writehdf5_part2d'
   integer :: ierr

   call write_dbg(cls_name, sname, cls_level, 'starts')

   ! call pwpart(this%pp,file,this%part,this%npp,1,ierr)
   call pwpart(this%pp,file,this%x, this%p, this%q, this%npp,1,ierr)

   call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine writehdf5_part2d
!
end module part2d_class
! part3d class for QPAD

module part3d_class

use param
use sysutil
use parallel_pipe_class
use grid_class
use field_class
use ufield_class
use fdist3d_class
use hdf5io_class
use part3d_comm
use mpi
use interpolation

implicit none

private

public :: part3d

type part3d

   private

! qbm = particle charge/mass ratio
! dt = time interval between successive calculations

   class(parallel_pipe), pointer, public :: pp => null()
   real :: qbm, dt, dr, dz
   real :: z0
   ! nbmax = size of buffer for passing particles between processors
   ! npp = number of particles in current partition
   ! npmax = maximum number of particles in each partition
   integer(kind=LG) :: npmax, nbmax, npp
   integer :: part_dim ! dimension of particle coordinates
   ! particle buffer for particle manager
   real, dimension(:,:), pointer :: pbuff => null()
   ! array for particle position
   real, dimension(:,:), pointer :: x => null()
   ! array for particle momenta
   real, dimension(:,:), pointer :: p => null()
   ! array for particle charge
   real, dimension(:), pointer :: q => null()
   ! array for particle spin
   real, dimension(:,:), pointer :: s => null()
   ! particle upper boundaries
   real, dimension(2) :: edge

   logical :: has_spin = .false.

   ! anomalous magnet moment for spin push
   real :: anom_mag_moment

   contains

   procedure :: new          => init_part3d
   procedure :: del          => end_part3d
   procedure :: qdeposit     => qdeposit_part3d
   procedure :: push_boris   => push_boris_part3d
   procedure :: push_reduced => push_reduced_part3d
   procedure :: update_bound => update_bound_part3d
   procedure, private :: push_spin => push_spin_part3d

   procedure :: pmv  => pmove
   procedure :: wr   => writehdf5_part3d
   procedure :: wrst => writerst_part3d
   procedure :: rrst => readrst_part3d

   procedure :: getnpp

end type

save

character(len=20), parameter :: cls_name = "part3d"
integer, parameter :: cls_level = 2

! buffer data for particle managers
real, dimension(:,:), allocatable :: sbufl, sbufr, rbufl, rbufr
integer(kind=LG), dimension(:), allocatable :: ihole

contains
!
subroutine init_part3d(this,pp,gp,pf,qbm,dt,has_spin,amm)

   implicit none

   class(part3d), intent(inout) :: this
   class(parallel_pipe), intent(in), pointer :: pp
   class(grid), intent(in), pointer :: gp
   class(fdist3d), intent(inout) :: pf
   real, intent(in) :: qbm, dt
   logical, intent(in) :: has_spin
   real, intent(in) :: amm
! local data
   character(len=18), save :: sname = 'init_part3d'
   integer :: prof, npmax, nbmax

   call write_dbg(cls_name, sname, cls_level, 'starts')

   this%has_spin = has_spin
   this%anom_mag_moment = amm

   this%pp => pp
   this%qbm = qbm
   this%dt = dt
   this%part_dim = p_x_dim + p_p_dim + 1
   if ( this%has_spin ) this%part_dim = this%part_dim + p_s_dim
   npmax = pf%getnpmax()
   this%dr = pf%getdx()
   this%dz = pf%getdz()
   this%npmax = npmax

   this%edge(1) = gp%get_nd(1) * this%dr
   this%edge(2) = (gp%get_nd(2) - 1) * this%dz

   ! *TODO* nbmax needs to be dynamically changed, otherwise it has the risk to overflow
   nbmax = int(0.01*this%npmax)
   this%nbmax = nbmax
   this%npp = 0
   prof = pf%getnpf()
   this%z0 = pf%getz0()

   allocate( this%x( p_x_dim, npmax ) )
   allocate( this%p( p_p_dim, npmax ) )
   allocate( this%q( npmax ) )
   if ( this%has_spin ) then
      allocate( this%s( p_s_dim, npmax ) )
   endif

   allocate( this%pbuff( this%part_dim, nbmax ) )

   ! initialize particle coordinates according to profile
   if ( this%has_spin ) then
      call pf%dist( this%x, this%p, this%q, this%npp, gp%get_noff(), &
         gp%get_ndp(), this%s )
   else
      call pf%dist( this%x, this%p, this%q, this%npp, gp%get_noff(), &
         gp%get_ndp() )
   endif

   if ( .not. allocated(sbufl) ) then
      allocate( sbufl( this%part_dim, nbmax ), sbufr( this%part_dim, nbmax ) )
      allocate( rbufl( this%part_dim, nbmax ), rbufr( this%part_dim, nbmax ) )
      allocate( ihole( nbmax*2 ) )
   endif
   call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine init_part3d
!
subroutine end_part3d(this)

   implicit none

   class(part3d), intent(inout) :: this
   character(len=18), save :: sname = 'end_part3d'

   call write_dbg(cls_name, sname, cls_level, 'starts')

   deallocate( this%x, this%p, this%q, this%pbuff )
   if ( this%has_spin ) then
      deallocate( this%s )
   endif

   call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine end_part3d
!
subroutine qdeposit_part3d( this, q )

   implicit none

   class(part3d), intent(inout) :: this
   class(field), intent(inout) :: q
   ! local data
   class(ufield), dimension(:), pointer :: q_re => null(), q_im => null()
   complex(kind=DB), dimension(p_cache_size) :: phase0
   complex(kind=DB) :: phase
   real, dimension(:,:,:), pointer :: q0 => null(), qr => null(), qi => null()
   real, dimension(0:1) :: wtr, wtz ! interpolation weight
   real, dimension(p_cache_size) :: pos_r, pos_z ! normalized position
   real :: idr, idz, r
   integer(kind=LG) :: ptrcur, pp
   integer :: i, j, k, nn, mm, noff1, noff2, n1p, n2p, np, mode, max_mode
   character(len=32), save :: sname = "qdeposit_part3d"

   call write_dbg(cls_name, sname, cls_level, 'starts')
   call start_tprof( 'deposit 3D particles' )

   idr = 1.0 / this%dr
   idz = 1.0 / this%dz

   q_re => q%get_rf_re()
   q_im => q%get_rf_im()

   max_mode = q%get_num_modes()

   noff1 = q_re(0)%get_noff(1)
   noff2 = q_re(0)%get_noff(2)
   n1p   = q_re(0)%get_ndp(1)
   n2p   = q_re(0)%get_ndp(2)
   q0    => q_re(0)%get_f2()

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
         pos_r(i) = sqrt( this%x(1,pp)**2 + this%x(2,pp)**2 ) * idr
         pos_z(i) = this%x(3,pp) * idz
         phase0(i) = cmplx( this%x(1,pp), -this%x(2,pp) ) / pos_r(i) * idr
         pp = pp + 1
      enddo

      pp = ptrcur
      do i = 1, np
         nn = int( pos_r(i) )
         mm = int( pos_z(i) )

         ! in-cell position
         pos_r(i) = pos_r(i) - real(nn)
         pos_z(i) = pos_z(i) - real(mm)

         nn = nn - noff1 + 1
         mm = mm - noff2 + 1

         call spline_linear( pos_r(i), wtr )
         call spline_linear( pos_z(i), wtz )

         phase = cmplx( 1.0, 0.0 ) * this%q(pp)
         ! deposit m=0 mode
         do k = 0, 1
         do j = 0, 1
            q0(1,nn+j,mm+k) = q0(1,nn+j,mm+k) + wtr(j) * wtz(k) * real(phase)
         enddo
         enddo

         ! deposit m>0 mode
         do mode = 1, max_mode
            qr => q_re(mode)%get_f2()
            qi => q_im(mode)%get_f2()
            phase = phase * phase0(i)

            do k = 0, 1
            do j = 0, 1
               qr(1,nn+j,mm+k) = qr(1,nn+j,mm+k) + wtr(j) * wtz(k) * real(phase)
               qi(1,nn+j,mm+k) = qi(1,nn+j,mm+k) + wtr(j) * wtz(k) * aimag(phase)
            enddo
            enddo
         enddo
         pp = pp + 1
      enddo
      
   enddo

   if (noff1 == 0) then

      q0(1,0,:) = 0.0
      q0(1,1,:) = 8.0 * q0(1,1,:)
      do j = 2, n1p+1
         r = j + noff1 - 1
         q0(1,j,:) = q0(1,j,:) / r
      end do

      do mode = 1, max_mode
         qr => q_re(mode)%get_f2()
         qi => q_im(mode)%get_f2()
         qr(1,0,:) = 0.0
         qi(1,0,:) = 0.0
         qr(1,1,:) = 0.0
         qi(1,1,:) = 0.0
         do j = 2, n1p+1
            r = j + noff1 - 1
            qr(1,j,:) = qr(1,j,:) / r
            qi(1,j,:) = qi(1,j,:) / r
         enddo
      enddo

   else

      do j = 0, n1p+1
         r = j + noff1 - 1
         q0(1,j,:) = q0(1,j,:) / r
      end do

      do mode = 1, max_mode
         qr => q_re(mode)%get_f2()
         qi => q_im(mode)%get_f2()
         do j = 0, n1p+1
            r = j + noff1 - 1
            qr(1,j,:) = qr(1,j,:) / r
            qi(1,j,:) = qi(1,j,:) / r
         enddo
      enddo

   endif

   call stop_tprof( 'deposit 3D particles' )
   call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine qdeposit_part3d

subroutine push_boris_part3d( this, ef, bf )

   implicit none

   class(part3d), intent(inout) :: this
   class(field), intent(in) :: ef, bf

   ! local
   class(ufield), dimension(:), pointer :: ef_re, ef_im, bf_re, bf_im
   integer :: i, np, max_mode
   integer(kind=LG) :: ptrcur, pp
   real :: qtmh, ostq, u2, gam_qtmh
   real, dimension(p_cache_size) :: gam
   real, dimension(p_p_dim, p_cache_size) :: bp, ep, utmp, p_old
   character(len=32), save :: sname = "push_boris_part3d"

   call write_dbg(cls_name, sname, cls_level, 'starts')
   call start_tprof( 'push 3D particles' )

   qtmh = 0.5 * this%qbm * this%dt
   max_mode = ef%get_num_modes()

   ef_re => ef%get_rf_re()
   ef_im => ef%get_rf_im()
   bf_re => bf%get_rf_re()
   bf_im => bf%get_rf_im()

   do ptrcur = 1, this%npp, p_cache_size

      ! check if last copy of table and set np
      if( ptrcur + p_cache_size > this%npp ) then
        np = this%npp - ptrcur + 1
      else
        np = p_cache_size
      endif

      ! interpolate fields to particles
      call interp_emf( ef_re, ef_im, bf_re, bf_im, max_mode, this%x, &
         this%dr, this%dz, bp, ep, np, ptrcur )

      ! store old momenta for spin push
      if ( this%has_spin ) then
         pp = ptrcur
         do i = 1, np
            p_old(1,i) = this%p(1,pp)
            p_old(2,i) = this%p(2,pp)
            p_old(3,i) = this%p(3,pp)
            pp = pp + 1
         enddo
      endif

      do i = 1, np
         ep(1,i) = ep(1,i) * qtmh
         ep(2,i) = ep(2,i) * qtmh
         ep(3,i) = ep(3,i) * qtmh
      enddo

      pp = ptrcur
      do i = 1, np

         ! first half of electric field acc.
         utmp(1,i) = this%p(1,pp) + ep(1,i)
         utmp(2,i) = this%p(2,pp) + ep(2,i)
         utmp(3,i) = this%p(3,pp) + ep(3,i)

         ! time-centered momenta and gamma
         u2 = utmp(1,i)*utmp(1,i) + utmp(2,i)*utmp(2,i) + utmp(3,i)*utmp(3,i)
         gam(i) = sqrt( 1.0 + u2 )
         pp = pp + 1
      enddo

      do i = 1, np
         gam_qtmh = qtmh / gam(i)
         bp(1,i) = bp(1,i) * gam_qtmh
         bp(2,i) = bp(2,i) * gam_qtmh
         bp(3,i) = bp(3,i) * gam_qtmh
      enddo

      ! push spin
      if ( this%has_spin ) then
         call this%push_spin( ep, bp, p_old, gam, ptrcur, np )
      endif

      pp = ptrcur
      do i = 1, np
         this%p(1,pp) = utmp(1,i) + utmp(2,i) * bp(3,i)
         this%p(2,pp) = utmp(2,i) + utmp(3,i) * bp(1,i)
         this%p(3,pp) = utmp(3,i) + utmp(1,i) * bp(2,i)

         this%p(1,pp) = this%p(1,i) - utmp(3,i) * bp(2,i)
         this%p(2,pp) = this%p(2,i) - utmp(1,i) * bp(3,i)
         this%p(3,pp) = this%p(3,i) - utmp(2,i) * bp(1,i)

         pp = pp + 1
      enddo

      do i = 1, np
         ostq = 2.0 / ( 1.0 + bp(1,i)**2 + bp(2,i)**2 + bp(3,i)**2 )
         bp(1,i) = bp(1,i) * ostq
         bp(2,i) = bp(2,i) * ostq
         bp(3,i) = bp(3,i) * ostq
      enddo

      pp = ptrcur
      do i = 1, np
         utmp(1,i) = utmp(1,i) + this%p(2,pp) * bp(3,i) - this%p(3,pp) * bp(2,i)
         utmp(2,i) = utmp(2,i) + this%p(3,pp) * bp(1,i) - this%p(1,pp) * bp(3,i)
         utmp(3,i) = utmp(3,i) + this%p(1,pp) * bp(2,i) - this%p(2,pp) * bp(1,i)
         pp = pp + 1
      enddo

      ! second half of electric field acc.
      pp = ptrcur
      do i = 1, np
         this%p(1,pp) = utmp(1,i) + ep(1,i)
         this%p(2,pp) = utmp(2,i) + ep(2,i)
         this%p(3,pp) = utmp(3,i) + ep(3,i)
         pp = pp + 1
      enddo

      ! advance the particle positions
      ! note the x(3) is xi instead of z
      pp = ptrcur
      do i = 1, np
         gam_qtmh = this%dt / &
            sqrt( 1.0 + this%p(1,pp)**2 + this%p(2,pp)**2 + this%p(3,pp)**2 )
         this%x(1,pp) = this%x(1,pp) + this%p(1,pp) * gam_qtmh
         this%x(2,pp) = this%x(2,pp) + this%p(2,pp) * gam_qtmh
         this%x(3,pp) = this%x(3,pp) - this%p(3,pp) * gam_qtmh + this%dt
         pp = pp + 1
      enddo

   enddo

   call stop_tprof( 'push 3D particles' )
   call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine push_boris_part3d

subroutine push_reduced_part3d( this, ef, bf )

   implicit none

   class(part3d), intent(inout) :: this
   class(field), intent(in) :: ef, bf
   ! local
   class(ufield), dimension(:), pointer :: ef_re, ef_im, bf_re, bf_im
   
   integer :: i, np, max_mode
   real :: qtmh, dt_gam
   real, dimension(p_p_dim, p_cache_size) :: bp, ep, p_old
   real, dimension(p_cache_size) :: gam
   integer(kind=LG) :: ptrcur, pp
   character(len=32), save :: sname = "push_reduced_part3d"

   call write_dbg(cls_name, sname, cls_level, 'starts')
   call start_tprof( 'push 3D particles' )

   qtmh = this%qbm * this%dt * 0.5
   max_mode = ef%get_num_modes()

   ef_re => ef%get_rf_re()
   ef_im => ef%get_rf_im()
   bf_re => bf%get_rf_re()
   bf_im => bf%get_rf_im()

   do ptrcur = 1, this%npp, p_cache_size

      ! check if last copy of table and set np
      if( ptrcur + p_cache_size > this%npp ) then
        np = this%npp - ptrcur + 1
      else
        np = p_cache_size
      endif

      ! interpolate fields to particles
      call interp_emf( ef_re, ef_im, bf_re, bf_im, max_mode, this%x, this%dr, &
         this%dz, bp, ep, np, ptrcur )

      ! store old momenta for spin push
      if ( this%has_spin ) then
         pp = ptrcur
         do i = 1, np
            p_old(1,i) = this%p(1,pp)
            p_old(2,i) = this%p(2,pp)
            p_old(3,i) = this%p(3,pp)
            pp = pp + 1
         enddo
      endif

      do i = 1, np
         ep(1,i) = ep(1,i) * qtmh
         ep(2,i) = ep(2,i) * qtmh
         ep(3,i) = ep(3,i) * qtmh
         bp(1,i) = bp(1,i) * qtmh
         bp(2,i) = bp(2,i) * qtmh
         bp(3,i) = bp(3,i) * qtmh
         ! calculate transverse force
         ep(1,i) = ep(1,i) - bp(2,i)
         ep(2,i) = ep(2,i) + bp(1,i)
      enddo

      pp = ptrcur
      do i = 1, np
         ! half advance momenta
         this%p(1,pp) = this%p(1,pp) + ep(1,i)
         this%p(2,pp) = this%p(2,pp) + ep(2,i)
         this%p(3,pp) = this%p(3,pp) + ep(3,i)
         gam(i) = sqrt( 1.0 + this%p(1,pp)**2 + this%p(2,pp)**2 + this%p(3,pp)**2 )
         pp = pp + 1
      enddo

      pp = ptrcur
      do i = 1, np
         ! half advance momenta
         this%p(1,pp) = this%p(1,pp) + ep(1,i)
         this%p(2,pp) = this%p(2,pp) + ep(2,i)
         this%p(3,pp) = this%p(3,pp) + ep(3,i)
         pp = pp + 1
      enddo

      ! push spin
      if ( this%has_spin ) then
         call this%push_spin( ep, bp, p_old, gam, ptrcur, np )
      endif

      ! advance the particle positions
      ! note the x(3) is xi instead of z
      pp = ptrcur
      do i = 1, np
         dt_gam = this%dt / &
            sqrt( 1.0 + this%p(1,pp)**2 + this%p(2,pp)**2 + this%p(3,pp)**2 )
         this%x(1,pp) = this%x(1,pp) + this%p(1,pp) * dt_gam
         this%x(2,pp) = this%x(2,pp) + this%p(2,pp) * dt_gam
         this%x(3,pp) = this%x(3,pp) - this%p(3,pp) * dt_gam + this%dt
         pp = pp + 1
      enddo

   enddo

   call stop_tprof( 'push 3D particles' )
   call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine push_reduced_part3d

subroutine push_spin_part3d( this, ep, bp, p_old, gam, ptrcur, np )

   implicit none

   class(part3d), intent(inout) :: this
   real, intent(in), dimension(:,:) :: ep, bp, p_old
   real, intent(in), dimension(:) :: gam
   integer, intent(in) :: np
   integer(kind=LG), intent(in) :: ptrcur

   ! local data
   integer :: i, pp
   real, dimension(p_p_dim, np) :: omega, vtemp, stemp
   real :: a, coef, vdotb
   character(len=32), save :: sname = "push_spin_part3d"

   call write_dbg(cls_name, sname, cls_level, 'starts')

   a = this%anom_mag_moment

   pp = ptrcur
   do i = 1, np

      ! calculate the time-centered velocity
      vtemp(1,i) = 0.5 * ( p_old(1,i) + this%p(1,pp) ) / gam(i)
      vtemp(2,i) = 0.5 * ( p_old(2,i) + this%p(2,pp) ) / gam(i)
      vtemp(3,i) = 0.5 * ( p_old(3,i) + this%p(3,pp) ) / gam(i)

      ! Now calculate the precession frequency Omega
      ! calculate contribution from B
      coef = a + 1.0 / gam(i)
      omega(1,i) = coef * bp(1,i) * gam(i)
      omega(2,i) = coef * bp(2,i) * gam(i)
      omega(3,i) = coef * bp(3,i) * gam(i)

      ! calculate (v cross E) contribution
      coef = -1.0 * ( a + 1.0 / ( 1.0 + gam(i) ) )
      omega(1,i) = omega(1,i) + coef * ( vtemp(2,i) * ep(3,i) - vtemp(3,i) * ep(2,i) )
      omega(2,i) = omega(2,i) + coef * ( vtemp(3,i) * ep(1,i) - vtemp(1,i) * ep(3,i) )
      omega(3,i) = omega(3,i) + coef * ( vtemp(1,i) * ep(2,i) - vtemp(2,i) * ep(1,i) )

      ! calculate (v dot B) contribution
      ! note that omega should be multiplied by 2 because ep and bp are already
      ! divided by 2 in dudt. But rotation vector is obtained via d = omega * dt/2,
      ! therefore the factor 2 is canceled.
      vdotb = vtemp(1,i) * bp(1,i) + vtemp(2,i) * bp(2,i) + vtemp(3,i) * bp(3,i)
      coef = -1.0 * ( a * gam(i)**2 / ( 1.0 + gam(i) ) * vdotb )
      omega(1,i) = omega(1,i) + coef * vtemp(1,i)
      omega(2,i) = omega(2,i) + coef * vtemp(2,i)
      omega(3,i) = omega(3,i) + coef * vtemp(3,i)

      ! calculate s_prime
      stemp(1,i) = this%s(1,pp) + ( this%s(2,pp) * omega(3,i) - this%s(3,pp) * omega(2,i) )
      stemp(2,i) = this%s(2,pp) + ( this%s(3,pp) * omega(1,i) - this%s(1,pp) * omega(3,i) )
      stemp(3,i) = this%s(3,pp) + ( this%s(1,pp) * omega(2,i) - this%s(2,pp) * omega(1,i) )

      coef = 2.0 / ( 1.0 + omega(1,i)**2 + omega(2,i)**2 + omega(3,i)**2 )
      this%s(1,pp) = this%s(1,pp) + coef * ( stemp(2,i) * omega(3,i) - stemp(3,i) * omega(2,i) )
      this%s(2,pp) = this%s(2,pp) + coef * ( stemp(3,i) * omega(1,i) - stemp(1,i) * omega(3,i) )
      this%s(3,pp) = this%s(3,pp) + coef * ( stemp(1,i) * omega(2,i) - stemp(2,i) * omega(1,i) )

      pp = pp + 1
   enddo

   call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine push_spin_part3d

subroutine update_bound_part3d( this )

   implicit none

   class(part3d), intent(inout) :: this
   ! local
   integer(kind=LG) :: i = 1
   real :: pos_r, pos_z
   character(len=32), save :: sname = "update_bound_part3d"

   call write_dbg(cls_name, sname, cls_level, 'starts')
   call start_tprof( 'push 3D particles' )

   if ( this%npp == 0 ) return

   do while ( i < this%npp )

      pos_r = sqrt( this%x(1,i)**2 + this%x(2,i)**2 )
      pos_z = this%x(3,i)

      ! check if particle goes out of r-edge
      if ( pos_r > this%edge(1) ) then
         this%x(:,i) = this%x(:, this%npp)
         this%p(:,i) = this%p(:, this%npp)
         this%q(i)   = this%q(this%npp)
         if ( this%has_spin ) then
            this%s(:,i) = this%s(:, this%npp)
         endif
         this%npp = this%npp - 1
         cycle
      endif

      ! check if particle goes out of z-edge
      if ( pos_z > this%edge(2) ) then
         this%x(:,i) = this%x(:, this%npp)
         this%p(:,i) = this%p(:, this%npp)
         this%q(i)   = this%q(this%npp)
         if ( this%has_spin ) then
            this%s(:,i) = this%s(:, this%npp)
         endif
         this%npp = this%npp - 1
         cycle
      endif

      i = i + 1
   enddo

   ! deal with the last particle
   pos_r = sqrt( this%x(1,this%npp)**2 + this%x(2,this%npp)**2 )
   pos_z = this%x(3,this%npp)

   if ( pos_r > this%edge(1) .or. pos_z > this%edge(2) ) then
      this%npp = this%npp - 1
   endif

   call stop_tprof( 'push 3D particles' )
   call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine update_bound_part3d

subroutine interp_emf( ef_re, ef_im, bf_re, bf_im, num_modes, x, dr, dz, bp, ep, np, ptrcur )

   implicit none

   class(ufield), dimension(:), pointer, intent(in) :: ef_re, ef_im, bf_re, bf_im
   integer, intent(in) :: num_modes, np
   real, intent(in) :: dr, dz
   real, dimension(:,:), intent(in) :: x
   real, dimension(:,:), intent(inout) :: bp, ep
   integer(kind=LG), intent(in) :: ptrcur

   real, dimension(:,:,:), pointer :: e0, b0, er, ei, br, bi
   integer :: noff1, noff2, n1p, n2p, n1, n2
   integer :: i, j, k, nn, mm, mode
   integer(kind=LG) :: pp
   real :: pos_r, pos_z, idr, idz, wt, ph_r, ph_i, cc, ss
   real, dimension(0:1) :: wtr, wtz
   complex(kind=DB) :: phase0, phase

   idr = 1.0 / dr
   idz = 1.0 / dz

   noff1 = ef_re(0)%get_noff(1)
   n1p   = ef_re(0)%get_ndp(1)
   n1    = ef_re(0)%get_nd(1)
   noff2 = ef_re(0)%get_noff(2)
   n2p   = ef_re(0)%get_ndp(2)
   n2    = ef_re(0)%get_nd(2)

   e0 => ef_re(0)%get_f2()
   b0 => bf_re(0)%get_f2()

   ep = 0.0
   bp = 0.0

   pp = ptrcur
   do i = 1, np
      pos_r = sqrt( x(1,pp)**2 + x(2,pp)**2 ) * idr
      pos_z = x(3,pp) * idz
      ! cosine and sine
      cc = x(1,pp) / pos_r * idr
      ss = x(2,pp) / pos_r * idr
      phase0 = cmplx( cc, ss )

      nn = int( pos_r )
      mm = int( pos_z )

      ! in-cell position
      pos_r = pos_r - real(nn)
      pos_z = pos_z - real(mm)

      ! cell index
      nn = nn - noff1 + 1
      mm = mm - noff2 + 1

      call spline_linear( pos_r, wtr )
      call spline_linear( pos_z, wtz )

      phase = cmplx( 1.0, 0.0 )

      ! interpolate m=0 mode
      do k = 0, 1
      do j = 0, 1
         wt = wtr(j) * wtz(k)
         ep(1,i) = ep(1,i) + e0(1,nn+j,mm+k) * wt
         ep(2,i) = ep(2,i) + e0(2,nn+j,mm+k) * wt
         ep(3,i) = ep(3,i) + e0(3,nn+j,mm+k) * wt
         bp(1,i) = bp(1,i) + b0(1,nn+j,mm+k) * wt
         bp(2,i) = bp(2,i) + b0(2,nn+j,mm+k) * wt
         bp(3,i) = bp(3,i) + b0(3,nn+j,mm+k) * wt
      enddo
      enddo

      ! interpolate m>0 modes
      do mode = 1, num_modes
         phase = phase * phase0
         ph_r = 2.0 * real(phase)
         ph_i = 2.0 * aimag(phase)

         er => ef_re(mode)%get_f2()
         ei => ef_im(mode)%get_f2()
         br => bf_re(mode)%get_f2()
         bi => bf_im(mode)%get_f2()

         do k = 0, 1
         do j = 0, 1
            wt = wtr(j) * wtz(k)
            ep(1,i) = ep(1,i) + ( er(1,nn+j,mm+k)*ph_r - ei(1,nn+j,mm+k)*ph_i ) * wt
            ep(2,i) = ep(2,i) + ( er(2,nn+j,mm+k)*ph_r - ei(2,nn+j,mm+k)*ph_i ) * wt
            ep(3,i) = ep(3,i) + ( er(3,nn+j,mm+k)*ph_r - ei(3,nn+j,mm+k)*ph_i ) * wt
            bp(1,i) = bp(1,i) + ( br(1,nn+j,mm+k)*ph_r - bi(1,nn+j,mm+k)*ph_i ) * wt
            bp(2,i) = bp(2,i) + ( br(2,nn+j,mm+k)*ph_r - bi(2,nn+j,mm+k)*ph_i ) * wt
            bp(3,i) = bp(3,i) + ( br(3,nn+j,mm+k)*ph_r - bi(3,nn+j,mm+k)*ph_i ) * wt
         enddo
         enddo
      enddo

      ! transform from cylindrical geometry to Cartesian geometry
      ! ph_r, ph_i are temporary variables here
      ph_r = ep(1,i) * cc - ep(2,i) * ss
      ph_i = ep(1,i) * ss + ep(2,i) * cc
      ep(1,i) = ph_r
      ep(2,i) = ph_i

      ph_r = bp(1,i) * cc - bp(2,i) * ss
      ph_i = bp(1,i) * ss + bp(2,i) * cc
      bp(1,i) = ph_r
      bp(2,i) = ph_i

      pp = pp + 1
   enddo

end subroutine interp_emf
!
subroutine pmove(this,fd,rtag,stag,sid)

   implicit none

   class(part3d), intent(inout) :: this
   class(field), intent(in) :: fd
   integer, intent(in) :: rtag, stag
   integer, intent(inout) :: sid
! local data
   character(len=18), save :: sname = 'pmove:'
   class(ufield), pointer :: ud
   integer, dimension(9) :: info


   call write_dbg(cls_name, sname, cls_level, 'starts')

   ud => fd%get_rf_re(0)

   if ( this%has_spin ) then
      call pmove_part3d(this%x, this%p, this%q, this%pp,ud,this%npp,this%dr,this%dz,sbufr,sbufl,&
      &rbufr,rbufl,ihole,this%pbuff,this%part_dim,this%npmax,this%nbmax,rtag,stag,sid,info, this%s)
   else
      call pmove_part3d(this%x, this%p, this%q, this%pp,ud,this%npp,this%dr,this%dz,sbufr,sbufl,&
      &rbufr,rbufl,ihole,this%pbuff,this%part_dim,this%npmax,this%nbmax,rtag,stag,sid,info)
   endif


   if (info(1) /= 0) then
      call write_err(cls_name//sname//'pmove_part3d error')
   endif

   if (this%pp%getstageid() == this%pp%getnstage() - 1) sid = MPI_REQUEST_NULL

   call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine pmove
!
function getnpp(this)

   implicit none

   class(part3d), intent(in) :: this
   integer(kind=LG) :: getnpp

   getnpp = this%npp

end function getnpp
!
subroutine writehdf5_part3d(this,file,dspl,rtag,stag,id)

   implicit none

   class(part3d), intent(inout) :: this
   class(hdf5file), intent(in) :: file
   integer, intent(in) :: dspl, rtag, stag
   integer, intent(inout) :: id
! local data
   character(len=18), save :: sname = 'writehdf5_part3d'
   integer :: ierr = 0

   call write_dbg(cls_name, sname, cls_level, 'starts')
   if ( this%has_spin ) then
      call pwpart_pipe(this%pp,file,this%x, this%p, this%q,this%npp,dspl,&
      &this%z0,rtag,stag,id,ierr, this%s)
   else
      call pwpart_pipe(this%pp,file,this%x, this%p, this%q,this%npp,dspl,&
      &this%z0,rtag,stag,id,ierr)
   endif
   call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine writehdf5_part3d
!
subroutine writerst_part3d(this,file)

   implicit none

   class(part3d), intent(inout) :: this
   class(hdf5file), intent(in) :: file
! local data
   character(len=18), save :: sname = 'writerst_part3d'
   integer :: ierr = 0

   call write_dbg(cls_name, sname, cls_level, 'starts')
   if ( this%has_spin ) then
      call wpart(this%pp,file,this%x, this%p, this%q,this%npp,1,ierr, this%s)
   else
      call wpart(this%pp,file,this%x, this%p, this%q,this%npp,1,ierr)
   endif
   call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine writerst_part3d
!
subroutine readrst_part3d(this,file)

   implicit none

   class(part3d), intent(inout) :: this
   class(hdf5file), intent(in) :: file
! local data
   character(len=18), save :: sname = 'readrst_part3d'
   integer :: ierr = 0

   call write_dbg(cls_name, sname, cls_level, 'starts')
   if ( this%has_spin ) then
      call rpart(this%pp,file,this%x, this%p, this%q,this%npp,ierr, this%s)
   else
      call rpart(this%pp,file,this%x, this%p, this%q,this%npp,ierr)
   endif
   call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine readrst_part3d
!
end module part3d_class
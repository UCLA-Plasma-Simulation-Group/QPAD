module part3d_tmplt_class

use param
use sysutil_module
use field_class
use ufield_class
use part3d_class
use interpolation

implicit none

private

public :: part3d_tmplt

type, extends( part3d ) :: part3d_tmplt

   ! private

   contains

   ! new procedure
   procedure :: push_tmplt => push_tmplt_part3d

end type

save

character(len=32), parameter :: cls_name = "part3d_tmplt"
integer, parameter :: cls_level = 2

contains

! this pusher is exactly the Boris pusher
subroutine push_tmplt_part3d( this, ef, bf )

   implicit none

   class(part3d_tmplt), intent(inout) :: this
   class(field), intent(in) :: ef, bf

   ! local
   type(ufield), dimension(:), pointer :: ef_re, ef_im, bf_re, bf_im
   integer :: i, np, max_mode
   integer(kind=LG) :: ptrcur, pp
   real :: qtmh, ostq, u2, gam_qtmh
   real, dimension(p_cache_size) :: gam
   real, dimension(p_p_dim, p_cache_size) :: bp, ep, utmp, p_old
   character(len=32), save :: sname = "push_tmplt_part3d"

   call write_dbg(cls_name, sname, cls_level, 'starts')
   call start_tprof( 'push 3D particles' )

   qtmh = 0.5 * this%qbm * this%dt
   max_mode = ef%get_max_mode()

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
      call interp_emf_part3d( ef_re, ef_im, bf_re, bf_im, max_mode, this%x, &
         this%dr, this%dz, bp, ep, np, ptrcur )

      ! store old momenta for spin push
      if ( this%has_spin ) then
         pp = ptrcur
         do i = 1, np
            p_old(:,i) = this%p(:,pp)
            pp = pp + 1
         enddo
      endif

      pp = ptrcur
      do i = 1, np
         ep(:,i) = ep(:,i) * qtmh

         ! first half of electric field acc.
         utmp(:,i) = this%p(:,pp) + ep(:,i)

         ! time-centered momenta and gamma
         u2 = utmp(1,i)*utmp(1,i) + utmp(2,i)*utmp(2,i) + utmp(3,i)*utmp(3,i)
         gam(i) = sqrt( 1.0 + u2 )
         pp = pp + 1
      enddo

      do i = 1, np
         gam_qtmh = qtmh / gam(i)
         bp(:,i) = bp(:,i) * gam_qtmh
      enddo

      ! push spin
      if ( this%has_spin ) then
         call this%push_spin( ep, bp, p_old, gam, ptrcur, np )
      endif

      pp = ptrcur
      do i = 1, np
         this%p(1,pp) = utmp(1,i) + utmp(2,i) * bp(3,i) - utmp(3,i) * bp(2,i)
         this%p(2,pp) = utmp(2,i) + utmp(3,i) * bp(1,i) - utmp(1,i) * bp(3,i)
         this%p(3,pp) = utmp(3,i) + utmp(1,i) * bp(2,i) - utmp(2,i) * bp(1,i)
         pp = pp + 1
      enddo

      do i = 1, np
         ostq = 2.0 / ( 1.0 + bp(1,i)**2 + bp(2,i)**2 + bp(3,i)**2 )
         bp(:,i) = bp(:,i) * ostq
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
         this%p(:,pp) = utmp(:,i) + ep(:,i)
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

end subroutine push_tmplt_part3d

subroutine interp_emf_part3d( ef_re, ef_im, bf_re, bf_im, max_mode, x, dr, dz, bp, ep, np, ptrcur )

   implicit none

   type(ufield), dimension(:), pointer, intent(in) :: ef_re, ef_im, bf_re, bf_im
   integer, intent(in) :: max_mode, np
   real, intent(in) :: dr, dz
   real, dimension(:,:), intent(in) :: x
   real, dimension(:,:), intent(inout) :: bp, ep
   integer(kind=LG), intent(in) :: ptrcur

   real, dimension(:,:,:), pointer :: e0, b0, er, ei, br, bi
   integer :: noff1, noff2, i, j, k, nn, mm, mode
   integer(kind=LG) :: pp
   real :: pos_r, pos_z, idr, idz, wt, ph_r, ph_i, cc, ss
   real, dimension(0:1) :: wtr, wtz
   complex(kind=DB) :: phase0, phase

   idr = 1.0 / dr
   idz = 1.0 / dz

   noff1 = ef_re(0)%get_noff(1)
   noff2 = ef_re(0)%get_noff(2)

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
         ep(:,i) = ep(:,i) + e0(:,nn+j,mm+k) * wt
         bp(:,i) = bp(:,i) + b0(:,nn+j,mm+k) * wt
      enddo
      enddo

      ! interpolate m>0 modes
      do mode = 1, max_mode
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
            ep(:,i) = ep(:,i) + ( er(:,nn+j,mm+k) * ph_r - ei(:,nn+j,mm+k) * ph_i ) * wt
            bp(:,i) = bp(:,i) + ( br(:,nn+j,mm+k) * ph_r - bi(:,nn+j,mm+k) * ph_i ) * wt
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

end subroutine interp_emf_part3d

end module part3d_tmplt_class

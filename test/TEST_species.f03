program test_species

use parallel_pipe_class
use grid_class
use field_e_class
use field_class
use field_src_class
use fdist2d_class
use system
use param
use species2d_class
use hdf5io_class
use input_class
use mpi
use ufield_class
use debug_tool

implicit none

type(species2d) :: spe
type(parallel_pipe), pointer :: pp => null()
type(grid), pointer :: gp => null()
type(hdf5file) :: file
type(input_json), pointer :: input => null()
integer :: nr, nz, nrp, noff, xdim, npf
integer :: num_modes, part_shape, i, id
character(len=:), allocatable :: shape
character(len=2) :: s1
real :: dr, dxi, rmin, rmax, zmin, zmax
type(field_rho), pointer :: rho
type fdist2d_wrap
   class(fdist2d), allocatable :: p
end type fdist2d_wrap
type(fdist2d_wrap) :: pf2d
type(hdf5file) :: file_spe
type(hdf5file), dimension(:), allocatable :: file_q

! type( field_e ) :: e
! type( field_jay ) :: jay

! real :: dr, dxi, r

! type( ufield ), dimension(:), pointer :: ujay_re => null(), ujay_im => null()
! type( ufield ), dimension(:), pointer :: ue_re => null(), ue_im => null()
! real, dimension(:,:), pointer :: p
! integer :: ierr, i, mode
! character(len=32) :: filename

allocate(pp,gp,input,rho)
call input%new()
pp => input%pp

call write_dbg( 'main', 'test_species', 0, 'starts' )

call input%get('simulation.max_mode',num_modes)
allocate(file_q(2*num_modes+1))
call input%get('simulation.interpolation',shape)
select case (trim(shape))
case('linear')
   part_shape = p_ps_linear
end select

call input%get('simulation.grid(1)',nr)
call input%get('simulation.grid(2)',nz)
call gp%new( pp, nr, nz )

call input%get('simulation.box.r(1)',rmin)
call input%get('simulation.box.r(2)',rmax)
dr = (rmax-rmin)/nr
call input%get('simulation.box.z(1)',zmin)
call input%get('simulation.box.z(2)',zmax)
dxi = (zmax-zmin)/nz

nrp = gp%get_ndp(1)
noff = gp%get_noff(1)

xdim = 8

call input%get('species(1).profile',npf)
select case (npf)
case (0)
   allocate(fdist2d_000::pf2d%p)
   call pf2d%p%new(input,1)
end select

call rho%new(pp, gp, dr, dxi, num_modes, part_shape)
call spe%new(pp,gp,pf2d%p,part_shape,dr,dxi,num_modes,-1.0,dxi,xdim,0.0)

call file_spe%new(&
&timeunit = '1 / \omega_p',&
&dt = dxi,&
&ty = 'particles',&
&filename = './spe_',&
&dataname = 'raw',&
&units = '',&
&label = 'Plasma Raw',&
&n = 0,&
&t = 0.0)
call spe%wr(file_spe)
do i = 2, nz
   call spe%cbq(i)
end do

call file_q(1)%new(&
&axismin = (/rmin,zmin,0.0/),&
&axismax = (/rmax,zmax,1.0/),&
&axisname  = (/'r  ','\xi','z  '/),&
&axislabel = (/'r  ','\xi','z  '/),&
&timeunit = '1 / \omega_p',&
&dt = dxi,&
&axisunits = (/'c / \omega_p','c / \omega_p','c / \omega_p'/),&
&rank = 2,&
&filename = './',&
&dataname = 'q0',&
&units = 'n_0',&
&label = 'Charge Density',&
&n = 0,&
&t = 0.0)


do i = 1, num_modes
   write (s1, '(I1.1)') i
   call file_q(2*i)%new(&
   &axismin = (/rmin,zmin,0.0/),&
   &axismax = (/rmax,zmax,1.0/),&
   &axisname  = (/'r  ','\xi','z  '/),&
   &axislabel = (/'r  ','\xi','z  '/),&
   &timeunit = '1 / \omega_p',&
   &dt = dxi,&
   &axisunits = (/'c / \omega_p','c / \omega_p','c / \omega_p'/),&
   &rank = 2,&
   &filename = './',&
   &dataname = 'qr'//trim(s1),&
   &units = 'n_0',&
   &label = 'Charge Density',&
   &n = 0,&
   &t = 0.0)
   call file_q(2*i+1)%new(&
   &axismin = (/rmin,zmin,0.0/),&
   &axismax = (/rmax,zmax,1.0/),&
   &axisname  = (/'r  ','\xi','z  '/),&
   &axislabel = (/'r  ','\xi','z  '/),&
   &timeunit = '1 / \omega_p',&
   &dt = dxi,&
   &axisunits = (/'c / \omega_p','c / \omega_p','c / \omega_p'/),&
   &rank = 2,&
   &filename = './',&
   &dataname = 'qi'//trim(s1),&
   &units = 'n_0',&
   &label = 'Charge Density',&
   &n = 0,&
   &t = 0.0)
end do
call spe%wrq(file_q,5,5,id)


! call jay%new( pp, gp, dr, dxi, num_modes, part_shape )
! call e%new( pp, gp, dr, dxi, num_modes, part_shape, entity=p_entity_plasma )


! ujay_re => jay%get_rf_re()
! ujay_im => jay%get_rf_im()

! ! set the charge
! do mode = 0, num_modes
  
!   do i = 1, nrp
!     r = (real(i+noff)-0.5)*dr
!     ujay_re(mode)%f1(1,i) = exp( -((r-0.5)/0.05)**2 )
!     ujay_re(mode)%f1(2,i) = exp( -((r-0.5)/0.05)**2 )
!     ujay_re(mode)%f1(3,i) = exp( -((r-0.5)/0.05)**2 )
!   enddo

!   if ( mode == 0 ) cycle

!   do i = 1, nrp
!     r = (real(i+noff)-0.5)*dr
!     ujay_im(mode)%f1(1,i) = exp( -((r-0.5)/0.05)**2 )
!     ujay_im(mode)%f1(2,i) = exp( -((r-0.5)/0.05)**2 )
!     ujay_im(mode)%f1(3,i) = exp( -((r-0.5)/0.05)**2 )
!   enddo

! enddo

! call e%solve(jay)

! ue_re => e%get_rf_re()
! ue_im => e%get_rf_im()

! p => ue_re(0)%get_f1()
! write( filename, '(A,I0.3,A)' ) 'ez-re-0-', pp%getlidproc(), '.txt'
! call write_data( p, trim(filename), 3 )
! p => ue_re(1)%get_f1()
! write( filename, '(A,I0.3,A)' ) 'ez-re-1-', pp%getlidproc(), '.txt'
! call write_data( p, trim(filename), 3 )
! p => ue_re(2)%get_f1()
! write( filename, '(A,I0.3,A)' ) 'ez-re-2-', pp%getlidproc(), '.txt'
! call write_data( p, trim(filename), 3 )
! p => ue_im(1)%get_f1()
! write( filename, '(A,I0.3,A)' ) 'ez-im-1-', pp%getlidproc(), '.txt'
! call write_data( p, trim(filename), 3 )
! p => ue_im(2)%get_f1()
! write( filename, '(A,I0.3,A)' ) 'ez-im-2-', pp%getlidproc(), '.txt'
! call write_data( p, trim(filename), 3 )

! call jay%del()
! call e%del()

call write_dbg( 'main', 'test_species', 0, 'ends' )

call end_errors()

end program test_species
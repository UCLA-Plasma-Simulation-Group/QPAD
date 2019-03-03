program test_species

use parallel_pipe_class
use grid_class
use field_e_class
use field_src_class
use field_b_class
use field_e_class
use field_psi_class
use fdist2d_class
use fdist3d_class
use system
use param
use species2d_class
use beam3d_class
use hdf5io_class
use input_class
use mpi
use debug_tool

implicit none

type(species2d) :: spe
type(beam3d) :: beam
type(parallel_pipe), pointer :: pp => null()
type(grid), pointer :: gp => null()
type(hdf5file) :: file
type(input_json), pointer :: input => null()
integer :: nr, nz, nrp, noff, xdim, npf, ierr
integer :: num_modes, part_shape, i, id
character(len=:), allocatable :: shape
character(len=2) :: s1
real :: dr, dxi, rmin, rmax, zmin, zmax, dt
type(field_rho), pointer :: rho, qb
type fdist2d_wrap
   class(fdist2d), allocatable :: p
end type fdist2d_wrap
type(fdist2d_wrap) :: pf2d
type fdist3d_wrap
   class(fdist3d), allocatable :: p
end type fdist3d_wrap
type(fdist3d_wrap) :: pf3d
type(hdf5file) :: file_spe
type(hdf5file), dimension(:), allocatable :: file_q
type(hdf5file) :: file_beam
type(hdf5file), dimension(:), allocatable :: file_qb

type( field_b ) :: b
type( field_e ) :: e
type( field_psi ) :: psi

! Initialization
allocate(pp,gp,input,rho,qb)
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

! call b%new( pp, gp, dr, dxi, num_modes, part_shape, entity=p_entity_plasma )
! call e%new( pp, gp, dr, dxi, num_modes, part_shape, entity=p_entity_plasma )
call b%new( pp, gp, dr, dxi, num_modes, part_shape, entity=p_entity_beam )
call e%new( pp, gp, dr, dxi, num_modes, part_shape, entity=p_entity_beam )
call psi%new( pp, gp, dr, dxi, num_modes, part_shape )

! do i = 1, nz
!    rho = 0.0
!    e = 0.0
!    psi = 0.0
!    b = 0.0
!    call spe%qdp(rho)
!    call psi%solve(rho)
!    call e%solve(b,psi)
!    call rho%copy_slice(i, p_copy_1to2)
!    call e%copy_slice(i, p_copy_1to2)
!    call b%copy_slice(i, p_copy_1to2)
!    call psi%copy_slice(i, p_copy_1to2)
! end do

call input%get('beam(1).profile',npf)
select case (npf)
case (0)
   allocate(fdist3d_000::pf3d%p)
   call pf3d%p%new(input,1)
end select

call input%get('simulation.dt',dt)
call qb%new(pp, gp, dr, dxi, num_modes, part_shape)
call qb%as(0.0)
call beam%new(pp,qb,gp,part_shape,pf3d%p,-1.0,dt,7)
call beam%qdp(qb)

do i = 1, nz
   write(2,*) "Step #", i
   call file_spe%new(&
   &timeunit = '1 / \omega_p',&
   &dt = dxi,&
   &ty = 'particles',&
   &filename = './spe_',&
   &dataname = 'raw',&
   &units = '',&
   &label = 'Plasma Raw',&
   &n = i,&
   &t = real(i))
   call spe%wr(file_spe)
   call qb%copy_slice(i,p_copy_2to1)
   rho = 0.0
   ! psi = 0.0
   call b%solve(qb)
   call e%solve(b)
   call spe%qdp(rho)
   call spe%push(e,b)
   ! call psi%solve(rho)
   ! call e%solve(b,psi)
   call rho%copy_slice(i, p_copy_1to2)
   call e%copy_slice(i, p_copy_1to2)
   call b%copy_slice(i, p_copy_1to2)
   ! call psi%copy_slice(i, p_copy_1to2)
end do

call file_beam%new(&
&timeunit = '1 / \omega_p',&
&dt = dt,&
&ty = 'particles',&
&filename = './beam_',&
&dataname = 'raw',&
&units = '',&
&label = 'Beam Raw',&
&n = 0,&
&t = 0.0)
call beam%wr(file_beam,1,6,6,id)

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
&dataname = 'q',&
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
! call spe%wrq(file_q,5,5,id)
call rho%write_hdf5(file_q,1,5,5,id)

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
&dataname = 'qb0',&
&units = 'n_0',&
&label = 'Beam Charge Density',&
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
   &dataname = 'qbr'//trim(s1),&
   &units = 'n_0',&
   &label = 'Beam Density',&
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
   &dataname = 'qbi'//trim(s1),&
   &units = 'n_0',&
   &label = 'Beam Density',&
   &n = 0,&
   &t = 0.0)
end do
call beam%wrq(file_q,8,8,id)


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
&dataname = 'psi',&
&units = 'mc^2/e',&
&label = '\Psi (0 Mode)',&
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
   &dataname = 'psir'//trim(s1),&
   &units = 'mc^2/e',&
   &label = 'Real Part of \Psi ('//trim(s1)//' Mode',&
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
   &dataname = 'psii'//trim(s1),&
   &units = 'mc^2/e',&
   &label = 'Imaginary Part of \Psi ('//trim(s1)//' Mode',&
   &n = 0,&
   &t = 0.0)
end do

call psi%write_hdf5(file_q,1,8,8,id)


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
&dataname = 'ex',&
&units = 'mc\omega_p/e',&
&label = 'E_x (0 Mode)',&
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
   &dataname = 'exr'//trim(s1),&
   &units = 'mc\omega_p/e',&
   &label = 'Real Part of E_x ('//trim(s1)//' Mode',&
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
   &dataname = 'exi'//trim(s1),&
   &units = 'mc\omega_p/e',&
   &label = 'Imaginary Part of E_x ('//trim(s1)//' Mode',&
   &n = 0,&
   &t = 0.0)
end do

call e%write_hdf5(file_q,1,8,8,id)

call file_q(1)%new(&
&dataname = 'ey',&
&label = 'E_y (0 Mode)')

do i = 1, num_modes
   write (s1, '(I1.1)') i
   call file_q(2*i)%new(&
   &dataname = 'eyr'//trim(s1),&
   &label = 'Real Part of E_y ('//trim(s1)//' Mode')
   call file_q(2*i+1)%new(&
   &dataname = 'eyi'//trim(s1),&
   &label = 'Imaginary Part of E_y ('//trim(s1)//' Mode')
end do

call e%write_hdf5(file_q,2,8,8,id)

call file_q(1)%new(&
&dataname = 'ez',&
&label = 'E_z (0 Mode)')

do i = 1, num_modes
   write (s1, '(I1.1)') i
   call file_q(2*i)%new(&
   &dataname = 'ezr'//trim(s1),&
   &label = 'Real Part of E_z ('//trim(s1)//' Mode')
   call file_q(2*i+1)%new(&
   &dataname = 'ezi'//trim(s1),&
   &label = 'Imaginary Part of E_z ('//trim(s1)//' Mode')
end do

call e%write_hdf5(file_q,3,8,8,id)

call file_q(1)%new(&
&dataname = 'bx',&
&label = 'B_x (0 Mode)')

do i = 1, num_modes
   write (s1, '(I1.1)') i
   call file_q(2*i)%new(&
   &dataname = 'bxr'//trim(s1),&
   &label = 'Real Part of B_x ('//trim(s1)//' Mode')
   call file_q(2*i+1)%new(&
   &dataname = 'bxi'//trim(s1),&
   &label = 'Imaginary Part of B_x ('//trim(s1)//' Mode')
end do

call b%write_hdf5(file_q,1,8,8,id)

call file_q(1)%new(&
&dataname = 'by',&
&label = 'B_y (0 Mode)')

do i = 1, num_modes
   write (s1, '(I1.1)') i
   call file_q(2*i)%new(&
   &dataname = 'byr'//trim(s1),&
   &label = 'Real Part of B_y ('//trim(s1)//' Mode')
   call file_q(2*i+1)%new(&
   &dataname = 'byi'//trim(s1),&
   &label = 'Imaginary Part of B_y ('//trim(s1)//' Mode')
end do

call b%write_hdf5(file_q,2,8,8,id)

call file_q(1)%new(&
&dataname = 'bz',&
&label = 'B_z (0 Mode)')

do i = 1, num_modes
   write (s1, '(I1.1)') i
   call file_q(2*i)%new(&
   &dataname = 'bzr'//trim(s1),&
   &label = 'Real Part of B_z ('//trim(s1)//' Mode')
   call file_q(2*i+1)%new(&
   &dataname = 'bzi'//trim(s1),&
   &label = 'Imaginary Part of B_z ('//trim(s1)//' Mode')
end do

call b%write_hdf5(file_q,3,8,8,id)


call write_dbg( 'main', 'test_species', 0, 'ends' )

call end_errors()

call MPI_FINALIZE(ierr)


end program test_species
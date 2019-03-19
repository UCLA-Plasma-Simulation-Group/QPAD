program test_main

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
type(input_json), pointer :: input => null()
integer :: nr, nz, nrp, noff, xdim, npf, ierr, iter
integer :: num_modes, part_shape, i, id, j, k, nt
integer :: st, so ! smooth
character(len=:), allocatable :: shape, st_str
character(len=2) :: s1
real :: dr, dxi, rmin, rmax, zmin, zmax, dt, tt, prec
type fdist2d_wrap
   class(fdist2d), allocatable :: p
end type fdist2d_wrap
type(fdist2d_wrap) :: pf2d
type fdist3d_wrap
   class(fdist3d), allocatable :: p
end type fdist3d_wrap
type(fdist3d_wrap) :: pf3d
type(hdf5file) :: file_spe
type(hdf5file), dimension(:), allocatable :: file_q, file_qe, file_psi,&
&file_er,file_eth,file_ez,file_br,file_bth,file_bz
type(hdf5file) :: file_beam

type(field_rho), pointer :: pqb => null()
type(field_rho), target :: qe,qb
type(field_jay) :: cu
type(field_djdxi) :: dcu,acu,amu
type(field_b) :: b,bb,bt
type(field_e) :: e
type(field_psi) :: psi


! Initialization
allocate(pp,gp,input)
call input%new()
pp => input%pp

call write_dbg( 'main', 'test_species', 0, 'starts' )

call input%get('simulation.max_mode',num_modes)
allocate(file_q(2*num_modes+1))
allocate(file_qe(2*num_modes+1))
allocate(file_psi(2*num_modes+1))
allocate(file_er(2*num_modes+1))
allocate(file_eth(2*num_modes+1))
allocate(file_ez(2*num_modes+1))
allocate(file_br(2*num_modes+1))
allocate(file_bth(2*num_modes+1))
allocate(file_bz(2*num_modes+1))
call input%get('simulation.interpolation',shape)
call input%get('simulation.iter',iter)
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

call input%get('simulation.solver_precision',prec)
call input%get('simulation.smooth_type',st_str)
call input%get('simulation.smooth_order',so)

select case ( trim(st_str) )
case ( 'binomial' )
   st = p_smooth_binomial
case ( 'compensated' )
   st = p_smooth_compensated
case default
   st = p_smooth_none
end select

nrp = gp%get_ndp(1)
noff = gp%get_noff(1)

call qb%new(pp, gp, dr, dxi, num_modes, part_shape, st, so)
call qe%new(pp, gp, dr, dxi, num_modes, part_shape, st, so)
call cu%new(pp, gp, dr, dxi, num_modes, part_shape, st, so)
call amu%new(pp, gp, dr, dxi, num_modes, part_shape, st, so)
call dcu%new(pp, gp, dr, dxi, num_modes, part_shape, st, so)
call acu%new(pp, gp, dr, dxi, num_modes, part_shape, st, so)
call b%new( pp, gp, dr, dxi, num_modes, part_shape, entity=p_entity_plasma, iter_tol=prec )
call e%new( pp, gp, dr, dxi, num_modes, part_shape, entity=p_entity_plasma, iter_tol=prec )
call bb%new( pp, gp, dr, dxi, num_modes, part_shape, entity=p_entity_beam_old, iter_tol=prec )
call bt%new( pp, gp, dr, dxi, num_modes, part_shape, entity=p_entity_beam_old, iter_tol=prec )
call psi%new( pp, gp, dr, dxi, num_modes, part_shape, iter_tol=prec )
pqb => qb

xdim = 8

call input%get('species(1).profile',npf)
select case (npf)
case (0)
   allocate(fdist2d_000::pf2d%p)
   call pf2d%p%new(input,1)
end select
call spe%new(pp,gp,pf2d%p,part_shape,dr,dxi,num_modes,-1.0,dxi,xdim,0.0,st,so)

call input%get('beam(1).profile',npf)
select case (npf)
case (0)
   allocate(fdist3d_000::pf3d%p)
   call pf3d%p%new(input,1)
end select

call input%get('simulation.dt',dt)
call input%get('simulation.time',tt)
nt = tt/dt
call beam%new(pp,pqb,gp,part_shape,pf3d%p,-1.0,dt,7,st,so)

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

call file_qe(1)%new(&
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
   call file_qe(2*i)%new(&
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
   call file_qe(2*i+1)%new(&
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

call file_psi(1)%new(&
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
   call file_psi(2*i)%new(&
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
   call file_psi(2*i+1)%new(&
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

call file_er(1)%new(&
&axismin = (/rmin,zmin,0.0/),&
&axismax = (/rmax,zmax,1.0/),&
&axisname  = (/'r  ','\xi','z  '/),&
&axislabel = (/'r  ','\xi','z  '/),&
&timeunit = '1 / \omega_p',&
&dt = dxi,&
&axisunits = (/'c / \omega_p','c / \omega_p','c / \omega_p'/),&
&rank = 2,&
&filename = './',&
&dataname = 'Er',&
&units = 'mc\omega_p/e',&
&label = 'E_r (0 Mode)',&
&n = 0,&
&t = 0.0)


do i = 1, num_modes
   write (s1, '(I1.1)') i
   call file_er(2*i)%new(&
   &axismin = (/rmin,zmin,0.0/),&
   &axismax = (/rmax,zmax,1.0/),&
   &axisname  = (/'r  ','\xi','z  '/),&
   &axislabel = (/'r  ','\xi','z  '/),&
   &timeunit = '1 / \omega_p',&
   &dt = dxi,&
   &axisunits = (/'c / \omega_p','c / \omega_p','c / \omega_p'/),&
   &rank = 2,&
   &filename = './',&
   &dataname = 'Er_r'//trim(s1),&
   &units = 'mc\omega_p/e',&
   &label = 'Real Part of E_r ('//trim(s1)//' Mode',&
   &n = 0,&
   &t = 0.0)
   call file_er(2*i+1)%new(&
   &axismin = (/rmin,zmin,0.0/),&
   &axismax = (/rmax,zmax,1.0/),&
   &axisname  = (/'r  ','\xi','z  '/),&
   &axislabel = (/'r  ','\xi','z  '/),&
   &timeunit = '1 / \omega_p',&
   &dt = dxi,&
   &axisunits = (/'c / \omega_p','c / \omega_p','c / \omega_p'/),&
   &rank = 2,&
   &filename = './',&
   &dataname = 'Er_i'//trim(s1),&
   &units = 'mc\omega_p/e',&
   &label = 'Imaginary Part of E_r ('//trim(s1)//' Mode',&
   &n = 0,&
   &t = 0.0)
end do

call file_eth(1)%new(&
&axismin = (/rmin,zmin,0.0/),&
&axismax = (/rmax,zmax,1.0/),&
&axisname  = (/'r  ','\xi','z  '/),&
&axislabel = (/'r  ','\xi','z  '/),&
&timeunit = '1 / \omega_p',&
&dt = dxi,&
&axisunits = (/'c / \omega_p','c / \omega_p','c / \omega_p'/),&
&rank = 2,&
&filename = './',&
&dataname = 'Eth',&
&units = 'mc\omega_p/e',&
&label = 'E_\theta (0 Mode)',&
&n = 0,&
&t = 0.0)


do i = 1, num_modes
   write (s1, '(I1.1)') i
   call file_eth(2*i)%new(&
   &axismin = (/rmin,zmin,0.0/),&
   &axismax = (/rmax,zmax,1.0/),&
   &axisname  = (/'r  ','\xi','z  '/),&
   &axislabel = (/'r  ','\xi','z  '/),&
   &timeunit = '1 / \omega_p',&
   &dt = dxi,&
   &axisunits = (/'c / \omega_p','c / \omega_p','c / \omega_p'/),&
   &rank = 2,&
   &filename = './',&
   &dataname = 'Eth_r'//trim(s1),&
   &units = 'mc\omega_p/e',&
   &label = 'Real Part of E_\theta ('//trim(s1)//' Mode',&
   &n = 0,&
   &t = 0.0)
   call file_eth(2*i+1)%new(&
   &axismin = (/rmin,zmin,0.0/),&
   &axismax = (/rmax,zmax,1.0/),&
   &axisname  = (/'r  ','\xi','z  '/),&
   &axislabel = (/'r  ','\xi','z  '/),&
   &timeunit = '1 / \omega_p',&
   &dt = dxi,&
   &axisunits = (/'c / \omega_p','c / \omega_p','c / \omega_p'/),&
   &rank = 2,&
   &filename = './',&
   &dataname = 'Eth_i'//trim(s1),&
   &units = 'mc\omega_p/e',&
   &label = 'Imaginary Part of E_th ('//trim(s1)//' Mode',&
   &n = 0,&
   &t = 0.0)
end do

call file_ez(1)%new(&
&axismin = (/rmin,zmin,0.0/),&
&axismax = (/rmax,zmax,1.0/),&
&axisname  = (/'r  ','\xi','z  '/),&
&axislabel = (/'r  ','\xi','z  '/),&
&timeunit = '1 / \omega_p',&
&dt = dxi,&
&axisunits = (/'c / \omega_p','c / \omega_p','c / \omega_p'/),&
&rank = 2,&
&filename = './',&
&dataname = 'Ez',&
&units = 'mc\omega_p/e',&
&label = 'E_z (0 Mode)',&
&n = 0,&
&t = 0.0)


do i = 1, num_modes
   write (s1, '(I1.1)') i
   call file_ez(2*i)%new(&
   &axismin = (/rmin,zmin,0.0/),&
   &axismax = (/rmax,zmax,1.0/),&
   &axisname  = (/'r  ','\xi','z  '/),&
   &axislabel = (/'r  ','\xi','z  '/),&
   &timeunit = '1 / \omega_p',&
   &dt = dxi,&
   &axisunits = (/'c / \omega_p','c / \omega_p','c / \omega_p'/),&
   &rank = 2,&
   &filename = './',&
   &dataname = 'Ez_r'//trim(s1),&
   &units = 'mc\omega_p/e',&
   &label = 'Real Part of E_z ('//trim(s1)//' Mode',&
   &n = 0,&
   &t = 0.0)
   call file_ez(2*i+1)%new(&
   &axismin = (/rmin,zmin,0.0/),&
   &axismax = (/rmax,zmax,1.0/),&
   &axisname  = (/'r  ','\xi','z  '/),&
   &axislabel = (/'r  ','\xi','z  '/),&
   &timeunit = '1 / \omega_p',&
   &dt = dxi,&
   &axisunits = (/'c / \omega_p','c / \omega_p','c / \omega_p'/),&
   &rank = 2,&
   &filename = './',&
   &dataname = 'Ez_i'//trim(s1),&
   &units = 'mc\omega_p/e',&
   &label = 'Imaginary Part of E_z ('//trim(s1)//' Mode',&
   &n = 0,&
   &t = 0.0)
end do

call file_br(1)%new(&
&axismin = (/rmin,zmin,0.0/),&
&axismax = (/rmax,zmax,1.0/),&
&axisname  = (/'r  ','\xi','z  '/),&
&axislabel = (/'r  ','\xi','z  '/),&
&timeunit = '1 / \omega_p',&
&dt = dxi,&
&axisunits = (/'c / \omega_p','c / \omega_p','c / \omega_p'/),&
&rank = 2,&
&filename = './',&
&dataname = 'Br',&
&units = 'mc\omega_p/e',&
&label = 'B_r (0 Mode)',&
&n = 0,&
&t = 0.0)


do i = 1, num_modes
   write (s1, '(I1.1)') i
   call file_br(2*i)%new(&
   &axismin = (/rmin,zmin,0.0/),&
   &axismax = (/rmax,zmax,1.0/),&
   &axisname  = (/'r  ','\xi','z  '/),&
   &axislabel = (/'r  ','\xi','z  '/),&
   &timeunit = '1 / \omega_p',&
   &dt = dxi,&
   &axisunits = (/'c / \omega_p','c / \omega_p','c / \omega_p'/),&
   &rank = 2,&
   &filename = './',&
   &dataname = 'Br_r'//trim(s1),&
   &units = 'mc\omega_p/e',&
   &label = 'Real Part of B_r ('//trim(s1)//' Mode',&
   &n = 0,&
   &t = 0.0)
   call file_br(2*i+1)%new(&
   &axismin = (/rmin,zmin,0.0/),&
   &axismax = (/rmax,zmax,1.0/),&
   &axisname  = (/'r  ','\xi','z  '/),&
   &axislabel = (/'r  ','\xi','z  '/),&
   &timeunit = '1 / \omega_p',&
   &dt = dxi,&
   &axisunits = (/'c / \omega_p','c / \omega_p','c / \omega_p'/),&
   &rank = 2,&
   &filename = './',&
   &dataname = 'Br_i'//trim(s1),&
   &units = 'mc\omega_p/e',&
   &label = 'Imaginary Part of B_r ('//trim(s1)//' Mode',&
   &n = 0,&
   &t = 0.0)
end do

call file_bth(1)%new(&
&axismin = (/rmin,zmin,0.0/),&
&axismax = (/rmax,zmax,1.0/),&
&axisname  = (/'r  ','\xi','z  '/),&
&axislabel = (/'r  ','\xi','z  '/),&
&timeunit = '1 / \omega_p',&
&dt = dxi,&
&axisunits = (/'c / \omega_p','c / \omega_p','c / \omega_p'/),&
&rank = 2,&
&filename = './',&
&dataname = 'Bth',&
&units = 'mc\omega_p/e',&
&label = 'B_\theta (0 Mode)',&
&n = 0,&
&t = 0.0)


do i = 1, num_modes
   write (s1, '(I1.1)') i
   call file_bth(2*i)%new(&
   &axismin = (/rmin,zmin,0.0/),&
   &axismax = (/rmax,zmax,1.0/),&
   &axisname  = (/'r  ','\xi','z  '/),&
   &axislabel = (/'r  ','\xi','z  '/),&
   &timeunit = '1 / \omega_p',&
   &dt = dxi,&
   &axisunits = (/'c / \omega_p','c / \omega_p','c / \omega_p'/),&
   &rank = 2,&
   &filename = './',&
   &dataname = 'Bth_r'//trim(s1),&
   &units = 'mc\omega_p/e',&
   &label = 'Real Part of B_\theta ('//trim(s1)//' Mode',&
   &n = 0,&
   &t = 0.0)
   call file_bth(2*i+1)%new(&
   &axismin = (/rmin,zmin,0.0/),&
   &axismax = (/rmax,zmax,1.0/),&
   &axisname  = (/'r  ','\xi','z  '/),&
   &axislabel = (/'r  ','\xi','z  '/),&
   &timeunit = '1 / \omega_p',&
   &dt = dxi,&
   &axisunits = (/'c / \omega_p','c / \omega_p','c / \omega_p'/),&
   &rank = 2,&
   &filename = './',&
   &dataname = 'Bth_i'//trim(s1),&
   &units = 'mc\omega_p/e',&
   &label = 'Imaginary Part of B_th ('//trim(s1)//' Mode',&
   &n = 0,&
   &t = 0.0)
end do

call file_bz(1)%new(&
&axismin = (/rmin,zmin,0.0/),&
&axismax = (/rmax,zmax,1.0/),&
&axisname  = (/'r  ','\xi','z  '/),&
&axislabel = (/'r  ','\xi','z  '/),&
&timeunit = '1 / \omega_p',&
&dt = dxi,&
&axisunits = (/'c / \omega_p','c / \omega_p','c / \omega_p'/),&
&rank = 2,&
&filename = './',&
&dataname = 'Bz',&
&units = 'mc\omega_p/e',&
&label = 'B_z (0 Mode)',&
&n = 0,&
&t = 0.0)


do i = 1, num_modes
   write (s1, '(I1.1)') i
   call file_bz(2*i)%new(&
   &axismin = (/rmin,zmin,0.0/),&
   &axismax = (/rmax,zmax,1.0/),&
   &axisname  = (/'r  ','\xi','z  '/),&
   &axislabel = (/'r  ','\xi','z  '/),&
   &timeunit = '1 / \omega_p',&
   &dt = dxi,&
   &axisunits = (/'c / \omega_p','c / \omega_p','c / \omega_p'/),&
   &rank = 2,&
   &filename = './',&
   &dataname = 'Bz_r'//trim(s1),&
   &units = 'mc\omega_p/e',&
   &label = 'Real Part of B_z ('//trim(s1)//' Mode',&
   &n = 0,&
   &t = 0.0)
   call file_bz(2*i+1)%new(&
   &axismin = (/rmin,zmin,0.0/),&
   &axismax = (/rmax,zmax,1.0/),&
   &axisname  = (/'r  ','\xi','z  '/),&
   &axislabel = (/'r  ','\xi','z  '/),&
   &timeunit = '1 / \omega_p',&
   &dt = dxi,&
   &axisunits = (/'c / \omega_p','c / \omega_p','c / \omega_p'/),&
   &rank = 2,&
   &filename = './',&
   &dataname = 'Bz_i'//trim(s1),&
   &units = 'mc\omega_p/e',&
   &label = 'Imaginary Part of B_z ('//trim(s1)//' Mode',&
   &n = 0,&
   &t = 0.0)
end do

do i = 1, nt

   call start_tprof( 'total simulation time' )
   call qb%as(0.0)
   call beam%qdp(qb)
      do j = 1, nz
      ! call file_spe%new(&
      ! &timeunit = '1 / \omega_p',&
      ! &dt = dxi,&
      ! &ty = 'particles',&
      ! &filename = './spe_',&
      ! &dataname = 'raw',&
      ! &units = '',&
      ! &label = 'Plasma Raw',&
      ! &n = j,&
      ! &t = real(j))
      ! call spe%wr(file_spe)

      call qb%copy_slice(j, p_copy_2to1)
      call qb%smooth_f1()
      call bb%solve_old(qb)
      qe = 0.0
      call spe%qdp(qe)
      call qe%smooth_f1()
      call psi%solve(qe)
      call spe%extpsi(psi)
      call cu%smooth_f1()
      call e%solve(cu)
      call b%solve(cu)
      do k = 1, iter
         bt = b + bb
         call e%solve(bt,psi)
         cu = 0.0
         acu = 0.0
         amu = 0.0
         call spe%amjdp(e,bt,cu,amu,acu)
         call acu%smooth_f1()
         call amu%smooth_f1()
         call cu%smooth_f1()
         call dcu%solve(acu,amu)
         call b%solve(dcu,cu)
         call e%solve(cu)
         call b%solve(cu)
         if (k == iter) then
            call spe%cbq(j+1)
         end if
      end do
      bb = b + bb
      call e%solve(bb,psi)
      cu = cu + dcu
      call spe%push(e,bb)
      call e%copy_slice(j+1, p_copy_1to2)
      call bb%copy_slice(j+1, p_copy_1to2)
      call psi%copy_slice(j+1, p_copy_1to2)
      write (2,*) "Step", j
   end do

   call beam%push(e,bb,7,7,id)
   call spe%renew(i*dt)

   call stop_tprof( 'total simulation time' )

   call file_q(1)%new(&
   &n = i,&
   &t = i*dt)
   
   do j = 1, num_modes
      call file_q(2*j)%new(&
      &n = i,&
      &t = i*dt)
      call file_q(2*j+1)%new(&
      &n = i,&
      &t = i*dt)
   end do   
   call beam%wrq(file_q,8,8,id)

   call file_beam%new(&
   &timeunit = '1 / \omega_p',&
   &dt = dt,&
   &ty = 'particles',&
   &filename = './beam_',&
   &dataname = 'raw',&
   &units = '',&
   &label = 'Beam Raw',&
   &n = i,&
   &t = i*dt)
   call beam%wr(file_beam,1,6,6,id)

   call file_qe(1)%new(&
   &n = i,&
   &t = i*dt)
   
   do j = 1, num_modes
      call file_qe(2*j)%new(&
      &n = i,&
      &t = i*dt)
      call file_qe(2*j+1)%new(&
      &n = i,&
      &t = i*dt)
   end do   
   call spe%wrq(file_qe,5,5,id)

   call file_psi(1)%new(&
   &n = i,&
   &t = i*dt)
   
   do j = 1, num_modes
      call file_psi(2*j)%new(&
      &n = i,&
      &t = i*dt)
      call file_psi(2*j+1)%new(&
      &n = i,&
      &t = i*dt)
   end do   
   call psi%write_hdf5(file_psi,1,8,8,id)

   call file_er(1)%new(&
   &n = i,&
   &t = i*dt)
   
   do j = 1, num_modes
      call file_er(2*j)%new(&
      &n = i,&
      &t = i*dt)
      call file_er(2*j+1)%new(&
      &n = i,&
      &t = i*dt)
   end do   
   call e%write_hdf5(file_er,1,8,8,id)

   call file_eth(1)%new(&
   &n = i,&
   &t = i*dt)
   
   do j = 1, num_modes
      call file_eth(2*j)%new(&
      &n = i,&
      &t = i*dt)
      call file_eth(2*j+1)%new(&
      &n = i,&
      &t = i*dt)
   end do   
   call e%write_hdf5(file_eth,2,8,8,id)

   call file_ez(1)%new(&
   &n = i,&
   &t = i*dt)
   
   do j = 1, num_modes
      call file_ez(2*j)%new(&
      &n = i,&
      &t = i*dt)
      call file_ez(2*j+1)%new(&
      &n = i,&
      &t = i*dt)
   end do   
   call e%write_hdf5(file_ez,3,8,8,id)

   call file_br(1)%new(&
   &n = i,&
   &t = i*dt)
   
   do j = 1, num_modes
      call file_br(2*j)%new(&
      &n = i,&
      &t = i*dt)
      call file_br(2*j+1)%new(&
      &n = i,&
      &t = i*dt)
   end do   
   call bb%write_hdf5(file_br,1,8,8,id)

   call file_bth(1)%new(&
   &n = i,&
   &t = i*dt)
   
   do j = 1, num_modes
      call file_bth(2*j)%new(&
      &n = i,&
      &t = i*dt)
      call file_bth(2*j+1)%new(&
      &n = i,&
      &t = i*dt)
   end do   
   call bb%write_hdf5(file_bth,2,8,8,id)

   call file_bz(1)%new(&
   &n = i,&
   &t = i*dt)
   
   do j = 1, num_modes
      call file_bz(2*j)%new(&
      &n = i,&
      &t = i*dt)
      call file_bz(2*j+1)%new(&
      &n = i,&
      &t = i*dt)
   end do   
   call bb%write_hdf5(file_bz,3,8,8,id)

end do

call write_tprof()
call write_dbg( 'main', 'test_main', 0, 'ends' )

call end_errors()

call MPI_FINALIZE(ierr)


end program test_main
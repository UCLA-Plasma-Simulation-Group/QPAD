program test_main

use parallel_pipe_class
use grid_class
use field_class
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
type(beam3d) :: beam1, beam2
class(parallel_pipe), pointer :: pp => null()
class(grid), pointer :: gp => null()
type(input_json), pointer :: input => null()
integer :: nr, nz, nrp, noff, xdim, npf, ierr, iter
integer :: num_modes, part_shape, fld_bnd, i, id, j, k, nt
integer :: st, so ! smooth
character(len=:), allocatable :: shape, st_str, fld_bnd_str
character(len=2) :: s1
real :: dr, dxi, rmin, rmax, zmin, zmax, dt, tt, prec
type fdist2d_wrap
   class(fdist2d), allocatable :: p
end type fdist2d_wrap
type(fdist2d_wrap) :: pf2d
type fdist3d_wrap
   class(fdist3d), allocatable :: p
end type fdist3d_wrap
type(fdist3d_wrap) :: pf3d1, pf3d2
type(hdf5file) :: file_spe
type(hdf5file), dimension(:), allocatable :: file_rjz, file_q1, file_qe, file_psi,&
&file_er,file_eth,file_ez,file_br,file_bth,file_bz, file_q2,file_jr,file_jth,file_jz
type(hdf5file) :: file_beam1, file_beam2

class(field), pointer :: pqb => null()
type(field_rho), target :: qe,qb
type(field_jay) :: cu, amu
type(field_djdxi) :: dcu,acu
type(field_b) :: b,bb,bt
type(field_e) :: e
type(field_psi) :: psi
integer :: ndump = 1


! Initialization
allocate(pp,gp,input)
call input%new()
pp => input%pp

call write_dbg( 'main', 'test_species', 0, 'starts' )

call input%get('simulation.max_mode',num_modes)
allocate(file_q1(2*num_modes+1))
allocate(file_q2(2*num_modes+1))
allocate(file_qe(2*num_modes+1))
allocate(file_rjz(2*num_modes+1))
allocate(file_psi(2*num_modes+1))
allocate(file_er(2*num_modes+1))
allocate(file_eth(2*num_modes+1))
allocate(file_ez(2*num_modes+1))
allocate(file_br(2*num_modes+1))
allocate(file_bth(2*num_modes+1))
allocate(file_bz(2*num_modes+1))
allocate(file_jr(2*num_modes+1))
allocate(file_jth(2*num_modes+1))
allocate(file_jz(2*num_modes+1))
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

call input%get( 'simulation.field_boundary', fld_bnd_str )
  select case ( trim(fld_bnd_str) )
  case ( 'zero' )
    fld_bnd = p_bnd_zero
  case ( 'conduct' )
    fld_bnd = p_bnd_conduct
  case ( 'open' )
    fld_bnd = p_bnd_open
  case default
    call write_err( 'Invalid field boundary type!' )
  end select

nrp = gp%get_ndp(1)
noff = gp%get_noff(1)

call qb%new(pp, gp, dr, dxi, num_modes, part_shape, st, so)
call qe%new(pp, gp, dr, dxi, num_modes, part_shape, st, so)
call cu%new(pp, gp, dr, dxi, num_modes, part_shape, st, so)
call amu%new(pp, gp, dr, dxi, num_modes, part_shape, st, so)
call dcu%new(pp, gp, dr, dxi, num_modes, part_shape, st, so)
call acu%new(pp, gp, dr, dxi, num_modes, part_shape, st, so)
call b%new( pp, gp, dr, dxi, num_modes, part_shape, fld_bnd, entity=p_entity_plasma, iter_tol=prec )
call e%new( pp, gp, dr, dxi, num_modes, part_shape, fld_bnd, entity=p_entity_plasma, iter_tol=prec )
call bb%new( pp, gp, dr, dxi, num_modes, part_shape, fld_bnd, entity=p_entity_beam, iter_tol=prec )
call bt%new( pp, gp, dr, dxi, num_modes, part_shape, fld_bnd, entity=p_entity_plasma, iter_tol=prec )
call psi%new( pp, gp, dr, dxi, num_modes, part_shape, fld_bnd, iter_tol=prec )

pqb => qb

xdim = 8

call input%get('species(1).profile',npf)
select case (npf)
case (0)
   allocate(fdist2d_000::pf2d%p)
   call pf2d%p%new(input,1)
case (12)
   allocate(fdist2d_012::pf2d%p)
   call pf2d%p%new(input,1)
end select
call spe%new(pp,gp,pf2d%p,part_shape,dr,dxi,num_modes,-1.0,dxi,xdim,0.0,st,so)

call input%get('beam(1).profile',npf)
select case (npf)
case (0)
   allocate(fdist3d_000::pf3d1%p)
   call pf3d1%p%new(input,1)
end select

call input%get('beam(2).profile',npf)
select case (npf)
case (0)
   allocate(fdist3d_000::pf3d2%p)
   call pf3d2%p%new(input,2)
end select

call input%get('simulation.dt',dt)
call input%get('simulation.time',tt)
nt = tt/dt
call beam1%new(pp,pqb,gp,part_shape,pf3d1%p,-1.0,dt,7,st,so)
call beam2%new(pp,pqb,gp,part_shape,pf3d2%p,-1.0,dt,7,st,so)

if (pp%getidproc() == 0) then
   call system('mkdir -p ./Beam01/Charge/')
   call system('mkdir -p ./Beam02/Charge/')
   call system('mkdir -p ./Beam01/Raw/')
   call system('mkdir -p ./Beam02/Raw/')
   call system('mkdir -p ./Species01/Charge/')
   call system('mkdir -p ./Species01/Raw/')
   call system('mkdir -p ./Fields/Rhojz/')
   call system('mkdir -p ./Fields/Psi/')
   call system('mkdir -p ./Fields/Er/')
   call system('mkdir -p ./Fields/Eth/')
   call system('mkdir -p ./Fields/Ez/')
   call system('mkdir -p ./Fields/Br/')
   call system('mkdir -p ./Fields/Bth/')
   call system('mkdir -p ./Fields/Bz/')
   call system('mkdir -p ./Fields/Jr/')
   call system('mkdir -p ./Fields/Jth/')
   call system('mkdir -p ./Fields/Jz/')
end if

call MPI_BARRIER(pp%getlworld(),ierr)


call file_q1(1)%new(&
&axismin = (/rmin,zmin,0.0/),&
&axismax = (/rmax,zmax,1.0/),&
&axisname  = (/'r  ','\xi','z  '/),&
&axislabel = (/'r  ','\xi','z  '/),&
&timeunit = '1 / \omega_p',&
&dt = dxi,&
&axisunits = (/'c / \omega_p','c / \omega_p','c / \omega_p'/),&
&rank = 2,&
&filename = './Beam01/Charge/',&
&dataname = 'qb0',&
&units = 'n_0',&
&label = 'Beam Charge Density',&
&n = 0,&
&t = 0.0)

do i = 1, num_modes
   write (s1, '(I1.1)') i
   call file_q1(2*i)%new(&
   &axismin = (/rmin,zmin,0.0/),&
   &axismax = (/rmax,zmax,1.0/),&
   &axisname  = (/'r  ','\xi','z  '/),&
   &axislabel = (/'r  ','\xi','z  '/),&
   &timeunit = '1 / \omega_p',&
   &dt = dxi,&
   &axisunits = (/'c / \omega_p','c / \omega_p','c / \omega_p'/),&
   &rank = 2,&
   &filename = './Beam01/Charge/',&
   &dataname = 'qbr'//trim(s1),&
   &units = 'n_0',&
   &label = 'Beam Density',&
   &n = 0,&
   &t = 0.0)
   call file_q1(2*i+1)%new(&
   &axismin = (/rmin,zmin,0.0/),&
   &axismax = (/rmax,zmax,1.0/),&
   &axisname  = (/'r  ','\xi','z  '/),&
   &axislabel = (/'r  ','\xi','z  '/),&
   &timeunit = '1 / \omega_p',&
   &dt = dxi,&
   &axisunits = (/'c / \omega_p','c / \omega_p','c / \omega_p'/),&
   &rank = 2,&
   &filename = './Beam01/Charge/',&
   &dataname = 'qbi'//trim(s1),&
   &units = 'n_0',&
   &label = 'Beam Density',&
   &n = 0,&
   &t = 0.0)
end do

call file_q2(1)%new(&
&axismin = (/rmin,zmin,0.0/),&
&axismax = (/rmax,zmax,1.0/),&
&axisname  = (/'r  ','\xi','z  '/),&
&axislabel = (/'r  ','\xi','z  '/),&
&timeunit = '1 / \omega_p',&
&dt = dxi,&
&axisunits = (/'c / \omega_p','c / \omega_p','c / \omega_p'/),&
&rank = 2,&
&filename = './Beam02/Charge/',&
&dataname = 'qb0',&
&units = 'n_0',&
&label = 'Beam Charge Density',&
&n = 0,&
&t = 0.0)

do i = 1, num_modes
   write (s1, '(I1.1)') i
   call file_q2(2*i)%new(&
   &axismin = (/rmin,zmin,0.0/),&
   &axismax = (/rmax,zmax,1.0/),&
   &axisname  = (/'r  ','\xi','z  '/),&
   &axislabel = (/'r  ','\xi','z  '/),&
   &timeunit = '1 / \omega_p',&
   &dt = dxi,&
   &axisunits = (/'c / \omega_p','c / \omega_p','c / \omega_p'/),&
   &rank = 2,&
   &filename = './Beam02/Charge/',&
   &dataname = 'qbr'//trim(s1),&
   &units = 'n_0',&
   &label = 'Beam Density',&
   &n = 0,&
   &t = 0.0)
   call file_q2(2*i+1)%new(&
   &axismin = (/rmin,zmin,0.0/),&
   &axismax = (/rmax,zmax,1.0/),&
   &axisname  = (/'r  ','\xi','z  '/),&
   &axislabel = (/'r  ','\xi','z  '/),&
   &timeunit = '1 / \omega_p',&
   &dt = dxi,&
   &axisunits = (/'c / \omega_p','c / \omega_p','c / \omega_p'/),&
   &rank = 2,&
   &filename = './Beam02/Charge/',&
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
&filename = './Species01/Charge/',&
&dataname = 'q0',&
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
   &filename = './Species01/Charge/',&
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
   &filename = './Species01/Charge/',&
   &dataname = 'qi'//trim(s1),&
   &units = 'n_0',&
   &label = 'Charge Density',&
   &n = 0,&
   &t = 0.0)
end do

call file_rjz(1)%new(&
&axismin = (/rmin,zmin,0.0/),&
&axismax = (/rmax,zmax,1.0/),&
&axisname  = (/'r  ','\xi','z  '/),&
&axislabel = (/'r  ','\xi','z  '/),&
&timeunit = '1 / \omega_p',&
&dt = dxi,&
&axisunits = (/'c / \omega_p','c / \omega_p','c / \omega_p'/),&
&rank = 2,&
&filename = './Fields/Rhojz/',&
&dataname = 'rjz0',&
&units = 'n_0',&
&label = '\rho - J_z',&
&n = 0,&
&t = 0.0)


do i = 1, num_modes
   write (s1, '(I1.1)') i
   call file_rjz(2*i)%new(&
   &axismin = (/rmin,zmin,0.0/),&
   &axismax = (/rmax,zmax,1.0/),&
   &axisname  = (/'r  ','\xi','z  '/),&
   &axislabel = (/'r  ','\xi','z  '/),&
   &timeunit = '1 / \omega_p',&
   &dt = dxi,&
   &axisunits = (/'c / \omega_p','c / \omega_p','c / \omega_p'/),&
   &rank = 2,&
   &filename = './Fields/Rhojz/',&
   &dataname = 'rjz'//trim(s1),&
   &units = 'n_0',&
   &label = '\rho - J_z',&
   &n = 0,&
   &t = 0.0)
   call file_rjz(2*i+1)%new(&
   &axismin = (/rmin,zmin,0.0/),&
   &axismax = (/rmax,zmax,1.0/),&
   &axisname  = (/'r  ','\xi','z  '/),&
   &axislabel = (/'r  ','\xi','z  '/),&
   &timeunit = '1 / \omega_p',&
   &dt = dxi,&
   &axisunits = (/'c / \omega_p','c / \omega_p','c / \omega_p'/),&
   &rank = 2,&
   &filename = './Fields/Rhojz/',&
   &dataname = 'rjz'//trim(s1),&
   &units = 'n_0',&
   &label = '\rho - J_z',&
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
&filename = './Fields/Psi/',&
&dataname = 'psi0',&
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
   &filename = './Fields/Psi/',&
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
   &filename = './Fields/Psi/',&
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
&filename = './Fields/Er/',&
&dataname = 'Er0',&
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
   &filename = './Fields/Er/',&
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
   &filename = './Fields/Er/',&
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
&filename = './Fields/Eth/',&
&dataname = 'Eth0',&
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
   &filename = './Fields/Eth/',&
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
   &filename = './Fields/Eth/',&
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
&filename = './Fields/Ez/',&
&dataname = 'Ez0',&
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
   &filename = './Fields/Ez/',&
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
   &filename = './Fields/Ez/',&
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
&filename = './Fields/Br/',&
&dataname = 'Br0',&
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
   &filename = './Fields/Br/',&
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
   &filename = './Fields/Br/',&
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
&filename = './Fields/Bth/',&
&dataname = 'Bth0',&
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
   &filename = './Fields/Bth/',&
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
   &filename = './Fields/Bth/',&
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
&filename = './Fields/Bz/',&
&dataname = 'Bz0',&
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
   &filename = './Fields/Bz/',&
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
   &filename = './Fields/Bz/',&
   &dataname = 'Bz_i'//trim(s1),&
   &units = 'mc\omega_p/e',&
   &label = 'Imaginary Part of B_z ('//trim(s1)//' Mode',&
   &n = 0,&
   &t = 0.0)
end do

call file_jr(1)%new(&
&axismin = (/rmin,zmin,0.0/),&
&axismax = (/rmax,zmax,1.0/),&
&axisname  = (/'r  ','\xi','z  '/),&
&axislabel = (/'r  ','\xi','z  '/),&
&timeunit = '1 / \omega_p',&
&dt = dxi,&
&axisunits = (/'c / \omega_p','c / \omega_p','c / \omega_p'/),&
&rank = 2,&
&filename = './Fields/Jr/',&
&dataname = 'Jr0',&
&units = 'n_pec',&
&label = 'J_r (0 Mode)',&
&n = 0,&
&t = 0.0)


do i = 1, num_modes
   write (s1, '(I1.1)') i
   call file_jr(2*i)%new(&
   &axismin = (/rmin,zmin,0.0/),&
   &axismax = (/rmax,zmax,1.0/),&
   &axisname  = (/'r  ','\xi','z  '/),&
   &axislabel = (/'r  ','\xi','z  '/),&
   &timeunit = '1 / \omega_p',&
   &dt = dxi,&
   &axisunits = (/'c / \omega_p','c / \omega_p','c / \omega_p'/),&
   &rank = 2,&
   &filename = './Fields/Jr/',&
   &dataname = 'Jr_r'//trim(s1),&
   &units = 'n_pec',&
   &label = 'Real Part of J_r ('//trim(s1)//' Mode',&
   &n = 0,&
   &t = 0.0)
   call file_jr(2*i+1)%new(&
   &axismin = (/rmin,zmin,0.0/),&
   &axismax = (/rmax,zmax,1.0/),&
   &axisname  = (/'r  ','\xi','z  '/),&
   &axislabel = (/'r  ','\xi','z  '/),&
   &timeunit = '1 / \omega_p',&
   &dt = dxi,&
   &axisunits = (/'c / \omega_p','c / \omega_p','c / \omega_p'/),&
   &rank = 2,&
   &filename = './Fields/Jr/',&
   &dataname = 'Jr_i'//trim(s1),&
   &units = 'n_pec',&
   &label = 'Imaginary Part of J_r ('//trim(s1)//' Mode',&
   &n = 0,&
   &t = 0.0)
end do

call file_jth(1)%new(&
&axismin = (/rmin,zmin,0.0/),&
&axismax = (/rmax,zmax,1.0/),&
&axisname  = (/'r  ','\xi','z  '/),&
&axislabel = (/'r  ','\xi','z  '/),&
&timeunit = '1 / \omega_p',&
&dt = dxi,&
&axisunits = (/'c / \omega_p','c / \omega_p','c / \omega_p'/),&
&rank = 2,&
&filename = './Fields/Jth/',&
&dataname = 'Jth0',&
&units = 'n_pec',&
&label = 'J_\theta (0 Mode)',&
&n = 0,&
&t = 0.0)


do i = 1, num_modes
   write (s1, '(I1.1)') i
   call file_jth(2*i)%new(&
   &axismin = (/rmin,zmin,0.0/),&
   &axismax = (/rmax,zmax,1.0/),&
   &axisname  = (/'r  ','\xi','z  '/),&
   &axislabel = (/'r  ','\xi','z  '/),&
   &timeunit = '1 / \omega_p',&
   &dt = dxi,&
   &axisunits = (/'c / \omega_p','c / \omega_p','c / \omega_p'/),&
   &rank = 2,&
   &filename = './Fields/Jth/',&
   &dataname = 'Jth_r'//trim(s1),&
   &units = 'n_pec',&
   &label = 'Real Part of J_\theta ('//trim(s1)//' Mode',&
   &n = 0,&
   &t = 0.0)
   call file_jth(2*i+1)%new(&
   &axismin = (/rmin,zmin,0.0/),&
   &axismax = (/rmax,zmax,1.0/),&
   &axisname  = (/'r  ','\xi','z  '/),&
   &axislabel = (/'r  ','\xi','z  '/),&
   &timeunit = '1 / \omega_p',&
   &dt = dxi,&
   &axisunits = (/'c / \omega_p','c / \omega_p','c / \omega_p'/),&
   &rank = 2,&
   &filename = './Fields/Jth/',&
   &dataname = 'Jth_i'//trim(s1),&
   &units = 'n_pec',&
   &label = 'Imaginary Part of J_\theta ('//trim(s1)//' Mode',&
   &n = 0,&
   &t = 0.0)
end do

call file_jz(1)%new(&
&axismin = (/rmin,zmin,0.0/),&
&axismax = (/rmax,zmax,1.0/),&
&axisname  = (/'r  ','\xi','z  '/),&
&axislabel = (/'r  ','\xi','z  '/),&
&timeunit = '1 / \omega_p',&
&dt = dxi,&
&axisunits = (/'c / \omega_p','c / \omega_p','c / \omega_p'/),&
&rank = 2,&
&filename = './Fields/Jz/',&
&dataname = 'Jz0',&
&units = 'n_pec',&
&label = 'J_z (0 Mode)',&
&n = 0,&
&t = 0.0)

do i = 1, num_modes
   write (s1, '(I1.1)') i
   call file_jz(2*i)%new(&
   &axismin = (/rmin,zmin,0.0/),&
   &axismax = (/rmax,zmax,1.0/),&
   &axisname  = (/'r  ','\xi','z  '/),&
   &axislabel = (/'r  ','\xi','z  '/),&
   &timeunit = '1 / \omega_p',&
   &dt = dxi,&
   &axisunits = (/'c / \omega_p','c / \omega_p','c / \omega_p'/),&
   &rank = 2,&
   &filename = './Fields/Jz/',&
   &dataname = 'Jz_r'//trim(s1),&
   &units = 'n_pec',&
   &label = 'Real Part of J_z ('//trim(s1)//' Mode',&
   &n = 0,&
   &t = 0.0)
   call file_jz(2*i+1)%new(&
   &axismin = (/rmin,zmin,0.0/),&
   &axismax = (/rmax,zmax,1.0/),&
   &axisname  = (/'r  ','\xi','z  '/),&
   &axislabel = (/'r  ','\xi','z  '/),&
   &timeunit = '1 / \omega_p',&
   &dt = dxi,&
   &axisunits = (/'c / \omega_p','c / \omega_p','c / \omega_p'/),&
   &rank = 2,&
   &filename = './Fields/Jz/',&
   &dataname = 'Jz_i'//trim(s1),&
   &units = 'n_pec',&
   &label = 'Imaginary Part of J_z ('//trim(s1)//' Mode',&
   &n = 0,&
   &t = 0.0)
end do

do i = 1, nt

   call start_tprof( 'total simulation time' )
   call qb%as(0.0)
   call qe%as(0.0)
   call beam1%qdp(qb)
   call beam2%qdp(qb)
   cu = 0.0
   b = 0.0
   e = 0.0
   bt = 0.0
   psi = 0.0
   do j = 1, nz

      call file_spe%new(&
      &timeunit = '1 / \omega_p',&
      &dt = dxi,&
      &ty = 'particles',&
      &filename = './Species01/Raw/',&
      &dataname = 'raw',&
      &units = '',&
      &label = 'Plasma Raw',&
      &n = j,&
      &t = real(j))
      call spe%wr(file_spe)

      call qb%copy_slice(j+1, p_copy_2to1)
      call qb%smooth_f1()
      call bb%solve(qb)
      qe = 0.0
      call spe%qdp(qe)
      call qe%smooth_f1()
      call qe%copy_slice(j+1, p_copy_1to2)
      call psi%solve(qe)
      call spe%extpsi(psi)
      call e%solve(psi,j+1)
      call bt%solve(cu)
      do k = 1, iter
         ! b = bt + bb
         call add_f1( bt, bb, b )
         call e%solve(b,psi)
         cu = 0.0
         acu = 0.0
         amu = 0.0
         call spe%amjdp(e,b,cu,amu,acu)
         call acu%smooth_f1()
         call amu%smooth_f1()
         call cu%smooth_f1()
         call dcu%solve(acu,amu)
         call bt%solve(dcu,cu)
         ! call e%solve(cu)
         call bt%solve(cu)
         if (k == iter) then
            call spe%cbq(j+1)
            call cu%copy_slice(j+1, p_copy_1to2)
         end if
      end do
      ! b = bt + bb
      call add_f1( bt, bb, b )
      call e%solve(b,psi)
      ! cu = cu + dcu*dxi
      call dot_f1( dxi, dcu )
      call add_f1( dcu, cu, (/1,2/), (/1,2/) )
      call spe%push(e,b)
      call e%copy_slice(j+1, p_copy_1to2)
      call b%copy_slice(j+1, p_copy_1to2)
      call psi%copy_slice(j+1, p_copy_1to2)
      write (2,*) "Step", j
   end do

   call beam1%push(e,b,7,7,id)
   call beam2%push(e,b,7,7,id)
   call spe%renew(i*dt)

   call stop_tprof( 'total simulation time' )

   if (mod ((i-1), ndump) == 0) then

      call file_rjz(1)%new(&
      &n = i,&
      &t = i*dt)
      
      do j = 1, num_modes
         call file_rjz(2*j)%new(&
         &n = i,&
         &t = i*dt)
         call file_rjz(2*j+1)%new(&
         &n = i,&
         &t = i*dt)
      end do   
      call qe%write_hdf5(file_rjz,1,8,8,id)

      call file_q1(1)%new(&
      &n = i,&
      &t = i*dt)
      
      do j = 1, num_modes
         call file_q1(2*j)%new(&
         &n = i,&
         &t = i*dt)
         call file_q1(2*j+1)%new(&
         &n = i,&
         &t = i*dt)
      end do   
      call beam1%wrq(file_q1,8,8,id)
   
      call file_q2(1)%new(&
      &n = i,&
      &t = i*dt)
      
      do j = 1, num_modes
         call file_q2(2*j)%new(&
         &n = i,&
         &t = i*dt)
         call file_q2(2*j+1)%new(&
         &n = i,&
         &t = i*dt)
      end do   
      call beam2%wrq(file_q2,8,8,id)
   
      call file_beam1%new(&
      &timeunit = '1 / \omega_p',&
      &dt = dt,&
      &ty = 'particles',&
      &filename = './Beam01/Raw/',&
      &dataname = 'raw01',&
      &units = '',&
      &label = 'Beam Raw',&
      &n = i,&
      &t = i*dt)
      call beam1%wr(file_beam1,1,6,6,id)
   
      call file_beam2%new(&
      &timeunit = '1 / \omega_p',&
      &dt = dt,&
      &ty = 'particles',&
      &filename = './Beam02/Raw/',&
      &dataname = 'raw02',&
      &units = '',&
      &label = 'Beam Raw',&
      &n = i,&
      &t = i*dt)
      call beam2%wr(file_beam2,1,6,6,id)
   
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
      call b%write_hdf5(file_br,1,8,8,id)
   
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
      call b%write_hdf5(file_bth,2,8,8,id)
   
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
      call b%write_hdf5(file_bz,3,8,8,id)

      call file_jr(1)%new(&
      &n = i,&
      &t = i*dt)
      
      do j = 1, num_modes
         call file_jr(2*j)%new(&
         &n = i,&
         &t = i*dt)
         call file_jr(2*j+1)%new(&
         &n = i,&
         &t = i*dt)
      end do   
      call cu%write_hdf5(file_jr,1,8,8,id)
   
      call file_jth(1)%new(&
      &n = i,&
      &t = i*dt)
      
      do j = 1, num_modes
         call file_jth(2*j)%new(&
         &n = i,&
         &t = i*dt)
         call file_jth(2*j+1)%new(&
         &n = i,&
         &t = i*dt)
      end do   
      call cu%write_hdf5(file_jth,2,8,8,id)
   
      call file_jz(1)%new(&
      &n = i,&
      &t = i*dt)
      
      do j = 1, num_modes
         call file_jz(2*j)%new(&
         &n = i,&
         &t = i*dt)
         call file_jz(2*j+1)%new(&
         &n = i,&
         &t = i*dt)
      end do   
      call cu%write_hdf5(file_jz,3,8,8,id)

   end if

end do

call write_tprof()
call write_dbg( 'main', 'test_main', 0, 'ends' )

call end_errors()

call MPI_FINALIZE(ierr)


end program test_main
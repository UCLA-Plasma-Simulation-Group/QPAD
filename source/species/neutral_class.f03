module neutral_class

use parallel_module
use param
use sysutil_module
use options_class
use fdist2d_class
use part2d_class
use part2d_comm
use field_class
use field_e_class
use field_b_class
use field_psi_class
use field_src_class
use hdf5io_class

implicit none

private

double precision, parameter, dimension(3,1) :: H_param = &
   reshape( (/8.50168d19, 3.42554d2, 1.0d0 /), &
                             (/ 3,1 /) )

double precision, parameter, dimension(3,2) :: He_param = &
   reshape( (/ 7.42273d18, 8.28266d2, 4.92d-1, & 
             2.7303d21, 2.74043d3, 1.00d0 /), &
                              (/ 3,2 /) )

double precision, parameter, dimension(3,3) :: Li_param = &
   reshape( (/ 3.51032d21, 8.54681d1, 2.18d0, & 
               3.64489d20, 4.49125d3, 6.9669d-1, &
               2.07194d22, 9.25111d3, 1.00025d0 /), &
                                (/ 3,3 /) )    

double precision, parameter, dimension(3,2) :: Cs_param = &
   reshape( (/ 1.00952125d22, 5.24016673d1, 2.74d0, & 
            4.36300414d23, 8.58877619d2, 1.94d0 /), &
                                (/ 3,2 /) )    

double precision, parameter, dimension(3,2) :: Rb_param = &
   reshape( (/ 8.1707423d21, 5.8422794d1, 2.6085756d0, & 
            2.5593350d23, 9.7516547d2, 1.8240161d0 /), &
                                (/ 3,2 /) )    

double precision, parameter, dimension(3,2) :: K_param = &
   reshape( (/ 7.19187444d21, 6.17526072d1, 2.54d0, & 
            9.15936246d22, 1.21498101d3, 1.62d0 /), &
                                (/ 3,2 /) )    

double precision, parameter, dimension(3,8) :: Ar_param = &
   reshape( (/ 4.58598d19, 4.27322d2, 8.58275d-1, & 
            2.29786d23, 9.96806d2, 1.80234d0, &
            1.86503d26, 1.78651d3, 2.46057d0, &
            2.21776d28, 3.15765d3, 2.81622d0, &
            3.38223d30, 4.43622d3, 3.25919d0, &
            2.96575d32, 5.95836d3, 3.63238d0, &
            2.23773d33, 9.4309d3,  3.63741d0, &
            9.00240d34, 1.17359d4, 3.92734d0 /), &
                                (/ 3,8 /) ) 

double precision, parameter, dimension(3,7) :: N_param = &
   reshape( (/ 6.40943d19, 3.79223d2, 9.33734d-1, & 
            1.45791d23, 1.10225d3, 1.70997d0, &
            4.90775d25, 2.23679d3, 2.21076d0, &
            1.41151d27, 4.66681d3, 2.35028d0, &
            1.28999d29, 6.62329d3, 2.72655d0, &
            2.11583d23, 8.87737d4, 8.82571d-1, &
            1.40334d24, 1.17903d5, 9.98104d-1/), &
                                (/ 3,7 /) ) 

double precision, parameter, dimension(3,8) :: O_param = &
   reshape( (/ 8.43069d19, 3.43915d2, 9.97766d-1, & 
            4.62554d22, 1.42425d3, 1.48807d0, & 
            1.36111d25, 2.78673d3, 1.98390d0, & 
            1.42268d27, 4.66149d3, 2.35156d0, & 
            2.14982d28, 8.31918d3, 2.45386d0, & 
            1.06826d30, 1.11093d4, 2.76372d0, & 
            5.07400d23, 1.37570d5, 8.97953d-1, & 
            2.72798d24, 1.76051d5, 9.97894d-1/), & 
                                (/ 3,8 /) ) 

double precision, parameter, dimension(3,6) :: C_param = &
   reshape( (/ 1.88271d20, 2.58065d2, 1.19846d0, & 
            5.51917d23, 8.22289d2, 1.98802d0, & 
            4.58202d25, 2.26304d3, 2.19830d0, & 
            9.85066d27, 3.53736d3, 2.67447d0, & 
            7.54229d22, 5.30261d4, 8.62810d-1, & 
            6.59436d23, 7.40778d4, 9.99632d-1/), & 
                                (/ 3,6 /) ) 

double precision, parameter, dimension(3,3) :: Xe_param = &
   reshape( (/ 1.38594d20, 2.88544d2, 1.11898d0, & 
                  1.45d24, 6.67163d2, 2.20491d0, &
               1.69d27, 1.24332d3, 2.90652d0 /), &
                              (/ 3,3 /) ) 

public :: neutral

type neutral

  private

  class(fdist2d), pointer :: pf => null()
  class(part2d), pointer :: part => null()

  class(field_rho), allocatable :: q
  class(field_jay), allocatable :: cu, amu
  class(field_djdxi), allocatable :: dcu

  integer :: multi_max
  ! qm is the 8th coordinate for every particle.
  ! If the longitudinal density of the plasma changes,
  ! the qm will change correspondingly
  real :: wp, dt, qm, density
  ! wp is plasma frequency in SI unit
  double precision, dimension(:,:), pointer :: rate_param
  class(field_rho), dimension(:), allocatable :: multi_ion_c
  ! rho_i_discrete is the ion charge density for diagnostic 
  class(field_rho), allocatable :: rho_i_discrete_old_c, rho_i_discrete_c, rho_i_discrete, rho_i_exact_c

  contains

  procedure :: new    => init_neutral
  procedure :: renew  => renew_neutral
  procedure :: update => update_neutral
  procedure :: del    => end_neutral
  procedure :: qdp    => qdeposit_neutral
  procedure :: amjdp  => amjdeposit_neutral
  procedure :: push   => push_neutral
  procedure :: psend  => psend_neutral
  procedure :: precv  => precv_neutral
  procedure :: wr     => writehdf5_neutral
  procedure :: wrq    => writeq_neutral
  procedure :: cbq    => cbq_neutral
  procedure :: get_multi_max

end type neutral

save

character(len=10) :: cls_name = 'neutral'
integer, parameter :: cls_level = 2

contains

subroutine init_neutral( this, opts, pf, max_mode, elem, max_e, qbm, wp, density, &
  s, smth_type, smth_ord )
! element is atomic number. e.g.: For Li, element = 3
! max_e is the maximum number of electrons that the programmer allow the atom
! to lose due to the ionization. It should be less or equal to element
   
  implicit none

  class(neutral), intent(inout) :: this
  type(options), intent(in) :: opts
  class(fdist2d), intent(inout), target :: pf
  integer, intent(in) :: max_mode, elem, max_e
  real, intent(in) :: qbm, wp, density, s
  integer, intent(in), optional :: smth_type, smth_ord

  ! local data
  character(len=18), save :: sname = 'init_neutral'
  integer :: i

  call write_dbg(cls_name, sname, cls_level, 'starts')

  this%pf => pf
  this%wp = wp
  this%dt = opts%get_dxi()

  allocate( this%q, this%cu, this%amu, this%dcu )
  allocate( this%rho_i_discrete_old_c, this%rho_i_discrete_c, this%rho_i_discrete, this%rho_i_exact_c )

  if ( present(smth_type) .and. present(smth_ord) ) then
    call this%q%new( opts, max_mode, p_ps_linear, smth_type, smth_ord )
    call this%cu%new( opts, max_mode, p_ps_linear, smth_type, smth_ord )
    call this%dcu%new( opts, max_mode, p_ps_linear, smth_type, smth_ord )
    call this%amu%new( opts, max_mode, p_ps_linear, smth_type, smth_ord )
  else
    call this%q%new( opts, max_mode, p_ps_linear )
    call this%cu%new( opts, max_mode, p_ps_linear )
    call this%dcu%new( opts, max_mode, p_ps_linear )
    call this%amu%new( opts, max_mode, p_ps_linear )
  endif

  call this%part%new( opts, pf, qbm, this%dt, s, if_empty=.true. )

  call this%rho_i_discrete_old_c%new( opts, max_mode, p_ps_linear )
  call this%rho_i_discrete_c%new( opts, max_mode, p_ps_linear )
  call this%rho_i_discrete%new( opts, max_mode, p_ps_linear )
  call this%rho_i_exact_c%new( opts, max_mode, p_ps_linear )

  this%q  = 0.0
  this%cu = 0.0
  this%rho_i_discrete_old_c = 0.0
  this%rho_i_discrete_c = 0.0
  this%rho_i_discrete = 0.0
  this%rho_i_exact_c = 0.0

  ! call this%q%ag()

  ! TODO: WHY NEED THIS?
  if ( id_stage() == 0 ) then
    call this%q%copy_slice( 1, p_copy_1to2 )      
  endif

  this%multi_max = max_e

  select case ( elem )
    
  case( p_neut_H )

    this%multi_max = min( this%multi_max, size( H_param, 2 ) )   
    allocate( this%rate_param( 3, this%multi_max ) )
    this%rate_param = H_param( :, 1:this%multi_max )

  case( p_neut_He )

    this%multi_max = min( this%multi_max, size( He_param, 2 ) )   
    allocate( this%rate_param( 3, this%multi_max ) )
    this%rate_param = He_param( :, 1:this%multi_max )

  case( p_neut_Li )

    this%multi_max = min( this%multi_max, size( Li_param, 2 ) )   
    allocate( this%rate_param( 3, this%multi_max ) )
    this%rate_param = Li_param( :, 1:this%multi_max )

  case( p_neut_C )

    this%multi_max = min( this%multi_max, size( C_param, 2 ) )   
    allocate( this%rate_param( 3, this%multi_max ) )
    this%rate_param = C_param( :, 1:this%multi_max )

  case( p_neut_N )

    this%multi_max = min( this%multi_max, size( N_param, 2 ) )   
    allocate( this%rate_param( 3, this%multi_max ) )
    this%rate_param = N_param( :, 1:this%multi_max )

  case( p_neut_O )

    this%multi_max = min( this%multi_max, size( O_param, 2 ) )   
    allocate( this%rate_param( 3, this%multi_max ) )
    this%rate_param = O_param( :, 1:this%multi_max )

  case( p_neut_Ar )

    this%multi_max = min( this%multi_max, size( Ar_param, 2 ) )   
    allocate( this%rate_param( 3, this%multi_max ) )
    this%rate_param = Ar_param( :, 1:this%multi_max )

  case( p_neut_K )

    this%multi_max = min( this%multi_max, size( K_param, 2 ) )   
    allocate( this%rate_param( 3, this%multi_max ) )
    this%rate_param = K_param( :, 1:this%multi_max )

  case( p_neut_Rb )

    this%multi_max = min( this%multi_max, size( Rb_param, 2 ) )   
    allocate( this%rate_param( 3, this%multi_max ) )
    this%rate_param = Rb_param( :, 1:this%multi_max )

  case( p_neut_Xe )

    this%multi_max = min( this%multi_max, size( Xe_param, 2 ) )   
    allocate( this%rate_param( 3, this%multi_max ) )
    this%rate_param = Xe_param( :, 1:this%multi_max )

  case( p_neut_Cs )

    this%multi_max = min( this%multi_max, size( Cs_param, 2 ) )   
    allocate( this%rate_param( 3, this%multi_max ) )
    this%rate_param = Cs_param( :, 1:this%multi_max )

  case default
     
    call write_err( 'Invalid neutral gas species!' )

  end select

  allocate( this%multi_ion_c( this%multi_max + 1 ) )

  do i = 1, this%multi_max + 1
    call this%multi_ion_c(i)%new( opts, max_mode, p_ps_linear )
    this%multi_ion_c(i) = 0.0
  enddo

  ! the first one is the neutral gas
  this%multi_ion_c(1) = 1.0

  ! normalized longitudinal density
  this%density = this%pf%get_density(0.0)
  
  ! normalized charge (coordinate of the particle array) 
  this%qm = this%multi_max * this%density * this%pf%den * this%pf%qm / &
    ( abs(this%pf%qm) * real(this%pf%ppc1) * real(this%pf%ppc2) )

  call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine init_neutral

subroutine update_neutral( this, e, psi )

  implicit none
  class(neutral), intent(inout) :: this
  class(field_e), intent(in) :: e
  class(field_psi), intent(in) :: psi
  ! local data
  character(len=18), save :: sname = 'update_neutral'
  integer :: i

  call write_dbg(cls_name, sname, cls_level, 'starts')

  this%rho_i_discrete_old_c = this%rho_i_discrete_c

  call ionize_neutral( this%multi_ion_c, this%rho_i_exact_c, this%rho_i_discrete_c, &
     this%rate_param, e, this%wp, this%dt, (/ this%pf%ppc1, this%pf%ppc2 /), this%multi_max )

  call add_particles( this%part, this%rho_i_discrete_c, this%rho_i_discrete_old_c, psi, &
    (/ this%pf%ppc1, this%pf%ppc2 /), this%multi_max, this%qm )

  call half_offset( this%rho_i_discrete_c, this%rho_i_discrete )

  call this%rho_i_discrete%acopy_gc_f1( dir=p_mpi_forward )

  ! account for the longitudinal varying plasma density
  call dot_f1( this%density, this%rho_i_discrete )

  call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine update_neutral

subroutine ionize_neutral( multi_ion_c, rho_i_exact_c, rho_i_discrete_c, adk_coef, &
  e, wp, dt, ppc, multi_max )

  implicit none

  class(field_rho), intent(inout), dimension(:), pointer :: multi_ion_c
  class(field_rho), intent(inout) :: rho_i_exact_c, rho_i_discrete_c
  double precision, dimension(:,:), intent(in) :: adk_coef
  class(field_e), intent(in) :: e
  real, intent(in) :: wp, dt
  integer, intent(in), dimension(2) :: ppc
  integer, intent(in) :: multi_max

  ! local data
  character(len=18), save :: sname = 'ionize_neutral'
  real, dimension(:,:), pointer :: multi_ion_re => null(), multi_ion_im => null()
  real, dimension(:,:), pointer :: rho_i_ext_re => null(), rho_i_ext_im => null()
  real, dimension(:,:), pointer :: rho_i_disc_re => null(), rho_i_disc_im => null()
  real, dimension(:,:), pointer :: e_re => null(), e_im => null()

  integer :: j, m, nrp, noff, ppc_tot
  real :: eij, e1, e2, e3, cons, inj, dens_temp
  real :: w_ion(20)
  logical :: shoot

  call write_dbg( cls_name, sname, cls_level, 'starts' )
  
  ! the current algorithm can only model the cylindrically symmetrical case.
  e_re          => e%rf_re(0)%get_f1()
  rho_i_ext_re  => rho_i_exact_c%rf_re(0)%get_f1()
  rho_i_disc_re => rho_i_discrete_c%rf_re(0)%get_f1()
  multi_ion_re  => multi_ion_c(1)%rf_re(0)%get_f1()

  nrp     = e%rf_re(0)%get_ndp(1)
  noff    = e%rf_re(0)%get_noff(1)
  ppc_tot = ppc(1) * ppc(2)

  do j = 1, nrp

    ! calculate the cell values
    e1 = 0.5 * ( e_re(1,j) + e_re(1,j+1) )
    e2 = 0.5 * ( e_re(2,j) + e_re(2,j+1) )
    e3 = 0.5 * ( e_re(3,j) + e_re(3,j+1) )

    ! total electric field in unit of GV/m
    eij = sqrt( e1 * e1 + e2 * e2 + e3 * e3 ) * wp * 1.708d-12

    if ( eij > 1.0d-6 ) then

      if ( rho_i_ext_re(1,j) < multi_max ) then

        shoot = .false.

        ! w = r1 * eij^-r3 * EXP(-r2/eij) * 1/wp
        w_ion(1:multi_max) = adk_coef(1,1:multi_max) * eij ** (-adk_coef(3,1:multi_max)) * &
          exp( -adk_coef(2,1:multi_max) / eij ) / wp ! w_ion is in normalized unit now

        ! Use the wrong 2nd order Runge Kutta method in OSIRIS
        ! TODO: check with Yujian
        cons = multi_ion_re(1,j) * w_ion(1) * dt * ( 1.0 + 0.5 * w_ion(1) * dt )

        ! overshoot
        if( cons > multi_ion_re(1,j) ) then
           shoot = .true.
           cons = multi_ion_re(1,j)
        endif
        ! subtract from 0.
        multi_ion_re(1,j) = multi_ion_re(1,j) - cons

        ! For the levels in the middle
        do m = 2, multi_max

          multi_ion_re => multi_ion_c(m)%rf_re(0)%get_f1()

          ! remember how much charge will be added to m-1. level
          inj = cons

          ! overshoot in previous level (we go to time centered densities)
          ! TODO: is dens_temp correct?
          if (shoot) then
            dens_temp = inj * 0.5
            shoot = .false.
          else
            dens_temp = 0.0
          endif

          ! ionizing m-1. level (index m) adding to m. level
          ! 2nd order Runge-Kutta (wrong)
          cons = ( multi_ion_re(1,j) + dens_temp ) * w_ion(m) * dt * ( 1.0 + 0.5 * w_ion(m) * dt )

          ! overshoot in current level (temp might come from previous overshoot)
          if ( cons > multi_ion_re(1,j) + dens_temp ) then
            shoot = .true.
            cons = multi_ion_re(1,j) + dens_temp
          endif

          ! update density for m-1. level (index m)
          ! min only due to round off
          ! Physically, multi_ion_re(1,j) - cons + inj will never exceed 1.0
          ! min is only used to round off the tiny error introduced by a double number (~10^-16)
          multi_ion_re(1,j) = min( multi_ion_re(1,j) - cons + inj, 1.0 )

        enddo ! all levels in the middle

        ! update the last level
        multi_ion_re => multi_ion_c(multi_max+1)%rf_re(0)%get_f1()
        multi_ion_re(1,j) = multi_ion_re(1,j) + cons

        ! End use Runge Kutta method

        ! calculate the total ion charge density by summing over the charge density from all 
        ! types of ions 
        rho_i_ext_re(1,j) = 0.0
        do m = 2, multi_max + 1
          multi_ion_re => multi_ion_c(m)%rf_re(0)%get_f1()
          rho_i_ext_re(1,j) = rho_i_ext_re(1,j) + real(m-1) * multi_ion_re(1,j)
        enddo

        ! We want rho_i_discrete_c that corresponds to release an electron at 'half' step (like OSIRIS and qpic2.0)
        rho_i_disc_re(1,j) = real(multi_max) / real(ppc_tot) * &
          int( rho_i_ext_re(1,j) * real(ppc_tot) / real(multi_max) + 0.5 )

      endif 

    endif ! eij > 1e-6

  enddo

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine ionize_neutral

subroutine add_particles( part, rho_i_disc, rho_i_disc_old, psi, ppc, multi_max, qm )

  implicit none

  class(part2d), intent(inout) :: part
  class(field_rho), intent(in) :: rho_i_disc, rho_i_disc_old
  class(field_psi), intent(in) :: psi
  integer, intent(in), dimension(2) :: ppc
  integer, intent(in) :: multi_max
  real, intent(in) :: qm

  ! local
  character(len=18), save :: sname = 'add_particles'
  real, dimension(:,:), pointer :: rho_i_new => null(), rho_i_old => null(), ppsi => null()
  integer :: nrp, noff, ppc_tot, i, j, k, ppc_add
  real :: r_k, r, dr, theta

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  dr      = psi%get_dr()
  nrp     = psi%rf_re(0)%get_ndp(1)
  noff    = psi%rf_re(0)%get_noff(1)
  ppc_tot = ppc(1) * ppc(2)

  rho_i_new => rho_i_disc%rf_re(0)%get_f1()
  rho_i_old => rho_i_disc_old%rf_re(0)%get_f1()
  ppsi      => psi%rf_re(0)%get_f1()

  do j = 1, nrp

    ! calculate the # of particles to inject based on the difference in ion density between 
    ! the previous time step and the current time step. 
    ppc_add = int( ( rho_i_new(1,j) - rho_i_old(1,j) ) / multi_max * real(ppc_tot) + 0.5 )

    do k = 1, ppc_add

      part%npp = part%npp + 1

      ! r_k is the in-cell position ranging from 0 to 1
      ! r_k = ( real(k) - 0.5 ) / ppc_add
      r_k = 0.5

      i = part%npp
      r = ( r_k + j + noff - 1 ) * dr
      theta = rand() * 2.0 * pi
      part%x(1,i)   = r * cos(theta)
      part%x(2,i)   = r * sin(theta)
      part%p(1,i)   = 0.0
      part%p(2,i)   = 0.0
      part%p(3,i)   = 0.0
      part%gamma(i) = 1.0
      part%q(i)     = qm
      part%psi(i)   = 1.0

    enddo

  enddo

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine add_particles

subroutine half_offset( rho_old, rho_new )

  implicit none

  class(field_rho), intent(in) :: rho_old
  class(field_rho), intent(inout) :: rho_new

  character(len=18), save :: sname = 'half_offset'
  integer :: nrp, j
  real, dimension(:,:), pointer :: f1_old => null(), f1_new => null()

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  nrp    = rho_old%rf_re(0)%get_ndp(1)
  f1_old => rho_old%rf_re(0)%get_f1()
  f1_new => rho_new%rf_re(0)%get_f1()

  f1_new = 0.0

  ! move the in-cell values onto the nodes (by interpolation)
  do j = 1, nrp
    f1_new(1,j)   = f1_new(1,j)   + 0.5 * f1_old(1,j)
    f1_new(1,j+1) = f1_new(1,j+1) + 0.5 * f1_old(1,j)
  enddo

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine half_offset

subroutine renew_neutral( this, s )

  implicit none

  class(neutral), intent(inout) :: this
  real, intent(in) :: s

  ! local data
  character(len=18), save :: sname = 'renew_neutral'
  integer :: i
        
  call write_dbg( cls_name, sname, cls_level, 'starts' )

  call this%part%renew( this%pf, s, if_empty=.true. )

  this%q = 0.0
  this%cu = 0.0
  this%rho_i_discrete_old_c = 0.0
  this%rho_i_discrete_c = 0.0
  this%rho_i_discrete = 0.0
  this%rho_i_exact_c = 0.0

  this%multi_ion_c(1) = 1.0

  do i = 2, this%multi_max + 1
    this%multi_ion_c(i) = 0.0
  enddo

  ! TODO: WHY NEED THIS
  if ( id_stage() == 0 ) then
    call this%q%copy_slice( 1, p_copy_1to2 )       
  endif

  ! normalized longitudinal density
  this%density = this%pf%get_density(s)
  
  ! normalized charge (coordinate of the particle array) 
  this%qm = this%multi_max * this%density * this%pf%den * this%pf%qm &
    / ( abs(this%pf%qm) * real(this%pf%ppc1) * real(this%pf%ppc2) )

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine renew_neutral

subroutine qdeposit_neutral( this, q )

  implicit none
  
  class(neutral), intent(inout) :: this
  class(field_rho), intent(inout) :: q
! local data
  character(len=18), save :: sname = 'qdeposit_neutral'
           
  call write_dbg( cls_name, sname, cls_level, 'starts' )

  this%q = 0.0

  call this%part%qdeposit( this%q )
  call this%q%acopy_gc_f1( dir=p_mpi_forward )
  call this%q%smooth_f1()
  call this%q%copy_gc_f1()
  
  call add_f1( this%q, q )
  call add_f1( this%rho_i_discrete, q )
           
  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine qdeposit_neutral

subroutine amjdeposit_neutral( this, e, b, cu, amu, dcu )

  implicit none
  
  class(neutral), intent(inout) :: this
  class(field_jay), intent(inout) :: cu, amu
  class(field_djdxi), intent(inout) :: dcu
  class(field_e), intent(in) :: e
  class(field_b), intent(in) :: b
  ! local data
  character(len=18), save :: sname = 'amjdeposit_neutral'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  this%cu = 0.0
  this%dcu = 0.0
  this%amu = 0.0
  call this%part%amjdeposit( e, b, this%cu, this%amu, this%dcu )  
  call this%cu%acopy_gc_f1( dir=p_mpi_forward )
  call this%dcu%acopy_gc_f1( dir=p_mpi_forward )
  call this%amu%acopy_gc_f1( dir=p_mpi_forward )
  call this%cu%smooth_f1()
  call this%dcu%smooth_f1()
  call this%amu%smooth_f1()
  call this%cu%copy_gc_f1()
  call this%dcu%copy_gc_f1()
  call this%amu%copy_gc_f1()
  call add_f1( this%cu, cu )
  call add_f1( this%dcu, dcu )
  call add_f1( this%amu, amu )

  call write_dbg( cls_name, sname, cls_level, 'ends' )
  
end subroutine amjdeposit_neutral

subroutine push_neutral( this, e, b )
  
  implicit none
  
  class(neutral), intent(inout) :: this
  class(field_e), intent(in) :: e
  class(field_b), intent(in) :: b
! local data
  character(len=18), save :: sname = 'push_neutral'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  call this%part%push( e, b )
  call this%part%update_bound()
  call move_part2d_comm( this%part )
  
  call write_dbg( cls_name, sname, cls_level, 'ends' )
     
end subroutine push_neutral

subroutine end_neutral( this )

  implicit none
  class(neutral), intent(inout) :: this

end subroutine end_neutral

subroutine psend_neutral( this, tag, id, tag_multi_ion_c, id_multi_ion_c )

  implicit none

  class(neutral), intent(inout) :: this
  integer, intent(in) :: tag
  integer, intent(in), dimension(:) :: tag_multi_ion_c
  integer, intent(inout) :: id
  integer, intent(inout), dimension(:) :: id_multi_ion_c
  ! local data
  character(len=18), save :: sname = 'psend_neutral'
  integer :: i
           
  call write_dbg( cls_name, sname, cls_level, 'starts' )

  call this%part%pipesend( tag, id )

  do i = 1, this%multi_max + 1
    call this%multi_ion_c(i)%pipe_send( tag_multi_ion_c(i), id_multi_ion_c(i), 'forward' )
  enddo 

  call this%rho_i_discrete_c%pipe_send( tag_multi_ion_c(this%multi_max + 2), &
    id_multi_ion_c(this%multi_max + 2), 'forward' )
      
  call this%rho_i_discrete%pipe_send( tag_multi_ion_c(this%multi_max + 3), &
    id_multi_ion_c(this%multi_max + 3), 'forward' )

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine psend_neutral

subroutine precv_neutral( this, tag, tag_multi_ion_c )

  implicit none

  class(neutral), intent(inout) :: this
  integer, intent(in) :: tag
  integer, intent(in), dimension(:) :: tag_multi_ion_c
  ! local data
  character(len=18), save :: sname = 'precv_neutral'
  integer :: i

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  call this%part%piperecv( tag )

  do i = 1, this%multi_max + 1
     call this%multi_ion_c(i)%pipe_recv( tag_multi_ion_c(i), 'forward', 'replace' )
  end do

  call this%rho_i_discrete_c%pipe_recv( tag_multi_ion_c(this%multi_max + 2), &
    'forward', 'replace' )

  call this%rho_i_discrete%pipe_recv( tag_multi_ion_c(this%multi_max + 3), &
    'forward', 'replace' )
           
  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine precv_neutral

subroutine writehdf5_neutral( this, file )

  implicit none

  class(neutral), intent(inout) :: this
  class(hdf5file), intent(in) :: file
  ! local data
  character(len=18), save :: sname = 'writehdf5_neutral'

  call write_dbg(cls_name, sname, cls_level, 'starts')

  call this%part%wr(file)

  call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine writehdf5_neutral

subroutine writeq_neutral( this, files, rtag, stag, id )

  implicit none
  
  class(neutral), intent(inout) :: this
  class(hdf5file), intent(in), dimension(:) :: files
  integer, intent(in) :: rtag, stag
  integer, intent(inout) :: id
! local data
  character(len=18), save :: sname = 'writeq_neutral'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  call this%q%write_hdf5( files, 1, rtag, stag, id ) 

  call write_dbg( cls_name, sname, cls_level, 'ends' )  

end subroutine writeq_neutral
   
subroutine cbq_neutral( this, pos )

  implicit none
  
  class(neutral), intent(inout) :: this
  integer, intent(in) :: pos
  ! local data
  character(len=18), save :: sname = 'cpq_neutral'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  call add_f1( this%cu, this%q, (/3/), (/1/) )
  call this%q%smooth_f1()
  call this%q%copy_slice( pos, p_copy_1to2 )

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine cbq_neutral

function get_multi_max(this)

  implicit none

  class(neutral), intent(in) :: this
  integer :: get_multi_max

  get_multi_max = this%multi_max

end function get_multi_max

end module neutral_class
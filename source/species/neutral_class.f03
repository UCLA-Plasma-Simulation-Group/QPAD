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
  class(part2d_buf), pointer :: part_buf => null()

  class(field_rho), allocatable :: q
  class(field_jay), allocatable :: cu, amu
  class(field_djdxi), allocatable :: dcu

  ! max ionization level
  integer :: multi_max
  ! index of neutral gas residue in multi_ion array
  integer :: idx_neut
  ! index of total ion in multi_ion array
  integer ::idx_ion

  ! qm is the 8th coordinate for every particle.
  ! If the longitudinal density of the plasma changes,
  ! the qm will change correspondingly
  real :: wp, dt, qm, density
  ! wp is plasma frequency in SI unit
  double precision, dimension(:,:), pointer :: rate_param

  ! in-cell ionization level. There are multi_max + 2 levels, the multi_max+1 level 
  ! is for the neutral, the first to multi_max levels are for ions, and the last
  ! level is for the total ion.
  real, dimension(:,:,:), allocatable :: multi_ion
  real, dimension(:,:), allocatable :: ion_old

  ! the ion charge density for diagnostic 
  class(field_rho), allocatable :: rho_ion

  contains

  procedure :: alloc  => alloc_neutral
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
  procedure :: wr_ion => write_ion_neutral
  procedure :: cbq    => cbq_neutral
  procedure :: get_multi_max
  procedure :: ion_deposit => ion_deposit_neutral

end type neutral

! This extended type only store the position for neutral deposition
type, extends( part2d ) :: part2d_buf

   contains

   procedure :: new => init_part2d_buf
   procedure :: end => end_part2d_buf

end type part2d_buf

save

character(len=10) :: cls_name = 'neutral'
integer, parameter :: cls_level = 2

integer, parameter :: p_num_theta = 32

! tables of sine and cosine for azimuthal Fourier decomposition
! real, dimension(:), allocatable, save :: sin_tab, cos_tab

contains

subroutine alloc_neutral( this )

   implicit none

   class(neutral), intent(inout) :: this

   if ( .not. associated( this%part ) ) then
      allocate( part2d :: this%part )
   endif

end subroutine alloc_neutral

subroutine init_neutral( this, opts, pf, max_mode, elem, max_e, qbm, wp, s, &
  smth_type, smth_ord )
! element is atomic number. e.g.: For Li, element = 3
! max_e is the maximum number of electrons that the programmer allow the atom
! to lose due to the ionization. It should be less or equal to element
   
  implicit none

  class(neutral), intent(inout) :: this
  type(options), intent(in) :: opts
  class(fdist2d), intent(inout), target :: pf
  integer, intent(in) :: max_mode, elem, max_e
  real, intent(in) :: qbm, wp, s
  integer, intent(in), optional :: smth_type, smth_ord

  ! local data
  character(len=18), save :: sname = 'init_neutral'
  integer :: i, nrp
  real :: pi2_ntheta

  call write_dbg(cls_name, sname, cls_level, 'starts')

  this%pf => pf
  this%wp = wp
  this%dt = opts%get_dxi()

  allocate( this%q, this%cu, this%amu, this%dcu, this%rho_ion )

  if ( present(smth_type) .and. present(smth_ord) ) then
    call this%q%new( opts, max_mode, p_ps_linear, smth_type, smth_ord )
    call this%rho_ion%new( opts, max_mode, p_ps_linear, smth_type, smth_ord )
    call this%cu%new( opts, max_mode, p_ps_linear, smth_type, smth_ord )
    call this%dcu%new( opts, max_mode, p_ps_linear, smth_type, smth_ord )
    call this%amu%new( opts, max_mode, p_ps_linear, smth_type, smth_ord )
  else
    call this%q%new( opts, max_mode, p_ps_linear )
    call this%rho_ion%new( opts, max_mode, p_ps_linear )
    call this%cu%new( opts, max_mode, p_ps_linear )
    call this%dcu%new( opts, max_mode, p_ps_linear )
    call this%amu%new( opts, max_mode, p_ps_linear )
  endif

  call this%part%new( opts, pf, qbm, this%dt, s, if_empty=.true. )
  call this%part_buf%new( opts, pf, qbm, this%dt, s )

  this%q  = 0.0
  this%cu = 0.0
  this%rho_ion = 0.0

  ! call this%q%ag()

  ! TODO: WHY NEED THIS?
  if ( id_stage() == 0 ) then
    call this%q%copy_slice( 1, p_copy_1to2 )      
  endif

  this%multi_max = max_e
  this%idx_neut  = this%multi_max + 1
  this%idx_ion   = this%multi_max + 2

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

  ! initialize the ionization level array
  nrp = opts%get_ndp(1)
  allocate( this%multi_ion( nrp, p_num_theta, this%multi_max + 2 ) )
  allocate( this%ion_old( nrp, p_num_theta ) )
  this%multi_ion = 0.0
  this%multi_ion( :, :, this%idx_neut ) = 1.0

  ! normalized longitudinal density
  this%density = this%pf%get_density(0.0)
  
  ! normalized charge (coordinate of the particle array) 
  this%qm = this%multi_max * this%density * this%pf%den * this%pf%qm / &
    ( abs(this%pf%qm) * real(this%pf%ppc1) * real(this%pf%ppc2) )

  ! initialize trigonometric function tables
  ! pi2_ntheta = 2.0 * pi / real( p_num_theta )
  ! do i = 1, p_num_theta
  !   sin_tab(i) = sin( pi2_ntheta * (i-1) )
  !   cos_tab(i) = cos( pi2_ntheta * (i-1) )
  ! enddo

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

  this%ion_old = this%multi_ion( :, :, this%idx_ion )

  call ionize_neutral( this%multi_ion, this%rate_param, e, this%wp, this%dt, &
    (/ this%pf%ppc1, this%pf%ppc2 /) )

  ! here multi_ion(:,:,idx_ion) and ion_old should be discrete.
  call add_particles( this%part, this%part_buf, this%multi_ion, this%ion_old, psi, &
    (/ this%pf%ppc1, this%pf%ppc2 /), this%qm )

  call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine update_neutral

subroutine ionize_neutral( multi_ion, adk_coef, e, wp, dt, ppc )

  implicit none

  real, intent(inout), dimension(:,:,:) :: multi_ion
  double precision, dimension(:,:), intent(in) :: adk_coef
  class(field_e), intent(in) :: e
  real, intent(in) :: wp, dt
  integer, intent(in), dimension(2) :: ppc

  ! local data
  character(len=18), save :: sname = 'ionize_neutral'
  real, dimension(:,:), pointer :: e_re => null(), e_im => null()
  integer :: i, j, k, m, nrp, noff, ppc_tot, max_mode, multi_max, idx_neut, idx_ion
  real :: eij, cons, inj, dens_temp, theta, pi2_ntheta
  real :: w_ion(20)
  real, dimension(p_num_theta) :: e1, e2, e3
  complex(kind=DB) :: phase, phase_inc
  logical :: shoot

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  nrp        = e%rf_re(0)%get_ndp(1)
  noff       = e%rf_re(0)%get_noff(1)
  max_mode   = e%get_max_mode()
  ppc_tot    = ppc(1) * ppc(2)
  pi2_ntheta = 2.0 * pi / real( p_num_theta )
  multi_max  = size( multi_ion, 3 ) - 2
  idx_neut   = multi_max + 1
  idx_ion    = multi_max + 2

  do k = 1, p_num_theta

    theta = pi2_ntheta * ( k - 1 )
    phase_inc = cmplx( cos(theta), sin(theta) )

    ! calculate the in-cell values of electric fields
    ! m = 0 mode
    e_re => e%rf_re(0)%get_f1()
    phase = cmplx( 1.0, 0.0 )
    do j = 1, nrp
      e1(j) = 0.5 * ( e_re(1,j) + e_re(1,j+1) )
      e2(j) = 0.5 * ( e_re(2,j) + e_re(2,j+1) )
      e3(j) = 0.5 * ( e_re(3,j) + e_re(3,j+1) )
    enddo

    ! m > 0 mode
    do m = 1, max_mode

      e_re => e%rf_re(m)%get_f1()
      e_im => e%rf_re(m)%get_f1()
      phase = phase * phase_inc

      do j = 1, nrp
        e1(j) = e1(j) + ( e_re(1,j) + e_re(1,j+1) ) * real(phase) - ( e_im(1,j) + e_im(1,j+1) ) * aimag(phase)
        e2(j) = e2(j) + ( e_re(2,j) + e_re(2,j+1) ) * real(phase) - ( e_im(2,j) + e_im(2,j+1) ) * aimag(phase)
        e3(j) = e3(j) + ( e_re(3,j) + e_re(3,j+1) ) * real(phase) - ( e_im(3,j) + e_im(3,j+1) ) * aimag(phase)
      enddo

    enddo

    do j = 1, nrp

      ! total electric field in unit of GV/m
      eij = sqrt( e1(j)**2 + e2(j)**2 + e3(j)**2 ) * wp * 1.708d-12

      if ( eij > 1.0d-6 ) then

        if ( multi_ion( j, k, idx_ion ) < real(multi_max) ) then

          shoot = .false.

          ! w = r1 * eij^-r3 * EXP(-r2/eij) * 1/wp
          ! w_ion is in normalized unit now
          do i = 1, multi_max
            w_ion(i) = adk_coef(1,i) * eij ** (-adk_coef(3,i)) * exp( -adk_coef(2,i) / eij ) / wp
          enddo

          ! Use the wrong 2nd order Runge Kutta method in OSIRIS
          ! TODO: check with Yujian
          cons = multi_ion(j, k, idx_neut) * w_ion(1) * dt * ( 1.0 + 0.5 * w_ion(1) * dt )

          ! overshoot
          if( cons > multi_ion(j, k, idx_neut) ) then
             shoot = .true.
             cons = multi_ion(j, k, idx_neut)
          endif
          ! subtract from 0.
          multi_ion(j, k, idx_neut) = multi_ion(j, k, idx_neut) - cons

          ! For the levels in the middle
          do i = 1, multi_max - 1

            ! remember how much charge will be added to i-th level
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
            cons = ( multi_ion(j,k,i) + dens_temp ) * w_ion(i+1) * dt * ( 1.0 + 0.5 * w_ion(i+1) * dt )

            ! overshoot in current level (temp might come from previous overshoot)
            if ( cons > multi_ion(j,k,i) + dens_temp ) then
              shoot = .true.
              cons = multi_ion(j,k,i) + dens_temp
            endif

            ! update density for i-th level
            ! Physically, multi_ion(j,k,i) - cons + inj will never exceed 1.0
            ! min is only used to round off the tiny error introduced by a double number (~10^-16)
            multi_ion(j,k,i) = min( multi_ion(j,k,i) - cons + inj, 1.0 )

          enddo ! all levels in the middle

          ! update the last level
          multi_ion(j, k, multi_max) = multi_ion(j, k, multi_max) + cons

          ! End use Runge Kutta method

          ! calculate the total ion charge density by summing over the charge density from all 
          ! types of ions 
          multi_ion(j, k, idx_ion) = 0.0
          do i = 1, multi_max
            multi_ion(j, k, idx_ion) = multi_ion(j, k, idx_ion) + real(i) * multi_ion(j,k,i)
          enddo

          ! We want the total ion level that corresponds to release an electron at 'half' step (like OSIRIS and qpic2.0)
          ! Now the multi_ion(:,:,idx_ion) is discrete
          multi_ion(j, k, idx_ion) = real(multi_max) / real(ppc_tot) * &
            int( multi_ion(j, k, idx_ion) * real(ppc_tot) / real(multi_max) + 0.5 )

        endif 

      endif ! eij > 1e-6

    enddo
  enddo

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine ionize_neutral

subroutine add_particles( part, part_buf, multi_ion, ion_old, psi, ppc, qm )

  implicit none

  class(part2d), intent(inout) :: part
  class(part2d_buf), intent(inout) :: part_buf
  real, intent(in), dimension(:,:,:) :: multi_ion
  real, intent(in), dimension(:,:) :: ion_old
  integer, intent(in), dimension(2) :: ppc
  class(field_psi), intent(in) :: psi
  real, intent(in) :: qm

  ! local
  character(len=18), save :: sname = 'add_particles'
  integer :: nrp, noff, ppc_tot, i, j, k, pp1, pp2, ppc_add, multi_max, idx_ion
  real :: r, dr, theta, dtheta

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  dr      = psi%get_dr()
  dtheta  = 2.0 * pi / real( p_num_theta )
  nrp     = psi%rf_re(0)%get_ndp(1)
  noff    = psi%rf_re(0)%get_noff(1)
  ppc_tot = ppc(1) * ppc(2)
  multi_max = size( multi_ion, 3 ) - 2
  idx_ion = multi_max + 2
  part_buf%npp = 0

  do k = 1, p_num_theta
    do j = 1, nrp

      ! calculate the # of particles to inject based on the difference in ion density between 
      ! the previous time step and the current time step.
      ! here multi_ion(:,:,idx_ion) and ion_old should both be discrete.
      ppc_add = int( ( multi_ion(j, k, idx_ion) - ion_old(j, k) ) / multi_max * real(ppc_tot) + 0.5 )

      pp1 = part%npp
      pp2 = part_buf%npp
      do i = 1, ppc_add

        pp1 = pp1 + 1
        pp2 = pp2 + 1

        ! randomly generate injection position
        r     = ( rand() + j + noff - 1 ) * dr
        theta = ( rand() + k - 1.5 ) * dtheta

        part%x(1,pp1)   = r * cos(theta)
        part%x(2,pp1)   = r * sin(theta)
        part%p(1,pp1)   = 0.0
        part%p(2,pp1)   = 0.0
        part%p(3,pp1)   = 0.0
        part%gamma(pp1) = 1.0
        part%q(pp1)     = qm
        part%psi(pp1)   = 1.0

        ! store the position of ionized particles for ion charge deposition
        part_buf%x(1,pp2) = part%x(1,pp1)
        part_buf%x(2,pp2) = part%x(2,pp1)
        part_buf%q(pp2)   = -part%q(pp1) ! note the sign

      enddo
      part%npp = part%npp + ppc_add
      part_buf%npp = part_buf%npp + ppc_add

    enddo
  enddo

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine add_particles

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
  this%rho_ion = 0.0

  this%multi_ion = 0.0
  this%multi_ion( :, :, this%idx_neut ) = 1.0

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

subroutine qdeposit_neutral( this, q_tot )

  implicit none
  
  class(neutral), intent(inout) :: this
  class(field_rho), intent(inout) :: q_tot
! local data
  character(len=18), save :: sname = 'qdeposit_neutral'
           
  call write_dbg( cls_name, sname, cls_level, 'starts' )

  this%q = 0.0

  call this%part%qdeposit( this%q )
  call this%q%acopy_gc_f1( dir=p_mpi_forward )
  call this%q%smooth_f1()
  call this%q%copy_gc_f1()
  
  call add_f1( this%q, q_tot )
           
  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine qdeposit_neutral

subroutine ion_deposit_neutral( this, q_tot )

  implicit none
  
  class(neutral), intent(inout) :: this
  class(field_rho), intent(inout) :: q_tot
  ! local data
  character(len=18), save :: sname = 'ion_deposit_neutral'
           
  call write_dbg( cls_name, sname, cls_level, 'starts' )

  this%rho_ion = 0.0

  call this%part_buf%qdeposit( this%rho_ion )
  call this%rho_ion%acopy_gc_f1( dir=p_mpi_forward )
  call this%rho_ion%smooth_f1()
  call this%rho_ion%copy_gc_f1()
  
  call add_f1( this%rho_ion, q_tot )

  ! clear up the buffer
  this%part_buf%npp = 0
           
  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine ion_deposit_neutral

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

subroutine psend_neutral( this, tag_q, id_q, tag_ion, id_ion )

  implicit none

  class(neutral), intent(inout) :: this
  integer, intent(in) :: tag_q, tag_ion
  integer, intent(inout) :: id_q, id_ion
  ! local data
  character(len=18), save :: sname = 'psend_neutral'
  integer :: i, idproc_des, ierr, count
           
  call write_dbg( cls_name, sname, cls_level, 'starts' )
  call start_tprof( 'pipeline' )

  call this%part%pipesend( tag_q, id_q )

  if ( id_stage() == num_stages() - 1 ) then
    id_ion = MPI_REQUEST_NULL
    call stop_tprof( 'pipeline' )
    call write_dbg( cls_name, sname, cls_level, 'ends' )
    return
  endif

  idproc_des = id_proc() + num_procs_loc()
  count = size( this%multi_ion )

  call mpi_isend( this%multi_ion, count, p_dtype_real, idproc_des, tag_ion, &
    comm_world(), id_ion, ierr )
  ! check for error
  if ( ierr /= 0 ) then
    call write_err( 'MPI_ISEND failed.' )
  endif

  call stop_tprof( 'pipeline' )
  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine psend_neutral

subroutine precv_neutral( this, tag_q, tag_ion )

  implicit none

  class(neutral), intent(inout) :: this
  integer, intent(in) :: tag_q, tag_ion
  ! local data
  character(len=18), save :: sname = 'precv_neutral'
  integer, dimension(MPI_STATUS_SIZE) :: stat
  integer :: i, idproc_src, ierr, count

  call write_dbg( cls_name, sname, cls_level, 'starts' )
  call start_tprof( 'pipeline' )

  call this%part%piperecv( tag_q )

  if ( id_stage() == 0 ) then
    call stop_tprof( 'pipeline' )
    call write_dbg( cls_name, sname, cls_level, 'ends' )
    return
  endif

  idproc_src = id_proc() - num_procs_loc()
  count = size( this%multi_ion)

  call mpi_recv( this%multi_ion, count, p_dtype_real, idproc_src, tag_ion, &
    comm_world(), stat, ierr )
  ! check for error
  if ( ierr /= 0 ) then
    call write_err( 'MPI_RECV failed.' )
  endif
           
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

subroutine write_ion_neutral( this, files, rtag, stag, id )

  implicit none
  
  class(neutral), intent(inout) :: this
  class(hdf5file), intent(in), dimension(:) :: files
  integer, intent(in) :: rtag, stag
  integer, intent(inout) :: id
  ! local data
  character(len=18), save :: sname = 'write_ion_neutral'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  call this%rho_ion%write_hdf5( files, 1, rtag, stag, id ) 

  call write_dbg( cls_name, sname, cls_level, 'ends' )  

end subroutine write_ion_neutral
   
subroutine cbq_neutral( this, pos )

  implicit none
  
  class(neutral), intent(inout) :: this
  integer, intent(in) :: pos
  ! local data
  character(len=18), save :: sname = 'cbq_neutral'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  call add_f1( this%cu, this%q, (/3/), (/1/) )
  call this%q%smooth_f1()
  call this%q%copy_slice( pos, p_copy_1to2 )

  call this%rho_ion%smooth_f1()
  call this%rho_ion%copy_slice( pos, p_copy_1to2 )

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine cbq_neutral

function get_multi_max(this)

  implicit none

  class(neutral), intent(in) :: this
  integer :: get_multi_max

  get_multi_max = this%multi_max

end function get_multi_max

subroutine init_part2d_buf( this, opts, pf, qbm, dt, s, if_empty )

  implicit none

  class(part2d_buf), intent(inout) :: this
  type(options), intent(in) :: opts
  class(fdist2d), intent(inout) :: pf
  real, intent(in) :: qbm, dt, s
  logical, intent(in), optional :: if_empty

  ! local data
  character(len=18), save :: sname = 'init_part2d_buf'
  integer :: npmax, nbmax

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  this%qbm = qbm
  this%dt  = dt
  this%part_dim = 3

  ! this is max number of particles to be ionized
  npmax      = pf%ppc1 * pf%ppc2 * p_num_theta
  nbmax      = max(int(0.1*npmax),100)
  this%npmax = npmax
  this%nbmax = nbmax
  this%npp   = 0

  this%dr   = pf%getdex()
  this%edge = opts%get_nd(1) * this%dr

  ! only allocate position and charge, no need to initialize with profile
  allocate( this%x( 2, npmax ), this%q( npmax ) )

  call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine init_part2d_buf

subroutine end_part2d_buf(this)

   implicit none

   class(part2d_buf), intent(inout) :: this
   character(len=18), save :: sname = 'end_part2d_buf'

   call write_dbg(cls_name, sname, cls_level, 'starts')

   deallocate( this%x, this%q )

   call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine end_part2d_buf

end module neutral_class
module diagnostics_class

use parallel_module
use options_class
use hdf5io_class
use sim_fields_class
use sim_beams_class
use sim_plasma_class
use field_class
use beam3d_class
use species2d_class
use neutral_class
use neutral2_class
use sysutil_module
use param
use input_class
use mpi

implicit none

private

public :: sim_diag, diag_node

type diag_node
  
  class( diag_node ), pointer :: next => null()
  ! HDF5 files
  type( hdf5file ), dimension(:), allocatable :: files
  ! Number of HDF5 files
  integer :: num_files
  ! Target object to be associated
  class(*), pointer :: obj => null()
  ! Text diagnostic file identifier
  integer :: fid
  ! Dimension of the grid data
  integer :: dim = 0
  ! Dump frequency
  integer :: dump_freq = 0
  ! Sampling percentage
  integer :: psample = 0
  ! Diagnostic type
  character(len=32) :: type_label = ''
  ! Tag and identifier for pipeline data communication
  integer :: tag = -1, id = MPI_REQUEST_NULL


  contains

  procedure :: new => init_diag_node
  procedure :: del => end_diag_node
  procedure :: set_sim_time

end type diag_node

type sim_diag

  class( diag_node ), pointer :: head => null()
  class( diag_node ), pointer :: diag => null()
  integer :: num_diag = 0
  integer :: ndump_gcd = 0 ! greatest common divisor of all the values of ndump
  logical :: has_vpotz = .false.
  logical :: has_vpott = .false.
  integer :: n ! the number of ion
  integer :: m ! 2*max_mode + 1

  contains

  procedure :: alloc => alloc_sim_diag
  procedure :: new => init_sim_diag
  procedure :: del => end_sim_diag
  procedure :: run => run_sim_diag
  procedure :: set_ndump_gcd
  procedure :: to_next, to_head, to_tail, is_tail
  generic :: add_diag => add_diag_cym, add_diag_raw, add_diag_rst, add_diag_ion

  procedure, private :: add_diag_cym, add_diag_raw, add_diag_rst, add_diag_ion
  procedure, private :: init_diag_beams
  procedure, private :: init_diag_plasma
  procedure, private :: init_diag_fields
  procedure, private :: init_diag_rst

end type sim_diag

contains

subroutine init_diag_node( this, obj, dump_freq, num_files, dim, psample, id )

  implicit none

  class( diag_node ), intent(inout) :: this
  class(*), intent(in), target :: obj
  integer, intent(in) :: dump_freq
  integer, intent(in), optional :: num_files, dim, psample, id

  integer, save :: cls_level = 3
  character(len=32), save :: cls_name = 'diag_node'
  character(len=32), save :: sname = 'init_diag_node'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  this%obj => obj
  this%dump_freq = dump_freq
  if ( present(num_files) ) this%num_files = num_files
  if ( present(dim) ) this%dim = dim
  if ( present(psample) ) this%psample = psample
  if ( present(id) ) this%id = id

  if ( this%num_files > 0 ) then
    allocate( this%files( this%num_files ) )
  else
    call write_err( 'Number of diagnostic files must be greater than 0' )
  endif

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine init_diag_node

subroutine end_diag_node( this )

  implicit none

  class( diag_node ), intent(inout) :: this

  integer, save :: cls_level = 3
  character(len=32), save :: cls_name = 'diag_node'
  character(len=32), save :: sname = 'end_diag_node'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  deallocate( this%files )

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine end_diag_node

subroutine set_sim_time( this, tstep, time )

  implicit none

  class( diag_node ), intent(inout) :: this
  integer, intent(in) :: tstep
  real, intent(in) :: time

  integer :: i
  integer, save :: cls_level = 3
  character(len=32), save :: cls_name = 'diag_node'
  character(len=32), save :: sname = 'set_sim_time'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  do i = 1, this%num_files
    call this%files(i)%new( n = tstep, t = time )
  enddo

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine set_sim_time

subroutine alloc_sim_diag( this, input )

  implicit none

  class( sim_diag ), intent(inout) :: this
  type( input_json ), intent(inout) :: input

  ! placeholder, do nothing currently

end subroutine alloc_sim_diag

subroutine init_diag_beams( this, input, beams )

  implicit none

  class( sim_diag ), intent(inout) :: this
  type( input_json ), intent(inout) :: input
  class( sim_beams ), intent(in), target :: beams

  ! local data
  integer :: nbeams, max_mode, ndump, psample
  integer :: i, j, k, m, n
  real :: rmin, rmax, zmin, zmax, dt
  character(len=:), allocatable :: ss

  call input%get( 'simulation.max_mode', max_mode )
  call input%get( 'simulation.nbeams', nbeams )

  call input%get( 'simulation.box.r(1)', rmin )
  call input%get( 'simulation.box.r(2)', rmax )
  call input%get( 'simulation.box.z(1)', zmin )
  call input%get( 'simulation.box.z(2)', zmax )
  call input%get( 'simulation.dt', dt )

  ! add beam diagnostics
  do i = 1, nbeams
    call input%info( 'beam('//num2str(i)//').diag', n_children=m )
    do j = 1, m
      call input%get( 'beam('//num2str(i)//').diag'//'('//num2str(j)//').ndump', ndump )
      if ( ndump > 0 ) then
        call input%info( 'beam('//num2str(i)//').diag'//'('//num2str(j)//').name', n_children=n )
        do k = 1, n
          if ( allocated(ss) ) deallocate(ss)
          call input%get( 'beam('//num2str(i)//').diag'//'('//num2str(j)//').name'//'('//num2str(k)//')', ss )
          select case ( trim(ss) )
          case ( 'charge_cyl_m' )
            call this%add_diag( &
              obj       = beams%beam(i), &
              max_mode  = max_mode, &
              dump_freq = ndump, &
              dim       = 1, &
              type_label= 'charge_cyl_m', &
              filename  = './Beam'//num2str(i)//'/Charge/', &
              dataname  = 'charge', &
              timeunit  = '1 / \omega_p', &
              dt        = dt, &
              axisname  = (/'r  ', '\xi', '   '/), &
              axislabel = (/'r  ', '\xi', '   '/), &
              axisunits = (/'c / \omega_p', 'c / \omega_p', '            '/), &
              axismax   = (/rmax, zmax, 0.0/), &
              axismin   = (/rmin, zmin, 0.0/), &
              units     = 'n_0', &
              label     = '\rho', &
              rank      = 2 )
          ! case ( 'charge' )
          ! this is the diagnostics for complete combined field, to be implemented
          case ( 'raw' )
            call input%get( 'beam('//num2str(i)//').diag'//'('//num2str(j)//').psample', psample )
            call this%add_diag( &
              obj       = beams%beam(i), &
              dump_freq = ndump, &
              psample   = psample, &
              type_label= 'raw', &
              filename  = './Beam'//num2str(i)//'/Raw/', &
              dataname  = 'raw', &
              timeunit  = '1 / \omega_p', &
              dt        = dt, &
              units     = '', &
              label     = 'Beam Raw' )
          end select
        enddo ! end of k
      endif
    enddo ! end of j
  enddo ! end of i

end subroutine init_diag_beams

subroutine init_diag_plasma( this, input, plasma )

  implicit none

  class( sim_diag ), intent(inout) :: this
  type( input_json ), intent(inout) :: input
  class( sim_plasma ), intent(in), target :: plasma
  ! local data
  integer :: nspecies, nneutrals, nneutral2s, max_mode, ndump, psample
  integer :: i, j, k, m, n, v, l,imax
  real :: rmin, rmax, zmin, zmax, dt
  character(len=:), allocatable :: ss

  call input%get( 'simulation.max_mode', max_mode )
  call input%get( 'simulation.nspecies', nspecies )
  call input%get( 'simulation.nneutrals', nneutrals )
  call input%get( 'simulation.nneutral2s', nneutral2s )
  call input%get( 'simulation.box.r(1)', rmin )
  call input%get( 'simulation.box.r(2)', rmax )
  call input%get( 'simulation.box.z(1)', zmin )
  call input%get( 'simulation.box.z(2)', zmax )
  call input%get( 'simulation.dt', dt )

  ! add species diagnostics
  do i = 1, nspecies
    call input%info( 'species('//num2str(i)//').diag', n_children=m )
    do j = 1, m
      call input%get( 'species('//num2str(i)//').diag'//'('//num2str(j)//').ndump', ndump )
      if ( ndump > 0 ) then
        call input%info( 'species('//num2str(i)//').diag'//'('//num2str(j)//').name', n_children=n )
        do k = 1, n
          if ( allocated(ss) ) deallocate(ss)
          call input%get( 'species('//num2str(i)//').diag'//'('//num2str(j)//').name'//'('//num2str(k)//')', ss )
          select case ( trim(ss) )
          case ( 'charge_cyl_m' )
            call this%add_diag( &
              obj       = plasma%spe(i), &
              max_mode  = max_mode, &
              dump_freq = ndump, &
              dim       = 1, &
              type_label= 'charge_cyl_m', &
              filename  = './Species'//num2str(i)//'/'//'Charge'//'/', &
              dataname  = 'charge', &
              timeunit  = '1 / \omega_p', &
              dt        = dt, &
              axisname  = (/'r  ', '\xi', '   '/), &
              axislabel = (/'r  ', '\xi', '   '/), &
              axisunits = (/'c / \omega_p', 'c / \omega_p', '            '/), &
              axismax   = (/rmax, zmax, 0.0/), &
              axismin   = (/rmin, zmin, 0.0/), &
              units     = 'n_0', &
              label     = '\rho', &
              rank      = 2 )
          case ( 'raw' )
            call input%get( 'species('//num2str(i)//').diag'//'('//num2str(j)//').psample', psample )
            call this%add_diag( &
              obj       = plasma%spe(i), &
              dump_freq = ndump, &
              psample   = psample, &
              type_label= 'raw', &
              filename  = './Species'//num2str(i)//'/Raw/', &
              dataname  = 'raw', &
              timeunit  = '1 / \omega_p', &
              dt        = dt, &
              units     = '', &
              label     = 'Species Raw' )
          end select
        enddo ! end of k
      endif
    enddo ! end of j
  enddo ! end of i

  ! add neutral diagnostics
  do i = 1, nneutrals
    call input%info( 'neutrals('//num2str(i)//').diag', n_children=m )
    do j = 1, m
      call input%get( 'neutrals('//num2str(i)//').diag'//'('//num2str(j)//').ndump', ndump )
      if ( ndump > 0 ) then
        call input%info( 'neutrals('//num2str(i)//').diag'//'('//num2str(j)//').name', n_children=n )
        do k = 1, n
          if ( allocated(ss) ) deallocate(ss)
          call input%get( 'neutrals('//num2str(i)//').diag'//'('//num2str(j)//').name'//'('//num2str(k)//')', ss )
          select case ( trim(ss) )
          case ( 'charge_cyl_m' )
            call this%add_diag( &
              obj       = plasma%neut(i), &
              max_mode  = max_mode, &
              dump_freq = ndump, &
              dim       = 1, &
              type_label= 'charge_cyl_m', &
              filename  = './Neutral'//num2str(i)//'/'//'Charge'//'/', &
              dataname  = 'charge', &
              timeunit  = '1 / \omega_p', &
              dt        = dt, &
              axisname  = (/'r  ', '\xi', '   '/), &
              axislabel = (/'r  ', '\xi', '   '/), &
              axisunits = (/'c / \omega_p', 'c / \omega_p', '            '/), &
              axismax   = (/rmax, zmax, 0.0/), &
              axismin   = (/rmin, zmin, 0.0/), &
              units     = 'n_0', &
              label     = '\rho', &
              rank      = 2 )
          case ( 'ion_cyl_m' )
            call this%add_diag( &
              obj       = plasma%neut(i), &
              max_mode  = max_mode, &
              dump_freq = ndump, &
              dim       = 1, &
              type_label= 'ion_cyl_m', &
              filename  = './Neutral'//num2str(i)//'/'//'Ion_charge'//'/', &
              dataname  = 'ion_charge', &
              timeunit  = '1 / \omega_p', &
              dt        = dt, &
              axisname  = (/'r  ', '\xi', '   '/), &
              axislabel = (/'r  ', '\xi', '   '/), &
              axisunits = (/'c / \omega_p', 'c / \omega_p', '            '/), &
              axismax   = (/rmax, zmax, 0.0/), &
              axismin   = (/rmin, zmin, 0.0/), &
              units     = 'n_0', &
              label     = '\rho', &
              rank      = 2 )
          case ( 'raw' )
            call input%get( 'neutrals('//num2str(i)//').diag'//'('//num2str(j)//').psample', psample )
            call this%add_diag( &
              obj       = plasma%neut(i), &
              dump_freq = ndump, &
              psample   = psample, &
              type_label= 'raw', &
              filename  = './Neutral'//num2str(i)//'/Raw/', &
              dataname  = 'raw', &
              timeunit  = '1 / \omega_p', &
              dt        = dt, &
              units     = '', &
              label     = 'Neutral Raw' )
          end select
        enddo ! end of k
      endif
    enddo ! end of j
  enddo ! end of i

  ! add neutral2 diagnostics
  do i = 1, nneutral2s
    call input%info( 'neutral2s('//num2str(i)//').diag', n_children=m )
    call input%get( 'neutral2s('//num2str(i)//').v', v )
    call input%get('neutral2s('//num2str(i)//').ion_max',imax)
    if(v == 0) then
      imax = imax + 1
    else
      imax = imax - v + 1
    endif
    write(2,*) v,imax, "neutral2 diagnostics_v"
    this%n = imax
    do j = 1, m
      call input%get( 'neutral2s('//num2str(i)//').diag'//'('//num2str(j)//').ndump', ndump )
      if ( ndump > 0 ) then
        call input%info( 'neutral2s('//num2str(i)//').diag'//'('//num2str(j)//').name', n_children=n )
        do k = 1, n
          if ( allocated(ss) ) deallocate(ss)
          call input%get( 'neutral2s('//num2str(i)//').diag'//'('//num2str(j)//').name'//'('//num2str(k)//')', ss )
          select case ( trim(ss) )
          case ( 'charge_cyl_m' )
            call this%add_diag( &
              obj       = plasma%neut2(i), &
              max_mode  = max_mode, &
              dump_freq = ndump, &
              dim       = 1, &
              type_label= 'charge_cyl_m', &
              filename  = './Neutral2'//num2str(i)//'/'//'Charge'//'/', &
              dataname  = 'charge', &
              timeunit  = '1 / \omega_p', &
              dt        = dt, &
              axisname  = (/'r  ', '\xi', '   '/), &
              axislabel = (/'r  ', '\xi', '   '/), &
              axisunits = (/'c / \omega_p', 'c / \omega_p', '            '/), &
              axismax   = (/rmax, zmax, 0.0/), &
              axismin   = (/rmin, zmin, 0.0/), &
              units     = 'n_0', &
              label     = '\rho', &
              rank      = 2 )
          case ( 'ion_cyl_m' )
!             do l = 1,imax
!               case ( 'ion_cyl_m'//num2str(l) )
              call this%add_diag( &
                obj       = plasma%neut2(i), &
                max_mode  = max_mode, &
                dump_freq = ndump, &
                dim       = 1, &
                type_label= 'ion_cyl_m', &
                filename  = './Neutral2'//num2str(i)//'/'//'Ion_charge'//'/', &
                dataname  = 'ion_charge', &
                timeunit  = '1 / \omega_p', &
                dt        = dt, &
                imax      = imax,&
                v         = v,&
                axisname  = (/'r  ', '\xi', '   '/), &
                axislabel = (/'r  ', '\xi', '   '/), &
                axisunits = (/'c / \omega_p', 'c / \omega_p', '            '/), &
                axismax   = (/rmax, zmax, 0.0/), &
                axismin   = (/rmin, zmin, 0.0/), &
                units     = 'n_0', &
                label     = '\rho', &
                rank      = 2 )
!             enddo
          case ( 'raw' )
            call input%get( 'neutral2s('//num2str(i)//').diag'//'('//num2str(j)//').psample', psample )
            call this%add_diag( &
              obj       = plasma%neut2(i), &
              dump_freq = ndump, &
              psample   = psample, &
              type_label= 'raw', &
              filename  = './Neutral2'//num2str(i)//'/Raw/', &
              dataname  = 'raw', &
              timeunit  = '1 / \omega_p', &
              dt        = dt, &
              units     = '', &
              label     = 'Neutral2 Raw' )
          end select
        enddo ! end of k
      endif
    enddo ! end of j
  enddo ! end of i

end subroutine init_diag_plasma

subroutine init_diag_fields( this, input, opts, fields )

  implicit none

  class( sim_diag ), intent(inout) :: this
  type( input_json ), intent(inout) :: input
  type( options ), intent(in) :: opts
  class( sim_fields ), intent(inout), target :: fields
  ! local data
  integer :: max_mode, ndump, dim
  integer :: j, k, m, n
  real :: rmin, rmax, zmin, zmax, dt
  character(len=32) :: sn1, sn2, sn3, sn4
  character(len=:), allocatable :: ss
  class(*), pointer :: obj => null()

  call input%get( 'simulation.max_mode', max_mode )
  call input%get( 'simulation.box.r(1)', rmin )
  call input%get( 'simulation.box.r(2)', rmax )
  call input%get( 'simulation.box.z(1)', zmin )
  call input%get( 'simulation.box.z(2)', zmax )
  call input%get( 'simulation.dt', dt )

  call input%info( 'field.diag', n_children=m )
  do j = 1, m
    call input%get( 'field.diag'//'('//num2str(j)//').ndump', ndump )
    if ( ndump > 0 ) then
      call input%info( 'field.diag'//'('//num2str(j)//').name', n_children=n )
      do k = 1, n
        if ( allocated(ss) ) deallocate(ss)
        call input%get( 'field.diag'//'('//num2str(j)//').name'//'('//num2str(k)//')', ss )
        select case ( trim(ss) )
        case ( 'er_cyl_m' )
          sn1 = 'Er'
          sn2 = 'er'
          sn3 = 'm_ec\omega_p/e'
          sn4 = 'E_r'
          dim = 1
          obj => fields%e
        case ( 'ephi_cyl_m' )
          sn1 = 'Ephi'
          sn2 = 'ephi'
          sn3 = 'm_ec\omega_p/e'
          sn4 = 'E_\phi'
          dim = 2
          obj => fields%e
        case ( 'ez_cyl_m' )
          sn1 = 'Ez'
          sn2 = 'ez'
          sn3 = 'm_ec\omega_p/e'
          sn4 = 'E_z'
          dim = 3
          obj => fields%e
        case ( 'br_cyl_m' )
          sn1 = 'Br'
          sn2 = 'br'
          sn3 = 'm_e\omega_p/e'
          sn4 = 'B_r'
          dim = 1
          obj => fields%b
        case ( 'bphi_cyl_m' )
          sn1 = 'Bphi'
          sn2 = 'bphi'
          sn3 = 'm_e\omega_p/e'
          sn4 = 'B_\phi'
          dim = 2
          obj => fields%b
        case ( 'bz_cyl_m' )
          sn1 = 'Bz'
          sn2 = 'bz'
          sn3 = 'm_e\omega_p/e'
          sn4 = 'B_z'
          dim = 3
          obj => fields%b
        case ( 'spec_er_cyl_m' )
          sn1 = 'Er_spec'
          sn2 = 'er'
          sn3 = 'm_ec\omega_p/e'
          sn4 = 'E_r'
          dim = 1
          obj => fields%e_spe
        case ( 'spec_ephi_cyl_m' )
          sn1 = 'Ephi_spec'
          sn2 = 'ephi'
          sn3 = 'm_ec\omega_p/e'
          sn4 = 'E_\phi'
          dim = 2
          obj => fields%e_spe
        case ( 'spec_br_cyl_m' )
          sn1 = 'Br_spec'
          sn2 = 'br'
          sn3 = 'm_e\omega_p/e'
          sn4 = 'B_r'
          dim = 1
          obj => fields%b_spe
        case ( 'spec_bphi_cyl_m' )
          sn1 = 'Bphi_spec'
          sn2 = 'bphi'
          sn3 = 'm_e\omega_p/e'
          sn4 = 'B_\phi'
          dim = 2
          obj => fields%b_spe
        case ( 'spec_bz_cyl_m' )
          sn1 = 'Bz_spec'
          sn2 = 'bz'
          sn3 = 'm_e\omega_p/e'
          sn4 = 'B_z'
          dim = 3
          obj => fields%b_spe
        case ( 'charge_cyl_m' )
          sn1 = 'Charge'
          sn2 = 'charge'
          sn3 = 'n_0'
          sn4 = '\rho'
          dim = 1
          obj => fields%q_spe
        case ( 'jr_cyl_m' )
          sn1 = 'Jr'
          sn2 = 'jr'
          sn3 = 'n_0 c'
          sn4 = 'J_r'
          dim = 1
          obj => fields%cu
        case ( 'jphi_cyl_m' )
          sn1 = 'Jphi'
          sn2 = 'jphi'
          sn3 = 'n_0 c'
          sn4 = 'J_\phi'
          dim = 2
          obj => fields%cu
        case ( 'jz_cyl_m' )
          sn1 = 'Jz'
          sn2 = 'jz'
          sn3 = 'n_0 c'
          sn4 = 'J_z'
          dim = 3
          obj => fields%cu
        case ( 'psi_cyl_m' )
          sn1 = 'Psi'
          sn2 = 'psi'
          sn3 = 'm_ec^2/e'
          sn4 = '\Psi'
          dim = 1
          obj => fields%psi
        case ( 'az_cyl_m' )
          sn1 = 'Az'
          sn2 = 'az'
          sn3 = 'm_ec/e'
          sn4 = 'A_z'
          dim = 3
          if ( .not. associated(fields%vpot) ) then
            allocate( fields%vpot )
            call fields%vpot%alloc( max_mode )
            call fields%vpot%new( opts, max_mode, p_ps_linear, &
              p_bnd_open )
          endif
          obj => fields%vpot
          this%has_vpotz = .true.
        case ( 'ar_cyl_m' )
          sn1 = 'Ar'
          sn2 = 'ar'
          sn3 = 'm_ec/e'
          sn4 = 'A_r'
          dim = 1
          if ( .not. associated(fields%vpot) ) then
            allocate( fields%vpot )
            call fields%vpot%alloc( max_mode )
            call fields%vpot%new( opts, max_mode, p_ps_linear, &
              p_bnd_open )
          endif
          obj => fields%vpot
          this%has_vpott = .true.
        case ( 'aphi_cyl_m' )
          sn1 = 'Aphi'
          sn2 = 'aphi'
          sn3 = 'm_ec/e'
          sn4 = 'A_\phi'
          dim = 2
          if ( .not. associated(fields%vpot) ) then
            allocate( fields%vpot )
            call fields%vpot%alloc( max_mode )
            call fields%vpot%new( opts, max_mode, p_ps_linear, &
              p_bnd_open )
          endif
          obj => fields%vpot
          this%has_vpott = .true.
        end select
        call this%add_diag( &
          obj       = obj, &
          max_mode  = max_mode, &
          dump_freq = ndump, &
          dim       = dim, &
          type_label= trim(ss), &
          filename  = './Fields/'//trim(sn1)//'/', &
          dataname  = trim(sn2), &
          timeunit  = '1 / \omega_p', &
          dt        = dt, &
          axisname  = (/'r  ', '\xi', '   '/), &
          axislabel = (/'r  ', '\xi', '   '/), &
          axisunits = (/'c / \omega_p', 'c / \omega_p', '            '/), &
          axismax   = (/rmax, zmax, 0.0/), &
          axismin   = (/rmin, zmin, 0.0/), &
          units     = trim(sn3), &
          label     = trim(sn4), &
          rank      = 2 )
      enddo ! end of k
    endif
  enddo ! end of j

end subroutine init_diag_fields

subroutine init_diag_rst( this, input, beams )

  implicit none

  class( sim_diag ), intent(inout) :: this
  type( input_json ), intent(inout) :: input
  class( sim_beams ), intent(in), target :: beams

  integer :: i, ndump, nbeams

  call input%get( 'simulation.ndump_restart', ndump )
  call input%get( 'simulation.nbeams', nbeams )
  do i = 1, nbeams
    call this%add_diag( &
      obj        = beams%beam(i), &
      dump_freq  = ndump, &
      type_label = 'rst', &
      filename   = './RST/Beam'//num2str(i,2)//'/', &
      dataname   = 'RST-beam'//num2str(i,2)//'-'//num2str(id_proc(),6) )
  enddo

end subroutine init_diag_rst

subroutine init_sim_diag( this, input, opts, fields, beams, plasma )

  implicit none

  class( sim_diag ), intent(inout) :: this
  type( input_json ), intent(inout) :: input
  type( options ), intent(in) :: opts
  class( sim_fields ), intent(inout) :: fields
  class( sim_beams ), intent(in) :: beams
  class( sim_plasma ), intent(in) :: plasma

  ! local data
  integer :: ierr
  logical :: rst

  integer, save :: cls_level = 2
  character(len=32), save :: cls_name = 'sim_diag'
  character(len=32), save :: sname = 'init_sim_diag'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  this%num_diag = 0

  call input%get( 'simulation.dump_restart', rst )

  ! initialize beam diagnostics
  call this%init_diag_beams( input, beams )

  ! initialize species diagnostics
  call this%init_diag_plasma( input, plasma )

  ! initialize field diagnostics
  call this%init_diag_fields( input, opts, fields )

  ! initialize restart file diagnostics
  if (rst) then
    call this%init_diag_rst( input, beams )
  endif

  call this%set_ndump_gcd()

  call mpi_barrier( comm_world(), ierr )

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine init_sim_diag

subroutine end_sim_diag( this )

  implicit none

  class( sim_diag ), intent(inout) :: this
  class( diag_node ), pointer :: to_be_del

  integer, save :: cls_level = 2
  character(len=32), save :: cls_name = 'sim_diag'
  character(len=32), save :: sname = 'end_sim_diag'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  call this%to_head()
  do while ( .not. this%is_tail() )
    to_be_del => this%diag
    call this%to_next()
    deallocate( to_be_del )
  enddo

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine end_sim_diag

subroutine run_sim_diag( this, tstep, dt )

  implicit none

  class( sim_diag ), intent(inout) :: this
  integer, intent(in) :: tstep
  real, intent(in) :: dt

  ! local data
  integer :: stag, rtag, ierr
  integer, dimension(MPI_STATUS_SIZE) :: istat
  integer :: i, n

  integer, save :: cls_level = 3
  character(len=32), save :: cls_name = 'sim_diag'
  character(len=32), save :: sname = 'run_sim_diag'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  if ( mod( tstep, this%ndump_gcd ) /= 0 ) return

  call this%to_head()
  if ( .not. associated(this%diag) ) return

  do
    if ( mod( tstep, this%diag%dump_freq ) == 0 ) then

      call this%diag%set_sim_time( tstep, tstep*dt )

      select type ( obj => this%diag%obj )

      class is ( field )

        rtag = ntag(); stag = rtag
        call mpi_wait( this%diag%id, istat, ierr )
        call obj%write_hdf5( this%diag%files, this%diag%dim, rtag, stag, this%diag%id )
      
      class is ( beam3d )

        select case ( trim( this%diag%type_label ) )

        case ( 'raw' )

          rtag = ntag(); stag = rtag
          call mpi_wait( this%diag%id, istat, ierr )
          call obj%wr( this%diag%files(1), this%diag%psample, rtag, stag, this%diag%id )

        case ( 'charge_cyl_m' )

          rtag = ntag(); stag = rtag
          call mpi_wait( this%diag%id, istat, ierr )
          call obj%wrq( this%diag%files, rtag, stag, this%diag%id )

        case ( 'rst' )

          call obj%wrst( this%diag%files(1) )

        end select

      class is ( species2d )

        select case ( trim( this%diag%type_label ) )

        case ( 'raw' )
          call obj%wr( this%diag%files(1) )

        case ( 'charge_cyl_m' )

          rtag = ntag(); stag = rtag
          call obj%wrq( this%diag%files, rtag, stag, this%diag%id )

        end select

      class is ( neutral )

        select case ( trim( this%diag%type_label ) )

        case ( 'raw' )
          call obj%wr( this%diag%files(1) )

        case ( 'charge_cyl_m' )

          rtag = ntag(); stag = rtag
          call obj%wrq( this%diag%files, rtag, stag, this%diag%id )

        case ( 'ion_cyl_m' )

          rtag = ntag(); stag = rtag
            call obj%wr_ion( this%diag%files, rtag, stag, this%diag%id )

        end select
        
      class is ( neutral2 )

        select case ( trim( this%diag%type_label ) )

        case ( 'raw' )

          rtag = ntag(); stag = rtag
          call obj%wr( this%diag%files(1) )

        case ( 'charge_cyl_m' )

          rtag = ntag(); stag = rtag
          call obj%wrq( this%diag%files, rtag, stag, this%diag%id )
          write(2,*) "true", "charge_cyl_m"

        case ( 'ion_cyl_m' )
          do i = 1,this%n
              rtag = ntag(); stag = rtag;
              call obj%wr_ion( this%diag%files((i-1)*this%m+1:i*this%m), rtag, stag, this%diag%id, i )
               write(2,*) i, "ion_cyl_m"
          enddo
        end select

      end select

    endif

    if ( .not. this%is_tail() ) then
      call this%to_next()
    else
      exit
    endif

  enddo

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine run_sim_diag

! add 2d-plot diagnostics of cylindrical modes
subroutine add_diag_cym( this, obj, max_mode, dump_freq, dim, type_label, filename, &
  dataname, timeunit, dt, axisname, axislabel, axisunits, axismax, axismin, units, label, rank )

  implicit none

  class( sim_diag ), intent(inout) :: this
  class(*), intent(in) :: obj
  integer, intent(in) :: max_mode, dump_freq, dim
  character(len=*), intent(in) :: filename, timeunit, dataname, units, label, type_label
  integer, intent(in) :: rank
  real, intent(in) :: dt
  character(len=*), dimension(3), intent(in) :: axisname, axislabel, axisunits
  real, dimension(3), intent(in) :: axismax, axismin

  character(len=64) :: cym_str
  integer :: i

  integer, save :: cls_level = 2
  character(len=32), save :: cls_name = 'sim_diag'
  character(len=32), save :: sname = 'add_diag_cym'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  if ( .not. associated( this%head ) ) then
    allocate( this%head )
    call this%head%new( obj, dump_freq, num_files=2*max_mode+1, dim=dim, id=MPI_REQUEST_NULL )
    this%diag => this%head
  else
    call this%to_tail()
    allocate( this%diag%next )
    this%diag => this%diag%next
    call this%diag%new( obj, dump_freq, num_files=2*max_mode+1, dim=dim, id=MPI_REQUEST_NULL )
  endif

  this%diag%type_label = trim(type_label)

  do i = 1, 2*max_mode+1

    if ( i == 1 ) then
      cym_str = 'Re0'
    elseif ( mod(i,2)==0 ) then
      cym_str = 'Re'//num2str(i/2)
    else
      cym_str = 'Im'//num2str(i/2)
    endif
    call system( 'mkdir -p '//trim(filename)//trim(cym_str)//'/' )
    call this%diag%files(i)%new( &
      timeunit  = timeunit, &
      dt        = dt, &
      axisname  = axisname, &
      axislabel = axislabel, &
      axismin   = axismin, &
      axismax   = axismax, &
      rank      = rank, &
      filename  = trim(filename)//trim(cym_str)//'/', &
      dataname  = dataname, &
      units     = units, &
      label     = label )

  enddo
  this%num_diag = this%num_diag + 1

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine add_diag_cym

! add ion
subroutine add_diag_ion( this, obj, max_mode, dump_freq, dim, type_label, filename, &
  dataname, timeunit, dt, imax, v, axisname, axislabel, axisunits, axismax, axismin, units, label, rank )

  implicit none

  class( sim_diag ), intent(inout) :: this
  class(*), intent(in) :: obj
  integer, intent(in) :: max_mode, dump_freq, dim, imax, v
  character(len=*), intent(in) :: filename, timeunit, dataname, units, label, type_label
  integer, intent(in) :: rank
  real, intent(in) :: dt
  character(len=*), dimension(3), intent(in) :: axisname, axislabel, axisunits
  real, dimension(3), intent(in) :: axismax, axismin

  character(len=64) :: ion_str, cym_str
  integer :: i, j

  integer, save :: cls_level = 2
  character(len=32), save :: cls_name = 'sim_diag'
  character(len=32), save :: sname = 'add_diag_ion'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  if ( .not. associated( this%head ) ) then
    allocate( this%head )
    call this%head%new( obj, dump_freq, num_files=imax*(2*max_mode+1), dim=dim, id=MPI_REQUEST_NULL )
    this%diag => this%head
  else
    call this%to_tail()
    allocate( this%diag%next )
    this%diag => this%diag%next
    call this%diag%new( obj, dump_freq, num_files=imax*(2*max_mode+1), dim=dim, id=MPI_REQUEST_NULL )
  endif

  write(2,*) imax, "imax"
  this%diag%type_label = trim(type_label)
  this%n = imax
  this%m = 2*max_mode + 1
  do j = 1, imax

    ion_str = 'ion'//num2str(v+j-1)
!     call system( 'mkdir -p '//trim(filename)//trim(ion_str)//'/' )
    do i = 1, 2*max_mode+1

      if ( i == 1 ) then
        cym_str = 'Re0'
      elseif ( mod(i,2)==0 ) then
        cym_str = 'Re'//num2str(i/2)
      else
        cym_str = 'Im'//num2str(i/2)
      endif
      call system( 'mkdir -p '//trim(filename)//trim(ion_str)//'/'//trim(cym_str)//'/' )
      call this%diag%files((2*max_mode+1)*(j-1)+i)%new( &
        timeunit  = timeunit, &
        dt        = dt, &
        axisname  = axisname, &
        axislabel = axislabel, &
        axismin   = axismin, &
        axismax   = axismax, &
        rank      = rank, &
        filename  = trim(filename)//trim(ion_str)//'/'//trim(cym_str)//'/', &
        dataname  = dataname, &
        units     = units, &
        label     = label )

    enddo
  this%num_diag = this%num_diag + 1
  enddo
  write(2,*) this%n, "add_diag_ion"
  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine add_diag_ion


! add restart file for particles
subroutine add_diag_rst( this, obj, dump_freq, type_label, filename, dataname )

  implicit none

  class( sim_diag ), intent(inout) :: this
  class(*), intent(in) :: obj
  integer, intent(in) :: dump_freq
  character(len=*), intent(in) :: filename, dataname, type_label

  integer, save :: cls_level = 2
  character(len=32), save :: cls_name = 'sim_diag'
  character(len=32), save :: sname = 'add_diag_rst'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  if ( .not. associated( this%head ) ) then
    allocate( this%head )
    call this%head%new( obj, dump_freq, num_files=1 )
    this%diag => this%head
  else
    call this%to_tail()
    allocate( this%diag%next )
    this%diag => this%diag%next
    call this%diag%new( obj, dump_freq, num_files=1 )
  endif

  this%diag%type_label = trim(type_label)

  call system( 'mkdir -p '//trim(filename) )
  call this%diag%files(1)%new( &
    ty        = 'restart', &
    filename  = filename, &
    dataname  = dataname )

  this%num_diag = this%num_diag + 1

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine add_diag_rst

! add raw data for particles
subroutine add_diag_raw( this, obj, dump_freq, psample, type_label, filename, &
  dataname, timeunit, dt, units, label )

  implicit none

  class( sim_diag ), intent(inout) :: this
  class(*), intent(in) :: obj
  integer, intent(in) :: dump_freq, psample
  character(len=*), intent(in) :: filename, timeunit, dataname, units, label, type_label
  real, intent(in) :: dt

  integer, save :: cls_level = 2
  character(len=32), save :: cls_name = 'sim_diag'
  character(len=32), save :: sname = 'add_diag_raw'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  if ( .not. associated( this%head ) ) then
    allocate( this%head )
    call this%head%new( obj, dump_freq, num_files=1, psample=psample, id=MPI_REQUEST_NULL )
    this%diag => this%head
  else
    call this%to_tail()
    allocate( this%diag%next )
    this%diag => this%diag%next
    call this%diag%new( obj, dump_freq, num_files=1, psample=psample, id=MPI_REQUEST_NULL )
  endif

  this%diag%type_label = trim(type_label)

  call system( 'mkdir -p '//trim(filename) )
  call this%diag%files(1)%new( &
    timeunit  = timeunit, &
    dt        = dt, &
    ty        = 'particles', &
    filename  = filename, &
    dataname  = dataname, &
    units     = units, &
    label     = label )

  this%num_diag = this%num_diag + 1

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine add_diag_raw

subroutine set_ndump_gcd( this )

  implicit none

  class( sim_diag ), intent(inout) :: this

  ! local data
  integer :: a, b, temp

  integer, save :: cls_level = 2
  character(len=32), save :: cls_name = 'sim_diag'
  character(len=32), save :: sname = 'set_ndump_gcd'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  call this%to_head()
  this%ndump_gcd = this%diag%dump_freq
  do
    ! use Euclidean algorithm to calculate GCD
    a = this%ndump_gcd
    b = this%diag%dump_freq
    if ( b == 0 ) b = 1
    do while ( mod(a,b) /= 0 )
      temp = a
      a = b
      b = mod(temp,b)
    enddo
    this%ndump_gcd = min( this%ndump_gcd, b )
    if ( this%is_tail() ) then
      exit
    else
      call this%to_next()
    endif
  enddo

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine set_ndump_gcd

! move to the next node of link table
subroutine to_next( this )

  implicit none
  class( sim_diag ), intent(inout) :: this

  if ( associated( this%head ) ) then
    if ( .not. this%is_tail() ) then
      this%diag => this%diag%next
    endif
  endif

end subroutine to_next

! move to the first node of link table
subroutine to_head( this )

  implicit none
  class( sim_diag ), intent(inout) :: this

  if ( associated( this%head ) ) then
    this%diag => this%head
  endif

end subroutine to_head

! move to the last node of link table
subroutine to_tail( this )

  implicit none
  class( sim_diag ), intent(inout) :: this

  if ( associated( this%head ) ) then
    do while ( .not. this%is_tail() )
      this%diag => this%diag%next
    enddo
  endif

end subroutine to_tail

function is_tail( this )

  implicit none
  class( sim_diag ), intent(in) :: this
  logical :: is_tail

  is_tail = ( .not. associated( this%diag%next ) )

end function is_tail

end module diagnostics_class
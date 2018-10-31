module diagnostics_class

use parallel_pipe_class
use hdf5io_class
use sim_fields_class
use sim_beams_class
use sim_species_class
use field_class
use system
use input_class

implicit none

private

public :: sim_diag

type diag_node

  ! private

  class( diag_node ), pointer :: next => null()
  type( hdf5file ), dimension(:), allocatable :: files
  class(*), pointer :: obj => null()
  integer :: num_files=1, dim = 0, df=0, psample=0

  contains

  generic :: new => init_diag_node
  generic :: del => end_diag_node
  procedure :: set_sim_time

  procedure, private :: init_diag_node, end_diag_node

end type diag_node

type sim_diag

  private

  class( parallel_pipe ), pointer :: pp => null()
  class( diag_node ), pointer :: head => null()
  class( diag_node ), pointer :: diag => null()
  integer :: num_diag = 0
  integer :: ndump_gcd = 0 ! greatest common divisor of all the values of ndump

  contains

  generic :: new => init_sim_diag
  generic :: del => end_sim_diag
  generic :: run => run_sim_diag
  generic :: add_diag => add_diag_cym, add_diag_raw, add_diag_rst

  procedure, private :: init_sim_diag, end_sim_diag, run_sim_diag
  procedure, private :: add_diag_cym, add_diag_raw, add_diag_rst
  procedure, private :: to_next, to_head, to_tail, is_tail
  procedure, private :: set_ndump_gcd

end type sim_diag

contains

subroutine init_diag_node( this, obj, df, num_files, dim, psample )

  implicit none

  class( diag_node ), intent(inout) :: this
  class(*), intent(in), target :: obj
  integer, intent(in) :: df
  integer, intent(in), optional :: num_files, dim, psample

  integer, save :: cls_level = 3
  character(len=32), save :: cls_name = 'diag_node'
  character(len=32), save :: sname = 'init_diag_node'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  this%obj => obj
  this%df = df
  if ( present(num_files) ) this%num_files = num_files
  if ( present(dim) ) this%dim = dim
  if ( present(psample) ) this%psample = psample

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

subroutine init_sim_diag( this, pp, input, fields, beams, species )

  implicit none

  class( sim_diag ), intent(inout) :: this
  class( parallel_pipe ), intent(in), pointer :: pp
  class( input_json ), intent(inout) :: input
  class( sim_fields ), intent(in), target :: fields
  class( sim_beams ), intent(in), target :: beams
  class( sim_species ), intent(in), target :: species

  ! local data
  integer :: nbeams, nspecies, num_modes, ndump, psample, dim
  integer :: i, j, k, m, n, ierr
  real :: rmin, rmax, zmin, zmax, dt
  logical :: rst
  character(len=32) :: sn1, sn2, sn3, sn4
  character(len=:), allocatable :: ss
  class(*), pointer :: obj => null()

  integer, save :: cls_level = 2
  character(len=32), save :: cls_name = 'sim_diag'
  character(len=32), save :: sname = 'init_sim_diag'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  this%pp => pp
  this%num_diag = 0

  call input%get( 'simulation.num_modes', num_modes )
  call input%get( 'simulation.nbeams', nbeams )
  call input%get( 'simulation.nspecies', nspecies )

  call input%get( 'simulation.box.r(1)', rmin )
  call input%get( 'simulation.box.r(2)', rmax )
  call input%get( 'simulation.box.z(1)', zmin )
  call input%get( 'simulation.box.z(2)', zmax )
  call input%get( 'simulation.dt', dt )

  call input%get( 'simulation.dump_restart', rst )

! ===============================================================================
! THIS PART IS TO BE FINISHED
! ===============================================================================
  ! ! add beam diagnostics
  ! do i = 1, nbeams
  !   call input%info( 'beam('//num2str(i),').diag', n_children=m )
  !   do j = 1, m
  !     call input%get( 'beam('//num2str(i),').diag'//'('//num2str(j)//').ndump', ndump )
  !     if ( ndump > 0 ) then
  !       call input%info( 'beam('//num2str(i),').diag'//'('//num2str(j)//').name', n_children=n )
  !       do k = 1, n
  !         if ( allocated(ss) ) deallocate(ss)
  !         call input%get( 'beam('//num2str(i),').diag'//'('//num2str(j)//').name'//'('//num2str(k)//')', ss )
  !         select case ( trim(ss) )
  !         case ( 'charge_cyl_m' )
  !           call this%add_diag( &
  !             obj       = beams%beam(i), &
  !             num_modes = num_modes, &
  !             df        = ndump, &
  !             dim       = 1, &
  !             filename  = './Beam'//num2str(i)//'/Charge/', &
  !             dataname  = 'charge', &
  !             timeunit  = '1 / \omega_p', &
  !             dt        = dt, &
  !             axisname  = (/'r', '\xi', ''/), &
  !             axislabel = (/'r', '\xi', ''/), &
  !             axisunits = (/'c / \omega_p', 'c / \omega_p', ''/), &
  !             axismax   = (/rmax, zmax, 0.0/), &
  !             axismin   = (/rmin, zmin, 0.0/), &
  !             units     = 'n_0', &
  !             label     = 'Charge Density', &
  !             rank      = 2 )
  !         ! case ( 'charge' )
  !         ! this is the diagnostics for complete combined field, to be implemented
  !         case ( 'raw' )
  !           call input%get( 'beam('//num2str(i),').diag'//'('//num2str(j)//').psample', psample )
  !           call this%add_diag( &
  !             obj       = beams%beam(i), &
  !             df        = ndump, &
  !             psample   = psample, &
  !             filename  = './Beam'//num2str(i)//'/Raw/', &
  !             dataname  = 'raw', &
  !             timeunit  = '1 / \omega_p', &
  !             dt        = dt, &
  !             ty        = 'particles', &
  !             units     = '', &
  !             label     = 'Beam Raw' )
  !         end select
  !       enddo ! end of k
  !     endif
  !   enddo ! end of j
  ! enddo ! end of i

  ! ! add species diagnostics
  ! do i = 1, nspecies
  !   call input%info( 'species('//num2str(i),').diag', n_children=m )
  !   do j = 1, m
  !     call input%get( 'species('//num2str(i),').diag'//'('//num2str(j)//').ndump', ndump )
  !     if ( ndump > 0 ) then
  !       call input%info( 'species('//num2str(i),').diag'//'('//num2str(j)//').name', n_children=n )
  !       do k = 1, n
  !         if ( allocated(ss) ) deallocate(ss)
  !         call input%get( 'species('//num2str(i),').diag'//'('//num2str(j)//').name'//'('//num2str(k)//')', ss )
  !         select case ( trim(ss) )
  !         case ( 'charge_cyl_m' )
  !           sn1 = 'Charge'
  !           sn2 = 'charge'
  !           sn3 = 'n_0'
  !           sn4 = 'Charge Density'
  !           dim = 1
  !           obj => species%spe(i)
  !         case ( 'jr_cyl_m' )
  !           sn1 = 'Jr'
  !           sn2 = 'jr'
  !           sn3 = 'n_0 c'
  !           sn4 = 'J_r'
  !           dim = 1
  !           obj => fields%jay
  !         case ( 'jphi_cyl_m' )
  !           sn1 = 'Jphi'
  !           sn2 = 'jphi'
  !           sn3 = 'n_0 c'
  !           sn4 = 'J_\phi'
  !           dim = 2
  !           obj => field%jay
  !         case ( 'jz_cyl_m' )
  !           sn1 = 'Jz'
  !           sn2 = 'jz'
  !           sn3 = 'n_0 c'
  !           sn4 = 'J_z'
  !           dim = 2
  !           obj => field%jay
  !         end select
  !         call this%add_diag( &
  !           obj       = obj, &
  !           num_modes = num_modes, &
  !           df        = ndump, &
  !           dim       = dim, &
  !           filename  = './Species'//num2str(i)//trim(sn1)//'/', &
  !           dataname  = trim(sn2), &
  !           timeunit  = '1 / \omega_p', &
  !           dt        = dt, &
  !           axisname  = (/'r', '\xi', ''/), &
  !           axislabel = (/'r', '\xi', ''/), &
  !           axisunits = (/'c / \omega_p', 'c / \omega_p', ''/), &
  !           axismax   = (/rmax, zmax, 0.0/), &
  !           axismin   = (/rmin, zmin, 0.0/), &
  !           units     = trim(sn3), &
  !           label     = trim(sn4), &
  !           rank      = 2 )
  !       enddo ! end of k
  !     endif
  !   enddo ! end of j
  ! enddo ! end of i

! ===============================================================================

  ! add field diagnostics
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
          sn3 = 'mc\omega_p/e'
          sn4 = 'Electric Field'
          dim = 1
          obj => fields%e_spe
        case ( 'ephi_cyl_m' )
          sn1 = 'Ephi'
          sn2 = 'ephi'
          sn3 = 'mc\omega_p/e'
          sn4 = 'Electric Field'
          dim = 2
          obj => fields%e_spe
        case ( 'ez_cyl_m' )
          sn1 = 'Ez'
          sn2 = 'ez'
          sn3 = 'mc\omega_p/e'
          sn4 = 'Electric Field'
          dim = 3
          obj => fields%e_spe
        case ( 'br_cyl_m' )
          sn1 = 'Br'
          sn2 = 'br'
          sn3 = 'mc\omega_p/e'
          sn4 = 'Magnetic Field'
          dim = 1
          obj => fields%b_spe
        case ( 'bphi_cyl_m' )
          sn1 = 'Bphi'
          sn2 = 'bphi'
          sn3 = 'mc\omega_p/e'
          sn4 = 'Magnetic Field'
          dim = 2
          obj => fields%b_spe
        case ( 'bz_cyl_m' )
          sn1 = 'Bz'
          sn2 = 'bz'
          sn3 = 'mc\omega_p/e'
          sn4 = 'Magnetic Field'
          dim = 3
          obj => fields%b_spe
        case ( 'psi_cyl_m' )
          sn1 = 'Psi'
          sn2 = 'psi'
          sn3 = 'mc^2'
          sn4 = '\Psi'
          dim = 1
          obj => fields%psi
        end select
        call this%add_diag( &
          obj       = obj, &
          num_modes = num_modes, &
          df        = ndump, &
          dim       = dim, &
          filename  = './Fields'//num2str(i)//trim(sn1)//'/', &
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

! ===============================================================================
! THIS PART IS TO BE FINISHED
! ===============================================================================
  ! if (rst) then
  !   call input%get( 'simulation.ndump_restart', ndump ) 
  !   do i = 1, nbeams
  !      call this%add_diag( &
  !       obj      = beams%beam(i), &
  !       df       = ndump, &
  !       filename = './RST/Beam-'//num2str(i)//'/', &
  !       dataname = 'RST-beam'//num2str(i)//'-'//num2str(pp%getidproc(),10) )
  !   enddo
  ! endif
! ===============================================================================

  call this%set_ndump_gcd()

  call MPI_BARRIER( pp%getlworld(), ierr )

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
  class(*), pointer :: obj => null()

  integer, save :: cls_level = 3
  character(len=32), save :: cls_name = 'sim_diag'
  character(len=32), save :: sname = 'run_sim_diag'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  if ( .not. ( mod( tstep-1, this%ndump_gcd ) == 0 ) ) return

  call this%to_head()
  do while ( .not. this%is_tail() )
    call this%diag%set_sim_time( tstep, tstep*dt )
    select type ( obj => this%diag%obj )
    class is ( field )
      call obj%write_hdf5( this%diag%files, this%diag%dim )
    ! class is ( beam )
    ! ! to be implemented
    ! class is ( species )
    ! ! to be implemented
    end select
    call this%to_next()
  enddo

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine run_sim_diag

! add 2d-plot diagnostics of cylindrical modes
subroutine add_diag_cym( this, obj, num_modes, df, dim, filename, dataname, timeunit, dt, &
  axisname, axislabel, axisunits, axismax, axismin, units, label, rank )

  implicit none

  class( sim_diag ), intent(inout) :: this
  class(*), intent(in) :: obj
  integer, intent(in) :: num_modes, df, dim
  character(len=*), intent(in) :: filename, timeunit, dataname, units, label
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
    call this%head%new( obj, df, num_files=2*num_modes+1, dim=dim )
    this%diag => this%head
  else
    call this%to_tail()
    allocate( this%diag%next )
    call this%diag%next%new( obj, df, num_files=2*num_modes+1, dim=dim )
    call this%to_next()
  endif

  do i = 1, 2*num_modes+1

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

! add restart file for particles
subroutine add_diag_rst( this, obj, df, filename, dataname, ty )

  implicit none

  class( sim_diag ), intent(inout) :: this
  class(*), intent(in) :: obj
  integer, intent(in) :: df
  character(len=*), intent(in) :: filename, dataname, ty

  integer, save :: cls_level = 2
  character(len=32), save :: cls_name = 'sim_diag'
  character(len=32), save :: sname = 'add_diag_rst'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  if ( .not. associated( this%head ) ) then
    allocate( this%head )
    call this%head%new( obj, df, num_files=1 )
    this%diag => this%head
  else
    call this%to_tail()
    allocate( this%diag%next )
    call this%diag%next%new( obj, df, num_files=1 )
    call this%to_next()
  endif

  call this%diag%files(1)%new( &
    ty        = ty, &
    filename  = filename, &
    dataname  = dataname )

  this%num_diag = this%num_diag + 1

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine add_diag_rst

! add raw data for particles
subroutine add_diag_raw( this, obj, df, psample, filename, dataname, timeunit, dt, &
  ty, units, label )

  implicit none

  class( sim_diag ), intent(inout) :: this
  class(*), intent(in) :: obj
  integer, intent(in) :: df, psample
  character(len=*), intent(in) :: filename, timeunit, dataname, units, label, ty
  real, intent(in) :: dt

  integer, save :: cls_level = 2
  character(len=32), save :: cls_name = 'sim_diag'
  character(len=32), save :: sname = 'add_diag_raw'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  if ( .not. associated( this%head ) ) then
    allocate( this%head )
    call this%head%new( obj, df, num_files=1, psample=psample )
    this%diag => this%head
  else
    call this%to_tail()
    allocate( this%diag%next )
    call this%diag%next%new( obj, df, num_files=1, psample=psample )
    call this%to_next()
  endif

  call this%diag%files(1)%new( &
    timeunit  = timeunit, &
    dt        = dt, &
    ty        = ty, &
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
  integer :: a, b

  integer, save :: cls_level = 2
  character(len=32), save :: cls_name = 'sim_diag'
  character(len=32), save :: sname = 'set_ndump_gcd'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  call this%to_head()
  this%ndump_gcd = this%diag%df
  do while ( .not. this%is_tail() )
    ! use Euclidean algorithm to calculate GCD
    a = this%diag%df
    b = this%diag%next%df
    do while ( .not. ( mod(a,b) == 0 ) )
      a = b
      b = mod(a,b)
    enddo
    this%ndump_gcd = min( this%ndump_gcd, b )
    call this%to_next()
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
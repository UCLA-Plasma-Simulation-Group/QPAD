module simulation_test_class

use simulation_class
use field_src_class
use field_psi_class
use ufield_class
use field_class
use system
use param
use debug_tool

implicit none
private

public :: simulation_test

type, extends( simulation ) :: simulation_test

  private

  contains

  generic :: run_test => run_simulation_test
  procedure, private :: run_simulation_test

end type simulation_test

character(len=18), save :: cls_name = 'simulation_test'
integer, save :: cls_level = 1

contains

subroutine run_simulation_test( this )

  implicit none

  class( simulation_test ), intent(inout) :: this

  integer :: i
  character(len=32), save :: sname = 'run_simulation_test'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  call set_rho( this%fields%q_spe )

  do i = 1, this%nstep1d

    call this%fields%q_spe%copy_slice( i, dir=p_copy_2to1 )
    call this%fields%psi%solve( this%fields%q_spe )
    call this%fields%psi%copy_slice( i, dir=p_copy_1to2 )

  enddo

  call this%diag%run( tstep=1, dt=this%dt )
  ! call write_field2d( this%fields%psi, 'psi', 1 )
  ! call write_field2d( this%fields%q_spe, 'q', 1 )

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine run_simulation_test

subroutine set_rho( q )

  implicit none

  class( field_rho ), intent(inout) :: q

  class( ufield ), dimension(:), pointer :: uq_re, uq_im
  integer :: num_modes, nrp, nr, nz, noff, i, j, mode
  real :: dr, dxi, r, xi

  uq_re     => q%get_rf_re()
  uq_im     => q%get_rf_im()
  num_modes = q%get_num_modes()
  dr        = q%get_dr()
  dxi       = q%get_dxi()
  nrp       = uq_re(0)%get_ndp(1)
  nr        = uq_re(0)%get_nd(1)
  nz        = uq_re(0)%get_nd(2)
  noff      = uq_re(0)%get_noff(1)

  do mode = 0, num_modes

    do j = 1, nz
      xi = (j-1)*dxi
      do i = 1, nrp
        r = ( real(i+noff)-0.5 )*dr
        uq_re(mode)%f2(1,i,j) = 0.5 * exp( -0.5*(r/0.2)**2 - 0.5*((xi-0.5)/0.2)**2 )
      enddo
    enddo

    if ( mode == 0 ) cycle

    do j = 1, nz
      xi = (j-1)*dxi
      do i = 1, nrp
        r = ( real(i+noff)-0.5 )*dr
        uq_im(mode)%f2(1,i,j) = 0.5 * exp( -0.5*(r/0.2)**2 - 0.5*((xi-0.5)/0.2)**2 )
      enddo
    enddo

  enddo

end subroutine set_rho

subroutine write_field2d( f, name, dim )

  implicit none

  class( field ), intent(in) :: f
  character(len=*), intent(in) :: name
  integer, intent(in) :: dim
  
  class( ufield ), dimension(:), pointer :: uf_re => null(), uf_im => null()
  real, dimension(:,:,:), pointer :: p => null()

  integer :: num_modes, i, mode

  uf_re => f%get_rf_re()
  uf_im => f%get_rf_im()

  num_modes = f%get_num_modes()

  do mode = 0, num_modes
    p => uf_re(mode)%get_f2()
    call write_data( p, trim(name)//'-re-'//num2str(mode)//'.txt', dim )
    if ( mode == 0 ) cycle
    p => uf_im(mode)%get_f2()
    call write_data( p, trim(name)//'-im-'//num2str(mode)//'.txt', dim )
  enddo

end subroutine write_field2d

end module simulation_test_class


program test_simulation

use simulation_test_class

implicit none

type( simulation_test ) :: sim

call sim%new()
call sim%run_test()
call sim%del()

end program test_simulation
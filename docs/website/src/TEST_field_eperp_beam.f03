program test_field_eperp

use field_b_class
use field_e_class
use field_class
use field_src_class
use system
use param
use mpi
use ufield_class
use debug_tool

implicit none

type( field_b ) :: b
type( field_e ) :: e
type( field_rho ) :: rho

integer :: num_modes = 3, dim = 3, order = p_fs_2order, part_shape = p_ps_linear
integer, dimension(2) :: nd = (/128, 1/), nvp = (/1, 1/)
integer, dimension(2,2) :: gc_num
real :: dr, dxi, r

type( ufield ), dimension(:), pointer :: uq_re => null(), uq_im => null()
type( ufield ), dimension(:), pointer :: ue_re => null(), ue_im => null()
real, dimension(:,:), pointer :: p
integer :: ierr, i, mode

call MPI_INIT( ierr )
call init_errors( 2, 3 )

call write_dbg( 'main', 'test_field_eperp', 0, 'starts' )

dr = 1.0 / (nd(1)-0.5)
dxi = 1.0

gc_num(:,1) = (/0,0/)
gc_num(:,2) = (/0,0/)

! Test beam field
call rho%new( num_modes, dr, dxi, nd, nvp, part_shape )
call b%new( num_modes, dr, dxi, nd, nvp, part_shape, p_entity_beam )
call e%new( num_modes, dr, dxi, nd, nvp, part_shape, p_entity_beam )

uq_re => rho%get_rf_re()
uq_im => rho%get_rf_im()

! set the charge
do mode = 0, num_modes
  
  do i = 1, nd(1)
    r = (i-0.5)*dr
    uq_re(mode)%f1(1,i) = exp( -((r-0.5)/0.05)**2 )
  enddo

  if ( mode == 0 ) cycle

  do i = 1, nd(1)
    r = (i-0.5)*dr
    uq_im(mode)%f1(1,i) = exp( -((r-0.5)/0.05)**2 )
  enddo

enddo

call b%solve(rho)
call e%solve(b)

ue_re => e%get_rf_re()
ue_im => e%get_rf_im()

p => ue_re(0)%get_f1()
call write_data( p, 'er-re-beam-0.txt', 1 )
call write_data( p, 'ephi-re-beam-0.txt', 2 )
p => ue_re(1)%get_f1()
call write_data( p, 'er-re-beam-1.txt', 1 )
call write_data( p, 'ephi-re-beam-1.txt', 2 )
p => ue_re(2)%get_f1()
call write_data( p, 'er-re-beam-2.txt', 1 )
call write_data( p, 'ephi-re-beam-2.txt', 2 )
p => ue_im(1)%get_f1()
call write_data( p, 'er-im-beam-1.txt', 1 )
call write_data( p, 'ephi-im-beam-1.txt', 2 )
p => ue_im(2)%get_f1()
call write_data( p, 'er-im-beam-2.txt', 1 )
call write_data( p, 'ephi-im-beam-2.txt', 2 )

call rho%del()
call b%del()
call e%del()

call write_dbg( 'main', 'test_field_eperp', 0, 'ends' )

call end_errors()
call MPI_FINALIZE( ierr )

end program test_field_eperp
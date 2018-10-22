program test_field_bperp_iter

use field_b_class
use field_class
use field_src_class
use system
use param
use mpi
use ufield_class
use debug_tool

implicit none

type( field_b ) :: b
type( field_jay ) :: jay
type( field_djdxi ) :: djdxi

integer :: num_modes = 2, dim = 3, part_shape = p_ps_linear
integer, dimension(2) :: nd = (/128, 1/), nvp = (/1, 1/)
integer, dimension(2,2) :: gc_num
real :: dr, dxi, r

type( ufield ), dimension(:), pointer :: ujay_re => null(), ujay_im => null()
type( ufield ), dimension(:), pointer :: udjdxi_re => null(), udjdxi_im => null()
type( ufield ), dimension(:), pointer :: ub_re_new => null(), ub_im_new => null()
type( ufield ), dimension(:), pointer :: ub_re_old => null(), ub_im_old => null()
type( ufield ) :: ub_res
real, dimension(:,:), pointer :: p
integer :: ierr, i, mode, iter
integer :: num_iter = 10
real, dimension(:,:), pointer :: res => null()

call MPI_INIT( ierr )
call init_errors( 2, 3 )

call write_dbg( 'main', 'test_field_bperp_iter', 0, 'starts' )

dr = 1.0 / (nd(1)-0.5)
dxi = 1.0

gc_num(:,1) = (/0,0/)
gc_num(:,2) = (/0,0/)

call jay%new( num_modes, dr, dxi, nd, nvp, part_shape )
call djdxi%new( num_modes, dr, dxi, nd, nvp, part_shape )
call b%new( num_modes, dr, dxi, nd, nvp, part_shape, p_entity_plasma )


ujay_re => jay%get_rf_re()
ujay_im => jay%get_rf_im()
udjdxi_re => djdxi%get_rf_re()
udjdxi_im => djdxi%get_rf_im()

! set the charge
do mode = 0, num_modes
  
  do i = 1, nd(1)
    r = (i-0.5)*dr
    ujay_re(mode)%f1(1,i) = 0.0
    ujay_re(mode)%f1(2,i) = 0.0
    ujay_re(mode)%f1(3,i) = 0.0
    udjdxi_re(mode)%f1(1,i) = 0.5*exp( -((r-0.4)/0.05)**2 )
    udjdxi_re(mode)%f1(2,i) = 0.5*exp( -((r-0.4)/0.05)**2 )
    udjdxi_re(mode)%f1(3,i) = 0.5*exp( -((r-0.4)/0.05)**2 )
  enddo

  if ( mode == 0 ) cycle

  do i = 1, nd(1)
    r = (i-0.5)*dr
    ujay_im(mode)%f1(1,i) = 0.0
    ujay_im(mode)%f1(2,i) = 0.0
    ujay_im(mode)%f1(3,i) = 0.0
    udjdxi_im(mode)%f1(1,i) = 0.5*exp( -((r-0.4)/0.05)**2 )
    udjdxi_im(mode)%f1(2,i) = 0.5*exp( -((r-0.4)/0.05)**2 )
    udjdxi_im(mode)%f1(3,i) = 0.5*exp( -((r-0.4)/0.05)**2 )
  enddo

enddo

! call HYPRE_StructMatrixGetBoxValues( psi%solver(0)%A, 1, nd(1), 3, &
!   psi%solver(0)%stencil_idx, values, ierr )

! print *, "A = ", values(1:nd(1)*3)

! call HYPRE_StructVectorGetBoxValues( psi%solver(0)%b )

allocate( res(num_iter, (2*num_modes+1)*2) )
allocate( ub_re_old(0:num_modes), ub_im_old(num_modes) )
call ub_res%new( dim, nd, nvp, gc_num, has_2d=.true. )

do i = 0, num_modes
  call ub_re_old(i)%new( dim, nd, nvp, gc_num, has_2d=.true. )
  if (i==0) cycle
  call ub_im_old(i)%new( dim, nd, nvp, gc_num, has_2d=.true. )
enddo

do iter = 1, num_iter

  print *, 'iter = ', iter

  call b%solve( djdxi, jay )
  ub_re_new => b%get_rf_re()
  ub_im_new => b%get_rf_im()

  do i = 0, num_modes
    ub_res = ub_re_new(i) - ub_re_old(i)
    ub_re_old(i) = ub_re_new(i)
    p => ub_res%get_f1()
    if (i==0) then
      res(iter,1) = maxval(abs(p(1,:)))
      res(iter,2) = maxval(abs(p(2,:)))
      cycle
    endif
    res(iter,4*i-1) = maxval(abs(p(1,:)))
    res(iter,4*i)   = maxval(abs(p(2,:)))

    ub_res = ub_im_new(i) - ub_im_old(i)
    ub_im_old(i) = ub_im_new(i)
    p => ub_res%get_f1()
    res(iter,4*i+1) = maxval(abs(p(1,:)))
    res(iter,4*i+2) = maxval(abs(p(2,:)))
  enddo

enddo

write(*,*) "Residue of each iteration"
write(*, '(10A16)') "Re(Br0)","Re(Bp0)","Re(Br1)","Re(Bp1)",&
"Im(Br1)","Im(Bp1)","Re(Br2)","Re(Bp2)","Im(Br2)","Im(Bp2)"
do iter = 1, num_iter
  write(*,'(10E16.4)') res(iter,:)
enddo


call jay%del()
call djdxi%del()
call b%del()
deallocate( ub_re_old, ub_im_old )

call write_dbg( 'main', 'test_field_bperp_iter', 0, 'ends' )

call end_errors()
call MPI_FINALIZE( ierr )

end program test_field_bperp_iter
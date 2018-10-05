module field_b_class

use field_class
use field_solver_class
use ufield_class
use param
use system
use debug_tool

implicit none

private

character(len=20), parameter :: cls_name = "field_b"
integer, parameter :: cls_level = 1

public :: field_b

type, extends( field ) :: field_b

  ! private

  class( field_solver ), dimension(:), pointer :: solver_bz => null()
  class( field_solver ), dimension(:), pointer :: solver_bperp => null()
  class( field_solver ), dimension(:), pointer :: solver_bperp_iter => null()

  real, dimension(:), pointer :: buf_re => null(), buf_im => null()
  real, dimension(:), pointer :: buf => null()

  contains

  generic :: new => init_field_b
  procedure :: del => end_field_b
  ! generic :: read_input => read_input_field_b
  generic :: solve => solve_field_bz, solve_field_bperp_iter

  procedure, private :: init_field_b
  procedure, private :: end_field_b
  procedure, private :: set_source_bz
  ! procedure, private :: set_source_bperp
  procedure, private :: set_source_bperp_iter
  procedure, private :: get_solution_bz
  ! procedure, private :: get_solution_bperp
  procedure, private :: get_solution_bperp_iter
  procedure, private :: solve_field_bz
  ! procedure, private :: solve_field_bperp
  procedure, private :: solve_field_bperp_iter

end type field_b

contains

subroutine init_field_b( this, num_modes, dr, dxi, nd, nvp, part_shape )

  implicit none

  class( field_b ), intent(inout) :: this
  integer, intent(in) :: num_modes, part_shape
  real, intent(in) :: dr, dxi
  integer, intent(in), dimension(2) :: nd, nvp

  integer, dimension(2,2) :: gc_num
  integer :: dim, i
  integer, dimension(2) :: ndp, noff
  real :: tol
  character(len=20), save :: sname = "init_field_b"

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  ndp = nd / nvp
  noff = (/0,0/)
  tol = 1.0d-6

  select case ( part_shape )
  
  case ( p_ps_linear )
  
    gc_num(:,1) = (/0, 1/)
    gc_num(:,2) = (/0, 1/)
  
  case ( p_ps_quadratic )

    print *, "Quadratic particle shape not implemented."
    stop

  case default

    print *, "Invalid particle shape."
    stop

  end select

  dim = 3
  ! call initialization routine of the parent class
  call this%field%new( num_modes, dim, dr, dxi, nd, nvp, gc_num )

  ! initialize solver
  allocate( this%solver_bz( 0:num_modes ) )
  allocate( this%solver_bperp_iter( 0:num_modes ) )
  do i = 0, num_modes
    call this%solver_bz(i)%new( nd, ndp, noff, p_fk_bz, i, dr, p_hypre_cycred, tol )
    call this%solver_bperp_iter(i)%new( nd, ndp, noff, p_fk_bperp_iter, i, dr, p_hypre_amg, tol )
  enddo 

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine init_field_b

subroutine end_field_b( this )

  implicit none

  class( field_b ), intent(inout) :: this

  integer :: i
  character(len=20), save :: sname = 'end_field_b'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  do i = 0, this%num_modes
    call this%solver_bz(i)%del()
    call this%solver_bperp_iter(i)%del()
  enddo
  deallocate( this%solver_bz )
  deallocate( this%solver_bperp_iter )

  if ( associated( this%buf_re ) ) deallocate( this%buf_re )
  if ( associated( this%buf_im ) ) deallocate( this%buf_im )
  if ( associated( this%buf ) ) deallocate( this%buf )

  call this%field%del()

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine end_field_b

subroutine set_source_bz( this, mode, jay_re, jay_im )

  implicit none

  class( field_b ), intent(inout) :: this
  class( ufield ), intent(in) :: jay_re
  class( ufield ), intent(in), optional :: jay_im
  integer, intent(in) :: mode

  integer :: i, nd1p
  real, dimension(:,:), pointer :: f1_re => null(), f1_im => null()
  real :: idrh, idr, a1, a2, a3, b
  character(len=20), save :: sname = 'set_source_bz'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  nd1p = jay_re%get_ndp(1)
  idr = 1.0 / this%dr
  idrh = 0.5 * idr
  
  f1_re => jay_re%get_f1()
  if ( .not. associated( this%buf_re ) ) then
    allocate( this%buf_re( nd1p ) )
  elseif ( size(this%buf_re) < nd1p ) then
    deallocate( this%buf_re )
    allocate( this%buf_re( nd1p ) )
  endif

  if ( present(jay_im) ) then
    f1_im => jay_im%get_f1()
    if ( .not. associated( this%buf_im ) ) then
      allocate( this%buf_im( nd1p ) )
    elseif ( size(this%buf_im) < nd1p ) then
      deallocate( this%buf_im )
      allocate( this%buf_im( nd1p ) )
    endif
  endif

  if ( mode == 0 ) then
    
    do i = 2, nd1p-1

      a1 = idrh * (i-1.0) / (i-0.5)
      a2 = -idrh / (i-0.5)
      a3 = -idrh * i / (i-0.5)

      this%buf_re(i) = a1 * f1_re(2,i-1) + a2 * f1_re(2,i) + a3 * f1_re(2,i+1)

    enddo

    ! calculate the derivatives at the boundary and axis
    this%buf_re(1) = -idr * ( f1_re(2,1) + f1_re(2,2) )
    a2 = idr * (nd1p+0.5) / (nd1p-0.5)
    this%buf_re(nd1p) = -idr * f1_re(2,nd1p-1) + a2 * f1_re(2,nd1p)

    ! call write_data( this%buf_re, 'bsource-re-0.txt' )

  elseif ( mode > 0 .and. present( jay_im ) ) then
    
    do i = 2, nd1p-1

      a1 = idrh * (i-1.0) / (i-0.5)
      a2 = -idrh / (i-0.5)
      a3 = -idrh * i / (i-0.5)
      b  = idr * real(mode) / (i-0.5)

      this%buf_re(i) = a1 * f1_re(2,i-1) + a2 * f1_re(2,i) + a3 * f1_re(2,i+1) - &
                        b * f1_im(1,i)
      this%buf_im(i) = a1 * f1_im(2,i-1) + a2 * f1_im(2,i) + a3 * f1_im(2,i+1) + &
                        b * f1_re(1,i)

    enddo

    ! calculate the derivatives at the boundary and axis
    this%buf_re(1) = -idr * ( f1_re(2,1) + f1_re(2,2) + 2.0 * real(mode) * f1_im(1,1) )
    this%buf_im(1) = -idr * ( f1_im(2,1) + f1_im(2,2) - 2.0 * real(mode) * f1_re(1,1) )
    a2 = idr * (nd1p+0.5) / (nd1p-0.5)
    b  = idr * real(mode) / (nd1p-0.5)
    this%buf_re(nd1p) = -idr * f1_re(2,nd1p-1) + a2 * f1_re(2,nd1p) - b * f1_im(1,nd1p)
    this%buf_im(nd1p) = -idr * f1_im(2,nd1p-1) + a2 * f1_im(2,nd1p) + b * f1_re(1,nd1p)

  else

    call write_err( 'Invalid input arguments!' )

  endif

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine set_source_bz

subroutine set_source_bperp_iter( this, mode, djdxi_re, jay_re, djdxi_im, jay_im )

  implicit none

  class( field_b ), intent(inout) :: this
  class( ufield ), intent(in) :: djdxi_re, jay_re
  class( ufield ), intent(in), optional :: djdxi_im, jay_im
  integer, intent(in) :: mode

  integer :: i, nd1p
  real, dimension(:,:), pointer :: f1_re => null(), f1_im => null()
  real, dimension(:,:), pointer :: f2_re => null(), f2_im => null()
  real, dimension(:,:), pointer :: f3_re => null(), f3_im => null()
  real :: idrh, idr, a1, a2, a3, b, ir
  character(len=20), save :: sname = 'set_source_bperp_iter'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  nd1p = jay_re%get_ndp(1)
  idr = 1.0 / this%dr
  idrh = 0.5 * idr
  
  f1_re => djdxi_re%get_f1()
  f2_re => jay_re%get_f1()
  f3_re => this%rf_re(mode)%get_f1()
  if ( .not. associated( this%buf ) ) then
    allocate( this%buf( nd1p*4 ) )
  endif

  if ( present(djdxi_im) .and. present(jay_im) ) then
    f1_im => djdxi_im%get_f1()
    f2_im => jay_im%get_f1()
    f3_im => this%rf_im(mode)%get_f1()
  endif

  this%buf = 0.0
  if ( mode == 0 ) then
    
    do i = 2, nd1p-1

      ir = idr / (i-0.5)
      ! ! source for Re(Br)
      ! this%buf(i) = -f1_re(2,i) - f3_re(1,i)
      ! ! source for Im(Bphi)
      ! this%buf(i+nd1p) = 0.0
      ! ! source for Im(Br)
      ! this%buf(i+2*nd1p) = 0.0
      ! ! source for Re(Bphi)
      ! this%buf(i+3*nd1p) = f1_re(1,i) - idrh * ( f2_re(3,i+1)-f2_re(3,i-1) ) - f3_re(2,i)


      ! Re(Br)
      this%buf(4*i-3) = -f1_re(2,i) - f3_re(1,i)
      ! Im(Br)
      this%buf(4*i-2) = 0.0
      ! Re(Bphi)
      this%buf(4*i-1) = f1_re(1,i) + idrh * ( f2_re(3,i+1)-f2_re(3,i-1) ) - f3_re(2,i)
      ! Im(Bphi)
      this%buf(4*i) = 0.0

    enddo

    ! calculate the derivatives at the boundary and axis
    this%buf(1) = -f1_re(2,1) - f3_re(1,1)
    this%buf(2) = 0.0
    this%buf(3) = f1_re(1,1) + idr * ( f2_re(3,2)-f2_re(3,1) ) - f3_re(2,1)
    this%buf(4) = 0.0

    this%buf(4*nd1p-3) = -f1_re(2,nd1p) - f3_re(1,nd1p)
    this%buf(4*nd1p-2) = 0.0
    this%buf(4*nd1p-1) = f1_re(1,nd1p) + idr * ( f2_re(3,nd1p)-f2_re(3,nd1p-1) ) - f3_re(2,nd1p)
    this%buf(4*nd1p)   = 0.0

    ! call write_data( this%buf, 'bsource-re-0.txt' )

  elseif ( mode > 0 .and. present( jay_im ) .and. present( djdxi_im ) ) then
    
    do i = 2, nd1p-1

      ir = idr / (i-0.5)
      this%buf(4*i-3) = -f1_re(2,i) + mode * f2_im(3,i) * ir - f3_re(1,i)
      this%buf(4*i-2) = -f1_im(2,i) - mode * f2_re(3,i) * ir - f3_im(1,i)
      this%buf(4*i-1) = f1_re(1,i) + idrh * ( f2_re(3,i+1)-f2_re(3,i-1) ) - f3_re(2,i)
      this%buf(4*i)   = f1_im(1,i) + idrh * ( f2_im(3,i+1)-f2_im(3,i-1) ) - f3_im(2,i)
      
    enddo

    ! calculate the derivatives at the boundary and axis
    this%buf(1) = -f1_re(2,1) + mode * f2_im(3,1) * ir - f3_re(1,1)
    this%buf(2) = -f1_im(2,1) - mode * f2_re(3,1) * ir - f3_im(1,1)
    this%buf(3) = f1_re(1,1) + idr * ( f2_re(3,2)-f2_re(3,1) ) - f3_re(2,1)
    this%buf(4) = f1_im(1,1) + idr * ( f2_im(3,2)-f2_im(3,1) ) - f3_im(2,1)

    this%buf(4*nd1p-3) = -f1_re(2,nd1p) + mode * f2_im(3,nd1p) * ir - f3_re(1,nd1p)
    this%buf(4*nd1p-2) = -f1_im(2,nd1p) - mode * f2_re(3,nd1p) * ir - f3_im(1,nd1p)
    this%buf(4*nd1p-1) = f1_re(1,nd1p) + idr * ( f2_re(3,nd1p)-f2_re(3,nd1p-1) ) - f3_re(2,nd1p)
    this%buf(4*nd1p)   = f1_im(1,nd1p) + idr * ( f2_im(3,nd1p)-f2_im(3,nd1p-1) ) - f3_im(2,nd1p)

  else

    call write_err( 'Invalid input arguments!' )

  endif

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine set_source_bperp_iter

subroutine get_solution_bz( this, mode )

  implicit none

  class( field_b ), intent(inout) :: this
  integer, intent(in) :: mode

  integer :: i, nd1p
  real, dimension(:,:), pointer :: f1_re => null(), f1_im => null()
  character(len=20), save :: sname = 'get_solution_bz'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  nd1p = this%rf_re(mode)%get_ndp(1)

  f1_re => this%rf_re(mode)%get_f1()
  do i = 1, nd1p
    f1_re(3,i) = this%buf_re(i)
  enddo

  if ( mode > 0 ) then
    f1_im => this%rf_im(mode)%get_f1()
    do i = 1, nd1p
      f1_im(3,i) = this%buf_im(i)
    enddo
  endif

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine get_solution_bz

subroutine get_solution_bperp_iter( this, mode )

  implicit none

  class( field_b ), intent(inout) :: this
  integer, intent(in) :: mode

  integer :: i, nd1p
  real, dimension(:,:), pointer :: f1_re => null(), f1_im => null()
  character(len=20), save :: sname = 'get_solution_bperp_iter'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  nd1p = this%rf_re(mode)%get_ndp(1)

  f1_re => this%rf_re(mode)%get_f1()
  do i = 1, nd1p
    f1_re(1,i) = this%buf(4*i-3)
    f1_re(2,i) = this%buf(4*i-1)
  enddo

  if ( mode > 0 ) then
    f1_im => this%rf_im(mode)%get_f1()
    do i = 1, nd1p
      f1_im(1,i) = this%buf(4*i-2)
      f1_im(2,i) = this%buf(4*i)
    enddo
  endif

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine get_solution_bperp_iter

subroutine solve_field_bz( this, jay )

  implicit none

  class( field_b ), intent(inout) :: this
  class( field ), intent(in) :: jay

  type( ufield ), dimension(:), pointer :: jay_re => null(), jay_im => null()
  integer :: i
  character(len=20), save :: sname = 'solve_field_bz'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  jay_re => jay%get_rf_re()
  jay_im => jay%get_rf_im()

  do i = 0, this%num_modes

    if ( i == 0 ) then
      call this%set_source_bz( i, jay_re(i) )
      call this%solver_bz(i)%solve( this%buf_re )
      call this%get_solution_bz(i)
      cycle
    endif

    call this%set_source_bz( i, jay_re(i), jay_im(i) )
    call this%solver_bz(i)%solve( this%buf_re )
    call this%solver_bz(i)%solve( this%buf_im )
    call this%get_solution_bz(i)

  enddo

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine solve_field_bz

subroutine solve_field_bperp_iter( this, djdxi, jay )

  implicit none

  class( field_b ), intent(inout) :: this
  class( field ), intent(in) :: djdxi, jay

  type( ufield ), dimension(:), pointer :: jay_re => null(), jay_im => null()
  type( ufield ), dimension(:), pointer :: djdxi_re => null(), djdxi_im => null()
  integer :: i
  character(len=20), save :: sname = 'solve_field_bperp_iter'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  djdxi_re => djdxi%get_rf_re()
  djdxi_im => djdxi%get_rf_im()
  jay_re => jay%get_rf_re()
  jay_im => jay%get_rf_im()

  do i = 0, this%num_modes

    if ( i == 0 ) then
      call this%set_source_bperp_iter( i, djdxi_re(i), jay_re(i) )
      call this%solver_bperp_iter(i)%solve( this%buf )
      call this%get_solution_bperp_iter(i)
      cycle
    endif

    call this%set_source_bperp_iter( i, djdxi_re(i), jay_re(i), djdxi_im(i), jay_im(i) )
    call this%solver_bperp_iter(i)%solve( this%buf )
    call this%get_solution_bperp_iter(i)

  enddo

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine solve_field_bperp_iter

end module field_b_class
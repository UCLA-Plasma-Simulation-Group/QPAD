module field_solver_all_modes_class

use options_class 
use parallel_module
use mpi
use param
use sysutil_module
use debug_tool
use ufield_class

implicit none

private

public :: field_solver_all_modes
public :: HYPRE_BUF

character(len=32), parameter :: cls_name = "field_solver_all_modes"
integer, parameter :: cls_level = 4

real, dimension(:), pointer, save :: HYPRE_BUF => null()
real, dimension(:), pointer :: src_sol => null()

type :: field_solver_all_modes ! class for HYPRE solver

  ! HYPRE parameters
  integer, dimension(:), pointer :: offsets => null()
  integer, dimension(:), pointer :: stencil_idx => null()
  integer :: num_stencil
  integer :: stype ! currently there's only cyclic reduction
  integer :: kind
  integer :: mode
  integer :: bnd

  integer(HYPRE_TYPE) :: A, b, x, grid, stencil, solver
  integer :: iupper, ilower

  contains

  procedure :: new => init_field_solver
  procedure :: del => end_field_solver 
  procedure :: solve => solve_equation
  procedure, private :: set_struct_solver
  procedure, private :: set_struct_grid
  procedure, private :: set_struct_stencil
  procedure, private :: set_struct_matrix 
 
end type field_solver_all_modes

integer, save :: nofff
real, save :: drr
integer :: nr
 
contains

! =====================================================================
! Class field_solver implementation
! =====================================================================

subroutine init_field_solver( this, opts, max_mode, dr, kind, bnd, stype )

  implicit none

  class( field_solver_all_modes ), intent(inout) :: this
  type( options ), intent(in) :: opts
  integer, intent(in) :: kind, stype, max_mode, bnd
  real, intent(in) :: dr

  integer :: ierr, comm
  character(len=32), save :: sname = "init_field_solver"

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  this%stype = stype
  this%mode  = max_mode
  this%kind  = kind
  this%bnd   = bnd
  nofff = opts%get_noff(1)
  drr = dr
  nr = opts%get_ndp(1)

  ! setup HYPRE grid
  comm = comm_loc()

  call this%set_struct_grid( opts )
  call this%set_struct_stencil()
!   call this%set_struct_matrix( opts, dr )

  call HYPRE_StructVectorCreate( comm, this%grid, this%b, ierr )
  call HYPRE_StructVectorInitialize( this%b, ierr )

  call HYPRE_StructVectorCreate( comm, this%grid, this%x, ierr )
  call HYPRE_StructVectorInitialize( this%x, ierr )

!   call this%set_struct_solver()

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine init_field_solver

subroutine end_field_solver( this )

  implicit none

  class( field_solver_all_modes ), intent(inout) :: this

  integer :: ierr
  character(len=32), save :: sname = "end_field_solver"

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  if ( associated(HYPRE_BUF) ) deallocate( HYPRE_BUF )
  if ( associated(src_sol) ) deallocate( src_sol )

  call HYPRE_StructGridDestroy( this%grid, ierr )
  call HYPRE_StructStencilDestroy( this%stencil, ierr )
  call HYPRE_StructVectorDestroy( this%b, ierr )
  call HYPRE_StructVectorDestroy( this%x, ierr )
  call HYPRE_StructMatrixDestroy( this%A, ierr )

  call HYPRE_StructCycRedDestroy( this%solver, ierr )

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine end_field_solver

subroutine set_struct_solver( this )

  implicit none

  class( field_solver_all_modes ), intent(inout) :: this

  integer :: ierr, comm
  character(len=32), save :: sname = "set_struct_solver"

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  comm = comm_loc()

  write(2,*) 'solver Initializing 2'
  ! call write_stdout( 'mode '//num2str(this%mode)//': Using Cyclic Reduction solver' )
  call HYPRE_StructCycRedCreate( comm, this%solver, ierr )
  write(2,*) 'solver Initializing 3'
  call HYPRE_StructCycRedSetup( this%solver, this%A, this%b, this%x, ierr )
  write(2,*) 'solver Initializing 4'

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine set_struct_solver

subroutine solve_equation( this, src, psi_re, q_re, psi_im, q_im, u, qbm)

  implicit none

  class( field_solver_all_modes ), intent(inout) :: this
  type(ufield), intent(inout), dimension(:), optional, pointer :: psi_re
  type(ufield), intent(inout), dimension(:), optional, pointer :: q_re
  type(ufield), intent(inout), dimension(:), optional, pointer :: psi_im
  type(ufield), intent(inout), dimension(:), optional, pointer :: q_im
  real, intent(inout), dimension(:,:,:), optional, pointer :: src
!   real, intent(inout), dimension(:), optional, pointer :: src1
!   real, intent(inout), dimension(:), optional, pointer :: src
  real, intent(inout), dimension(:,:,:), optional :: u
  real, intent(inout), optional :: qbm

  integer :: ierr, a, im, k, l, kk
  real :: jj
  integer :: size_sol
  real, dimension(:,:), pointer :: f1_re => null(), f1_im => null()
  real, dimension(:,:), pointer :: f2_re => null(), f2_im => null()
  real, dimension(:,:), pointer :: f3_re => null(), f3_im => null()
  character(len=32), save :: sname = "solve_equation"

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  size_sol = nr * (4 * this%mode + 1)  

  if (.not. associated( src_sol)) then
    allocate( src_sol(size_sol) )
  end if

  select case ( this%kind )
  case ( p_fk_coef )
    call set_struct_matrix(this, psi_re = psi_re, psi_im = psi_im, qbm=qbm)
    call this%set_struct_solver()
    write(2,*) 'this%solver_coef%solve||HYPRE_StructVectorSetBoxValues 0' 
    a = 1
    do k = 1, nr
      do l = 0, 2*this%mode
        im = l - this%mode
!         write(2,*) im,"im"
        if (im > 0 ) then
          f1_re => q_re(im)%get_f1()
          f1_im => q_im(im)%get_f1()
!           write(2,*) f1_re(1,k),"f1_re(1,k) im>0"
!           write(2,*) f1_im(1,k),"f1_im(1,k) im>0"
          src_sol(a)  = -qbm*f1_re(1,k)
          src_sol(a + 1)  = -qbm*f1_im(1,k)
          a=a + 2
        elseif (im == 0) then
          f1_re => q_re(0)%get_f1()
!           write(2,*) f1_re(1,k),"f1_re(1,k) im=0"
          src_sol(a)  = -qbm*f1_re(1,k)
          a = a + 1
        else
          im = abs(im)
          f1_re => q_re(im)%get_f1()
          f1_im => q_im(im)%get_f1()
          write(2,*) f1_re(1,k),"f1_re(1,k) im<0"
          write(2,*) f1_im(1,k),"f1_im(1,k) im<0"
          write(2,*) -qbm*f1_re(1,k),"src_sol(a) im<0"
          write(2,*) a,"a0"
          src_sol(a) = -qbm*f1_re(1,k)
          src_sol(a + 1)  = qbm*f1_im(1,k)
!           write(2,*) qbm*f1_im(1,k),"src_sol(a+1) 0  im<0"
!           write(2,*) src_sol(a + 1),"src_sol(a+1) 1  im<0"
          a=a + 2 
!           write(2,*) a,"a"         
        endif
      enddo
    enddo
    write(2,*) a,"a"
    write(2,*) 'this%solver_coef%solve||HYPRE_StructVectorSetBoxValues 1'     
    call HYPRE_StructVectorSetBoxValues( this%b, this%ilower, this%iupper, src_sol, ierr )
    write(2,*) 'this%solver_coef%solve||HYPRE_StructVectorAssembl'   
    call HYPRE_StructVectorAssemble( this%b, ierr )
    write(2,*) 'this%solver_coef%solve||HYPRE_StructCycRedSolve'  
    call HYPRE_StructCycRedSolve( this%solver, this%A, this%b, this%x, ierr )
    write(2,*) 'this%solver_coef%solve||HYPRE_StructVectorGetBoxValues 0' 
    call HYPRE_StructVectorGetBoxValues( this%x, this%ilower, this%iupper, src_sol, ierr )
    a = 1
    do k = 1, nr
      do l = 0, 2*this%mode
        im = l - this%mode
        write(2,*) im,"im get"
        write(2,*) a,"a get"
          if (im == 0 ) then
            if (a > size(src_sol)) then
              print *, "Error: a is out of bounds for src_sol."
              stop
            end if
            write(2,*) src_sol(a),"src_sol(a) im=0 0"
            src(l,k,0) = src_sol(a)
            write(2,*) src_sol(a),"src_sol(a) im=0 1"
            a=a+1
          else
            if (a > size(src_sol)) then
              print *, "Error: a is out of bounds for src_sol."
              stop
            end if
            write(2,*) src_sol(a),"src_sol(a) im><0 0"
            src(l,k,0) = src_sol(a)
            write(2,*) src_sol(a),"src_sol(a) im><0 1"
            src(l,k,1) = src_sol(a + 1)
            a=a+2          
        endif
      enddo 
    enddo
    write(2,*) 'this%solver_coef%solve||HYPRE_StructVectorGetBoxValues 0' 
  case ( p_fk_all_B_minus, p_fk_all_B_plus )
    call set_struct_matrix(this, u=u)
    call this%set_struct_solver() 
    a = 1
    do k = 1, nr
      do l = 0, 2*this%mode
          im = l - this%mode
          if (im > 0 ) then
            src_sol(a)  = src(l,k,0)
            src_sol(a+1)  = src(l,k,1)
            a=a+2
          elseif (im == 0) then
            src_sol(a)  = src(l,k,0)
            a = a+1
          else
            src_sol(a)  = src(l,k,0)
            src_sol(a+1)  = -src(l,k,1)
            a=a+2          
        endif
      enddo 
    enddo
    call HYPRE_StructVectorSetBoxValues( this%b, this%ilower, this%iupper, src_sol, ierr )
    call HYPRE_StructVectorAssemble( this%b, ierr )
    call HYPRE_StructCycRedSolve( this%solver, this%A, this%b, this%x, ierr )
    call HYPRE_StructVectorGetBoxValues( this%x, this%ilower, this%iupper, src_sol, ierr )
    a = 1
    do k = 1, nr
      do l = 0, 2*this%mode
        im = l - this%mode
        if (im == 0 ) then
          src(l,k,0) = src_sol(a)
          a=a+1
        else
          src(l,k,0) = src_sol(a)
          src(l,k,1) = src_sol(a + 1)
          a=a+2          
        endif
      enddo 
    enddo
  end select 

  select case ( this%kind )
  case ( p_fk_coef)
    call stop_tprof( 'solve coef' )
  case ( p_fk_all_B_minus)
    call stop_tprof( 'solve B_minus' )
  case ( p_fk_all_B_plus)
    call stop_tprof( 'solve B_plus' )

  end select

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine solve_equation

subroutine set_struct_grid( this, opts )

  implicit none

  class( field_solver_all_modes ), intent(inout) :: this
  type( options ), intent(in) :: opts

  integer :: comm, ierr
  character(len=32), save :: sname = "set_struct_grid"

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  comm = comm_loc()
  this%ilower = opts%get_noff(1) + 1
  this%iupper = opts%get_noff(1) + opts%get_ndp(1)*(4*this%mode+1) 

  call HYPRE_StructGridCreate( comm, 1, this%grid, ierr )
  call HYPRE_StructGridSetExtents( this%grid, this%ilower, this%iupper, ierr )
  call HYPRE_StructGridAssemble( this%grid, ierr )

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine set_struct_grid

subroutine set_struct_stencil( this )

  implicit none

  class( field_solver_all_modes ), intent(inout) :: this

  integer :: i, j, k, l, ierr
  character(len=32), save :: sname = "set_struct_stencil"

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  select case(this%kind)

  case(p_fk_coef)
    
    this%num_stencil = 4*this%mode + 1
    if ( .not. associated( this%offsets ) ) then
      allocate( this%offsets(this%num_stencil ) )
    endif

    if ( .not. associated( this%stencil_idx ) ) then
      allocate( this%stencil_idx( this%num_stencil ) )
      do i = 1, this%num_stencil
        this%stencil_idx(i) = i-1
      enddo
    endif

    k = 1
    do j = 0, 4*this%mode
      this%offsets(k) = j-2*this%mode
      k=k+1
    enddo

    call HYPRE_StructStencilCreate( 1, this%num_stencil, this%stencil, ierr )
    do i = 1, this%num_stencil
      call HYPRE_StructStencilSetElement( this%stencil, i-1, this%offsets(i), ierr )
    enddo
  
  case( p_fk_all_B_minus, p_fk_all_B_plus )

    this%num_stencil = (4*this%mode + 1)*3
    if ( .not. associated( this%offsets ) ) then
      allocate( this%offsets( this%num_stencil ) )
    endif

    if ( .not. associated( this%stencil_idx ) ) then
      allocate( this%stencil_idx( this%num_stencil ) )
      do i = 1, this%num_stencil
        this%stencil_idx(i) = i-1
      enddo
    endif

    k = 1
    do j = 0, (4*this%mode+1)*3-1
      this%offsets(k) = j-(6*this%mode+1)
      k=k+1
    enddo

    call HYPRE_StructStencilCreate( 1, this%num_stencil, this%stencil, ierr )
    do i = 1, this%num_stencil
      call HYPRE_StructStencilSetElement( this%stencil, i-1, this%offsets(i), ierr )
    enddo

  end select

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine set_struct_stencil

subroutine set_struct_matrix( this, psi_re, psi_im, u, qbm)

  implicit none

  class( field_solver_all_modes ), intent(inout) :: this
  type(ufield), intent(inout), dimension(:), optional, pointer :: psi_re
  type(ufield), intent(inout), dimension(:), optional, pointer :: psi_im
  real, intent(inout), dimension(:,:,:), optional :: u
  real, intent(inout), optional :: qbm

  integer :: i, ierr, local_vol, noff, m, aa, bb, nn, kk
  integer :: k, n, d, bb_offeset, aa_offeset,l, demode, Anumber
  integer :: comm, lidproc, lnvp
  real :: dr2, m2, j, jmax, dr
  real, dimension(:,:), pointer :: f1_re => null(), f1_im => null()
  real, dimension(:,:), pointer :: f2_re => null(), f2_im => null()
  character(len=32), save :: sname = "set_struct_matrix"

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  comm    = comm_loc()
  lidproc = id_proc_loc()
  lnvp    = num_procs_loc()
  noff   = nofff
  dr     = drr

  dr2 = drr*drr
  m = this%mode
  m2 = real(m*m)
  local_vol = (this%iupper - this%ilower + 1)* this%num_stencil

  if ( .not. associated( HYPRE_BUF ) ) then
    allocate( HYPRE_BUF( local_vol ) )
  elseif ( size(HYPRE_BUF) < local_vol ) then
    deallocate( HYPRE_BUF )
    allocate( HYPRE_BUF( local_vol ) )
  endif

  call HYPRE_StructMatrixCreate( comm, this%grid, this%stencil, this%A, ierr )
  call HYPRE_StructMatrixInitialize( this%A, ierr )

  ! set the matrix element and lower boundary
  select case ( this%kind )

  case ( p_fk_coef )

    j = real(noff)
    i = 1
    do nn = 1, nr
      kk = int(j)
      do aa = 0, 4*this%mode
        if( aa - 2*this%mode > 0 ) then
          k = (aa + 1)/2 - this%mode
        else
          k = aa/2 - this%mode
        endif
        do bb = aa - 2*this%mode, aa + 2*this%mode
          if(  (bb >= 0) .and. (bb <= 4*this%mode) ) then
            if( bb - 2*this%mode > 0 ) then
              n = (bb + 1)/2 - this%mode 
            else
              n = bb/2 - this%mode 
            endif
            demode = abs(k-n)
            if (aa == bb) then
              f1_re => psi_re(0)%get_f1()
              HYPRE_BUF(i)   = 1.0 - qbm * f1_re(1,kk)
              i = i + 1
            else
              if (demode > this%mode) then
                HYPRE_BUF(i) = 0.0
                i = i + 1
              else
                if (demode == 0) then
                  HYPRE_BUF(i) = 0.0
                  i = i + 1
                else
                  if (aa - 2*this%mode > 0) then
                    aa_offeset = aa + 1 - (k+this%mode)*2
                  else
                    aa_offeset = aa - (k+this%mode)*2
                  endif
                  if( bb - 2*this%mode > 0 ) then
                    bb_offeset = bb + 1 - (n+this%mode)*2
                  else
                    bb_offeset = bb - (n+this%mode)*2
                  endif
                  if ( aa_offeset == 1  .and. bb_offeset == 1) then
                    f1_re => psi_re(demode)%get_f1()
                    HYPRE_BUF(i) = -qbm * f1_re(1,kk)
                    i = i + 1
                  elseif ( aa_offeset == 1  .and. bb_offeset == 0) then
                    f1_im => psi_im(demode)%get_f1()
                    HYPRE_BUF(i) = -sign(1,k-n)*qbm * f1_im(1,kk)
                    i = i + 1
                  elseif ( aa_offeset == 0  .and. bb_offeset == 1) then
                    f1_im => psi_im(demode)%get_f1()
                    HYPRE_BUF(i) = sign(1,k-n)*qbm * f1_im(1,kk)
                    i = i + 1
                  elseif ( aa_offeset == 0  .and. bb_offeset == 0) then
                    f1_re => psi_re(demode)%get_f1()
                    HYPRE_BUF(i) = -qbm * f1_re(1,kk)
                    i = i + 1
                  endif
                endif
              endif
            endif
          else            
            HYPRE_BUF(i) = 0.0
            i = i + 1
          endif
        enddo
      enddo
      j = j + 1.0 
    enddo
          
    call HYPRE_StructMatrixSetBoxValues( this%A, this%ilower, this%iupper, this%num_stencil, &
      this%stencil_idx, HYPRE_BUF, ierr )

    call HYPRE_StructMatrixAssemble( this%A, ierr )

  case( p_fk_all_B_minus )

    ! set the first grid point of each partition
    i = ((4*this%mode + 1)*3) * (4*this%mode + 1) + 1
    j = real(noff) 
     do nn = 2, nr
      j = j + 1.0
      kk = int(j)
      do aa = 0, 4*this%mode
        if( aa - 2*this%mode > 0 ) then
          k = (aa + 1)/2 - this%mode
        else
          k = aa/2 - this%mode
        endif
        do bb = aa - 6*this%mode - 1, aa + 6*this%mode + 1
          if(  (bb >= 0) .and. (bb <= 4*this%mode) ) then
            if( bb - 2*this%mode > 0 ) then
              n = (bb + 1)/2 - this%mode 
            else
              n = bb/2 - this%mode 
            endif
            demode = abs(k-n)
            if (aa == bb) then
              HYPRE_BUF(i)   = (-2.0 - ((n-1)/j)**2)/dr2 + u(this%mode,kk,0)
              i = i + 1
            else
              if (demode > this%mode) then
                HYPRE_BUF(i) = 0.0
                i = i + 1
              else
                if (demode == 0) then
                  HYPRE_BUF(i) = 0.0
                  i = i + 1
                else
                  if (aa - 2*this%mode > 0) then
                    aa_offeset = aa + 1 - (k+this%mode)*2
                  else
                    aa_offeset = aa - (k+this%mode)*2
                  endif
                  if( bb - 2*this%mode > 0 ) then
                    bb_offeset = bb + 1 - (n+this%mode)*2
                  else
                    bb_offeset = bb - (n+this%mode)*2
                  endif
                  if ( aa_offeset == 1  .and. bb_offeset == 1) then
                    HYPRE_BUF(i) = u(this%mode+k-n,kk,0)
                    i = i + 1
                  elseif ( aa_offeset == 1  .and. bb_offeset == 0) then
                    HYPRE_BUF(i) = u(this%mode+k-n,kk,1)
                    i = i + 1
                  elseif ( aa_offeset == 0  .and. bb_offeset == 1) then
                    HYPRE_BUF(i) = -u(this%mode+k-n,kk,1)
                    i = i + 1
                  elseif ( aa_offeset == 0  .and. bb_offeset == 0) then
                    HYPRE_BUF(i) = u(this%mode+k-n,kk,0)
                    i = i + 1
                  endif
                endif
              endif
            endif
          elseif (bb == aa - 6*this%mode + 1) then
            HYPRE_BUF(i) = 1/dr2 - 1/(j*dr2)
            i = i + 1
          elseif (bb == aa + 6*this%mode + 1) then
            HYPRE_BUF(i) = 1/dr2 - 1/(j*dr2)
            i = i + 1  
          else               
            HYPRE_BUF(i) = 0.0
            i = i + 1
          endif
        enddo
      enddo 
    enddo
! inner boundary
    i = 1
    if (lidproc == 0) then
      do aa = 0, 4*this%mode
        do bb = aa - 6*this%mode - 1, aa + 6*this%mode + 1
          if (aa == bb) then
!             m=-1
            if (aa == 2*this%mode-1) then
              HYPRE_BUF(i) = u(this%mode,0,0) - u(this%mode + 2,0,0) - 4.0/dr2
              i = i + 1
            elseif(aa == 2*this%mode-2) then
              HYPRE_BUF(i) = u(this%mode,0,0) + u(this%mode + 2,0,0) - 4.0/dr2
              i = i + 1
!               m=1
            elseif(aa == 2*this%mode+1) then
              HYPRE_BUF(i) = u(this%mode,0,0) - 4.0/dr2
              i = i + 1
            elseif(aa == 2*this%mode+2) then
              HYPRE_BUF(i) = u(this%mode,0,0) - 4.0/dr2
              i = i + 1
            else
              HYPRE_BUF(i) = 1.0
              i = i + 1
            endif
          else
            HYPRE_BUF(i) = 0.0
            i = i + 1
          endif
        enddo
      enddo 
    else
      j = real(noff) 
      kk = int(j)
      do aa = 0, 4*this%mode
        if( aa - 2*this%mode > 0 ) then
          k = (aa + 1)/2 - this%mode
        else
          k = aa/2 - this%mode
        endif
        do bb = aa - 6*this%mode - 1, aa + 6*this%mode + 1
          if(  (bb >= 0) .and. (bb <= 4*this%mode) ) then
            if( bb - 2*this%mode > 0 ) then
              n = (bb + 1)/2 - this%mode 
            else
              n = bb/2 - this%mode 
            endif
            demode = abs(k-n)
            if (aa == bb) then
              HYPRE_BUF(i)   = (-2.0 - ((n-1)/j)**2)/dr2 + u(this%mode,kk,0)
              i = i + 1
            else
              if (demode > this%mode) then
                HYPRE_BUF(i) = 0.0
                i = i + 1
              else
                if (demode == 0) then
                  HYPRE_BUF(i) = 0.0
                  i = i + 1
                else
                  if (aa - 2*this%mode > 0) then
                    aa_offeset = aa + 1 - (k+this%mode)*2
                  else
                    aa_offeset = aa - (k+this%mode)*2
                  endif
                  if( bb - 2*this%mode > 0 ) then
                    bb_offeset = bb + 1 - (n+this%mode)*2
                  else
                    bb_offeset = bb - (n+this%mode)*2
                  endif
                  if ( aa_offeset == 1  .and. bb_offeset == 1) then
                    HYPRE_BUF(i) = u(this%mode+k-n,kk,0)
                    i = i + 1
                  elseif ( aa_offeset == 1  .and. bb_offeset == 0) then
                    HYPRE_BUF(i) = u(this%mode+k-n,kk,1)
                    i = i + 1
                  elseif ( aa_offeset == 0  .and. bb_offeset == 1) then
                    HYPRE_BUF(i) = -u(this%mode+k-n,kk,1)
                    i = i + 1
                  elseif ( aa_offeset == 0  .and. bb_offeset == 0) then
                    HYPRE_BUF(i) = u(this%mode+k-n,kk,0)
                    i = i + 1
                  endif
                endif
              endif
            endif
          elseif (bb == aa - 6*this%mode + 1) then
            HYPRE_BUF(i) = 1/dr2 - 1/(j*dr2)
            i = i + 1
          elseif (bb == aa + 6*this%mode + 1) then
            HYPRE_BUF(i) = 1/dr2 - 1/(j*dr2)
            i = i + 1  
          else               
            HYPRE_BUF(i) = 0.0
            i = i + 1
          endif
        enddo
      enddo
    endif 

  case( p_fk_all_B_plus )

    ! set the first grid point of each partition
    i = ((4*this%mode + 1)*3) * (4*this%mode + 1) + 1
    j = real(noff) 
     do nn = 2, nr
      j = j + 1.0
      kk = int(j)
      do aa = 0, 4*this%mode
        if( aa - 2*this%mode > 0 ) then
          k = (aa + 1)/2 - this%mode
        else
          k = aa/2 - this%mode
        endif
        do bb = aa - 6*this%mode - 1, aa + 6*this%mode + 1
          if(  (bb >= 0) .and. (bb <= 4*this%mode) ) then
            if( bb - 2*this%mode > 0 ) then
              n = (bb + 1)/2 - this%mode 
            else
              n = bb/2 - this%mode 
            endif
            demode = abs(k-n)
            if (aa == bb) then
              HYPRE_BUF(i)   = (-2.0 - ((n+1)/j)**2)/dr2 + u(this%mode,kk,0)
              i = i + 1
            else
              if (demode > this%mode) then
                HYPRE_BUF(i) = 0.0
                i = i + 1
              else
                if (demode == 0) then
                  HYPRE_BUF(i) = 0.0
                  i = i + 1
                else
                  if (aa - 2*this%mode > 0) then
                    aa_offeset = aa + 1 - (k+this%mode)*2
                  else
                    aa_offeset = aa - (k+this%mode)*2
                  endif
                  if( bb - 2*this%mode > 0 ) then
                    bb_offeset = bb + 1 - (n+this%mode)*2
                  else
                    bb_offeset = bb - (n+this%mode)*2
                  endif
                  if ( aa_offeset == 1  .and. bb_offeset == 1) then
                    HYPRE_BUF(i) = u(this%mode+k-n,kk,0)
                    i = i + 1
                  elseif ( aa_offeset == 1  .and. bb_offeset == 0) then
                    HYPRE_BUF(i) = u(this%mode+k-n,kk,1)
                    i = i + 1
                  elseif ( aa_offeset == 0  .and. bb_offeset == 1) then
                    HYPRE_BUF(i) = -u(this%mode+k-n,kk,1)
                    i = i + 1
                  elseif ( aa_offeset == 0  .and. bb_offeset == 0) then
                    HYPRE_BUF(i) = u(this%mode+k-n,kk,0)
                    i = i + 1
                  endif
                endif
              endif
            endif
          elseif (bb == aa -  6*this%mode + 1) then
            HYPRE_BUF(i) = 1/dr2 - 1/(j*dr2)
            i = i + 1
          elseif (bb == aa +  6*this%mode + 1) then
            HYPRE_BUF(i) = 1/dr2 - 1/(j*dr2)
            i = i + 1  
          else               
            HYPRE_BUF(i) = 0.0
            i = i + 1
          endif
        enddo
      enddo 
    enddo
! inner boundary
    i = 1
    if (lidproc == 0) then
      do aa = 0, 4*this%mode
        do bb = aa - 6*this%mode - 1, aa + 6*this%mode + 1
          if (aa == bb) then
            HYPRE_BUF(i)   = 1.0
            i = i + 1
          else
            HYPRE_BUF(i) = 0.0
            i = i + 1
          endif
        enddo
      enddo 
    else
      j = real(noff) 
      kk = int(j)
      do aa = 0, 4*this%mode
        if( aa - 2*this%mode > 0 ) then
          k = (aa + 1)/2 - this%mode
        else
          k = aa/2 - this%mode
        endif
        do bb = aa - 6*this%mode - 1, aa + 6*this%mode + 1
          if(  (bb >= 0) .and. (bb <= 4*this%mode) ) then
            if( bb - 2*this%mode > 0 ) then
              n = (bb + 1)/2 - this%mode 
            else
              n = bb/2 - this%mode 
            endif
            demode = abs(k-n)
            if (aa == bb) then
              HYPRE_BUF(i)   = (-2.0 - ((n+1)/j)**2)/dr2 + u(this%mode,kk,0)
              i = i + 1
            else
              if (demode > this%mode) then
                HYPRE_BUF(i) = 0.0
                i = i + 1
              else
                if (demode == 0) then
                  HYPRE_BUF(i) = 0.0
                  i = i + 1
                else
                  if (aa - 2*this%mode > 0) then
                    aa_offeset = aa + 1 - (k+this%mode)*2
                  else
                    aa_offeset = aa - (k+this%mode)*2
                  endif
                  if( bb - 2*this%mode > 0 ) then
                    bb_offeset = bb + 1 - (n+this%mode)*2
                  else
                    bb_offeset = bb - (n+this%mode)*2
                  endif
                  if ( aa_offeset == 1  .and. bb_offeset == 1) then
                    HYPRE_BUF(i) = u(this%mode+k-n,kk,0)
                    i = i + 1
                  elseif ( aa_offeset == 1  .and. bb_offeset == 0) then
                    HYPRE_BUF(i) = u(this%mode+k-n,kk,1)
                    i = i + 1
                  elseif ( aa_offeset == 0  .and. bb_offeset == 1) then
                    HYPRE_BUF(i) = -u(this%mode+k-n,kk,1)
                    i = i + 1
                  elseif ( aa_offeset == 0  .and. bb_offeset == 0) then
                    HYPRE_BUF(i) = u(this%mode+k-n,kk,0)
                    i = i + 1
                  endif
                endif
              endif
            endif
          elseif (bb == aa - 6*this%mode + 1) then
            HYPRE_BUF(i) = 1/dr2 - 1/(j*dr2)
            i = i + 1
          elseif (bb == aa + 6*this%mode + 1) then
            HYPRE_BUF(i) = 1/dr2 - 1/(j*dr2)
            i = i + 1  
          else               
            HYPRE_BUF(i) = 0.0
            i = i + 1
          endif
        enddo
      enddo
    endif 

  end select

  ! set the upper boundary
  if ( lidproc == lnvp-1 ) then

    select case ( this%kind )

    case ( p_fk_all_B_minus, p_fk_all_B_plus )
      
      i = local_vol - ((4*this%mode + 1)*3) * (4*this%mode + 1) + 1
      jmax = real(noff + nr)
      do aa = 0, 4*this%mode
        do bb = aa - 6*this%mode - 1, aa + 6*this%mode + 1
!           beta+ becomes 0.0
          if ( bb == aa + 6*this%mode + 1 ) then
            HYPRE_BUF(i) = 0.0
            i = i + 1
          endif
        enddo
      enddo 

    end select

  endif

  call HYPRE_StructMatrixSetBoxValues( this%A, this%ilower, this%iupper, this%num_stencil, &
    this%stencil_idx, HYPRE_BUF, ierr )

  call HYPRE_StructMatrixAssemble( this%A, ierr )

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine set_struct_matrix

end module field_solver_all_modes_class
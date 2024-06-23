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

public :: field_solver
public :: HYPRE_BUF

character(len=32), parameter :: cls_name = "field_solver"
integer, parameter :: cls_level = 4

real, dimension(:,:,:), pointer, save :: HYPRE_BUF => null()

type :: field_solver ! class for HYPRE solver

  ! HYPRE parameters
  integer, dimension(:), pointer :: offsets => null()
  integer, dimension(:), pointer :: stencil_idx => null()
  integer :: num_stencil
  integer :: stype ! currently there's only cyclic reduction
  integer :: kind
  integer :: mode
  integer :: bnd

  integer(HYPRE_TYPE) :: A, b, x, grid, stencil, solver
  integer, dimension(2) :: iupper, ilower

  contains

  procedure :: new => init_field_solver
  procedure :: del => end_field_solver 
  procedure :: solve => solve_equation
  procedure, private :: set_struct_solver
  procedure, private :: set_struct_grid
  procedure, private :: set_struct_stencil
  procedure, private :: set_struct_matrix 
  procedure, private :: set_struct_matrixv1 
 
end type field_solver

integer, save :: nofff
real, save :: drr
 
contains

! =====================================================================
! Class field_solver implementation
! =====================================================================

subroutine init_field_solver( this, opts, max_mode, dr, kind, bnd, stype )

  implicit none

  class( field_solver ), intent(inout) :: this
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

  ! setup HYPRE grid
  comm = comm_loc()

  call this%set_struct_grid( opts )
  call this%set_struct_stencil()
  call this%set_struct_matrix( opts, dr )

  call HYPRE_StructVectorCreate( comm, this%grid, this%b, ierr )
  call HYPRE_StructVectorInitialize( this%b, ierr )

  call HYPRE_StructVectorCreate( comm, this%grid, this%x, ierr )
  call HYPRE_StructVectorInitialize( this%x, ierr )

  call this%set_struct_solver()

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine init_field_solver

subroutine end_field_solver( this )

  implicit none

  class( field_solver ), intent(inout) :: this

  integer :: ierr
  character(len=32), save :: sname = "end_field_solver"

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  if ( associated(HYPRE_BUF) ) deallocate( HYPRE_BUF )

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

  class( field_solver ), intent(inout) :: this

  integer :: ierr, comm
  character(len=32), save :: sname = "set_struct_solver"

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  comm = comm_loc()

  ! call write_stdout( 'mode '//num2str(this%mode)//': Using Cyclic Reduction solver' )
  call HYPRE_StructCycRedCreate( comm, this%solver, ierr )
  call HYPRE_StructCycRedSetup( this%solver, this%A, this%b, this%x, ierr )

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine set_struct_solver

subroutine solve_equation( this, src, psi_re, q_re, psi_im, q_im，u_re, u_im, qbm)

  implicit none

  class( field_solver ), intent(inout) :: this
  type(ufield), intent(inout), pointer, optional :: psi_re
  type(ufield), intent(inout), pointer, optional :: q_re
  type(ufield), intent(inout), pointer, optional :: psi_im
  type(ufield), intent(inout), pointer, optional :: q_im
  real, intent(inout), dimension(:,:,:), pointer, optional :: src
  real, intent(inout), optional :: u_re
  real, intent(inout), optional :: u_im
  real, intent(inout), optional :: qbm
!   real, intent(inout), optional:: psisum
!   real, intent(inout), optional :: qsum

  integer :: ierr, a, im
  real, dimension(:), pointer :: src_sol
  real, dimension(:,:), pointer :: f1_re => null(), f1_im => null()
  real, dimension(:,:), pointer :: f2_re => null(), f2_im => null()
  real, dimension(:,:), pointer :: f3_re => null(), f3_im => null()
  character(len=32), save :: sname = "solve_equation"

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  select case ( this%kind )
  case ( p_fk_coef )
    call set_struct_matrix(this, opts, dr, psi_re, psi_im, qbm)
    call this%set_struct_solver() 
    a = 1
    do k = this%ilower(2),this%iupper(2)
      do l = this%ilower(1),this%iupper(1)
          im = l - this%mode
          if (im > 0 ) then
            f1_re => q_re(im)%get_f1()
            f1_re => q_im(im)%get_f1()
            src_sol(a)  = -qbm*f1_re(1,k)
            src_sol(a+1)  = -qbm*f1_im(1,k)
            a=a+2
          elseif (im = 0) then
            f1_re => q_re(0)%get_f1()
            src_sol(a)  = -qbm*f1_re(1,k)
            a = a+1
          else
            f1_re => q_re(im)%get_f1()
            f1_re => q_im(im)%get_f1()
            src_sol(a)  = -qbm*f1_re(1,k)
            src_sol(a+1)  = qbm*f1_im(1,k)
            a=a+2          
        endif
      enddo 
    enddo
    call HYPRE_StructVectorSetBoxValues( this%b, this%ilower, this%iupper, src_sol, ierr )
    call HYPRE_StructVectorAssemble( this%b, ierr )
    call HYPRE_StructCycRedSolve( this%solver, this%A, this%b, this%x, ierr )
    call HYPRE_StructVectorGetBoxValues( this%x, this%ilower, this%iupper, src_sol, ierr )
    a = 1
    do k = this%ilower(2), this%iupper(2)
      do l = 0, 4*this%mode
          if (l == 2*this%mode ) then
            src(im,k,0) = src_sol(a)
            a=a+1
          else
            src(im,k,0) = src_sol(a)
            src(im,k,1) = src_sol(a + 1)
            a=a+2          
        endif
      enddo 
    enddo

  case (p_fk_all_B_minus,p_fk_all_B_bplus)
    call set_struct_matrix(this, opts, dr, u)
    call this%set_struct_solver() 
    a = 1
    do k = this%ilower(2),this%iupper(2)
      do l = this%ilower(1),this%iupper(1)
          im = l - this%mode
          if (im > 0 ) then
            src_sol(a)  = src(im,k,0)
            src_sol(a+1)  = src(im,k,1)
            a=a+2
          elseif (im = 0) then
            src_sol(a)  = src(im,k,0)
            a = a+1
          else
            src_sol(a)  = src(im,k,0)
            src_sol(a+1)  = -src(im,k,1)
            a=a+2          
        endif
      enddo 
    enddo
    call HYPRE_StructVectorSetBoxValues( this%b, this%ilower, this%iupper, src_sol, ierr )
    call HYPRE_StructVectorAssemble( this%b, ierr )
    call HYPRE_StructCycRedSolve( this%solver, this%A, this%b, this%x, ierr )
    call HYPRE_StructVectorGetBoxValues( this%x, this%ilower, this%iupper, src_sol, ierr )
    a = 1
    do k = this%ilower(2), this%iupper(2)
      do l = 0, 4*this%mode
          if (l == 2*this%mode ) then
            src(im,k,0) = src_sol(a)
            a=a+1
          else
            src(im,k,0) = src_sol(a)
            src(im,k,1) = src_sol(a + 1)
            a=a+2          
        endif
      enddo 
    enddo
  end select 

  select case ( this%kind )
  case ( p_fk_coef)
    call stop_tprof( 'solve coef' )
  end select

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine solve_equation

subroutine set_struct_grid( this, opts )

  implicit none

  class( field_solver ), intent(inout) :: this
  type( options ), intent(in) :: opts

  integer :: comm, ierr
  character(len=32), save :: sname = "set_struct_grid"

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  comm = comm_loc()
  this%ilower(1) = 0
  this%ilower(2) = opts%get_noff(1) + 1
  this%ilower(3) = 0
  this%iupper(1) = 2*this%mode
  this%iupper(2) = opts%get_noff(1) + opts%get_ndp(1)
  this%iupper(3) = 1

  call HYPRE_StructGridCreate( comm, 3, this%grid, ierr )
  call HYPRE_StructGridSetExtents( this%grid, this%ilower, this%iupper, ierr )
  call HYPRE_StructGridAssemble( this%grid, ierr )

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine set_struct_grid

subroutine set_struct_stencil( this )

  implicit none

  class( field_solver ), intent(inout) :: this

  integer :: i, j, k, l, ierr
  character(len=32), save :: sname = "set_struct_stencil"

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  select case(this%kind)

  case(p_fk_coef)
    
    this%num_stencil = 4*this%mode + 1
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
    do j = 0, 2*this%mode
      if (j == this%mode) then
        this%offsets(k) = (/0, 0, 0/)
        k=k+1
      else
        this%offsets(k) = (/j-this%mode, 0, 0/)
        this%offsets(k+1) = (/j-this%mode, 0, 1/)
        k=k+2
      endif
    enddo

    call HYPRE_StructStencilCreate( 3, this%num_stencil, this%stencil, ierr )
    do i = 1, this%num_stencil
      call HYPRE_StructStencilSetElement( this%stencil, i-1, this%offsets(i), ierr )
    enddo
  
  case(p_fk_all_B_bplus,p_fk_all_B_minus)

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
    do l = -1, 1
      do j = 0, 2*this%mode
        if (j == this%mode) then
          this%offsets(k) = (/0, l, 0/)
          k=k+1
        else
          this%offsets(k) = (/j-this%mode, l, 0/)
          this%offsets(k+1) = (/j-this%mode, l, 1/)
          k=k+2
        endif
      enddo
    enddo

    call HYPRE_StructStencilCreate( 3, this%num_stencil, this%stencil, ierr )
    do i = 1, this%num_stencil
      call HYPRE_StructStencilSetElement( this%stencil, i-1, this%offsets(i), ierr )
    enddo

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine set_struct_stencil

subroutine set_struct_matrix( this, opts, dr, psi_re, psi_im, u, qbm)

  implicit none

  class( field_solver ), intent(inout) :: this
  type( options ), intent(in) :: opts
  real, intent(in) :: dr
  type(ufield), intent(inout), pointer, optional :: psi_re
  type(ufield), intent(inout), pointer, optional :: q_re
  type(ufield), intent(inout), pointer, optional :: psi_im
  type(ufield), intent(inout), pointer, optional :: q_im
  real, intent(inout), dimension(:,:,:), pointer, optional :: src
  real, intent(inout), optional :: u
  real, intent(inout), optional :: qbm

  integer :: i, ierr, local_vol, nr, noff, m, aa, bb, cc
  integer :: k, n, d, bb_offeset, cc_offeset,l
  integer :: comm, lidproc, lnvp
  real :: dr2, m2, j, jmax
  real, intent(inout), dimension(:,:,:), pointer :: src_sol
  real, dimension(:,:), pointer :: f1_re => null(), f1_im => null()
  real, dimension(:,:), pointer :: f2_re => null(), f2_im => null()
  real, dimension(:,:), pointer :: f3_re => null(), f3_im => null()
  character(len=32), save :: sname = "set_struct_matrix"

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  comm    = comm_loc()
  lidproc = id_proc_loc()
  lnvp    = num_procs_loc()
  noff    = opts%get_noff(1)
  nofff   = noff
  drr     = dr

  dr2 = dr*dr
  m = this%mode
  m2 = real(m*m)
  nr = this%iupper(2) - this%ilower(2) + 1
  local_vol = nr * this%num_stencil * (4 * this%mode +1)

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

    ! set from the first grid point of each partition，mode is close.
!     i = 4*this%mode+1
      do aa = 1, nr
        do bb = 0, 4*this%mode
          do cc = 0, 4*this%mode
            if( cc - 2*this%mode > 0 ) then
              n = (cc + 1)/2 - this%mode
            else
              n = cc/2 - this%mode
            endif
            if( bb - 2*this%mode > 0 ) then
              k = (bb + 1)/2 - this%mode
            else
              k = bb/2 - this%mode
            endif
            dmode = abs(k-n)
            if (bb == cc) then
              f1_re => psi_re(0)%get_f1()
              HYPRE_BUF(i)   = 1.0 - qbm * f1(1,aa)
              i = i + 1
            else
              if (dmode > this%mode) then
                HYPRE_BUF(i) = 0
                i = i + 1
              else
                if (demode = 0) then
                  HYPRE_BUF(i) = 0
                  i = i + 1
                else
                  bb_offeset = bb - (k+this%mode)*2
                  cc_offeset = cc - (n+this%mode)*2
                  if ( bb_offeset == 1 ) .and. (cc_offeset == 1) then
                    f1_re => psi_re(demode)%get_f1()
                    HYPRE_BUF(i) = -qbm * f1(1,aa)
                    i = i + 1
                  elseif ( bb_offeset == 1 ) .and. (cc_offeset == 0) then
                    f1_im => psi_im(demode)%get_f1()
                    HYPRE_BUF(i) = -sign(1,k-n)*qbm * f1(1,aa)
                    i = i + 1
                  elseif ( bb_offeset == 0 ) .and. (cc_offeset == 1) then
                    f1_im => psi_im(demode)%get_f1()
                    HYPRE_BUF(i) = sign(1,k-n)*qbm * f1(1,aa)
                    i = i + 1
                  elseif ( bb_offeset == 0 ) .and. (cc_offeset == 0) then
                    f1_re => psi_re(demode)%get_f1()
                    HYPRE_BUF(i) = -qbm * f1(1,aa)
                    i = i + 1
                  endif
                endif
              endif
            endif
          enddo
        enddo
      enddo
    call HYPRE_StructMatrixSetBoxValues( this%A, this%ilower, this%iupper, this%num_stencil, &
      this%stencil_idx, HYPRE_BUF, ierr )

    call HYPRE_StructMatrixAssemble( this%A, ierr )

  case(p_fk_all_B_minus)

    ! set the first grid point of each partition
    i=1
    if (lidproc == 0) then
      do bb = 0, 4*this%mode
        do l = -1, 1
          do cc = 0, 4*this%mode
            if( cc - 2*this%mode > 0 ) then
              n = (cc + 1)/2 - this%mode
            else
              n = cc/2 - this%mode
            endif
            if (n == -1) .or. (n == 1) then
              if (l = -1) then
                HYPRE_BUF(i) = 0.0
                i = i + 1
              elseif (l = 0 ) then
                HYPRE_BUF(i) = -4.0 - dr2
                i = i + 1
              elseif (l = 0 ) then
                HYPRE_BUF(i) = 0.0
                i = i + 1
              endif
            else
              if (l = -1) then
                HYPRE_BUF(i) = 0.0
                i = i + 1
              elseif (l = 0 ) then
                HYPRE_BUF(i) = 1.0
                i = i + 1
              elseif (l = 0 ) then
                HYPRE_BUF(i) = 0.0
                i = i + 1
              endif           
            endif
          enddo
        enddo
      enddo
    else
      j = real(noff)
      do aa = 1, nr
        j = j + 1.0
        do bb = 0, 4*this%mode
          do l = -1, 1
            do cc = 0, 4*this%mode
              if( cc - 2*this%mode > 0 ) then
                n = (cc + 1)/2 - this%mode
              else
                n = cc/2 - this%mode
              endif
              if( bb - 2*this%mode > 0 ) then
                k = (bb + 1)/2 - this%mode
              else
                k = bb/2 - this%mode
              endif
              dmode = abs(k-n)
              if (l == -1) then
                if( bb == cc) then
                  HYPRE_BUF(i) = ( 1.0 - 0.5 / j )/dr2
                  i = i + 1
                else
                  HYPRE_BUF(i) = 0
                  i = i + 1 
                endif                 
              elseif( l == 0) then  
                if (bb == cc) then
                  HYPRE_BUF(i)   = (-2.0 - ((n-1)/j)**2)/dr2 + u(this%mode,aa,0)
                  i = i + 1
                else
                  if (dmode > this%mode) then
                    HYPRE_BUF(i) = 0
                    i = i + 1
                  else
                    if (demode = 0) then
                      HYPRE_BUF(i) = 0
                      i = i + 1
                    else
                      bb_offeset = bb - (k+this%mode)*2
                      cc_offeset = cc - (n+this%mode)*2
                      if ( bb_offeset == 1 ) .and. (cc_offeset == 1) then
                        HYPRE_BUF(i) = u(k-n+this%mode,aa,0)
                        i = i + 1
                      elseif ( bb_offeset == 1 ) .and. (cc_offeset == 0) then)
                        HYPRE_BUF(i) = u(k-n+this%mode,aa,1)
                        i = i + 1
                      elseif ( bb_offeset == 0 ) .and. (cc_offeset == 1) then
                        HYPRE_BUF(i) = -u(k-n+this%mode,aa,0)
                        i = i + 1
                      elseif ( bb_offeset == 0 ) .and. (cc_offeset == 0) then
                        HYPRE_BUF(i) = u(k-n+this%mode,aa,0)
                        i = i + 1
                      endif
                    endif
                  endif
                elseif (l== 1) then
                  if( bb == cc) then
                    HYPRE_BUF(i) = ( 1.0 + 0.5 / j ) / dr2
                    i = i + 1
                  else
                    HYPRE_BUF(i) = 0
                    i = i + 1 
                  endif 
                endif 
              endif
            enddo
          enddo
        enddo
      enddo
  case(p_fk_all_B_bplus)

    ! set the first grid point of each partition
    i=1
    if (lidproc == 0) then
      do bb = 0, 4*this%mode
        do l = -1, 1
          do cc = 0, 4*this%mode
            if( cc - 2*this%mode > 0 ) then
              n = (cc + 1)/2 - this%mode
            else
              n = cc/2 - this%mode
            endif
            if (n == -1) .or. (n == 1) then
              if (l = -1) then
                HYPRE_BUF(i) = 0.0
                i = i + 1
              elseif (l = 0 ) then
                HYPRE_BUF(i) = -4.0 - dr2
                i = i + 1
              elseif (l = 0 ) then
                HYPRE_BUF(i) = 0.0
                i = i + 1
              endif
            else
              if (l = -1) then
                HYPRE_BUF(i) = 0.0
                i = i + 1
              elseif (l = 0 ) then
                HYPRE_BUF(i) = 1.0
                i = i + 1
              elseif (l = 0 ) then
                HYPRE_BUF(i) = 0.0
                i = i + 1
              endif           
            endif
          enddo
        enddo
      enddo
    else
      j = real(noff)
      do aa = 1, nr
        j = j + 1.0
        do bb = 0, 4*this%mode
          do l = -1, 1
            do cc = 0, 4*this%mode
              if( cc - 2*this%mode > 0 ) then
                n = (cc + 1)/2 - this%mode
              else
                n = cc/2 - this%mode
              endif
              if( bb - 2*this%mode > 0 ) then
                k = (bb + 1)/2 - this%mode
              else
                k = bb/2 - this%mode
              endif
              dmode = abs(k-n)
              if (l == -1) then
                if( bb == cc) then
                  HYPRE_BUF(i) = ( 1.0 - 0.5 / j )/dr2
                  i = i + 1
                else
                  HYPRE_BUF(i) = 0
                  i = i + 1 
                endif                 
              elseif( l == 0) then  
                if (bb == cc) then
                  HYPRE_BUF(i)   = (-2.0 - ((n+1)/j)**2)/dr2 + u(this%mode,aa,0)
                  i = i + 1
                else
                  if (dmode > this%mode) then
                    HYPRE_BUF(i) = 0
                    i = i + 1
                  else
                    if (demode = 0) then
                      HYPRE_BUF(i) = 0
                      i = i + 1
                    else
                      bb_offeset = bb - (k+this%mode)*2
                      cc_offeset = cc - (n+this%mode)*2
                      if ( bb_offeset == 1 ) .and. (cc_offeset == 1) then
                        HYPRE_BUF(i) = u(k-n+this%mode,aa,0)
                        i = i + 1
                      elseif ( bb_offeset == 1 ) .and. (cc_offeset == 0) then)
                        HYPRE_BUF(i) = u(k-n+this%mode,aa,1)
                        i = i + 1
                      elseif ( bb_offeset == 0 ) .and. (cc_offeset == 1) then
                        HYPRE_BUF(i) = -u(k-n+this%mode,aa,0)
                        i = i + 1
                      elseif ( bb_offeset == 0 ) .and. (cc_offeset == 0) then
                        HYPRE_BUF(i) = u(k-n+this%mode,aa,0)
                        i = i + 1
                      endif
                    endif
                  endif
                elseif (l== 1) then
                  if( bb == cc) then
                    HYPRE_BUF(i) = ( 1.0 + 0.5 / j ) / dr2
                    i = i + 1
                  else
                    HYPRE_BUF(i) = 0
                    i = i + 1 
                  endif 
                endif 
              endif
            enddo
          enddo
        enddo
      enddo

  ! set the upper boundary
  if ( lidproc == lnvp-1 ) then

    select case ( this%bnd )

    ! to be deleted
    case ( p_bnd_zero )
      
      i = local_vol - (4*this%mode+1)^2*3*2
      do bb = 0, 4*this%mode
        do l = 0, 1
          do cc = 0, 4*this%mode
            HYPRE_BUF(i) = 0.0
            i = i + 1
          enddo
        enddo
      enddo

    case ( p_bnd_open )

      jmax = real(noff + nr)

      select case ( this%kind )

      case ( p_fk_all_B_minus )
        
        i = local_vol - (4*this%mode+1)^2*3*2
        do bb = 0, 4*this%mode
          do l = 0, 1
            do cc = 0, 4*this%mode
              if( cc - 2*this%mode > 0 ) then
                n = (cc + 1)/2 - this%mode
              else
                n = cc/2 - this%mode
              endif
              if (l = 0 ) then
                HYPRE_BUF(i)=HYPRE_BUF(i) + (1.0-(n+this%mode-1)/jmax) * HYPRE_BUF(i)
                i = i + 1
              elseif (l = 1 ) then
                HYPRE_BUF(i) = 0.0
                i = i + 1
              endif          
            enddo
          enddo
        enddo

      case ( p_fk_all_B_bplus )
        
        i = local_vol - (4*this%mode+1)^2*3*2
        do bb = 0, 4*this%mode
          do l = 0, 1
            do cc = 0, 4*this%mode
              if( cc - 2*this%mode > 0 ) then
                n = (cc + 1)/2 - this%mode
              else
                n = cc/2 - this%mode
              endif
              if (l = 0 ) then
                HYPRE_BUF(i)=HYPRE_BUF(i) + (1.0-(n+this%mode-1)/jmax) * HYPRE_BUF(i)
                i = i + 1
              elseif (l = 1 ) then
                HYPRE_BUF(i) = 0.0
                i = i + 1
              endif          
            enddo
          enddo
        enddo

      case default
        call write_err( 'Invalid field type!' )
      end select

    case default
      call write_err( 'Invalid boundary condition!' )
    end select

  endif

  call HYPRE_StructMatrixSetBoxValues( this%A, this%ilower, this%iupper, this%num_stencil, &
    this%stencil_idx, HYPRE_BUF, ierr )

  call HYPRE_StructMatrixAssemble( this%A, ierr )

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine set_struct_matrix

end module field_solver_all_modes_class

module field_solver_all_modesv1_class

use options_class 
use parallel_module
use mpi
use param
use sysutil_module
use debug_tool
use ufield_class

implicit none

private

public :: field_solver_all_modesv1
public :: HYPRE_BUF

character(len=32), parameter :: cls_name = "field_solver_all_modesv1"
integer, parameter :: cls_level = 4

real, dimension(:), pointer, save :: HYPRE_BUF => null()
real, dimension(:), pointer :: src_sol => null()

type :: field_solver_all_modesv1 ! class for HYPRE solver

  ! HYPRE parameters
  integer, dimension(:), pointer :: offsets => null()
  integer, dimension(:), pointer :: stencil_idx => null()
  integer, dimension(:), pointer :: indices => null()
  integer :: num_stencil
  integer :: stype ! currently there's only grmes reduction
  integer :: kind
  integer :: mode
  integer :: num_rows
  integer :: bnd

!   integer(HYPRE_TYPE) :: A, b, x, grid, stencil, solver
  integer(HYPRE_TYPE) :: A, A_internal, b, b_internal, x, x_internal, solver, precond
  integer :: iupper, ilower

  contains

  procedure :: new => init_field_solver
  procedure :: del => end_field_solver 
  procedure :: solve => solve_equation
  procedure, private :: set_IJ_solver
  procedure, private :: set_IJ_matrix 
  procedure, private :: set_hypre_src_coef 
  procedure, private :: set_hypre_src_b 
  procedure, private :: get_hypre_src_coef 
  procedure, private :: get_hypre_src_b 
 
end type field_solver_all_modesv1

integer, save :: nofff
real, save :: drr
integer :: nr
 
contains

! =====================================================================
! Class field_solver implementation
! =====================================================================

subroutine init_field_solver( this, opts, max_mode, dr, kind, bnd, stype )

  implicit none

  class( field_solver_all_modesv1 ), intent(inout) :: this
  type( options ), intent(in) :: opts
  integer, intent(in) :: kind, stype, max_mode, bnd
  real, intent(in) :: dr

  integer :: ierr, comm, i
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
  this%ilower = opts%get_noff(1) + 1
  this%iupper = opts%get_noff(1) + opts%get_ndp(1)*(4*this%mode+1)
  this%num_rows = this%iupper - this%ilower + 1
  ! 分配并初始化 indices 数组
  allocate(this%indices(this%num_rows))
  do i = 1, this%num_rows
      this%indices(i) = this%ilower + i - 1
  end do

  call HYPRE_IJVectorCreate(comm, this%ilower, this%iupper, this%b, ierr)
  call HYPRE_IJVectorSetObjectType(this%b, HYPRE_PARCSR, ierr)
  call HYPRE_IJVectorInitialize(this%b, ierr)

  call HYPRE_IJVectorCreate( comm, this%ilower, this%iupper, this%x, ierr)
  call HYPRE_IJVectorSetObjectType(this%x, HYPRE_PARCSR, ierr)
  call HYPRE_IJVectorInitialize(this%x, ierr)

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine init_field_solver

subroutine end_field_solver( this )

  implicit none

  class( field_solver_all_modesv1 ), intent(inout) :: this

  integer :: ierr
  character(len=32), save :: sname = "end_field_solver"

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  call HYPRE_IJMatrixDestroy(this%A)
  call HYPRE_IJVectorDestroy(this%b)
  call HYPRE_IJVectorDestroy(this%x)
  call HYPRE_ParCSRGMRESDestroy(this%solver)
  call HYPRE_BoomerAMGDestroy(this%precond) 

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine end_field_solver

subroutine set_IJ_solver( this )

  implicit none

  class( field_solver_all_modesv1 ), intent(inout) :: this

  integer :: ierr, comm, max_iter
  real :: tol
  external HYPRE_BoomerAMGSolve, HYPRE_BoomerAMGSetup
  character(len=32), save :: sname = "set_struct_solver"

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  comm = comm_loc()

  tol = 1.0e-7
  max_iter = 100
  call HYPRE_ParCSRGMRESCreate(comm, this%solver)
!   call HYPRE_GMRESSetTol(this%solver, tol)
!   call HYPRE_GMRESSetMaxIter(this%solver, max_iter)

  ! 设置预处理器（BoomerAMG 作为示例）
  call HYPRE_BoomerAMGCreate(this%precond)
  call HYPRE_BoomerAMGSetCoarsenType(this%precond, 6)
  call HYPRE_BoomerAMGSetRelaxType(this%precond, 3)
!   call HYPRE_GMRESSetPrecond(this%solver, HYPRE_BoomerAMGSolve, HYPRE_BoomerAMGSetup, this%precond)
!   call HYPRE_ParCSRGMRESSetup(this%solver, this%A, this%b, this%x)

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine set_IJ_solver

subroutine solve_equation( this, src, psi_re, qe_re, qn1_re, psi_im, qe_im, qn1_im, u, qbm)

  implicit none

  class( field_solver_all_modesv1 ), intent(inout) :: this
  type(ufield), intent(inout), dimension(:), optional, pointer :: psi_re
  type(ufield), intent(inout), dimension(:), optional, pointer :: qe_re
  type(ufield), intent(inout), dimension(:), optional, pointer :: qn1_re
  type(ufield), intent(inout), dimension(:), optional, pointer :: psi_im
  type(ufield), intent(inout), dimension(:), optional, pointer :: qe_im
  type(ufield), intent(inout), dimension(:), optional, pointer :: qn1_im
  real, intent(inout), dimension(:,:,:), optional, pointer :: src
!   real, intent(inout), dimension(:), optional, pointer :: src1
!   real, intent(inout), dimension(:), optional, pointer :: src
  real, intent(inout), dimension(:,:,:), optional :: u
  real, intent(inout), optional :: qbm

  integer :: ierr, a, im, k, l, kk
  real :: jj, j
!   real, dimension(:), pointer :: value=> null()
  real(8), dimension(this%num_rows) :: value
  integer :: size_sol
  real, dimension(:,:), pointer :: f1_re => null(), f1_im => null()
  real, dimension(:,:), pointer :: f2_re => null(), f2_im => null()
  real, dimension(:,:), pointer :: f3_re => null(), f3_im => null()
  character(len=32), save :: sname = "solve_equation"

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  size_sol = nr* (4 * this%mode + 1)
!   write(2,*) size_sol, 'size_sol' 

  if (.not. associated( src_sol)) then
    allocate( src_sol(size_sol) )
  end if
  src_sol(:) = 0.0

  value(:) = 0.0

!   write(2,*) this%kind, 'kind' 
  select case ( this%kind )
  case ( p_fk_coef )

    write(2,*) 'set_IJ_matrix'
    call set_IJ_matrix(this, psi_re = psi_re, psi_im = psi_im, qbm=qbm)
    write(2,*) 'set_IJ_solver'
    call this%set_IJ_solver()
    write(2,*) 'set_hypre_src_coef'
    call this%set_hypre_src_coef( qe_re, qn1_re, qe_im, qn1_im, qbm )
    call HYPRE_IJVectorAssemble(this%b)
    write(2,*) 'HYPRE_IJVectorSetValues'
    call HYPRE_IJVectorSetValues(this%x, this%num_rows, this%indices, value)
    call HYPRE_IJVectorAssemble(this%x)
    ! 获取 HYPRE 内部的 ParCSR 向量对象
    call HYPRE_IJVectorGetObject(this%x, this%x_internal)
    call HYPRE_IJVectorGetObject(this%b, this%b_internal)
    write(2,*) 'ParCSRGMRESSolve' 
    call HYPRE_ParCSRGMRESSetup(this%solver, this%A, this%b_internal, this%x_internal)
    call HYPRE_ParCSRGMRESSolve(this%solver, this%A, this%b_internal, this%x_internal)
    call this%get_hypre_src_coef( src )
 
  case ( p_fk_all_B_minus, p_fk_all_B_plus )
    call set_IJ_matrix(this, u=u)
    call this%set_IJ_solver() 
    call this%set_hypre_src_b(src)
    call HYPRE_IJVectorAssemble(this%b)
    ! 设置初始解向量 x 的值为 0.0
    call HYPRE_IJVectorSetValues(this%x, this%num_rows, this%indices, value)
    call HYPRE_IJVectorAssemble(this%x)
    ! 获取 HYPRE 内部的 ParCSR 向量对象
    call HYPRE_IJVectorGetObject(this%x, this%x_internal)
    call HYPRE_IJVectorGetObject(this%b, this%b_internal)
    call HYPRE_ParCSRGMRESSetup(this%solver, this%A, this%b_internal, this%x_internal)
    call HYPRE_ParCSRGMRESSolve(this%solver, this%A, this%b_internal, this%x_internal)
    call this%get_hypre_src_b(src)
 
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

subroutine set_IJ_matrix( this, psi_re, psi_im, u, qbm)

  implicit none

  class( field_solver_all_modesv1 ), intent(inout) :: this
  type(ufield), intent(inout), dimension(:), optional, pointer :: psi_re
  type(ufield), intent(inout), dimension(:), optional, pointer :: psi_im
  real, intent(inout), dimension(:,:,:), optional :: u
  real, intent(inout), optional :: qbm

  integer :: i, ii, ierr, local_vol, noff, m, aa, bb, nn, kk, mm, g, jj
  integer :: k, n, d, bb_offeset, aa_offeset,l, demode, Anumber
  integer :: comm, lidproc, lnvp
  real :: dr2, m2, j, jmax, dr
  real, dimension(1) :: value
  integer, dimension(1) :: row, col
  integer, dimension(:),allocatable :: rows, cols, ncols
  real, dimension(:,:), pointer :: f1_re => null(), f1_im => null()
  real, dimension(:,:), pointer :: f2_re => null(), f2_im => null()
  character(len=32), save :: sname = "set_ij_matrix"

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  write(2,*) 'set_IJ_matrix 0'
  comm    = comm_loc()
  lidproc = id_proc_loc()
  lnvp    = num_procs_loc()
  noff   = nofff
  dr     = drr

  write(2,*) 'set_IJ_matrix 01'
  dr2 = drr*drr
  m = this%mode
  m2 = real(m*m)

  ! 创建 HYPRE_IJMatrix
  write(2,*) 'set_IJ_matrix 1'
  call HYPRE_IJMatrixCreate(comm, this%ilower, this%iupper, this%ilower, this%iupper, this%A, ierr)
  write(2,*) 'set_IJ_matrix 2'
  call HYPRE_IJMatrixSetObjectType(this%A, HYPRE_PARCSR,ierr)
  write(2,*) 'set_IJ_matrix 3'
  call HYPRE_IJMatrixInitialize(this%A,ierr)

  ! set the matrix element and lower boundary
  select case ( this%kind )

  case ( p_fk_coef )

    local_vol = (4*m + 1)*nr
    write(2,*) 'set_IJ_matrix 4'

    if ( .not. associated( HYPRE_BUF ) ) then
      allocate( HYPRE_BUF( local_vol*(4*m +1) ) )
    elseif ( size(HYPRE_BUF) < local_vol*(4*m +1) ) then
      deallocate( HYPRE_BUF )
      allocate( HYPRE_BUF( local_vol*(4*m +1) ) )
    endif

    write(2,*) 'set_IJ_matrix 5'
    allocate(rows(local_vol))
    allocate(cols(local_vol*(4*m +1)))
    allocate(ncols(local_vol))
      
    i = 1 
    ii = 1
    jj = 1
    do nn = 1, nr
      do aa = 0, 4*m
        if( aa - 2*m > 0 ) then
          k = (aa + 1)/2 - m
        else
          k = aa/2 - m
        endif
        rows(ii) = (nn - 1) *(4*m +1) + aa + 1
        do bb = 0, 4*m
          cols(jj) = (nn - 1) *(4*m +1) + aa + bb + 1
          if( bb - 2*m > 0 ) then
            n = (bb + 1)/2 - m
          else
            n = bb/2 - m
          endif
          demode = abs(k-n)
          if (aa == bb) then
            f1_re => psi_re(0)%get_f1()
            HYPRE_BUF(i)   = 1.0 - qbm * f1_re(1,nn)
            i = i + 1
          else
            if (demode > m) then
              HYPRE_BUF(i) = 0.0
              i = i + 1
            else
              if (demode == 0) then
                HYPRE_BUF(i) = 0.0
                i = i + 1
              else
                if (aa - 2*m > 0) then
                  aa_offeset = aa + 1 - (k+m)*2
                else
                  aa_offeset = aa - (k+m)*2
                endif
                if( bb - 2*m > 0 ) then
                  bb_offeset = bb + 1 - (n+m)*2
                else
                  bb_offeset = bb - (n+m)*2
                endif
                if ( aa_offeset == 1  .and. bb_offeset == 1) then
                  f1_re => psi_re(demode)%get_f1()
                  HYPRE_BUF(i) = -qbm * f1_re(1,nn)
                  i = i + 1
                elseif ( aa_offeset == 1  .and. bb_offeset == 0) then
                  f1_im => psi_im(demode)%get_f1()
                  HYPRE_BUF(i) = -sign(1,k-n)*qbm * f1_im(1,nn)
                  i = i + 1
                elseif ( aa_offeset == 0  .and. bb_offeset == 1) then
                  f1_im => psi_im(demode)%get_f1()
                  HYPRE_BUF(i) = sign(1,k-n)*qbm * f1_im(1,nn)
                  i = i + 1
                elseif ( aa_offeset == 0  .and. bb_offeset == 0) then
                  f1_re => psi_re(demode)%get_f1()
                  HYPRE_BUF(i) = -qbm * f1_re(1,nn)
                  i = i + 1
                endif
              endif
            endif
          endif
          jj = jj +1
        enddo
        ii = ii + 1
      enddo 
    enddo

  case( p_fk_all_B_minus )

    local_vol = (4*m + 1)*nr

    if ( .not. associated( HYPRE_BUF ) ) then
      allocate( HYPRE_BUF( local_vol*(12*m +3) ) )
    elseif ( size(HYPRE_BUF) < local_vol ) then
      deallocate( HYPRE_BUF )
      allocate( HYPRE_BUF( local_vol ) )
    endif

    allocate(rows(local_vol))
    allocate(cols(local_vol*(12*m +3)))
    allocate(ncols(local_vol))

    ! set the first grid point of each partition
    j = real(noff)
    i = 1 + (12*m + 3)*(4*m + 1)
    ii = 1 + (4*m + 1)
    jj = 1 + (12*m + 3)*(4*m + 1)
     do nn = 2, nr
      j = j + 1.0
      do aa = 0, 4*m
        if( aa - 2*m > 0 ) then
          k = (aa + 1)/2 - m
        else
          k = aa/2 - m
        endif
        rows(ii) = (nn - 1) *(4*m +1) + aa + 1
        do bb = - 4*m - 1, 8*m + 1
          cols(jj) = (nn - 1) *(4*m +1) + aa + bb + 1
          if(  (bb >= 0) .and. (bb <= 4*m) ) then
            if( bb - 2*m > 0 ) then
              n = (bb + 1)/2 - m
            else
              n = bb/2 - m 
            endif
            demode = abs(k-n)
            if (aa == bb) then
              HYPRE_BUF(i)   = (-2.0 - ((n -1)/j)**2)/dr2 + u(1+m,nn,1)
              i = i + 1
            else
              if (demode > m) then
                HYPRE_BUF(i) = 0.0
                i = i + 1
              else
                if (demode == 0) then
                  HYPRE_BUF(i) = 0.0
                  i = i + 1
                else
                  if (aa - m > 0) then
                    aa_offeset = aa + 1 - (k+m)*2
                  else
                    aa_offeset = aa - (k+m)*2
                  endif
                  if( bb - 2*m > 0 ) then
                    bb_offeset = bb + 1 - (n+m)*2
                  else
                    bb_offeset = bb - (n+m)*2
                  endif
                  if ( aa_offeset == 1  .and. bb_offeset == 1) then
                    HYPRE_BUF(i) = u(1+m+k-n,nn,1)
!                     write(2,*) u(1+this%mode+k-n,kk,1),"u(1+this%mode+k-n,kk,1)"
                    i = i + 1
                  elseif ( aa_offeset == 1  .and. bb_offeset == 0) then
                    HYPRE_BUF(i) = u(1+m+k-n,nn,2)
!                     write(2,*) u(1+this%mode+k-n,kk,2),"u(1+this%mode+k-n,kk,2)"
                    i = i + 1
                  elseif ( aa_offeset == 0  .and. bb_offeset == 1) then
                    HYPRE_BUF(i) = -u(1+m+k-n,nn,2)
                    i = i + 1
                  elseif ( aa_offeset == 0  .and. bb_offeset == 0) then
                    HYPRE_BUF(i) = u(1+m+k-n,nn,1)
                    i = i + 1
                  endif
                endif
              endif
            endif
          elseif (bb == aa - 4*m - 1) then
            HYPRE_BUF(i) = 1.0/dr2 - 0.5/(j*dr2)
            i = i + 1
          elseif (bb == aa + 4*m + 1) then
            HYPRE_BUF(i) = 1.0/dr2 + 0.5/(j*dr2)
            i = i + 1  
          else               
            HYPRE_BUF(i) = 0.0
            i = i + 1
          endif
          jj = jj + 1
        enddo
        ii = ii + 1
      enddo
    enddo
! inner boundary
    i = 1
    ii = 1
    jj = 1
    if (lidproc == 0) then
      do aa = 0, 4*m
        rows(ii) = (nn - 1) *(4*m +1) + aa + 1
        do bb = 0, 12*m +3
          cols(jj) = (nn - 1) *(4*m +1) + aa + bb + 1
          if (aa == bb) then
!             m=-1
            if (aa == 2*m-1) then
              if(m > 1) then
                HYPRE_BUF(i) = u(1 + m,1,1) - u(m + 3,1,1) + 4.0/dr2
                i = i + 1
              else
                HYPRE_BUF(i) = u(1 + m,1,1) + 4.0/dr2
!                 HYPRE_BUF(i) = -4.0/dr2
                i = i + 1
              endif
            elseif(aa == 2*m -2) then
              if(m > 1) then
                HYPRE_BUF(i) = u(1+m,1,1) + u(m + 3,1,1) - 4.0/dr2
                i = i + 1
              else
                HYPRE_BUF(i) = u(1 + m,1,1) - 4.0/dr2
!                 HYPRE_BUF(i) = -4.0/dr2
                i = i + 1
              endif
!           m=1
!             elseif(aa == 2*this%mode+1) then
            elseif(aa == 2*m+ 1 .or. aa == 2*m+2) then
              HYPRE_BUF(i) = u(1 + m,1,1) - 4.0/dr2
!               HYPRE_BUF(i) = -4.0/dr2
              i = i + 1
!             elseif(aa == 2*this%mode+2) then
!               HYPRE_BUF(i) = u(1+this%mode,1,1) - 4.0/dr2
! !               HYPRE_BUF(i) = -4.0/dr2
!               i = i + 1
            else
              HYPRE_BUF(i) = 1.0/dr2
              i = i + 1
            endif
          else
            if( this%mode > 1) then
              if( bb == aa + 4*m + 1 ) then
                if (aa /= 2*m - 1 .and. aa /= 2*m - 2 .and. aa /= 2*m + 1 .and. aa /= 2*m + 2) then
                  HYPRE_BUF(i) = 0.0
                  i = i + 1
                elseif(aa == 2*m - 1) then
                  HYPRE_BUF(i) = -4.0/dr2
                  i = i + 1
                else
                  HYPRE_BUF(i) = 4.0/dr2
                  i = i + 1
                endif
              elseif( bb == aa - 1) then
                if(aa == 2*m - 1) then
                  HYPRE_BUF(i) = -u(m + 3,1,2)
                  i = i + 1
                else
                  HYPRE_BUF(i) = 0.0
                  i = i + 1
                endif
              elseif( bb == aa + 1) then
                if(aa == 2*m - 2) then
                  HYPRE_BUF(i) = -u(m + 3,1,2)
                  i = i + 1
                else
                  HYPRE_BUF(i) = 0.0
                  i = i + 1
                endif
              elseif( bb == aa - 3) then
                if(aa == 2*m + 1) then
                  HYPRE_BUF(i) = u(m + 3,1,1)
                  i = i + 1
                elseif(aa == 2*m + 2) then
                  HYPRE_BUF(i) = u(m + 3,1,2)
                  i = i + 1
                else
                  HYPRE_BUF(i) = 0.0
                  i = i + 1
                endif
              elseif( bb == aa - 2) then
                if(aa == 2*m + 1) then
                  HYPRE_BUF(i) = -u(m + 3,1,2)
                  i = i + 1
                elseif(aa == 2*m + 2) then
                  HYPRE_BUF(i) = u(m + 3,1,1)
                  i = i + 1
                else
                  HYPRE_BUF(i) = 0.0
                  i = i + 1
                endif
              else
                HYPRE_BUF(i) = 0.0
                i = i + 1
              endif
            elseif(m == 1) then
              if( bb == aa + 4*m + 1 ) then
                if (aa /= 2*this%mode - 1 .and. aa /= 2*this%mode - 2 .and. aa /= 2*this%mode + 1 .and. aa /= 2*this%mode + 2) then
                  HYPRE_BUF(i) = 0.0
                  i = i + 1
                elseif(aa == 2*this%mode - 1) then
                  HYPRE_BUF(i) = -4.0/dr2
!                   HYPRE_BUF(i) = 4.0/dr2
                  i = i + 1
                else
                  HYPRE_BUF(i) = 4.0/dr2
!                   HYPRE_BUF(i) = 0.0
                  i = i + 1
                endif
              else
                HYPRE_BUF(i) = 0.0
                i = i + 1
              endif
            else
              HYPRE_BUF(i) = 0.0
              i = i + 1
            endif
          endif
          jj = jj + 1
        enddo
        ii = ii + 1
      enddo
      i = (12*m +3)*(4*m + 1)+1
      ii = (4*m + 1)+1
      jj = (12*m +3)*(4*m + 1)+1
      do aa = 0, 4*m
        rows(ii + aa) = ii +aa
        if (aa /= 2*m - 1 .and. aa /= 2*m - 2 .and. aa /= 2*m + 1 .and. aa /= 2*m + 2) then
          cols(jj + aa*15) = aa + 1
          HYPRE_BUF(i + aa*15) = 0.0
        endif
      enddo 
    else
      i = 1
      ii = 1
      jj = 1
      j = real(noff) 
      kk = int(j)
      do aa = 0, 4*m
        if( aa - 2*m > 0 ) then
          k = (aa + 1)/2 - m
        else
          k = aa/2 - m
        endif
        rows(ii) = (nn - 1)*(4*m + 1)+aa +1
        do bb = 0, 12*this%mode + 3
          cols(jj) = (nn - 1)*(4*m + 1)+aa+bb +1
          if(  (bb >= 0) .and. (bb <= 4*m) ) then
            if( bb - 2*m > 0 ) then
              n = (bb + 1)/2 - m 
            else
              n = bb/2 - m 
            endif
            demode = abs(k-n)
            if (aa == bb) then
              HYPRE_BUF(i)   = (-2.0 - ((n-1)/j)**2)/dr2 + u(1+m,1,1)
              i = i + 1
            else
              if (demode > m) then
                HYPRE_BUF(i) = 0.0
                i = i + 1
              else
                if (demode == 0) then
                  HYPRE_BUF(i) = 0.0
                  i = i + 1
                else
                  if (aa - 2*m > 0) then
                    aa_offeset = aa + 1 - (k+m)*2
                  else
                    aa_offeset = aa - (k+m)*2
                  endif
                  if( bb - 2*m > 0 ) then
                    bb_offeset = bb + 1 - (n+m)*2
                  else
                    bb_offeset = bb - (n+m)*2
                  endif
                  if ( aa_offeset == 1  .and. bb_offeset == 1) then
                    HYPRE_BUF(i) = u(1+m+k-n,1,1)
                    i = i + 1
                  elseif ( aa_offeset == 1  .and. bb_offeset == 0) then
                    HYPRE_BUF(i) = u(1+m+k-n,1,2)
                    i = i + 1
                  elseif ( aa_offeset == 0  .and. bb_offeset == 1) then
                    HYPRE_BUF(i) = -u(1+m+k-n,1,2)
                    i = i + 1
                  elseif ( aa_offeset == 0  .and. bb_offeset == 0) then
                    HYPRE_BUF(i) = u(1+m+k-n,1,1)
                    i = i + 1
                  endif
                endif
              endif
            endif
          elseif (bb == aa - 4*m - 1) then
            HYPRE_BUF(i) = 1.0/dr2 - 0.5/(j*dr2)
            i = i + 1
          elseif (bb == aa + 4*m + 1) then
            HYPRE_BUF(i) = 1.0/dr2 + 0.5/(j*dr2)
            i = i + 1  
          else               
            HYPRE_BUF(i) = 0.0
            i = i + 1
          endif
          jj = jj + 1
        enddo
        ii = ii + 1
      enddo
    endif 

  case( p_fk_all_B_plus )

    ! set the first grid point of each partition
!     write(2,*) u(:,:,:),"u"
    j = real(noff)
    i = 1 + (12*m + 3)*(4*m + 1)
    ii = 1 + (4*m + 1)
    jj = 1 + (12*m + 3)*(4*m + 1)
     do nn = 2, nr
      j = j + 1.0
      kk = int(j)
      do aa = 0, 4*m
        if( aa - 2*m > 0 ) then
          k = (aa + 1)/2 - m
        else
          k = aa/2 - m
        endif
        rows(ii) = (nn - 1)*(4*m + 1)+aa +1
        do bb = - 4*m - 1, 8*this%mode + 1
          cols(jj) = (nn - 1)*(4*m + 1)+aa+bb +1
          if(  (bb >= 0) .and. (bb <= 4*m) ) then
            if( bb - 2*m > 0 ) then
              n = (bb + 1)/2 - m 
            else
              n = bb/2 - m 
            endif
            demode = abs(k-n)
            if (aa == bb) then
!               HYPRE_BUF(i) = (-2.0 - ((n-1)/j)**2)/dr2 + u(1+this%mode,kk,1)
              HYPRE_BUF(i) = (-2.0 - ((n + 1)/j)**2)/dr2 + u(1+m,nn,1)
              i = i + 1
            else
              if (demode > m) then
                HYPRE_BUF(i) = 0.0
                i = i + 1
              else
                if (demode == 0) then
                  HYPRE_BUF(i) = 0.0
                  i = i + 1
                else
                  if (aa - 2*m > 0) then
                    aa_offeset = aa + 1 - (k+m)*2
                  else
                    aa_offeset = aa - (k+m)*2
                  endif
                  if( bb - 2*m > 0 ) then
                    bb_offeset = bb + 1 - (n+m)*2
                  else
                    bb_offeset = bb - (n+m)*2
                  endif
                  if ( aa_offeset == 1  .and. bb_offeset == 1) then
                    HYPRE_BUF(i) = u(1+m+k-n,nn,1)
                    i = i + 1
                  elseif ( aa_offeset == 1  .and. bb_offeset == 0) then
                    HYPRE_BUF(i) = u(1+m+k-n,nn,2)
                    i = i + 1
                  elseif ( aa_offeset == 0  .and. bb_offeset == 1) then
                    HYPRE_BUF(i) = -u(1+m+k-n,nn,2)
                    i = i + 1
                  elseif ( aa_offeset == 0  .and. bb_offeset == 0) then
                    HYPRE_BUF(i) = u(1+m+k-n,nn,1)
                    i = i + 1
                  endif
                endif
              endif
            endif
          elseif (bb == aa -  4*m - 1) then
            HYPRE_BUF(i) = 1.0/dr2 - 0.5/(j*dr2)
            i = i + 1
          elseif (bb == aa +  4*m + 1) then
            HYPRE_BUF(i) = 1.0/dr2 + 0.5/(j*dr2)
            i = i + 1  
          else               
            HYPRE_BUF(i) = 0.0
            i = i + 1
          endif
          jj = jj + 1
        enddo
        ii = ii + 1
      enddo
    enddo
! inner boundary
    i = 1
    ii = 1
    jj = 1
    if (lidproc == 0) then
      do aa = 0, 4*m
        rows(ii) = (nn - 1) *(4*m +1) + aa + 1
        do bb = 0, 12*m + 3
          cols(jj) = (nn - 1) *(4*m +1) + aa + bb + 1
          if (aa == bb) then
            HYPRE_BUF(i)   = 1.0/dr2
            i = i + 1
          else
            HYPRE_BUF(i) = 0.0
            i = i + 1
          endif
        enddo
        ii = ii + 1
      enddo
      i = (12*m +3)*(4*m + 1)+1
      ii = (4*m + 1)+1
      jj = (12*m +3)*(4*m + 1)+1
      do aa = 0, 4*m
        rows(ii + aa) = ii +aa
        if (aa /= 2*m - 1 .and. aa /= 2*m - 2 .and. aa /= 2*m + 1 .and. aa /= 2*m + 2) then
          cols(jj + aa*15) = aa + 1
          HYPRE_BUF(i + aa*15) = 0.0
        endif
      enddo
    else
      i = 1
      ii = 1
      jj = 1
      j = real(noff) 
      kk = int(j)
      do aa = 0, 4*this%mode
        if( aa - 2*this%mode > 0 ) then
          k = (aa + 1)/2 - this%mode
        else
          k = aa/2 - this%mode
        endif
        rows(ii) = (nn - 1) *(4*m +1) + aa + 1
        do bb = 0, 12*this%mode + 3
          cols(jj) = (nn - 1) *(4*m +1) + aa + bb + 1
          if(  (bb >= 0) .and. (bb <= 4*this%mode) ) then
            if( bb - 2*this%mode > 0 ) then
              n = (bb + 1)/2 - this%mode 
            else
              n = bb/2 - this%mode 
            endif
            demode = abs(k-n)
            if (aa == bb) then
              HYPRE_BUF(i) = (-2.0 - ((n + 1)/j)**2)/dr2 + u(1+m,1,1)
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
                    aa_offeset = aa + 1 - (k+m)*2
                  else
                    aa_offeset = aa - (k+m)*2
                  endif
                  if( bb - 2*this%mode > 0 ) then
                    bb_offeset = bb + 1 - (n+m)*2
                  else
                    bb_offeset = bb - (n+m)*2
                  endif
                  if ( aa_offeset == 1  .and. bb_offeset == 1) then
                    HYPRE_BUF(i) = u(1+m+k-n,1,1)
                    i = i + 1
                  elseif ( aa_offeset == 1  .and. bb_offeset == 0) then
                    HYPRE_BUF(i) = u(1+m+k-n,1,2)
                    i = i + 1
                  elseif ( aa_offeset == 0  .and. bb_offeset == 1) then
                    HYPRE_BUF(i) = -u(1+m+k-n,1,2)
                    i = i + 1
                  elseif ( aa_offeset == 0  .and. bb_offeset == 0) then
                    HYPRE_BUF(i) = u(1+m+k-n,1,1)
                    i = i + 1
                  endif
                endif
              endif
            endif
          elseif (bb == aa - 4*m - 1) then
            HYPRE_BUF(i) = 1.0/dr2 - 0.5/(j*dr2)
            i = i + 1
          elseif (bb == aa + 4*m + 1) then
            HYPRE_BUF(i) = 1.0/dr2 + 0.5/(j*dr2)
            i = i + 1  
          else               
            HYPRE_BUF(i) = 0.0
            i = i + 1
          endif
          jj =jj +1
        enddo
        ii = ii +1
      enddo
    endif 

  end select

  ! set the upper boundary
  if ( lidproc == lnvp-1 ) then

    select case ( this%kind )

    case ( p_fk_all_B_minus, p_fk_all_B_plus )
      
      select case ( this%bnd )

      case ( p_bnd_zero )

        i = local_vol - (12*this%mode+3)*(4*this%mode + 1) + 1
        do aa = 0, 4*this%mode
          do bb = aa - 6*this%mode - 1, aa + 6*this%mode + 1
            if ( bb == aa + 4*this%mode + 1) then
              HYPRE_BUF(i) = 0.0
              i = i + 1
            else
              i = i + 1
            endif
          enddo
        enddo 
      
      case ( p_bnd_open )
        i = nr*(4*m +1)*(12*m +3) - (12*m +3)*(4*m + 1) + 1
        ii = 1+(4*m +1)*(nr - 1)
        jj = nr*(4*m +1)*(12*m +3) - (12*m +3)*(4*m + 1) + 1
        jmax = real(noff + nr)
        do aa = 0, 4*m
          rows(ii + aa) = ii + aa 
          do bb = -4*m - 1, 8*this%mode + 1
            cols(jj) = ii + aa + bb + 1
            if (aa == bb) then
              if( aa - 2*m > 0 ) then
                g = (aa + 1)/2 - m
!                 write(2,*) m,"m"
              else
                g = aa/2 - m
!                 write(2,*) m,"m"
              endif
!               g = abs(g)
              HYPRE_BUF(i) = HYPRE_BUF(i) + (1.0-(g + 1)/jmax) * HYPRE_BUF(i + 4*m + 1)
              i = i + 1
            elseif ( bb == aa + 4*m + 1) then
!             elseif ( bb > aa ) then
              HYPRE_BUF(i) = 0.0
              i = i + 1
            else
              i = i + 1
            endif
            jj = jj + 1
          enddo
        enddo 
      end select

    end select

  endif
        
  call HYPRE_IJMatrixSetValues(this%A, local_vol, ncols, rows, cols, HYPRE_BUF, ierr)      
  call HYPRE_IJMatrixAssemble(this%A)
  call HYPRE_IJMatrixGetObject(this%A, this%A_internal)  ! 获取内部 ParCSR 矩阵对象

  deallocate(rows, cols, ncols)

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine set_IJ_matrix

subroutine set_hypre_src_coef(this, qe_re, qn1_re, qe_im, qn1_im, qbm )

  implicit none

  class( field_solver_all_modesv1 ), intent(inout) :: this
  type(ufield), intent(inout), dimension(:), optional, pointer :: qe_re
  type(ufield), intent(inout), dimension(:), optional, pointer :: qn1_re
  type(ufield), intent(inout), dimension(:), optional, pointer :: qe_im
  type(ufield), intent(inout), dimension(:), optional, pointer :: qn1_im
  real, intent(inout), optional :: qbm

  integer :: ierr, a, im, k, l, kk
  real :: jj, j
  real, dimension(:,:), pointer :: f1_re => null(), f1_im => null()
  real, dimension(:,:), pointer :: f2_re => null(), f2_im => null()
  real, dimension(:,:), pointer :: f3_re => null(), f3_im => null()
  character(len=32), save :: sname = "set_hypre_src_coef"

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  a = 1
!   j = real(nofff) + 1.0
!   write(2,*) j,'j'
  do k = 1, nr
!     kk = int(j)
    do l = 0, 2*this%mode
      im = l - this%mode
!         write(2,*) im,"im"
      if (im > 0 ) then
        f1_re => qe_re(im)%get_f1()
        f1_im => qe_im(im)%get_f1()
        if ( present(qe_im) .and. present(qn1_im) ) then
!             write(2,*) "present(qn_im)"
          f2_re => qn1_re(im)%get_f1()
          f2_im => qn1_im(im)%get_f1()
          f1_re(1,k) = f1_re (1,k)- f2_re(1,k)
          f1_im(1,k) = f1_im (1,k)- f2_im(1,k)
!           write(2,*) f1_im(1,k),'f1_im(1,k) im>0'
        endif
        src_sol(a)  = -qbm*f1_re(1,k)
        src_sol(a + 1)  = -qbm*f1_im(1,k)
!         write(2,*) src_sol(a),'f1_re(1,k) im>0'
!         write(2,*) src_sol(a + 1),'f1_im(1,k) im>0'
        a=a + 2
!         write(2,*) a,'a'
      elseif (im == 0) then
        f1_re => qe_re(0)%get_f1()
!         write(2,*) qe_re(0)%get_f1(),'size(qe_re(0)%get_f1())'
        if ( present(qe_re) .and. present(qn1_re) ) then
!             write(2,*) "present(qn_re)"
          f2_re => qn1_re(0)%get_f1()
!             write(2,*) f2_re(1,k),"f2_re(1,k) im=0"
          f1_re(1,k) = f1_re(1,k) - f2_re(1,k)
!             write(2,*) f1_re(1,kk),"f1_re(1,kk) im=0"
!           write(2,*) f1_im(1,k),'f1_im(1,k) im=0'
        endif
        src_sol(a)  = -qbm*f1_re(1,k)
!         write(2,*) src_sol(a),"f1_re(1,k) im=0"
        a = a + 1
!         write(2,*) a,'a'
      else
        im = abs(im)
        f1_re => qe_re(im)%get_f1()
        f1_im => qe_im(im)%get_f1()
        if ( present(qe_im) .and. present(qn1_im) ) then
!             write(2,*) "present(qn_im)"
          f2_re => qn1_re(im)%get_f1()
          f2_im => qn1_im(im)%get_f1()
          f1_re(1,k)= f1_re(1,k) - f2_re(1,k)
          f1_im(1,k) = f1_im(1,k) - f2_im(1,k)
!           write(2,*) f1_im(1,k),'f1_im(1,k) im<0'
        endif
        src_sol(a) = -qbm*f1_re(1,k)
        src_sol(a + 1)  = qbm*f1_im(1,k)
!         write(2,*) src_sol(a),'f1_re(1,k) im<0'
!         write(2,*) src_sol(a + 1),'f1_im(1,k) im<0'
        a=a + 2 
!         write(2,*) a,"a"         
      endif
    enddo
!     j = j + 1.0
  enddo
!     write(2,*) a,"a"
!     write(2,*) 'this%solver_coef%solve||HYPRE_StructVectorSetBoxValues 1'    
  call HYPRE_IJVectorSetValues(this%b, this%num_rows, this%indices, src_sol)
  write(2,*) size(src_sol),'size(src_sol)'
!   write(2,*) src_sol,'src_sol'

end subroutine set_hypre_src_coef

subroutine set_hypre_src_b(this, src)

  implicit none

  class( field_solver_all_modesv1 ), intent(inout) :: this
  real, intent(inout), dimension(:,:,:), optional, pointer :: src

  integer :: ierr, a, im, k, l, kk, j
  character(len=32), save :: sname = "set_hypre_src_b"

  a = 1
  j = real(nofff) + 1.0
  write(2,*) j,'j'
  do k = 1, nr
!     kk = int(j)
    do l = 0, 2*this%mode
        im = l - this%mode
!         if (im > 0 ) then
!           src_sol(a)  = src(1+l,k,1)
! !             src_sol(a)  = 0.0
!           src_sol(a+1)  = src(1+l,k,2)
! !             src_sol(a+1)  = 0.0
!           a=a+2
        if (im == 0) then
          src_sol(a)  = src(1+l,k,1)
!             write(2,*) src_sol(a),"src"
          a = a+1
        else
          src_sol(a)  = src(1+l,k,1)
!             src_sol(a)  = 0.0
          src_sol(a+1)  = src(1+l,k,2)
!             src_sol(a+1)  = 0.0
          a=a+2          
      endif
    enddo
    j = j + 1.0 
  enddo
  
  call HYPRE_IJVectorSetValues(this%b, this%num_rows, this%indices, src_sol)
  write(2,*) size(src_sol),'size(src_sol)'
  write(2,*) src_sol,'src_sol'

end subroutine set_hypre_src_b

subroutine get_hypre_src_coef(this, src)

  implicit none

  class( field_solver_all_modesv1 ), intent(inout) :: this
  real, intent(inout), dimension(:,:,:), optional, pointer :: src

  integer :: ierr, a, im, k, l
  character(len=32), save :: sname = "get_hypre_src_coef"

  call HYPRE_IJVectorGetValues(this%x, this%num_rows, this%indices, src_sol)
!   call HYPRE_IJVectorPrint( this%x,"struct_x_output.out")
!   call HYPRE_StructVectorPrint( "struct_x_output.out", this%x)
  write(2,*) size(src_sol),'size(src_sol)'
  write(2,*) src_sol(:),'coef src_sol'
    a = 1
    do k = 1, nr
      do l = 0, 2*this%mode
        im = l - this%mode
!         write(2,*) im,"im get"
!         write(2,*) a,"a get"
          if (im == 0 ) then
!             write(2,*) src_sol(a),"src_sol(a) im=0 0"
            src(1+l,k,1) = src_sol(a)
!             write(2,*) src_sol(a),"src_sol(a) im=0 1"
            a=a+1
          else
            src(1+l,k,1) = src_sol(a)
            src(1+l,k,2) = src_sol(a + 1)
            a=a+2          
        endif
      enddo 
    enddo
  write(2,*) src(:,:,:),'coef src'

end subroutine get_hypre_src_coef

subroutine get_hypre_src_b(this, src)

  implicit none

  class( field_solver_all_modesv1 ), intent(inout) :: this
  real, intent(inout), dimension(:,:,:), optional, pointer :: src

  integer :: ierr, a, im, k, l
  character(len=32), save :: sname = "get_hypre_src_b"

  call HYPRE_IJVectorGetValues(this%x, this%num_rows, this%indices, src_sol)
!   call HYPRE_IJVectorPrint( this%x,"vector_x.out")
  write(2,*) size(src_sol),'size(src_sol)'
!     write(2,*) src_sol,"src_sol 3" 
  a = 1
  do k = 1, nr
    do l = 0, 2*this%mode
      im = l - this%mode
      if (im == 0 ) then
        src(1+l,k,1) = src_sol(a)
        a=a+1
      else
        src(1+l,k,1) = src_sol(a)
        src(1+l,k,2) = src_sol(a + 1)
        a=a+2          
      endif
    enddo 
  enddo

end subroutine get_hypre_src_b

end module field_solver_all_modesv1_class
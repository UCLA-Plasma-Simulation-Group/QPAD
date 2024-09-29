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

character(len=32), parameter :: cls_name = "field_solver_all_modes"
integer, parameter :: cls_level = 4

real, dimension(:), pointer, save :: HYPRE_BUF => null()
real, dimension(:), pointer :: src_sol => null()

type :: field_solver_all_modesv1 ! class for HYPRE solver

  ! HYPRE parameters
  integer, dimension(:), pointer :: offsets => null()
!   integer, dimension(:), pointer :: stencil_idx => null()
  integer :: num_stencil
  integer :: stype 
  integer :: kind
  integer :: mode
  integer :: bnd

!   integer(HYPRE_TYPE) :: A, b, x, grid, stencil, solver, precond
  integer(HYPRE_TYPE) :: solver, precond
  integer(HYPRE_TYPE) :: A_parcsr, b_par, x_par
  integer(HYPRE_TYPE) :: A_ij, b_ij, x_ij
  integer :: iupper, ilower

  contains

  procedure :: new => init_field_solver
  procedure :: del => end_field_solver 
  procedure :: solve => solve_equation
  procedure, private :: set_struct_solver
!   procedure, private :: set_struct_grid
!   procedure, private :: set_struct_stencil
  procedure, private :: set_struct_matrix 
 
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
  write(2,*) comm,'comm_loc'
  this%ilower = opts%get_noff(1) 
  this%iupper = opts%get_noff(1) + opts%get_ndp(1)*(4*this%mode+1) - 1

!   call this%set_struct_grid( opts )
!   call this%set_struct_stencil()
!   call this%set_struct_matrix( opts, dr )

!   call HYPRE_StructVectorCreate( comm, this%grid, this%b, ierr )
!   call HYPRE_StructVectorInitialize( this%b, ierr )
  call HYPRE_IJVectorCreate(comm, this%ilower, this%iupper, this%b_ij, ierr)
  call HYPRE_IJVectorSetObjectType(this%b_ij, HYPRE_PARCSR, ierr)
  call HYPRE_IJVectorInitialize(this%b_ij, ierr)
  call HYPRE_IJVectorGetObject(this%b_ij, this%b_par, ierr)

!   call HYPRE_StructVectorCreate( comm, this%grid, this%x, ierr )
!   call HYPRE_StructVectorInitialize( this%x, ierr )
  call HYPRE_IJVectorCreate(comm, this%ilower, this%iupper, this%x_ij, ierr)
  call HYPRE_IJVectorSetObjectType(this%x_ij, HYPRE_PARCSR, ierr)
  call HYPRE_IJVectorInitialize(this%x_ij, ierr)
  call HYPRE_IJVectorGetObject(this%x_ij, this%x_par, ierr) 

!   call this%set_struct_solver()

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine init_field_solver

subroutine end_field_solver( this )

  implicit none

  class( field_solver_all_modesv1 ), intent(inout) :: this

  integer :: ierr
  character(len=32), save :: sname = "end_field_solver"

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  if ( associated(HYPRE_BUF) ) deallocate( HYPRE_BUF )
  if ( associated(src_sol) ) deallocate( src_sol )

!   call HYPRE_StructGridDestroy( this%grid, ierr )
!   call HYPRE_StructStencilDestroy( this%stencil, ierr )
!   call HYPRE_StructVectorDestroy( this%b, ierr )
!   call HYPRE_StructVectorDestroy( this%x, ierr )
!   call HYPRE_StructMatrixDestroy( this%A, ierr )
  call HYPRE_ParCSRMatrixDestroy( this%A_parcsr, ierr )
  call HYPRE_ParVectorDestroy( this%b_par, ierr )
  call HYPRE_ParVectorDestroy( this%x_par, ierr )
  call HYPRE_IJVectorDestroy(this%A_ij, ierr)
  call HYPRE_IJVectorDestroy(this%b_ij, ierr)
  call HYPRE_IJVectorDestroy(this%x_ij, ierr)

  call HYPRE_ParCSRGMRESDestroy( this%solver, ierr )

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine end_field_solver

subroutine set_struct_solver( this )

  implicit none

  class( field_solver_all_modesv1 ), intent(inout) :: this

  ! External functions
!   external HYPRE_BoomerAMGSetup, HYPRE_BoomerAMGSolve

  integer :: ierr, comm
  character(len=32), save :: sname = "set_struct_solver"

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  comm = comm_loc()

  write(2,*) comm,'comm_loc'
  write(2,*) 'solver Initializing 2'
  ! 使用 BoomerAMG 作为预处理器
!   call HYPRE_BoomerAMGCreate(this%precond, ierr)
!   call HYPRE_BoomerAMGSetMaxIter(this%precond, 1, ierr)
!   call HYPRE_BoomerAMGSetTol(this%precond, 1e-6, ierr)

  call HYPRE_ParCSRGMRESCreate(comm, this%solver, ierr)
  call HYPRE_ParCSRGMRESSetMaxIter(this%solver, 1000, ierr)
  call HYPRE_ParCSRGMRESSetTol(this%solver, 1e-6, ierr)

!   call HYPRE_ParCSRGMRESSetPrecond(this%solver, HYPRE_BoomerAMGSolve, HYPRE_BoomerAMGSetup, this%precond, ierr)
  write(2,*) 'solver Initializing 3'
!   write(2,*) this%A,'A'
!   write(2,*) this%b,'b'
!   write(2,*) this%x,'x'
  call HYPRE_ParCSRGMRESSetup( this%solver, this%A_parcsr, this%b_par, this%x_par, ierr )
  write(2,*) 'solver Initializing 4'

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine set_struct_solver

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

  integer :: ierr, a, im, k, l, kk, i
  real :: jj, j
  integer :: size_sol
  integer :: first_local_index, last_local_index
  integer :: local_size
  real(8), allocatable :: local_values(:)
  integer, allocatable :: indices(:)
  integer, dimension(:), pointer :: row_indices => null()
  integer, dimension(:), pointer :: col_indices => null()
  real, dimension(:), pointer :: values => null()
  real, dimension(:,:), pointer :: f1_re => null(), f1_im => null()
  real, dimension(:,:), pointer :: f2_re => null(), f2_im => null()
  real, dimension(:,:), pointer :: f3_re => null(), f3_im => null()
  character(len=32), save :: sname = "solve_equation"

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  size_sol = nr * (4 * this%mode + 1)  

  if (.not. associated( src_sol)) then
    allocate( src_sol(size_sol) )
  end if
  if (.not. associated( row_indices)) then
    allocate( row_indices(size_sol) )
  end if
  if (.not. associated( col_indices)) then
    allocate( col_indices(size_sol) )
  end if
  if (.not. associated( values)) then
    allocate( values(size_sol) )
  end if

  write(2,*) this%kind, 'kind' 
  select case ( this%kind )
  case ( p_fk_coef )
    call set_struct_matrix(this, psi_re = psi_re, psi_im = psi_im, qbm=qbm)
!     call this%set_struct_solver()
    write(2,*) 'this%solver_coef%solve||HYPRE_StructVectorSetBoxValues 0' 
    a = 1
    j = real(nofff)
    do k = 1, nr
      kk = int(j)
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
            f1_re(1,kk) = f1_re (1,kk)- f2_re(1,kk)
            f1_im(1,kk) = f1_im (1,kk)- f2_im(1,kk)
          endif
          src_sol(a)  = -qbm*f1_re(1,kk)
          src_sol(a + 1)  = -qbm*f1_im(1,kk)
          a=a + 2
        elseif (im == 0) then
          f1_re => qe_re(0)%get_f1()
          if ( present(qe_re) .and. present(qn1_re) ) then
!             write(2,*) "present(qn_re)"
            f2_re => qn1_re(0)%get_f1()
!             write(2,*) f2_re(1,kk),"f2_re(1,kk) im=0"
            f1_re(1,kk) = f1_re(1,kk) - f2_re(1,kk)
!             write(2,*) f1_re(1,kk),"f1_re(1,kk) im=0"
          endif
!           write(2,*) f1_re(1,k),"f1_re(1,k) im=0"
          src_sol(a)  = -qbm*f1_re(1,kk)
!           src_sol(a)  = 0.0
          a = a + 1
        else
          im = abs(im)
          f1_re => qe_re(im)%get_f1()
          f1_im => qe_im(im)%get_f1()
          if ( present(qe_im) .and. present(qn1_im) ) then
!             write(2,*) "present(qn_im)"
            f2_re => qn1_re(im)%get_f1()
            f2_im => qn1_im(im)%get_f1()
            f1_re(1,kk)= f1_re(1,kk) - f2_re(1,kk)
            f1_im(1,kk) = f1_im(1,kk) - f2_im(1,kk)
          endif
          src_sol(a) = -qbm*f1_re(1,kk)
          src_sol(a + 1)  = qbm*f1_im(1,kk)
!           write(2,*) qbm*f1_im(1,k),"src_sol(a+1) 0  im<0"
!           write(2,*) src_sol(a + 1),"src_sol(a+1) 1  im<0"
          a=a + 2 
!           write(2,*) a,"a"         
        endif
      enddo
      j = j + 1.0
    enddo
!     write(2,*) a,"a"
    write(2,*) 'this%solver_coef%solve||HYPRE_StructVectorSetBoxValues 1'     
!     call HYPRE_StructVectorSetBoxValues( this%b, this%ilower, this%iupper, src_sol, ierr )
!     write(2,*) 'this%solver_coef%solve||HYPRE_StructVectorAssembl'   
!     call HYPRE_StructVectorAssemble( this%b, ierr )
    a = 1
    do i = 1, size_sol
       row_indices(i) = i - 1         ! Adjust as necessary, e.g., row i
       col_indices(i) = 0             ! Adjust as necessary, e.g., column 0
       values(i) = src_sol(a)         ! Set values from src_sol
       a = a + 1
    end do

    ! Set the values in b_ij
    call HYPRE_IJMatrixSetValues(this%b_ij, size_sol, row_indices, col_indices, values, ierr)

    ! Finalize the matrix
    call HYPRE_IJMatrixAssemble(this%b_ij, ierr)
    call this%set_struct_solver()
    write(2,*) 'this%solver_coef%solve||HYPRE_ParCSRGMRESSolve 0'  
    ! 获取 ParCSR 格式的矩阵和向量
    call HYPRE_ParCSRGMRESSolve( this%solver, this%A_parcsr, this%b_par, this%x_par, ierr )
    write(2,*) 'this%solver_coef%solve||HYPRE_ParCSRGMRESSolve 1' 
    ! Retrieve the solution from the x_par ParCSR vector
    call HYPRE_IJVectorGetLocalRange(this%x_ij, first_local_index, last_local_index, ierr)
    write(2,*) 'this%solver_coef%solve||HYPRE_ParCSRGMRESSolve 2'
    ! 提取本地解

    local_size = last_local_index - first_local_index + 1

    allocate(local_values(local_size))

    ! 使用 HYPRE_IJVectorGetValues 来获取解

    allocate(indices(local_size))
    do i = 1, local_size
      indices(i) = first_local_index + (i - 1)
    end do

    ! 获取局部向量的值
    call HYPRE_IJVectorGetValues(this%x_ij, local_size, indices, local_values, ierr)
    ! Map solution from src_sol to 'src'
    a = 1
    do k = 1, nr
      do l = 0, 2*this%mode
        im = l - this%mode
        if (im == 0) then
            src(1+l,k,1) = local_values(a)
            a = a + 1
        else
            src(1+l,k,1) = local_values(a)
            src(1+l,k,2) = local_values(a + 1)
            a = a + 2          
        endif
      end do
    end do
  case ( p_fk_all_B_minus, p_fk_all_B_plus )
    call set_struct_matrix(this, u=u)
!     call this%set_struct_solver()
    write(2,*) HYPRE_BUF(:),"HYPRE_BUF"
    write(2,*) src,"src"  
    a = 1
    j = real(nofff) + 1.0
    do k = 0, nr
      kk = int(j)
      do l = 0, 2*this%mode
          im = l - this%mode
          if (im > 0 ) then
!             src_sol(a)  = src(1+l,kk,1)
            src_sol(a)  = 0.0
!             src_sol(a+1)  = src(1+l,kk,2)
            src_sol(a+1)  = 0.0
            a=a+2
          elseif (im == 0) then
            src_sol(a)  = src(1+l,kk,1)
!             write(2,*) src_sol(a),"src"
            a = a+1
          else
!             src_sol(a)  = src(1+l,kk,1)
            src_sol(a)  = 0.0
!             src_sol(a+1)  = -src(1+l,kk,2)
            src_sol(a+1)  = 0.0
            a=a+2          
        endif
      enddo
      j = j + 1.0 
    enddo
!     write(2,*) src_sol,"src_sol 1"
!     call HYPRE_StructVectorSetBoxValues( this%b, this%ilower, this%iupper, src_sol, ierr )
    write(2,*) src_sol,"src_sol 1"
!     call HYPRE_StructVectorAssemble( this%b, ierr )
    a = 1
    do i = 1, size_sol
       row_indices(i) = i - 1         ! Adjust as necessary, e.g., row i
       col_indices(i) = 0             ! Adjust as necessary, e.g., column 0
       values(i) = src_sol(a)         ! Set values from src_sol
       a = a + 1
    end do

    ! Set the values in b_ij
    call HYPRE_IJMatrixSetValues(this%b_ij, size_sol, row_indices, col_indices, values, ierr)

    ! Finalize the matrix
    call HYPRE_IJMatrixAssemble(this%b_ij, ierr)
    call this%set_struct_solver()
      ! 转换结构化矩阵和向量为 ParCSR 格式
    ! 获取 ParCSR 格式的矩阵和向量
    call HYPRE_IJVectorGetObject(this%b_ij, this%b_par, ierr)
    call HYPRE_IJVectorGetObject(this%x_ij, this%x_par, ierr)
    call HYPRE_ParCSRGMRESSolve( this%solver, this%A_parcsr, this%b_par, this%x_par, ierr )
    ! Retrieve the solution from the x_par ParCSR vector
    call HYPRE_IJVectorGetLocalRange(this%x_ij, first_local_index, last_local_index, ierr)
    ! 提取本地解

    local_size = last_local_index - first_local_index + 1

    allocate(local_values(local_size))

    ! 使用 HYPRE_IJVectorGetValues 来获取解

    allocate(indices(local_size))
    do i = 1, local_size
      indices(i) = first_local_index + (i - 1)
    end do

    ! 获取局部向量的值
    call HYPRE_IJVectorGetValues(this%x_ij, local_size, indices, local_values, ierr)

    ! Map solution from src_sol to 'src'
    a = 1
    do k = 1, nr
      do l = 0, 2*this%mode
        im = l - this%mode
        if (im == 0) then
            src(1+l,k,1) = local_values(a)
            a = a + 1
        else
            src(1+l,k,1) = local_values(a)
            src(1+l,k,2) = local_values(a + 1)
            a = a + 2          
        endif
      end do
    end do
    write(2,*) src,"src end" 
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

! subroutine set_struct_grid( this, opts )

!   implicit none

!   class( field_solver_all_modesv1 ), intent(inout) :: this
!   type( options ), intent(in) :: opts

!   integer :: comm, ierr
!   character(len=32), save :: sname = "set_struct_grid"

!   call write_dbg( cls_name, sname, cls_level, 'starts' )

!   comm = comm_loc()
!   this%ilower = opts%get_noff(1) + 1
!   this%iupper = opts%get_noff(1) + opts%get_ndp(1)*(4*this%mode+1) 

!   call HYPRE_StructGridCreate( comm, 1, this%grid, ierr )
!   call HYPRE_StructGridSetExtents( this%grid, this%ilower, this%iupper, ierr )
!   call HYPRE_StructGridAssemble( this%grid, ierr )

!   call write_dbg( cls_name, sname, cls_level, 'ends' )

! end subroutine set_struct_grid

! subroutine set_struct_stencil( this )

!   implicit none

!   class( field_solver_all_modesv1 ), intent(inout) :: this

!   integer :: i, j, k, l, ierr
!   character(len=32), save :: sname = "set_struct_stencil"

!   call write_dbg( cls_name, sname, cls_level, 'starts' )

!   select case(this%kind)

!   case(p_fk_coef)
    
!     this%num_stencil = 4*this%mode + 1
!     if ( .not. associated( this%offsets ) ) then
!       allocate( this%offsets(this%num_stencil ) )
!     endif

!     if ( .not. associated( this%stencil_idx ) ) then
!       allocate( this%stencil_idx( this%num_stencil ) )
!       do i = 1, this%num_stencil
!         this%stencil_idx(i) = i-1
!       enddo
!     endif

!     k = 1
!     do j = 0, 4*this%mode
!       this%offsets(k) = j-2*this%mode
!       k=k+1
!     enddo

!     call HYPRE_StructStencilCreate( 1, this%num_stencil, this%stencil, ierr )
!     do i = 1, this%num_stencil
!       call HYPRE_StructStencilSetElement( this%stencil, i-1, this%offsets(i), ierr )
!     enddo
  
!   case( p_fk_all_B_minus, p_fk_all_B_plus )

!     this%num_stencil = (4*this%mode + 1)*3
!     if ( .not. associated( this%offsets ) ) then
!       allocate( this%offsets( this%num_stencil ) )
!     endif

!     if ( .not. associated( this%stencil_idx ) ) then
!       allocate( this%stencil_idx( this%num_stencil ) )
!       do i = 1, this%num_stencil
!         this%stencil_idx(i) = i-1
!       enddo
!     endif

!     k = 1
!     do j = 0, (4*this%mode+1)*3-1
!       this%offsets(k) = j-(6*this%mode+1)
!       k=k+1
!     enddo

!     call HYPRE_StructStencilCreate( 1, this%num_stencil, this%stencil, ierr )
!     do i = 1, this%num_stencil
!       call HYPRE_StructStencilSetElement( this%stencil, i-1, this%offsets(i), ierr )
!     enddo

!   end select

!   call write_dbg( cls_name, sname, cls_level, 'ends' )

! end subroutine set_struct_stencil

subroutine set_struct_matrix( this, psi_re, psi_im, u, qbm)

  implicit none

  class( field_solver_all_modesv1 ), intent(inout) :: this
  type(ufield), intent(inout), dimension(:), optional, pointer :: psi_re
  type(ufield), intent(inout), dimension(:), optional, pointer :: psi_im
  real, intent(inout), dimension(:,:,:), optional :: u
  real, intent(inout), optional :: qbm

  integer :: i, ierr, local_vol, noff, m, aa, bb, nn, kk, mm, g, ii
  integer :: k, n, d, bb_offeset, aa_offeset,l, demode, Anumber
  integer :: comm, lidproc, lnvp
  integer, dimension(:), pointer :: row_indices => null()
  integer, dimension(:), pointer :: col_indices => null()
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
  write(2,*) comm,'comm_loc'
  write(2,*) lidproc,'id_proc_loc()'
  write(2,*) lnvp,'num_procs_loc()'
  select case ( this%kind )

  case ( p_fk_coef )
    this%num_stencil = 4*this%mode + 1
!     do ii = 1, this%num_stencil
!       this%stencil_idx(ii) = ii  ! 计算列索引
!     end do
  case( p_fk_all_B_plus, p_fk_all_B_minus )
    this%num_stencil = (4*this%mode + 1)*3
!     do ii = 1, this%num_stencil
!       this%stencil_idx(ii) = ii  ! 计算列索引
!     end do
  end select

  local_vol = (this%iupper - this%ilower + 1)* this%num_stencil
  write(2,*) 'this%solver_coef%solve||HYPRE_IJMatrixInitialize 3'

  if ( .not. associated( HYPRE_BUF ) ) then
    allocate( HYPRE_BUF( local_vol ) )
  elseif ( size(HYPRE_BUF) < local_vol ) then
    deallocate( HYPRE_BUF )
    allocate( HYPRE_BUF( local_vol ) )
  endif
  allocate( row_indices( local_vol ) )
  allocate( col_indices( local_vol ) )

  write(2,*) 'this%solver_coef%solve||HYPRE_IJMatrixInitialize 1'
!   call HYPRE_StructMatrixCreate( comm, this%grid, this%stencil, this%A, ierr )
!   call HYPRE_StructMatrixInitialize( this%A, ierr )

  call HYPRE_IJMatrixCreate(comm, this%ilower, this%iupper, this%ilower, this%iupper, this%A_ij, ierr)
  write(2,*) 'this%solver_coef%solve||HYPRE_IJMatrixInitialize 2'
  call HYPRE_IJMatrixSetObjectType(this%A_ij, HYPRE_PARCSR, ierr)
  write(2,*) 'this%solver_coef%solve||HYPRE_IJMatrixInitialize 3'
  call HYPRE_IJMatrixInitialize(this%A_ij, ierr)
  write(2,*) 'this%solver_coef%solve||HYPRE_IJMatrixInitialize 4'
  call HYPRE_IJMatrixGetObject(this%A_ij, this%A_parcsr, ierr)
  ! set the matrix element and lower boundary
  select case ( this%kind )

  case ( p_fk_coef )

    write(2,*) 'this%solver_coef%solve||set_struct_matrix' 
    j = real(noff)
    write(2,*) j,"j"
    write(2,*) qbm,"qbm"
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
!               write(2,*) aa,"aa"
!               write(2,*) n,"n"
              f1_re => psi_re(0)%get_f1()
!               write(2,*) f1_re(1,kk),"f1_re(1,kk)"
!               write(2,*) kk,"kk"
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
    write(2,*) 'this%solver_coef%solve||set_struct_matrix 1111' 
!     call HYPRE_StructMatrixSetBoxValues( this%A, this%ilower, this%iupper, this%num_stencil, &
!       this%stencil_idx, HYPRE_BUF, ierr )

!       ! 遍历每一行（网格点）
!     do i = this%ilower, this%iupper
!      ! 行索引
!      row_indices(1) = i

!      ! 对于每个行索引，设置相应的列索引 (根据 stencil 结构)
!       do ii = 1, this%num_stencil
!         col_indices(ii) = i + ii - 1 ! 计算列索引，邻接点根据 stencil 偏移
!         ! 插入非零值到 IJMatrix 中
!         write(2,*) 'this%solver_coef%solve||set_struct_matrix iiii'
!         call HYPRE_IJMatrixSetValues(this%A_ij, 1, 1, row_indices, col_indices(ii), HYPRE_BUF(ii), ierr)
!       end do
!     end do
!     write(2,*) 'this%solver_coef%solve||set_struct_matrix 222'
!     ! 完成矩阵的值插入后，进行矩阵组装
!     call HYPRE_IJMatrixAssemble(this%A_ij, ierr)
!     write(2,*) 'this%solver_coef%solve||set_struct_matrix 333'
!     ! 将 IJMatrix 转换为 ParCSRMatrix
!     call HYPRE_IJMatrixGetObject(this%A_ij, this%A_parcsr, ierr)

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
              HYPRE_BUF(i)   = (-2.0 - ((n-1)/j)**2)/dr2 + u(1+this%mode,kk,1)
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
                    HYPRE_BUF(i) = u(1+this%mode+k-n,kk,1)
!                     write(2,*) u(1+this%mode+k-n,kk,1),"u(1+this%mode+k-n,kk,1)"
                    i = i + 1
                  elseif ( aa_offeset == 1  .and. bb_offeset == 0) then
                    HYPRE_BUF(i) = u(1+this%mode+k-n,kk,2)
!                     write(2,*) u(1+this%mode+k-n,kk,2),"u(1+this%mode+k-n,kk,2)"
                    i = i + 1
                  elseif ( aa_offeset == 0  .and. bb_offeset == 1) then
                    HYPRE_BUF(i) = -u(1+this%mode+k-n,kk,2)
                    i = i + 1
                  elseif ( aa_offeset == 0  .and. bb_offeset == 0) then
                    HYPRE_BUF(i) = u(1+this%mode+k-n,kk,1)
                    i = i + 1
                  endif
                endif
              endif
            endif
          elseif (bb == aa - 4*this%mode - 1) then
            HYPRE_BUF(i) = 1.0/dr2 - 0.5/(j*dr2)
            i = i + 1
          elseif (bb == aa + 4*this%mode + 1) then
            HYPRE_BUF(i) = 1.0/dr2 + 0.5/(j*dr2)
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
              if(this%mode > 1) then
                HYPRE_BUF(i) = u(1 + this%mode,0,1) - u(this%mode + 3,0,1) + 4.0/dr2
                i = i + 1
              else
                HYPRE_BUF(i) = u(1 + this%mode,0,1) + 4.0/dr2
!                 HYPRE_BUF(i) = -4.0/dr2
                i = i + 1
              endif
            elseif(aa == 2*this%mode-2) then
              if(this%mode > 1) then
                HYPRE_BUF(i) = u(1+this%mode,1,1) + u(this%mode + 3,1,1) - 4.0/dr2
                i = i + 1
              else
                HYPRE_BUF(i) = u(1 + this%mode,0,1) - 4.0/dr2
!                 HYPRE_BUF(i) = -4.0/dr2
                i = i + 1
              endif
!           m=1
!             elseif(aa == 2*this%mode+1) then
            elseif(aa == 2*this%mode+1 .or. aa == 2*this%mode+2) then
              HYPRE_BUF(i) = u(1 + this%mode,0,1) - 4.0/dr2
!               HYPRE_BUF(i) = -4.0/dr2
              i = i + 1
!             elseif(aa == 2*this%mode+2) then
!               HYPRE_BUF(i) = u(1+this%mode,1,1) - 4.0/dr2
! !               HYPRE_BUF(i) = -4.0/dr2
!               i = i + 1
            else
              HYPRE_BUF(i) = 1.0
              i = i + 1
            endif
          else
            if( this%mode > 1) then
              if( bb == aa + 4*this%mode + 1 ) then
                if (aa /= 2*this%mode - 1 .and. aa /= 2*this%mode - 2 .and. aa /= 2*this%mode + 1 .and. aa /= 2*this%mode + 2) then
                  HYPRE_BUF(i) = 0.0
                  i = i + 1
                elseif(aa == 2*this%mode - 1) then
                  HYPRE_BUF(i) = -4.0/dr2
                  i = i + 1
                else
                  HYPRE_BUF(i) = 4.0/dr2
                  i = i + 1
                endif
              elseif( bb == aa - 1) then
                if(aa == 2*this%mode - 1) then
                  HYPRE_BUF(i) = -u(this%mode + 3,0,2)
                  i = i + 1
                else
                  HYPRE_BUF(i) = 0.0
                  i = i + 1
                endif
              elseif( bb == aa + 1) then
                if(aa == 2*this%mode - 2) then
                  HYPRE_BUF(i) = -u(this%mode + 3,0,2)
                  i = i + 1
                else
                  HYPRE_BUF(i) = 0.0
                  i = i + 1
                endif
              elseif( bb == aa - 3) then
                if(aa == 2*this%mode + 1) then
                  HYPRE_BUF(i) = u(this%mode + 3,0,1)
                  i = i + 1
                elseif(aa == 2*this%mode + 2) then
                  HYPRE_BUF(i) = u(this%mode + 3,0,2)
                  i = i + 1
                else
                  HYPRE_BUF(i) = 0.0
                  i = i + 1
                endif
              elseif( bb == aa - 2) then
                if(aa == 2*this%mode + 1) then
                  HYPRE_BUF(i) = -u(this%mode + 3,0,2)
                  i = i + 1
                elseif(aa == 2*this%mode + 2) then
                  HYPRE_BUF(i) = u(this%mode + 3,0,1)
                  i = i + 1
                else
                  HYPRE_BUF(i) = 0.0
                  i = i + 1
                endif
              else
                HYPRE_BUF(i) = 0.0
                i = i + 1
              endif
            elseif(this%mode == 1) then
              if( bb == aa + 4*this%mode + 1 ) then
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
        enddo
      enddo
      mm = ((4*this%mode + 1)*3)*(4*this%mode + 1) + 2*this%mode + 1
      do aa = 0, 4*this%mode
        if (aa /= 2*this%mode - 1 .and. aa /= 2*this%mode - 2 .and. aa /= 2*this%mode + 1 .and. aa /= 2*this%mode + 2) then
          i = mm + ((4*this%mode + 1)*3)*aa 
          HYPRE_BUF(i) = 0.0
        endif
      enddo 
    else
!       i = 1
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
              HYPRE_BUF(i)   = (-2.0 - ((n-1)/j)**2)/dr2 + u(1+this%mode,kk,1)
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
                    HYPRE_BUF(i) = u(1+this%mode+k-n,kk,1)
                    i = i + 1
                  elseif ( aa_offeset == 1  .and. bb_offeset == 0) then
                    HYPRE_BUF(i) = u(1+this%mode+k-n,kk,2)
                    i = i + 1
                  elseif ( aa_offeset == 0  .and. bb_offeset == 1) then
                    HYPRE_BUF(i) = -u(1+this%mode+k-n,kk,2)
                    i = i + 1
                  elseif ( aa_offeset == 0  .and. bb_offeset == 0) then
                    HYPRE_BUF(i) = u(1+this%mode+k-n,kk,1)
                    i = i + 1
                  endif
                endif
              endif
            endif
          elseif (bb == aa - 4*this%mode - 1) then
            HYPRE_BUF(i) = 1.0/dr2 - 0.5/(j*dr2)
            i = i + 1
          elseif (bb == aa + 4*this%mode + 1) then
            HYPRE_BUF(i) = 1.0/dr2 + 0.5/(j*dr2)
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
    write(2,*) u(:,:,:),"u"
    i = ((4*this%mode + 1)*3) * (4*this%mode + 1) + 1
    j = real(nofff) 
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
!               HYPRE_BUF(i) = (-2.0 - ((n-1)/j)**2)/dr2 + u(1+this%mode,kk,1)
              HYPRE_BUF(i) = (-2.0 - ((n + 1)/j)**2)/dr2 + u(1+this%mode,kk,1)
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
                    HYPRE_BUF(i) = u(1+this%mode+k-n,kk,1)
                    i = i + 1
                  elseif ( aa_offeset == 1  .and. bb_offeset == 0) then
                    HYPRE_BUF(i) = u(1+this%mode+k-n,kk,2)
                    i = i + 1
                  elseif ( aa_offeset == 0  .and. bb_offeset == 1) then
                    HYPRE_BUF(i) = -u(1+this%mode+k-n,kk,2)
                    i = i + 1
                  elseif ( aa_offeset == 0  .and. bb_offeset == 0) then
                    HYPRE_BUF(i) = u(1+this%mode+k-n,kk,1)
                    i = i + 1
                  endif
                endif
              endif
            endif
          elseif (bb == aa -  4*this%mode - 1) then
            HYPRE_BUF(i) = 1.0/dr2 - 0.5/(j*dr2)
            i = i + 1
          elseif (bb == aa +  4*this%mode + 1) then
            HYPRE_BUF(i) = 1.0/dr2 + 0.5/(j*dr2)
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
      mm = ((4*this%mode + 1)*3)*(4*this%mode + 1) + 2*this%mode + 1
      do aa = 0, 4*this%mode
        i = mm + ((4*this%mode + 1)*3)*aa 
        HYPRE_BUF(i) = 0.0
      enddo 
    else
!       i = 1
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
              HYPRE_BUF(i) = (-2.0 - ((n + 1)/j)**2)/dr2 + u(1+this%mode,kk,1)
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
                    HYPRE_BUF(i) = u(1+this%mode+k-n,kk,1)
                    i = i + 1
                  elseif ( aa_offeset == 1  .and. bb_offeset == 0) then
                    HYPRE_BUF(i) = u(1+this%mode+k-n,kk,2)
                    i = i + 1
                  elseif ( aa_offeset == 0  .and. bb_offeset == 1) then
                    HYPRE_BUF(i) = -u(1+this%mode+k-n,kk,2)
                    i = i + 1
                  elseif ( aa_offeset == 0  .and. bb_offeset == 0) then
                    HYPRE_BUF(i) = u(1+this%mode+k-n,kk,1)
                    i = i + 1
                  endif
                endif
              endif
            endif
          elseif (bb == aa - 4*this%mode - 1) then
            HYPRE_BUF(i) = 1.0/dr2 - 0.5/(j*dr2)
            i = i + 1
          elseif (bb == aa + 4*this%mode + 1) then
            HYPRE_BUF(i) = 1.0/dr2 + 0.5/(j*dr2)
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
        i = local_vol - (12*this%mode+3)*(4*this%mode + 1) + 1
        jmax = real(noff + nr)
        do aa = 0, 4*this%mode
          do bb = aa - 6*this%mode - 1, aa + 6*this%mode + 1
            if (aa == bb) then
              if( aa - 2*this%mode > 0 ) then
                g = (aa + 1)/2 - this%mode
                write(2,*) m,"m"
              else
                g = aa/2 - this%mode
                write(2,*) m,"m"
              endif
!               g = abs(g)
              HYPRE_BUF(i) = HYPRE_BUF(i) + (1.0-(g + 1)/jmax) * HYPRE_BUF(i + 4*this%mode + 1)
              i = i + 1
            elseif ( bb == aa + 4*this%mode + 1) then
!             elseif ( bb > aa ) then
              HYPRE_BUF(i) = 0.0
              i = i + 1
            else
              i = i + 1
            endif
          enddo
        enddo 
!         HYPRE_BUF(local_vol-1) = HYPRE_BUF(local_vol-1) + (1.0-(m+1)/jmax) * HYPRE_BUF(local_vol)
!         HYPRE_BUF(local_vol) = 0.0
      end select

    end select

  endif

!   call HYPRE_StructMatrixSetBoxValues( this%A, this%ilower, this%iupper, this%num_stencil, &
!     this%stencil_idx, HYPRE_BUF, ierr )
  ! 遍历每一行（网格点）
  write(2,*) 'this%solver_coef%solve||set_struct_matrix 0000'
  do i = this%ilower, this%iupper
   ! 行索引
   row_indices(1) = i

   ! 对于每个行索引，设置相应的列索引 (根据 stencil 结构)
    do ii = 1, this%num_stencil
      col_indices(ii) = i + ii - 1  ! 计算列索引，邻接点根据 stencil 偏移
      ! 插入非零值到 IJMatrix 中
!       write(2,*) 'this%solver_coef%solve||set_struct_matrix iiii'
      call HYPRE_IJMatrixSetValues(this%A_ij, 1, 1, row_indices, col_indices(ii), HYPRE_BUF(ii), ierr)
    end do
  end do
  write(2,*) 'this%solver_coef%solve||set_struct_matrix 2222'
  ! 完成矩阵的值插入后，进行矩阵组装
  call HYPRE_IJMatrixAssemble(this%A_ij, ierr)
  write(2,*) 'this%solver_coef%solve||set_struct_matrix 333'
  ! 将 IJMatrix 转换为 ParCSRMatrix

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine set_struct_matrix

end module field_solver_all_modesv1_class
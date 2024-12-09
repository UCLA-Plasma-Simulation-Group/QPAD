module field_solver_all_modesv4_class

use options_class 
use parallel_module
use mpi
use param
use sysutil_module
use debug_tool
use ufield_class

implicit none

private

public :: field_solver_all_modesv4
public :: HYPRE_BUF

character(len=32), parameter :: cls_name = "field_solver_all_modesv2"
integer, parameter :: cls_level = 4

real, dimension(:), pointer, save :: HYPRE_BUF => null()
real, dimension(:), pointer :: src_sol => null()

type :: field_solver_all_modesv4 ! class for HYPRE solver

  ! HYPRE parameters
  integer, dimension(:,:), pointer :: offsets => null()
  integer, dimension(:), pointer :: stencil_idx => null()
  integer :: num_stencil
  integer :: stype ! currently there's only cyclic reduction
  integer :: kind
  integer :: mode
  integer :: bnd

!   integer(HYPRE_TYPE) :: A, b, x, grid, stencil, solver
  integer(HYPRE_TYPE) :: A, b, x, grid, stencil, solver, precond
  integer, dimension(2) :: iupper, ilower

  contains

  procedure :: new => init_field_solver
  procedure :: del => end_field_solver 
  procedure :: solve => solve_equation
  procedure, private :: set_struct_solver
  procedure, private :: set_struct_grid
  procedure, private :: set_struct_stencil
  procedure, private :: set_struct_matrix 
  procedure, private :: set_hypre_src_coef 
  procedure, private :: set_hypre_src_b 
  procedure, private :: get_hypre_src_coef 
  procedure, private :: get_hypre_src_b 
 
end type field_solver_all_modesv4

integer, save :: nofff
real, save :: drr
integer :: nr
 
contains

! =====================================================================
! Class field_solver implementation
! =====================================================================

subroutine init_field_solver( this, opts, max_mode, dr, kind, bnd, stype )

  implicit none

  class( field_solver_all_modesv4 ), intent(inout) :: this
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
  write(2,*) nr,'nr'

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

  class( field_solver_all_modesv4 ), intent(inout) :: this

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

  select case(this%kind)
  case(p_fk_coef)
    call HYPRE_StructBiCGSTABDestroy( this%solver, ierr )
  case(p_fk_all_B_plus,p_fk_all_B_minus)
    call HYPRE_StructBiCGSTABDestroy( this%solver, ierr )
  end select  

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine end_field_solver

subroutine set_struct_solver( this )

  implicit none

  class( field_solver_all_modesv4 ), intent(inout) :: this

  integer :: ierr, comm
  external HYPRE_BoomerAMGSolve, HYPRE_BoomerAMGSetup
  character(len=32), save :: sname = "set_struct_solver"

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  comm = comm_loc()

  select case(this%kind)
  case( p_fk_coef )
    call HYPRE_StructBiCGSTABCreate( comm, this%solver, ierr )
    call HYPRE_StructBiCGSTABSetup( this%solver, this%A, this%b, this%x, ierr )
  case( p_fk_all_B_minus, p_fk_all_B_plus)
  !   write(2,*) 'solver Initializing 2'
    ! call write_stdout( 'mode '//num2str(this%mode)//': Using Cyclic Reduction solver' )
    call HYPRE_StructBiCGSTABCreate( comm, this%solver, ierr )
!     call HYPRE_StructGMRESSetTol(this%solver, 1.0e-15, ierr)
! !     call HYPRE_StructGMRESSetMaxIter(this%solver, 100000, ierr)
!   !  创建 BoomerAMG 预处理器
!     call HYPRE_BoomerAMGCreate(this%precond, ierr)
!     call HYPRE_BoomerAMGSetCycleType(this%precond, 2, ierr)  ! 2 W-cycle 1 V-cycle
!     call HYPRE_BoomerAMGSetMaxLevels(this%precond, 50, ierr)  ! 设置最大网格层数
!     call HYPRE_BoomerAMGSetRelaxType(this%precond, 6, ierr)  ! Schwarz 或 Symmetric Gauss-Seidel
! !   !   设置 BoomerAMG 参数，例如最大迭代次数和容忍度
!     call HYPRE_BoomerAMGSetMaxIter(this%precond, 10000, ierr)  ! 预处理器的最大迭代次数
!     call HYPRE_BoomerAMGSetTol(this%precond, 1.0e-6, ierr)    ! 预处理器的容忍度（GMRES 控制收敛）
!   !   设置 BoomerAMG 作为 GMRES 的预处理器
!     call HYPRE_StructGMRESSetPrecond(this%solver, HYPRE_BoomerAMGSolve, HYPRE_BoomerAMGSetup, this%precond, ierr)
  !   write(2,*) 'solver Initializing 3'
  !   write(2,*) this%A,'A'
  !   write(2,*) this%b,'b'
  !   write(2,*) this%x,'x'
    call HYPRE_StructBiCGSTABSetup( this%solver, this%A, this%b, this%x, ierr )
!   write(2,*) 'solver Initializing 4'
  end select

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine set_struct_solver

subroutine solve_equation( this, src, psi_re, qe_re, qn1_re, psi_im, qe_im, qn1_im, u, qbm)

  implicit none

  class( field_solver_all_modesv4 ), intent(inout) :: this
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
  integer :: size_sol
  real, dimension(:,:), pointer :: f1_re => null(), f1_im => null()
  character(len=32), save :: sname = "solve_equation"

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  size_sol = nr* (4 * this%mode + 1)
!   write(2,*) size_sol, 'size_sol' 

  if (.not. associated( src_sol)) then
    allocate( src_sol(size_sol) )
  end if
  src_sol(:) = 0.0

!   write(2,*) this%kind, 'kind' 
  select case ( this%kind )
  case ( p_fk_coef )
    call set_struct_matrix(this, psi_re = psi_re, psi_im = psi_im, qbm=qbm)
    call this%set_struct_solver()
    call this%set_hypre_src_coef( qe_re, qn1_re, qe_im, qn1_im, qbm )
!     write(2,*) src_sol,"src coef"
!     write(2,*) HYPRE_BUF(:),"HYPRE_BUF"
!     write(2,*) 'this%solver_coef%solve||HYPRE_StructVectorAssembl'   
    call HYPRE_StructVectorAssemble( this%b, ierr )
!     write(2,*) 'this%solver_coef%solve||HYPRE_StructCycRedSolve'  
    call HYPRE_StructBiCGSTABSolve( this%solver, this%A, this%b, this%x, ierr )
!     call HYPRE_HYPRE_StructBiCGSTABGetNumIterations(this%solver,this%numiterations)
!     write(2,*) 'this%numiterations'
!     write(2,*) 'this%solver_coef%solve||HYPRE_StructVectorGetBoxValues 0' 
    call this%get_hypre_src_coef( src )
!     write(2,*) sum(src),"src sum(src)"
!     write(2,*) src,"src sum(src) coef end"
  case ( p_fk_all_B_minus, p_fk_all_B_plus )
    call set_struct_matrix(this, u=u)
    call this%set_struct_solver()
! !     write(2,*) HYPRE_BUF(:),"HYPRE_BUF"
!     write(2,*) "src"  
    call this%set_hypre_src_b(src)
!     write(2,*) "src_sol 1"
    call HYPRE_StructVectorAssemble( this%b, ierr )
!     write(2,*) "src_sol 2"
    call HYPRE_StructBiCGSTABSolve( this%solver, this%A, this%b, this%x, ierr )
!     call HYPRE_HYPRE_StructBiCGSTABGetNumIterations(this%solver,this%numiterations)
!     write(2,*) 'this%numiterations'
    call this%get_hypre_src_b(src)
!     write(2,*) sum(src),"src sum(src)" 
!     write(2,*) src,"src end" 
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

  class( field_solver_all_modesv4 ), intent(inout) :: this
  type( options ), intent(in) :: opts

  integer :: comm, ierr
  character(len=32), save :: sname = "set_struct_grid"

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  comm = comm_loc()
  this%ilower = (/1,opts%get_noff(1) + 1/)
  this%iupper = (/1 + 4*this%mode,opts%get_noff(1) + opts%get_ndp(1)/) 

  call HYPRE_StructGridCreate( comm, 2, this%grid, ierr )
  call HYPRE_StructGridSetExtents( this%grid, this%ilower, this%iupper, ierr )
  call HYPRE_StructGridAssemble( this%grid, ierr )

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine set_struct_grid

subroutine set_struct_stencil( this )

  implicit none

  class( field_solver_all_modesv4 ), intent(inout) :: this

  integer :: i, j, k, l, ierr
  character(len=32), save :: sname = "set_struct_stencil"

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  select case(this%kind)

  case(p_fk_coef)
    
    this%num_stencil = 4*this%mode + 1
    if ( .not. associated( this%offsets ) ) then
      allocate( this%offsets(2,this%num_stencil ) )
    endif

    if ( .not. associated( this%stencil_idx ) ) then
      allocate( this%stencil_idx( this%num_stencil ) )
      do i = 1, this%num_stencil
        this%stencil_idx(i) = i-1
      enddo
    endif

    k = 1
    do j = 0, 4*this%mode
      this%offsets(:,k) = (/j - 2*this%mode,0/)
      k=k + 1
    enddo

    call HYPRE_StructStencilCreate( 2, this%num_stencil, this%stencil, ierr )
    do i = 1, this%num_stencil
      call HYPRE_StructStencilSetElement( this%stencil, i-1, this%offsets(1,i), ierr )
    enddo
  
  case( p_fk_all_B_minus, p_fk_all_B_plus )

    this%num_stencil = 4*this%mode + 3
    if ( .not. associated( this%offsets ) ) then
      allocate( this%offsets( 2,this%num_stencil ) )
    endif

    if ( .not. associated( this%stencil_idx ) ) then
      allocate( this%stencil_idx( this%num_stencil ) )
      do i = 1, this%num_stencil
        this%stencil_idx(i) = i-1
      enddo
    endif

    k = 3
    this%offsets(:,1) = (/0,-1/)
    this%offsets(:,2) = (/0,1/)
    do j = 0, 4*this%mode
      this%offsets(:,k) = (/j - 2*this%mode,0/)
      k=k + 1
    enddo

    call HYPRE_StructStencilCreate( 2, this%num_stencil, this%stencil, ierr )
    do i = 1, this%num_stencil
      call HYPRE_StructStencilSetElement( this%stencil, i-1, this%offsets(1,i), ierr )
    enddo

  end select

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine set_struct_stencil

subroutine set_struct_matrix( this, psi_re, psi_im, u, qbm)

  implicit none

  class( field_solver_all_modesv4 ), intent(inout) :: this
  type(ufield), intent(inout), dimension(:), optional, pointer :: psi_re
  type(ufield), intent(inout), dimension(:), optional, pointer :: psi_im
  real, intent(inout), dimension(:,:,:), optional :: u
  real, intent(inout), optional :: qbm

  integer :: i, ierr, local_vol, noff, m, aa, bb, nn, kk, mm, g
  integer :: k, n, d, bb_offeset, aa_offeset,l, demode, Anumber
  integer :: comm, lidproc, lnvp
  real :: dr2, m2, j, jmax, dr
  real, dimension(:,:), pointer :: f1_re => null(), f1_im => null()
  real, dimension(:,:), pointer :: f2_re => null(), f2_im => null()
  character(len=32), save :: sname = "set_struct_matrix"
!   write(2,*) '-1'
  call write_dbg( cls_name, sname, cls_level, 'starts' )

!   write(2,*) '0'
  comm    = comm_loc()
  lidproc = id_proc_loc()
  lnvp    = num_procs_loc()
  noff   = nofff
  dr     = drr
!   write(2,*) '1'
  dr2 = drr*drr
  m = this%mode
  m2 = real(m*m)
  local_vol = (this%iupper(1) - this%ilower(1) + 1)*(this%iupper(2) - this%ilower(2) + 1)*this%num_stencil

  if ( .not. associated( HYPRE_BUF ) ) then
    allocate( HYPRE_BUF( local_vol ) )
  elseif ( size(HYPRE_BUF) < local_vol ) then
    deallocate( HYPRE_BUF )
    allocate( HYPRE_BUF( local_vol ) )
  endif

  call HYPRE_StructMatrixCreate( comm, this%grid, this%stencil, this%A, ierr )
  call HYPRE_StructMatrixInitialize( this%A, ierr )
  write(2,*) size(HYPRE_BUF),'size(HYPRE_BUF)'

  ! set the matrix element and lower boundary
  select case ( this%kind )

  case ( p_fk_coef )

    do nn = 1, nr  
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
              HYPRE_BUF(i) = 1.0 - qbm * f1_re(1,nn)
            else
              if (demode > this%mode) then
                HYPRE_BUF(i) = 0.0
              else
                if (demode == 0) then
                  HYPRE_BUF(i) = 0.0
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
                    HYPRE_BUF(i) = -qbm * f1_re(1,nn)
                  elseif ( aa_offeset == 1  .and. bb_offeset == 0) then
                    f1_im => psi_im(demode)%get_f1()
                    HYPRE_BUF(i) = -sign(1,k-n)*qbm * f1_im(1,nn)
                  elseif ( aa_offeset == 0  .and. bb_offeset == 1) then
                    f1_im => psi_im(demode)%get_f1()
                    HYPRE_BUF(i) = sign(1,k-n)*qbm * f1_im(1,nn)
                  elseif ( aa_offeset == 0  .and. bb_offeset == 0) then
                    f1_re => psi_re(demode)%get_f1()
                    HYPRE_BUF(i) = -qbm * f1_re(1,nn)
                  endif
                endif
              endif
            endif
          else            
            HYPRE_BUF(i) = 0.0
          endif
          i = i + 1
        enddo
      enddo
    enddo

  case( p_fk_all_B_minus )

    ! set the first grid point of each partition
    i = (4*this%mode + 3) * (4*this%mode + 1) + 1
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
        HYPRE_BUF(i) = 1.0/dr2 - 0.5/(j*dr2)
        i = i + 1
        HYPRE_BUF(i) = 1.0/dr2 + 0.5/(j*dr2)
        i = i + 1
        do bb = aa - 2*this%mode, aa + 2*this%mode
          if(  (bb >= 0) .and. (bb <= 4*this%mode) ) then
            if( bb - 2*this%mode > 0 ) then
              n = (bb + 1)/2 - this%mode 
            else
              n = bb/2 - this%mode 
            endif
            demode = abs(k-n)
            if (aa == bb) then
              HYPRE_BUF(i)   = (-2.0 - ((n-1)/j)**2)/dr2 + u(1+this%mode,nn,1)
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
                    HYPRE_BUF(i) = u(1+this%mode+k-n,nn,1)
!                     write(2,*) u(1+this%mode+k-n,kk,1),"u(1+this%mode+k-n,kk,1)"
                    i = i + 1
                  elseif ( aa_offeset == 1  .and. bb_offeset == 0) then
                    HYPRE_BUF(i) = u(1+this%mode+k-n,nn,2)
!                     write(2,*) u(1+this%mode+k-n,kk,2),"u(1+this%mode+k-n,kk,2)"
                    i = i + 1
                  elseif ( aa_offeset == 0  .and. bb_offeset == 1) then
                    HYPRE_BUF(i) = -u(1+this%mode+k-n,nn,2)
                    i = i + 1
                  elseif ( aa_offeset == 0  .and. bb_offeset == 0) then
                    HYPRE_BUF(i) = u(1+this%mode+k-n,nn,1)
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
    enddo
! inner boundary
    i = 1
    if (lidproc == 0) then
      do aa = 0, 4*this%mode
        HYPRE_BUF(i) = 0.0
        i = i + 1
        if (aa /= 2*this%mode - 1 .and. aa /= 2*this%mode - 2 .and. aa /= 2*this%mode + 1 .and. aa /= 2*this%mode + 2) then
          HYPRE_BUF(i) = 0.0
        elseif(aa == 2*this%mode - 1) then
          HYPRE_BUF(i) = -4.0/dr2
        else
          HYPRE_BUF(i) = 4.0/dr2
        endif
        i = i + 1 
        do bb = aa - 2*this%mode, aa + 2*this%mode
          if (aa == bb) then
!             m=-1
            if (aa == 2*this%mode-1) then
              if(this%mode > 1) then
                HYPRE_BUF(i) = u(1 + this%mode,1,1) - u(this%mode + 3,1,1) + 4.0/dr2
                i = i + 1
              else
                HYPRE_BUF(i) = u(1 + this%mode,1,1) + 4.0/dr2
!                 HYPRE_BUF(i) = -4.0/dr2
                i = i + 1
              endif
            elseif(aa == 2*this%mode-2) then
              if(this%mode > 1) then
                HYPRE_BUF(i) = u(1+this%mode,1,1) + u(this%mode + 3,1,1) - 4.0/dr2
                i = i + 1
              else
                HYPRE_BUF(i) = u(1 + this%mode,1,1) - 4.0/dr2
!                 HYPRE_BUF(i) = -4.0/dr2
                i = i + 1
              endif
!           m=1
!             elseif(aa == 2*this%mode+1) then
            elseif(aa == 2*this%mode+1 .or. aa == 2*this%mode+2) then
              HYPRE_BUF(i) = u(1 + this%mode,1,1) - 4.0/dr2
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
              if( bb == aa - 1) then
                if(aa == 2*this%mode - 1) then
                  HYPRE_BUF(i) = -u(this%mode + 3,1,2)
                  i = i + 1
                else
                  HYPRE_BUF(i) = 0.0
                  i = i + 1
                endif
              elseif( bb == aa + 1) then
                if(aa == 2*this%mode - 2) then
                  HYPRE_BUF(i) = -u(this%mode + 3,1,2)
                  i = i + 1
                else
                  HYPRE_BUF(i) = 0.0
                  i = i + 1
                endif
              elseif( bb == aa - 3) then
                if(aa == 2*this%mode + 1) then
                  HYPRE_BUF(i) = u(this%mode + 3,1,1)
                  i = i + 1
                elseif(aa == 2*this%mode + 2) then
                  HYPRE_BUF(i) = u(this%mode + 3,1,2)
                  i = i + 1
                else
                  HYPRE_BUF(i) = 0.0
                  i = i + 1
                endif
              elseif( bb == aa - 2) then
                if(aa == 2*this%mode + 1) then
                  HYPRE_BUF(i) = -u(this%mode + 3,1,2)
                  i = i + 1
                elseif(aa == 2*this%mode + 2) then
                  HYPRE_BUF(i) = u(this%mode + 3,1,1)
                  i = i + 1
                else
                  HYPRE_BUF(i) = 0.0
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
      mm = (4*this%mode + 3)*(4*this%mode + 1) + 1
      do aa = 0, 4*this%mode
        if (aa /= 2*this%mode - 1 .and. aa /= 2*this%mode - 2 .and. aa /= 2*this%mode + 1 .and. aa /= 2*this%mode + 2) then
          i = mm + (4*this%mode + 3)*aa 
          HYPRE_BUF(i) = 0.0
        endif
      enddo 
    else
      i = 1
      j = real(noff) 
      kk = int(j)
      do aa = 0, 4*this%mode
        if( aa - 2*this%mode > 0 ) then
          k = (aa + 1)/2 - this%mode
        else
          k = aa/2 - this%mode
        endif
        HYPRE_BUF(i) = 1.0/dr2 - 0.5/(j*dr2)
        i = i + 1
        HYPRE_BUF(i) = 1.0/dr2 + 0.5/(j*dr2)
        i = i + 1 
        do bb = aa - 2*this%mode, aa + 2*this%mode
          if(  (bb >= 0) .and. (bb <= 4*this%mode) ) then
            if( bb - 2*this%mode > 0 ) then
              n = (bb + 1)/2 - this%mode 
            else
              n = bb/2 - this%mode 
            endif
            demode = abs(k-n)
            if (aa == bb) then
              HYPRE_BUF(i)   = (-2.0 - ((n-1)/j)**2)/dr2 + u(1+this%mode,1,1)
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
                    HYPRE_BUF(i) = u(1+this%mode+k-n,1,1)
                    i = i + 1
                  elseif ( aa_offeset == 1  .and. bb_offeset == 0) then
                    HYPRE_BUF(i) = u(1+this%mode+k-n,1,2)
                    i = i + 1
                  elseif ( aa_offeset == 0  .and. bb_offeset == 1) then
                    HYPRE_BUF(i) = -u(1+this%mode+k-n,1,2)
                    i = i + 1
                  elseif ( aa_offeset == 0  .and. bb_offeset == 0) then
                    HYPRE_BUF(i) = u(1+this%mode+k-n,1,1)
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
    endif 

  case( p_fk_all_B_plus )

    ! set the first grid point of each partition
!     write(2,*) u(:,:,:),"u"
    i = (4*this%mode + 3) * (4*this%mode + 1) + 1
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
        HYPRE_BUF(i) = 1.0/dr2 - 0.5/(j*dr2)
        i = i + 1
        HYPRE_BUF(i) = 1.0/dr2 + 0.5/(j*dr2)
        i = i + 1 
        do bb = aa - 2*this%mode, aa + 2*this%mode
          if(  (bb >= 0) .and. (bb <= 4*this%mode) ) then
            if( bb - 2*this%mode > 0 ) then
              n = (bb + 1)/2 - this%mode 
            else
              n = bb/2 - this%mode 
            endif
            demode = abs(k-n)
            if (aa == bb) then
!               HYPRE_BUF(i) = (-2.0 - ((n-1)/j)**2)/dr2 + u(1+this%mode,kk,1)
              HYPRE_BUF(i) = (-2.0 - ((n + 1)/j)**2)/dr2 + u(1+this%mode,nn,1)
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
                    HYPRE_BUF(i) = u(1+this%mode+k-n,nn,1)
                    i = i + 1
                  elseif ( aa_offeset == 1  .and. bb_offeset == 0) then
                    HYPRE_BUF(i) = u(1+this%mode+k-n,nn,2)
                    i = i + 1
                  elseif ( aa_offeset == 0  .and. bb_offeset == 1) then
                    HYPRE_BUF(i) = -u(1+this%mode+k-n,nn,2)
                    i = i + 1
                  elseif ( aa_offeset == 0  .and. bb_offeset == 0) then
                    HYPRE_BUF(i) = u(1+this%mode+k-n,nn,1)
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
    enddo
! inner boundary
    i = 1
    if (lidproc == 0) then
      do aa = 0, 4*this%mode
        HYPRE_BUF(i) = 0.0
        i = i + 1
        HYPRE_BUF(i) = 0.0
        i = i + 1
        do bb = aa - 2*this%mode, aa + 2*this%mode
          if (aa == bb) then
            HYPRE_BUF(i)   = 1.0/dr2
            i = i + 1
          else
            HYPRE_BUF(i) = 0.0
            i = i + 1
          endif
        enddo
      enddo
      mm = (4*this%mode + 3)*(4*this%mode + 1) + 1
      do aa = 0, 4*this%mode
        i = mm + (4*this%mode + 3)*aa 
        HYPRE_BUF(i) = 0.0
      enddo 
    else
      i = 1
      j = real(noff) 
      kk = int(j)
      do aa = 0, 4*this%mode
        if( aa - 2*this%mode > 0 ) then
          k = (aa + 1)/2 - this%mode
        else
          k = aa/2 - this%mode
        endif
        HYPRE_BUF(i) = 1.0/dr2 - 0.5/(j*dr2)
        i = i + 1
        HYPRE_BUF(i) = 1.0/dr2 + 0.5/(j*dr2)
        i = i + 1 
        do bb = aa - 2*this%mode, aa + 2*this%mode
          if(  (bb >= 0) .and. (bb <= 4*this%mode) ) then
            if( bb - 2*this%mode > 0 ) then
              n = (bb + 1)/2 - this%mode 
            else
              n = bb/2 - this%mode 
            endif
            demode = abs(k-n)
            if (aa == bb) then
              HYPRE_BUF(i) = (-2.0 - ((n + 1)/j)**2)/dr2 + u(1+this%mode,1,1)
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
                    HYPRE_BUF(i) = u(1+this%mode+k-n,1,1)
                    i = i + 1
                  elseif ( aa_offeset == 1  .and. bb_offeset == 0) then
                    HYPRE_BUF(i) = u(1+this%mode+k-n,1,2)
                    i = i + 1
                  elseif ( aa_offeset == 0  .and. bb_offeset == 1) then
                    HYPRE_BUF(i) = -u(1+this%mode+k-n,1,2)
                    i = i + 1
                  elseif ( aa_offeset == 0  .and. bb_offeset == 0) then
                    HYPRE_BUF(i) = u(1+this%mode+k-n,1,1)
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
    endif 

  end select

  ! set the upper boundary
  if ( lidproc == lnvp-1 ) then

    select case ( this%kind )

    case ( p_fk_all_B_minus, p_fk_all_B_plus )
      
      select case ( this%bnd )

      case ( p_bnd_zero )

        i = local_vol - (4*this%mode+3)*(4*this%mode + 1) + 1
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
        i = local_vol - (4*this%mode +3)*(4*this%mode + 1) + 1
        jmax = real(noff + nr)
        do aa = 0, 4*this%mode
          i = i + 2
          do bb = aa - 2*this%mode, aa + 2*this%mode
            if (aa == bb) then
              if( aa - 2*this%mode > 0 ) then
                g = (aa + 1)/2 - this%mode
!                 write(2,*) m,"m"
              else
                g = aa/2 - this%mode
!                 write(2,*) m,"m"
              endif
!               g = abs(g)
              HYPRE_BUF(i) = HYPRE_BUF(i) + (1.0-(g + 1)/jmax) * HYPRE_BUF(i - 2*this%mode - 1)
              HYPRE_BUF(i - 2*this%mode - 1) = 0.0
              i = i + 1
            else
              i = i + 1
            endif
          enddo
        enddo 

      end select

    end select

  endif

  call HYPRE_StructMatrixSetBoxValues( this%A, this%ilower, this%iupper, this%num_stencil, &
    this%stencil_idx, HYPRE_BUF, ierr )

  call HYPRE_StructMatrixAssemble( this%A, ierr )

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine set_struct_matrix

subroutine set_hypre_src_coef(this, qe_re, qn1_re, qe_im, qn1_im, qbm )

  implicit none

  class( field_solver_all_modesv4 ), intent(inout) :: this
  type(ufield), intent(inout), dimension(:), optional, pointer :: qe_re
  type(ufield), intent(inout), dimension(:), optional, pointer :: qn1_re
  type(ufield), intent(inout), dimension(:), optional, pointer :: qe_im
  type(ufield), intent(inout), dimension(:), optional, pointer :: qn1_im
  real, intent(inout), optional :: qbm

  integer :: ierr, a, im, k, l, kk, m
  real :: jj, j
  real, dimension(:,:), pointer :: f1_re => null(), f1_im => null()
  real, dimension(:,:), pointer :: f2_re => null(), f2_im => null()
  real, dimension(:,:), pointer :: f3_re => null(), f3_im => null()
  character(len=32), save :: sname = "set_hypre_src_coef"

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  a = 1
  m = this%mode
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

  call HYPRE_StructVectorSetBoxValues( this%b, this%ilower, this%iupper, src_sol, ierr )
!   write(2,*) size(src_sol),'size(src_sol)'

end subroutine set_hypre_src_coef

subroutine set_hypre_src_b(this, src)

  implicit none

  class( field_solver_all_modesv4 ), intent(inout) :: this
  real, intent(inout), dimension(:,:,:), optional, pointer :: src

  integer :: ierr, a, m, k, l, kk, j, im
  character(len=32), save :: sname = "set_hypre_src_b"

  a = 1
  j = real(nofff) + 1.0
!   write(2,*) j,'j'
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
  
  call HYPRE_StructVectorSetBoxValues( this%b, this%ilower, this%iupper, src_sol, ierr )
  write(2,*) size(src_sol),'size(src_sol)'
!   write(2,*) src_sol,'src_sol'

end subroutine set_hypre_src_b

subroutine get_hypre_src_coef(this, src)

  implicit none

  class( field_solver_all_modesv4 ), intent(inout) :: this
  real, intent(inout), dimension(:,:,:), optional, pointer :: src

  integer :: ierr, a, m, k, l, im
  character(len=32), save :: sname = "get_hypre_src_coef"

  call HYPRE_StructVectorGetBoxValues( this%x, this%ilower, this%iupper, src_sol, ierr )
!   call HYPRE_IJVectorPrint( this%x,"struct_x_output.out")
!   call HYPRE_StructVectorPrint( "struct_x_output.out", this%x)
!   write(2,*) size(src_sol),'size(src_sol)'
!   write(2,*) src_sol(:),'coef src_sol'
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
!   write(2,*) src(:,:,:),'coef src'

end subroutine get_hypre_src_coef

subroutine get_hypre_src_b(this, src)

  implicit none

  class( field_solver_all_modesv4 ), intent(inout) :: this
  real, intent(inout), dimension(:,:,:), optional, pointer :: src

  integer :: ierr, a, m, k, l, im
  character(len=32), save :: sname = "get_hypre_src_b"

  call HYPRE_StructVectorGetBoxValues( this%x, this%ilower, this%iupper, src_sol, ierr )
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

end module field_solver_all_modesv4_class
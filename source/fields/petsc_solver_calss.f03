module petsc_solver_class
  use petscksp
  use mpi
  implicit none

  type :: PetscSolver
     type(Mat) :: A        ! 矩阵 A
     type(Vec) :: b, x     ! 向量 b 和 x
     type(KSP) :: ksp      ! 求解器
     type(PC) :: pc        ! 预条件器
     integer :: n          ! 方程规模
  contains
     procedure :: initialize
     procedure :: set_matrix
     procedure :: solve
  end type PetscSolver

contains

  ! 初始化 PETSc 环境
  subroutine initialize(this, n)
    class(PetscSolver), intent(inout) :: this
    integer, intent(in) :: n
    integer :: ierr

    call MPI_Init(ierr)
    call KSPCreate(PETSC_COMM_WORLD, this%ksp)  ! 确保ksp正确初始化
    this%n = n
    call VecCreateSeq(PETSC_COMM_SELF, n, this%b)
    call VecCreateSeq(PETSC_COMM_SELF, n, this%x)
    call MatCreateSeqAIJ(PETSC_COMM_SELF, n, n, 3, NULL(), this%A)
  end subroutine initialize

  ! 设置矩阵 A 和向量 b
  subroutine set_matrix(this, A_values, b_values)
    class(PetscSolver), intent(inout) :: this
    real(8), intent(in) :: A_values(:), b_values(:)
    integer :: ierr, i, j

    ! 设置矩阵 A 和向量 b 的值
    do i = 1, this%n
       call VecSetValue(this%b, i-1, b_values(i), INSERT_VALUES)
    end do
    do i = 1, this%n
       do j = 1, 3
          call MatSetValue(this%A, i-1, j-1, A_values(i), INSERT_VALUES)
       end do
    end do
    call VecAssemblyBegin(this%b)
    call VecAssemblyEnd(this%b)
    call MatAssemblyBegin(this%A, MAT_FINAL_ASSEMBLY)
    call MatAssemblyEnd(this%A, MAT_FINAL_ASSEMBLY)
  end subroutine set_matrix

  ! 求解线性方程组
  subroutine solve(this)
    class(PetscSolver), intent(inout) :: this
    integer :: ierr

    ! 设置求解器和预条件器
    call KSPSetOperators(this%ksp, this%A, this%A)
    call KSPSetTolerances(this%ksp, 1.0d-5, 1.0d-50, 1000, 0)
    call KSPSetType(this%ksp, KSPGMRES)
    call KSPGetPC(this%ksp, this%pc)
    call PCSetType(this%pc, PCJACOBI)

    ! 求解
    call KSPSolve(this%ksp, this%b, this%x)

    ! 输出结果
    call VecView(this%x, PETSC_VIEWER_STDOUT_WORLD)
  end subroutine solve

end module petsc_solver_class
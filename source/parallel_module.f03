module parallel_module

use mpi
use omp_lib

implicit none

private

public :: init_parallel, end_parallel, init_pipeline
public :: num_procs, id_proc, comm_world
public :: num_procs_loc, id_proc_loc, id_stage, num_stages, comm_loc, root_loc
public :: ntag

interface init_parallel
  module procedure init_parallel
end interface init_parallel

interface init_pipeline
  module procedure init_pipeline
end interface init_pipeline

interface end_parallel
  module procedure end_parallel
end interface end_parallel

interface num_procs
  module procedure num_procs
end interface num_procs

interface id_proc
  module procedure id_proc
end interface id_proc

interface comm_world
  module procedure comm_world
end interface comm_world

interface num_procs_loc
  module procedure num_procs_loc
end interface num_procs_loc

interface id_proc_loc
  module procedure id_proc_loc
end interface id_proc_loc

interface root_loc
  module procedure root_loc
end interface root_loc

interface id_stage
  module procedure id_stage
end interface id_stage

interface num_stages
  module procedure num_stages
end interface num_stages

interface comm_loc
  module procedure comm_loc
end interface comm_loc

interface ntag
  module procedure ntag
end interface ntag

! module variables
! mpi data types
integer, public :: p_dtype_real, p_dtype_int, p_dtype_double, p_dtype_cplx, p_dtype_char
! total number of processors
integer :: num_procs_
! processor id
integer :: id_proc_
! mpi global communicator
integer :: comm_world_
! number of pipeline stages
integer :: num_stages_
! number of processors in local stage
integer :: num_procs_loc_
! processor id in local stage
integer :: id_proc_loc_
! stage id
integer :: id_stage_
! local mpi communicator
integer :: comm_loc_

contains

function num_procs()

  implicit none

  integer :: num_procs

  num_procs = num_procs_

end function num_procs

function id_proc()

  implicit none

  integer :: id_proc

  id_proc = id_proc_

end function id_proc

function comm_world()

  implicit none

  integer :: comm_world

  comm_world = comm_world_

end function comm_world

function num_stages()

  implicit none

  integer :: num_stages

  num_stages = num_stages_

end function num_stages

function num_procs_loc()

  implicit none

  integer :: num_procs_loc

  num_procs_loc = num_procs_loc_

end function num_procs_loc

function id_proc_loc()

  implicit none

  integer :: id_proc_loc

  id_proc_loc = id_proc_loc_

end function id_proc_loc

function root_loc()

  implicit none

  integer :: root_loc

  root_loc = num_procs_loc_ * id_stage_

end function root_loc

function id_stage()

  implicit none

  integer :: id_stage

  id_stage = id_stage_

end function id_stage

function comm_loc()

  implicit none

  integer :: comm_loc

  comm_loc = comm_loc_

end function comm_loc

subroutine init_parallel( n_threads )

  implicit none

  integer, intent(in) :: n_threads

  integer :: ierr
  logical :: flag
  real :: prec

  ! initialize for shared memory parallel processing using openmp
  call init_omp(n_threads)

  ! initialize basic mpi environment
  call mpi_initialized( flag, ierr )
  if (.not. flag) then
     call mpi_init(ierr)
     if (ierr /= 0) stop
  endif
  comm_world_ = MPI_COMM_WORLD
  ! determine the rank of the calling process in the communicator
  call mpi_comm_rank( comm_world_, id_proc_, ierr )
  ! determine the size of the group associated with a communicator
  call mpi_comm_size( comm_world_, num_procs_, ierr )

  ! set default data types
  p_dtype_int = MPI_INTEGER
  p_dtype_double  = MPI_DOUBLE_PRECISION
  p_dtype_char    = MPI_CHARACTER
  
  if ( digits(prec) > 24 ) then
    ! double precision real
    p_dtype_real    = MPI_DOUBLE_PRECISION
    p_dtype_cplx = MPI_DOUBLE_COMPLEX
  else
    ! single precision real
    p_dtype_real    = MPI_REAL
    p_dtype_cplx = MPI_COMPLEX
  endif

end subroutine init_parallel

subroutine init_pipeline( n_stages )

  implicit none

  integer, intent(in) :: n_stages

  integer :: ierr

  ! initialize pipeline environment
  num_stages_     = n_stages
  num_procs_loc_ = num_procs_ / n_stages
  id_proc_loc_   = mod( id_proc_, num_procs_loc_ )
  id_stage_      = id_proc_ / num_procs_loc_

  call mpi_comm_split( comm_world_, id_stage_, id_proc_loc_, comm_loc_, ierr )
  call mpi_comm_rank( comm_loc_, id_proc_loc_, ierr )
  call mpi_comm_size( comm_loc_, num_procs_loc_, ierr )

  ! initialize pseudo-random number sequence
  call random_seed()

end subroutine init_pipeline

subroutine init_omp(nth)
  ! initialize openmp library
  ! use nth threads if nth > 0; otherwise, use the number found
  implicit none
  integer, intent(in) :: nth
  ! local data
  integer :: ncpus, nthreads, tid
  ! determine how many processors are available
  ncpus = omp_get_num_procs()
  nthreads = omp_get_max_threads()
  if (nth > 0) nthreads = nth
  call omp_set_num_threads(nthreads)

  !$omp parallel private(tid)
  tid = omp_get_thread_num()
  !$omp end parallel

end subroutine init_omp

! this subroutine terminates parallel processing
subroutine end_parallel()
  
  implicit none

  integer :: ierr
  logical :: flag

  ! indicate whether mpi_init has been called
  call mpi_initialized( flag, ierr )
  if (flag) then
    ! synchronize processes
    call mpi_barrier( comm_world_, ierr )
    ! terminate mpi execution environment
    call mpi_finalize( ierr )
  endif

end subroutine end_parallel

function ntag()

   implicit none
   integer, save :: tag = 0
   integer :: ntag

   ntag = tag
   tag = tag + 1

end function ntag

end module parallel_module
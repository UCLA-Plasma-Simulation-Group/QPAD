module parallel_class

use mpi
use omp_lib

implicit none

private

public :: parallel
public :: ntag

interface ntag
  module procedure ntag
end interface ntag

type parallel

  private
  ! nvp: number of MPI nodes
  ! idproc: processor id
  ! kstrt: idproc+1
  ! mreal = default datatype for reals
  ! mint = default datatype for integers
  ! mcplx = default datatype for complex type
  ! mdouble = default double precision type
  ! lworld = MPI_COMM_WORLD communicator
  integer :: nvp
  integer :: idproc
  integer :: kstrt
  integer :: mreal, mint, mcplx, mdouble, mchar, lworld

  ! these datatypes are for huge array communication
  integer :: mreal_pack, mint_pack

  contains

  generic :: new => init_parallel
  generic :: del => end_parallel
  procedure :: getnvp
  procedure :: getidproc
  procedure :: getkstrt
  procedure :: getlworld
  procedure :: getmreal
  procedure :: getmint
  procedure :: getmdouble
  procedure :: getmcplx
  procedure :: getmchar
  procedure, private :: init_parallel
  procedure, private :: end_parallel

end type parallel

contains
!
function getnvp(this)

  implicit none

  class(parallel), intent(in) :: this
  integer :: getnvp

  getnvp = this%nvp

end function getnvp
!
function getidproc(this)

  implicit none

  class(parallel), intent(in) :: this
  integer :: getidproc

  getidproc = this%idproc

end function getidproc
!
function getkstrt(this)

  implicit none

  class(parallel), intent(in) :: this
  integer :: getkstrt

  getkstrt = this%kstrt

end function getkstrt
!
function getlworld(this)

  implicit none

  class(parallel), intent(in) :: this
  integer :: getlworld

  getlworld = this%lworld

end function getlworld
!
function getmint(this)

  implicit none

  class(parallel), intent(in) :: this
  integer :: getmint

  getmint = this%mint

end function getmint
!
function getmreal(this)

  implicit none

  class(parallel), intent(in) :: this
  integer :: getmreal

  getmreal = this%mreal

end function getmreal
!
function getmdouble(this)

  implicit none

  class(parallel), intent(in) :: this
  integer :: getmdouble

  getmdouble = this%mdouble

end function getmdouble
!
function getmcplx(this)

  implicit none

  class(parallel), intent(in) :: this
  integer :: getmcplx

  getmcplx = this%mcplx

end function getmcplx
!
function getmchar(this)

  implicit none

  class(parallel), intent(in) :: this
  integer :: getmchar

  getmchar = this%mchar

end function getmchar
!
function getmint_pack(this)

  implicit none

  class(parallel), intent(in) :: this
  integer :: getmint_pack

  getmint_pack = this%mint_pack

end function getmint_pack
!
function getmreal_pack(this)

  implicit none

  class(parallel), intent(in) :: this
  integer :: getmreal_pack

  getmreal_pack = this%mreal_pack

end function getmreal_pack
!
subroutine init_parallel(this)

  implicit none

  class(parallel), intent(inout) :: this
  ! nvpp = number of shared memory threads (0=default)
  integer :: nvpp = 0

  ! initialize for shared memory parallel processing using openmp
  call init_omp(nvpp)

  ! initialize for distributed memory parallel processing using mpi
  call ppinit2(this%idproc,this%nvp,this%lworld,&
  &this%mint,this%mreal,this%mdouble,this%mcplx,this%mchar)
  this%kstrt = this%idproc + 1

end subroutine init_parallel
!
subroutine ppinit2(idproc,nvp,lworld,mint,mreal,mdouble,mcplx,mchar)
  ! this subroutine initializes parallel processing using mpi
  implicit none

  integer, intent(inout) :: idproc, nvp
  integer, intent(inout) :: lworld,mint,mreal,mdouble,mcplx,mchar
  ! nproc = number of real or virtual processors obtained
  ! mreal = default datatype for reals
  ! mint = default datatype for integers
  ! mcplx = default datatype for complex type
  ! mdouble = default double precision type
  ! mchar = default datatype for character type
  ! lworld = MPI_COMM_WORLD communicator
  ! local data
  integer :: ierror, ndprec, idprec
  integer :: iprec
  logical :: flag
  real :: prec
  integer :: nproc
  ! ndprec = (0,1) = (no,yes) use (normal,autodouble) precision
  if (digits(prec) > 24) then
     ndprec = 1
  else
     ndprec = 0
  endif
  ! idprec = (0,1) = (no,yes) use (normal,autodouble) integer precision
  if (digits(iprec) > 31) then
     idprec = 1
  else
     idprec = 0
  endif
  ! this segment is used for mpi computers
  ! indicate whether MPI_INIT has been called
  call MPI_INITIALIZED(flag,ierror)
  if (.not.flag) then
  ! initialize the MPI execution environment
     call MPI_INIT(ierror)
     if (ierror /= 0) stop
  endif
  lworld = MPI_COMM_WORLD
  ! determine the rank of the calling process in the communicator
  call MPI_COMM_RANK(lworld,idproc,ierror)
  ! determine the size of the group associated with a communicator
  call MPI_COMM_SIZE(lworld,nproc,ierror)
  ! set default datatypes
  mint = MPI_INTEGER
  mdouble = MPI_DOUBLE_PRECISION
  mchar = MPI_CHARACTER
  ! single precision real
  if (ndprec==0) then
     mreal = MPI_REAL
     mcplx = MPI_COMPLEX
  ! double precision real
  else
     mreal = MPI_DOUBLE_PRECISION
     mcplx = MPI_DOUBLE_COMPLEX
  endif

  ! call MPI_TYPE_CONTIGUOUS( huge(1), mint, mint_pack )
  ! call MPI_TYPE_CONTIGUOUS( huge(1), mreal, mreal_pack )

  nvp = nproc

end subroutine ppinit2
!
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

  !$OMP PARALLEL PRIVATE(tid)
  tid = OMP_GET_THREAD_NUM()
  !$OMP END PARALLEL

end subroutine init_omp
!
subroutine end_parallel(this)
  ! this subroutine terminates parallel processing
  implicit none

  class(parallel), intent(inout) :: this
  ! lworld = MPI_COMM_WORLD communicator
  ! local data
  integer :: ierror
  logical :: flag
  ! indicate whether MPI_INIT has been called
  call MPI_INITIALIZED(flag,ierror)
  if (flag) then
  ! synchronize processes
     call MPI_BARRIER(this%getlworld(),ierror)
  ! terminate MPI execution environment
     call MPI_FINALIZE(ierror)
  endif

end subroutine end_parallel
!
function ntag()

   implicit none
   integer, save :: tag = 0
   integer :: ntag

   ntag = tag
   tag = tag + 1

end function ntag

end module parallel_class

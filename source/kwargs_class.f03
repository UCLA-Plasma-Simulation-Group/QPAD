module kwargs_class

implicit none

private

interface kw_arg
  module procedure :: kw_arg_constructor
end interface kw_arg

type :: kw_arg
  private
  type( kw_arg ), pointer :: next => null()
  character(len=:), allocatable :: key
  class(*), allocatable :: value
  contains
  procedure :: get_key => get_key_kw_arg
  generic :: get_value => get_value_real_kw_arg, &
                          get_value_integer_kw_arg, &
                          get_value_logical_kw_arg, &
                          get_value_complex_kw_arg, &
                          get_value_character_kw_arg
  procedure, private :: get_value_real_kw_arg
  procedure, private :: get_value_integer_kw_arg
  procedure, private :: get_value_logical_kw_arg
  procedure, private :: get_value_complex_kw_arg
  procedure, private :: get_value_character_kw_arg
  generic :: set_value => set_value_real_kw_arg, &
                          set_value_integer_kw_arg, &
                          set_value_logical_kw_arg, &
                          set_value_complex_kw_arg, &
                          set_value_character_kw_arg
  procedure, private :: set_value_real_kw_arg
  procedure, private :: set_value_integer_kw_arg
  procedure, private :: set_value_logical_kw_arg
  procedure, private :: set_value_complex_kw_arg
  procedure, private :: set_value_character_kw_arg
  final :: kw_arg_finalizer
end type kw_arg

type :: kw_list
  private
  type( kw_arg ), pointer :: head => null(), tail => null()
  integer :: num_kw_arg = 0
  contains
  procedure :: append => append_kw_list
  generic   :: get => get_real_kw_list, &
                      get_integer_kw_list, &
                      get_complex_kw_list, &
                      get_logical_kw_list, &
                      get_character_kw_list
  procedure, private :: get_real_kw_list
  procedure, private :: get_integer_kw_list
  procedure, private :: get_complex_kw_list
  procedure, private :: get_logical_kw_list
  procedure, private :: get_character_kw_list
  generic   :: set => set_real_kw_list, &
                      set_integer_kw_list, &
                      set_complex_kw_list, &
                      set_logical_kw_list, &
                      set_character_kw_list
  procedure, private :: set_real_kw_list
  procedure, private :: set_integer_kw_list
  procedure, private :: set_complex_kw_list
  procedure, private :: set_logical_kw_list
  procedure, private :: set_character_kw_list
  final :: kw_list_finalizer
end type kw_list

public :: kw_list, kw_arg

integer, parameter :: KEY_MAX_LEN = 64

contains

function kw_arg_constructor( key, value ) result(res)
  implicit none
  character(len=*), intent(in), optional :: key
  class(*), intent(in), optional :: value
  type(kw_arg) :: res
  allocate( res%value, source=value )
  res%key = trim(key)
end function kw_arg_constructor

subroutine kw_arg_finalizer( this )
  implicit none
  type( kw_arg ), intent(inout) :: this
  if ( associated( this%next ) ) nullify( this%next )
  if ( allocated( this%key ) ) deallocate( this%key )
  if ( allocated( this%value ) ) deallocate( this%value )
end subroutine kw_arg_finalizer

subroutine get_key_kw_arg( this, key, ierr )

  implicit none
  class( kw_arg ), intent(in) :: this
  character(*), intent(out) :: key
  integer, intent(out) :: ierr

  if ( allocated( this%key ) ) then
    key = this%key
    ierr = 0
  else
    ierr = 1
  endif

end subroutine get_key_kw_arg

subroutine get_value_real_kw_arg( this, value, ierr )

  implicit none
  class( kw_arg ), intent(in) :: this
  real, intent(out) :: value
  integer, intent(out) :: ierr

  if ( allocated( this%value ) ) then
    select type ( val => this%value )
      type is (real)
        value = val
        ierr = 0
      class default
        ierr = 1
    end select
  else
    ierr = 2
  endif

end subroutine get_value_real_kw_arg

subroutine get_value_integer_kw_arg( this, value, ierr )

  implicit none
  class( kw_arg ), intent(in) :: this
  integer, intent(out) :: value
  integer, intent(out) :: ierr

  if ( allocated( this%value ) ) then
    select type ( val => this%value )
      type is (integer)
        value = val
        ierr = 0
      class default
        ierr = 1
    end select
  else
    ierr = 2
  endif

end subroutine get_value_integer_kw_arg

subroutine get_value_logical_kw_arg( this, value, ierr )
  implicit none
  class( kw_arg ), intent(in) :: this
  logical, intent(out) :: value
  integer, intent(out) :: ierr

  if ( allocated( this%value ) ) then
    select type ( val => this%value )
      type is (logical)
        value = val
        ierr = 0
      class default
        ierr = 1
    end select
  else
    ierr = 2
  endif

end subroutine get_value_logical_kw_arg

subroutine get_value_complex_kw_arg( this, value, ierr )
  implicit none
  class( kw_arg ), intent(in) :: this
  complex, intent(out) :: value
  integer, intent(out) :: ierr

  if ( allocated( this%value ) ) then
    select type ( val => this%value )
      type is (complex)
        value = val
        ierr = 0
      class default
        ierr = 1
    end select
  else
    ierr = 2
  endif

end subroutine get_value_complex_kw_arg

subroutine get_value_character_kw_arg( this, value, ierr )
  implicit none
  class( kw_arg ), intent(in) :: this
  character(*), intent(out) :: value
  integer, intent(out) :: ierr

  if ( allocated( this%value ) ) then
    select type ( val => this%value )
      type is (character(*))
        value = val
        ierr = 0
      class default
        ierr = 1
    end select
  else
    ierr = 2
  endif

end subroutine get_value_character_kw_arg

subroutine set_value_real_kw_arg( this, value, ierr )

  implicit none
  class( kw_arg ), intent(inout) :: this
  real, intent(in) :: value
  integer, intent(out) :: ierr

  if ( allocated( this%value ) ) then
    select type ( val => this%value )
      type is (real)
        val = value
        ierr = 0
      class default
        ierr = 1
    end select
  else
    ierr = 2
  endif

end subroutine set_value_real_kw_arg

subroutine set_value_integer_kw_arg( this, value, ierr )

  implicit none
  class( kw_arg ), intent(inout) :: this
  integer, intent(in) :: value
  integer, intent(out) :: ierr

  if ( allocated( this%value ) ) then
    select type ( val => this%value )
      type is (integer)
        val = value
        ierr = 0
      class default
        ierr = 1
    end select
  else
    ierr = 2
  endif

end subroutine set_value_integer_kw_arg

subroutine set_value_logical_kw_arg( this, value, ierr )
  implicit none
  class( kw_arg ), intent(inout) :: this
  logical, intent(in) :: value
  integer, intent(out) :: ierr

  if ( allocated( this%value ) ) then
    select type ( val => this%value )
      type is (logical)
        val = value
        ierr = 0
      class default
        ierr = 1
    end select
  else
    ierr = 2
  endif

end subroutine set_value_logical_kw_arg

subroutine set_value_complex_kw_arg( this, value, ierr )
  implicit none
  class( kw_arg ), intent(inout) :: this
  complex, intent(in) :: value
  integer, intent(out) :: ierr

  if ( allocated( this%value ) ) then
    select type ( val => this%value )
      type is (complex)
        val = value
        ierr = 0
      class default
        ierr = 1
    end select
  else
    ierr = 2
  endif

end subroutine set_value_complex_kw_arg

subroutine set_value_character_kw_arg( this, value, ierr )
  implicit none
  class( kw_arg ), intent(inout) :: this
  character(*), intent(in) :: value
  integer, intent(out) :: ierr

  if ( allocated( this%value ) ) then
    select type ( val => this%value )
      type is (character(*))
        val = value
        ierr = 0
      class default
        ierr = 1
    end select
  else
    ierr = 2
  endif

end subroutine set_value_character_kw_arg

subroutine append_kw_list( this, key, value )

  implicit none
  class( kw_list ), intent(inout) :: this
  character(len=*), intent(in) :: key
  class(*), intent(in) :: value

  if ( associated( this%tail ) ) then
    allocate( this%tail%next, source=kw_arg( key, value ) )
    this%tail => this%tail%next
  else
    allocate( this%head, source=kw_arg( key, value ) )
    this%tail => this%head
  endif
  this%num_kw_arg = this%num_kw_arg + 1

end subroutine append_kw_list

subroutine get_real_kw_list( this, key, value, err )

  implicit none
  class( kw_list ), intent(in) :: this
  character(len=*), intent(in) :: key
  real, intent(out) :: value
  integer, intent(out), optional :: err

  character(len=KEY_MAX_LEN) :: key_tmp
  type( kw_arg ), pointer :: current => null()
  integer :: ierr, err_tmp

  current => this%head
  do while ( associated(current) )
    call current%get_key( key_tmp, ierr )
    if ( trim(key_tmp) == trim(key) ) then
      call current%get_value( value, ierr )
      err_tmp = ierr
      current => null()
      return
    endif
    current => current%next
  enddo
  err_tmp = 3
  current => null()
  if ( present(err) ) err = err_tmp

end subroutine get_real_kw_list

subroutine get_integer_kw_list( this, key, value, err )

  implicit none
  class( kw_list ), intent(in) :: this
  character(len=*), intent(in) :: key
  integer, intent(out) :: value
  integer, intent(out), optional :: err

  character(len=KEY_MAX_LEN) :: key_tmp
  type( kw_arg ), pointer :: current => null()
  integer :: ierr, err_tmp

  current => this%head
  do while ( associated(current) )
    call current%get_key( key_tmp, ierr )
    if ( trim(key_tmp) == trim(key) ) then
      call current%get_value( value, ierr )
      err_tmp = ierr
      current => null()
      return
    endif
    current => current%next
  enddo
  err_tmp = 3
  current => null()
  if ( present(err) ) err = err_tmp

end subroutine get_integer_kw_list

subroutine get_logical_kw_list( this, key, value, err )

  implicit none
  class( kw_list ), intent(in) :: this
  character(len=*), intent(in) :: key
  logical, intent(out) :: value
  integer, intent(out), optional :: err

  character(len=KEY_MAX_LEN) :: key_tmp
  type( kw_arg ), pointer :: current => null()
  integer :: ierr, err_tmp

  current => this%head
  do while ( associated(current) )
    call current%get_key( key_tmp, ierr )
    if ( trim(key_tmp) == trim(key) ) then
      call current%get_value( value, ierr )
      err_tmp = ierr
      current => null()
      return
    endif
    current => current%next
  enddo
  err_tmp = 3
  current => null()
  if ( present(err) ) err = err_tmp

end subroutine get_logical_kw_list

subroutine get_complex_kw_list( this, key, value, err )

  implicit none
  class( kw_list ), intent(in) :: this
  character(len=*), intent(in) :: key
  complex, intent(out) :: value
  integer, intent(out), optional :: err

  character(len=KEY_MAX_LEN) :: key_tmp
  type( kw_arg ), pointer :: current => null()
  integer :: ierr, err_tmp

  current => this%head
  do while ( associated(current) )
    call current%get_key( key_tmp, ierr )
    if ( trim(key_tmp) == trim(key) ) then
      call current%get_value( value, ierr )
      err_tmp = ierr
      current => null()
      return
    endif
    current => current%next
  enddo
  err_tmp = 3
  current => null()
  if ( present(err) ) err = err_tmp

end subroutine get_complex_kw_list

subroutine get_character_kw_list( this, key, value, err )

  implicit none
  class( kw_list ), intent(in) :: this
  character(len=*), intent(in) :: key
  character(*), intent(out) :: value
  integer, intent(out), optional :: err

  character(len=KEY_MAX_LEN) :: key_tmp
  type( kw_arg ), pointer :: current => null()
  integer :: ierr, err_tmp

  current => this%head
  do while ( associated(current) )
    call current%get_key( key_tmp, ierr )
    if ( trim(key_tmp) == trim(key) ) then
      call current%get_value( value, ierr )
      err_tmp = ierr
      current => null()
      return
    endif
    current => current%next
  enddo
  err_tmp = 3
  current => null()
  if ( present(err) ) err = err_tmp

end subroutine get_character_kw_list

subroutine set_real_kw_list( this, key, value, err )

  implicit none
  class( kw_list ), intent(in) :: this
  character(len=*), intent(in) :: key
  real, intent(in) :: value
  integer, intent(out), optional :: err

  character(len=KEY_MAX_LEN) :: key_tmp
  type( kw_arg ), pointer :: current => null()
  integer :: ierr, err_tmp

  current => this%head
  do while ( associated(current) )
    call current%get_key( key_tmp, ierr )
    if ( trim(key_tmp) == trim(key) ) then
      call current%set_value( value, ierr )
      err_tmp = ierr
      current => null()
      return
    endif
    current => current%next
  enddo
  err_tmp = 3
  current => null()
  if ( present(err) ) err = err_tmp

end subroutine set_real_kw_list

subroutine set_integer_kw_list( this, key, value, err )

  implicit none
  class( kw_list ), intent(in) :: this
  character(len=*), intent(in) :: key
  integer, intent(in) :: value
  integer, intent(out), optional :: err

  character(len=KEY_MAX_LEN) :: key_tmp
  type( kw_arg ), pointer :: current => null()
  integer :: ierr, err_tmp

  current => this%head
  do while ( associated(current) )
    call current%get_key( key_tmp, ierr )
    if ( trim(key_tmp) == trim(key) ) then
      call current%set_value( value, ierr )
      err_tmp = ierr
      current => null()
      return
    endif
    current => current%next
  enddo
  err_tmp = 3
  current => null()
  if ( present(err) ) err = err_tmp

end subroutine set_integer_kw_list

subroutine set_logical_kw_list( this, key, value, err )

  implicit none
  class( kw_list ), intent(in) :: this
  character(len=*), intent(in) :: key
  logical, intent(in) :: value
  integer, intent(out), optional :: err

  character(len=KEY_MAX_LEN) :: key_tmp
  type( kw_arg ), pointer :: current => null()
  integer :: ierr, err_tmp

  current => this%head
  do while ( associated(current) )
    call current%get_key( key_tmp, ierr )
    if ( trim(key_tmp) == trim(key) ) then
      call current%set_value( value, ierr )
      err_tmp = ierr
      current => null()
      return
    endif
    current => current%next
  enddo
  err_tmp = 3
  current => null()
  if ( present(err) ) err = err_tmp

end subroutine set_logical_kw_list

subroutine set_complex_kw_list( this, key, value, err )

  implicit none
  class( kw_list ), intent(in) :: this
  character(len=*), intent(in) :: key
  complex, intent(in) :: value
  integer, intent(out), optional :: err

  character(len=KEY_MAX_LEN) :: key_tmp
  type( kw_arg ), pointer :: current => null()
  integer :: ierr, err_tmp

  current => this%head
  do while ( associated(current) )
    call current%get_key( key_tmp, ierr )
    if ( trim(key_tmp) == trim(key) ) then
      call current%set_value( value, ierr )
      err_tmp = ierr
      current => null()
      return
    endif
    current => current%next
  enddo
  err_tmp = 3
  current => null()
  if ( present(err) ) err = err_tmp

end subroutine set_complex_kw_list

subroutine set_character_kw_list( this, key, value, err )

  implicit none
  class( kw_list ), intent(in) :: this
  character(len=*), intent(in) :: key
  character(*), intent(in) :: value
  integer, intent(out), optional :: err

  character(len=KEY_MAX_LEN) :: key_tmp
  type( kw_arg ), pointer :: current => null()
  integer :: ierr, err_tmp

  current => this%head
  do while ( associated(current) )
    call current%get_key( key_tmp, ierr )
    if ( trim(key_tmp) == trim(key) ) then
      call current%set_value( value, ierr )
      err_tmp = ierr
      current => null()
      return
    endif
    current => current%next
  enddo
  err_tmp = 3
  current => null()
  if ( present(err) ) err = err_tmp

end subroutine set_character_kw_list

subroutine kw_list_finalizer( this )
  implicit none
  type(kw_list), intent(inout) :: this
  this%num_kw_arg = 0
  if ( associated(this%head) ) deallocate(this%head)
end subroutine kw_list_finalizer

end module kwargs_class
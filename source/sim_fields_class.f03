module sim_fields_class

use parallel_pipe_class
use grid_class
use field_psi_class
use field_e_class
use field_b_class
use field_src_class

use sys
use param

implicit none

private

public :: sim_fields

type sim_fields

  ! private

  class( parallel_pipe ), pointer :: pp => null()
  class( grid ), pointer :: gp => null()

  type( field_psi ), allocatable :: psi
  type( field_b ), allocatable :: b_spe, b_beam
  type( field_e ), allocatable :: e_spe, e_beam
  type( field_jay ), allocatable :: jay
  type( field_rho ), allocatable :: q_spe, q_beam
  type( field_djdxi ), allocatable :: djdxi

  contains

  generic :: new => init_sim_fields
  generic :: del => end_sim_fields

  procedure, private :: init_sim_fields, end_sim_fields

end type sim_fields

character(len=18), save :: cls_name = 'sim_fields'
integer, save :: cls_level = 2

contains

subroutine init_sim_fields( this, pp, gp, dr, dxi, num_modes, part_shape, boundary, iter_tol )

  implicit none

  class( sim_fields ), intent(inout) :: this
  class( parallel_pipe ), intent(in), pointer :: pp
  class( grid ), intent(in), pointer :: gp
  real, intent(in) :: dr, dxi, iter_tol
  integer, intent(in) :: num_modes, part_shape, boundary

  ! local data
  character(len=18), save :: sname = 'init_sim_fields'
  character(len=20) :: s1, s2, s3
  character(len=:), allocatable :: ff
  integer :: i,n,ndump,j,k,l,m
  integer :: entity

  this%gp => gp
  this%pp => pp

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  allocate( this%psi, this%e_spe, this%b_spe, this%e_beam, this%b_beam, &
    this%jay, this%q_spe, this%q_beam, this%djdxi )

  call this%psi%new( this%pp, this%gp, dr, dxi, num_modes, part_shape, boundary, iter_tol )
  call this%jay%new( this%pp, this%gp, dr, dxi, num_modes, part_shape )
  call this%q_spe%new( this%pp, this%gp, dr, dxi, num_modes, part_shape )
  call this%q_beam%new( this%pp, this%gp, dr, dxi, num_modes, part_shape )
  call this%djdxi%new( this%pp, this%gp, dr, dxi, num_modes, part_shape )
  entity = p_entity_plasma
  call this%e_spe%new( this%pp, this%gp, dr, dxi, num_modes, part_shape, boundary, entity, iter_tol )
  call this%b_spe%new( this%pp, this%gp, dr, dxi, num_modes, part_shape, boundary, entity, iter_tol )
  entity = p_entity_beam
  call this%e_beam%new( this%pp, this%gp, dr, dxi, num_modes, part_shape, boundary, entity, iter_tol )
  call this%b_beam%new( this%pp, this%gp, dr, dxi, num_modes, part_shape, boundary, entity, iter_tol )

  ! call input%get('simulation.nspecies',n)

  ! loop1: do i = 1, n
  !   write (s1, '(I4.4)') i
  !   call input%info('species('//trim(s1)//').diag',n_children=m)
  !   do j = 1, m
  !      write (s2, '(I4.4)') j
  !      call input%get('species('//trim(s1)//').diag'//'('//trim(s2)//').ndump',ndump)
  !      if (ndump>0) then
  !         call input%info('species('//trim(s1)//').diag'//'('//trim(s2)//').name',n_children=l)
  !         do k = 1, l
  !            write (s3, '(I4.4)') k
  !            if(allocated(ff)) deallocate(ff)
  !            call input%get('species('//trim(s1)//').diag'//'('//trim(s2)//').name'&
  !            &//'('//trim(s3)//')',ff)
  !            if (ff == 'jx' .or. ff == 'jy' .or. ff == 'jz') then
  !               allocate(this%cu3d)
  !               call this%cu3d%new(this%p,this%err,this%sp3,dim=1)
  !               exit loop1
  !            end if
  !         end do
  !      end if
  !   end do
  ! end do loop1

  ! call input%info('field.diag',n_children=n)

  ! loop2: do i = 1, n
  !   write (s1,'(I4.4)') i
  !   call input%get('field.diag('//trim(s1)//').ndump',ndump)
  !   if (ndump > 0) then
  !      call input%info('field.diag('//trim(s1)//').name',n_children=m)
  !      do j = 1, m
  !         write (s2,'(I4.4)') j
  !         if(allocated(ff)) deallocate(ff)
  !         call input%get('field.diag('//trim(s1)//').name('//trim(s2)//')',ff)
  !         if (ff == 'psi') then
  !            allocate(this%psi3d)
  !            call this%psi3d%new(this%p,this%err,this%sp3,dim=1)
  !            exit loop2
  !         end if
  !      end do
  !   end if
  ! end do loop2

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine init_sim_fields

subroutine end_sim_fields( this )

  implicit none

  class(sim_fields), intent(inout) :: this
  ! local data
  character(len=18), save :: sname = 'end_sim_fields'
  integer :: i, n

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  call this%psi%del()
  call this%e_spe%del()
  call this%b_spe%del()
  call this%e_beam%del()
  call this%b_beam%del()
  call this%jay%del()
  call this%q_spe%del()
  call this%q_beam%del()
  call this%djdxi%del()

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine end_sim_fields

end module sim_fields_class
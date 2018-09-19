module param

implicit none

public

! ================================================================
! mathematics
! ================================================================
real, parameter :: pi = 3.14159265359
integer, parameter :: p_real = 0, p_imag = 1

! ================================================================
! dimension related
! ================================================================
integer, parameter :: p_max_xdim = 3
integer, parameter :: p_lower = 1, p_upper = 2

! ================================================================
! particle shape
! ================================================================
integer, parameter :: p_ps_linear = 1, p_ps_quadratic = 2

! ================================================================
! electromagnetic field solver
! ================================================================
! field solver order
integer, parameter :: p_fs_2order = 1, p_fs_4order = 2
! boundary condition
integer, parameter :: p_bnd_axial = 0, p_bnd_conduct = 1
! field kind
integer, parameter :: p_fk_psi = 0, p_fk_ez = 1, p_fk_bz = 2, &
                      p_fk_bperp = 3, p_fk_br_iter = 4, p_fk_bphi_iter = 5

! ================================================================
! HYPRE parameters
! ================================================================
integer, parameter :: HYPRE_TYPE = 8
integer, parameter :: p_hypre_cycred = 1, p_hypre_pcg = 2

! ================================================================
! copy field, used by ufield class
! ================================================================
integer, parameter :: p_copy_1to2 = 1, p_copy_2to1 = -1

end module param
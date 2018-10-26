module param

implicit none

public

! ================================================================
! data kind values
! ================================================================

integer, parameter :: SG = kind(0.0)
integer, parameter :: DB = kind(0.0d0)
integer, parameter :: LG = selected_int_kind(12)

! ================================================================
! mathematics
! ================================================================
real, parameter :: pi = 4*atan(1.0_DB)
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
                      p_fk_bperp = 3, p_fk_br_iter = 4, p_fk_bphi_iter = 5, &
                      p_fk_bperp_iter = 6

! ================================================================
! entity of electromagnetic field
! ================================================================
integer, parameter :: p_entity_none = 0, p_entity_beam = 1, p_entity_plasma = 2

! ================================================================
! HYPRE parameters
! ================================================================
integer, parameter :: HYPRE_TYPE = 8, HYPRE_PARCSR = 5555
integer, parameter :: p_hypre_cycred = 1, p_hypre_pcg = 2, p_hypre_smg = 3, &
                      p_hypre_gmres = 4, p_hypre_amg = 5, p_hypre_parpcg = 6

! ================================================================
! copy field, used by ufield class
! ================================================================
integer, parameter :: p_copy_1to2 = 1, p_copy_2to1 = -1

! ================================================================
! parameters for copy guard cells
! ================================================================
integer, parameter :: p_mpi_forward = 1, p_mpi_backward = -1, p_mpi_bothway = 0

end module param
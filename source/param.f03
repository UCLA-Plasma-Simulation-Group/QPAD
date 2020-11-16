module param

implicit none

public

! ================================================================
! hardware related
! ================================================================
integer, parameter :: p_cache_size = 1024

! ================================================================
! simulation module
! ================================================================
integer, parameter :: p_sim_standard = 0, &
                      p_sim_popas = 1, &
                      p_sim_tmplt = 99

! ================================================================
! geometry related
! ================================================================
integer, parameter :: p_cylindrical = 0, p_cartesian = 1

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
integer, parameter :: p_x_dim = 3, p_p_dim = 3, p_s_dim = 3

! ================================================================
! particle shape
! ================================================================
integer, parameter :: p_ps_linear = 1, p_ps_quadratic = 2

! ================================================================
! particle pusher
! ================================================================
integer, parameter :: p_push_reduced = 1, p_push_boris = 2

! ================================================================
! electromagnetic field solver
! ================================================================
! field solver order
integer, parameter :: p_fs_2order = 1, p_fs_4order = 2
! boundary condition
integer, parameter :: p_bnd_axial = 0, &
                      p_bnd_zero  = 2, &
                      p_bnd_open  = 3
! field kind
integer, parameter :: p_fk_psi     = 0, &
                      p_fk_ez      = 1, &
                      p_fk_bz      = 2, &
                      p_fk_bt      = 3, &
                      p_fk_bplus   = 4, &
                      p_fk_bminus  = 5, &
                      p_fk_vpotz   = 6, &
                      p_fk_vpotp   = 7, &
                      p_fk_vpotm   = 8

! ================================================================
! entity of electromagnetic field
! ================================================================
integer, parameter :: p_entity_none       = 0, &
                      p_entity_beam       = 1, &
                      p_entity_plasma     = 2

! ================================================================
! HYPRE parameters
! ================================================================
integer, parameter :: HYPRE_TYPE = 8, HYPRE_PARCSR = 5555
! TODO: now all the solvers are using cyclic reduction
! other options may be deleted in the future
integer, parameter :: p_hypre_cycred = 1, &
                      p_hypre_pcg    = 2, &
                      p_hypre_smg    = 3, &
                      p_hypre_gmres  = 4, &
                      p_hypre_amg    = 5, &
                      p_hypre_parpcg = 6

! ================================================================
! copy field, used by ufield class
! ================================================================
integer, parameter :: p_copy_1to2 = 1, p_copy_2to1 = -1

! ================================================================
! parameters for data communication
! ================================================================
integer, parameter :: p_mpi_forward  = 1, &
                      p_mpi_backward = -1, &
                      p_mpi_bothway  = 0

! ================================================================
! parameters for cell types
! ================================================================
integer, parameter :: p_cell_guard = 0, &
                      p_cell_inner = 1

! ================================================================
! parameters for smooth
! ================================================================
integer, parameter :: p_smooth_none        = 0, &
                      p_smooth_binomial    = 1, &
                      p_smooth_compensated = 2

! ================================================================
! parameters for diagnostic type
! ================================================================
integer, parameter :: p_tdiag_grid = 0, &
                      p_tdiag_raw  = 1, &
                      p_tdiag_rst  = 2

! ================================================================
! others
! ================================================================
integer, parameter :: MAX_LEN_LABEL = 32

end module param
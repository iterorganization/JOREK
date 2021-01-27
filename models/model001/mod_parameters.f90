module mod_parameters

  implicit none

  integer, parameter :: jorek_model    = 001       !< JOREK physics model

  integer, parameter :: var_A3   = 0                       ! place of variable psi/mag pot 3               (ps or A3)
  integer, parameter :: var_AR   = 0                       ! place of variable mag pot  1                  (AR)
  integer, parameter :: var_AZ   = 0                       ! place of variable mag pot  2                  (AZ)
  integer, parameter :: var_uR   = 0                       ! place of variable velocity 1                  (UR)
  integer, parameter :: var_uZ   = 0                       ! place of variable velocity 2                  (UZ)
  integer, parameter :: var_up   = 0                       ! place of variable velocity 3                  (Up)
  integer, parameter :: var_rho  = 0                       ! place of variable density                     (r or rho)
  integer, parameter :: var_T    = 0                       ! place of variable temperature                 (T)
  integer, parameter :: var_psi  = 0                       ! place of variable psi/mag pot 3               (ps or A3)  
  integer, parameter :: var_u    = 1                       ! place of variable velocity stream function    (u)
  integer, parameter :: var_zj   = 0                       ! place of variable toroidal current density    (zj)
  integer, parameter :: var_w    = 2                       ! place of variable vorticity                   (w)
  integer, parameter :: var_Vpar = 0                       ! place of variable parallel velocity           (Vpar)
  integer, parameter :: var_rhon = 0                       ! place of variable neutral or impurity density (rn)
  integer, parameter :: var_Ti   = 0                       ! place of variable ion temperature             (Ti)
  integer, parameter :: var_Te   = 0                       ! place of variable electron temperature        (Te)
  integer, parameter :: var_jec  = 0                       ! place of variable ECCD current                (jec)
  integer, parameter :: var_jec1 = 0                       ! place of variable ECCD current #1             (jec1)
  integer, parameter :: var_jec2 = 0                       ! place of variable ECCD current #2             (jec2)

  integer, parameter :: n_var          = 2         !< number of variables
  integer, parameter :: n_dim          = 2         !< number of dimensions
  integer, parameter :: n_order        = 3         !< order of the polynomial basis
  integer, parameter :: n_tor          = 1         !< number of toroidal harmonics in physics variables
  integer, parameter :: n_coord_tor    = 1         !< number of toroidal harmonics in (R, Z) coordinates
  integer, parameter :: n_period       = 1         !< periodicity in toroidal direction
  integer, parameter :: n_plane        = 1         !< number of toroidal angles
  integer, parameter :: n_vertex_max   = 4         !< maximum number of corners of an element
  integer, parameter :: n_nodes_max    = 80001     !< maximum number of nodes
  integer, parameter :: n_elements_max = 80001     !< maximum number of elements
  integer, parameter :: n_boundary_max = 1001      !< maximum number of boundary elements
  integer, parameter :: n_pieces_max   = 6001      !< maximum number of line pieces describing a flux surface
  integer, parameter :: n_degrees      = n_order+1 !< degrees of freedom per variable per node
  integer, parameter :: nref_max       = 10000     !< (refinement)
  integer, parameter :: n_ref_list     = 10000     !< (refinement)
  character*11 :: variable_names(n_var)
 
  !> Names of the physical variables
  parameter (variable_names = (/ 'Potential  ','Vorticity  ' /))
  
  !> element_matrix and element_matrix_fft combined into a single one?
  logical, parameter :: unified_element_matrix = .false.

endmodule

module mod_settings

implicit none

integer, parameter :: n_tor             = 3      !< number of toroidal harmonics
integer, parameter :: n_coord_tor       = 1      !< number of toroidal harmonics in (R, Z) coordinates
integer, parameter :: l_pol_domm        = 0      !< highest poloidal mode in the Dommaschk potentials        
integer, parameter :: n_period          = 1      !< periodicity in toroidal direction
integer, parameter :: n_coord_period    = 1      !< periodicity of the device in toroidal direction: equivalent to number of field periods
integer, parameter :: n_plane           = 4      !< number of toroidal angles
integer, parameter :: n_order           = 3      !< order of the polynomial basis
integer, parameter :: n_nodes_max       = 60001  !< maximum number of nodes
integer, parameter :: n_elements_max    = 60001  !< maximum number of elements
integer, parameter :: n_boundary_max    = 1001   !< maximum number of boundary elements
integer, parameter :: n_pieces_max      = 6001   !< maximum number of line pieces describing a flux surface

! ##################################################################################################
! ####  @USERS: Please do not change below this line ###############################################
! ##################################################################################################

! The following line is needed by ./util/config.sh:
! #SETTINGS# n_tor n_coord_tor l_pol_domm n_period n_coord_period n_plane n_order n_nodes_max n_elements_max n_boundary_max n_pieces_max

! --- a few constants that should not be touched
integer, parameter :: n_dim             = 2                !< number of dimensions
integer, parameter :: n_vertex_max      = 4                !< maximum number of corners of an element
integer, parameter :: n_degrees_1d      = (n_order+1)/2    !< degrees of freedom per variable per node in 1D (used for boundary conditions)
integer, parameter :: n_degrees         = n_degrees_1d**2  !< degrees of freedom per variable per node in 2D
integer, parameter :: nref_max          = 10               !< (refinement; not presently working)
integer, parameter :: n_ref_list        = 10               !< (refinement; not presently working)

end module mod_settings

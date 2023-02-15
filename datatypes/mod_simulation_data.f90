!> Definitions of derived data types for grid nodes and elements, boundary nodes and elements,
!! and flux surface elements, as well as the shattered pellets
module mod_simulation_data

  use equil_info,       only: t_equil_state
  use data_structure,   only: type_node_list, type_element_list, type_bnd_node_list, type_bnd_element_list
  use equil_info,       only: t_equil_state
  use mod_integer_types
  
  implicit none
  
  private
  public type_MHD_SIM

!> MHD simulation type, containing all variables pertaining to a simulation.
  type :: type_MHD_SIM
    real(kind=8)                                  :: time = 0.d0 !< time of the simulation
    type(type_node_list), pointer                 :: node_list => null() !< Current node list
    type(type_element_list), pointer              :: element_list => null() !< Current element list
    type(type_bnd_node_list), pointer             :: bnd_node_list => null()
    type(type_bnd_element_list), pointer          :: bnd_elm_list => null()
    type(t_equil_state), pointer                  :: es => null() !< Pointer to the equilibrium state
    
    integer, dimension(:), pointer                :: local_elms => null()
    integer                                       :: n_local_elms

    !< MPI settings
    integer                                       :: my_id = 0
    integer                                       :: n_cpu = 1 ! if not initialized, act as if there is no mpi
    
    integer                                       :: n_tor = 1
    
    logical                                       :: freeboundary
    logical                                       :: restart
    
    integer                                       :: sr_n_tor !< to pass sr%n_tor for direct construction
  end type type_MHD_SIM
  
  contains
  
end module mod_simulation_data  
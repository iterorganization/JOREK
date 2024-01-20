!>  Module that contais functions to fill IDSs from JOREK data
module mod_jorek2IMAS

#ifdef USE_IMAS
  use ids_schemas !, only: ids_equilibrium
  use ids_routines, only: imas_open_env, &
     imas_create_env, imas_close, ids_get, ids_put, ids_put_slice

  use mod_parameters 
  use data_structure
  use nodes_elements
  use constants
  
  implicit none
  
 
  public
 
  contains


  subroutine fill_mhd_IDS(first_step, idx)  

    use phys_module, only : t_start, F0, central_density, sqrt_mu0_rho0, &
                           sqrt_mu0_over_rho0, central_mass

    implicit none

    ! --- External parameters
    logical,            intent(in) :: first_step   ! is this the first step?
    integer,            intent(in) :: idx          ! IMAS identifier

   
    ! --- Local parameters 
    integer    :: i, j, k, m, etype, irst, int, i_var, i_tor, index, index_node, my_id, ierr
    real*8     :: fact_T, fact_time, fact_v, fact_zj, fact_psi, rho0, fact_phi, fact_rho, fact_w
    
    
    ! **********************************************************************************
    ! ******************************* IMAS **********************************************
    ! **********************************************************************************
    type(ids_mhd),                      target  :: mhd_ids
    type(ids_generic_grid_scalar),      pointer :: ggd_scalar
    type(ids_generic_grid_aos3_root),   pointer :: grid
    
    integer:: num_nodes, stat
    
    integer :: n_slice, i_slice, grid_ind, grid_sub_ind, n_grid_sub, n_grid
    ! **********************************************************************************
  
    ! --- Number of grids and grid subsets
    n_grid       = 1
    n_grid_sub   = 1
    grid_ind     = 1  ! Index
    grid_sub_ind = 1  ! Index
  
    if (first_step) then
      ! --- Put the grid in GGD
      allocate( mhd_ids%grid_ggd(n_grid) )
      grid => mhd_ids%grid_ggd(grid_ind)
      call grid2ggd( grid, node_list, element_list )
    endif

    ! --- Normalization factors for IMAS
    rho0               = central_density * 1.d20 * central_mass * mass_proton
    sqrt_mu0_rho0      = sqrt( mu_zero * rho0 )
    sqrt_mu0_over_rho0 = sqrt( mu_zero / rho0 )
  
    fact_psi  = -2.d0 * PI                  ! Transform to COCOS convention 8 --> 11
    fact_time =  sqrt_mu0_rho0 
    fact_v    =  1.d0 /  sqrt_mu0_rho0 
    fact_w    = -1.d0 /  sqrt_mu0_rho0      ! Transform for COCOS convention of toroidal direction (anti-clockwise) 
    fact_phi  = -1.d0 /  sqrt_mu0_rho0 * F0 ! COCOS convection: F0 depends on phi direction
    fact_zj   = -1.d0 / mu_zero * (-1.d0)   ! Last sign due to COCOS transformation
    fact_rho  =  rho0 
    fact_T    =  1.d0 / ( EL_CHG * mu_zero * central_density * 1.d20 )   
  
  
    ! --- Set times
    n_slice = 1  
    i_slice = 1
    allocate(  mhd_ids%time(n_slice) )
    allocate(  mhd_ids%ggd(n_slice ) )

    mhd_ids%ids_properties%homogeneous_time = 1
  
    mhd_ids%time(i_slice)     = t_start * fact_time 
    mhd_ids%ggd(i_slice)%time = t_start * fact_time

    ! --- Fill MHD data
    do i=1, n_var 
  
      ! --- Poloidal magnetic flux
      if (variable_names(i) == 'Psi') then      
        allocate( mhd_ids%ggd(i_slice)%psi(n_grid_sub))
        ggd_scalar => mhd_ids%ggd(i_slice)%psi(grid_sub_ind)
        call fill_Bezier_coefficients( ggd_scalar, node_list, var_psi, grid_ind, grid_sub_ind, fact_psi )
      endif
  
      ! --- Electrostatic potential 
      if (variable_names(i) == 'u') then      
        allocate( mhd_ids%ggd(i_slice)%phi_potential(n_grid_sub))
        ggd_scalar => mhd_ids%ggd(i_slice)%phi_potential(grid_sub_ind)
        call fill_Bezier_coefficients( ggd_scalar, node_list, var_u, grid_ind, grid_sub_ind, fact_phi )
      endif
  
      ! --- Toroidal current density * R
      if (variable_names(i) == 'zj') then      
        allocate( mhd_ids%ggd(i_slice)%j_tor_r(n_grid_sub))
        ggd_scalar => mhd_ids%ggd(i_slice)%j_tor_r(grid_sub_ind)
        call fill_Bezier_coefficients( ggd_scalar, node_list, var_zj, grid_ind, grid_sub_ind, fact_zj )
      endif
  
      ! --- Toroidal vorticity / R 
      if (variable_names(i) == 'omega') then      
        allocate( mhd_ids%ggd(i_slice)%vorticity_over_r(n_grid_sub))
        ggd_scalar => mhd_ids%ggd(i_slice)%vorticity_over_r(grid_sub_ind)
        call fill_Bezier_coefficients( ggd_scalar, node_list, var_w, grid_ind, grid_sub_ind, fact_w )
      endif
  
      ! --- Mass density
      if (variable_names(i) == 'rho') then      
        allocate( mhd_ids%ggd(i_slice)%mass_density(n_grid_sub))
        ggd_scalar => mhd_ids%ggd(i_slice)%mass_density(grid_sub_ind)
        call fill_Bezier_coefficients( ggd_scalar, node_list, var_rho, grid_ind, grid_sub_ind, fact_rho )
      endif
  
      ! --- Total temperature 
      if (variable_names(i) == 'T') then      
        ! --- Te
        allocate( mhd_ids%ggd(i_slice)%electrons%temperature(n_grid_sub))
        ggd_scalar => mhd_ids%ggd(i_slice)%electrons%temperature(grid_sub_ind)
        call fill_Bezier_coefficients( ggd_scalar, node_list, var_T, grid_ind, grid_sub_ind, fact_T*0.5d0 )
  
        ! --- Ti
        allocate( mhd_ids%ggd(i_slice)%t_i_average(n_grid_sub))
        ggd_scalar => mhd_ids%ggd(i_slice)%t_i_average(grid_sub_ind)
        call fill_Bezier_coefficients( ggd_scalar, node_list, var_T, grid_ind, grid_sub_ind, fact_T*0.5d0 )
      endif
  
      ! --- Ion temperature
      if (variable_names(i) == 'T_i') then      
        allocate( mhd_ids%ggd(i_slice)%t_i_average(n_grid_sub))
        ggd_scalar => mhd_ids%ggd(i_slice)%t_i_average(grid_sub_ind)
        call fill_Bezier_coefficients( ggd_scalar, node_list, var_Ti, grid_ind, grid_sub_ind, fact_T )
      endif
  
      ! --- Electron temperature
      if (variable_names(i) == 'T_e') then      
        allocate( mhd_ids%ggd(i_slice)%electrons%temperature(n_grid_sub))
        ggd_scalar => mhd_ids%ggd(i_slice)%electrons%temperature(grid_sub_ind)
        call fill_Bezier_coefficients( ggd_scalar, node_list, var_Te, grid_ind, grid_sub_ind, fact_T )
      endif
  
      ! --- Parallel velocity
      if (variable_names(i) == 'v_par') then      
        allocate( mhd_ids%ggd(i_slice)%velocity_parallel_over_b_field(n_grid_sub))
        ggd_scalar => mhd_ids%ggd(i_slice)%velocity_parallel_over_b_field(grid_sub_ind)
        call fill_Bezier_coefficients( ggd_scalar, node_list, var_vpar, grid_ind, grid_sub_ind, fact_v )
      endif
  
    enddo
  
    ! --- Put data into local database
    if (first_step) then  
      call ids_put(idx,'mhd',mhd_ids,stat)
    else
      call ids_put_slice(idx,'mhd',mhd_ids,stat)
    endif

    if (stat==0) then
       write(*,*) '    MHD IDS exported'
    else
       write(*,*) '    Something went wrong writting the MHD IDS!'
    endif

  end subroutine fill_mhd_IDS

 






  subroutine fill_radiation_IDS(first_step, idx)  

    use phys_module, only : t_start, F0, central_density, sqrt_mu0_rho0, &
                           sqrt_mu0_over_rho0, central_mass, imp_type, &
                           gamma, index_main_imp
    implicit none

    ! --- External parameters
    logical,      intent(in) :: first_step   ! is this the first step?
    integer,      intent(in) :: idx          ! IMAS identifier
   
    ! --- Local parameters 
    integer    :: i, j, k, m, var_rad, i_var, i_tor, index, index_node, my_id, ierr
    real*8     :: fact_time, rho0, fact_rad
    
    
    ! **********************************************************************************
    ! ******************************* IMAS **********************************************
    ! **********************************************************************************
    type(ids_radiation),                target  :: radiation_ids
    type(ids_generic_grid_scalar),      pointer :: ggd_scalar
    type(ids_generic_grid_aos3_root),   pointer :: grid
    
    integer:: num_nodes, stat
    
    integer :: n_slice, i_slice, grid_ind, grid_sub_ind, n_grid_sub, n_grid
    ! **********************************************************************************
  
    ! --- Number of grids and grid subsets
    n_grid       = 1
    n_grid_sub   = 1
    grid_ind     = 1  ! Index
    grid_sub_ind = 1  ! Index
  
    if (first_step) then
      ! --- Put the grid in GGD
      allocate( radiation_ids%grid_ggd(n_grid) )
      grid => radiation_ids%grid_ggd(grid_ind)
      call grid2ggd( grid, node_list, element_list )
    endif
 
    ! --- Normalization factors for IMAS
    rho0               = central_density * 1.d20 * central_mass * mass_proton
    sqrt_mu0_rho0      = sqrt( mu_zero * rho0 )
    sqrt_mu0_over_rho0 = sqrt( mu_zero / rho0 )

    fact_rad = 1.d0 / ( (gamma-1.d0) * MU_ZERO * sqrt_mu0_rho0 )
    fact_time =  sqrt_mu0_rho0 

    ! --- Set times
    n_slice = 1  
    i_slice = 1
    allocate(  radiation_ids%time(n_slice) )

    radiation_ids%ids_properties%homogeneous_time = 1
    allocate( radiation_ids%process(1))   ! --- 1 type of radiation
    allocate( radiation_ids%process(1)%ggd(n_slice) )
  
    radiation_ids%time(i_slice)                = t_start * fact_time 
    radiation_ids%process(1)%ggd(i_slice)%time = t_start * fact_time


 
    ! --- Fill radiation data 
    var_rad = 2
  
    allocate( radiation_ids%process(1)%ggd(i_slice)%ion(1))
    allocate( radiation_ids%process(1)%ggd(i_slice)%ion(1)%emissivity(n_grid_sub))
    allocate( radiation_ids%process(1)%ggd(i_slice)%ion(1)%label(1) )  
    allocate( radiation_ids%process(1)%identifier%name(1) )

    radiation_ids%process(1)%identifier%name  = "Line radiation"
    radiation_ids%process(1)%identifier%index = 10

    radiation_ids%process(1)%ggd(i_slice)%ion(1)%label = imp_type(index_main_imp) 
  
    ggd_scalar => radiation_ids%process(1)%ggd(i_slice)%ion(1)%emissivity(grid_sub_ind)
    call fill_Bezier_coefficients( ggd_scalar, aux_node_list, var_rad, grid_ind, grid_sub_ind, fact_rad )

    ! --- Put data into local database
    if (first_step) then  
      call ids_put(idx,'radiation',radiation_ids,stat)
    else
      call ids_put_slice(idx,'radiation',radiation_ids,stat)
    endif

    if (stat==0) then
       write(*,*) '    Radiation IDS exported'
    else
       write(*,*) '    Something went wrong writting the radiation IDS!'
    endif

  end subroutine fill_radiation_IDS




 


  ! --- Fills Bezier coefficients in GGD
  subroutine fill_Bezier_coefficients( ggd_scalar, node_list, var_index, grid_ind, grid_sub_ind, res_fact )
  
    implicit none
  
    ! --- External parameters
    type(ids_generic_grid_scalar),  intent(inout) ::  ggd_scalar
    type (type_node_list), intent(in)  :: node_list
    integer,               intent(in)  :: var_index, grid_ind, grid_sub_ind
    real*8,                intent(in)  :: res_fact
  
    ! --- Local parameters
    integer :: num_nodes, idof, inode, icoeff
    
    num_nodes = node_list%n_nodes
  
    if ( associated(ggd_scalar%coefficients) ) then
      deallocate( ggd_scalar%coefficients )
      allocate( ggd_scalar%coefficients(n_tor, num_nodes*(n_order+1)) ) 
    else
      allocate( ggd_scalar%coefficients(n_tor, num_nodes*(n_order+1)) ) 
    endif
  
    do inode=1, num_nodes
      do idof=1, n_order+1
  
        icoeff = inode + (idof-1)*num_nodes
  
        ggd_scalar%coefficients(:,icoeff)=node_list%node(inode)%values(:,idof,var_index)
  
      enddo
    enddo
  
    ggd_scalar%grid_index = grid_ind
    ggd_scalar%grid_subset_index = grid_sub_ind
 
   ! --- Re-scale coefficients
    ggd_scalar%coefficients = ggd_scalar%coefficients * res_fact
  
  end subroutine fill_Bezier_coefficients
  
  
  
  
  !< Fills JOREK grid into GGD
  subroutine grid2ggd( grid, node_list, element_list )
  
    implicit none
  
    ! --- External parameters
    type(ids_generic_grid_aos3_root),     pointer   :: grid
    type (type_node_list),    intent(in)            :: node_list
    type (type_element_list), intent(in)            :: element_list
  
  
    ! --- Local parameters
    type(ids_generic_grid_dynamic_space), pointer   ::  space_RZ
    type(ids_generic_grid_dynamic_space), pointer   ::  space_fourier
    type(ids_generic_grid_dynamic_space_dimension), pointer :: ids_cells
    
    integer:: idx, shot_number, run_number, num_nodes, num_cells   
    integer :: gs_index, i, j
    integer, allocatable :: vertex_elm_array(:,:)
    real*8,  allocatable :: RZ(:,:)
    
    integer :: itor, idof, n_slice, i_slice, grid_ind, grid_sub_ind, n_grid_sub
  
    ! --- create vertex - elements array
    allocate(  vertex_elm_array(n_vertex_max, element_list%n_elements)  )
    do i=1, element_list%n_elements
      vertex_elm_array(:,i) = element_list%element(i)%vertex
    enddo 
  
    ! Get values of R and Z at the nodes
    allocate( RZ(node_list%n_nodes, 2) )
    do i=1, node_list%n_nodes
      RZ(i,:) = node_list%node(i)%x(1,1,:)  
    enddo
  
    ! Write grid geometry
    allocate(  grid%space(2)                           )
    allocate(  grid%space(1)%objects_per_dimension(4)  )
    allocate(  grid%space(1)%coordinates_type(2)       )
  
    ! Set coordinates type to [R, Z]
    grid%space(1)%coordinates_type = (/ 4, 3 /)
  
    allocate(grid%identifier%description(1))
    allocate(grid%identifier%name(1))
    grid%identifier%description(1) = "Mesh JOREK output HDF5 file grid with quantities"
    grid%identifier%name = "JOREK output HDF5 file grid with quantities"
    grid%identifier%index = 1
  
    num_nodes = size(RZ,1)
    space_RZ  => grid%space(1)
  
    ! Fill simplified grid nodes (uses geometry instead of geometry_2D and misses Bezier representation) 
    allocate( space_RZ%objects_per_dimension(1)%object(num_nodes) )  ! Allocate to number of nodes
    do i=1, num_nodes 
      allocate( space_RZ%objects_per_dimension(1)%object(i)%geometry(2) ) ! Allocate dimensions per each node
      space_RZ%objects_per_dimension(1)%object(i)%geometry(:) = RZ(i,:)
    enddo
  
    ! Fill dummy variables for 1D elements (edges)
    allocate( space_RZ%objects_per_dimension(2)%object(1) )          ! Allocate just one edge    
    allocate( space_RZ%objects_per_dimension(2)%object(1)%nodes(1))  
    space_RZ%objects_per_dimension(2)%object(1)%nodes(1) = 0
  
    ! Fill JOREK 2D elements (or cells)
    num_cells = element_list%n_elements
    ids_cells => space_RZ%objects_per_dimension(3)
    allocate(    ids_cells%object(num_cells) )
    do i=1, num_cells
      allocate(  ids_cells%object(i)%nodes(n_vertex_max)  )
      ids_cells%object(i)%nodes(:) = vertex_elm_array(:,i) 
    enddo
  
    ! Writing grid_subsets
    allocate(grid%grid_subset(2))  ! 2 grid subsets 
  
    ! Subset for points
    gs_index = 1
  
    allocate( grid%grid_subset(gs_index)%identifier%name(1)         )
    allocate( grid%grid_subset(gs_index)%identifier%description(1)  )
    grid%grid_subset(gs_index)%identifier%name(1)        = "Nodes"
    grid%grid_subset(gs_index)%identifier%index          = gs_index 
    grid%grid_subset(gs_index)%identifier%description(1) = "All points/nodes/vertices/0D objects in the domain."
    grid%grid_subset(gs_index)%dimension                 = 1
   
    allocate( grid%grid_subset(gs_index)%element(num_nodes) )
   
    do i=1, num_nodes  
      allocate(  grid%grid_subset(gs_index)%element(i)%object(1)   )
      grid%grid_subset(gs_index)%element(i)%object(1)%space     = 1
      grid%grid_subset(gs_index)%element(i)%object(1)%index     = i 
      grid%grid_subset(gs_index)%element(i)%object(1)%dimension = 1
    enddo 
  
    ! Subset for cells
    gs_index = 2
  
    allocate(grid%grid_subset(gs_index)%identifier%name(1))
    allocate(grid%grid_subset(gs_index)%identifier%description(1))
  
    grid%grid_subset(gs_index)%identifier%name(1)        = "2D cells "
    grid%grid_subset(gs_index)%identifier%index          = gs_index 
    grid%grid_subset(gs_index)%identifier%description(1) = "All points/nodes/vertices/0D objects in the domain."
    grid%grid_subset(gs_index)%dimension                 = 3   ! A bit confusing, but 1 means point, 2 line and 3 surface
   
    allocate( grid%grid_subset(gs_index)%element(num_cells) )
  
    do i=1, num_cells
      allocate(  grid%grid_subset(gs_index)%element(i)%object(1)   )
      grid%grid_subset(gs_index)%element(i)%object(1)%space     = 1
      grid%grid_subset(gs_index)%element(i)%object(1)%index     = i 
      grid%grid_subset(gs_index)%element(i)%object(1)%dimension = 3
    enddo 
  
    ! Fill toroidal space (must be adapted for multiple time slices?)
    space_fourier  => grid%space(2)
    allocate(    space_fourier%coordinates_type(1)    )
    allocate(    space_fourier%identifier%description(1)  )
    space_fourier%coordinates_type(1) = 5          ! The coordinate type is 5, phi angle
    space_fourier%geometry_type%index = n_period   ! Fourier periodicity
    space_fourier%identifier%description(1) = "Toroidal space"             
  
    allocate(  space_fourier%objects_per_dimension(1)               )  ! We have only one dimension of
    allocate(  space_fourier%objects_per_dimension(1)%object(n_tor) )  ! toroidal harmonics
  
    do i=1, n_tor
      allocate(  space_fourier%objects_per_dimension(1)%object(i)%geometry(1) )  ! toroidal harmonics
      space_fourier%objects_per_dimension(1)%object(i)%geometry(1) = i
    enddo
  
    ! Fill in grid Bezier coefficients
    ! Needs generalization for JOREK 3D STELLERATOR EXTENSION!!
    space_RZ%geometry_type%index = 0  ! Standard geometry (non Fourier)
    do i=1, num_nodes
      allocate( space_RZ%objects_per_dimension(1)%object(i)%geometry_2d(2,n_order+1) )
      space_RZ%objects_per_dimension(1)%object(i)%geometry_2d(1,:) = node_list%node(i)%x(1,:,1 )   ! R dofs
      space_RZ%objects_per_dimension(1)%object(i)%geometry_2d(2,:) = node_list%node(i)%x(1,:,2 )   ! Z dofs
    enddo
  
    ! JOREK element sizes
    do i=1, num_cells
      allocate( space_RZ%objects_per_dimension(3)%object(i)%geometry_2d(n_order+1, n_vertex_max) )
      do j=1, n_vertex_max
        space_RZ%objects_per_dimension(3)%object(i)%geometry_2d(:,j) = element_list%element(i)%size(j,:)   
      enddo
    enddo
  
  end subroutine grid2ggd




  ! --- Checks if a restart file exists in the current directory
  logical function restart_file_exists(i_step)

    use phys_module, only : rst_hdf5

    implicit none

    integer,  intent(in) :: i_step
    character(len=64)    :: file_name

    write(file_name,'(a,i5.5)') 'jorek', i_step
    if ( rst_hdf5 .ne. 0 ) then
      inquire (file=trim(file_name)//'.h5', exist=restart_file_exists)
    else
      inquire (file=trim(file_name)//'.rst', exist=restart_file_exists)
    end if

  end function restart_file_exists



#endif

end module mod_jorek2IMAS

!> Routines to import a restart file written out by a routine in [[export_restart]].
module mod_import_restart
implicit none
contains
!> Imports a restart file written out by the routine export_restart.

subroutine import_restart(node_list, element_list, filename, format_rst, ierr, no_perturbations)

  use tr_module
  use data_structure
  use phys_module
  use pellet_module
  use equil_info
  use mod_boundary
  use basis_at_gaussian

  implicit none
  
  ! --- Routine parameters
  type(type_node_list),    intent(inout) :: node_list
  type(type_element_list), intent(inout) :: element_list
  character*(*)          , intent(in)    :: filename
  integer,                 intent(out)   :: ierr
  integer,                 intent(in)    :: format_rst  ! format of restart file 
  logical, optional,       intent(in)    :: no_perturbations ! don't initialize new harmonics
 
  ! --- Local parameters
  type (type_bnd_element_list)           :: bnd_elm_list    
  type (type_bnd_node_list)              :: bnd_node_list 

  if ( rst_hdf5 == 0 ) then
    write(*,*) " Restart from BINARY file " // trim(filename) // '.rst'
    call import_binary_restart(node_list, element_list, trim(filename)//'.rst', &
            format_rst, ierr, no_perturbations)
  else if ( rst_hdf5 == 1 ) then
    write(*,*) " Restart from HDF5 file " // trim(filename) // '.h5'
    call import_hdf5_restart(node_list, element_list, trim(filename)//'.h5', &
            format_rst,ierr, no_perturbations)
  end if
  
  ! --- Required initializations to update equilibrium state
  call initialise_basis
  call boundary_from_grid(node_list, element_list, bnd_node_list, bnd_elm_list, .false.)
  call update_equil_state(0,node_list, element_list, bnd_elm_list, xpoint, xcase)
  write(*,*) " "
  write(*,*) " The equilibrium state has been updated "
  
end subroutine import_restart


!
! Import a binary restart file
subroutine import_binary_restart(node_list, element_list, filename, format_rst, error, no_perturbations)

  use tr_module 
  use data_structure
  use phys_module
  use pellet_module
  use vacuum, only: import_restart_vacuum, current_FB_fact
  use mod_element_rtree, only: populate_element_rtree
  
  implicit none
  
  ! --- Routine parameters
  type(type_node_list),    intent(inout) :: node_list
  type(type_element_list), intent(inout) :: element_list
  character(len=*),        intent(in)    :: filename
  integer,                 intent(out)   :: error
  integer,                 intent(in)    :: format_rst  ! format of restart file
  logical, optional,       intent(in)    :: no_perturbations ! don't initialize new harmonics
  
  ! --- Local variables
  integer              :: i, j, m, k, n_tor_tmp
  real*8               :: growth_mag, growth_kin, amplitude
  integer, allocatable :: mode_tmp(:)
  real*8,  allocatable :: values_tmp(:,:,:), deltas_tmp(:,:,:)
  real*8,  allocatable :: spi_R_arr (:)
  real*8,  allocatable :: spi_Z_arr (:)
  real*8,  allocatable :: spi_phi_arr (:)
  real*8,  allocatable :: spi_phi_init_arr (:)
  real*8,  allocatable :: spi_Vel_R_arr (:)
  real*8,  allocatable :: spi_Vel_Z_arr (:)
  real*8,  allocatable :: spi_Vel_RxZ_arr (:)
  real*8,  allocatable :: spi_radius_arr (:)
  real*8,  allocatable :: spi_abl_arr (:)
  real*8, allocatable  :: spi_species_arr (:)

  integer              :: n_spi_check, n_inj_check
  logical              :: modes_changed
 
  real*8, allocatable :: t_energies(:,:,:)   !< Magnetic and kinetic mode energies at previous timesteps.
  real*8, allocatable :: t_energies2(:,:,:)  !< Magnetic and kinetic mode energies at previous timesteps.
  real*8, allocatable :: t_energies3(:,:,:)  !< Magnetic and kinetic mode energies at previous timesteps.
  real*8, allocatable :: t_energies4(:,:,:)  !< Magnetic and kinetic mode energies at previous timesteps.

  ! --- Perturbation-Import variables
  type (type_node_list)   , pointer	:: node_list_perturbation
  type (type_element_list), pointer	:: element_list_perturbation
  integer              			:: n_tor_tmp_perturbation
  integer, allocatable 			:: mode_tmp_perturbation(:)
  real*8,  allocatable 			:: values_tmp_perturbation(:,:,:), deltas_tmp_perturbation(:,:,:)
  logical, parameter   			:: import_perturbation = .false.
  logical                               :: no_pert
  
  no_pert = .false.
  if ( present(no_perturbations) ) no_pert = no_perturbations

  error = 0

  write(*,*) 'Importing restart file "', trim(filename), '".'
  write(*,*) ' Using format : ',rst_format

  open(21,file=trim(filename), form='unformatted', status='old', action='read', iostat=error)
  if ( error /= 0 ) then
    write(*,*) '...failed!'
    return
  end if

  read(21) n_tor_tmp

  allocate(mode_tmp(n_tor_tmp), values_tmp(n_tor_tmp,n_order+1,n_var), deltas_tmp(n_tor_tmp,n_order+1,n_var))

  if (format_rst == 1) then
    read(21) mode_tmp
    write(*,*) ' NEW format (1) : ',mode_tmp
  elseif (format_rst == 0) then
    write(*,*) ' mode : ',mode
    modes_changed = .false.
    if (n_tor .eq. n_tor_tmp) then 
       mode_tmp = mode
    else
       mode_tmp(1:min(n_tor,n_tor_tmp)) = mode(1:min(n_tor,n_tor_tmp))
       modes_changed = .true.
    endif
    if (modes_changed) then
      write(*,*) 'OLD format (0) : '
      write(*,'(A,999i4)') '  previous modenumbers : ',mode_tmp
      write(*,'(A,999i4)') '  new mode numbers     : ',mode
    else
      write(*,*) 'OLD format (0)'
    endif
  elseif ( format_rst > 2 ) then
    write(*,'(A,i3)') ' restart file format not supported : ',format_rst
    stop
  endif

  if (n_tor_tmp .gt. n_tor) write(*,'(3(a,i4))') &
       ' Warning: Reducing number of harmonics from', n_tor_tmp, ' to', n_tor, '!'
  if (n_tor_tmp .lt. n_tor) write(*,'(3(a,i4))') &
       ' Warning: Increasing number of harmonics from', n_tor_tmp, ' to', n_tor, '!'

  write(*,'(A,i5,A)') ' Importing ',n_tor_tmp,' harmonics'

  read(21) node_list%n_nodes,element_list%n_elements
  read(21) node_list%n_dof

  do i=1,node_list%n_nodes
     
    read(21) node_list%node(i)%x
    read(21) values_tmp
    read(21) deltas_tmp

#ifdef fullmhd
    read(21) node_list%node(i)%psi_eq               !< equilibrium flux at the nodes
    read(21) node_list%node(i)%Fprof_eq             !< equilibrium profile R*B_phi at the nodes
#elif altcs
    read(21) node_list%node(i)%psi_eq               !< equilibrium flux at the nodes
#endif
    read(21) node_list%node(i)%index
    read(21) node_list%node(i)%boundary
    read(21) node_list%node(i)%axis_node
    read(21) node_list%node(i)%parents
    read(21) node_list%node(i)%parent_elem
    read(21) node_list%node(i)%ref_lambda
    read(21) node_list%node(i)%ref_mu
    read(21) node_list%node(i)%constrained

    node_list%node(i)%values = 0.d0
    node_list%node(i)%deltas = 0.d0

    do m=1,n_tor_tmp,2
      do k=1, n_tor,2
        if (mode_tmp(m) .eq. mode(k)) then
          if ((m .eq. 1) .and. (k.eq.1)) then
            node_list%node(i)%values(k,:,:) = values_tmp(m,:,:)
            node_list%node(i)%deltas(k,:,:) = deltas_tmp(m,:,:)
          else
            node_list%node(i)%values(k-1,:,:) = values_tmp(m-1,:,:)
            node_list%node(i)%deltas(k-1,:,:) = deltas_tmp(m-1,:,:)
            node_list%node(i)%values(k,:,:)   = values_tmp(m,:,:)
            node_list%node(i)%deltas(k,:,:)   = deltas_tmp(m,:,:)
          endif
        endif
      enddo
    enddo
  enddo

  read(21) element_list%element(1:element_list%n_elements)
  read(21) tstep,eta_rst,visco_rst,visco_par_rst
  read(21) index_start
  read(21) t_start
  
  write(*,*) 'CHECK (1): allocating energies in import_restart : ',index_start,index_start+nstep

  if (index_start .ge. 1) then

    write(*,*) 'CHECK (2): allocating energies in import_restart : ',index_start,index_start+nstep

    if (allocated(xtime)) call tr_deallocate(xtime,"xtime",CAT_UNKNOWN)
    call tr_allocate(xtime,1,index_start+nstep,"xtime",CAT_UNKNOWN)
    
    if (allocated(t_energies))   call tr_deallocate(t_energies,"t_energies",CAT_UNKNOWN)
    call tr_allocate(t_energies,1,n_tor_tmp,1,2,1,index_start+nstep,"t_energies",CAT_UNKNOWN)
    t_energies = 0.d0
    
    if (allocated(energies)) call tr_deallocate(energies,"energies",CAT_UNKNOWN)
    call tr_allocate(energies,1,n_tor,1,2,1,index_start+nstep,"energies",CAT_UNKNOWN)
    energies = 0.d0
    
    if (allocated(R_axis_t)) call tr_deallocate(R_axis_t,"R_axis_t",CAT_UNKNOWN)
    call tr_allocate(R_axis_t,1,index_start+nstep,"R_axis_t",CAT_UNKNOWN)
    R_axis_t = 0.d0
    
    if (allocated(Z_axis_t)) call tr_deallocate(Z_axis_t,"Z_axis_t",CAT_UNKNOWN)
    call tr_allocate(Z_axis_t,1,index_start+nstep,"Z_axis_t",CAT_UNKNOWN)
    Z_axis_t = 0.d0
    
    if (allocated(psi_axis_t)) call tr_deallocate(psi_axis_t,"psi_axis_t",CAT_UNKNOWN)
    call tr_allocate(psi_axis_t,1,index_start+nstep,"psi_axis_t",CAT_UNKNOWN)
    psi_axis_t = 0.d0
    
    if (allocated(R_xpoint_t)) call tr_deallocate(R_xpoint_t,"R_xpoint_t",CAT_UNKNOWN)
    call tr_allocate(R_xpoint_t,1,index_start+nstep,1,2,"R_xpoint_t",CAT_UNKNOWN)
    R_xpoint_t = 0.d0
    
    if (allocated(Z_xpoint_t)) call tr_deallocate(Z_xpoint_t,"Z_xpoint_t",CAT_UNKNOWN)
    call tr_allocate(Z_xpoint_t,1,index_start+nstep,1,2,"Z_xpoint_t",CAT_UNKNOWN)
    Z_xpoint_t = 0.d0

    if (allocated(psi_xpoint_t)) call tr_deallocate(psi_xpoint_t,"psi_xpoint_t",CAT_UNKNOWN)
    call tr_allocate(psi_xpoint_t,1,index_start+nstep,1,2,"psi_xpoint_t",CAT_UNKNOWN)
    psi_xpoint_t = 0.d0

    if (allocated(R_bnd_t)) call tr_deallocate(R_bnd_t,"R_bnd_t",CAT_UNKNOWN)
    call tr_allocate(R_bnd_t,1,index_start+nstep,"R_bnd_t",CAT_UNKNOWN)
    R_bnd_t = 0.d0
    
    if (allocated(Z_bnd_t)) call tr_deallocate(Z_bnd_t,"Z_bnd_t",CAT_UNKNOWN)
    call tr_allocate(Z_bnd_t,1,index_start+nstep,"Z_bnd_t",CAT_UNKNOWN)
    Z_bnd_t = 0.d0
    
    if (allocated(psi_bnd_t)) call tr_deallocate(psi_bnd_t,"psi_bnd_t",CAT_UNKNOWN)
    call tr_allocate(psi_bnd_t,1,index_start+nstep,"psi_bnd_t",CAT_UNKNOWN)
    psi_bnd_t = 0.d0
    
    if (allocated(current_t)) call tr_deallocate(current_t,"current_t",CAT_UNKNOWN)
    call tr_allocate(current_t,1,index_start+nstep,"current_t",CAT_UNKNOWN)
    current_t = 0.d0
    
    if (allocated(beta_p_t)) call tr_deallocate(beta_p_t,"beta_p_t",CAT_UNKNOWN)
    call tr_allocate(beta_p_t,1,index_start+nstep,"beta_p_t",CAT_UNKNOWN)
    beta_p_t = 0.d0
    
    if (allocated(beta_t_t)) call tr_deallocate(beta_t_t,"beta_t_t",CAT_UNKNOWN)
    call tr_allocate(beta_t_t,1,index_start+nstep,"beta_t_t",CAT_UNKNOWN)
    beta_t_t = 0.d0
    
    if (allocated(beta_n_t)) call tr_deallocate(beta_n_t,"beta_n_t",CAT_UNKNOWN)
    call tr_allocate(beta_n_t,1,index_start+nstep,"beta_n_t",CAT_UNKNOWN)
    beta_n_t = 0.d0
    
    if (allocated(density_in_t)) call tr_deallocate(density_in_t,"density_in_t",CAT_UNKNOWN)
    call tr_allocate(density_in_t,1,index_start+nstep,"density_in_t",CAT_UNKNOWN)
    density_in_t = 0.d0
    
    if (allocated(density_out_t)) call tr_deallocate(density_out_t,"density_out_t",CAT_UNKNOWN)
    call tr_allocate(density_out_t,1,index_start+nstep,"density_out_t",CAT_UNKNOWN)
    density_out_t = 0.d0
    
    if (allocated(pressure_in_t)) call tr_deallocate(pressure_in_t,"pressure_in_t",CAT_UNKNOWN)
    call tr_allocate(pressure_in_t,1,index_start+nstep,"pressure_in_t",CAT_UNKNOWN)
    pressure_in_t = 0.d0
    
    if (allocated(pressure_out_t)) call tr_deallocate(pressure_out_t,"pressure_out_t",CAT_UNKNOWN)
    call tr_allocate(pressure_out_t,1,index_start+nstep,"pressure_out_t",CAT_UNKNOWN)
    pressure_out_t = 0.d0
    
    if (allocated(heat_src_in_t)) call tr_deallocate(heat_src_in_t,"heating_power_t",CAT_UNKNOWN)
    call tr_allocate(heat_src_in_t,1,index_start+nstep,"heat_src_in_t",CAT_UNKNOWN)
    heat_src_in_t = 0.d0
    
    if (allocated(heat_src_out_t)) call tr_deallocate(heat_src_out_t,"heating_power_t",CAT_UNKNOWN)
    call tr_allocate(heat_src_out_t,1,index_start+nstep,"heat_src_out_t",CAT_UNKNOWN)
    heat_src_out_t = 0.d0
    
    if (allocated(part_src_in_t)) call tr_deallocate(part_src_in_t,"parting_power_t",CAT_UNKNOWN)
    call tr_allocate(part_src_in_t,1,index_start+nstep,"part_src_in_t",CAT_UNKNOWN)
    part_src_in_t = 0.d0
    
    if (allocated(part_src_out_t)) call tr_deallocate(part_src_out_t,"parting_power_t",CAT_UNKNOWN)
    call tr_allocate(part_src_out_t,1,index_start+nstep,"part_src_out_t",CAT_UNKNOWN)
    part_src_out_t = 0.d0
    
    if (allocated(E_tot_t)) call tr_deallocate(E_tot_t,"E_tot_t",CAT_UNKNOWN)
    call tr_allocate(E_tot_t,1,index_start+nstep,"E_tot_t",CAT_UNKNOWN)
    E_tot_t = 0.d0

    if (allocated(helicity_tot_t)) call tr_deallocate(helicity_tot_t,"helicity_tot_t",CAT_UNKNOWN)
    call tr_allocate(helicity_tot_t,1,index_start+nstep,"helicity_tot_t",CAT_UNKNOWN)
    helicity_tot_t = 0.d0

    if (allocated(thermal_tot_t)) call tr_deallocate(thermal_tot_t,"thermal_tot_t",CAT_UNKNOWN)
    call tr_allocate(thermal_tot_t,1,index_start+nstep,"thermal_tot_t",CAT_UNKNOWN)
    thermal_tot_t = 0.d0

    if (allocated(thermal_e_tot_t)) call tr_deallocate(thermal_e_tot_t,"thermal_e_tot_t",CAT_UNKNOWN)
    call tr_allocate(thermal_e_tot_t,1,index_start+nstep,"thermal_e_tot_t",CAT_UNKNOWN)
    thermal_e_tot_t = 0.d0

    if (allocated(thermal_i_tot_t)) call tr_deallocate(thermal_i_tot_t,"thermal_i_tot_t",CAT_UNKNOWN)
    call tr_allocate(thermal_i_tot_t,1,index_start+nstep,"thermal_i_tot_t",CAT_UNKNOWN)
    thermal_i_tot_t = 0.d0

    if (allocated(kin_par_tot_t)) call tr_deallocate(kin_par_tot_t,"kin_par_tot_t",CAT_UNKNOWN)
    call tr_allocate(kin_par_tot_t,1,index_start+nstep,"kin_par_tot_t",CAT_UNKNOWN)
    kin_par_tot_t = 0.d0

    if (allocated(kin_perp_tot_t)) call tr_deallocate(kin_perp_tot_t,"kin_perp_tot_t",CAT_UNKNOWN)
    call tr_allocate(kin_perp_tot_t,1,index_start+nstep,"kin_perp_tot_t",CAT_UNKNOWN)
    kin_perp_tot_t = 0.d0

    if (allocated(Ip_tot_t)) call tr_deallocate(Ip_tot_t,"Ip_tot_t",CAT_UNKNOWN)
    call tr_allocate(Ip_tot_t,1,index_start+nstep,"Ip_tot_t",CAT_UNKNOWN)
    Ip_tot_t = 0.d0

    if (allocated(ohmic_tot_t)) call tr_deallocate(ohmic_tot_t,"ohmic_tot_t",CAT_UNKNOWN)
    call tr_allocate(ohmic_tot_t,1,index_start+nstep,"ohmic_tot_t",CAT_UNKNOWN)
    ohmic_tot_t = 0.d0

    if (allocated(Wmag_tot_t)) call tr_deallocate(Wmag_tot_t,"Wmag_tot_t",CAT_UNKNOWN)
    call tr_allocate(Wmag_tot_t,1,index_start+nstep,"Wmag_tot_t",CAT_UNKNOWN)
    Wmag_tot_t = 0.d0

    if (allocated(Magwork_tot_t)) call tr_deallocate(Magwork_tot_t,"Magwork_tot_t",CAT_UNKNOWN)
    call tr_allocate(Magwork_tot_t,1,index_start+nstep,"Magwork_tot_t",CAT_UNKNOWN)
    Magwork_tot_t = 0.d0

    if (allocated(flux_qpar_t)) call tr_deallocate(flux_qpar_t,"flux_qpar_t",CAT_UNKNOWN)
    call tr_allocate(flux_qpar_t,1,index_start+nstep,"flux_qpar_t",CAT_UNKNOWN)
    flux_qpar_t = 0.d0

    if (allocated(flux_qperp_t)) call tr_deallocate(flux_qperp_t,"flux_qperp_t",CAT_UNKNOWN)
    call tr_allocate(flux_qperp_t,1,index_start+nstep,"flux_qperp_t",CAT_UNKNOWN)
    flux_qperp_t = 0.d0

    if (allocated(flux_kinpar_t)) call tr_deallocate(flux_kinpar_t,"flux_kinpar_t",CAT_UNKNOWN)
    call tr_allocate(flux_kinpar_t,1,index_start+nstep,"flux_kinpar_t",CAT_UNKNOWN)
    flux_kinpar_t = 0.d0

    if (allocated(flux_poynting_t)) call tr_deallocate(flux_poynting_t,"flux_poynting_t",CAT_UNKNOWN)
    call tr_allocate(flux_poynting_t,1,index_start+nstep,"flux_poynting_t",CAT_UNKNOWN)
    flux_poynting_t = 0.d0

    if (allocated(flux_Pvn_t)) call tr_deallocate(flux_Pvn_t,"flux_Pvn_t",CAT_UNKNOWN)
    call tr_allocate(flux_Pvn_t,1,index_start+nstep,"flux_Pvn_t",CAT_UNKNOWN)
    flux_Pvn_t = 0.d0

    if (allocated(dE_tot_dt)) call tr_deallocate(dE_tot_dt,"dE_tot_dt",CAT_UNKNOWN)
    call tr_allocate(dE_tot_dt,1,index_start+nstep,"dE_tot_dt",CAT_UNKNOWN)
    dE_tot_dt = 0.d0

    if (allocated(dWmag_tot_dt)) call tr_deallocate(dWmag_tot_dt,"dWmag_tot_dt",CAT_UNKNOWN)
    call tr_allocate(dWmag_tot_dt,1,index_start+nstep,"dWmag_tot_dt",CAT_UNKNOWN)
    dWmag_tot_dt = 0.d0

    if (allocated(dthermal_tot_dt)) call tr_deallocate(dthermal_tot_dt,"dthermal_tot_dt",CAT_UNKNOWN)
    call tr_allocate(dthermal_tot_dt,1,index_start+nstep,"dthermal_tot_dt",CAT_UNKNOWN)
    dthermal_tot_dt = 0.d0

    if (allocated(dkinperp_tot_dt)) call tr_deallocate(dkinperp_tot_dt,"dkinperp_tot_dt",CAT_UNKNOWN)
    call tr_allocate(dkinperp_tot_dt,1,index_start+nstep,"dkinperp_tot_dt",CAT_UNKNOWN)
    dkinperp_tot_dt = 0.d0

    if (allocated(dkinpar_tot_dt)) call tr_deallocate(dkinpar_tot_dt,"dkinpar_tot_dt",CAT_UNKNOWN)
    call tr_allocate(dkinpar_tot_dt,1,index_start+nstep,"dkinpar_tot_dt",CAT_UNKNOWN)
    dkinpar_tot_dt = 0.d0

    if (allocated(heat_src_tot_t)) call tr_deallocate(heat_src_tot_t,"heat_src_tot_t",CAT_UNKNOWN)
    call tr_allocate(heat_src_tot_t,1,index_start+nstep,"heat_src_tot_t",CAT_UNKNOWN)
    heat_src_tot_t = 0.d0

    if (allocated(part_src_tot_t)) call tr_deallocate(part_src_tot_t,"part_src_tot_t",CAT_UNKNOWN)
    call tr_allocate(part_src_tot_t,1,index_start+nstep,"part_src_tot_t",CAT_UNKNOWN)
    part_src_tot_t = 0.d0

    if (allocated(li3_tot_t)) call tr_deallocate(li3_tot_t,"li3_tot_t",CAT_UNKNOWN)
    call tr_allocate(li3_tot_t,1,index_start+nstep,"li3_tot_t",CAT_UNKNOWN)
    li3_tot_t = 0.d0

    if (allocated(li3_t)) call tr_deallocate(li3_t,"li3_t",CAT_UNKNOWN)
    call tr_allocate(li3_t,1,index_start+nstep,"li3_t",CAT_UNKNOWN)
    li3_t = 0.d0

    if (allocated(viscopar_flux_t)) call tr_deallocate(viscopar_flux_t,"viscopar_flux_t",CAT_UNKNOWN)
    call tr_allocate(viscopar_flux_t,1,index_start+nstep,"viscopar_flux_t",CAT_UNKNOWN)
    viscopar_flux_t = 0.d0

    if (allocated(viscopar_dissip_tot_t)) call tr_deallocate(viscopar_dissip_tot_t,"viscopar_dissip_tot_t",CAT_UNKNOWN)
    call tr_allocate(viscopar_dissip_tot_t,1,index_start+nstep,"viscopar_dissip_tot_t",CAT_UNKNOWN)
    viscopar_dissip_tot_t = 0.d0

    if (allocated(friction_dissip_tot_t)) call tr_deallocate(friction_dissip_tot_t,"friction_dissip_tot_t",CAT_UNKNOWN)
    call tr_allocate(friction_dissip_tot_t,1,index_start+nstep,"friction_dissip_tot_t",CAT_UNKNOWN)
    friction_dissip_tot_t = 0.d0

    if (allocated(thmwork_tot_t)) call tr_deallocate(thmwork_tot_t,"thmwork_tot_t",CAT_UNKNOWN)
    call tr_allocate(thmwork_tot_t,1,index_start+nstep,"thmwork_tot_t",CAT_UNKNOWN)
    thmwork_tot_t = 0.d0

    if (allocated(volume_t)) call tr_deallocate(volume_t,"volume_t",CAT_UNKNOWN)
    call tr_allocate(volume_t,1,index_start+nstep,"volume_t",CAT_UNKNOWN)
    volume_t = 0.d0

    if (allocated(area_t)) call tr_deallocate(area_t,"area_t",CAT_UNKNOWN)
    call tr_allocate(area_t,1,index_start+nstep,"area_t",CAT_UNKNOWN)
    area_t = 0.d0

    if (allocated(mag_ener_src_tot)) call tr_deallocate(mag_ener_src_tot,"mag_ener_src_tot",CAT_UNKNOWN)
    call tr_allocate(mag_ener_src_tot,1,index_start+nstep,"mag_ener_src_tot",CAT_UNKNOWN)
    mag_ener_src_tot = 0.d0

#ifdef JECCD
    if (allocated(energies2)) call tr_deallocate(energies2,"energies2",CAT_UNKNOWN)
    call tr_allocate(energies2,1,n_tor,1,2,1,index_start+nstep,"energies2",CAT_UNKNOWN)

    if (allocated(energies3)) call tr_deallocate(energies3,"energies3",CAT_UNKNOWN)
    call tr_allocate(energies3,1,n_tor,1,2,1,index_start+nstep,"energies3",CAT_UNKNOWN)

#ifdef JEC2DIAG
    if (allocated(energies4)) call tr_deallocate(energies4,"energies4",CAT_UNKNOWN)
    call tr_allocate(energies4,1,n_tor,1,2,1,index_start+nstep,"energies4",CAT_UNKNOWN)
#endif

    energies2 = 0.d0
    energies3 = 0.d0
#ifdef JEC2DIAG
    energies4 = 0.d0
#endif
#endif

  read(21) xtime(1:index_start)
  read(21) energies(1:n_tor_tmp,:,1:index_start)

#ifdef JECCD
  read(21) energies2(1:n_tor_tmp,:,1:index_start)
  read(21) energies3(1:n_tor_tmp,:,1:index_start)
#ifdef JEC2DIAG
  read(21) energies4(1:n_tor_tmp,:,1:index_start)
#endif
#endif
endif

  call import_restart_vacuum(21, freeboundary, resistive_wall)  
  
  !--- Some parameters need to be scaled when importing a free-boundary equilibrium
  T_0  = T_0  * current_FB_fact / prev_FB_fact
  T_1  = T_1  * current_FB_fact / prev_FB_fact
  FF_0 = FF_0 * current_FB_fact / prev_FB_fact
  FF_1 = FF_1 * current_FB_fact / prev_FB_fact
  prev_FB_fact = current_FB_fact

  if (use_pellet) then
    if (index_start .ge. 1) then
      if (allocated(xtime_pellet_R)) call tr_deallocate(xtime_pellet_R,"xtime_pellet_R",CAT_UNKNOWN)
      call tr_allocate(xtime_pellet_R,1,index_start+nstep,"xtime_pellet_R",CAT_UNKNOWN)
      if (allocated(xtime_pellet_Z)) call tr_deallocate(xtime_pellet_Z,"xtime_pellet_Z",CAT_UNKNOWN)
      call tr_allocate(xtime_pellet_Z,1,index_start+nstep,"xtime_pellet_Z",CAT_UNKNOWN)
      if (allocated(xtime_pellet_psi)) call tr_deallocate(xtime_pellet_psi,"xtime_pellet_psi",CAT_UNKNOWN)
      call tr_allocate(xtime_pellet_psi,1,index_start+nstep,"xtime_pellet_psi",CAT_UNKNOWN)
      if (allocated(xtime_pellet_particles)) &
           call tr_deallocate(xtime_pellet_particles,"xtime_pellet_particles",CAT_UNKNOWN)
      call tr_allocate(xtime_pellet_particles,1,index_start+nstep,"xtime_pellet_particles",CAT_UNKNOWN)
      if (allocated(xtime_phys_ablation)) &
           call tr_deallocate(xtime_phys_ablation,"xtime_phys_ablation",CAT_UNKNOWN)
      call tr_allocate(xtime_phys_ablation,1,index_start+nstep,"xtime_phys_ablation",CAT_UNKNOWN)

      read(21,err=999, end=999)  xtime_pellet_R(1:index_start)
      read(21)  xtime_pellet_Z(1:index_start)
      read(21)  xtime_pellet_psi(1:index_start)
      read(21)  xtime_pellet_particles(1:index_start)
      read(21)  xtime_phys_ablation(1:index_start)
    endif
    read(21,err=999, end=999)  pellet_particles, pellet_R, pellet_Z
    write(*,'(A,e12.4,2f10.5)') ' *** PELLET PARAMETERS : ',pellet_particles, pellet_R, pellet_Z
  endif

#if (defined WITH_Neutrals) || (defined WITH_Impurities)
  if (index_start >= 1) then
    if (allocated(xtime_radiation)) &
      call tr_deallocate(xtime_radiation,"xtime_radiation",CAT_UNKNOWN)
    call tr_allocate(xtime_radiation,1,index_start+nstep,"xtime_radiation",CAT_UNKNOWN)
    read(21)  xtime_radiation(1:index_start)
    if (allocated(xtime_rad_power)) &
      call tr_deallocate(xtime_rad_power,"xtime_rad_power",CAT_UNKNOWN)
    call tr_allocate(xtime_rad_power,1,index_start+nstep,"xtime_rad_power",CAT_UNKNOWN)
    read(21)  xtime_rad_power(1:index_start)
    if (allocated(xtime_E_ion)) &
      call tr_deallocate(xtime_E_ion,"xtime_E_ion",CAT_UNKNOWN)
    call tr_allocate(xtime_E_ion,1,index_start+nstep,"xtime_E_ion",CAT_UNKNOWN)
    read(21)  xtime_E_ion(1:index_start)
    if (allocated(xtime_E_ion_power)) &
      call tr_deallocate(xtime_E_ion_power,"xtime_E_ion_power",CAT_UNKNOWN)
    call tr_allocate(xtime_E_ion_power,1,index_start+nstep,"xtime_E_ion_power",CAT_UNKNOWN)
    read(21)  xtime_E_ion_power(1:index_start)
    if (allocated(xtime_P_ei)) &
      call tr_deallocate(xtime_P_ei,"xtime_P_ei",CAT_UNKNOWN)
    call tr_allocate(xtime_P_ei,1,index_start+nstep,"xtime_P_ei",CAT_UNKNOWN)
    read(21)  xtime_P_ei(1:index_start)
  end if
#endif

  if (using_spi) then
    if (n_spi_tot >= 1) then

      if (index_start >= 1) then

        if (allocated(xtime_spi_ablation)) &
          call tr_deallocate(xtime_spi_ablation,"xtime_spi_ablation",CAT_UNKNOWN)
        call tr_allocate(xtime_spi_ablation,1,n_spi_tot,1,index_start+nstep,"xtime_spi_ablation",CAT_UNKNOWN)
        if (allocated(xtime_spi_ablation_rate)) &
          call tr_deallocate(xtime_spi_ablation_rate,"xtime_spi_ablation_rate",CAT_UNKNOWN)
        call tr_allocate(xtime_spi_ablation_rate,1,n_spi_tot,1,index_start+nstep,"xtime_spi_ablation_rate",CAT_UNKNOWN)
        if (allocated(xtime_spi_ablation_bg)) &
          call tr_deallocate(xtime_spi_ablation_bg,"xtime_spi_ablation_bg",CAT_UNKNOWN)
        call tr_allocate(xtime_spi_ablation_bg,1,n_spi_tot,1,index_start+nstep,"xtime_spi_ablation_bg",CAT_UNKNOWN)
        if (allocated(xtime_spi_ablation_bg_rate)) &
          call tr_deallocate(xtime_spi_ablation_bg_rate,"xtime_spi_ablation_bg_rate",CAT_UNKNOWN)
        call tr_allocate(xtime_spi_ablation_bg_rate,1,n_spi_tot,1,index_start+nstep,"xtime_spi_ablation_bg_rate",CAT_UNKNOWN)

        read(21)  xtime_spi_ablation(1:n_spi_tot,1:index_start)
        read(21)  xtime_spi_ablation_rate(1:n_spi_tot,1:index_start)
        read(21)  xtime_spi_ablation_bg(1:n_spi_tot,1:index_start)
        read(21)  xtime_spi_ablation_bg_rate(1:n_spi_tot,1:index_start)
      end if

      read(21,err=999, end=999) n_spi_check

      if (n_spi_check /= n_spi_tot) then
        write(*,*) "Inconsistency in n_spi_tot detected, exiting!"
        stop
      end if

      read(21,err=999, end=999) n_inj_check

      if (n_inj_check /= n_inj) then
        write(*,*) "Inconsistency in n_inj detected, exiting!"
        stop
      end if      

      allocate (spi_R_arr(n_spi_tot))
      allocate (spi_Z_arr(n_spi_tot))
      allocate (spi_phi_arr(n_spi_tot))
      allocate (spi_phi_init_arr(n_spi_tot))
      allocate (spi_Vel_R_arr(n_spi_tot))
      allocate (spi_Vel_Z_arr(n_spi_tot))
      allocate (spi_Vel_RxZ_arr(n_spi_tot))
      allocate (spi_radius_arr(n_spi_tot))
      allocate (spi_abl_arr(n_spi_tot))
      allocate (spi_species_arr(n_spi_tot))
    
      read(21,err=999, end=999)  spi_R_arr(1:n_spi_tot)
      read(21,err=999, end=999)  spi_Z_arr(1:n_spi_tot)
      read(21,err=999, end=999)  spi_phi_arr(1:n_spi_tot)
      read(21,err=999, end=999)  spi_phi_init_arr(1:n_spi_tot)
      read(21,err=999, end=999)  spi_Vel_R_arr(1:n_spi_tot)
      read(21,err=999, end=999)  spi_Vel_Z_arr(1:n_spi_tot)
      read(21,err=999, end=999)  spi_Vel_RxZ_arr(1:n_spi_tot)
      read(21,err=999, end=999)  spi_radius_arr(1:n_spi_tot)
      read(21,err=999, end=999)  spi_abl_arr(1:n_spi_tot)
      read(21,err=999, end=999)  spi_species_arr(1:n_spi_tot)

      do i=1, n_spi_tot
        pellets(i)%spi_R       = spi_R_arr(i)
        pellets(i)%spi_Z       = spi_Z_arr(i)
        pellets(i)%spi_phi     = spi_phi_arr(i)
        pellets(i)%spi_phi_init= spi_phi_init_arr(i)
        pellets(i)%spi_Vel_R   = spi_Vel_R_arr(i)
        pellets(i)%spi_Vel_Z   = spi_Vel_Z_arr(i)
        pellets(i)%spi_Vel_RxZ = spi_Vel_RxZ_arr(i)
        pellets(i)%spi_radius  = spi_radius_arr(i)
        pellets(i)%spi_abl     = spi_abl_arr(i)
        pellets(i)%spi_species = spi_species_arr(i)

        write(*,'(A,I5,6ES10.2)') ' *** SHATTERED PELLET PARAMETERS : ',i, pellets(i)%spi_R, pellets(i)%spi_Z, &
                        pellets(i)%spi_phi, pellets(i)%spi_Vel_R, pellets(i)%spi_Vel_Z, pellets(i)%spi_radius
      end do

      deallocate (spi_R_arr)
      deallocate (spi_Z_arr)
      deallocate (spi_phi_arr)
      deallocate (spi_phi_init_arr)
      deallocate (spi_Vel_R_arr)
      deallocate (spi_Vel_Z_arr)
      deallocate (spi_Vel_RxZ_arr)
      deallocate (spi_radius_arr)
      deallocate (spi_abl_arr)
      deallocate (spi_species_arr)

      if (spi_tor_rot) then
        read(21,err=999, end=999) ns_phi_rotate 
      end if

    end if
  end if

999 continue
  
  close(21)
  
  write(*,*) '************* restart ******************'
  write(*,'(A,I6,F14.6,A)') ' *  restart time       : ',index_start,t_start,' *'
  write(*,*) '****************************************'

  do i=2,index_start
    if ( (energies(n_tor,1,i).ne.0.) .and. (energies(n_tor,1,i-1).ne.0.)) then
      Growth_mag  = 0.5d0*log(abs(energies(n_tor,1,i)/energies(n_tor,1,i-1))) &
    	    / (xtime(i)-xtime(i-1))
    else
       Growth_mag  = 0.
    endif
    if ( (energies(n_tor,2,i).ne.0.) .and. (energies(n_tor,2,i-1).ne.0.)) then
      Growth_kin  = 0.5d0*log(abs(energies(n_tor,2,i)/energies(n_tor,2,i-1))) &
    	    / (xtime(i)-xtime(i-1))
    else
      Growth_kin  = 0.
    endif

!     write(*,'(i7,f10.3,200e14.6)') i,xtime(i),energies(1:n_tor,:,i),growth_mag,growth_kin
!     write(*,'(i7,f10.3,200e14.6)') i,xtime(i),energies(1:n_tor,:,i)

  enddo
  
  ! --- initialise new harmonics (only density and temperature, to be improved)
  if ( (.not. no_pert) .and. (n_tor_tmp .lt. n_tor) ) then
    ! --- Using an already computated mode
    if ( (import_perturbation) .and. (n_tor .gt. 1) ) then
      
      write(*,*)'Importing perturbation from jorek_perturbation.rst file...'

      open(21,file='jorek_perturbation.rst', form='unformatted', status='old', action='read', iostat=error)
      if (error .ne. 0) write(*,*) '...failed to open file jorek_perturbation.rst !'
      if (error .ne. 0) return
      
      read(21) n_tor_tmp_perturbation
      
      if (n_tor_tmp_perturbation .lt. 3) then
   	write(*,*)'The jorek_perturbation.rst file does not have n_tor>=3, required...'
   	return
      endif

      allocate( mode_tmp_perturbation  (n_tor_tmp_perturbation                ) )
      allocate( values_tmp_perturbation(n_tor_tmp_perturbation,n_order+1,n_var) )
      allocate( deltas_tmp_perturbation(n_tor_tmp_perturbation,n_order+1,n_var) )

      if (format_rst == 1) then
   	read(21) mode_tmp_perturbation
   	write(*,*) ' NEW format (1) : ',mode_tmp_perturbation
      elseif (format_rst == 0) then
   	write(*,*) ' mode : ',mode
   	if (n_tor .eq. n_tor_tmp) then 
   	  mode_tmp_perturbation = mode
   	else
   	  mode_tmp_perturbation(1:min(n_tor,n_tor_tmp_perturbation)) = mode(1:min(n_tor,n_tor_tmp_perturbation))
   	endif
   	write(*,*) ' OLD format (0) : '
   	write(*,'(A,999i4)') ' previous modenumbers : ',mode_tmp_perturbation
   	write(*,'(A,999i4)') ' new mode numbers     : ',mode
      elseif (format_rst > 2 ) then
   	write(*,'(A,i3)') ' restart file format not supported : ',format_rst
        stop
      endif

      allocate(node_list_perturbation, element_list_perturbation)
      
      read(21) node_list_perturbation%n_nodes,element_list_perturbation%n_elements
      read(21) node_list_perturbation%n_dof

      do i=1,node_list_perturbation%n_nodes

   	read(21) node_list_perturbation%node(i)%x

   	read(21) values_tmp_perturbation
   	read(21) deltas_tmp_perturbation

#ifdef fullmhd
        read(21) node_list_perturbation%node(i)%psi_eq               !< equilibrium flux at the nodes
        read(21) node_list_perturbation%node(i)%Fprof_eq             !< equilibrium profile R*B_phi at the nodes
#elif altcs
        read(21) node_list_perturbation%node(i)%psi_eq               !< equilibrium flux at the nodes
#endif
   	read(21) node_list_perturbation%node(i)%index
   	read(21) node_list_perturbation%node(i)%boundary
   	read(21) node_list_perturbation%node(i)%axis_node
   	read(21) node_list_perturbation%node(i)%parents
   	read(21) node_list_perturbation%node(i)%parent_elem
   	read(21) node_list_perturbation%node(i)%ref_lambda
   	read(21) node_list_perturbation%node(i)%ref_mu
   	read(21) node_list_perturbation%node(i)%constrained

   	node_list_perturbation%node(i)%values = 0.d0
   	node_list_perturbation%node(i)%deltas = 0.d0

   	do m=1,n_tor_tmp_perturbation,2
   	  do k=1, n_tor,2
   	    if (mode_tmp_perturbation(m) .eq. mode(k)) then
   	      if ((m .eq. 1) .and. (k.eq.1)) then
   		node_list_perturbation%node(i)%values(k,:,:)   = values_tmp_perturbation(m,:,:)
   		node_list_perturbation%node(i)%deltas(k,:,:)   = deltas_tmp_perturbation(m,:,:)
   	      else
   		node_list_perturbation%node(i)%values(k-1,:,:) = values_tmp_perturbation(m-1,:,:)
   		node_list_perturbation%node(i)%deltas(k-1,:,:) = deltas_tmp_perturbation(m-1,:,:)
   		node_list_perturbation%node(i)%values(k,:,:)   = values_tmp_perturbation(m,:,:)
   		node_list_perturbation%node(i)%deltas(k,:,:)   = deltas_tmp_perturbation(m,:,:)
   	      endif
   	    endif
   	  enddo
   	enddo

      enddo

      ! --- Import (n_tor,n_period) = (3,XX) mode only to another (n_tor,n_period) = (3,XX) equilibrium
      amplitude = 1.d0
      write(*,*)'Copying perturbation into node-structure using amplitude',amplitude
      do i=1,node_list%n_nodes
   	do j=1,n_var
   	  do k = 1, 4
   	    do m = 2, n_tor_tmp_perturbation
 	      node_list%node(i)%values(m,k,j) = amplitude*node_list_perturbation%node(i)%values(m,k,j)
   	      node_list%node(i)%deltas(m,k,j) = amplitude*node_list_perturbation%node(i)%deltas(m,k,j)
    	    enddo
   	  enddo
   	enddo
      enddo
      
      ! --- Deallocate temporary nodes/elements
      deallocate(node_list_perturbation)
      deallocate(element_list_perturbation)
      deallocate(mode_tmp_perturbation  )
      deallocate(values_tmp_perturbation)
      deallocate(deltas_tmp_perturbation)

  
    ! --- Using just noise
    else

      amplitude = 1.d-10
      do i=1,node_list%n_nodes
        node_list%node(i)%values(n_tor_tmp+1:n_tor,:,:)= 0.d0
        do j=n_tor_tmp+1, n_tor
          node_list%node(i)%values(j,:,var_rho)= amplitude * node_list%node(i)%values(1,:,var_rho)
#ifdef WITH_TiTe
          node_list%node(i)%values(j,:,var_Ti)= amplitude * node_list%node(i)%values(1,:,var_Ti)
          node_list%node(i)%values(j,:,var_Te)= amplitude * node_list%node(i)%values(1,:,var_Te)
#else
          node_list%node(i)%values(j,:,var_T)  = amplitude * node_list%node(i)%values(1,:,var_T)
#endif
#ifdef fullmhd
          node_list%node(i)%values(j,:,var_AR)= amplitude * node_list%node(i)%values(1,:,var_AR)
          node_list%node(i)%values(j,:,var_AZ)= amplitude * node_list%node(i)%values(1,:,var_AZ)
          node_list%node(i)%values(j,:,var_A3)= amplitude * node_list%node(i)%values(1,:,var_A3)
#endif
        enddo
      enddo
    endif
  endif

  ! End reading binary restart file  
  write(*,*) '********* end restart ******************'

  !call add_pellet(node_list,element_list,25.d0,0.06d0,0.03d0,3.78d0,0.14d0)

  ! -> Deallocate temporary arrays 
  if (allocated(mode_tmp))   call tr_deallocate(mode_tmp,"mode_tmp",CAT_UNKNOWN)
  if (allocated(values_tmp)) call tr_deallocate(values_tmp,"values_tmp",CAT_UNKNOWN)
  if (allocated(deltas_tmp)) call tr_deallocate(deltas_tmp,"deltas_tmp",CAT_UNKNOWN)

  call populate_element_rtree(node_list, element_list)

  return
end subroutine import_binary_restart


!
! Import an HDF5 restart file
subroutine import_hdf5_restart(node_list, element_list, filename, format_rst, error, no_perturbations)

#include "version.h"

  use tr_module 
  use data_structure
  use phys_module
  use pellet_module
  use vacuum, only: import_HDF5_restart_vacuum, current_FB_fact
  use mod_element_rtree, only: populate_element_rtree
#ifdef USE_HDF5
  use hdf5
  use hdf5_io_module
  use mod_parameters 
#endif
  
  implicit none
  
  ! --- Routine parameters
  type(type_node_list),    intent(inout) :: node_list
  type(type_element_list), intent(inout) :: element_list
  character(len=*),        intent(in)    :: filename
  integer,                 intent(in)    :: format_rst  ! format of restart file
  integer,                 intent(out)   :: error
  logical, optional,       intent(in)    :: no_perturbations ! don't initialize new harmonics
  
  ! --- Perturbation-Import variables
  type (type_node_list)   , pointer	:: node_list_perturbation
  type (type_element_list), pointer	:: element_list_perturbation
  integer              			:: n_tor_tmp_perturbation
  integer, allocatable 			:: mode_tmp_perturbation(:)
  real*8,  allocatable 			:: values_tmp_perturbation(:,:,:), deltas_tmp_perturbation(:,:,:)
  logical, parameter   			:: import_perturbation = .false.

  ! --- Local variables
  integer              :: i, j, m, k, n_tor_tmp, n_coord_tor_tmp, jorek_model_tmp, n_var_tmp, n_order_tmp, n_period_tmp, rst_hdf5_version_tmp
  integer              :: n_plane_tmp, n_vertex_max_tmp, n_nodes_max_tmp, n_elements_max_tmp,n_boundary_max_tmp
  integer              :: n_pieces_max_tmp, n_degrees_tmp, nref_max_tmp, n_ref_list_tmp, n_new_modes
  real*8               :: growth_mag, growth_kin, amplitude
  integer, allocatable :: mode_tmp(:), new_mode(:)
  real*8,  allocatable :: values_tmp(:,:,:), deltas_tmp(:,:,:)
  character*50         :: version_control, version_control_tmp
  logical              :: kept, modes_changed, import_3xx_4xx
  
#ifdef USE_HDF5
  integer(HID_T)     :: file_id, datatype, dataset
  integer            :: ind, n_spi_check, n_inj_check
  
  real(RKIND), allocatable :: t_x(:,:,:,:)
  real(RKIND), allocatable :: t_values(:,:,:,:)
  real(RKIND), allocatable :: t_deltas(:,:,:,:)

  real(RKIND), allocatable :: t_psi_eq(:,:)
  real(RKIND), allocatable :: t_Fprof_eq(:,:)

  integer,     allocatable :: t_index(:,:)
  integer,     allocatable :: t_boundary(:)
  character,   allocatable :: t_axis_node(:)     
  integer,     allocatable :: t_parents(:,:)
  integer,     allocatable :: t_parent_elem(:)
  real(RKIND), allocatable :: t_ref_lambda(:)
  real(RKIND), allocatable :: t_ref_mu(:)
  character,   allocatable :: t_constrained(:)     

  integer,     allocatable :: t_vertex(:,:)
  integer,     allocatable :: t_neighbours(:,:)
  real(RKIND), allocatable :: t_size(:,:,:)
  integer,     allocatable :: t_father(:)
  integer,     allocatable :: t_n_sons(:)
  integer,     allocatable :: t_n_gen(:)
  integer,     allocatable :: t_sons(:,:)
  integer,     allocatable :: t_contain_node(:,:)
  integer,     allocatable :: t_nref(:)

! local variables

  real*8, allocatable :: spi_R_arr (:)
  real*8, allocatable :: spi_Z_arr (:)
  real*8, allocatable :: spi_phi_arr (:)
  real*8, allocatable :: spi_phi_init_arr (:)
  real*8, allocatable :: spi_Vel_R_arr (:)
  real*8, allocatable :: spi_Vel_Z_arr (:)
  real*8, allocatable :: spi_Vel_RxZ_arr (:)
  real*8, allocatable :: spi_radius_arr (:)
  real*8, allocatable :: spi_abl_arr (:)
  real*8, allocatable :: spi_species_arr (:)
  integer, allocatable :: spi_species_arr_old (:)  !< For backward compatibility only

  integer :: err_exists, dterr, n_spi_begin, i_inj
  logical :: flag_exists, type_match

  real*8, allocatable :: t_energies(:,:,:)   !< Magnetic and kinetic mode energies at previous timesteps.
  real*8, allocatable :: t_energies2(:,:,:)  !< Magnetic and kinetic mode energies at previous timesteps.
  real*8, allocatable :: t_energies3(:,:,:)  !< Magnetic and kinetic mode energies at previous timesteps.
  real*8, allocatable :: t_energies4(:,:,:)  !< Magnetic and kinetic mode energies at previous timesteps.
  logical                               :: no_pert
  
  no_pert = .false.
  if ( present(no_perturbations) ) no_pert = no_perturbations
  
  !
#endif
  error = 0
#ifdef USE_HDF5

  ! ->  Reading HDF5 file
  write(*,*) 'Importing HDF5 restart file "', trim(filename), '".' 
  
  ! -> Open HDF5 file
  call HDF5_open(trim(filename),file_id,error)
  if ( error /= 0 ) then
    write(*,*) '...failed!'
    return
  end if

  ! Restart file version
  rst_hdf5_version_tmp = 0
  call HDF5_integer_reading(file_id,rst_hdf5_version_tmp,"rst_hdf5_version")
  write(*,*) 'Restart file has rst_hdf5_version=', rst_hdf5_version_tmp
  if ( rst_hdf5_version_tmp > rst_hdf5_version_supported ) then
    write(*,*) 'ERROR: Cannot read the hdf5 restart file "', trim(filename), '" since it was created with a more recent code version.'
    write(*,*) '* rst_hdf5_version of the restart file: ', rst_hdf5_version_tmp
    write(*,*) '* rst_hdf5_version of the code version: ', rst_hdf5_version_supported
    write(*,*) '* Note that you can export an older restart file version with your newer code by'
    write(*,*) '    explicitly setting rst_hdf5_version in the namelist input file.'
    stop
  end if

  call HDF5_char_reading(file_id,version_control_tmp, "RCS_version")
  version_control = trim(adjustl(RCS_VERSION))

  call HDF5_integer_reading(file_id,jorek_model_tmp,"jorek_model")
  call HDF5_integer_reading(file_id,n_var_tmp,"n_var")
  import_3xx_4xx = .false.
  if ( (jorek_model >= 400) .and. (jorek_model <= 499) .and. (jorek_model_tmp >= 300) .and. (jorek_model_tmp <= 399) ) then
    import_3xx_4xx = .true. ! Import a JOREK model 3XX restart file into a 4XX binary
    write(*,*) 'WARNING: Restarting a JOREK model 3XX simulation with model 4XX.'
  else if ( n_var /= n_var_tmp ) then
    write(*,*) 'ERROR: The number of variables in the restart file and the compiled JOREK binary does not agree.'
    write(*,*) '  Restarting normally works only with the same JOREK model.'
    write(*,*) '  As an exception, importing a 3XX restart file into model 4XX has been implemented.'
    stop
  end if
  call HDF5_integer_reading(file_id,n_order_tmp,"n_order")
  call HDF5_integer_reading(file_id,n_tor_tmp, "n_tor")
  if (rst_hdf5_version_tmp .eq. 2) then
    call HDF5_integer_reading(file_id,n_coord_tor_tmp, "n_coord_tor")
  else
    n_coord_tor_tmp = 1
  endif  
  call HDF5_integer_reading(file_id,n_period_tmp, "n_period")
  call HDF5_integer_reading(file_id,n_plane_tmp, "n_plane")
  call HDF5_integer_reading(file_id,n_vertex_max_tmp, "n_vertex_max")
  call HDF5_integer_reading(file_id,n_nodes_max_tmp, "n_nodes_max")
  call HDF5_integer_reading(file_id,n_elements_max_tmp, "n_elements_max")
  call HDF5_integer_reading(file_id,n_boundary_max_tmp, "n_boundary_max")
  call HDF5_integer_reading(file_id,n_pieces_max_tmp, "n_pieces_max")
  call HDF5_integer_reading(file_id,n_degrees_tmp, "n_degrees")
  call HDF5_integer_reading(file_id,nref_max_tmp, "nref_max")
  call HDF5_integer_reading(file_id,n_ref_list_tmp, "n_ref_list")


  if (allocated(mode_tmp))   call tr_deallocate(mode_tmp,"mode_tmp",CAT_UNKNOWN)
  allocate(mode_tmp(n_tor_tmp))
  mode_tmp = -1 ! unset

  if (format_rst == 1) then
    call HDF5_array1D_reading_int(file_id,mode_tmp,"mode_tmp")
    write(*,*) " import_restart, HDF5 file : n_var     = ",mode_tmp
    write(*,*) ' NEW format (1) : ',mode_tmp
  elseif (format_rst == 0) then
    do i=1, n_tor_tmp
       mode_tmp(i) = int(i / 2) * n_period_tmp
    end do
    modes_changed = .false.
    if (n_tor_tmp .ne. n_tor) then
      modes_changed = .true.
    elseif (sum(abs(mode_tmp-mode)) .gt. 0) then
      modes_changed = .true.
    end if
    
    write(*,*) ' OLD format (0) : '
    write(*,'(A,999i4)') ' previous modenumbers : ',mode_tmp
    write(*,'(A,999i4)') ' new mode numbers     : ',mode
    do i = 1, n_tor_tmp, 2
      kept = .false.
      do j = 1, n_tor, 2
        if ( mode_tmp(i) == mode(j) ) kept = .true.
      end do
      if ( .not. kept ) write (*,'(1x,a,i5,a)') 'Warning: The mode n=', mode_tmp(i), ' is being dropped!'
    end do
  elseif ( format_rst > 2 ) then
    write(*,'(A,i3)') ' restart file format not supported : ',format_rst
    stop
  endif

  if (n_tor_tmp .gt. n_tor) write(*,'(3(a,i5))') &
       ' Warning: Reducing number of harmonics from', n_tor_tmp, ' to', n_tor, '!'
  if (n_tor_tmp .lt. n_tor) write(*,'(3(a,i5))') &
       ' Warning: Increasing number of harmonics from', n_tor_tmp, ' to', n_tor, '!'
  if (n_period_tmp .ne. n_period) write(*,'(3(a,i5))') &
       ' Warning: n_period has changed from', n_period_tmp, ' to', n_period
  if (n_coord_tor_tmp .ne. n_coord_tor) then
    write(*,'(3(a,i5))') "Error: The number of toroidal harmonics in the grid representation has changed from ", &
                         n_coord_tor_tmp, " to ", n_coord_tor, "!"
    stop
  endif
    
  !write(*,'(2(A,i5))') ' Importing ',n_tor_tmp,' harmonics with n_period=', n_period_tmp 

  call HDF5_integer_reading(file_id,node_list%n_nodes,"n_nodes")
  call HDF5_integer_reading(file_id,element_list%n_elements,"n_elements")
  call HDF5_integer_reading(file_id,node_list%n_dof,"n_dof")

  ! -> Allocate temporary arrays 
  call tr_allocate(t_x,     1,node_list%n_nodes,1,n_coord_tor_tmp,1,n_order+1,1,n_dim,         "node_list%x",     CAT_UNKNOWN)
  call tr_allocate(t_values,1,node_list%n_nodes,1,      n_tor_tmp,1,n_order+1,1,n_var_tmp, "node_list%values",CAT_UNKNOWN)
  call tr_allocate(t_deltas,1,node_list%n_nodes,1,      n_tor_tmp,1,n_order+1,1,n_var_tmp, "node_list%deltas",CAT_UNKNOWN)
 
#ifdef fullmhd
  call tr_allocate(t_psi_eq,  1,node_list%n_nodes,1,n_order+1, "node_list%psi_eq",  CAT_UNKNOWN)
  call tr_allocate(t_Fprof_eq,1,node_list%n_nodes,1,n_order+1, "node_list%Fprof_eq",CAT_UNKNOWN)
#elif altcs
  call tr_allocate(t_psi_eq,  1,node_list%n_nodes,1,n_order+1, "node_list%psi_eq",  CAT_UNKNOWN)
#endif
 
  call tr_allocate(t_index,      1,node_list%n_nodes,1,n_order+1,"index",      CAT_UNKNOWN)
  call tr_allocate(t_boundary,   1,node_list%n_nodes,            "boundary",   CAT_UNKNOWN)
  call tr_allocate(t_axis_node,  1,node_list%n_nodes,            "axis_node",  CAT_UNKNOWN)
  call tr_allocate(t_parents,    1,node_list%n_nodes,1,2,        "parent",     CAT_UNKNOWN)
  call tr_allocate(t_parent_elem,1,node_list%n_nodes,            "parent_elem",CAT_UNKNOWN)
  call tr_allocate(t_ref_lambda, 1,node_list%n_nodes,            "ref_lambda" ,CAT_UNKNOWN)
  call tr_allocate(t_ref_mu,     1,node_list%n_nodes,            "ref_mu",     CAT_UNKNOWN)
  call tr_allocate(t_constrained,1,node_list%n_nodes,            "constrained",CAT_UNKNOWN)

  ! type_element, element_list%n_elements
  call tr_allocate(t_vertex,      1,element_list%n_elements,1,n_vertex_max,             "vertex",CAT_UNKNOWN)
  call tr_allocate(t_neighbours,  1,element_list%n_elements,1,n_vertex_max,             "neighbours",CAT_UNKNOWN)
  call tr_allocate(t_size,        1,element_list%n_elements,1,n_vertex_max,1,n_order+1, "size",CAT_UNKNOWN)
  call tr_allocate(t_father,      1,element_list%n_elements,                            "father",CAT_UNKNOWN)
  call tr_allocate(t_n_sons,      1,element_list%n_elements,                            "n_sons",CAT_UNKNOWN)
  call tr_allocate(t_n_gen,       1,element_list%n_elements,                            "n_gen",CAT_UNKNOWN)
  call tr_allocate(t_sons,        1,element_list%n_elements,1,4,                        "sons",CAT_UNKNOWN)
  call tr_allocate(t_contain_node,1,element_list%n_elements,1,5,                        "contain_node",CAT_UNKNOWN)
  call tr_allocate(t_nref,        1,element_list%n_elements,                            "nref",CAT_UNKNOWN)

  if (rst_hdf5_version .eq. 2) then
    call HDF5_array4D_reading(file_id,t_x,        'x')
  else
    call HDF5_array3D_reading(file_id,t_x(:,1,:,:),        'x')
  endif
  call HDF5_array4D_reading(file_id,t_values,   'values')
  call HDF5_array4D_reading(file_id,t_deltas,   'deltas')

#ifdef fullmhd
  call HDF5_array2D_reading(file_id,t_psi_eq,   'psi_eq')
  call HDF5_array2D_reading(file_id,t_Fprof_eq, 'Fprof_eq')
#elif altcs
  call HDF5_array2D_reading(file_id,t_psi_eq,   'psi_eq')
#endif

  call HDF5_array2D_reading_int (file_id,t_index,       'index')
  call HDF5_array1D_reading_int (file_id,t_boundary,    'boundary')
  call HDF5_array1D_reading_char(file_id,t_axis_node,   'axis_node')
  call HDF5_array2D_reading_int (file_id,t_parents,     'parents')
  call HDF5_array1D_reading_int (file_id,t_parent_elem, 'parent_elem')
  call HDF5_array1D_reading     (file_id,t_ref_lambda,  'ref_lambda')
  call HDF5_array1D_reading     (file_id,t_ref_mu,      'ref_mu')
  call HDF5_array1D_reading_char(file_id,t_constrained, 'constrained')

  ! --- Detect new modes that need to be initialized to noise level
  if (allocated(new_mode))   call tr_deallocate(new_mode,"new_mode",CAT_UNKNOWN)
  allocate(new_mode(n_tor))
  new_mode(:)=1

  do m=1,n_tor_tmp,2
    do k=1, n_tor,2 
      if (mode_tmp(m) .eq. mode(k)) then
        if ((m .eq. 1) .and. (k.eq.1)) then
          new_mode(k)=0
        else
          new_mode(k-1)=0
          new_mode(k)=0
        end if
      end if
    end do
  end do
  if (any(new_mode .ne. 0)) write(*,'(a,999i4)') ' need initialization  : ', new_mode
  
  do i=1,node_list%n_nodes
    node_list%node(i)%x = t_x(i,:,:,:) 

    node_list%node(i)%values = 0.d0 
    node_list%node(i)%deltas = 0.d0 

    do m=1,n_tor_tmp,2
      do k=1, n_tor,2 
        if (mode_tmp(m) .eq. mode(k)) then
          if ((m .eq. 1) .and. (k.eq.1)) then
            node_list%node(i)%values(k,:,1:n_var_tmp)   = t_values(i,m,:,1:n_var_tmp)
            node_list%node(i)%deltas(k,:,1:n_var_tmp)   = t_deltas(i,m,:,1:n_var_tmp)
          else
            node_list%node(i)%values(k-1,:,1:n_var_tmp) = t_values(i,m-1,:,1:n_var_tmp)
            node_list%node(i)%deltas(k-1,:,1:n_var_tmp) = t_deltas(i,m-1,:,1:n_var_tmp) 
            node_list%node(i)%values(k,:,1:n_var_tmp)   = t_values(i,m,:,1:n_var_tmp) 
            node_list%node(i)%deltas(k,:,1:n_var_tmp)   = t_deltas(i,m,:,1:n_var_tmp)
          end if
        end if
      end do
    end do

    ! --- Split "total" temperature into electron and ion temperature
    if ( import_3xx_4xx ) then
      node_list%node(i)%values(:,:,var_Te) = node_list%node(i)%values(:,:,6) / 2.d0
      node_list%node(i)%deltas(:,:,var_Te) = node_list%node(i)%deltas(:,:,6) / 2.d0
      node_list%node(i)%values(:,:,var_Ti) = node_list%node(i)%values(:,:,6) / 2.d0
      node_list%node(i)%deltas(:,:,var_Ti) = node_list%node(i)%deltas(:,:,6) / 2.d0
    end if


#ifdef fullmhd
    node_list%node(i)%psi_eq   = t_psi_eq(i,:)
    node_list%node(i)%Fprof_eq = t_Fprof_eq(i,:)
#elif altcs
    node_list%node(i)%psi_eq   = t_psi_eq(i,:)
#endif

    node_list%node(i)%index = t_index(i,:)
    node_list%node(i)%boundary = t_boundary(i)
    if (t_axis_node(i) == 'T') then
       node_list%node(i)%axis_node = .true.
    else
       node_list%node(i)%axis_node = .false.
    end if
    node_list%node(i)%parents = t_parents(i,:)
    node_list%node(i)%parent_elem = t_parent_elem(i)
    node_list%node(i)%ref_lambda = t_ref_lambda(i)
    node_list%node(i)%ref_mu = t_ref_mu(i)
    if (t_constrained(i) == 'T') then
       node_list%node(i)%constrained = .true.
    else
       node_list%node(i)%constrained = .false.
    end if
  end do

  call HDF5_array2D_reading_int(file_id,t_vertex,      'vertex')
  call HDF5_array2D_reading_int(file_id,t_neighbours,  'neighbours')
  call HDF5_array3D_reading    (file_id,t_size,        'size')
  call HDF5_array1D_reading_int(file_id,t_father,      'father')
  call HDF5_array1D_reading_int(file_id,t_n_sons,      'n_sons')
  call HDF5_array1D_reading_int(file_id,t_n_gen,       'n_gen')
  call HDF5_array2D_reading_int(file_id,t_sons,        'sons')
  call HDF5_array2D_reading_int(file_id,t_contain_node,'contain_node')
  call HDF5_array1D_reading_int(file_id,t_nref,        'nref')

  do i=1,element_list%n_elements
    element_list%element(i)%vertex	 = t_vertex(i,:)
    element_list%element(i)%neighbours   = t_neighbours(i,:)
    element_list%element(i)%size	 = t_size(i,:,:)
    element_list%element(i)%father	 = t_father(i)
    element_list%element(i)%n_sons	 = t_n_sons(i)
    element_list%element(i)%n_gen	 = t_n_gen(i)
    element_list%element(i)%sons	 = t_sons(i,:)
    element_list%element(i)%contain_node = t_contain_node(i,:)
    element_list%element(i)%nref	 = t_nref(i)
  end do
 
  call HDF5_real_reading(file_id,tstep,'tstep')
  call HDF5_real_reading(file_id,eta_rst,'eta')
  call HDF5_real_reading(file_id,visco_rst,'visco')
  call HDF5_real_reading(file_id,visco_par_rst,'visco_par')
  call HDF5_integer_reading(file_id,index_start,'index_now')
  call HDF5_real_reading(file_id,t_start,'t_now')
  
  if (index_start .ge. 1) then

    if (allocated(xtime)) call tr_deallocate(xtime,"xtime",CAT_UNKNOWN)
    call tr_allocate(xtime,1,index_start+nstep,"xtime",CAT_UNKNOWN)
    call HDF5_array1D_reading(file_id,xtime,'xtime')

    if (allocated(t_energies))   call tr_deallocate(t_energies,"t_energies",CAT_UNKNOWN)
    call tr_allocate(t_energies,1,n_tor_tmp,1,2,1,index_start+nstep,"t_energies",CAT_UNKNOWN)
    t_energies = 0.d0
    call HDF5_array3D_reading(file_id,t_energies,'energies')

    if (allocated(energies))   call tr_deallocate(energies,"energies",CAT_UNKNOWN)
    call tr_allocate(energies,1,n_tor,1,2,1,index_start+nstep,"energies",CAT_UNKNOWN)
    energies = 0.d0

    do m=1,n_tor_tmp,2
      do k=1, n_tor,2 
        if (mode_tmp(m) .eq. mode(k)) then
          if ((m .eq. 1) .and. (k.eq.1)) then
            energies(k,:,:) = t_energies(m,:,:)
          else
            energies(k-1:k,:,:) = t_energies(m-1:m,:,:)
          end if
        end if
      end do
    end do

    if (allocated(R_axis_t)) call tr_deallocate(R_axis_t,"R_axis_t",CAT_UNKNOWN)
    call tr_allocate(R_axis_t,1,index_start+nstep,"R_axis_t",CAT_UNKNOWN)
    R_axis_t = 0.d0
    call HDF5_array1D_reading(file_id,R_axis_t,'R_axis_t')
    
    if (allocated(Z_axis_t)) call tr_deallocate(Z_axis_t,"Z_axis_t",CAT_UNKNOWN)
    call tr_allocate(Z_axis_t,1,index_start+nstep,"Z_axis_t",CAT_UNKNOWN)
    Z_axis_t = 0.d0
    call HDF5_array1D_reading(file_id,Z_axis_t,'Z_axis_t')
    
    if (allocated(psi_axis_t)) call tr_deallocate(psi_axis_t,"psi_axis_t",CAT_UNKNOWN)
    call tr_allocate(psi_axis_t,1,index_start+nstep,"psi_axis_t",CAT_UNKNOWN)
    psi_axis_t = 0.d0
    call HDF5_array1D_reading(file_id,psi_axis_t,'psi_axis_t')
    
    if (allocated(R_xpoint_t)) call tr_deallocate(R_xpoint_t,"R_xpoint_t",CAT_UNKNOWN)
    call tr_allocate(R_xpoint_t,1,index_start+nstep,1,2,"R_xpoint_t",CAT_UNKNOWN)
    R_xpoint_t = 0.d0
    call HDF5_array2D_reading(file_id,R_xpoint_t,'R_xpoint_t')
    
    if (allocated(Z_xpoint_t)) call tr_deallocate(Z_xpoint_t,"Z_xpoint_t",CAT_UNKNOWN)
    call tr_allocate(Z_xpoint_t,1,index_start+nstep,1,2,"Z_xpoint_t",CAT_UNKNOWN)
    Z_xpoint_t = 0.d0
    call HDF5_array2D_reading(file_id,Z_xpoint_t,'Z_xpoint_t')

    if (allocated(psi_xpoint_t)) call tr_deallocate(psi_xpoint_t,"psi_xpoint_t",CAT_UNKNOWN)
    call tr_allocate(psi_xpoint_t,1,index_start+nstep,1,2,"psi_xpoint_t",CAT_UNKNOWN)
    psi_xpoint_t = 0.d0
    call HDF5_array2D_reading(file_id,psi_xpoint_t,'psi_xpoint_t')

    if (allocated(R_bnd_t)) call tr_deallocate(R_bnd_t,"R_bnd_t",CAT_UNKNOWN)
    call tr_allocate(R_bnd_t,1,index_start+nstep,"R_bnd_t",CAT_UNKNOWN)
    R_bnd_t = 0.d0
    call HDF5_array1D_reading(file_id,R_bnd_t,'R_bnd_t')
    
    if (allocated(Z_bnd_t)) call tr_deallocate(Z_bnd_t,"Z_bnd_t",CAT_UNKNOWN)
    call tr_allocate(Z_bnd_t,1,index_start+nstep,"Z_bnd_t",CAT_UNKNOWN)
    Z_bnd_t = 0.d0
    call HDF5_array1D_reading(file_id,Z_bnd_t,'Z_bnd_t')
    
    if (allocated(psi_bnd_t)) call tr_deallocate(psi_bnd_t,"psi_bnd_t",CAT_UNKNOWN)
    call tr_allocate(psi_bnd_t,1,index_start+nstep,"psi_bnd_t",CAT_UNKNOWN)
    psi_bnd_t = 0.d0
    call HDF5_array1D_reading(file_id,psi_bnd_t,'psi_bnd_t')
    
    if (allocated(current_t)) call tr_deallocate(current_t,"current_t",CAT_UNKNOWN)
    call tr_allocate(current_t,1,index_start+nstep,"current_t",CAT_UNKNOWN)
    current_t = 0.d0
    call HDF5_array1D_reading(file_id,current_t,'current_t')
    
    if (allocated(beta_p_t)) call tr_deallocate(beta_p_t,"beta_p_t",CAT_UNKNOWN)
    call tr_allocate(beta_p_t,1,index_start+nstep,"beta_p_t",CAT_UNKNOWN)
    beta_p_t = 0.d0
    call HDF5_array1D_reading(file_id,beta_p_t,'beta_p_t')
    
    if (allocated(beta_t_t)) call tr_deallocate(beta_t_t,"beta_t_t",CAT_UNKNOWN)
    call tr_allocate(beta_t_t,1,index_start+nstep,"beta_t_t",CAT_UNKNOWN)
    beta_t_t = 0.d0
    call HDF5_array1D_reading(file_id,beta_t_t,'beta_t_t')
    
    if (allocated(beta_n_t)) call tr_deallocate(beta_n_t,"beta_n_t",CAT_UNKNOWN)
    call tr_allocate(beta_n_t,1,index_start+nstep,"beta_n_t",CAT_UNKNOWN)
    beta_n_t = 0.d0
    call HDF5_array1D_reading(file_id,beta_n_t,'beta_n_t')
    
    if (allocated(density_in_t)) call tr_deallocate(density_in_t,"density_in_t",CAT_UNKNOWN)
    call tr_allocate(density_in_t,1,index_start+nstep,"density_in_t",CAT_UNKNOWN)
    density_in_t = 0.d0
    call HDF5_array1D_reading(file_id,density_in_t,'density_in_t')
    
    if (allocated(density_out_t)) call tr_deallocate(density_out_t,"density_out_t",CAT_UNKNOWN)
    call tr_allocate(density_out_t,1,index_start+nstep,"density_out_t",CAT_UNKNOWN)
    density_out_t = 0.d0
    call HDF5_array1D_reading(file_id,density_out_t,'density_out_t')
    
    if (allocated(pressure_in_t)) call tr_deallocate(pressure_in_t,"pressure_in_t",CAT_UNKNOWN)
    call tr_allocate(pressure_in_t,1,index_start+nstep,"pressure_in_t",CAT_UNKNOWN)
    pressure_in_t = 0.d0
    call HDF5_array1D_reading(file_id,pressure_in_t,'pressure_in_t')
    
    if (allocated(pressure_out_t)) call tr_deallocate(pressure_out_t,"pressure_out_t",CAT_UNKNOWN)
    call tr_allocate(pressure_out_t,1,index_start+nstep,"pressure_out_t",CAT_UNKNOWN)
    pressure_out_t = 0.d0
    call HDF5_array1D_reading(file_id,pressure_out_t,'pressure_out_t')
    
    if (allocated(heat_src_in_t)) call tr_deallocate(heat_src_in_t,"heating_power_t",CAT_UNKNOWN)
    call tr_allocate(heat_src_in_t,1,index_start+nstep,"heat_src_in_t",CAT_UNKNOWN)
    heat_src_in_t = 0.d0
    call HDF5_array1D_reading(file_id,heat_src_in_t,'heat_src_in_t')
    
    if (allocated(heat_src_out_t)) call tr_deallocate(heat_src_out_t,"heating_power_t",CAT_UNKNOWN)
    call tr_allocate(heat_src_out_t,1,index_start+nstep,"heat_src_out_t",CAT_UNKNOWN)
    heat_src_out_t = 0.d0
    call HDF5_array1D_reading(file_id,heat_src_out_t,'heat_src_out_t')
    
    if (allocated(part_src_in_t)) call tr_deallocate(part_src_in_t,"parting_power_t",CAT_UNKNOWN)
    call tr_allocate(part_src_in_t,1,index_start+nstep,"part_src_in_t",CAT_UNKNOWN)
    part_src_in_t = 0.d0
    call HDF5_array1D_reading(file_id,part_src_in_t,'part_src_in_t')
    
    if (allocated(part_src_out_t)) call tr_deallocate(part_src_out_t,"parting_power_t",CAT_UNKNOWN)
    call tr_allocate(part_src_out_t,1,index_start+nstep,"part_src_out_t",CAT_UNKNOWN)
    part_src_out_t = 0.d0
    call HDF5_array1D_reading(file_id,part_src_out_t,'part_src_out_t')
  
    if (allocated(E_tot_t)) call tr_deallocate(E_tot_t,"E_tot_t",CAT_UNKNOWN)
    call tr_allocate(E_tot_t,1,index_start+nstep,"E_tot_t",CAT_UNKNOWN)
    E_tot_t = 0.d0
    call HDF5_array1D_reading(file_id,E_tot_t,'E_tot_t')

    if (allocated(helicity_tot_t)) call tr_deallocate(helicity_tot_t,"helicity_tot_t",CAT_UNKNOWN)
    call tr_allocate(helicity_tot_t,1,index_start+nstep,"helicity_tot_t",CAT_UNKNOWN)
    helicity_tot_t = 0.d0
    call HDF5_array1D_reading(file_id,helicity_tot_t,'helicity_tot_t')

    if (allocated(thermal_tot_t)) call tr_deallocate(thermal_tot_t,"thermal_tot_t",CAT_UNKNOWN)
    call tr_allocate(thermal_tot_t,1,index_start+nstep,"thermal_tot_t",CAT_UNKNOWN)
    thermal_tot_t = 0.d0
    call HDF5_array1D_reading(file_id,thermal_tot_t,'thermal_tot_t')

    if (allocated(thermal_e_tot_t)) call tr_deallocate(thermal_e_tot_t,"thermal_e_tot_t",CAT_UNKNOWN)
    call tr_allocate(thermal_e_tot_t,1,index_start+nstep,"thermal_e_tot_t",CAT_UNKNOWN)
    thermal_e_tot_t = 0.d0
    call HDF5_array1D_reading(file_id,thermal_e_tot_t,'thermal_e_tot_t')

    if (allocated(thermal_i_tot_t)) call tr_deallocate(thermal_i_tot_t,"thermal_i_tot_t",CAT_UNKNOWN)
    call tr_allocate(thermal_i_tot_t,1,index_start+nstep,"thermal_i_tot_t",CAT_UNKNOWN)
    thermal_i_tot_t = 0.d0
    call HDF5_array1D_reading(file_id,thermal_i_tot_t,'thermal_i_tot_t')

    if (allocated(kin_par_tot_t)) call tr_deallocate(kin_par_tot_t,"kin_par_tot_t",CAT_UNKNOWN)
    call tr_allocate(kin_par_tot_t,1,index_start+nstep,"kin_par_tot_t",CAT_UNKNOWN)
    kin_par_tot_t = 0.d0
    call HDF5_array1D_reading(file_id,kin_par_tot_t,'kin_par_tot_t')

    if (allocated(kin_perp_tot_t)) call tr_deallocate(kin_perp_tot_t,"kin_perp_tot_t",CAT_UNKNOWN)
    call tr_allocate(kin_perp_tot_t,1,index_start+nstep,"kin_perp_tot_t",CAT_UNKNOWN)
    kin_perp_tot_t = 0.d0
    call HDF5_array1D_reading(file_id,kin_perp_tot_t,'kin_perp_tot_t')

    if (allocated(Ip_tot_t)) call tr_deallocate(Ip_tot_t,"Ip_tot_t",CAT_UNKNOWN)
    call tr_allocate(Ip_tot_t,1,index_start+nstep,"Ip_tot_t",CAT_UNKNOWN)
    Ip_tot_t = 0.d0
    call HDF5_array1D_reading(file_id,Ip_tot_t,'Ip_tot_t')

    if (allocated(ohmic_tot_t)) call tr_deallocate(ohmic_tot_t,"ohmic_tot_t",CAT_UNKNOWN)
    call tr_allocate(ohmic_tot_t,1,index_start+nstep,"ohmic_tot_t",CAT_UNKNOWN)
    ohmic_tot_t = 0.d0
    call HDF5_array1D_reading(file_id,ohmic_tot_t,'ohmic_tot_t')

    if (allocated(Wmag_tot_t)) call tr_deallocate(Wmag_tot_t,"Wmag_tot_t",CAT_UNKNOWN)
    call tr_allocate(Wmag_tot_t,1,index_start+nstep,"Wmag_tot_t",CAT_UNKNOWN)
    Wmag_tot_t = 0.d0
    call HDF5_array1D_reading(file_id,Wmag_tot_t,'Wmag_tot_t')

    if (allocated(Magwork_tot_t)) call tr_deallocate(Magwork_tot_t,"Magwork_tot_t",CAT_UNKNOWN)
    call tr_allocate(Magwork_tot_t,1,index_start+nstep,"Magwork_tot_t",CAT_UNKNOWN)
    Magwork_tot_t = 0.d0
    call HDF5_array1D_reading(file_id,Magwork_tot_t,'Magwork_tot_t')

    if (allocated(flux_qpar_t)) call tr_deallocate(flux_qpar_t,"flux_qpar_t",CAT_UNKNOWN)
    call tr_allocate(flux_qpar_t,1,index_start+nstep,"flux_qpar_t",CAT_UNKNOWN)
    flux_qpar_t = 0.d0
    call HDF5_array1D_reading(file_id,flux_qpar_t,'flux_qpar_t')

    if (allocated(flux_qperp_t)) call tr_deallocate(flux_qperp_t,"flux_qperp_t",CAT_UNKNOWN)
    call tr_allocate(flux_qperp_t,1,index_start+nstep,"flux_qperp_t",CAT_UNKNOWN)
    flux_qperp_t = 0.d0
    call HDF5_array1D_reading(file_id,flux_qperp_t,'flux_qperp_t')

    if (allocated(flux_kinpar_t)) call tr_deallocate(flux_kinpar_t,"flux_kinpar_t",CAT_UNKNOWN)
    call tr_allocate(flux_kinpar_t,1,index_start+nstep,"flux_kinpar_t",CAT_UNKNOWN)
    flux_kinpar_t = 0.d0
    call HDF5_array1D_reading(file_id,flux_kinpar_t,'flux_kinpar_t')

    if (allocated(flux_poynting_t)) call tr_deallocate(flux_poynting_t,"flux_poynting_t",CAT_UNKNOWN)
    call tr_allocate(flux_poynting_t,1,index_start+nstep,"flux_poynting_t",CAT_UNKNOWN)
    flux_poynting_t = 0.d0
    call HDF5_array1D_reading(file_id,flux_poynting_t,'flux_poynting_t')

    if (allocated(flux_Pvn_t)) call tr_deallocate(flux_Pvn_t,"flux_Pvn_t",CAT_UNKNOWN)
    call tr_allocate(flux_Pvn_t,1,index_start+nstep,"flux_Pvn_t",CAT_UNKNOWN)
    flux_Pvn_t = 0.d0
    call HDF5_array1D_reading(file_id,flux_Pvn_t,'flux_Pvn_t')

    if (allocated(dE_tot_dt)) call tr_deallocate(dE_tot_dt,"dE_tot_dt",CAT_UNKNOWN)
    call tr_allocate(dE_tot_dt,1,index_start+nstep,"dE_tot_dt",CAT_UNKNOWN)
    dE_tot_dt = 0.d0
    call HDF5_array1D_reading(file_id,dE_tot_dt,'dE_tot_dt')

    if (allocated(dWmag_tot_dt)) call tr_deallocate(dWmag_tot_dt,"dWmag_tot_dt",CAT_UNKNOWN)
    call tr_allocate(dWmag_tot_dt,1,index_start+nstep,"dWmag_tot_dt",CAT_UNKNOWN)
    dWmag_tot_dt = 0.d0
    call HDF5_array1D_reading(file_id,dWmag_tot_dt,'dWmag_tot_dt')

    if (allocated(dthermal_tot_dt)) call tr_deallocate(dthermal_tot_dt,"dthermal_tot_dt",CAT_UNKNOWN)
    call tr_allocate(dthermal_tot_dt,1,index_start+nstep,"dthermal_tot_dt",CAT_UNKNOWN)
    dthermal_tot_dt = 0.d0
    call HDF5_array1D_reading(file_id,dthermal_tot_dt,'dthermal_tot_dt')

    if (allocated(dkinperp_tot_dt)) call tr_deallocate(dkinperp_tot_dt,"dkinperp_tot_dt",CAT_UNKNOWN)
    call tr_allocate(dkinperp_tot_dt,1,index_start+nstep,"dkinperp_tot_dt",CAT_UNKNOWN)
    dkinperp_tot_dt = 0.d0
    call HDF5_array1D_reading(file_id,dkinperp_tot_dt,'dkinperp_tot_dt')

    if (allocated(dkinpar_tot_dt)) call tr_deallocate(dkinpar_tot_dt,"dkinpar_tot_dt",CAT_UNKNOWN)
    call tr_allocate(dkinpar_tot_dt,1,index_start+nstep,"dkinpar_tot_dt",CAT_UNKNOWN)
    dkinpar_tot_dt = 0.d0
    call HDF5_array1D_reading(file_id,dkinpar_tot_dt,'dkinpar_tot_dt')

    if (allocated(heat_src_tot_t)) call tr_deallocate(heat_src_tot_t,"heat_src_tot_t",CAT_UNKNOWN)
    call tr_allocate(heat_src_tot_t,1,index_start+nstep,"heat_src_tot_t",CAT_UNKNOWN)
    heat_src_tot_t = 0.d0
    call HDF5_array1D_reading(file_id,heat_src_tot_t,'heat_src_tot_t')

    if (allocated(part_src_tot_t)) call tr_deallocate(part_src_tot_t,"part_src_tot_t",CAT_UNKNOWN)
    call tr_allocate(part_src_tot_t,1,index_start+nstep,"part_src_tot_t",CAT_UNKNOWN)
    part_src_tot_t = 0.d0
    call HDF5_array1D_reading(file_id,part_src_tot_t,'part_src_tot_t')

    if (allocated(li3_tot_t)) call tr_deallocate(li3_tot_t,"li3_tot_t",CAT_UNKNOWN)
    call tr_allocate(li3_tot_t,1,index_start+nstep,"li3_tot_t",CAT_UNKNOWN)
    li3_tot_t = 0.d0
    call HDF5_array1D_reading(file_id,li3_tot_t,'li3_tot_t')

    if (allocated(li3_t)) call tr_deallocate(li3_t,"li3_t",CAT_UNKNOWN)
    call tr_allocate(li3_t,1,index_start+nstep,"li3_t",CAT_UNKNOWN)
    li3_t = 0.d0
    call HDF5_array1D_reading(file_id,li3_t,'li3_t')

    if (allocated(viscopar_flux_t)) call tr_deallocate(viscopar_flux_t,"viscopar_flux_t",CAT_UNKNOWN)
    call tr_allocate(viscopar_flux_t,1,index_start+nstep,"viscopar_flux_t",CAT_UNKNOWN)
    viscopar_flux_t = 0.d0
    call HDF5_array1D_reading(file_id,viscopar_flux_t,'viscopar_flux_t')

    if (allocated(viscopar_dissip_tot_t)) call tr_deallocate(viscopar_dissip_tot_t,"viscopar_dissip_tot_t",CAT_UNKNOWN)
    call tr_allocate(viscopar_dissip_tot_t,1,index_start+nstep,"viscopar_dissip_tot_t",CAT_UNKNOWN)
    viscopar_dissip_tot_t = 0.d0
    call HDF5_array1D_reading(file_id,viscopar_dissip_tot_t,'viscopar_dissip_tot_t')

    if (allocated(friction_dissip_tot_t)) call tr_deallocate(friction_dissip_tot_t,"friction_dissip_tot_t",CAT_UNKNOWN)
    call tr_allocate(friction_dissip_tot_t,1,index_start+nstep,"friction_dissip_tot_t",CAT_UNKNOWN)
    friction_dissip_tot_t = 0.d0
    call HDF5_array1D_reading(file_id,friction_dissip_tot_t,'friction_dissip_tot_t')

    if (allocated(thmwork_tot_t)) call tr_deallocate(thmwork_tot_t,"thmwork_tot_t",CAT_UNKNOWN)
    call tr_allocate(thmwork_tot_t,1,index_start+nstep,"thmwork_tot_t",CAT_UNKNOWN)
    thmwork_tot_t = 0.d0
    call HDF5_array1D_reading(file_id,thmwork_tot_t,'thmwork_tot_t')

    if (allocated(volume_t)) call tr_deallocate(volume_t,"volume_t",CAT_UNKNOWN)
    call tr_allocate(volume_t,1,index_start+nstep,"volume_t",CAT_UNKNOWN)
    volume_t = 0.d0
    call HDF5_array1D_reading(file_id,volume_t,'volume_t')

    if (allocated(area_t)) call tr_deallocate(area_t,"area_t",CAT_UNKNOWN)
    call tr_allocate(area_t,1,index_start+nstep,"area_t",CAT_UNKNOWN)
    area_t = 0.d0
    call HDF5_array1D_reading(file_id,area_t,'area_t')

    if (allocated(mag_ener_src_tot)) call tr_deallocate(mag_ener_src_tot,"mag_ener_src_tot",CAT_UNKNOWN)
    call tr_allocate(mag_ener_src_tot,1,index_start+nstep,"mag_ener_src_tot",CAT_UNKNOWN)
    mag_ener_src_tot = 0.d0
    call HDF5_array1D_reading(file_id,mag_ener_src_tot,'mag_ener_src_tot')

    if (allocated(part_flux_Dpar_t)) call tr_deallocate(part_flux_Dpar_t,"part_flux_Dpar_t",CAT_UNKNOWN)
    call tr_allocate(part_flux_Dpar_t,1,index_start+nstep,"part_flux_Dpar_t",CAT_UNKNOWN)
    part_flux_Dpar_t = 0.d0
    call HDF5_array1D_reading(file_id,part_flux_Dpar_t,'part_flux_Dpar_t')

    if (allocated(part_flux_Dperp_t)) call tr_deallocate(part_flux_Dperp_t,"part_flux_Dperp_t",CAT_UNKNOWN)
    call tr_allocate(part_flux_Dperp_t,1,index_start+nstep,"part_flux_Dperp_t",CAT_UNKNOWN)
    part_flux_Dperp_t = 0.d0
    call HDF5_array1D_reading(file_id,part_flux_Dperp_t,'part_flux_Dperp_t')

    if (allocated(part_flux_Vpar_t)) call tr_deallocate(part_flux_Vpar_t,"part_flux_Vpar_t",CAT_UNKNOWN)
    call tr_allocate(part_flux_Vpar_t,1,index_start+nstep,"part_flux_Vpar_t",CAT_UNKNOWN)
    part_flux_Vpar_t = 0.d0
    call HDF5_array1D_reading(file_id,part_flux_Vpar_t,'part_flux_Vpar_t')

    if (allocated(part_flux_Vperp_t)) call tr_deallocate(part_flux_Vperp_t,"part_flux_Vperp_t",CAT_UNKNOWN)
    call tr_allocate(part_flux_Vperp_t,1,index_start+nstep,"part_flux_Vperp_t",CAT_UNKNOWN)
    part_flux_Vperp_t = 0.d0
    call HDF5_array1D_reading(file_id,part_flux_Vperp_t,'part_flux_Vperp_t')

    if (allocated(npart_flux_t)) call tr_deallocate(npart_flux_t,"npart_flux_t",CAT_UNKNOWN)
    call tr_allocate(npart_flux_t,1,index_start+nstep,"npart_flux_t",CAT_UNKNOWN)
    npart_flux_t = 0.d0
    call HDF5_array1D_reading(file_id,npart_flux_t,'npart_flux_t')

    if (allocated(dpart_tot_dt)) call tr_deallocate(dpart_tot_dt,"dpart_tot_dt",CAT_UNKNOWN)
    call tr_allocate(dpart_tot_dt,1,index_start+nstep,"dpart_tot_dt",CAT_UNKNOWN)
    dpart_tot_dt = 0.d0
    call HDF5_array1D_reading(file_id,dpart_tot_dt,'dpart_tot_dt')

    if (allocated(dnpart_tot_dt)) call tr_deallocate(dnpart_tot_dt,"dnpart_tot_dt",CAT_UNKNOWN)
    call tr_allocate(dnpart_tot_dt,1,index_start+nstep,"dnpart_tot_dt",CAT_UNKNOWN)
    dnpart_tot_dt = 0.d0
    call HDF5_array1D_reading(file_id,dnpart_tot_dt,'dnpart_tot_dt')

    if (allocated(npart_tot_t)) call tr_deallocate(npart_tot_t,"npart_tot_t",CAT_UNKNOWN)
    call tr_allocate(npart_tot_t,1,index_start+nstep,"npart_tot_t",CAT_UNKNOWN)
    npart_tot_t = 0.d0
    call HDF5_array1D_reading(file_id,npart_tot_t,'npart_tot_t')

    if (allocated(density_tot_t)) call tr_deallocate(density_tot_t,"density_tot_t",CAT_UNKNOWN)
    call tr_allocate(density_tot_t,1,index_start+nstep,"density_tot_t",CAT_UNKNOWN)
    density_tot_t = 0.d0
    call HDF5_array1D_reading(file_id,density_tot_t,'density_tot_t')

#ifdef JECCD                   
    if (allocated(t_energies2))   call tr_deallocate(t_energies2,"t_energies2",CAT_UNKNOWN)
    call tr_allocate(t_energies2,1,n_tor_tmp,1,2,1,index_start+nstep, "t_energies2",CAT_UNKNOWN)
    if (allocated(t_energies3))   call tr_deallocate(t_energies3,"t_energies3",CAT_UNKNOWN)
    call tr_allocate(t_energies3,1,n_tor_tmp,1,2,1,index_start+nstep, "t_energies3",CAT_UNKNOWN)
    t_energies2 = 0.d0
    t_energies3 = 0.d0
    call HDF5_array3D_reading(file_id,t_energies2,'energies2')
    call HDF5_array3D_reading(file_id,t_energies3,'energies3')

    if (allocated(energies2))   call tr_deallocate(energies2,"energies2",CAT_UNKNOWN)
    call tr_allocate(energies2,1,n_tor,1,2,1,index_start+nstep, "energies2",CAT_UNKNOWN)
    if (allocated(energies3))   call tr_deallocate(energies3,"energies3",CAT_UNKNOWN)
    call tr_allocate(energies3,1,n_tor,1,2,1,index_start+nstep, "energies3",CAT_UNKNOWN)
    energies2 = 0.d0
    energies3 = 0.d0

    do m=1,n_tor_tmp,2
      do k=1, n_tor,2 
        if (mode_tmp(m) .eq. mode(k)) then
          if ((m .eq. 1) .and. (k.eq.1)) then
            energies2(k,:,:) = t_energies2(m,:,:)
            energies3(k,:,:) = t_energies3(m,:,:)
          else
            energies2(k-1:k,:,:) = t_energies2(m-1:m,:,:)
            energies3(k-1:k,:,:) = t_energies3(m-1:m,:,:)
          end if
        end if
      end do
    end do

#ifdef JEC2DIAG
    if (allocated(t_energies4))   call tr_deallocate(t_energies4,"t_energies4",CAT_UNKNOWN)
    call tr_allocate(t_energies4,1,n_tor_tmp,1,2,1,index_start+nstep, "t_energies4",CAT_UNKNOWN)
    t_energies4 = 0.d0
    call HDF5_array3D_reading(file_id,t_energies4,'energies4')

    if (allocated(energies4))   call tr_deallocate(energies4,"energies4",CAT_UNKNOWN)
    call tr_allocate(energies4,1,n_tor,1,2,1,index_start+nstep, "energies4",CAT_UNKNOWN)
    energies4 = 0.d0
    do m=1,n_tor_tmp,2
      do k=1, n_tor,2 
        if (mode_tmp(m) .eq. mode(k)) then
          if ((m .eq. 1) .and. (k.eq.1)) then
            energies4(k,:,:) = t_energies4(m,:,:)
          else
            energies4(k-1:k,:,:) = t_energies4(m-1:m,:,:)
          end if
        end if
      end do
    end do
#endif

#endif
  end if

  ! Import restart Vacuum 
  call import_HDF5_restart_vacuum(file_id, freeboundary, resistive_wall)
  
  !--- Some parameters need to be scaled when importing a free-boundary equilibrium
  T_0  = T_0  * current_FB_fact / prev_FB_fact
  T_1  = T_1  * current_FB_fact / prev_FB_fact
  FF_0 = FF_0 * current_FB_fact / prev_FB_fact
  FF_1 = FF_1 * current_FB_fact / prev_FB_fact
  prev_FB_fact = current_FB_fact
  
  if (use_pellet) then
     if (index_start .ge. 1) then
        if (allocated(xtime_pellet_R)) call tr_deallocate(xtime_pellet_R,"xtime_pellet_R",CAT_UNKNOWN)
        call tr_allocate(xtime_pellet_R,1,index_start+nstep,"xtime_pellet_R",CAT_UNKNOWN)
        if (allocated(xtime_pellet_Z)) call tr_deallocate(xtime_pellet_Z,"xtime_pellet_Z",CAT_UNKNOWN)
        call tr_allocate(xtime_pellet_Z,1,index_start+nstep,"xtime_pellet_Z",CAT_UNKNOWN)
        if (allocated(xtime_pellet_psi)) call tr_deallocate(xtime_pellet_psi,"xtime_pellet_psi",CAT_UNKNOWN)
        call tr_allocate(xtime_pellet_psi,1,index_start+nstep,"xtime_pellet_psi",CAT_UNKNOWN)
        if (allocated(xtime_pellet_particles)) &
             call tr_deallocate(xtime_pellet_particles,"xtime_pellet_particles",CAT_UNKNOWN)
        call tr_allocate(xtime_pellet_particles,1,index_start+nstep,"xtime_pellet_particles",CAT_UNKNOWN)
        if (allocated(xtime_phys_ablation)) &
             call tr_deallocate(xtime_phys_ablation,"xtime_phys_ablation",CAT_UNKNOWN)
        call tr_allocate(xtime_phys_ablation,1,index_start+nstep,"xtime_phys_ablation",CAT_UNKNOWN)

        call HDF5_array1D_reading(file_id,xtime_pellet_R,"xtime_pellet_R")
        call HDF5_array1D_reading(file_id,xtime_pellet_Z,"xtime_pellet_Z")
        call HDF5_array1D_reading(file_id,xtime_pellet_psi,"xtime_pellet_psi")
        call HDF5_array1D_reading(file_id,xtime_pellet_particles,"xtime_pellet_particles")
        call HDF5_array1D_reading(file_id,xtime_phys_ablation,"xtime_phys_ablation")
     end if
     call HDF5_real_reading(file_id,pellet_R,"pellet_R")
     call HDF5_real_reading(file_id,pellet_Z,"pellet_Z")
     call HDF5_real_reading(file_id,pellet_particles,"pellet_particles")
  endif

#if (defined WITH_Neutrals) || (defined WITH_Impurities)
  if (index_start >= 1) then
    if (allocated(xtime_radiation)) &
      call tr_deallocate(xtime_radiation,"xtime_radiation",CAT_UNKNOWN)
    call tr_allocate(xtime_radiation,1,index_start+nstep,"xtime_radiation",CAT_UNKNOWN)
    call HDF5_array1D_reading(file_id,xtime_radiation,"xtime_radiation")
    if (allocated(xtime_rad_power)) &
      call tr_deallocate(xtime_rad_power,"xtime_rad_power",CAT_UNKNOWN)
    call tr_allocate(xtime_rad_power,1,index_start+nstep,"xtime_rad_power",CAT_UNKNOWN)
    call HDF5_array1D_reading(file_id,xtime_rad_power,"xtime_rad_power")
    if (allocated(xtime_E_ion)) &
      call tr_deallocate(xtime_E_ion,"xtime_E_ion",CAT_UNKNOWN)
    call tr_allocate(xtime_E_ion,1,index_start+nstep,"xtime_E_ion",CAT_UNKNOWN)
    call HDF5_array1D_reading(file_id,xtime_E_ion,"xtime_E_ion")
    if (allocated(xtime_E_ion_power)) &
      call tr_deallocate(xtime_E_ion_power,"xtime_E_ion_power",CAT_UNKNOWN)
    call tr_allocate(xtime_E_ion_power,1,index_start+nstep,"xtime_E_ion_power",CAT_UNKNOWN)
    call HDF5_array1D_reading(file_id,xtime_E_ion_power,"xtime_E_ion_power")
    if (allocated(xtime_P_ei)) &
      call tr_deallocate(xtime_P_ei,"xtime_P_ei",CAT_UNKNOWN)
    call tr_allocate(xtime_P_ei,1,index_start+nstep,"xtime_P_ei",CAT_UNKNOWN)
    call HDF5_array1D_reading(file_id,xtime_P_ei,"xtime_P_ei")
  end if
#endif

  if (using_spi) then
    if (n_spi_tot >= 1) then

      if (index_start >= 1) then
        if (allocated(xtime_spi_ablation)) &
          call tr_deallocate(xtime_spi_ablation,"xtime_spi_ablation",CAT_UNKNOWN)
        call tr_allocate(xtime_spi_ablation,1,n_spi_tot,1,index_start+nstep,"xtime_spi_ablation",CAT_UNKNOWN)
        if (allocated(xtime_spi_ablation_rate)) &
          call tr_deallocate(xtime_spi_ablation_rate,"xtime_spi_ablation_rate",CAT_UNKNOWN)
        call tr_allocate(xtime_spi_ablation_rate,1,n_spi_tot,1,index_start+nstep,"xtime_spi_ablation_rate",CAT_UNKNOWN)
        if (allocated(xtime_spi_ablation_bg)) &
          call tr_deallocate(xtime_spi_ablation_bg,"xtime_spi_ablation_bg",CAT_UNKNOWN)
        call tr_allocate(xtime_spi_ablation_bg,1,n_spi_tot,1,index_start+nstep,"xtime_spi_ablation_bg",CAT_UNKNOWN)
        if (allocated(xtime_spi_ablation_bg_rate)) &
          call tr_deallocate(xtime_spi_ablation_bg_rate,"xtime_spi_ablation_bg_rate",CAT_UNKNOWN)
        call tr_allocate(xtime_spi_ablation_bg_rate,1,n_spi_tot,1,index_start+nstep,"xtime_spi_ablation_bg_rate",CAT_UNKNOWN)

        call HDF5_array2D_reading(file_id,xtime_spi_ablation,"xtime_spi_ablation")
        call HDF5_array2D_reading(file_id,xtime_spi_ablation_rate,"xtime_spi_ablation_rate")

        call H5Lexists_f(file_id,"xtime_spi_ablation_bg",flag_exists,err_exists) !Backward compatibility
        if (flag_exists .and. err_exists == 0) then
          call HDF5_array2D_reading(file_id,xtime_spi_ablation_bg,"xtime_spi_ablation_bg")
          call HDF5_array2D_reading(file_id,xtime_spi_ablation_bg_rate,"xtime_spi_ablation_bg_rate")
        else
          xtime_spi_ablation_bg = 0.
          xtime_spi_ablation_bg_rate = 0.
          write(*,*)"Backward Compatibility: No bg species ablation history information found, assuming none."
        end if
      end if

      call H5Lexists_f(file_id,"n_spi_tot",flag_exists,err_exists) !Backward compatibility
      if (flag_exists .and. err_exists == 0) then
        call HDF5_integer_reading(file_id,n_spi_check,"n_spi_tot")
        if (n_spi_check /= n_spi_tot) then
          write(*,*) "Inconsistency in n_spi_tot detected, exiting!"
          stop
        end if
      else if (n_spi_tot == n_spi(1)) then
        write(*,*)"Backward Compatibility: No n_spi_tot information found, assuming consistent."
      else
        write(*,*)"Backward Compatibility: No n_spi_tot information found, but n_spi_tot is not equal to n_spi(1)."
        stop
      end if

      call H5Lexists_f(file_id,"n_inj",flag_exists,err_exists) !Backward compatibility
      if (flag_exists .and. err_exists == 0) then
        call HDF5_integer_reading(file_id,n_inj_check,"n_inj")
        if (n_inj_check /= n_inj) then
          write(*,*) "Inconsistency in n_inj detected, exiting!"
          stop
        end if
      else if (n_inj == 1) then
        write(*,*)"Backward Compatibility: No n_inj information found, assuming consistent."
      else 
        write(*,*)"Backward Compatibility: No n_inj information found, but n_inj larger than 1, aborting."
        stop
      end if

      allocate (spi_R_arr(n_spi_tot))
      allocate (spi_Z_arr(n_spi_tot))
      allocate (spi_phi_arr(n_spi_tot))
      allocate (spi_phi_init_arr(n_spi_tot))
      allocate (spi_Vel_R_arr(n_spi_tot))
      allocate (spi_Vel_Z_arr(n_spi_tot))
      allocate (spi_Vel_RxZ_arr(n_spi_tot))
      allocate (spi_radius_arr(n_spi_tot))
      allocate (spi_abl_arr(n_spi_tot))
      allocate (spi_species_arr(n_spi_tot))

      call HDF5_array1D_reading(file_id,spi_R_arr,"spi_R_arr")
      call HDF5_array1D_reading(file_id,spi_Z_arr,"spi_Z_arr")
      call HDF5_array1D_reading(file_id,spi_phi_arr,"spi_phi_arr")

      call H5Lexists_f(file_id,"spi_phi_init_arr",flag_exists,err_exists)
      if (flag_exists .and. err_exists == 0) then
        call HDF5_array1D_reading(file_id,spi_phi_init_arr,"spi_phi_init_arr")
      else
        n_spi_begin = 1
        do i_inj = 1, n_inj
          if(n_spi(i_inj)>0) spi_phi_init_arr(n_spi_begin:(n_spi_begin+n_spi(i_inj)-1)) = ns_phi(i_inj)
          n_spi_begin = n_spi_begin + n_spi(i_inj)
        end do
        write(*,*)"Backward Compatibility: No spi_phi_init location found, assuming to be ns_phi."
      end if

      call HDF5_array1D_reading(file_id,spi_Vel_R_arr,"spi_Vel_R_arr")
      call HDF5_array1D_reading(file_id,spi_Vel_Z_arr,"spi_Vel_Z_arr")
      call HDF5_array1D_reading(file_id,spi_Vel_RxZ_arr,"spi_Vel_RxZ_arr")
      call HDF5_array1D_reading(file_id,spi_radius_arr,"spi_radius_arr")
      call HDF5_array1D_reading(file_id,spi_abl_arr,"spi_abl_arr")

      call H5Lexists_f(file_id,"spi_species_arr",flag_exists,err_exists)
      if (flag_exists .and. err_exists == 0) then
        call H5Dopen_f(file_id,"spi_species_arr",dataset,dterr)
        call H5Dget_type_f(dataset,datatype,dterr)
        call H5Tequal_f(datatype,H5T_NATIVE_INTEGER,type_match,dterr)
        call H5Tclose_f(datatype,dterr)
        call H5Dclose_f(dataset,dterr)
        if (type_match .and. dterr == 0) then
          write(*,*) "Backward Compatibility: Converting integer spi_species into double precision"
          allocate (spi_species_arr_old(n_spi_tot))
          call HDF5_array1D_reading_int(file_id,spi_species_arr_old,"spi_species_arr")
          spi_species_arr = REAL(spi_species_arr_old,8)
        else if (dterr == 0) then
          call HDF5_array1D_reading(file_id,spi_species_arr,"spi_species_arr")
        else
          write(*,*) "Error while trying to determine spi_species type, exiting!"
          stop
        end if
      else
#ifdef WITH_Impurities
        spi_species_arr = 1.0
        write(*,*)"Backward Compatibility: No species information found, assuming full impurity."
#endif
#ifdef WITH_Neutrals
        spi_species_arr = 0.0
        write(*,*)"Backward Compatibility: No species information found, assuming pure deuterium."
#endif
      end if

      do i=1, n_spi_tot
        pellets(i)%spi_R       = spi_R_arr(i)
        pellets(i)%spi_Z       = spi_Z_arr(i)
        pellets(i)%spi_phi     = spi_phi_arr(i)
        pellets(i)%spi_phi_init= spi_phi_init_arr(i)
        pellets(i)%spi_Vel_R   = spi_Vel_R_arr(i)
        pellets(i)%spi_Vel_Z   = spi_Vel_Z_arr(i)
        pellets(i)%spi_Vel_RxZ = spi_Vel_RxZ_arr(i)
        pellets(i)%spi_radius  = spi_radius_arr(i)
        pellets(i)%spi_abl     = spi_abl_arr(i)
        pellets(i)%spi_species = spi_species_arr(i)

        write(*,'(A,I5,6ES10.2)') ' *** SHATTERED PELLET PARAMETERS : ',i, pellets(i)%spi_R, pellets(i)%spi_Z, &
                        pellets(i)%spi_phi, pellets(i)%spi_Vel_R, pellets(i)%spi_Vel_Z, pellets(i)%spi_radius
      end do

      deallocate (spi_R_arr)
      deallocate (spi_Z_arr)
      deallocate (spi_phi_arr)
      deallocate (spi_phi_init_arr)
      deallocate (spi_Vel_R_arr)
      deallocate (spi_Vel_Z_arr)
      deallocate (spi_Vel_RxZ_arr)
      deallocate (spi_radius_arr)
      deallocate (spi_abl_arr)
      deallocate (spi_species_arr)
      if (allocated(spi_species_arr_old)) deallocate (spi_species_arr_old)

      if (spi_tor_rot) then
        call HDF5_real_reading(file_id,ns_phi_rotate,"ns_phi_rotate")
      end if


    end if
  end if

  call HDF5_close(file_id)
 
  write(*,*) '************* restart ******************'
  write(*,'(A19,i6,f14.6,A)') ' *  restart time : ',index_start,t_start,' *'
  write(*,*) '****************************************'
  
  do i=2,index_start
    if ( (energies(n_tor,1,i).ne.0.) .and. (energies(n_tor,1,i-1).ne.0.)) then
      Growth_mag  = 0.5d0*log(abs(energies(n_tor,1,i)/energies(n_tor,1,i-1))) &
            / (xtime(i)-xtime(i-1))
    else
      Growth_mag  = 0.
    endif
    if ( (energies(n_tor,2,i).ne.0.) .and. (energies(n_tor,2,i-1).ne.0.)) then
      Growth_kin  = 0.5d0*log(abs(energies(n_tor,2,i)/energies(n_tor,2,i-1))) &
            / (xtime(i)-xtime(i-1))
    else
      Growth_kin  = 0.
    endif

    ! write(*,'(i7,f10.3,200e14.6)') i,xtime(i),energies(1:n_tor,:,i),growth_mag,growth_kin
    ! write(*,'(i7,f10.3,200e14.6)') i,xtime(i),energies(1:n_tor,:,i)
  enddo
 
  ! --- initialise new harmonics (only density and temperature, to be improved)
  n_new_modes = sum(new_mode(1:n_tor))
  if ( (.not. no_pert) .and. (n_new_modes .gt. 0) ) then
    write(*,*), 'Warning:', n_new_modes, ' new modes initialized to noise level' 
    ! --- Using an already computed mode
    if ( (import_perturbation) .and. (n_tor .gt. 1) ) then
      write(*,*) 'ERROR: Importing perturbation from jorek_perturbation.rst file...'
      write(*,*) 'ERROR: Not yet implemeted!'
      stop
    ! --- Using just noise
    else
      amplitude = 1.d-10
      do i=1,node_list%n_nodes
        do m=2,n_tor
          if ( new_mode(m) .eq. 1 ) then
          node_list%node(i)%values(m,:,:) = 0.d0
          node_list%node(i)%values(m,:,var_rho)   = amplitude * node_list%node(i)%values(1,:,var_rho)
#ifdef WITH_TiTe
          node_list%node(i)%values(m,:,var_Ti)   = amplitude * node_list%node(i)%values(1,:,var_Ti)
          node_list%node(i)%values(m,:,var_Te)   = amplitude * node_list%node(i)%values(1,:,var_Te)
#else
          node_list%node(i)%values(m,:,var_T)    = amplitude * node_list%node(i)%values(1,:,var_T)
#endif
#ifdef fullmhd
          node_list%node(i)%values(m,:,var_AR)= amplitude * node_list%node(i)%values(1,:,var_AR)
          node_list%node(i)%values(m,:,var_AZ)= amplitude * node_list%node(i)%values(1,:,var_AZ)
          node_list%node(i)%values(m,:,var_A3)= amplitude * node_list%node(i)%values(1,:,var_A3)
#endif
          end if
        end do
      end do
    endif
  end if

  !call add_pellet(node_list,element_list,25.d0,0.06d0,0.03d0,3.78d0,0.14d0)
  
  ! -> Deallocate temporary arrays 
  call tr_deallocate(mode_tmp,"mode_tmp",CAT_UNKNOWN)
  call tr_deallocate(new_mode,"new_mode",CAT_UNKNOWN)
  
  call tr_deallocate(t_x,"t_x",CAT_UNKNOWN)
  call tr_deallocate(t_values,"t_values",CAT_UNKNOWN)
  call tr_deallocate(t_deltas,"t_deltas",CAT_UNKNOWN)
  call tr_deallocate(t_energies,"t_energies",CAT_UNKNOWN)

#ifdef JECCD                   
  call tr_deallocate(t_energies2,"t_energies2",CAT_UNKNOWN)
  call tr_deallocate(t_energies3,"t_energies3",CAT_UNKNOWN)
#ifdef JEC2DIAG
  call tr_deallocate(t_energies4,"t_energies4",CAT_UNKNOWN)
#endif
#endif
 
#ifdef fullmhd
  call tr_deallocate(t_psi_eq,"t_psi_eq",CAT_UNKNOWN)
  call tr_deallocate(t_Fprof_eq,"t_Fprof",CAT_UNKNOWN) 
#elif altcs
  call tr_deallocate(t_psi_eq,"t_psi_eq",CAT_UNKNOWN)
#endif
 
  call tr_deallocate(t_index,"index",CAT_UNKNOWN)
  call tr_deallocate(t_boundary,"boundary",CAT_UNKNOWN)
  call tr_deallocate(t_axis_node,"axis_node",CAT_UNKNOWN)
  call tr_deallocate(t_parents,"parents",CAT_UNKNOWN)
  call tr_deallocate(t_parent_elem,"parent_elem",CAT_UNKNOWN)
  call tr_deallocate(t_ref_lambda,"ref_lambda",CAT_UNKNOWN)
  call tr_deallocate(t_ref_mu,"ref_mu",CAT_UNKNOWN)
  call tr_deallocate(t_constrained,"constrained",CAT_UNKNOWN)
 
  call tr_deallocate(t_vertex,"t_vertex",CAT_UNKNOWN)
  call tr_deallocate(t_neighbours,"t_neighbours",CAT_UNKNOWN)
  call tr_deallocate(t_size,"t_size",CAT_UNKNOWN)
  call tr_deallocate(t_father,"t_father",CAT_UNKNOWN)
  call tr_deallocate(t_n_sons,"t_n_sons",CAT_UNKNOWN)
  call tr_deallocate(t_n_gen,"t_n_gen",CAT_UNKNOWN)
  call tr_deallocate(t_sons,"t_sons",CAT_UNKNOWN)
  call tr_deallocate(t_contain_node,"t_contain_node",CAT_UNKNOWN)
  call tr_deallocate(t_nref,"t_nref",CAT_UNKNOWN)

#else
  write (6,*) " ERROR: trying to import with hdf5 but USE_HDF5 was not set at compile-time"
#endif
  call populate_element_rtree(node_list, element_list)

  return
end subroutine import_hdf5_restart
end module mod_import_restart

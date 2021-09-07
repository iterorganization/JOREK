module mod_export_restart
implicit none
contains
!> Export the current simulation state as a restart file that can be read back into JOREK or into
!! a diagnostic program by the routine import_restart.
subroutine export_restart(node_list,element_list,filename)

  use mod_parameters
  use data_structure
  use phys_module
  use pellet_module

  implicit none

  ! --- Routine parameters
  type(type_node_list),    intent(in) :: node_list
  type(type_element_list), intent(in) :: element_list
  character(len=*)       , intent(in) :: filename

  character*17 :: fileout

  if ( rst_hdf5 == 0 ) then
    ! --- Write restart binary file
    fileout = trim(filename)//".rst"
    write (6,*) " =============>, jorek2, filename = ", fileout
    call export_binary_restart(node_list, element_list, fileout)
  elseif ( rst_hdf5 == 1 ) then
    ! --- Write restart HDF5 file
    fileout = trim(filename)//".h5"
    write (6,*) " =============>, jorek2, filename = ", fileout
    call export_hdf5_restart(node_list, element_list, fileout)
  end if

end subroutine export_restart

!
! Export in a binary restart file
subroutine export_binary_restart(node_list,element_list,filename)

  use mod_parameters
  use data_structure
  use phys_module
  use pellet_module
  use vacuum, only: export_restart_vacuum

  implicit none

#include "version.h"

  ! --- Routine parameters
  type(type_node_list),    intent(in) :: node_list
  type(type_element_list), intent(in) :: element_list
  character(len=*),        intent(in) :: filename

  ! --- Local variables
  integer :: i
  character*50 :: version_control
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

  ! -> Write binary restart file
  open(21, file=filename, form='unformatted', status='replace', action='write')

  write(21) n_tor
  write(21) node_list%n_nodes,element_list%n_elements
  write(21) node_list%n_dof

  do i=1,node_list%n_nodes
     write(21) node_list%node(i)%x
     write(21) node_list%node(i)%values
     write(21) node_list%node(i)%deltas
#ifdef fullmhd
     write(21) node_list%node(i)%psi_eq               !< equilibrium flux at the nodes
     write(21) node_list%node(i)%Fprof_eq             !< equilibrium profile R*B_phi at the nodes
#elif altcs
     write(21) node_list%node(i)%psi_eq               !< equilibrium flux at the nodes
#endif
     write(21) node_list%node(i)%index
     write(21) node_list%node(i)%boundary
     write(21) node_list%node(i)%axis_node
     write(21) node_list%node(i)%parents
     write(21) node_list%node(i)%parent_elem
     write(21) node_list%node(i)%ref_lambda
     write(21) node_list%node(i)%ref_mu
     write(21) node_list%node(i)%constrained
  enddo

  write(21) element_list%element(1:element_list%n_elements)
  write(21) tstep,eta,visco,visco_par
  write(21) index_now
  write(21) t_now

  if (index_now .gt. 0) then
     write(21) xtime(1:index_now)
     write(21) energies(:,:,1:index_now)
#ifdef JECCD
     write(21) energies2(:,:,1:index_now)
     write(21) energies3(:,:,1:index_now)
#ifdef JEC2DIAG
     write(21) energies4(:,:,1:index_now)
#endif
#endif
  endif

  call export_restart_vacuum(21, freeboundary, resistive_wall)

  if (use_pellet) then
     if (index_now .gt. 0) then
        write(21) xtime_pellet_R(1:index_now)
        write(21) xtime_pellet_Z(1:index_now)
        write(21) xtime_pellet_psi(1:index_now)
        write(21) xtime_pellet_particles(1:index_now)
        write(21) xtime_phys_ablation(1:index_now)
     endif
     write(21) pellet_particles, pellet_R, pellet_Z
  endif

  ! Radiation and ionization energy history
  if (index_now .gt. 0) write(21) xtime_radiation(1:index_now)
  if (index_now .gt. 0) write(21) xtime_rad_power(1:index_now)
  if (index_now .gt. 0) write(21) xtime_E_ion(1:index_now)
  if (index_now .gt. 0) write(21) xtime_E_ion_power(1:index_now)
  if (index_now .gt. 0) write(21) xtime_P_ei(1:index_now)

  ! Dynamically allocate memeries for temporary arrays in order to export
  if (using_spi .and. n_spi_tot >= 1) then

    if (index_now .gt. 0) then
      write(21) xtime_spi_ablation(:,1:index_now)
      write(21) xtime_spi_ablation_rate(:,1:index_now)
      write(21) xtime_spi_ablation_bg(:,1:index_now)
      write(21) xtime_spi_ablation_bg_rate(:,1:index_now)
    endif

    write(21) n_inj
    write(21) n_spi_tot

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

    do i=1, n_spi_tot
      spi_R_arr(i)       = pellets(i)%spi_R
      spi_Z_arr(i)       = pellets(i)%spi_Z
      spi_phi_arr(i)     = pellets(i)%spi_phi
      spi_phi_init_arr(i)= pellets(i)%spi_phi_init
      spi_Vel_R_arr(i)   = pellets(i)%spi_Vel_R
      spi_Vel_Z_arr(i)   = pellets(i)%spi_Vel_Z
      spi_Vel_RxZ_arr(i) = pellets(i)%spi_Vel_RxZ
      spi_radius_arr(i)  = pellets(i)%spi_radius
      spi_abl_arr(i)     = pellets(i)%spi_abl
      spi_species_arr(i) = pellets(i)%spi_species
    end do

    write(21) spi_R_arr(1:n_spi_tot)
    write(21) spi_Z_arr(1:n_spi_tot)
    write(21) spi_phi_arr(1:n_spi_tot)
    write(21) spi_phi_init_arr(1:n_spi_tot)
    write(21) spi_Vel_R_arr(1:n_spi_tot)
    write(21) spi_Vel_Z_arr(1:n_spi_tot)
    write(21) spi_Vel_RxZ_arr(1:n_spi_tot)
    write(21) spi_radius_arr(1:n_spi_tot)
    write(21) spi_abl_arr(1:n_spi_tot)
    write(21) spi_species_arr(1:n_spi_tot)

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
      write(21) ns_phi_rotate
    end if

  end if
   
  ! save Revision control
  write(version_control,'(A)') trim(adjustl(RCS_VERSION))
  write(21) version_control

  ! save parameters
  write(21) jorek_model
  
  write(21) n_var
  write(21) n_order
  write(21) n_tor
  write(21) n_period
  write(21) n_plane
  write(21) n_vertex_max
  write(21) n_nodes_max
  write(21) n_elements_max
  write(21) n_boundary_max
  write(21) n_pieces_max
  write(21) n_degrees
  write(21) nref_max
  write(21) n_ref_list

  close(21)

  return
end subroutine export_binary_restart

 ! 
 ! Export in a HDF5 binary restart file
subroutine export_hdf5_restart(node_list,element_list,filename)
 
  use data_structure
  use phys_module
  use pellet_module
  use vacuum, only : export_HDF5_restart_vacuum
  
#ifdef USE_HDF5
  use hdf5
  use hdf5_io_module
  use tr_module
  use mod_parameters
#endif
 
  implicit none
 
#include "version.h"
 
  ! --- Routine parameters
  type(type_node_list),    intent(in) :: node_list
  type(type_element_list), intent(in) :: element_list
  character*(*),           intent(in) :: filename

  ! --- Local variables
  integer :: i
  character(len=50)        :: version_control

#ifdef USE_HDF5
  integer(HID_T)     :: file_id
  integer            :: ind, ierr

  ! type_node, node_list%n_nodes
  real(RKIND), allocatable :: t_x(:,:,:,:)                 ! n_coord_tor, n_order+1, n_dim
  real(RKIND), allocatable :: t_values(:,:,:,:)            !       n_tor, n_order+1, n_var
  real(RKIND), allocatable :: t_deltas(:,:,:,:)            !       n_tor, n_order+1, n_var

  real(RKIND), allocatable :: t_psi_eq(:,:)                ! n_order+1
  real(RKIND), allocatable :: t_Fprof_eq(:,:)              ! n_order+1

  integer,     allocatable :: t_index(:,:)                 ! n_order+1
  integer,     allocatable :: t_boundary(:)                ! 
  character,   allocatable :: t_axis_node(:)     
  integer,     allocatable :: t_parents(:,:)               ! 2
  integer,     allocatable :: t_parent_elem(:)             ! 
  real(RKIND), allocatable :: t_ref_lambda(:)
  real(RKIND), allocatable :: t_ref_mu(:)
  character,   allocatable :: t_constrained(:)     

  ! element, element_list%n_elements
  integer,     allocatable :: t_vertex(:,:)                ! n_vertex_max
  integer,     allocatable :: t_neighbours(:,:)            ! n_vertex_max
  real(RKIND), allocatable :: t_size(:,:,:)                ! n_vertex_max,n_order+1
  integer,     allocatable :: t_father(:)
  integer,     allocatable :: t_n_sons(:)
  integer,     allocatable :: t_n_gen(:)
  integer,     allocatable :: t_sons(:,:)                  ! 4
  integer,     allocatable :: t_contain_node(:,:)          ! 5
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

  ! index_now+nstep
  real(RKIND), allocatable :: t_xtime(:)                   ! nstep
  real(RKIND), allocatable :: t_energies(:,:,:)            ! n_tor,2,index_start+nstep
#ifdef JECCD                                          
  real(RKIND), allocatable :: t_energies2(:,:,:)           ! n_tor,2,index_start+nstep
  real(RKIND), allocatable :: t_energies3(:,:,:)           ! n_tor,2,index_start+nstep
#ifdef JEC2DIAG                                       
  real(RKIND), allocatable :: t_energies4(:,:,:)           ! n_tor,2,index_start+nstep
#endif
#endif

  ! type_node, node_list%n_nodes
  call tr_allocate(t_x,1,node_list%n_nodes,1,n_coord_tor,1,n_order+1,1,n_dim, &
      "node_list%x",CAT_UNKNOWN)
  call tr_allocate(t_values,1,node_list%n_nodes,1,n_tor,1,n_order+1,1,n_var, &
       "node_list%values",CAT_UNKNOWN)
  call tr_allocate(t_deltas,1,node_list%n_nodes,1,n_tor,1,n_order+1,1,n_var, &
       "node_list%deltas",CAT_UNKNOWN)

#ifdef fullmhd
  call tr_allocate(t_psi_eq,1,node_list%n_nodes,1,n_order+1, &
       "node_list%psi_eq",CAT_UNKNOWN)
  call tr_allocate(t_Fprof_eq,1,node_list%n_nodes,1,n_order+1, &
       "node_list%Fprof_eq",CAT_UNKNOWN)
#elif altcs
  call tr_allocate(t_psi_eq,1,node_list%n_nodes,1,n_order+1, &
       "node_list%psi_eq",CAT_UNKNOWN)
#endif

  call tr_allocate(t_index,1,node_list%n_nodes,1,n_order+1,"index",CAT_UNKNOWN)
  call tr_allocate(t_boundary,1,node_list%n_nodes,"boundary",CAT_UNKNOWN)
  call tr_allocate(t_axis_node,1,node_list%n_nodes,"axis_node",CAT_UNKNOWN)
  call tr_allocate(t_parents,1,node_list%n_nodes,1,2,"parent",CAT_UNKNOWN)
  call tr_allocate(t_parent_elem,1,node_list%n_nodes,"parent_elem",CAT_UNKNOWN)
  call tr_allocate(t_ref_lambda,1,node_list%n_nodes,"ref_lambade",CAT_UNKNOWN)
  call tr_allocate(t_ref_mu,1,node_list%n_nodes,"ref_mu",CAT_UNKNOWN)
  call tr_allocate(t_constrained,1,node_list%n_nodes,"constrained",CAT_UNKNOWN)

  ! element_list%n_elements
  call tr_allocate(t_vertex,1,element_list%n_elements,1,n_vertex_max,"vertex",CAT_UNKNOWN)
  call tr_allocate(t_neighbours,1,element_list%n_elements,1,n_vertex_max,"neighbours",CAT_UNKNOWN)
  call tr_allocate(t_size,1,element_list%n_elements,1,n_vertex_max,1,n_order+1,"size",CAT_UNKNOWN)
  call tr_allocate(t_father,1,element_list%n_elements,"father",CAT_UNKNOWN)
  call tr_allocate(t_n_sons,1,element_list%n_elements,"n_sons",CAT_UNKNOWN)
  call tr_allocate(t_n_gen,1,element_list%n_elements,"n_gen",CAT_UNKNOWN)
  call tr_allocate(t_sons,1,element_list%n_elements,1,4,"sons",CAT_UNKNOWN)
  call tr_allocate(t_contain_node,1,element_list%n_elements,1,5,"contain_node",CAT_UNKNOWN)
  call tr_allocate(t_nref,1,element_list%n_elements,"nref",CAT_UNKNOWN)

  ! index_now+nstep
  if (index_now .gt. 0) then
     if (allocated(t_xtime)) call tr_deallocate(t_xtime,"xtime",CAT_UNKNOWN)
     call tr_allocate(t_xtime,1,index_now,"xtime",CAT_UNKNOWN)
     t_xtime(:) = xtime(1:index_now)

     if (allocated(t_energies)) call tr_deallocate(t_energies,"energies",CAT_UNKNOWN)
     call tr_allocate(t_energies,1,n_tor,1,2,1,index_now,"energies",CAT_UNKNOWN)
     t_energies(:,:,:) = energies(:,:,1:index_now)
#ifdef JECCD
     if (allocated(t_energies2)) call tr_deallocate(t_energies2,"energies2",CAT_UNKNOWN)
     call tr_allocate(t_energies2,1,n_tor,1,2,1,index_now,"energies2",CAT_UNKNOWN)
     t_energies2(:,:,:) = t_energies2(:,:,1:index_now)

     if (allocated(t_energies3)) call tr_deallocate(t_energies3,"energies3",CAT_UNKNOWN)
     call tr_allocate(t_energies3,1,n_tor,1,2,1,index_now, "energies3",CAT_UNKNOWN)
     t_energies3(:,:,:) = t_energies3(:,:,1:index_now)

#ifdef JEC2DIAG
     if (allocated(t_energies4)) call tr_deallocate(t_energies4,"energies4",CAT_UNKNOWN)
     call tr_allocate(t_energies4,1,n_tor,1,2,1,index_now, "energies4",CAT_UNKNOWN)
     t_energies4(:,:,:) = t_energies4(:,:,1:index_now)
#endif
#endif

  end if

  !
  do i=1,node_list%n_nodes
     t_x(i,:,:,:)        = node_list%node(i)%x
     t_values(i,:,:,:) = node_list%node(i)%values
     t_deltas(i,:,:,:) = node_list%node(i)%deltas

#ifdef fullmhd
     t_psi_eq(i,:)     = node_list%node(i)%psi_eq
     t_Fprof_eq(i,:)   = node_list%node(i)%Fprof_eq
#elif altcs
     t_psi_eq(i,:)     = node_list%node(i)%psi_eq
#endif

     t_index(i,:)      = node_list%node(i)%index
     t_boundary(i)     = node_list%node(i)%boundary
     if (node_list%node(i)%axis_node) then
        t_axis_node(i)  = 'T'
     else
        t_axis_node(i)  = 'F'
     end if
     t_parents(i,:)    = node_list%node(i)%parents(1:2)
     t_parent_elem(i)  = node_list%node(i)%parent_elem
     t_ref_lambda(i)   = node_list%node(i)%ref_lambda
     t_ref_mu(i)       = node_list%node(i)%ref_mu
     if (node_list%node(i)%constrained) then
        t_constrained(i)  = 'T'
     else
        t_constrained(i)  = 'F'
     end if
  end do

  do i=1,element_list%n_elements
     t_vertex(i,:)       = element_list%element(i)%vertex
     t_neighbours(i,:)   = element_list%element(i)%neighbours
     t_size(i,:,:)       = element_list%element(i)%size
     t_father(i)         = element_list%element(i)%father
     t_n_sons(i)         = element_list%element(i)%n_sons
     t_n_gen(i)          = element_list%element(i)%n_gen
     t_sons(i,:)         = element_list%element(i)%sons
     t_contain_node(i,:) = element_list%element(i)%contain_node
     t_nref(i)           = element_list%element(i)%nref
  end do

  ! -> Create and open HDF5 file
  write (6,*) " HDF5 file ", filename
  call HDF5_create(trim(filename),file_id,ierr)
  if (ierr.ne.0) then
     print*,'pglobal_id = ',pglobal_id, &
          ' ==> error for opening of HDF5 file',filename
  end if
  
  ! Store the hdf5 restart file version
  write(*,*) 'Exporting HDF5 restart file with rst_hdf5_version=', rst_hdf5_version
  if ( rst_hdf5_version > rst_hdf5_version_supported ) then
    write(*,*) 'ERROR: Cannot write HDF5 restart file with rst_hdf5_version=', rst_hdf5_version
    write(*,*) '  This code version supports rst_hdf5_version_supported=', rst_hdf5_version_supported
    write(*,*) '  Please select an appropriate value in the namelist or use the default value.'
    stop
  end if
  call HDF5_integer_saving(file_id,rst_hdf5_version,'rst_hdf5_version'//char(0)) 

  ! -> Save version of revision control system
  write(version_control,'(A)') trim(adjustl(RCS_VERSION))
  version_control = trim(adjustl(version_control))
  call HDF5_char_saving(file_id,version_control,"RCS_version"//char(0))

  ! -> Save parameters
  call HDF5_integer_saving(file_id,jorek_model,'jorek_model'//char(0))
  call HDF5_integer_saving(file_id,n_var,'n_var'//char(0))
  call HDF5_integer_saving(file_id,n_dim,'n_dim'//char(0))
  call HDF5_integer_saving(file_id,n_order,'n_order'//char(0))
  call HDF5_integer_saving(file_id,n_tor,'n_tor'//char(0))
  call HDF5_integer_saving(file_id,n_coord_tor,'n_coord_tor'//char(0))
  call HDF5_integer_saving(file_id,n_period,'n_period'//char(0))
  call HDF5_integer_saving(file_id,n_plane,'n_plane'//char(0))
  call HDF5_integer_saving(file_id,n_vertex_max,'n_vertex_max'//char(0))
  call HDF5_integer_saving(file_id,n_nodes_max,'n_nodes_max'//char(0))
  call HDF5_integer_saving(file_id,n_elements_max,'n_elements_max'//char(0))
  call HDF5_integer_saving(file_id,n_boundary_max,'n_boundary_max'//char(0))
  call HDF5_integer_saving(file_id,n_pieces_max,'n_pieces_max'//char(0))
  call HDF5_integer_saving(file_id,n_degrees,'n_degrees'//char(0))
  call HDF5_integer_saving(file_id,nref_max,'nref_max'//char(0))
  call HDF5_integer_saving(file_id,n_ref_list,'n_ref_list'//char(0))

  ! -> 
  call HDF5_integer_saving(file_id,node_list%n_nodes,'n_nodes'//char(0))
  call HDF5_integer_saving(file_id,element_list%n_elements,'n_elements'//char(0))
  call HDF5_integer_saving(file_id,node_list%n_dof,'n_dof'//char(0))
  
  if (rst_hdf5_version .eq. 2) then
    call HDF5_array4D_saving(file_id,t_x, &
         node_list%n_nodes,n_coord_tor,n_order+1,n_dim,'x'//char(0))
  else
    call HDF5_array3D_saving(file_id,t_x(:,1,:,:), &
         node_list%n_nodes,n_order+1,n_dim,'x'//char(0))
  endif
  call HDF5_array4D_saving(file_id,t_values, &
       node_list%n_nodes,n_tor,n_order+1,n_var,'values'//char(0))
  call HDF5_array4D_saving(file_id,t_deltas, &
       node_list%n_nodes,n_tor,n_order+1,n_var,'deltas'//char(0))

#ifdef fullmhd
  call HDF5_array2D_saving(file_id,t_psi_eq, &
       node_list%n_nodes,n_order+1,'psi_eq'//char(0))
  call HDF5_array2D_saving(file_id,t_Fprof_eq, &
       node_list%n_nodes,n_order+1,'Fprof_eq'//char(0))
#elif altcs
  call HDF5_array2D_saving(file_id,t_psi_eq, &
       node_list%n_nodes,n_order+1,'psi_eq'//char(0))
#endif

  call HDF5_array2D_saving_int(file_id,t_index, &
       node_list%n_nodes,n_order+1,'index'//char(0))
  call HDF5_array1D_saving_int(file_id,t_boundary, &
       node_list%n_nodes,'boundary'//char(0))
  call HDF5_array1D_saving_char(file_id,t_axis_node, &
       node_list%n_nodes,'axis_node'//char(0))
  call HDF5_array2D_saving_int(file_id,t_parents, &
       node_list%n_nodes,2,'parents'//char(0))
  call HDF5_array1D_saving_int(file_id,t_parent_elem, &
       node_list%n_nodes,'parent_elem'//char(0))
  call HDF5_array1D_saving(file_id,t_ref_lambda, &
       node_list%n_nodes,'ref_lambda'//char(0))
  call HDF5_array1D_saving(file_id,t_ref_mu, &
       node_list%n_nodes,'ref_mu'//char(0))
  call HDF5_array1D_saving_char(file_id,t_constrained, &
       node_list%n_nodes,'constrained'//char(0))

  call HDF5_array2D_saving_int(file_id,t_vertex, &
       element_list%n_elements,n_vertex_max,'vertex'//char(0))
  call HDF5_array2D_saving_int(file_id,t_neighbours, &
       element_list%n_elements,n_vertex_max,'neighbours'//char(0))
  call HDF5_array3D_saving(file_id,t_size, &
       element_list%n_elements,n_vertex_max,n_order+1,'size'//char(0))
  call HDF5_array1D_saving_int(file_id,t_father, &
       element_list%n_elements,'father'//char(0))
  call HDF5_array1D_saving_int(file_id,t_n_sons, &
       element_list%n_elements,'n_sons'//char(0))
  call HDF5_array1D_saving_int(file_id,t_n_gen, &
       element_list%n_elements,'n_gen'//char(0))
  call HDF5_array2D_saving_int(file_id,t_sons, &
       element_list%n_elements,4,'sons'//char(0))
  call HDF5_array2D_saving_int(file_id,t_contain_node, &
       element_list%n_elements,5,'contain_node'//char(0))
  call HDF5_array1D_saving_int(file_id,t_nref, &
       element_list%n_elements,'nref'//char(0))
  call HDF5_real_saving(file_id,tstep,'tstep'//char(0))
  call HDF5_real_saving(file_id,eta,'eta'//char(0))
  call HDF5_real_saving(file_id,visco,'visco'//char(0))
  call HDF5_real_saving(file_id,visco_par,'visco_par'//char(0))
  call HDF5_integer_saving(file_id,index_now,'index_now'//char(0)) 
  call HDF5_real_saving(file_id,t_now,'t_now'//char(0))
  call HDF5_real_saving(file_id,central_density,'central_density'//char(0))
  call HDF5_real_saving(file_id,central_mass,'central_mass'//char(0))
  call HDF5_real_saving(file_id,F0,'F0'//char(0))
  call HDF5_real_saving(file_id,sqrt_mu0_rho0,'sqrt_mu0_rho0'//char(0))
  call HDF5_real_saving(file_id,sqrt_mu0_rho0,'t_norm'//char(0))
  call HDF5_real_saving(file_id,sqrt_mu0_over_rho0,'sqrt_mu0_over_rho0'//char(0))

  if (index_now .gt. 0) then
     call HDF5_array1D_saving(file_id,t_xtime,index_now,'xtime'//char(0))
     call HDF5_array3D_saving(file_id,t_energies, &
          n_tor,2,index_now,'energies'//char(0))
     !           n_tor,2,index_now,'energies'//char(0))

     call HDF5_array1D_saving(file_id,R_axis_t(1:index_now),index_now,'R_axis_t'//char(0))
     call HDF5_array1D_saving(file_id,Z_axis_t(1:index_now),index_now,'Z_axis_t'//char(0))
     call HDF5_array1D_saving(file_id,psi_axis_t(1:index_now),index_now,'psi_axis_t'//char(0))
     call HDF5_array2D_saving(file_id,R_xpoint_t(1:index_now,1:2),index_now,2,'R_xpoint_t'//char(0))
     call HDF5_array2D_saving(file_id,Z_xpoint_t(1:index_now,1:2),index_now,2,'Z_xpoint_t'//char(0))
     call HDF5_array2D_saving(file_id,psi_xpoint_t(1:index_now,1:2),index_now,2,'psi_xpoint_t'//char(0))
     call HDF5_array1D_saving(file_id,R_bnd_t(1:index_now),index_now,'R_bnd_t'//char(0))
     call HDF5_array1D_saving(file_id,Z_bnd_t(1:index_now),index_now,'Z_bnd_t'//char(0))
     call HDF5_array1D_saving(file_id,psi_bnd_t(1:index_now),index_now,'psi_bnd_t'//char(0))
     call HDF5_array1D_saving(file_id,current_t(1:index_now),index_now,'current_t'//char(0))
     call HDF5_array1D_saving(file_id,beta_p_t(1:index_now),index_now,'beta_p_t'//char(0))
     call HDF5_array1D_saving(file_id,beta_t_t(1:index_now),index_now,'beta_t_t'//char(0))
     call HDF5_array1D_saving(file_id,beta_n_t(1:index_now),index_now,'beta_n_t'//char(0))
     call HDF5_array1D_saving(file_id,density_tot_t(1:index_now),index_now,'density_tot_t'//char(0))
     call HDF5_array1D_saving(file_id,density_in_t(1:index_now),index_now,'density_in_t'//char(0))
     call HDF5_array1D_saving(file_id,density_out_t(1:index_now),index_now,'density_out_t'//char(0))
     call HDF5_array1D_saving(file_id,pressure_in_t(1:index_now),index_now,'pressure_in_t'//char(0))
     call HDF5_array1D_saving(file_id,pressure_out_t(1:index_now),index_now,'pressure_out_t'//char(0))
     call HDF5_array1D_saving(file_id,heat_src_in_t(1:index_now),index_now,'heat_src_in_t'//char(0))
     call HDF5_array1D_saving(file_id,heat_src_out_t(1:index_now),index_now,'heat_src_out_t'//char(0))
     call HDF5_array1D_saving(file_id,heat_src_tot_t(1:index_now),index_now,'heat_src_tot_t'//char(0))
     call HDF5_array1D_saving(file_id,part_src_in_t(1:index_now),index_now,'part_src_in_t'//char(0))
     call HDF5_array1D_saving(file_id,part_src_out_t(1:index_now),index_now,'part_src_out_t'//char(0))
     call HDF5_array1D_saving(file_id,part_src_tot_t(1:index_now),index_now,'part_src_tot_t'//char(0))

     call HDF5_array1D_saving(file_id,E_tot_t(1:index_now),index_now,'E_tot_t'//char(0))
     call HDF5_array1D_saving(file_id,helicity_tot_t(1:index_now),index_now,'helicity_tot_t'//char(0))
     call HDF5_array1D_saving(file_id,Ip_tot_t(1:index_now),index_now,'Ip_tot_t'//char(0))
     call HDF5_array1D_saving(file_id,thermal_tot_t(1:index_now),index_now,'thermal_tot_t'//char(0))
     call HDF5_array1D_saving(file_id,thermal_e_tot_t(1:index_now),index_now,'thermal_e_tot_t'//char(0))
     call HDF5_array1D_saving(file_id,thermal_i_tot_t(1:index_now),index_now,'thermal_i_tot_t'//char(0))
     call HDF5_array1D_saving(file_id,kin_par_tot_t(1:index_now),index_now,'kin_par_tot_t'//char(0))
     call HDF5_array1D_saving(file_id,kin_perp_tot_t(1:index_now),index_now,'kin_perp_tot_t'//char(0))
     call HDF5_array1D_saving(file_id,Wmag_tot_t(1:index_now),index_now,'Wmag_tot_t'//char(0))
     call HDF5_array1D_saving(file_id,ohmic_tot_t(1:index_now),index_now,'ohmic_tot_t'//char(0))
     call HDF5_array1D_saving(file_id,Magwork_tot_t(1:index_now),index_now,'Magwork_tot_t'//char(0))
     call HDF5_array1D_saving(file_id,flux_qpar_t(1:index_now),index_now,'flux_qpar_t'//char(0))
     call HDF5_array1D_saving(file_id,flux_qperp_t(1:index_now),index_now,'flux_qperp_t'//char(0))
     call HDF5_array1D_saving(file_id,flux_kinpar_t(1:index_now),index_now,'flux_kinpar_t'//char(0))
     call HDF5_array1D_saving(file_id,flux_poynting_t(1:index_now),index_now,'flux_poynting_t'//char(0))
     call HDF5_array1D_saving(file_id,flux_pvn_t(1:index_now),index_now,'flux_Pvn_t'//char(0))
     call HDF5_array1D_saving(file_id,part_flux_Dpar_t(1:index_now),index_now,'part_flux_Dpar_t'//char(0))
     call HDF5_array1D_saving(file_id,part_flux_Dperp_t(1:index_now),index_now,'part_flux_Dperp_t'//char(0))
     call HDF5_array1D_saving(file_id,part_flux_Vpar_t(1:index_now),index_now,'part_flux_Vpar_t'//char(0))
     call HDF5_array1D_saving(file_id,part_flux_Vperp_t(1:index_now),index_now,'part_flux_Vperp_t'//char(0))
     call HDF5_array1D_saving(file_id,npart_flux_t(1:index_now),index_now,'npart_flux_t'//char(0))
     call HDF5_array1D_saving(file_id,dpart_tot_dt(1:index_now),index_now,'dpart_tot_dt'//char(0))
     call HDF5_array1D_saving(file_id,dnpart_tot_dt(1:index_now),index_now,'dnpart_tot_dt'//char(0))
     call HDF5_array1D_saving(file_id,npart_tot_t(1:index_now),index_now,'npart_tot_t'//char(0))
     call HDF5_array1D_saving(file_id,dE_tot_dt(1:index_now),index_now,'dE_tot_dt'//char(0))
     call HDF5_array1D_saving(file_id,dWmag_tot_dt(1:index_now),index_now,'dWmag_tot_dt'//char(0))
     call HDF5_array1D_saving(file_id,dthermal_tot_dt(1:index_now),index_now,'dthermal_tot_dt'//char(0))
     call HDF5_array1D_saving(file_id,dkinpar_tot_dt(1:index_now),index_now,'dkinpar_tot_dt'//char(0))
     call HDF5_array1D_saving(file_id,dkinperp_tot_dt(1:index_now),index_now,'dkinperp_tot_dt'//char(0))
     call HDF5_array1D_saving(file_id,thmwork_tot_t(1:index_now),index_now,'thmwork_tot_t'//char(0))
     call HDF5_array1D_saving(file_id,viscopar_dissip_tot_t(1:index_now),index_now,'viscopar_dissip_tot_t'//char(0))
     call HDF5_array1D_saving(file_id,friction_dissip_tot_t(1:index_now),index_now,'friction_dissip_tot_t'//char(0))
     call HDF5_array1D_saving(file_id,viscopar_flux_t(1:index_now),index_now,'viscopar_flux_t'//char(0))
     call HDF5_array1D_saving(file_id,li3_t(1:index_now),index_now,'li3_t'//char(0))
     call HDF5_array1D_saving(file_id,li3_tot_t(1:index_now),index_now,'li3_tot_t'//char(0))
     call HDF5_array1D_saving(file_id,area_t(1:index_now),index_now,'area_t'//char(0))
     call HDF5_array1D_saving(file_id,volume_t(1:index_now),index_now,'volume_t'//char(0))
     call HDF5_array1D_saving(file_id,mag_ener_src_tot(1:index_now),index_now,'mag_ener_src_tot'//char(0))

#ifdef JECCD                   
     call HDF5_array3D_saving(file_id,t_energies2, &
          n_tor,2,index_now,'energies2'//char(0))
     !           n_tor,2,index_now,'energies2'//char(0))
     call HDF5_array3D_saving(file_id,t_energies3, &
          n_tor,2,index_now,'energies3'//char(0))
     !           n_tor,2,index_now,'energies3'//char(0))
#ifdef JEC2DIAG
     call HDF5_array3D_saving(file_id,t_energies4, &
          n_tor,2,index_now,'energies4'//char(0))
     !           n_tor,2,index_now,'energies4'//char(0))
#endif
#endif
  end if

  if (use_pellet) then
     if (index_now .gt. 0) then
        call HDF5_array1D_saving(file_id,xtime_pellet_R, &
             index_now,'xtime_pellet_R'//char(0))
        call HDF5_array1D_saving(file_id,xtime_pellet_Z, &
             index_now,'xtime_pellet_Z'//char(0))
        call HDF5_array1D_saving(file_id,xtime_pellet_psi, &
             index_now,'xtime_pellet_psi'//char(0))
        call HDF5_array1D_saving(file_id,xtime_pellet_particles, &
             index_now,'xtime_pellet_particles'//char(0))
        call HDF5_array1D_saving(file_id,xtime_phys_ablation, &
             index_now,'xtime_phys_ablation'//char(0))
     end if
     call HDF5_real_saving(file_id,pellet_particles,"pellet_particles"//char(0))
     call HDF5_real_saving(file_id,pellet_R,"pellet_R"//char(0))
     call HDF5_real_saving(file_id,pellet_Z,"pellet_Z"//char(0))
  end if

  ! Radiation and ionization energy history
  if (index_now .gt. 0) then
    if ( allocated(xtime_radiation)   ) call HDF5_array1D_saving(file_id,xtime_radiation, index_now,'xtime_radiation'//char(0))
    if ( allocated(xtime_rad_power)   ) call HDF5_array1D_saving(file_id,xtime_rad_power, index_now,'xtime_rad_power'//char(0))
    if ( allocated(xtime_E_ion)       ) call HDF5_array1D_saving(file_id,xtime_E_ion, index_now,'xtime_E_ion'//char(0))
    if ( allocated(xtime_E_ion_power) ) call HDF5_array1D_saving(file_id,xtime_E_ion_power, index_now,'xtime_E_ion_power'//char(0))
    if ( allocated(xtime_P_ei)        ) call HDF5_array1D_saving(file_id,xtime_P_ei, index_now,'xtime_P_ei'//char(0))
  end if

  ! Dynamically allocate memeries for temporary arrays in order to export
  if (using_spi .and. n_spi_tot>=1) then
    if (index_now .gt. 0) then
      call HDF5_array2D_saving(file_id,xtime_spi_ablation, &
             n_spi_tot,index_now,'xtime_spi_ablation'//char(0))
      call HDF5_array2D_saving(file_id,xtime_spi_ablation_rate, &
             n_spi_tot,index_now,'xtime_spi_ablation_rate'//char(0))
      call HDF5_array2D_saving(file_id,xtime_spi_ablation_bg, &
             n_spi_tot,index_now,'xtime_spi_ablation_bg'//char(0))
      call HDF5_array2D_saving(file_id,xtime_spi_ablation_bg_rate, &
             n_spi_tot,index_now,'xtime_spi_ablation_bg_rate'//char(0))
    end if

    call HDF5_integer_saving(file_id,n_inj,"n_inj"//char(0))
    call HDF5_integer_saving(file_id,n_spi_tot,"n_spi_tot"//char(0))

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

    do i=1, n_spi_tot
      spi_R_arr(i)       = pellets(i)%spi_R
      spi_Z_arr(i)       = pellets(i)%spi_Z
      spi_phi_arr(i)     = pellets(i)%spi_phi
      spi_phi_init_arr(i)= pellets(i)%spi_phi_init
      spi_Vel_R_arr(i)   = pellets(i)%spi_Vel_R
      spi_Vel_Z_arr(i)   = pellets(i)%spi_Vel_Z
      spi_Vel_RxZ_arr(i) = pellets(i)%spi_Vel_RxZ
      spi_radius_arr(i)  = pellets(i)%spi_radius
      spi_abl_arr(i)     = pellets(i)%spi_abl
      spi_species_arr(i) = pellets(i)%spi_species
    end do

    call HDF5_array1D_saving(file_id,spi_R_arr, &
             n_spi_tot,'spi_R_arr'//char(0))
    call HDF5_array1D_saving(file_id,spi_Z_arr, &
             n_spi_tot,'spi_Z_arr'//char(0))
    call HDF5_array1D_saving(file_id,spi_phi_arr, &
             n_spi_tot,'spi_phi_arr'//char(0))
    call HDF5_array1D_saving(file_id,spi_phi_init_arr, &
             n_spi_tot,'spi_phi_init_arr'//char(0))
    call HDF5_array1D_saving(file_id,spi_Vel_R_arr, &
             n_spi_tot,'spi_Vel_R_arr'//char(0))
    call HDF5_array1D_saving(file_id,spi_Vel_Z_arr, &
             n_spi_tot,'spi_Vel_Z_arr'//char(0))
    call HDF5_array1D_saving(file_id,spi_Vel_RxZ_arr, &
             n_spi_tot,'spi_Vel_RxZ_arr'//char(0))
    call HDF5_array1D_saving(file_id,spi_radius_arr, &
             n_spi_tot,'spi_radius_arr'//char(0))
    call HDF5_array1D_saving(file_id,spi_abl_arr, &
             n_spi_tot,'spi_abl_arr'//char(0))
    call HDF5_array1D_saving(file_id,spi_species_arr, &
             n_spi_tot,'spi_species_arr'//char(0))

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

    if (spi_tor_rot) call HDF5_real_saving(file_id,ns_phi_rotate,"ns_phi_rotate"//char(0)) 
  else
    n_spi_tot = 0
    call HDF5_integer_saving(file_id,n_spi_tot,"n_spi_tot"//char(0))
  end if

  ! Export restart vacuum 
  call export_HDF5_restart_vacuum(file_id, freeboundary, resistive_wall)

  ! -> clode file
  call HDF5_close(file_id)

  ! -> Deallocate arrays
  call tr_deallocate(t_x,"x",CAT_UNKNOWN)
  call tr_deallocate(t_values,"values",CAT_UNKNOWN)
  call tr_deallocate(t_deltas,"deltas",CAT_UNKNOWN)
#ifdef fullmhd
  call tr_deallocate(t_psi_eq,"psi_eq",CAT_UNKNOWN)
  call tr_deallocate(t_Fprof_eq,"Fprof_eq",CAT_UNKNOWN)
#elif altcs
  call tr_deallocate(t_psi_eq,"psi_eq",CAT_UNKNOWN)
#endif
  call tr_deallocate(t_index,"index",CAT_UNKNOWN)
  call tr_deallocate(t_boundary,"boundary",CAT_UNKNOWN)
  call tr_deallocate(t_axis_node,"axis_node",CAT_UNKNOWN)
  call tr_deallocate(t_parents,"parents",CAT_UNKNOWN)
  call tr_deallocate(t_parent_elem,"parent_elem",CAT_UNKNOWN)
  call tr_deallocate(t_ref_lambda,"ref_lambda",CAT_UNKNOWN)
  call tr_deallocate(t_ref_mu,"ref_mu",CAT_UNKNOWN)
  call tr_deallocate(t_constrained,"constrained",CAT_UNKNOWN)

  call tr_deallocate(t_vertex,"vertex",CAT_UNKNOWN)
  call tr_deallocate(t_neighbours,"neighbours",CAT_UNKNOWN)
  call tr_deallocate(t_size,"size",CAT_UNKNOWN)
  call tr_deallocate(t_father,"father",CAT_UNKNOWN)
  call tr_deallocate(t_n_sons,"n_sons",CAT_UNKNOWN)
  call tr_deallocate(t_n_gen,"n_gen",CAT_UNKNOWN)
  call tr_deallocate(t_sons,"sons",CAT_UNKNOWN)
  call tr_deallocate(t_contain_node,"contain_node",CAT_UNKNOWN)
  call tr_deallocate(t_nref,"nref",CAT_UNKNOWN)

  if (index_now .gt. 0) then
     call tr_deallocate(t_xtime,"xtime",CAT_UNKNOWN)
     call tr_deallocate(t_energies,"energies",CAT_UNKNOWN)
#ifdef JECCD
     call tr_deallocate(t_energies2,"energies2",CAT_UNKNOWN)
     call tr_deallocate(t_energies3,"energies3",CAT_UNKNOWN)
#ifdef JEC2DIAG
     call tr_deallocate(t_energies4,"energies4",CAT_UNKNOWN)
#endif
#endif
  end if

#else
  write (6,*) " ERROR: trying to export with hdf5 but USE_HDF5 was not set at compile-time"
#endif

  return
end subroutine export_hdf5_restart
end module mod_export_restart

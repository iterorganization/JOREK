module mod_direct_construction
#ifdef DIRECT_CONSTRUCTION

  implicit none
  public update_pc_mat

contains

  subroutine update_pc_mat(pc, a_mat, mhd_sim)

    use mod_parameters, only : n_tor, n_var
    use mpi_mod
    use mod_integer_types
    use data_structure, only: type_SP_MATRIX, type_PRECOND, type_RHS
    use mod_simulation_data, only: type_MHD_SIM
    use construct_matrix_mod, only: construct_matrix
    
    implicit none
    
    type(type_PRECOND)                 :: pc
    type(type_SP_MATRIX)               :: a_mat
    type(type_MHD_SIM)                 :: mhd_sim
    
    integer                            :: ierr    
    
    if (.not.pc%structured) call set_pc_structure(pc, a_mat, mhd_sim)
    
    call construct_matrix(mhd_sim, pc%local_elms, pc%n_local_elms, pc%mat, pc%rhs, harmonic_matrix=.true.)
    
  end subroutine update_pc_mat
  
  subroutine set_pc_structure(pc, a_mat, mhd_sim)

    use tr_module
    use mod_integer_types
    use mod_parameters, only : n_tor, n_var
    use data_structure, only: type_SP_MATRIX, type_PRECOND
    use mod_simulation_data, only: type_MHD_SIM
    use mod_global_matrix_structure, only: global_matrix_structure
    use global_distributed_matrix, only: global_matrix_structure_vacuum
    
    implicit none
    
    type(type_PRECOND)                 :: pc
    type(type_MHD_SIM)                 :: mhd_sim
    type(type_SP_MATRIX)               :: a_mat
    integer                            :: i_tor_min, i_tor_max
    integer                            :: i
    integer(kind=int_all)              :: ind
    
    if (pc%my_id.eq.0) write(*,*) "Analyzing preconditioner"
    
    pc%mat%ng = (pc%mode_set_n)*(a_mat%ng)/n_tor ! rank of local PC matrix
    pc%mat%nr = pc%mat%ng
    pc%mat%nc = pc%mat%ng
    pc%n_glob = a_mat%ng
    
    i_tor_min = pc%mode_set(1)
    i_tor_max = pc%mode_set(pc%mode_set_n)
    
    call tr_allocatep(pc%local_elms,1,mhd_sim%element_list%n_elements,"local_elms_harm",CAT_FEM)
    
    call distribute_nodes_elements(mhd_sim%my_id, pc%n_cpu_n, mhd_sim%n_cpu, mhd_sim%node_list, mhd_sim%element_list, .true., pc%local_elms, & 
                                   pc%n_local_elms, mhd_sim%restart, mhd_sim%freeboundary, pc%mat)
                                   
    call global_matrix_structure(mhd_sim%node_list, mhd_sim%element_list, mhd_sim%bnd_elm_list, mhd_sim%freeboundary, &
                                 pc%local_elms, pc%n_local_elms, pc%mat, i_tor_min=i_tor_min, i_tor_max=i_tor_max)
                                 
    if (mhd_sim%freeboundary .and. (mhd_sim%sr_n_tor /= 0)) then 
      call global_matrix_structure_vacuum(mhd_sim%node_list, mhd_sim%bnd_node_list, pc%mat, i_tor_min=i_tor_min, i_tor_max=i_tor_max) 
    endif
    
    allocate(pc%rhs%val(pc%mat%ng))
    pc%rhs%val(1:pc%mat%ng) = 0.d0
    pc%rhs%n = pc%mat%ng

    allocate(pc%row_index(pc%rhs%n))
    
    do ind = 0, a_mat%ng/mhd_sim%n_tor - 1
      do i = 1, pc%mode_set_n
        pc%row_index(i + ind*pc%mode_set_n) =  pc%mode_set(i) + ind*mhd_sim%n_tor
      enddo
    enddo
    
    write(*,*) pc%my_id, "PC matrix: ng, nnz", pc%mat%ng, pc%mat%nnz
    
    pc%structured = .true.    
    
  end subroutine set_pc_structure

#endif
end module mod_direct_construction

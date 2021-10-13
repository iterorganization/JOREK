module construct_matrix_mod

use mod_parameters, only : n_var

implicit none

logical  :: difference_found, rhs_problem(n_var), elm_problem(n_var,n_var)

contains

  !> subroutine that will construct elementary matrices
  subroutine elementary_matrix_build(element, nodes, xpoint2, xcase2, R_axis,         &
       &                             Z_axis, psi_axis, psi_bnd, R_xpoint, Z_xpoint,   &
       &                             omp_tid, ife, n_local_elms, node_list, i_tor_min, i_tor_max, &
                                     aux_nodes)

    ! --- Modules
    use mod_parameters,           only : n_tor, jorek_model, n_vertex_max, n_order, unified_element_matrix
    use phys_module,              only : bc_natural_open, bc_natural_flux, n_tor_fft_thresh, grid_to_wall, n_wall_blocks, keep_n0_const
    USE data_structure,           only : type_element, type_node, type_node_list, thread_struct
    use mod_boundary_matrix_open, only : boundary_matrix_open
    use mod_elt_matrix,           only : element_matrix
    use mod_elt_matrix_fft,       only : element_matrix_fft
    use mod_locate_irn_jcn
    use mod_global_matrix_structure
    use mpi_mod

    ! --- Routine parameters
    type (type_element),              intent(in)     :: element
    type (type_node),                 intent(inout)  :: nodes(n_vertex_max)
    logical,                          intent(in)     :: xpoint2
    integer,                          intent(in)     :: xcase2
    real*8,                           intent(in)     :: R_axis
    real*8,                           intent(in)     :: Z_axis
    real*8,                           intent(in)     :: psi_axis
    real*8,                           intent(in)     :: psi_bnd
    real*8,                           intent(in)     :: R_xpoint(2)
    real*8,                           intent(in)     :: Z_xpoint(2)
    integer,                          intent(in)     :: omp_tid
    integer,                          intent(in)     :: ife
    integer,                          intent(in)     :: n_local_elms
    integer,                          intent(in)     :: i_tor_min   
    integer,                          intent(in)     :: i_tor_max   
    TYPE (type_node_list),            intent(in)     :: node_list
    type (type_node), optional,       intent(inout)  :: aux_nodes(n_vertex_max)
    
    ! -- internal parameters
    integer :: iv, iv2, iv3, iv4, inode1, inode2, inode3, inode4, i, j
    integer :: vertex(2), direction(2), bnd1, bnd2, side1, side2
    integer :: i_max   ! for keep_n0_const max index which should be updated
    integer :: n_tor_local

#ifdef COMPARE_ELEMENT_MATRIX
    integer  :: jvertex, jorder, jvar, jtor, ivertex, iorder, ivar, itor
    integer  :: my_id, rank, ierr

    ! --- Determine ID of each MPI proc
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
    my_id = rank
#endif

    ! --- Call element_matrix
    if ( ( (i_tor_min .eq. 1) .and. (i_tor_max .eq. n_tor) .and. (n_tor .ge. n_tor_fft_thresh) )   &
      .or. (unified_element_matrix) ) then
      ! (use the FFT element matrix construction or the unified one in case it has been combined
      ! for the respective model)
      call element_matrix_fft(element,nodes, xpoint2, xcase2, R_axis, Z_axis, psi_axis, psi_bnd,   &
        R_xpoint, Z_xpoint, thread_struct(omp_tid)%ELM, thread_struct(omp_tid)%RHS, omp_tid,       &
        thread_struct(omp_tid)%ELM_p, thread_struct(omp_tid)%ELM_n, thread_struct(omp_tid)%ELM_k,  &
        thread_struct(omp_tid)%ELM_kn, thread_struct(omp_tid)%RHS_p, thread_struct(omp_tid)%RHS_k, &
        thread_struct(omp_tid)%eq_g, thread_struct(omp_tid)%eq_s, thread_struct(omp_tid)%eq_t,     &
        thread_struct(omp_tid)%eq_p, thread_struct(omp_tid)%eq_ss, thread_struct(omp_tid)%eq_st,   &
        thread_struct(omp_tid)%eq_tt, thread_struct(omp_tid)%delta_g,                              &
        thread_struct(omp_tid)%delta_s, thread_struct(omp_tid)%delta_t, i_tor_min, i_tor_max,      &
        aux_nodes)
    else
      ! (use the element matrix by toroidal integration in case of very few harmonics or in case
      ! of direct construction of the harmonic matrices used in preconditioning)
      call element_matrix(element,nodes, xpoint2, xcase2, R_axis, Z_axis, psi_axis, psi_bnd,       &
        R_xpoint, Z_xpoint, thread_struct(omp_tid)%ELM, thread_struct(omp_tid)%RHS, omp_tid,       &
        i_tor_min, i_tor_max, aux_nodes)
    endif
    
    ! --- Apply sheath boundary conditions at the targets
    if (bc_natural_open) then
      ! --- Loop over the 4 nodes
      do iv = 1, n_vertex_max
        
        iv2  = mod(iv, n_vertex_max) + 1
        iv3 = mod(iv2, n_vertex_max) + 1
        iv4 = mod(iv3, n_vertex_max) + 1
        
        inode1 = element%vertex(iv)
        inode2 = element%vertex(iv2)
        inode3 = element%vertex(iv3)
        inode4 = element%vertex(iv4)
        
        bnd1 = node_list%node(inode1)%boundary
        bnd2 = node_list%node(inode2)%boundary
        
        ! --- carry on only if on boundary
        if ( (bnd1 .eq. 0) .or. (bnd2 .eq. 0)) cycle
        
        nodes(1) = node_list%node(inode1)
        nodes(2) = node_list%node(inode2)
        nodes(3) = node_list%node(inode3)
        nodes(4) = node_list%node(inode4)
        vertex    = (/ iv, iv2 /)
        
        if ( (grid_to_wall) .and. (n_wall_blocks .gt. 0) ) then
          side1 = 0                  ; side2 = 0
          if (bnd1 .eq. 1) side1 = 2 ; if (bnd2 .eq. 1) side2 = 2
          if (bnd1 .eq.11) side1 = 2 ; if (bnd2 .eq.11) side2 = 2
          if (bnd1 .eq. 5) side1 = 3 ; if (bnd2 .eq. 5) side2 = 3
          if (bnd1 .eq.15) side1 = 3 ; if (bnd2 .eq.15) side2 = 3
          if (bnd1 .eq. 2) side1 = 3 ; if (bnd2 .eq. 2) side2 = 3
          if (bnd1 .eq.12) side1 = 2 ; if (bnd2 .eq.12) side2 = 2
          if (bnd1 .eq. 4) side1 = 2 ; if (bnd2 .eq. 4) side2 = 2
          if     ( (side1 .eq. 2) .or. (side2 .eq. 2) ) then
            direction = (/  1, 2  /)
          elseif ( (side1 .eq. 3) .or. (side2 .eq. 3) ) then
            direction = (/  1, 3  /)
          endif
          ! --- This should never happen, but just in case...
          if (     ((side1 .eq. 2) .and. (side2 .eq. 3)) &
              .or. ((side1 .eq. 3) .and. (side2 .eq. 2)) ) then
            write(*,'(A,4i8)') 'WARNING: boundary_matrix_open, boundary element incoherent ',&
                               inode1,node_list%node(inode1)%boundary,inode2,node_list%node(inode2)%boundary  
            cycle
          endif
        else
          ! --- The target has boundary 1 or 3
          if (     (  ((bnd1 .eq. 1) .or. (bnd1 .eq. 3)) .and. ((bnd2 .eq. 1) .or. (bnd2 .eq. 3))  ) &
              .or. (  ((bnd1 .eq. 1) .or. (bnd1 .eq. 9)) .and. ((bnd2 .eq. 1) .or. (bnd2 .eq. 9))  ) &
              .or. (  ((bnd1 .eq. 4) .or. (bnd1 .eq. 9)) .and. ((bnd2 .eq. 4) .or. (bnd2 .eq. 9))  ) &
              .or. (  ((bnd1 .eq. 1) .or. (bnd1 .eq. 4)) .and. ((bnd2 .eq. 4) .or. (bnd2 .eq. 1))  ) ) then
            
            direction = (/  1, 2  /)
            
          elseif (  ((bnd1 .eq. 5) .or. (bnd1 .eq. 9)) .and. ((bnd2 .eq. 5) .or. (bnd2 .eq. 9)) ) then
            
            direction = (/  1, 3  /)
            
          elseif (  ((bnd1 .eq. 2) .or. (bnd1 .eq. 3)) .and. ((bnd2 .eq. 2) .or. (bnd2 .eq. 3)) ) then
            
            direction = (/  1, 3  /)
            
          else
            write(*,'(A,4i8)') 'WARNING: boundary_matrix_open, boundary element not included ',&
                               inode1,node_list%node(inode1)%boundary,inode2,node_list%node(inode2)%boundary  
            cycle
          endif
        endif
          

        ! --- Build matrix elements for boundary
        call boundary_matrix_open(vertex, direction, element, nodes, & 
                                  xpoint2, xcase2, R_axis, Z_axis, psi_axis, psi_bnd, R_xpoint, Z_xpoint, &
                                  thread_struct(omp_tid)%ELM, thread_struct(omp_tid)%RHS, i_tor_min, i_tor_max)
       
      enddo
    endif
    
   
    n_tor_local = i_tor_max - i_tor_min + 1
    ! --- If keep_n0_const then the n0 component should be frozen = diagonal entries high
    if ( keep_n0_const ) then
      i_max =  (n_order+1)*n_vertex_max*n_var*n_tor_local
#ifdef JECCD
      ! n0 component of eccd current should not be frozen when keep_n0_const=.t. (last variable)
      i_max = (n_order+1)*n_vertex_max*(n_var-1)*n_tor_local
#endif
      do i = 1, i_max, n_tor_local
        thread_struct(omp_tid)%ELM(i,i) = 1.d15
      enddo
    endif
    
    ! --- Compare the two element_matrix routines (error thresholds might need to be adapted!)
#ifdef COMPARE_ELEMENT_MATRIX
    ! --- Comparison is performed only for one finite element
    if (ife .eq. n_local_elms/2) then

      ! --- Call both routines
      if (     (jorek_model .eq. 303) &
          .or. (jorek_model .eq. 333) &
          .or. (jorek_model .eq. 500) &
          .or. (jorek_model .eq. 710) &
          .or. (jorek_model .eq. 711) &
          .or. (jorek_model .eq. 712) &
         ) n_tor_fft_thresh = 1
      call element_matrix_fft(element,nodes, xpoint2, xcase2, R_axis, Z_axis, psi_axis, psi_bnd, R_xpoint, Z_xpoint, &
        thread_struct(omp_tid)%ELM2, thread_struct(omp_tid)%RHS2, omp_tid, &
        thread_struct(omp_tid)%ELM_p, thread_struct(omp_tid)%ELM_n, thread_struct(omp_tid)%ELM_k, thread_struct(omp_tid)%ELM_kn, &
        thread_struct(omp_tid)%RHS_p, thread_struct(omp_tid)%RHS_k,  thread_struct(omp_tid)%eq_g, thread_struct(omp_tid)%eq_s, &
        thread_struct(omp_tid)%eq_t, thread_struct(omp_tid)%eq_p, thread_struct(omp_tid)%eq_ss, thread_struct(omp_tid)%eq_st, &
        thread_struct(omp_tid)%eq_tt, thread_struct(omp_tid)%delta_g, thread_struct(omp_tid)%delta_s, &
        thread_struct(omp_tid)%delta_t, i_tor_min, i_tor_max)
      if (     (jorek_model .eq. 303) &
          .or. (jorek_model .eq. 333) &
          .or. (jorek_model .eq. 500) &
          .or. (jorek_model .eq. 710) &
          .or. (jorek_model .eq. 711) &
          .or. (jorek_model .eq. 712) &
         ) n_tor_fft_thresh = 300
      call element_matrix    (element,nodes, xpoint2, xcase2, R_axis, Z_axis, psi_axis, psi_bnd, R_xpoint, Z_xpoint, &
        thread_struct(omp_tid)%ELM,  thread_struct(omp_tid)%RHS,  omp_tid, i_tor_min, i_tor_max)
      
      ! --- Compare right hand side
      write(*,*)
      write(*,*) 'Comparing rhs:'
      write(*,*)

      write(*,'(A)') '  #    my_id       i    ivtx iodr itor         ivar' // &
                     '                    RHS            RHS2        RHS-RHS2'
      do i = 1, n_tor_local*n_vertex_max*(n_order+1)*n_var
    	
    	if (abs(thread_struct(omp_tid)%RHS(i)-thread_struct(omp_tid)%RHS2(i)) / &
            (abs(thread_struct(omp_tid)%RHS(i))+abs(thread_struct(omp_tid)%RHS2(i))+1.d0) .gt. 1.d-12) then
    	  call decrypt_index(i, ivertex, iorder, itor, ivar)
    	  write(*,'(4x,2i8,4x,3i4,7x,1i8,7x,3es16.8)') my_id, i, ivertex, iorder, ivar, itor, thread_struct(omp_tid)%RHS(i), &
    	      thread_struct(omp_tid)%RHS2(i), thread_struct(omp_tid)%RHS(i)-thread_struct(omp_tid)%RHS2(i)
    	  rhs_problem(ivar) = .true.
    	  difference_found  = .true.
    	endif
    	
      enddo
      
      ! --- Compare matrix entries
      write(*,*)
      write(*,*) 'Comparing elm:'
      write(*,*)
      write(*,'(A)') '  #    my_id       i       j    ivtx iodr itor      ivar       jvtx jodr jtor      jvar' // &
                     '                    ELM            ELM2        ELM-ELM2'
      do i = 1, n_tor_local*n_vertex_max*(n_order+1)*n_var
    	do j = 1, n_tor_local*n_vertex_max*(n_order+1)*n_var
    	  
    	  if (abs(thread_struct(omp_tid)%ELM(i,j)-thread_struct(omp_tid)%ELM2(i,j))/  &
    	      (abs(thread_struct(omp_tid)%ELM(i,j))+abs(thread_struct(omp_tid)%ELM2(i,j))+1.d0) .gt. 1.d-10) then
    	    call decrypt_index(i, ivertex, iorder, ivar, itor)
    	    call decrypt_index(j, jvertex, jorder, jvar, jtor)
    	    write(*,'(4x,3i8,4x,3i4,4x,1i8,7x,3i4,4x,1i8,7x,3es16.8)') my_id, i, j, ivertex, iorder, itor, ivar, &
                                                                                    jvertex, jorder, jtor, jvar, &
                                                                       thread_struct(omp_tid)%ELM(i,j), thread_struct(omp_tid)%ELM2(i,j), &
    	      thread_struct(omp_tid)%ELM(i,j)-thread_struct(omp_tid)%ELM2(i,j)
    	    elm_problem(ivar,jvar) = .true.
    	    difference_found	   = .true.
    	  endif
    	  
        enddo
      enddo
      
    endif
#endif /* End of element_matrix comparison */
    
  end subroutine elementary_matrix_build



!> Construct the main matrix from the contributions of the Bezier elements.
!!
!! The element contributions are determined by element_matrix(_fft). Additional
!! contributions from boundary conditions and the free boundary extension are
!! added by external routine calls.
subroutine construct_matrix(my_id, MPI_COMM_N, my_id_n, MPI_COMM_MASTER, my_id_master, local_elms, n_local_elms, index_min, index_max,& 
                            xpoint2, xcase2, R_axis, Z_axis, psi_axis, psi_bnd, R_xpoint, Z_xpoint, psi_xpoint,i_tor_min, i_tor_max,  &
                            n, nz, ndof, n_matrix_block_size, A_mat, rhs, irn, jcn, ijA_index, ijA_size, irn_jcn, harmonic_matrix)
  
  use mumps_module 
  use tr_module 
  use mod_parameters
  use data_structure
  use phys_module
  use pellet_module
  use nodes_elements
  use vacuum, only: sr
  use vacuum_response, only: vacuum_boundary_integral
  use mod_ch_nod_rhs_elm
  use mod_boundary_matrix_open
  use mod_elt_matrix
  use mod_elt_matrix_fft
  use mpi_mod
  use mod_boundary_conditions, only : boundary_conditions
  use mod_fix_axis_nodes, only : fix_nodes_on_axis
  use mod_locate_irn_jcn
  use mod_integer_types
  !$ use omp_lib
  implicit none
  
#include "r3_info.h"

  ! --- Routine parameters
  integer,               intent(in) :: my_id
  integer,               intent(in) :: MPI_COMM_N
  integer,               intent(in) :: my_id_n
  integer,               intent(in) :: MPI_COMM_MASTER
  integer,               intent(in) :: my_id_master
  integer,               intent(in) :: local_elms(*)
  integer,               intent(in) :: n_local_elms
  integer,               intent(in) :: index_min
  integer,               intent(in) :: index_max
  integer,               intent(in) :: xcase2
  real*8,                intent(in) :: R_axis
  real*8,                intent(in) :: Z_axis
  real*8,                intent(in) :: psi_axis
  real*8,                intent(in) :: psi_bnd
  real*8,                intent(in) :: R_xpoint(2)
  real*8,                intent(in) :: Z_xpoint(2)
  real*8,                intent(in) :: psi_xpoint(2)
  logical,               intent(in) :: xpoint2
  integer,               intent(in) :: i_tor_min
  integer,               intent(in) :: i_tor_max
  integer(kind=int_all), intent(in) :: n
  integer(kind=int_all), intent(in) :: nz
  integer(kind=int_all), intent(in) :: ndof
  integer(kind=int_all), intent(in) :: n_matrix_block_size
  logical,               intent(in) :: harmonic_matrix
  real*8,                intent(inout), allocatable :: A_mat(:)
  real*8,                intent(inout), allocatable :: rhs(:)
  integer(kind=int_all), intent(inout), allocatable :: irn(:)
  integer(kind=int_all), intent(inout), allocatable :: jcn(:)
  integer(kind=int_all), intent(in),    allocatable :: ijA_index(:,:), ijA_size(:), irn_jcn(:,:)
  
  !--- Internal variables
  type (type_element)               :: element
  type (type_node)                  :: nodes(n_vertex_max), aux_nodes(n_vertex_max)
  type (type_element)               :: element_father
  type (type_node)                  :: nodes_father(n_vertex_max)
  real*8,              allocatable  :: rhs_local(:)
  integer                           :: i, ife, iv, iv2, inode, inode1, inode2, knode, j, k, l, index_ij, index_kl
  integer(kind=int_all), parameter  :: Int1=1
  integer                           :: index_node1, index_node2, index_min_loc, index_max_loc
  integer(kind=int_all)             :: ijA_position
  integer(kind=int_all)             :: index_large_i, index_large_k, ilarge2
  integer                           :: i_order, k_order, ielm, ierr
  integer                           :: vertex(2), direction(2)
  integer                           :: omp_nthreads, omp_tid, n_tor_local
  integer                           :: node_out(n_vertex_max)
  integer                           :: i_father,INODE_FATHER, ios
  integer                           :: ilarge_vp, in, ivertex, iorder, ivar, itor, jvertex, jorder, jvar, jtor
  integer                           :: random_element, n_var_reduced, v1, v2, im, index_ij_model400_e, index_kl_model400_e
  real*8                            :: tmp_rhs, tmp_elm, tmp_elm_v2_8
  CHARACTER(LEN=128)                :: fname

  ! --- Timing call
  call r3_info_begin (r3_info_index_0, 'construct_matrix')

  ! --- Printout
  if (my_id .eq. 0) then
    if(.NOT.harmonic_matrix) then
      write(*,*) '****************************************'
      write(*,*) '*  construct global matrix             *'
      write(*,*) '****************************************'
    else
      write(*,*) '****************************************'
      write(*,*) '*  construct harmonic matrix           *'
      write(*,*) '****************************************'
    endif
  endif
  
  ! --- Memory tracking
  call tr_print_memsize("DebConstM")

  if (.not. harmonic_matrix) then
    
    ! --- Local min-max indices for the nodes of our local elements (local in the MPI sense)
    do i = 1, n_local_elms
      ielm = local_elms(i)

      do iv = 1, n_vertex_max
        if (ielm > n_elements_max) then
          write(*,*) "WARNING: ielm, n_elements_max = ", ielm, n_elements_max
        end if

        inode = element_list%element(ielm)%vertex(iv)

        if (i == 1 .and. iv == 1) then
          index_min_loc = minval(node_list%node(iv)%index)
          index_max_loc = maxval(node_list%node(iv)%index)
        else
          index_min_loc = min(index_min_loc, minval(node_list%node(iv)%index))
          index_max_loc = max(index_max_loc, maxval(node_list%node(iv)%index))
        end if
      enddo

    enddo

    ! --- Initialise internal variables
    difference_found = .false.
    rhs_problem(:)   = .false.
    elm_problem(:,:) = .false.

  endif ! (.not. harmonic_matrix)

  ! --- Memory allocation
  if (allocated(A_mat))   call tr_deallocate(A_mat,"A_mat",CAT_DMATRIX) 
  call tr_allocate(A_mat,Int1,nz,"A_mat",  CAT_DMATRIX)
  A_mat = 0.0d0

  if (allocated(rhs)) call tr_deallocate(rhs,"rhs",CAT_DMATRIX) 
  call tr_allocate(rhs, Int1,ndof,"rhs", CAT_DMATRIX)
  rhs = 0.0d0 

  call tr_allocate(rhs_local, Int1,ndof,"rhs_local", CAT_DMATRIX)
  rhs_local  = 0.d0
  
  ! --- Declare shared and private variables for omp
  !$omp parallel default(none) &
  !$omp   shared(n_local_elms,local_elms,element_list,node_list, aux_node_list,                                &
  !$omp          index_min, index_max,xpoint2,xcase2,R_axis,Z_axis,psi_axis,psi_bnd,Z_xpoint,harmonic_matrix,  &
  !$omp          A_mat, rhs_local, rhs, irn, jcn,                                                              &
  !$omp          R_xpoint,my_id,bc_natural_open,bc_natural_flux,refinement,thread_struct,n_tor_fft_thresh,     &
  !$omp          difference_found,rhs_problem,elm_problem, i_tor_min, i_tor_max, ijA_index, ijA_size, irn_jcn) &
  !$omp   private(ife,ielm,iv,inode,element,nodes, aux_nodes, i,inode1,i_order,index_node1, n_tor_local,   &
  !$omp           index_large_i,j,index_ij,k,knode,k_order,index_node2,index_large_k,ijA_position,         &
  !$omp           l,index_kl,ilarge2,iv2,vertex,direction,inode2,omp_nthreads,omp_tid,                     &
  !$omp           i_father,element_father, nodes_father, inode_father, node_out, ivertex, iorder,          &
  !$omp           ivar, itor, jvertex, jorder, jvar, jtor, random_element, n_var_reduced, v1, v2, im,      &
  !$omp           index_ij_model400_e, index_kl_model400_e,  tmp_rhs, tmp_elm, tmp_elm_v2_8        )

! --- omp id
#ifdef _OPENMP
  omp_nthreads = omp_get_num_threads()
  omp_tid      = 1+omp_get_thread_num()
#else
  omp_nthreads = 1
  omp_tid      = 1
#endif
  
  n_tor_local = i_tor_max - i_tor_min + 1 
! --- Loop over local elements
  !$omp do schedule(runtime)
  do ife = 1, n_local_elms
    
    ! --- Get element
    ielm = local_elms(ife)
    element = element_list%element(ielm)
    
    ! --- Define nodes (this depends on whether our element has been refined)
    if (refinement .and. .not. harmonic_matrix) then

      i_father = element_list%element(ielm)%father

      if (i_father .ne. 0) then
      element_father = element_list%element(i_father)
      do iv = 1, n_vertex_max
        inode_father=element_father%vertex(iv)
        nodes_father(iv) = node_list%node(inode_father)
      enddo
     endif

    else
       
      do iv = 1, n_vertex_max
        inode         = element%vertex(iv)
        nodes(iv)     = node_list%node(inode)
        aux_nodes(iv) = aux_node_list%node(inode)
      enddo

    endif

    call elementary_matrix_build(element, nodes, xpoint2, xcase2, R_axis, Z_axis, psi_axis,        &
      psi_bnd, R_xpoint, Z_xpoint, omp_tid, ife, n_local_elms, node_list, i_tor_min, i_tor_max, aux_nodes)
    
#ifdef PRINT_ELM_RHS
    if (.not. harmonic_matrix) then
      ! --- Write out rhs and elm for one element to compare models.
      !     Switch this on by adding -DPRINT_ELM_RHS as compiler flag.
      if (ielm == element_list%n_elements/2) then
      
        open(unit = 387, file = 'comp_matrix_elements_elm.dat', status='REPLACE', action='WRITE')
        open(unit = 388, file = 'comp_matrix_elements_rhs.dat', status='REPLACE', action='WRITE')
        write(387, "( '#', A17, 8A4)" ), 'ELM', 'v1', 'v2', 'i', 'k', 'j', 'l', 'im', 'in'
        write(388, "( '#', A17, 4A4)" ), 'RHS', 'v1', 'i', 'j', 'im'
      
        n_var_reduced = n_var
#ifdef WITH_TiTe
        n_var_reduced = n_var - 1
#endif
      
        do v1 = 1, n_var_reduced
        
          do v2 = 1, n_var_reduced
          
            do i = 1, n_vertex_max
            
              do k = 1, n_vertex_max
              
                do j = 1, n_order+1
                
                  do l = 1, n_order+1
                  
                    do im = 1, n_tor_local
                    
                    
                      do in = 1, n_tor_local
                      
                        ! --- Indices for RHS (index_ij) and ELM matrix (index_ij, index_kl)
                        index_ij            = n_tor_local * n_var * (n_order+1) * (i-1) &
                                            + n_tor_local * n_var * (j-1) + n_tor_local * (v1-1) + im
                      
                        index_kl            = n_tor_local * n_var * (n_order+1) * (k-1) &
                                            + n_tor_local * n_var * (l-1) + n_tor_local * (v2-1) + in
                      
                        ! --- Indices for T_e (model400)
                        index_ij_model400_e = n_tor_local * n_var * (n_order+1) * (i-1) &
                                            + n_tor_local * n_var * (j-1) + n_tor_local * (var_Te -1) + im
                        index_kl_model400_e = n_tor_local * n_var * (n_order+1) * (k-1) &
                                            + n_tor_local * n_var * (l-1) + n_tor_local * (var_Te -1) + in

                        !--- RHS: simple output of vector element
                        if ( (k==1) .and. (l==1) .and. (v2==1) .and. (in==1) ) then
                        
                          tmp_rhs = thread_struct(omp_tid)%RHS(index_ij)
                        
#ifdef WITH_TiTe
                          !--- RHS: for model400, add T_e (v1=8) to T_i (v1=6)
                          if (v1 == var_Ti) tmp_rhs = tmp_rhs + thread_struct(omp_tid)%RHS(index_ij_model400_e)
#endif
                        
                          write(388, "( E18.6, 4I4 )" ) tmp_rhs, v1, i, j, im
                      
                        end if


                        ! --- ELM: simple output of matrix element
                        tmp_elm   = thread_struct(omp_tid)%ELM(index_ij,index_kl)

                        ! --- ELM: additional output for comparison with model400, output as v2=8, below v2=6
                        ! --- for model != model400: duplicate of ELM(v1, v2=6)
                        ! --- for model  = model400: contribution from ELM(v1, v2=8)
                        tmp_elm_v2_8 = tmp_elm

                        ! --- for model400, when v1==6, add the ELM(v1=8, v2) contribution to ELM(v1=6, v2), for both tmp_elm and tmp_elm_v2_8

#ifdef WITH_TiTe
                        if (v2 == var_Ti ) then
                          tmp_elm_v2_8 = thread_struct(omp_tid)%ELM(index_ij, index_kl_model400_e) 
                          if (v1 == var_Ti) then
                            tmp_elm_v2_8 = tmp_elm_v2_8 + thread_struct(omp_tid)%ELM(index_ij_model400_e, index_kl_model400_e) 
                          end if
                        end if
                      
                        if (v1 == var_Ti) then
                          tmp_elm = tmp_elm + thread_struct(omp_tid)%ELM(index_ij_model400_e, index_kl)
                        end if
#endif


                        write(387, "( E18.6, 8I4 )" )   tmp_elm,      v1, v2, i, k, j, l, im, in
                      
                        if (v2 == 6) then
                          write(387, "( E18.6, 8I4 )" ) tmp_elm_v2_8, v1, 8,  i, k, j, l, im, in
                        end if
                      
                      end do ! n_tor_local
                    
                    end do ! n_tor_local
                  
                  end do ! n_order+1
                
                end do ! n_order+1
              
              end do ! n_vertex_max
            
            end do ! n_vertex_max
          
          end do ! n_var_reduced
        
        end do ! n_var_reduced
      
        close(387)
        close(388)
      end if ! i_elm
    ! --- end: Write out rhs and elm for one element to compare models.
    endif ! .not. harmonic_matrix  
#endif
 
    ! --- Define element nodes (depends if it's refined)
    if (refinement .and. .not. harmonic_matrix) then   
      call ch_nod_rhs_elm(ielm,element,nodes,element_father,nodes_father, &
              thread_struct(omp_tid)%ELM, thread_struct(omp_tid)%RHS,node_out)
    else
      do i=1, n_vertex_max
        node_out(i) = element%vertex(i)   
      enddo 
    endif

    ! --- We only look at non-refined elements
    if ((.not. refinement) .or. (refinement .and. (element%n_sons .eq. 0))) then
    
      do i=1,n_vertex_max

        inode1 = node_out(i)

        do i_order = 1, n_order+1

          index_node1 = node_list%node(inode1)%index(i_order)

          index_large_i = n_tor_local * n_var * (index_node1 - 1)

          if ((index_node1 .ge. index_min) .and. (index_node1 .le. index_max)) then

            do j = 1, n_var * n_tor_local

              index_ij = n_tor_local * n_var * (n_order+1) * (i-1) + n_tor_local * n_var * (i_order-1) + j   ! index in the ELM matrix
             
              !$omp atomic
              rhs_local(index_large_i+j) = rhs_local(index_large_i+j) + thread_struct(omp_tid)%RHS(index_ij) 
              !$omp end atomic
            enddo

            do k=1,n_vertex_max

              knode = node_out(k)

              do k_order = 1, n_order+1

                index_node2 = node_list%node(knode)%index(k_order)

                index_large_k = n_tor_local * n_var * (index_node2 - 1)

                call locate_irn_jcn(index_node1,index_node2,index_min,index_max,ijA_position,ijA_index, ijA_size, irn_jcn)

                thread_struct(omp_tid)%synch_buff(1:n_var*n_tor_local*n_var*n_tor_local) = 0.d0
                do j = 1, n_var * n_tor_local
                  index_ij = n_tor_local * n_var * (n_order+1) * (i-1) + n_tor_local * n_var * (i_order-1) + j   ! index in the ELM matrix

                  do l = 1, n_var * n_tor_local

                    index_kl = n_tor_local * n_var * (n_order+1) * (k-1) +  n_tor_local * n_var * (k_order-1) + l   ! index in the ELM matrix

                    ilarge2 = ijA_position - 1 + (j-1) * n_var * n_tor_local + l

                    irn(ilarge2) = index_large_i	+ j
                    jcn(ilarge2) = index_large_k	+ l

                    thread_struct(omp_tid)%synch_buff((j-1)*n_var*n_tor_local+l) = &
                      thread_struct(omp_tid)%synch_buff((j-1)*n_var*n_tor_local+l) + thread_struct(omp_tid)%ELM(index_ij,index_kl)
                     
                  enddo ! n_var * n_tor_local

                enddo ! n_var * n_tor_local

                !$omp critical
                A_mat(ijA_position : ijA_position + n_var*n_tor_local*n_var*n_tor_local - 1) = &
                  A_mat(ijA_position : ijA_position + n_var*n_tor_local*n_var*n_tor_local - 1) +  &
                  thread_struct(omp_tid)%synch_buff(1:n_var*n_tor_local*n_var*n_tor_local)
                !$omp end critical 

              enddo ! n_order+1
            enddo ! n_vertex_max

          endif ! index_min < index < index_max

        enddo ! n_order+1

      enddo ! n_vertex_max
      
    end if

  end do
  !$omp end do
  !$omp end parallel
  
  ! --- Memory tracking
  call tr_vnorms("cm_A_bef_bc",A_mat,nz)
  
  ! --- Apply boundary conditions.
  call boundary_conditions(my_id, node_list, element_list,  bnd_node_list,local_elms, n_local_elms,  &
                           index_min, index_max,  rhs_local, xpoint2, xcase2, R_axis, Z_axis,        & 
                           psi_axis, psi_bnd, R_xpoint, Z_xpoint, psi_xpoint, .false., .false.,      & 
                           ijA_index, ijA_size, irn_jcn, irn, jcn, A_mat, i_tor_min, i_tor_max )

  if (fix_axis_nodes) then
    call fix_nodes_on_axis(node_list, element_list, local_elms, n_local_elms, index_min, index_max, & 
                           ijA_index, ijA_size, irn_jcn, irn, jcn, A_mat, i_tor_min, i_tor_max )
  endif

  ! --- Memory tracking
  call tr_vnorms("cm_A_aft_bc",A_mat,nz)


    ! --- Add vacuum response (boundary integral) for free boundary computations
    if ( freeboundary .and. ( sr%n_tor /= 0 ) ) then
      call vacuum_boundary_integral(my_id, bnd_node_list, node_list, bnd_elm_list, freeboundary_equil, & 
                                  resistive_wall, index_min, index_max, rhs_local, A_mat, tstep, index_now, & 
                                  irn, jcn, n_matrix_block_size, ijA_index, ijA_size, irn_jcn, i_tor_min, i_tor_max)
    end if
  
  if ( .not. harmonic_matrix ) then 
    ! --- Summarise element_matrix comparison
#ifdef COMPARE_ELEMENT_MATRIX
    call MPI_Barrier(MPI_COMM_WORLD,ierr)
    if ( difference_found ) then
      write(*,*)
      write(*,'(i3,a)') my_id, ' ERROR: DIFFERENCES BETWEEN ELEMENT_MATRIX AND ELEMENT_MATRIX_FFT!'
      do i = 1, n_var
        if ( rhs_problem(i) ) write(*,'(i5," rhs_ij_",i1)') my_id, i
      end do
      write(*,*)
      do i = 1, n_var
        do j = 1, n_var
          if ( elm_problem(i,j) ) write(*,'(i5," amat_",2i1)') my_id, i, j
        end do
      end do
      write(*,*)
      write(*,*)
    else
      write(*,*)
      write(*,'(i3,a)') my_id, ' BOTH ELEMENT_MATRIX ROUTINES SEEM TO BE CONSISTENT.'
      write(*,*)
    end if
    call MPI_Barrier(MPI_COMM_WORLD,ierr)
    if ( difference_found ) stop
#endif

#ifdef NORMTRACE
    ! --- For debugging purpose
    call MPI_Reduce(RHS_local,RHS,ndof,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    call tr_locvnorms("cm_Rhs",RHS,ndof)
    if (my_id .eq. 0) then
      write(fname,'(A,I6.6)')"rhs",index_now
      call tr_vdump(fname,RHS,ndof)
    end if
#endif
    call MPI_Reduce(RHS_local,RHS,ndof,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)

  else ! ( if harmonic_matrix)
  
    ! --- Form a global rhs from the rhss of the individual mpi tasks
    call MPI_Reduce(RHS_local,RHS,ndof,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_N,ierr)

  endif 

  call tr_deallocate(RHS_local,"RHS_local",CAT_DMATRIX)
     
  ! --- Memory tracking
  call tr_locvnorms("cm_BCRhs",RHS,ndof)
  call tr_debug_write("ndof",ndof)
  
  ! --- Timing
  call r3_info_end(r3_info_index_0)
  call tr_print_memsize("EndConstM")

end subroutine construct_matrix


!> Helps to interprete an element matrix index
subroutine  decrypt_index(ind, ivertex, iorder, ivar, itor)

  use mod_parameters,  only : n_tor, jorek_model, n_vertex_max, n_var, n_order 

  integer, intent(in)  :: ind     !< Element matrix index
  integer, intent(out) :: ivertex !< Vertex index
  integer, intent(out) :: iorder  !< Degree of freedom
  integer, intent(out) :: ivar    !< Variable index
  integer, intent(out) :: itor    !< Toroidal mode index

  integer :: ind2

  ind2 = ind

  ivertex = ( ind2 - 1 ) / ( n_tor*n_var*(n_order+1) ) + 1
  ind2 = ind2 - ( ivertex - 1 ) * ( n_tor*n_var*(n_order+1) )

  iorder = ( ind2 - 1 ) / ( n_tor*n_var ) + 1
  ind2 = ind2 - ( iorder - 1 ) * ( n_tor*n_var )

  ivar = ( ind2 - 1 ) / ( n_tor ) + 1
  ind2 = ind2 - ( ivar - 1 ) * ( n_tor )

  itor = ind2

end subroutine decrypt_index


end module construct_matrix_mod


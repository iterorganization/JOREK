!> Definitions of derived data types for grid nodes and elements, boundary nodes and elements,
!! and flux surface elements, as well as the shattered pellets
module data_structure
  use mod_parameters
  use tr_module
  use gauss
  use ISO_C_BINDING, ONLY : C_INT

  implicit none

  type type_node                                  !< type definition of a node (i.e. a vertex)
    real*8     :: x(n_coord_tor,n_order+1,n_dim)        !< x,y,z coordinates of points and additional nodal geometry
    real*8     :: values(n_tor,n_order+1,n_var)   !< Variable values and derivatives
    real*8     :: deltas(n_tor,n_order+1,n_var)   !< Change of variable values and derivatives in last timestep
#ifdef fullmhd
    real*8     :: psi_eq(n_order+1)               !< equilibrium flux at the nodes
    real*8     :: Fprof_eq(n_order+1)             !< equilibrium profile R*B_phi at the nodes
#elif altcs
    real*8     :: psi_eq(n_order+1)               !< equilibrium flux at the nodes
#endif
    integer    :: index(n_order+1)                !< index in the main matrix
    integer    :: boundary                        !< = 1, 2 or 3 for boundary nodes.
                                                  !< For wall-aligned grids, check routine update_boundary_types_final
                                                  !< in grids/grid_utils/update_boundary_types.f90
    integer    :: boundary_index                  !< index of the boundary node 
    logical    :: axis_node                       !< Flag nodes that are on the axis (and can/need-to-be be stabilised)
    integer    :: axis_dof                        !< which dof to enforce to zero
    integer    :: parents(2)                      !< Parent nodes (used if node is constrained)"refinement"
    integer    :: parent_elem                     !< which element do parent nodes belong to ? "refinement"
    real*8     :: ref_lambda, ref_mu              !< Local coordinates of node inside the parent element. "refinement"
    logical    :: constrained                     !< Constrained node or not..."refinement"
  end type type_node

  type type_node_list                             !< type definition of a list of nodes
    integer            :: n_nodes                 !< the number of nodes in the list
    integer            :: n_dof                   !< the total number of degrees of freedom
    type (type_node)   :: node(n_nodes_max)       !< an allocatable list of nodes
  end type type_node_list

  type type_element                               !< type definition for one elements
    integer :: vertex(n_vertex_max)               !< nodes of the corners
    integer :: neighbours(n_vertex_max)           !< neighbouring elements
    real*8  :: size(n_vertex_max,n_order+1)       !< size of vectors at each vertex of the element
    integer :: father                             !< index of father element (0 if no father)"refinement"
    integer :: n_sons                             !< Number of sons elements"refinement"
    integer :: n_gen                              !< Generation rank of the element"refinement"
    integer :: sons(4)                            !< Sons of the element (=0 if no son)"refinement"
    integer :: contain_node(5)                    !< nodes belonging to the element"refinement"
    integer :: nref                               !< How the element has been refined (if so)"refinement"
  end type type_element

  type type_element_list                          !< type definition for a list of elements
    integer :: n_elements                         !< number of elements in the list
    type (type_element) :: element(n_elements_max)!< list of elements
  end type type_element_list

  type type_bnd_element                           !< type definition for one boundary element (1D element)
    integer :: vertex(2)                          !< indices of the nodes of the corners in the node list
    integer :: bnd_vertex(2)                      !< indices of the nodes of the corners in the boundary node list
    integer :: direction(2,2)                     !< indicates which direction of the nodes is along the boundary (2 or 3)
    integer :: element                            !< boundary element is part of this element
    integer :: side                               !< boundary element corresponds to this side of the originating element
    real*8  :: size(2,2)                          !< size of vectors at each vertex of the element : size(vertex,order)
  end type type_bnd_element

  type type_bnd_element_list                      !< type definition for a list of boundary elements
    integer :: n_bnd_elements                     !< number of boundary elements in the list
    type (type_bnd_element) :: bnd_element(n_boundary_max) !< list of boundary elements
  end type type_bnd_element_list
  
  type type_bnd_node                              !< type definition for one boundary node
    integer :: index_jorek                        !< index of the node in the node_list
    integer :: index_starwall(2)                  !< index of the node in STARWALL numbering
    integer :: n_dof                              !< total number of degrees of freedom for this boundary node
    integer :: direction(2)                       !< which direction is along the boundary?
  end type type_bnd_node
  
  type type_bnd_node_list                         !< type definition for a list of boundary nodes
    integer :: n_bnd_nodes                        !< number of boundary nodes
    type (type_bnd_node) :: bnd_node(n_boundary_max)    !< list of the boundary nodes
  end type type_bnd_node_list

  type type_surface                               !< type definition for a fluxsurface (in 2D)
    integer :: flag                           	  !< Flag surface if you want (used for wall-grid)
    real*8  :: psi                           	  !< psi-value of the surface
    integer :: n_pieces                           !< total number of pieces (each piece is a 3rd order polynomial)
    integer :: n_parts                            !< number of surface parts (eg. one core surf + on private surf)
    integer :: parts_index(10)                    !< index of the first piece of each surface part (assuming 10 parts max)
    integer :: elm(n_pieces_max)                  !< element containg the current piece
    real*8  :: s(4,n_pieces_max), t(4,n_pieces_max)     !< 4 variables per line piece of the flux surface
   end type type_surface

  type  type_surface_list                         !< type definition for a list of surfaces
    integer                         :: n_psi      !< the number of surfaces
    real*8,allocatable              :: psi_values(:)    !< the values of the poloidal flux at the surfaces
    type (type_surface),allocatable :: flux_surfaces(:) !< the list of surfaces
  end type type_surface_list

  TYPE type_thread_buffer
     real*8, dimension (:,:,:), allocatable :: ELM_p
     real*8, dimension (:,:,:), allocatable:: ELM_n
     real*8, dimension (:,:,:), allocatable:: ELM_k
     real*8, dimension (:,:,:), allocatable:: ELM_kn
     real*8, dimension (:,:)  , allocatable:: RHS_p
     real*8, dimension (:,:)  , allocatable:: RHS_k
     real*8, dimension (:,:)  , allocatable :: ELM
     real*8, dimension (:,:)  , allocatable :: ELM2
     real*8, dimension (:)    , allocatable :: RHS
     real*8, dimension (:)    , allocatable :: RHS2

     real*8, dimension(:,:,:,:) , allocatable :: eq_g, eq_s, eq_t
     real*8, dimension(:,:,:,:) , allocatable :: eq_p, eq_pp
     real*8, dimension(:,:,:,:) , allocatable:: eq_ss, eq_st, eq_tt   
     real*8, dimension(:,:,:,:) , allocatable :: delta_g, delta_s, delta_t

     real*8, dimension(:), allocatable  :: synch_buff
  END TYPE type_thread_buffer

  !> One shard of a shattered pellet (or the complete pellet if unshattered)
  type type_SPI
    real*8  :: spi_R                 !< R coordinate of shard (m)
    real*8  :: spi_Z                 !< Z coordinate of shard (m)
    real*8  :: spi_phi               !< Phi coordinate of shard (radian)
    real*8  :: spi_phi_init          !< The initial phi coordinate of shard (radian) for trajectory calculation.
    real*8  :: spi_Vel_R             !< Velocity in R direction (m/s)
    real*8  :: spi_Vel_Z             !< Velocity in Z direction (m/s)
    real*8  :: spi_Vel_RxZ           !< Velocity in RxZ direction (m/s)
    real*8  :: spi_radius            !< Shard radius (assuming spherical shard) (m)
    real*8  :: spi_abl               !< Shard ablation rate (atom/s)
    real*8  :: spi_species           !< Fraction of impurity atoms relative to the total number of atoms (model501)
                                     !! 0.: pure background species
                                     !! 1.: pure impurity shard
    real*8  :: spi_vol               !< Numerically integrated volume of the gas source at the shard position
    real*8  :: spi_psi               !< Psi value at the shard position
    real*8  :: spi_grad_psi          !< Value of grad(Psi)=sqrt(PSI_R * PSI_R + PSI_Z * PSI_Z) at the shard position

    real*8  :: spi_vol_drift         !< Numerically integrated volume of the gas source depositing at the post-drift position
    real*8  :: spi_psi_drift         !< Psi value at the post-drift deposition position
    real*8  :: spi_grad_psi_drift    !< Value of grad(Psi)=sqrt(PSI_R * PSI_R + PSI_Z * PSI_Z) at the post-drift deposition position

  end type type_SPI
 
  integer                                         , public :: nbthreads
  TYPE(type_thread_buffer), dimension(:), pointer , public :: thread_struct => NULL()
  
contains

  subroutine init_threads()
#ifdef _OPENMP
    use omp_lib
    !$OMP PARALLEL shared(nbthreads)
    !$OMP master
    call omp_set_dynamic(.false.)
    nbthreads = omp_get_num_threads()
    !$OMP end master
    !$OMP barrier
    !$OMP end PARALLEL
#else
    nbthreads = 1
#endif
  end subroutine init_threads

  subroutine new_thread_buffers()
    integer i
    if (.not. associated(thread_struct)) then
       allocate(thread_struct(nbthreads))
       call tr_register_mem(sizeof(thread_struct),"thread_struct",CAT_MATELEM)
       do i = 1, nbthreads
          call tr_debug_write("Init thread_struct, thread_id=",i)
          call tr_allocate(thread_struct(i)%ELM_p, 1,n_plane,1,n_vertex_max*n_var*(n_order+1),1,n_vertex_max*n_var*(n_order+1),"ELM_p",CAT_MATELEM)
          call tr_allocate(thread_struct(i)%ELM_n, 1,n_plane,1,n_vertex_max*n_var*(n_order+1),1,n_vertex_max*n_var*(n_order+1),"ELM_n",CAT_MATELEM)
          call tr_allocate(thread_struct(i)%ELM_k, 1,n_plane,1,n_vertex_max*n_var*(n_order+1),1,n_vertex_max*n_var*(n_order+1),"ELM_k",CAT_MATELEM)
          call tr_allocate(thread_struct(i)%ELM_kn,1,n_plane,1,n_vertex_max*n_var*(n_order+1),1,n_vertex_max*n_var*(n_order+1),"ELM_kn",CAT_MATELEM)
          call tr_allocate(thread_struct(i)%RHS_p, 1,n_plane,1,n_vertex_max*n_var*(n_order+1),"RHS_p",CAT_MATELEM)                                     
          call tr_allocate(thread_struct(i)%RHS_k, 1,n_plane,1,n_vertex_max*n_var*(n_order+1),"RHS_k",CAT_MATELEM)                                     
          call tr_allocate(thread_struct(i)%ELM,   1,n_tor*n_vertex_max*(n_order+1)*n_var,1,n_tor*n_vertex_max*(n_order+1)*n_var,"ELM",CAT_MATELEM)       
          call tr_allocate(thread_struct(i)%RHS,   1,n_tor*n_vertex_max*(n_order+1)*n_var,"RHS",CAT_MATELEM)
          call tr_allocate(thread_struct(i)%synch_buff, 1,n_tor*n_var*n_tor*n_var,"synch_buff",CAT_MATELEM)
          thread_struct(i)%ELM_p   = 0.d0
          thread_struct(i)%ELM_n   = 0.d0
          thread_struct(i)%ELM_k   = 0.d0
          thread_struct(i)%ELM_kn  = 0.d0
          thread_struct(i)%RHS_p   = 0.d0
          thread_struct(i)%RHS_k   = 0.d0
          thread_struct(i)%ELM     = 0.d0
          thread_struct(i)%RHS     = 0.d0
          thread_struct(i)%synch_buff     = 0.d0
          call tr_allocate(thread_struct(i)%eq_g   ,1,n_plane,1,n_var,1,n_gauss,1,n_gauss,"eq_g",CAT_MATELEM)
          call tr_allocate(thread_struct(i)%eq_s   ,1,n_plane,1,n_var,1,n_gauss,1,n_gauss,"eq_s",CAT_MATELEM)
          call tr_allocate(thread_struct(i)%eq_t   ,1,n_plane,1,n_var,1,n_gauss,1,n_gauss,"eq_t",CAT_MATELEM)
          call tr_allocate(thread_struct(i)%eq_p   ,1,n_plane,1,n_var,1,n_gauss,1,n_gauss,"eq_p",CAT_MATELEM)
          call tr_allocate(thread_struct(i)%eq_ss  ,1,n_plane,1,n_var,1,n_gauss,1,n_gauss,"eq_ss",CAT_MATELEM)
          call tr_allocate(thread_struct(i)%eq_st  ,1,n_plane,1,n_var,1,n_gauss,1,n_gauss,"eq_st",CAT_MATELEM)
          call tr_allocate(thread_struct(i)%eq_tt  ,1,n_plane,1,n_var,1,n_gauss,1,n_gauss,"eq_tt",CAT_MATELEM)
          call tr_allocate(thread_struct(i)%eq_pp  ,1,n_plane,1,n_var,1,n_gauss,1,n_gauss,"eq_pp",CAT_MATELEM) 
          call tr_allocate(thread_struct(i)%delta_g,1,n_plane,1,n_var,1,n_gauss,1,n_gauss,"delta_g",CAT_MATELEM) 
          call tr_allocate(thread_struct(i)%delta_s,1,n_plane,1,n_var,1,n_gauss,1,n_gauss,"delta_s",CAT_MATELEM)
          call tr_allocate(thread_struct(i)%delta_t,1,n_plane,1,n_var,1,n_gauss,1,n_gauss,"delta_t",CAT_MATELEM)
          thread_struct(i)%eq_g    = 0.d0
          thread_struct(i)%eq_s    = 0.d0
          thread_struct(i)%eq_t    = 0.d0
          thread_struct(i)%eq_p    = 0.d0
          thread_struct(i)%eq_ss   = 0.d0
          thread_struct(i)%eq_st   = 0.d0
          thread_struct(i)%eq_tt   = 0.d0
          thread_struct(i)%eq_pp   = 0.d0
          thread_struct(i)%delta_g = 0.d0
          thread_struct(i)%delta_s = 0.d0
          thread_struct(i)%delta_t = 0.d0
#ifdef COMPARE_ELEMENT_MATRIX
          call tr_allocate(thread_struct(i)%ELM2,  1,n_tor*n_vertex_max*(n_order+1)*n_var,1,n_tor*n_vertex_max*(n_order+1)*n_var,"ELM2",CAT_MATELEM)
          call tr_allocate(thread_struct(i)%RHS2,  1,n_tor*n_vertex_max*(n_order+1)*n_var,"RHS2",CAT_MATELEM)
          thread_struct(i)%ELM2    = 0.d0
          thread_struct(i)%RHS2    = 0.d0
#endif
       end do
    end if
  end subroutine new_thread_buffers
  
  subroutine del_thread_buffers()
    integer i
    do i = 1, nbthreads
       call tr_deallocate(thread_struct(i)%ELM_p,"ELM_p",CAT_MATELEM)
       call tr_deallocate(thread_struct(i)%ELM_n,"ELM_n",CAT_MATELEM)
       call tr_deallocate(thread_struct(i)%ELM_k,"ELM_k",CAT_MATELEM)
       call tr_deallocate(thread_struct(i)%ELM_kn,"ELM_kn",CAT_MATELEM)
       call tr_deallocate(thread_struct(i)%RHS_p,"RHS_p",CAT_MATELEM)                                     
       call tr_deallocate(thread_struct(i)%RHS_k,"RHS_k",CAT_MATELEM)                                     
       call tr_deallocate(thread_struct(i)%ELM,"ELM",CAT_MATELEM)
       call tr_deallocate(thread_struct(i)%RHS,"RHS",CAT_MATELEM)
       call tr_deallocate(thread_struct(i)%synch_buff,"synch_buff",CAT_MATELEM)
       call tr_deallocate(thread_struct(i)%eq_g   ,"eq_g",CAT_MATELEM)
       call tr_deallocate(thread_struct(i)%eq_s   ,"eq_s",CAT_MATELEM)
       call tr_deallocate(thread_struct(i)%eq_t   ,"eq_t",CAT_MATELEM)
       call tr_deallocate(thread_struct(i)%eq_p   ,"eq_p",CAT_MATELEM)
       call tr_deallocate(thread_struct(i)%eq_ss  ,"eq_ss",CAT_MATELEM)
       call tr_deallocate(thread_struct(i)%eq_st  ,"eq_st",CAT_MATELEM)
       call tr_deallocate(thread_struct(i)%eq_tt  ,"eq_tt",CAT_MATELEM)
       call tr_deallocate(thread_struct(i)%eq_pp  ,"eq_pp",CAT_MATELEM) 
       call tr_deallocate(thread_struct(i)%delta_g,"delta_g",CAT_MATELEM) 
       call tr_deallocate(thread_struct(i)%delta_s,"delta_s",CAT_MATELEM)
       call tr_deallocate(thread_struct(i)%delta_t,"delta_t",CAT_MATELEM)
#ifdef COMPARE_ELEMENT_MATRIX
       call tr_deallocate(thread_struct(i)%ELM2,"ELM2",CAT_MATELEM)
       call tr_deallocate(thread_struct(i)%RHS2,"RHS2",CAT_MATELEM)
#endif
    end do
    call tr_unregister_mem(sizeof(thread_struct),"thread_struct",CAT_MATELEM)
    deallocate(thread_struct)
  end subroutine del_thread_buffers
  
end module data_structure



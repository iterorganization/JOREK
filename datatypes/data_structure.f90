!> Definitions of derived data types for grid nodes and elements, boundary nodes and elements,
!! and flux surface elements, as well as the shattered pellets
module data_structure
  use mod_parameters
  use mod_integer_types
  use tr_module
  use gauss
  use ISO_C_BINDING, ONLY : C_INT, C_DOUBLE

  
  implicit none

  type type_node                                  !< type definition of a node (i.e. a vertex)
    real*8     :: x(n_coord_tor,n_degrees,n_dim)        !< x,y,z coordinates of points and additional nodal geometry
    real*8     :: values(n_tor,n_degrees,n_var)   !< Variable values and derivatives
    real*8     :: deltas(n_tor,n_degrees,n_var)   !< Change of variable values and derivatives in last timestep
#ifdef fullmhd
    real*8     :: psi_eq(n_degrees)               !< equilibrium flux at the nodes
    real*8     :: Fprof_eq(n_degrees)             !< equilibrium profile R*B_phi at the nodes
#elif altcs
    real*8     :: psi_eq(n_degrees)               !< equilibrium flux at the nodes
#endif
    integer    :: index(n_degrees)                !< index in the main matrix
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
    real*8  :: size(n_vertex_max,n_degrees)       !< size of vectors at each vertex of the element
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
     real*8, dimension (:,:,:), allocatable:: ELM_pnn
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

     real*8, dimension(:,:,:,:) , allocatable :: eq_aux_g, eq_aux_s, eq_aux_t, eq_aux_p
     real*8, dimension(:,:,:,:) , allocatable :: delta_aux_g, delta_aux_s, delta_aux_t, delta_aux_p

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
    integer :: plasmoid_in_domain    !< Flag representing whether (post-teleportation) plasmoids are in computational domain
                                     !! this is only relevant if drift_distance /= 0
  end type type_SPI
 
  !> Sparse matrix type (generally distributed)
  type type_SP_MATRIX
    integer(kind=int_all), dimension(:), pointer   :: irn => Null()
    integer(kind=int_all), dimension(:), pointer   :: jcn => Null()
    real(kind=8), dimension(:), pointer            :: val => Null()
    integer(kind=int_all), dimension(:), pointer   :: ijA_size => Null()
    integer(kind=int_all), dimension(:,:), pointer :: ijA_index => Null()
    integer(kind=int_all), dimension(:,:), pointer :: irn_jcn => Null()
    
    real(kind=8), dimension(:), pointer          :: column_scaling => Null()    !< global column scaling, vector size of ng
    integer                                      :: indexing = 1         !< matrix indexing (1 is standart FORTRAN)
    integer(kind=int_all)                        :: ng = 0               !< matrix total rank
    integer(kind=int_all)                        :: nr = 0               !< number of local rows
    integer(kind=int_all)                        :: nc = 0               !< number of local cols
    integer(kind=int_all)                        :: nnz = 0              !< number of local nonzero entries
    integer, dimension(:), pointer               :: index_min => Null()  !< minimum node index in global range for all MPI ranks
    integer, dimension(:), pointer               :: index_max => Null()  !< maximum node index in global range for all MPI ranks
    integer                                      :: my_ind_min = 0
    integer                                      :: my_ind_max = 0
    integer                                      :: i_tor_min = 0        ! minimum toroidal Fourier number used in construction
    integer                                      :: i_tor_max = 0        ! maximum toroidal Fourier number used in construction
    integer                                      :: block_size = 1
    integer                                      :: comm = 0             !< communicator over which the matrix is distributed
    logical                                      :: scaled = .false.
    logical                                      :: row_distributed = .false.
    logical                                      :: col_distributed = .false.
    logical                                      :: reduced = .false. !< matrix is available on all comm ranks (not distribued)

  contains
    procedure :: copy_to
    procedure :: move_to
    procedure :: reset
  end type type_SP_MATRIX
  
  !> RHS vector type
  type type_RHS
    real(kind=8), dimension(:), pointer :: val => Null()
    integer(kind=int_all)               :: n                    !< vector length
  end type type_RHS  
  
  !> Preconditioner type  
  type type_PRECOND
    type(type_SP_MATRIX)                         :: mat                           !< PC matrix structure
    type(type_RHS)                               :: rhs                           !< PC rhs structure
    
    integer                                      :: n_mode_families               !< number of mode families (input)
    integer, dimension(:), pointer               :: modes_per_family => Null()    !< number of toroidal modes per mode family (input)    
    integer, dimension(:), pointer               :: ranks_per_family => Null()    !< number of MPI tasks per mode family (input)
    logical                                      :: autodistribute_modes          !< if true - use single mode par family (input)
    logical                                      :: autodistribute_ranks          !< if true - distribute MPI ranks equally between mode families (input)
    integer(kind=int_all), dimension(:), pointer :: row_index => Null()           !< Row indices of local mode family in global RHS
    real(kind=8)                                 :: row_factor                    !< Multiplying factor of current mode family in global RHS           

    integer                                      :: family_id                     !< family id (MPI private)
    integer                                      :: mode_set_n                    !< number of modes in current mode family    
    integer, dimension(:), pointer               :: mode_set => Null()            !< toroidal modes in current mode family
    integer, dimension(:,:), pointer             :: mode_families_ranks => Null() !< MPI ranks which belong to each mode family
    integer, dimension(:,:), pointer             :: mode_families_modes => Null() !< Toroidal modes which belong to each mode family
    
    integer, dimension(:), pointer               :: rank_range => Null()          !< range of MPI ranks which belong to mode families
    integer                                      :: my_id, n_cpu, comm    
    integer                                      :: my_id_n, n_cpu_n, MPI_COMM_N
    integer                                      :: my_id_master, n_masters, MPI_COMM_MASTER, MPI_COMM_TRANS, MPI_GROUP_WORLD, MPI_GROUP_MASTER
! the following variables are used in PC distribution (they are set only once to save computation time)
    integer, dimension(:,:), pointer             :: send_counts => Null()         !< number of entries sent to each other MPI ranks (PC distribution)
    integer, dimension(:,:), pointer             :: recv_counts => Null()         !< number of entries received from each other MPI ranks (PC distribution)
    integer, dimension(:,:), pointer             :: send_disp => Null()           !< send dispalcements for mpi_alltoallv (PC distribution)
    integer, dimension(:,:), pointer             :: splt_disp => Null()           !< receive displacement for split communication
    integer(kind=int_all), dimension(:,:), pointer :: recv_disp => Null()           !< receive dispalcements for mpi_alltoallv (PC distribution)
    integer(kind=int_all), dimension(:), pointer :: istart => Null()              !< start-index for split communication
    integer(kind=int_all), dimension(:), pointer :: ifinish => Null()             !< end-index for split communication
    integer                                      :: nsplit                        !< number of communication splits
    integer(kind=int_all), dimension(:), pointer :: n_per_rank => Null()          !< min number of rows/cols per MPI rank for each family

    logical                                      :: initialized = .false.
    logical                                      :: structured = .false.          !< flag indicating the allocation of PC matrix structure
    integer(kind=int_all)                        :: n_glob                        !< global number of unknowns
    
#ifdef DIRECT_CONSTRUCTION
    integer, dimension(:), pointer               :: local_elms => null()
    integer                                      :: n_local_elms
#endif    
  end type type_PRECOND

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
          call tr_allocate(thread_struct(i)%ELM_p, 1,n_plane,1,n_vertex_max*n_var*n_degrees,1,n_vertex_max*n_var*n_degrees,"ELM_p",CAT_MATELEM)
          call tr_allocate(thread_struct(i)%ELM_n, 1,n_plane,1,n_vertex_max*n_var*n_degrees,1,n_vertex_max*n_var*n_degrees,"ELM_n",CAT_MATELEM)
          call tr_allocate(thread_struct(i)%ELM_k, 1,n_plane,1,n_vertex_max*n_var*n_degrees,1,n_vertex_max*n_var*n_degrees,"ELM_k",CAT_MATELEM)
          call tr_allocate(thread_struct(i)%ELM_kn,1,n_plane,1,n_vertex_max*n_var*n_degrees,1,n_vertex_max*n_var*n_degrees,"ELM_kn",CAT_MATELEM)
          call tr_allocate(thread_struct(i)%ELM_pnn,1,n_plane,1,n_vertex_max*n_var*n_degrees,1,n_vertex_max*n_var*n_degrees,"ELM_pnn",CAT_MATELEM)
          call tr_allocate(thread_struct(i)%RHS_p, 1,n_plane,1,n_vertex_max*n_var*n_degrees,"RHS_p",CAT_MATELEM)                                     
          call tr_allocate(thread_struct(i)%RHS_k, 1,n_plane,1,n_vertex_max*n_var*n_degrees,"RHS_k",CAT_MATELEM)                                     
          call tr_allocate(thread_struct(i)%ELM,   1,n_tor*n_vertex_max*n_degrees*n_var,1,n_tor*n_vertex_max*n_degrees*n_var,"ELM",CAT_MATELEM)       
          call tr_allocate(thread_struct(i)%RHS,   1,n_tor*n_vertex_max*n_degrees*n_var,"RHS",CAT_MATELEM)
          call tr_allocate(thread_struct(i)%synch_buff, 1,n_tor*n_var*n_tor*n_var,"synch_buff",CAT_MATELEM)
          thread_struct(i)%ELM_p   = 0.d0
          thread_struct(i)%ELM_n   = 0.d0
          thread_struct(i)%ELM_k   = 0.d0
          thread_struct(i)%ELM_kn  = 0.d0
          thread_struct(i)%ELM_pnn = 0.d0
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
          call tr_allocate(thread_struct(i)%eq_aux_g,1,n_plane,1,n_var,1,n_gauss,1,n_gauss,"eq_aux_g",CAT_MATELEM)
          call tr_allocate(thread_struct(i)%eq_aux_s,1,n_plane,1,n_var,1,n_gauss,1,n_gauss,"eq_aux_s",CAT_MATELEM)
          call tr_allocate(thread_struct(i)%eq_aux_t,1,n_plane,1,n_var,1,n_gauss,1,n_gauss,"eq_aux_t",CAT_MATELEM)
          call tr_allocate(thread_struct(i)%eq_aux_p,1,n_plane,1,n_var,1,n_gauss,1,n_gauss,"eq_aux_p",CAT_MATELEM)
          call tr_allocate(thread_struct(i)%delta_aux_g,1,n_plane,1,n_var,1,n_gauss,1,n_gauss,"delta_aux_g",CAT_MATELEM)
          call tr_allocate(thread_struct(i)%delta_aux_s,1,n_plane,1,n_var,1,n_gauss,1,n_gauss,"delta_aux_s",CAT_MATELEM)
          call tr_allocate(thread_struct(i)%delta_aux_t,1,n_plane,1,n_var,1,n_gauss,1,n_gauss,"delta_aux_t",CAT_MATELEM)
          call tr_allocate(thread_struct(i)%delta_aux_p,1,n_plane,1,n_var,1,n_gauss,1,n_gauss,"delta_aux_p",CAT_MATELEM)
          thread_struct(i)%eq_g     = 0.d0
          thread_struct(i)%eq_s     = 0.d0
          thread_struct(i)%eq_t     = 0.d0
          thread_struct(i)%eq_p     = 0.d0
          thread_struct(i)%eq_ss    = 0.d0
          thread_struct(i)%eq_st    = 0.d0
          thread_struct(i)%eq_tt    = 0.d0
          thread_struct(i)%eq_pp    = 0.d0
          thread_struct(i)%delta_g  = 0.d0
          thread_struct(i)%delta_s  = 0.d0
          thread_struct(i)%delta_t  = 0.d0
          thread_struct(i)%eq_aux_g = 0.d0
          thread_struct(i)%eq_aux_s = 0.d0
          thread_struct(i)%eq_aux_t = 0.d0
          thread_struct(i)%eq_aux_p = 0.d0
          thread_struct(i)%delta_aux_g = 0.d0
          thread_struct(i)%delta_aux_s = 0.d0
          thread_struct(i)%delta_aux_t = 0.d0
          thread_struct(i)%delta_aux_p = 0.d0
#ifdef COMPARE_ELEMENT_MATRIX
          call tr_allocate(thread_struct(i)%ELM2,  1,n_tor*n_vertex_max*n_degrees*n_var,1,n_tor*n_vertex_max*n_degrees*n_var,"ELM2",CAT_MATELEM)
          call tr_allocate(thread_struct(i)%RHS2,  1,n_tor*n_vertex_max*n_degrees*n_var,"RHS2",CAT_MATELEM)
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
       call tr_deallocate(thread_struct(i)%ELM_pnn,"ELM_pnn",CAT_MATELEM)
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
       call tr_deallocate(thread_struct(i)%eq_aux_g,"eq_aux_g",CAT_MATELEM)
       call tr_deallocate(thread_struct(i)%eq_aux_s,"eq_aux_s",CAT_MATELEM)
       call tr_deallocate(thread_struct(i)%eq_aux_t,"eq_aux_t",CAT_MATELEM)
       call tr_deallocate(thread_struct(i)%eq_aux_p,"eq_aux_p",CAT_MATELEM)
       call tr_deallocate(thread_struct(i)%delta_aux_g,"eq_aux_g",CAT_MATELEM)
       call tr_deallocate(thread_struct(i)%delta_aux_s,"eq_aux_s",CAT_MATELEM)
       call tr_deallocate(thread_struct(i)%delta_aux_t,"eq_aux_t",CAT_MATELEM)
       call tr_deallocate(thread_struct(i)%delta_aux_p,"eq_aux_p",CAT_MATELEM)
#ifdef COMPARE_ELEMENT_MATRIX
       call tr_deallocate(thread_struct(i)%ELM2,"ELM2",CAT_MATELEM)
       call tr_deallocate(thread_struct(i)%RHS2,"RHS2",CAT_MATELEM)
#endif
    end do
    call tr_unregister_mem(sizeof(thread_struct),"thread_struct",CAT_MATELEM)
    deallocate(thread_struct)
  end subroutine del_thread_buffers

! move one matrix structure into another
  subroutine move_to(self, mat_a, with_data)
    class(type_SP_MATRIX), intent(inout)    :: self
    class(type_SP_MATRIX), intent(inout) :: mat_a
    logical                              :: with_data

    if (with_data) then
      mat_a%irn            => self%irn; self%irn => null()
      mat_a%jcn            => self%jcn; self%jcn => null()
      mat_a%val            => self%val; self%val => null()

      mat_a%ijA_size       => self%ijA_size
      mat_a%ijA_index      => self%ijA_index
      mat_a%irn_jcn        => self%irn_jcn

      mat_a%column_scaling => self%column_scaling
      mat_a%index_min      => self%index_min
      mat_a%index_max      => self%index_max
    endif

    mat_a%indexing        = self%indexing
    mat_a%ng              = self%ng
    mat_a%nr              = self%nr
    mat_a%nc              = self%nc
    mat_a%nnz             = self%nnz
    mat_a%my_ind_min      = self%my_ind_min
    mat_a%my_ind_max      = self%my_ind_max
    mat_a%i_tor_min       = self%i_tor_min
    mat_a%i_tor_max       = self%i_tor_max
    mat_a%block_size      = self%block_size
    mat_a%comm            = self%comm
    mat_a%scaled          = self%scaled
    mat_a%row_distributed = self%row_distributed
    mat_a%col_distributed = self%col_distributed
    mat_a%reduced         = self%reduced

    return
  end subroutine

! copy one matrix structure into another
  subroutine copy_to(self, mat_a)
    class(type_SP_MATRIX), intent(inout)    :: self
    class(type_SP_MATRIX), intent(inout)    :: mat_a

    call mat_a%reset()

    call self%move_to(mat_a, with_data=.false.)

    if (associated(self%irn)) then
      allocate(mat_a%irn(mat_a%nnz))
      mat_a%irn(1:mat_a%nnz) = self%irn(1:self%nnz)
    endif
    if (associated(self%jcn)) then
      allocate(mat_a%jcn(mat_a%nnz))
      mat_a%jcn(1:mat_a%nnz) = self%jcn(1:self%nnz)
    endif
    if (associated(self%val)) then
      allocate(mat_a%val(mat_a%nnz))
      mat_a%val(1:mat_a%nnz) = self%val(1:self%nnz)
    endif
    if (associated(self%column_scaling)) then
      allocate(mat_a%column_scaling(mat_a%ng))
      mat_a%column_scaling(1:mat_a%ng) = self%column_scaling(1:self%ng)
    endif
    !mat_a%ijA_size       => self%ijA_size
    !mat_a%ijA_index      => self%ijA_index
    !mat_a%irn_jcn        => self%irn_jcn
    !mat_a%index_min      => self%index_min
    !mat_a%index_max      => self%index_max

    return
  end subroutine

  subroutine reset(self)
    class(type_SP_MATRIX), intent(inout)    :: self

    if (associated(self%irn)) then
      deallocate(self%irn); self%irn => Null()
    endif
    if (associated(self%jcn)) then
      deallocate(self%jcn); self%jcn => Null()
    endif
    if (associated(self%val)) then
      deallocate(self%val); self%val => Null()
    endif
    if (associated(self%ijA_size)) then
      deallocate(self%ijA_size); self%ijA_size => Null()
    endif
    if (associated(self%ijA_index)) then
      deallocate(self%ijA_index); self%ijA_index => Null()
    endif
    if (associated(self%irn_jcn)) then
      deallocate(self%irn_jcn); self%irn_jcn => Null()
    endif
    if (associated(self%column_scaling)) then
      deallocate(self%column_scaling); self%column_scaling => Null()
    endif
    if (associated(self%index_min)) then
      deallocate(self%index_min); self%index_min => Null()
    endif
    if (associated(self%index_max)) then
      deallocate(self%index_max); self%index_max => Null()
    endif

    self%indexing = 1
    self%ng = 0
    self%nr = 0
    self%nc = 0
    self%nnz = 0
    self%my_ind_min = 0
    self%my_ind_max = 0
    self%i_tor_min = 0
    self%i_tor_max = 0
    self%block_size = 1
    self%comm = 0
    self%scaled = .false.
    self%row_distributed = .false.
    self%col_distributed = .false.
    self%reduced = .false.
  end subroutine

end module data_structure



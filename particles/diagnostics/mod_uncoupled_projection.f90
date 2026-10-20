module mod_uncoupled_projection

implicit none
private

public assemble_system, project_only
public assemble_projection_matrix, assemble_projection_matrix_n0 !< for testing purpouse

contains
  subroutine assemble_system(node_list, element_list, n_tor_local, i_tor_local,  &
                                  this_mpi_comm_world, this_mpi_comm_n, this_mpi_comm_master,  &
                                  a_mat, area, volume,                                         &
                                  filter, filter_hyper, filter_parallel,                       &
                                  filter_n0, filter_hyper_n0, filter_parallel_n0,              &
                                  integral_weights, do_zonal,                                  &
                                  apply_dirichlet_condition_in,                                &
                                  system_size, n_dof)

    use data_structure, only: type_node_list, type_element_list
    use mod_sparse_matrix, only: type_SP_MATRIX

    implicit none

    type (type_node_list), intent(in)    :: node_list 
    type (type_element_list), intent(in) :: element_list
    integer, intent(in)                  :: n_tor_local, i_tor_local
    integer, intent(in)                  :: this_mpi_comm_world, this_mpi_comm_n, this_mpi_comm_master
    real*8, intent(in)                   :: filter,    filter_hyper,    filter_parallel    !< normal, hyper and parallel smoothing
    real*8, intent(in)                   :: filter_n0, filter_hyper_n0, filter_parallel_n0 !< for n=0
    logical, intent(in), optional        :: do_zonal
    logical, intent(in), optional        :: apply_dirichlet_condition_in
    real*8, intent(out), dimension(:), allocatable :: integral_weights !< these will be filled with the volume of each basis function
    real*8, intent(out)                  :: area, volume

    integer,                 intent(out)   :: system_size
    integer,                 intent(out)   :: n_dof

    type(type_SP_MATRIX),    intent(inout) :: a_mat

    if (n_tor_local .eq. 1) then

      call assemble_projection_matrix_n0(node_list, element_list, n_tor_local, i_tor_local,  & 
                                this_mpi_comm_world, this_mpi_comm_n, this_mpi_comm_master,  &
                                a_mat, area, volume,                                         &
                                filter_n0, filter_hyper_n0, filter_parallel_n0,              &
                                integral_weights=integral_weights, do_zonal=do_zonal,        &
                                apply_dirichlet_condition_in=apply_dirichlet_condition_in) 
      
      system_size = a_mat%ng
      n_dof = system_size / 2.d0
    
    else

      call assemble_projection_matrix(node_list, element_list, n_tor_local, i_tor_local,  &
                            this_mpi_comm_world, this_mpi_comm_n, this_mpi_comm_master,   &
                            a_mat, filter, filter_hyper, filter_parallel,                 &
                            apply_dirichlet_condition_in=apply_dirichlet_condition_in)
      
      system_size = a_mat%ng
      n_dof = system_size / n_tor_local
    
    endif

  end subroutine assemble_system


  !> Gather all of the rhs-es into a single matrix and feed it to mumps, and then
  !> broadcast the result
  subroutine project_only(this, sim)

    use mod_projection_type, only: t_projection

    use mod_solve_sparse_projection, only: solve_sparse_projection_system
    use mod_particle_sim, only: particle_sim
    use mod_parameters, only: n_tor, n_vertex_max, n_degrees
    use phys_module, only: n_aux_var
    use data_structure, only: init_node

    use mpi_mod
    use mod_event

    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
    use, intrinsic :: ieee_exceptions

    !$ use omp_lib
    class(t_projection), intent(inout) :: this
    type(particle_sim), intent(inout)    :: sim
    integer :: my_id, ierr
    integer :: in, i_elm, i, j, k, i_var, i_start, i_tor
    integer :: index_large_i, inode, index
    real*8  :: t0, t1, ostart, oend, mmm(3), mmm2(3)
    real*8,  allocatable :: my_rhs(:,:), y_tmp(:)
    integer, allocatable :: recv_counts(:), recv_disp(:)
    integer :: n_rhs, n_rhs_f, i_rhs, n_tor_local, i_tor_local, n_loc_n
    integer :: in_local, in_global, index_n, id_master_in_world, offset, i_glob1, i_glob2
    logical :: halt(size(IEEE_USUAL,1)), found_nan

    call MPI_COMM_RANK(MPI_COMM_WORLD, my_id, ierr)
    call cpu_time(t0)
    !$ ostart = omp_get_wtime()

    ! Safety checks
    if (.not. allocated(sim%groups)) return
    
    if (.not. allocated(this%rhs)) then
      n_rhs = 0
    else
      n_rhs = size(this%rhs,5)
    end if

    if (.not. allocated(this%rhs_f)) then
      n_rhs_f = 0
    else 
      n_rhs_f = size(this%rhs_f,5)
    end if

    ! reinitialise the storage node_list to ensure all projections fit
    do i=1, this%node_list%n_nodes
      call init_node(this%node_list%node(i), n_rhs_f+n_rhs)
    enddo

    n_tor_local = this%n_tor_local
    i_tor_local = this%i_tor_local

    this%rhs_vec%nrhs = (n_rhs + n_rhs_f)
    this%rhs_vec%n = this%system_size

    if (2*this%n_dof .ne. this%rhs_vec%n) then
      write(*,*) 'FATAL : 2*this%n_dof .ne. this%rhs_vec%n'
      write(*,*) n_tor_local*this%n_dof,  this%rhs_vec%n
    endif

    if (associated(this%rhs_vec%val)) then
      deallocate(this%rhs_vec%val); this%rhs_vec%val => Null()
    endif
    allocate(this%rhs_vec%val(this%rhs_vec%n*this%rhs_vec%nrhs))
    this%rhs_vec%val = 0.d0

    allocate(my_rhs(this%rhs_vec%n * this%rhs_vec%nrhs, (n_tor+1)/2))
    
    my_rhs = 0.d0

    do i_rhs=1,n_rhs

      i_start =  this%rhs_vec%n * (i_rhs-1)
          
      do i_elm=1,this%element_list%n_elements

        do i=1,n_vertex_max

          inode = this%element_list%element(i_elm)%vertex(i)

          do j=1,n_degrees

            index_large_i = 2 * (this%node_list%node(inode)%index(j)-1) + 1 + i_start ! base index in the main matrix + rhs index

            my_rhs(index_large_i,1)   = my_rhs(index_large_i,1)   + this%rhs(j, i, i_elm, 1, i_rhs) !/ this%rhs_gather_time
            my_rhs(index_large_i+1,1) = 0.d0

            do in=2, n_tor, 2

              index_n = in / 2 + 1

              my_rhs(index_large_i,  index_n) = my_rhs(index_large_i,  index_n) + this%rhs(j, i, i_elm, in,   i_rhs) !/ this%rhs_gather_time
              my_rhs(index_large_i+1,index_n) = my_rhs(index_large_i+1,index_n) + this%rhs(j, i, i_elm, in+1, i_rhs) !/ this%rhs_gather_time

            enddo
          enddo

        enddo
      enddo
    enddo
    
    this%rhs_gather_time = 0.d0

    ! Reset the RHS for the next step
    if (allocated(this%rhs)) this%rhs = 0.d0

    do i_rhs=1,n_rhs_f
      ! Fill projection function part
      i_start =  this%rhs_vec%n * (n_rhs + i_rhs - 1)
      
      do i_elm=1,this%element_list%n_elements
          
        do i=1,n_vertex_max
        
          inode = this%element_list%element(i_elm)%vertex(i)
            
          do j=1,n_degrees

            index_large_i = 2*(this%node_list%node(inode)%index(j)-1) + 1 + i_start ! base index in the main matrix + rhs index

            my_rhs(index_large_i,1)   = my_rhs(index_large_i,1)   + this%rhs_f(j, i, i_elm, 1, i_rhs)
            my_rhs(index_large_i+1,1) = 0.

            do in=2, n_tor, 2

              index_n = in / 2 + 1
          
              my_rhs(index_large_i,  index_n) = my_rhs(index_large_i,  index_n) + this%rhs_f(j, i, i_elm, in,   i_rhs)
              my_rhs(index_large_i+1,index_n) = my_rhs(index_large_i+1,index_n) + this%rhs_f(j, i, i_elm, in+1, i_rhs)
              
            enddo
      
          enddo
        enddo

      enddo
    enddo

    ! Gather the RHS's to the root process
    ! results are gathered on the my_id_n=0 nodes
    call MPI_Reduce(my_rhs(:,1),this%rhs_vec%val,this%rhs_vec%n*this%rhs_vec%nrhs, MPI_REAL8, MPI_SUM, 0, this%mpi_comm_world, ierr)

    do in=2, n_tor, 2
      id_master_in_world = in/2 * this%m_mpi
      index_n = in/2 + 1
      call MPI_Reduce(my_rhs(:,index_n),this%rhs_vec%val,this%rhs_vec%n*this%rhs_vec%nrhs, MPI_REAL8, MPI_SUM, id_master_in_world, this%mpi_comm_world, ierr)
    enddo

    if ((this%my_id .eq. 0) .and. (this%do_zonal)) then   

      write(*,'(A,3e14.6)') 'check project :',maxval(this%rhs_vec%val(:)), maxval(this%scaling_integral_weights*this%integral_weights), &
                                              maxval(this%rhs_vec%val(:) - this%scaling_integral_weights * this%integral_weights)

      this%rhs_vec%val(1:this%rhs_vec%n) = this%rhs_vec%val(1:this%rhs_vec%n) - this%scaling_integral_weights * this%integral_weights(:)
    endif


    call ieee_get_halting_mode(IEEE_USUAL, halt)
    call ieee_set_halting_mode(IEEE_USUAL, [.false., .false., .false.])

    ! Compute the solution of Ax=B (B = RHSes)
    call solve_sparse_projection_system(this%a_mat, this%rhs_vec, this%solver)

    call ieee_set_halting_mode(IEEE_USUAL, halt)
    
    ! collect the solution of all the toroidal harmonics (ntor+1)/2 to my_id=0
    if (this%my_id_n .eq. 0) then

      n_loc_n = this%rhs_vec%n * this%rhs_vec%nrhs
    
      allocate(y_tmp(n_loc_n*(n_tor+1)))            ! allocate only on my_id=0???

      allocate(recv_counts(this%n_mpi/this%m_mpi))
      allocate(recv_disp(this%n_mpi/this%m_mpi))
      
      y_tmp = 0.d0

      recv_counts(1) = n_loc_n
      do i=2,(n_tor+1)/2
        recv_counts(i) = n_loc_n
      enddo

      recv_disp(1) = 0
      do i=2,(n_tor+1)/2
        recv_disp(i) = recv_disp(i-1) + recv_counts(i-1)
      enddo

      call mpi_gatherv(this%rhs_vec%val, this%rhs_vec%n * this%rhs_vec%nrhs, MPI_DOUBLE_PRECISION, &
                      y_tmp, recv_counts, recv_disp, MPI_DOUBLE_PRECISION, 0, this%mpi_comm_master,ierr)

    endif

    call MPI_BARRIER(this%mpi_comm_world, ierr)

    ! Write the solution to the node_list
    if (this%my_id .eq. 0) then

      do i_var=1,n_rhs+n_rhs_f
    
        found_nan = .false.
        
        do i=1,this%node_list%n_nodes

          do k=1,n_degrees
        
            index = this%node_list%node(i)%index(k)

            i_glob1 = 2*(index-1) + 1 + 2*this%n_dof*(i_var-1)

            if (this%do_zonal) then
              this%node_list%node(i)%deltas(1,k,i_var) = y_tmp(i_glob1) + y_tmp(i_glob1+1) - this%node_list%node(i)%values(1,k, i_var)
              this%node_list%node(i)%deltas(1,k,2)     = y_tmp(i_glob1+1)                  - this%node_list%node(i)%values(1,k,     2)
              
              this%node_list%node(i)%values(1,k,i_var) = y_tmp(i_glob1) + y_tmp(i_glob1+1)
              this%node_list%node(i)%values(1,k,2)     = y_tmp(i_glob1+1)
            else
              this%node_list%node(i)%deltas(1,k,i_var) = y_tmp(i_glob1) - this%node_list%node(i)%values(1,   k, i_var)
              this%node_list%node(i)%values(1,k,i_var) = y_tmp(i_glob1)
            endif

            offset = 2*this%n_dof * this%rhs_vec%nrhs
            
            do i_tor=2,n_tor,2

              i_glob2 = 2*(index-1) + 1 + offset + 2*this%n_dof*(i_var-1) + (i_tor-2)*this%n_dof * this%rhs_vec%nrhs

              this%node_list%node(i)%deltas(i_tor,  k,i_var) = y_tmp(i_glob2)   - this%node_list%node(i)%values(i_tor,   k, i_var)
              this%node_list%node(i)%deltas(i_tor+1,k,i_var) = y_tmp(i_glob2+1) - this%node_list%node(i)%values(i_tor+1, k, i_var)

              this%node_list%node(i)%values(i_tor,  k,i_var) = y_tmp(i_glob2)   
              this%node_list%node(i)%values(i_tor+1,k,i_var) = y_tmp(i_glob2+1)
            
            end do
          
          enddo    ! order
          
          ! Check for NaNs in the projection
          if (any(ieee_is_nan(this%node_list%node(i)%values(:,:,i_var))) .or. any(ieee_is_nan(this%node_list%node(i)%deltas(:,:,i_var)))) then
            found_nan = .true.
          end if
        
        enddo      ! nodes

        if (found_nan) then
          write(*,*) "Found NaNs in projection number ", i_var
        end if

      enddo
    endif

    call MPI_BARRIER(this%mpi_comm_world, ierr)

    call broadcast_nodes(this%my_id, this%node_list)

#ifdef DEBUG
      call cpu_time(t1)
      !$ oend = omp_get_wtime()
      !$ mmm = mpi_minmeanmax(t1-t0)
      !$ mmm2 = mpi_minmeanmax(oend-ostart)
      if (this%my_id .eq. 0) then
        write(*,"(A,3g12.5)") "projection cpu time", mmm
        !$ write(*,"(A,3g12.5)") "projection wall time", mmm2
      end if
#endif

  end subroutine project_only

  !> Project particles by weight onto the elements
  !> Saves output in node_list%values(1)
  !> The projection is done by solving
  !> $$\int [p(x) u(x) + \lambda p'(x) u'(x)] dV =
  !> \int \sum \delta(x_i-x) w_i u(x) dV $$
  !> where p(x) is in the bernstein representation, $u(x)$ are the test functions
  !> composed of two basis functions and $w_i$ is the particle weight.
  !> A smoothing factor lambda is included
  !> x is a vector (R,Z,phi) and dV is r dr dphi
  !> divide by 1 or 2pi on both sides (LHS gets 2pi for n=0 mode, 1pi for other modes)
  !>
  !> See also [project_particles] and [mod_elt_matrix] for reference of the integration method
  subroutine assemble_projection_matrix(node_list, element_list, n_tor_local, i_tor_local,           &
                              this_mpi_comm_world, this_mpi_comm_n, this_mpi_comm_master,  &
                              a_mat, filter, filter_hyper, filter_parallel,            &
                              skip_factorisation, apply_dirichlet_condition_in)
  use phys_module, only : F0, TWOPI, mode, fix_axis_nodes
  use data_structure
  use basis_at_gaussian
  use mod_basisfunctions
  use mpi_mod
  use, intrinsic :: ieee_exceptions
  implicit none

  type (type_node_list), intent(in)    :: node_list !< A copy of the node list which will be used to save variables
  type (type_element_list), intent(in) :: element_list
  integer, intent(in)                  :: n_tor_local, i_tor_local
  integer, intent(in)                  :: this_mpi_comm_world, this_mpi_comm_n, this_mpi_comm_master
  real*8, intent(in)                   :: filter
  real*8, intent(in)                   :: filter_hyper
  real*8, intent(in)                   :: filter_parallel
  logical, intent(in), optional        :: skip_factorisation
  logical,intent(in),optional          :: apply_dirichlet_condition_in

  type (type_element)      :: element
  type (type_node)         :: nodes(n_vertex_max)

  real*8, allocatable      :: ELM(:,:)
  real*8     :: wgauss2(n_gauss)
  real*8, dimension(n_gauss,n_gauss) :: x_g,   x_s,   x_t,   x_ss,   x_tt,   x_st
  real*8, dimension(n_gauss,n_gauss) :: y_g,   y_s,   y_t,   y_ss,   y_tt,   y_st
  real*8, dimension(n_gauss,n_gauss) :: psi_g, psi_s, psi_t, psi_ss, psi_tt, psi_st
  real*8     :: v, v_s, v_t, v_ss, v_st, v_tt, v_x, v_y, v_xx, v_yy, v_p
  real*8     :: p, p_s, p_t, p_ss, p_st, p_tt, p_x, p_y, p_xx, p_yy, p_p
  real*8     :: wst, area, volume, xjac, xjac_x, xjac_y, psi_x, psi_y
  real*8     :: Bgrad_p, Bgrad_v_star, BB2
  integer    :: i, j, k, l, m, in, im, ilarge, index_large_i, index_large_k, inode, knode
  integer    :: nz_AA, n_AA, nz_bnd, i_elm, index_ij, index_kl, im_index, in_index, index1
  integer    :: ms, mt, mp, my_id, my_id_n, my_id_master, ierr, MPI_COMM_MUMPS
  logical    :: apply_dirichlet_condition
  logical    :: do_facto

  type (type_SP_MATRIX) :: a_mat  !< Projection matrix

  ! We need a separate communicator to be able to run multiple MUMPSes
  call MPI_Comm_dup(this_mpi_comm_n, MPI_COMM_MUMPS, ierr)
  a_mat%comm = MPI_COMM_MUMPS

  call MPI_COMM_RANK(this_mpi_comm_world,  my_id  , ierr)
  call MPI_COMM_RANK(a_mat%comm,           my_id_n, ierr)

  nz_AA = 4 * element_list%n_elements * (n_vertex_max * n_degrees)**2
  n_AA  = 2 * maxval(node_list%node(1:node_list%n_nodes)%index(4))

  apply_dirichlet_condition = .true.
  if(present(apply_dirichlet_condition_in)) apply_dirichlet_condition = apply_dirichlet_condition_in

  nz_bnd = 0
  if (apply_dirichlet_condition) then
    do i=1,node_list%n_nodes
      if (node_list%node(i)%boundary .eq. 1) nz_bnd = nz_bnd + 4
      if (node_list%node(i)%boundary .eq. 2) nz_bnd = nz_bnd + 4
      if (node_list%node(i)%boundary .eq. 3) nz_bnd = nz_bnd + 8
      if (node_list%node(i)%boundary .eq. 4) nz_bnd = nz_bnd + 4
      if (node_list%node(i)%boundary .eq. 5) nz_bnd = nz_bnd + 4
      if (node_list%node(i)%boundary .eq. 9) nz_bnd = nz_bnd + 8
    enddo
  endif

  ! Only perform the construction of the matrix on the host
  if (my_id_n .eq. 0) then

    allocate(ELM(2*n_vertex_max*n_degrees,2*n_vertex_max*n_degrees))

    ! Allocate space for elements

    allocate(  a_mat%val(nz_AA+nz_bnd),    a_mat%irn(nz_AA+nz_bnd),    a_mat%jcn(nz_AA+nz_bnd))

    a_mat%irn = 0
    a_mat%jcn = 0
    a_mat%val = 0.d0

    ! Copy wgauss into wgauss2 to get around gfortran not recognizing it as a shared
    ! thing https://groups.google.com/forum/#!topic/comp.lang.fortran/VKhoAm8m9KE
    wgauss2 = wgauss

    write(*,*) '**************************************************'
    write(*,'(A,i2,A)') ' * constructing particle projection matrix (n=',mode(i_tor_local),') *'
    write(*,*) '**************************************************'
    write(*,*) ' n_AA = ',n_AA, n_tor_local
    write(*,'(2i3,A,3e12.4)') my_id,my_id_n,'  filters       ',filter, filter_hyper, filter_parallel

    if (apply_dirichlet_condition) write(*,*) 'applying Dirichlet conditions'

    !$omp parallel do default(none) &
    !$omp shared(element_list, node_list, n_tor_local, i_tor_local,                       &
    !$omp        H, H_s, H_t, H_ss, H_st, H_tt, Hz, Hz_p, a_mat, wgauss2,                 &
    !$omp        filter, filter_hyper, filter_parallel, F0, fix_axis_nodes, my_id_master) &
    !$omp private(ELM, i_elm, element, i, j, k, l, ms, mt, in, im, mp,           &
    !$omp         x_g, x_s, x_t, x_ss, x_st, x_tt,                                      &
    !$omp         y_g, y_s, y_t, y_ss, y_st, y_tt,                                      &
    !$omp         psi_g, psi_s, psi_t, psi_ss, psi_tt, psi_st,                          &
    !$omp         v, v_s, v_t, v_ss, v_st, v_tt, v_x, v_y, v_xx, v_yy, v_p,             &
    !$omp         p, p_s, p_t, p_ss, p_st, p_tt, p_x, p_y, p_xx, p_yy, p_p,             &
    !$omp         wst, xjac, xjac_x, xjac_y, psi_x, psi_y, BB2, Bgrad_p, Bgrad_v_star,  &
    !$omp         index_ij, index_kl, ilarge, in_index, im_index,                       &
    !$omp         inode, index_large_i, knode, index_large_k)                           &
    !$omp firstprivate(nodes)                                                           & 
    !$omp schedule(static) 
    do i_elm=1,element_list%n_elements
      
      ELM = 0.d0

      element = element_list%element(i_elm)
      do m=1,n_vertex_max
        call make_deep_copy_node(node_list%node(element%vertex(m)), nodes(m))
      enddo

      ! Set up gauss points in this element
      x_g = 0.d0;   x_s = 0.d0;   x_t = 0.d0;   x_ss = 0.d0;   x_st = 0.d0;   x_tt = 0.d0
      y_g = 0.d0;   y_s = 0.d0;   y_t = 0.d0;   y_ss = 0.d0;   y_st = 0.d0;   y_tt = 0.d0
      psi_g = 0.d0; psi_s = 0.d0; psi_t = 0.d0; psi_ss = 0.d0; psi_st = 0.d0; psi_tt = 0.d0

      do i=1,n_vertex_max
        do j=1,n_degrees
          do ms=1, n_gauss
            do mt=1, n_gauss
              x_g(ms,mt)  = x_g(ms,mt)  + nodes(i)%x(1,j,1) * element%size(i,j) * H(i,j,ms,mt)
              x_s(ms,mt)  = x_s(ms,mt)  + nodes(i)%x(1,j,1) * element%size(i,j) * H_s(i,j,ms,mt)
              x_t(ms,mt)  = x_t(ms,mt)  + nodes(i)%x(1,j,1) * element%size(i,j) * H_t(i,j,ms,mt)

              x_ss(ms,mt) = x_ss(ms,mt) + nodes(i)%x(1,j,1) * element%size(i,j) * H_ss(i,j,ms,mt)
              x_st(ms,mt) = x_st(ms,mt) + nodes(i)%x(1,j,1) * element%size(i,j) * H_st(i,j,ms,mt)
              x_tt(ms,mt) = x_tt(ms,mt) + nodes(i)%x(1,j,1) * element%size(i,j) * H_tt(i,j,ms,mt)

              y_g(ms,mt)  = y_g(ms,mt)  + nodes(i)%x(1,j,2) * element%size(i,j) * H(i,j,ms,mt)
              y_s(ms,mt)  = y_s(ms,mt)  + nodes(i)%x(1,j,2) * element%size(i,j) * H_s(i,j,ms,mt)
              y_t(ms,mt)  = y_t(ms,mt)  + nodes(i)%x(1,j,2) * element%size(i,j) * H_t(i,j,ms,mt)

              y_ss(ms,mt) = y_ss(ms,mt) + nodes(i)%x(1,j,2) * element%size(i,j) * H_ss(i,j,ms,mt)
              y_st(ms,mt) = y_st(ms,mt) + nodes(i)%x(1,j,2) * element%size(i,j) * H_st(i,j,ms,mt)
              y_tt(ms,mt) = y_tt(ms,mt) + nodes(i)%x(1,j,2) * element%size(i,j) * H_tt(i,j,ms,mt)

              psi_g(ms,mt)  = psi_g(ms,mt)  + nodes(i)%values(1,j,1) * element%size(i,j) * H(i,j,ms,mt)
              psi_s(ms,mt)  = psi_s(ms,mt)  + nodes(i)%values(1,j,1) * element%size(i,j) * H_s(i,j,ms,mt)
              psi_t(ms,mt)  = psi_t(ms,mt)  + nodes(i)%values(1,j,1) * element%size(i,j) * H_t(i,j,ms,mt)

              psi_ss(ms,mt) = psi_ss(ms,mt) + nodes(i)%values(1,j,1) * element%size(i,j) * H_ss(i,j,ms,mt)
              psi_st(ms,mt) = psi_st(ms,mt) + nodes(i)%values(1,j,1) * element%size(i,j) * H_st(i,j,ms,mt)
              psi_tt(ms,mt) = psi_tt(ms,mt) + nodes(i)%values(1,j,1) * element%size(i,j) * H_tt(i,j,ms,mt)
            enddo
          enddo
        enddo
      enddo

      do ms=1, n_gauss
        do mt=1, n_gauss

          wst = wgauss2(ms)*wgauss2(mt)
          xjac =  x_s(ms,mt)*y_t(ms,mt) - x_t(ms,mt)*y_s(ms,mt)

          xjac_x  = (x_ss(ms,mt)*y_t(ms,mt)**2 - y_ss(ms,mt)*x_t(ms,mt)*y_t(ms,mt) - 2.d0*x_st(ms,mt)*y_s(ms,mt)*y_t(ms,mt)   &
                  + y_st(ms,mt)*(x_s(ms,mt)*y_t(ms,mt) + x_t(ms,mt)*y_s(ms,mt))                                               &
                  + x_tt(ms,mt)*y_s(ms,mt)**2 - y_tt(ms,mt)*x_s(ms,mt)*y_s(ms,mt)) / xjac

          xjac_y  = (y_tt(ms,mt)*x_s(ms,mt)**2 - x_tt(ms,mt)*y_s(ms,mt)*x_s(ms,mt) - 2.d0*y_st(ms,mt)*x_t(ms,mt)*x_s(ms,mt)   &
                  + x_st(ms,mt)*(y_t(ms,mt)*x_s(ms,mt) + y_s(ms,mt)*x_t(ms,mt))                                               &
                  + y_ss(ms,mt)*x_t(ms,mt)**2 - x_ss(ms,mt)*y_t(ms,mt)*x_t(ms,mt)) / xjac

          psi_x = (  y_t(ms,mt) * psi_s(ms,mt) - y_s(ms,mt) * psi_t(ms,mt)) / xjac
          psi_y = (- x_t(ms,mt) * psi_s(ms,mt) + x_s(ms,mt) * psi_t(ms,mt)) / xjac

          BB2 = 1.d0
          if (filter_parallel .gt. 0.d0) BB2 = (F0*F0 + psi_x * psi_x + psi_y * psi_y )/x_g(ms,mt)**2

          do mp = 1, n_plane

            do i=1,n_vertex_max
              do j=1,n_degrees

                do im = 1, n_tor_local

                  im_index = i_tor_local + im - 1   ! i_tor_local is the starting index in HZ

                  index_ij = 2*n_degrees*(i-1) + 2 * (j-1) + im   ! index in the ELM matrix

                  v    = H(i,j,ms,mt)    * element%size(i,j) * HZ(im_index,mp)
                  v_s  = H_s(i,j,ms,mt)  * element%size(i,j) * HZ(im_index,mp)
                  v_t  = H_t(i,j,ms,mt)  * element%size(i,j) * HZ(im_index,mp)
                  v_p  = H(i,j,ms,mt)    * element%size(i,j) * HZ_p(im_index,mp)

                  v_ss = H_ss(i,j,ms,mt) * element%size(i,j) * HZ(im_index,mp)
                  v_tt = H_tt(i,j,ms,mt) * element%size(i,j) * HZ(im_index,mp)
                  v_st = H_st(i,j,ms,mt) * element%size(i,j) * HZ(im_index,mp)

                  v_x = (  y_t(ms,mt) * v_s - y_s(ms,mt) * v_t) / xjac
                  v_y = (- x_t(ms,mt) * v_s + x_s(ms,mt) * v_t) / xjac

                  v_xx = (v_ss * y_t(ms,mt)**2 - 2.d0*v_st * y_s(ms,mt)*y_t(ms,mt) + v_tt * y_s(ms,mt)**2  &
                      + v_s * (y_st(ms,mt)*y_t(ms,mt) - y_tt(ms,mt)*y_s(ms,mt) )                          &
                      + v_t * (y_st(ms,mt)*y_s(ms,mt) - y_ss(ms,mt)*y_t(ms,mt) ) )  / xjac**2             &
                      - xjac_x * (v_s * y_t(ms,mt) - v_t * y_s(ms,mt)) / xjac**2

                  v_yy = (v_ss * x_t(ms,mt)**2 - 2.d0*v_st * x_s(ms,mt)*x_t(ms,mt) + v_tt * x_s(ms,mt)**2  &
                      + v_s * (x_st(ms,mt)*x_t(ms,mt) - x_tt(ms,mt)*x_s(ms,mt) )                          &
                      + v_t * (x_st(ms,mt)*x_s(ms,mt) - x_ss(ms,mt)*x_t(ms,mt) ) )     / xjac**2          &
                      - xjac_y * (- v_s * x_t(ms,mt) + v_t * x_s(ms,mt) ) / xjac**2


                  Bgrad_v_star = 0.d0
                  if (filter_parallel .gt. 0.d0) Bgrad_v_star = ( F0 / x_g(ms,mt) * v_p  +  v_x  * psi_y - v_y  * psi_x ) / x_g(ms,mt)


                  do k=1,n_vertex_max
                    do l=1,n_degrees

                      do in = 1, n_tor_local

                        in_index = i_tor_local + in - 1

                        index_kl = 2*n_degrees*(k-1) + 2 * (l-1) + in   ! index in the ELM matrix

                        p   = h(k,l,ms,mt)     * element%size(k,l) * HZ(in_index,mp)
                        p_s = h_s(k,l,ms,mt)   * element%size(k,l) * HZ(in_index,mp)
                        p_t = h_t(k,l,ms,mt)   * element%size(k,l) * HZ(in_index,mp)
                        p_p = h(k,l,ms,mt)     * element%size(k,l) * HZ_p(in_index,mp)

                        p_ss = h_ss(k,l,ms,mt) * element%size(k,l) * HZ(in_index,mp)
                        p_tt = h_tt(k,l,ms,mt) * element%size(k,l) * HZ(in_index,mp)
                        p_st = h_st(k,l,ms,mt) * element%size(k,l) * HZ(in_index,mp)

                        p_x = (  y_t(ms,mt) * p_s - y_s(ms,mt) * p_t) / xjac
                        p_y = (- x_t(ms,mt) * p_s + x_s(ms,mt) * p_t) / xjac

                        p_xx = (p_ss * y_t(ms,mt)**2 - 2.d0*p_st * y_s(ms,mt)*y_t(ms,mt) + p_tt * y_s(ms,mt)**2  &
                            + p_s * (y_st(ms,mt)*y_t(ms,mt) - y_tt(ms,mt)*y_s(ms,mt) )                          &
                            + p_t * (y_st(ms,mt)*y_s(ms,mt) - y_ss(ms,mt)*y_t(ms,mt) ) )  / xjac**2             &
                            - xjac_x * (p_s * y_t(ms,mt) - p_t * y_s(ms,mt)) / xjac**2

                        p_yy = (p_ss * x_t(ms,mt)**2 - 2.d0*p_st * x_s(ms,mt)*x_t(ms,mt) + p_tt * x_s(ms,mt)**2  &
                            + p_s * (x_st(ms,mt)*x_t(ms,mt) - x_tt(ms,mt)*x_s(ms,mt) )                          &
                            + p_t * (x_st(ms,mt)*x_s(ms,mt) - x_ss(ms,mt)*x_t(ms,mt) ) )     / xjac**2          &
                            - xjac_y * (- p_s * x_t(ms,mt) + p_t * x_s(ms,mt) ) / xjac**2

                        Bgrad_p = 0.d0
                        if (filter_parallel .gt. 0.d0) Bgrad_p = ( F0 / x_g(ms,mt) * p_p +  p_x * psi_y - p_y * psi_x ) / x_g(ms,mt)

                        ELM(index_ij,index_kl) = ELM(index_ij,index_kl) &
                        
                                              + p * v * xjac * x_g(ms,mt) * wst &

                                              + filter          * (p_x * v_x + p_y * v_y) * xjac * x_g(ms,mt) * wst &

                                              + filter_hyper    * (v_xx + v_x/x_g(ms,mt) + v_yy)*(p_xx + p_x/x_g(ms,mt) + p_yy) * xjac * x_g(ms,mt) * wst &

                                              + filter_parallel * Bgrad_v_star * Bgrad_p / BB2 * xjac * x_g(ms,mt) * wst
                      enddo
                    enddo
                  enddo
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo

      ! Save contribution of this element in MUMPS format
      do i=1,n_vertex_max

        inode = element_list%element(i_elm)%vertex(i)
      
        do j=1,n_degrees

          do im =1, n_tor_local
        
            index_ij = 2*n_degrees*(i-1) + 2 * (j-1) + im   ! index in the ELM matrix

            index_large_i = 2*(node_list%node(inode)%index(j)-1) + im   ! base index in the main matrix

            do k=1,n_vertex_max
          
              knode = element_list%element(i_elm)%vertex(k)
            
              do l=1,n_degrees

                do in =1, n_tor_local
            
                  index_kl = 2*n_degrees*(k-1) + 2 * (l-1) + in   ! index in the ELM matrix

                  index_large_k = 2*(node_list%node(knode)%index(l)-1) + in   ! base index in the main matrix

                ! Explicitly calculate the index

                  ilarge = in + 2*(l-1) + 2*(k-1)*n_degrees      &
                          
                        + 2*(im-1)* n_vertex_max*n_degrees       &
                        
                        + 4*(j-1) * n_vertex_max*n_degrees       &
                        
                        + 4*(i-1) * n_vertex_max*n_degrees**2    &
                        
                        + 4*(i_elm-1)*((n_vertex_max*n_degrees)**2 )


    !$omp critical
                  a_mat%irn(ilarge)     = index_large_i
                  a_mat%jcn(ilarge)     = index_large_k

                  if( fix_axis_nodes .and.  (node_list%node(inode)%axis_node .and. (j .eq. 3 .or. j .eq. 4)) &
                    .and. (index_large_i .eq. index_large_k) ) then
                      a_mat%val(ilarge)   = 1.d12
                  else
                      a_mat%val(ilarge)     = ELM(index_ij,index_kl) * TWOPI / real(n_plane,8)
                  endif
    !$omp end critical

                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
      do m=1,n_vertex_max
        call dealloc_node(nodes(m))
      enddo
    enddo
    !$omp end parallel do
    ilarge = nz_AA

    if (apply_dirichlet_condition) then

      do i=1,node_list%n_nodes
        
        if ((node_list%node(i)%boundary .eq. 2) .or. (node_list%node(i)%boundary .eq. 3) .or. &
            (node_list%node(i)%boundary .eq. 5) .or. (node_list%node(i)%boundary .eq. 9)) then

          do j=1,3,2             ! order
            do k=1,2             ! variables

              index1 = node_list%node(i)%index(j)

              ilarge = ilarge + 1

              a_mat%irn(ilarge)     = 2*(index1-1) + k
              a_mat%jcn(ilarge)     = 2*(index1-1) + k
              a_mat%val(ilarge)     = 1.d12
            enddo
          enddo

        elseif ((node_list%node(i)%boundary .eq. 1) .or. (node_list%node(i)%boundary .eq. 3) .or. &
                (node_list%node(i)%boundary .eq. 4) .or. (node_list%node(i)%boundary .eq. 9)) then

          do j=1,2               ! order
            do k=1,2             ! variables

              index1 = node_list%node(i)%index(j)

              ilarge = ilarge + 1

              a_mat%irn(ilarge)     = 2*(index1-1) + k
              a_mat%jcn(ilarge)     = 2*(index1-1) + k
              a_mat%val(ilarge)     = 1.d12
            enddo
          enddo
      
        endif
      enddo

    endif

    nz_AA = ilarge

    a_mat%ng  = n_AA
    a_mat%nnz = nz_AA

  end if

  ! MPI broadcast the ng and nnz
  call MPI_Bcast(a_mat%ng, 1, MPI_INTEGER, 0, this_mpi_comm_world, ierr)
  call MPI_Bcast(a_mat%nnz, 1, MPI_INTEGER, 0, this_mpi_comm_world, ierr)

  end subroutine assemble_projection_matrix


  subroutine assemble_projection_matrix_n0(node_list, element_list, n_tor_local, i_tor_local,           &
                                  this_mpi_comm_world, this_mpi_comm_n, this_mpi_comm_master,  &
                                  a_mat, area, volume,                                     &
                                  filter, filter_hyper, filter_parallel,                       &
                                  integral_weights, skip_factorisation, do_zonal,              &
                                  apply_dirichlet_condition_in)
  use phys_module, only : F0, TWOPI, fix_axis_nodes
  use data_structure
  use basis_at_gaussian
  use mod_basisfunctions
  use mpi_mod
  use, intrinsic :: ieee_exceptions
  implicit none

  type (type_node_list), intent(in)    :: node_list !< A copy of the node list which will be used to save variables
  type (type_element_list), intent(in) :: element_list
  integer, intent(in)                  :: n_tor_local, i_tor_local
  integer, intent(in)                  :: this_mpi_comm_world, this_mpi_comm_n, this_mpi_comm_master
  real*8, intent(in)                   :: filter
  real*8, intent(in)                   :: filter_hyper
  real*8, intent(in)                   :: filter_parallel
  logical, intent(in), optional        :: do_zonal
  logical, intent(in), optional        :: skip_factorisation
  logical, intent(in), optional        :: apply_dirichlet_condition_in
  real*8, intent(out), dimension(:), allocatable :: integral_weights !< these will be filled with the volume of each basis function
  real*8, intent(out)                  :: area, volume
  type (type_element)      :: element
  type (type_node)         :: nodes(n_vertex_max)

  real*8, allocatable      :: ELM(:,:)
  real*8     :: wgauss2(n_gauss)
  real*8, dimension(n_gauss,n_gauss) :: x_g,   x_s,   x_t,   x_ss,   x_tt,   x_st
  real*8, dimension(n_gauss,n_gauss) :: y_g,   y_s,   y_t,   y_ss,   y_tt,   y_st
  real*8, dimension(n_gauss,n_gauss) :: psi_g, psi_s, psi_t, psi_ss, psi_tt, psi_st
  real*8     :: v, v_s, v_t, v_ss, v_st, v_tt, v_x, v_y, v_xx, v_yy, v_p
  real*8     :: p, p_s, p_t, p_ss, p_st, p_tt, p_x, p_y, p_xx, p_yy, p_p
  real*8     :: wst, xjac, xjac_x, xjac_y, psi_x, psi_y
  real*8     :: Bgrad_p, Bgrad_v_star, BB2
  real*8     :: filter_n0, filter_hyper_n0, filter_parallel_n0, zonal_factor
  integer    :: i, j, k, l, m, in, im, ilarge, index_large_i, index_large_k, inode, knode
  integer    :: nz_AA, n_AA, nz_bnd, i_elm, index_ij, index_kl, im_index, in_index, index1, index2, index_rhs
  integer    :: ms, mt, mp, my_id, my_id_n, my_id_master, ierr, MPI_COMM_MUMPS
  logical    :: halt(size(IEEE_USUAL,1)), do_facto
  logical    :: apply_dirichlet_condition, apply_zonal
  real*8, dimension(n_vertex_max,n_degrees) :: basisfunction_volume

  type (type_SP_MATRIX) :: a_mat  !< Projection matrix

  ! We need a separate communicator to be able to run multiple MUMPSes
  call MPI_Comm_dup(this_mpi_comm_n, MPI_COMM_MUMPS, ierr)
  a_mat%comm = MPI_COMM_MUMPS

  call MPI_COMM_RANK(this_mpi_comm_world,  my_id  , ierr)
  call MPI_COMM_RANK(a_mat%comm,           my_id_n, ierr)

  apply_zonal = .false.
  if (present(do_zonal)) apply_zonal = do_zonal

  apply_dirichlet_condition = .true.
  if (present(apply_dirichlet_condition_in)) apply_dirichlet_condition = apply_dirichlet_condition_in
  if (apply_zonal)  apply_dirichlet_condition = .true.


  nz_AA = 4 * element_list%n_elements * (n_vertex_max * n_degrees)**2
  n_AA =  2 * maxval(node_list%node(1:node_list%n_nodes)%index(4))

  nz_bnd = 0
  if (apply_dirichlet_condition) then
    do i=1,node_list%n_nodes
      if (node_list%node(i)%boundary .eq. 1) nz_bnd = nz_bnd + 4
      if (node_list%node(i)%boundary .eq. 2) nz_bnd = nz_bnd + 4
      if (node_list%node(i)%boundary .eq. 3) nz_bnd = nz_bnd + 8
      if (node_list%node(i)%boundary .eq. 4) nz_bnd = nz_bnd + 4
      if (node_list%node(i)%boundary .eq. 5) nz_bnd = nz_bnd + 4
      if (node_list%node(i)%boundary .eq. 9) nz_bnd = nz_bnd + 8
    enddo
  endif

  ! Only perform the construction of the matrix on the host
  if (my_id_n .eq. 0) then

      allocate(ELM(2*n_vertex_max*n_degrees,2*n_vertex_max*n_degrees))
      allocate(a_mat%val(nz_AA+nz_bnd),a_mat%irn(nz_AA+nz_bnd),a_mat%jcn(nz_AA+nz_bnd))
      
      a_mat%irn     = 0
      a_mat%jcn     = 0
      a_mat%val     = 0.d0

      allocate(integral_weights(n_AA))
    
      integral_weights = 0.d0
      
      ! Copy wgauss into wgauss2 to get around gfortran not recognizing it as a shared
      ! thing https://groups.google.com/forum/#!topic/comp.lang.fortran/VKhoAm8m9KE
      wgauss2 = wgauss

      area   = 0.
      volume = 0.

      filter_n0             = filter 
      filter_hyper_n0       = filter_hyper 
      filter_parallel_n0    = filter_parallel 
      zonal_factor          = 1.d0

      write(*,*) '*************************************************'
      write(*,*) '* constructing particle projection matrix (n=0) *'
      write(*,*) '*************************************************'
      write(*,*)  ' n_AA = ',n_AA
      write(*,'(2I3,A,3e12.4)') my_id, my_id_n,'  filters (n=0) : ',filter_n0, filter_hyper_n0, filter_parallel_n0
      
      if (apply_zonal)               write(*,*) 'using n=0 zonal flow equations'
      if (apply_dirichlet_condition) write(*,*) 'applying Dirichlet conditions'

    !$omp parallel do default(none) &
    !$omp shared(element_list, node_list, n_tor_local, i_tor_local,                       &
    !$omp        apply_dirichlet_condition, zonal_factor, apply_zonal,                    &
    !$omp        H, H_s, H_t, H_ss, H_st, H_tt, Hz, Hz_p, a_mat, wgauss2,                 &
    !$omp        filter_n0, filter_hyper_n0, filter_parallel_n0, integral_weights, F0,    &
    !$omp        fix_axis_nodes, my_id_master)                                            &
    !$omp private(ELM, i_elm, element, i, j, k, l, ms, mt, in, im, mp,                    &
    !$omp         x_g, x_s, x_t, x_ss, x_st, x_tt,                                        &
    !$omp         y_g, y_s, y_t, y_ss, y_st, y_tt,                                        &
    !$omp         psi_g, psi_s, psi_t, psi_ss, psi_tt, psi_st,                            &
    !$omp         v, v_s, v_t, v_ss, v_st, v_tt, v_x, v_y, v_xx, v_yy, v_p,               &
    !$omp         p, p_s, p_t, p_ss, p_st, p_tt, p_x, p_y, p_xx, p_yy, p_p,               &
    !$omp         wst, xjac, xjac_x, xjac_y, psi_x, psi_y, BB2, Bgrad_p, Bgrad_v_star,    &
    !$omp         index_ij, index_kl, ilarge, in_index, im_index, basisfunction_volume,   &
    !$omp         inode, index_large_i, knode, index_large_k, index_rhs)                  &
    !$omp firstprivate(nodes)                                                             & 
    !$omp reduction(+:area,volume) schedule(static)
    do i_elm=1,element_list%n_elements
      
      ELM = 0.d0

      element = element_list%element(i_elm)
      do m=1,n_vertex_max
        call make_deep_copy_node(node_list%node(element%vertex(m)), nodes(m))
      enddo

      ! Set up gauss points in this element
      x_g = 0.d0;   x_s = 0.d0;   x_t = 0.d0;   x_ss = 0.d0;   x_st = 0.d0;   x_tt = 0.d0
      y_g = 0.d0;   y_s = 0.d0;   y_t = 0.d0;   y_ss = 0.d0;   y_st = 0.d0;   y_tt = 0.d0
      psi_g = 0.d0; psi_s = 0.d0; psi_t = 0.d0; psi_ss = 0.d0; psi_st = 0.d0; psi_tt = 0.d0

      do i=1,n_vertex_max
        do j=1,n_degrees
          do ms=1, n_gauss
            do mt=1, n_gauss
              x_g(ms,mt)  = x_g(ms,mt)  + nodes(i)%x(1,j,1) * element%size(i,j) * H(i,j,ms,mt)
              x_s(ms,mt)  = x_s(ms,mt)  + nodes(i)%x(1,j,1) * element%size(i,j) * H_s(i,j,ms,mt)
              x_t(ms,mt)  = x_t(ms,mt)  + nodes(i)%x(1,j,1) * element%size(i,j) * H_t(i,j,ms,mt)

              x_ss(ms,mt) = x_ss(ms,mt) + nodes(i)%x(1,j,1) * element%size(i,j) * H_ss(i,j,ms,mt)
              x_st(ms,mt) = x_st(ms,mt) + nodes(i)%x(1,j,1) * element%size(i,j) * H_st(i,j,ms,mt)
              x_tt(ms,mt) = x_tt(ms,mt) + nodes(i)%x(1,j,1) * element%size(i,j) * H_tt(i,j,ms,mt)

              y_g(ms,mt)  = y_g(ms,mt)  + nodes(i)%x(1,j,2) * element%size(i,j) * H(i,j,ms,mt)
              y_s(ms,mt)  = y_s(ms,mt)  + nodes(i)%x(1,j,2) * element%size(i,j) * H_s(i,j,ms,mt)
              y_t(ms,mt)  = y_t(ms,mt)  + nodes(i)%x(1,j,2) * element%size(i,j) * H_t(i,j,ms,mt)

              y_ss(ms,mt) = y_ss(ms,mt) + nodes(i)%x(1,j,2) * element%size(i,j) * H_ss(i,j,ms,mt)
              y_st(ms,mt) = y_st(ms,mt) + nodes(i)%x(1,j,2) * element%size(i,j) * H_st(i,j,ms,mt)
              y_tt(ms,mt) = y_tt(ms,mt) + nodes(i)%x(1,j,2) * element%size(i,j) * H_tt(i,j,ms,mt)

              psi_g(ms,mt)  = psi_g(ms,mt)  + nodes(i)%values(1,j,1) * element%size(i,j) * H(i,j,ms,mt)
              psi_s(ms,mt)  = psi_s(ms,mt)  + nodes(i)%values(1,j,1) * element%size(i,j) * H_s(i,j,ms,mt)
              psi_t(ms,mt)  = psi_t(ms,mt)  + nodes(i)%values(1,j,1) * element%size(i,j) * H_t(i,j,ms,mt)

              psi_ss(ms,mt) = psi_ss(ms,mt) + nodes(i)%values(1,j,1) * element%size(i,j) * H_ss(i,j,ms,mt)
              psi_st(ms,mt) = psi_st(ms,mt) + nodes(i)%values(1,j,1) * element%size(i,j) * H_st(i,j,ms,mt)
              psi_tt(ms,mt) = psi_tt(ms,mt) + nodes(i)%values(1,j,1) * element%size(i,j) * H_tt(i,j,ms,mt)

            enddo
          enddo
        enddo
      enddo

      basisfunction_volume = 0.d0

      do ms=1, n_gauss
        do mt=1, n_gauss

          wst = wgauss2(ms)*wgauss2(mt)
          xjac =  x_s(ms,mt)*y_t(ms,mt) - x_t(ms,mt)*y_s(ms,mt)

          xjac_x  = (x_ss(ms,mt)*y_t(ms,mt)**2 - y_ss(ms,mt)*x_t(ms,mt)*y_t(ms,mt) - 2.d0*x_st(ms,mt)*y_s(ms,mt)*y_t(ms,mt)   &
                  + y_st(ms,mt)*(x_s(ms,mt)*y_t(ms,mt) + x_t(ms,mt)*y_s(ms,mt))                                               &
                  + x_tt(ms,mt)*y_s(ms,mt)**2 - y_tt(ms,mt)*x_s(ms,mt)*y_s(ms,mt)) / xjac

          xjac_y  = (y_tt(ms,mt)*x_s(ms,mt)**2 - x_tt(ms,mt)*y_s(ms,mt)*x_s(ms,mt) - 2.d0*y_st(ms,mt)*x_t(ms,mt)*x_s(ms,mt)   &
                  + x_st(ms,mt)*(y_t(ms,mt)*x_s(ms,mt) + y_s(ms,mt)*x_t(ms,mt))                                               &
                  + y_ss(ms,mt)*x_t(ms,mt)**2 - x_ss(ms,mt)*y_t(ms,mt)*x_t(ms,mt)) / xjac

          psi_x = (  y_t(ms,mt) * psi_s(ms,mt) - y_s(ms,mt) * psi_t(ms,mt)) / xjac
          psi_y = (- x_t(ms,mt) * psi_s(ms,mt) + x_s(ms,mt) * psi_t(ms,mt)) / xjac

          BB2 = 1.d0
          if (filter_parallel_n0 .gt. 0.d0) BB2 = (F0*F0 + psi_x * psi_x + psi_y * psi_y )/x_g(ms,mt)**2

          area   = area   + xjac * wst
          volume = volume + TWOPI * x_g(ms,mt) * xjac * wst

          do i=1,n_vertex_max
            do j=1,n_degrees

              index_ij = 2*n_degrees*(i-1) + 2*(j-1) + 1   ! index in the ELM matrix

              v    = H(i,j,ms,mt)    * element%size(i,j) 
              v_s  = H_s(i,j,ms,mt)  * element%size(i,j)
              v_t  = H_t(i,j,ms,mt)  * element%size(i,j)
              v_p  = 0.d0

              v_ss = H_ss(i,j,ms,mt) * element%size(i,j)
              v_tt = H_tt(i,j,ms,mt) * element%size(i,j)
              v_st = H_st(i,j,ms,mt) * element%size(i,j)

              v_x = (  y_t(ms,mt) * v_s - y_s(ms,mt) * v_t) / xjac
              v_y = (- x_t(ms,mt) * v_s + x_s(ms,mt) * v_t) / xjac

              v_xx = (v_ss * y_t(ms,mt)**2 - 2.d0*v_st * y_s(ms,mt)*y_t(ms,mt) + v_tt * y_s(ms,mt)**2  &
                  + v_s * (y_st(ms,mt)*y_t(ms,mt) - y_tt(ms,mt)*y_s(ms,mt) )                          &
                  + v_t * (y_st(ms,mt)*y_s(ms,mt) - y_ss(ms,mt)*y_t(ms,mt) ) )  / xjac**2             &
                  - xjac_x * (v_s * y_t(ms,mt) - v_t * y_s(ms,mt)) / xjac**2

              v_yy = (v_ss * x_t(ms,mt)**2 - 2.d0*v_st * x_s(ms,mt)*x_t(ms,mt) + v_tt * x_s(ms,mt)**2  &
                  + v_s * (x_st(ms,mt)*x_t(ms,mt) - x_tt(ms,mt)*x_s(ms,mt) )                          &
                  + v_t * (x_st(ms,mt)*x_s(ms,mt) - x_ss(ms,mt)*x_t(ms,mt) ) )     / xjac**2          &
                  - xjac_y * (- v_s * x_t(ms,mt) + v_t * x_s(ms,mt) ) / xjac**2

              Bgrad_v_star = 0.d0
              if (filter_parallel_n0 .gt. 0.d0) Bgrad_v_star = ( F0 / x_g(ms,mt) * v_p  +  v_x  * psi_y - v_y  * psi_x ) / x_g(ms,mt)

              basisfunction_volume(i,j) = basisfunction_volume(i,j) +  v * TWOPI * x_g(ms,mt) * xjac * wst

              do k=1,n_vertex_max
                do l=1,n_degrees

                  index_kl = 2*n_degrees*(k-1) + 2*(l-1) + 1   ! index in the ELM matrix

                  p   = h(k,l,ms,mt)     * element%size(k,l) 
                  p_s = h_s(k,l,ms,mt)   * element%size(k,l)
                  p_t = h_t(k,l,ms,mt)   * element%size(k,l)
                  p_p = 0.d0

                  p_ss = h_ss(k,l,ms,mt) * element%size(k,l)
                  p_tt = h_tt(k,l,ms,mt) * element%size(k,l)
                  p_st = h_st(k,l,ms,mt) * element%size(k,l)

                  p_x = (  y_t(ms,mt) * p_s - y_s(ms,mt) * p_t) / xjac
                  p_y = (- x_t(ms,mt) * p_s + x_s(ms,mt) * p_t) / xjac

                  p_xx = (p_ss * y_t(ms,mt)**2 - 2.d0*p_st * y_s(ms,mt)*y_t(ms,mt) + p_tt * y_s(ms,mt)**2  &
                      + p_s * (y_st(ms,mt)*y_t(ms,mt) - y_tt(ms,mt)*y_s(ms,mt) )                          &
                      + p_t * (y_st(ms,mt)*y_s(ms,mt) - y_ss(ms,mt)*y_t(ms,mt) ) )  / xjac**2             &
                      - xjac_x * (p_s * y_t(ms,mt) - p_t * y_s(ms,mt)) / xjac**2

                  p_yy = (p_ss * x_t(ms,mt)**2 - 2.d0*p_st * x_s(ms,mt)*x_t(ms,mt) + p_tt * x_s(ms,mt)**2  &
                      + p_s * (x_st(ms,mt)*x_t(ms,mt) - x_tt(ms,mt)*x_s(ms,mt) )                          &
                      + p_t * (x_st(ms,mt)*x_s(ms,mt) - x_ss(ms,mt)*x_t(ms,mt) ) )     / xjac**2          &
                      - xjac_y * (- p_s * x_t(ms,mt) + p_t * x_s(ms,mt) ) / xjac**2

                  Bgrad_p = 0.d0
                  if (filter_parallel_n0 .gt. 0.d0) Bgrad_p = ( F0 / x_g(ms,mt) * p_p +  p_x * psi_y - p_y * psi_x ) / x_g(ms,mt)


                  if (apply_zonal) then

                    ELM(index_ij,index_kl)     = ELM(index_ij,index_kl)     + p * v        * xjac * x_g(ms,mt) * wst &

                                              + filter_n0       * (p_x * v_x + p_y * v_y) * xjac * x_g(ms,mt) * wst &

                                              + filter_hyper_n0 * (v_xx + v_x/x_g(ms,mt) + v_yy)*(p_xx + p_x/x_g(ms,mt) + p_yy) * xjac * x_g(ms,mt) * wst

                    ELM(index_ij,index_kl+1)   = ELM(index_ij,index_kl+1) &

                                              + filter_n0  * (p_x * v_x + p_y * v_y) * xjac * x_g(ms,mt) * wst 


                    ELM(index_ij+1,index_kl)   = ELM(index_ij+1,index_kl)   + p * v   * xjac * x_g(ms,mt) * wst * zonal_factor

                    ELM(index_ij+1,index_kl+1) = ELM(index_ij+1,index_kl+1) &

                                              - filter_parallel_n0 * Bgrad_v_star * Bgrad_p / BB2 * xjac * x_g(ms,mt) * wst

                  else

                    ELM(index_ij,index_kl)     = ELM(index_ij,index_kl)     + p * v           * xjac * x_g(ms,mt) * wst &

                                              + filter_n0          * (p_x * v_x + p_y * v_y) * xjac * x_g(ms,mt) * wst &

                                              + filter_hyper_n0    * (v_xx + v_x/x_g(ms,mt) + v_yy)*(p_xx + p_x/x_g(ms,mt) + p_yy) * xjac * x_g(ms,mt) * wst &

                                              + filter_parallel_n0 * Bgrad_v_star * Bgrad_p / BB2 * xjac * x_g(ms,mt) * wst

                    ELM(index_ij+1,index_kl+1) = ELM(index_ij+1,index_kl+1)  + p * v      * xjac * x_g(ms,mt) * wst       ! (dummy equation)
                    
                  endif
                  
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo

      ! Save contribution of this element in MUMPS format
      do i=1,n_vertex_max

        inode = element_list%element(i_elm)%vertex(i)
      
        do j=1,n_degrees

          index_rhs = 2*(node_list%node(inode)%index(j)-1) + 1   ! base index in the main matrix

          integral_weights(index_rhs) = integral_weights(index_rhs) + basisfunction_volume(i,j)

          do im =1, 2
        
            index_ij = 2*n_degrees*(i-1) + 2 * (j-1) + im   ! index in the ELM matrix

            index_large_i = 2*(node_list%node(inode)%index(j)-1) + im   ! base index in the main matrix

            do k=1,n_vertex_max
          
              knode = element_list%element(i_elm)%vertex(k)
            
              do l=1,n_degrees

                do in =1, 2
            
                  index_kl = 2*n_degrees*(k-1) + 2 * (l-1) + in   ! index in the ELM matrix

                  index_large_k = 2*(node_list%node(knode)%index(l)-1) + in   ! base index in the main matrix

                ! Explicitly calculate the index

                  ilarge = in + 2*(l-1) + 2*(k-1)*n_degrees      &
                          
                        + 2*(im-1)* n_vertex_max*n_degrees       &
                        
                        + 4*(j-1) * n_vertex_max*n_degrees       &
                        
                        + 4*(i-1) * n_vertex_max*n_degrees**2    &
                        
                        + 4*(i_elm-1)*((n_vertex_max*n_degrees)**2 )

                  a_mat%irn(ilarge) = index_large_i
                  a_mat%jcn(ilarge) = index_large_k

                  if( fix_axis_nodes .and.  (node_list%node(inode)%axis_node .and. (j .eq. 3 .or. j .eq. 4)) .and. (index_large_i .eq. index_large_k) ) then
                      a_mat%val(ilarge) = 1.d12
                  else
                      a_mat%val(ilarge) = ELM(index_ij,index_kl) * TWOPI
                  endif

                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
      do m=1,n_vertex_max
        call dealloc_node(nodes(m))
      enddo
    enddo
    !$omp end parallel do

    ! boundary conditions 

    ilarge = nz_AA

    if (apply_dirichlet_condition) then

      do i=1,node_list%n_nodes
        
        if ((node_list%node(i)%boundary .eq. 2) .or. (node_list%node(i)%boundary .eq. 3) .or. &
            (node_list%node(i)%boundary .eq. 5) .or. (node_list%node(i)%boundary .eq. 9)) then

          do j=1,3,2         ! order
            do k=1,2         ! variables

              index1 = node_list%node(i)%index(j)

              ilarge = ilarge + 1

              a_mat%irn(ilarge) = 2*(index1-1) + k
              a_mat%jcn(ilarge) = 2*(index1-1) + k
              a_mat%val(ilarge) = 1.d12
            enddo
          enddo

        elseif ((node_list%node(i)%boundary .eq. 1) .or. (node_list%node(i)%boundary .eq. 3) .or. &
                (node_list%node(i)%boundary .eq. 4) .or. (node_list%node(i)%boundary .eq. 9)) then

          do j=1,2           ! order
            do k=1,2         ! variables

              index1 = node_list%node(i)%index(j)

              ilarge = ilarge + 1

              a_mat%irn(ilarge)     = 2*(index1-1) + k
              a_mat%jcn(ilarge)     = 2*(index1-1) + k
              a_mat%val(ilarge)     = 1.d12
            enddo
          enddo
      
        endif
      enddo

    endif

    nz_AA = ilarge

    write(*,'(A,2e16.8)') 'area volume : ',area, volume

    a_mat%ng              = n_AA
    a_mat%nnz             = nz_AA
  end if

  call MPI_Bcast(a_mat%ng, 1, MPI_INTEGER, 0, this_mpi_comm_world, ierr)
  call MPI_Bcast(a_mat%nnz, 1, MPI_INTEGER, 0, this_mpi_comm_world, ierr)

  end subroutine assemble_projection_matrix_n0

end module mod_uncoupled_projection

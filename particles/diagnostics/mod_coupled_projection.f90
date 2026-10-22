module mod_coupled_projection

implicit none
private

public assemble_system, project_only

contains

  !> Constructor for project_particles
  !> Be sure to use keyword arguments when initializing, to avoid confusion
  subroutine assemble_system(node_list, element_list, mpi_comm_world, a_mat, filter, filter_hyper, filter_parallel, apply_dirichlet, system_size, n_dof)
    use data_structure, only: type_element_list, type_node_list
    use mod_sparse_matrix, only: type_SP_MATRIX
  
    ! Declare the input variables
    type(type_node_list),    intent(in)    :: node_list
    type(type_element_list), intent(in)    :: element_list

    integer,                 intent(in)    :: mpi_comm_world

    real*8,                  intent(in)    :: filter, filter_hyper, filter_parallel

    logical,                 intent(in)    :: apply_dirichlet

    ! Declare the output variables
    integer,                 intent(out)   :: system_size
    integer,                 intent(out)   :: n_dof

    type(type_SP_MATRIX),    intent(inout) :: a_mat

    ! Assemble the projection matrix
    call assemble_projection_matrix(node_list, element_list,                       &
                                    mpi_comm_world, a_mat,                         &
                                    filter, filter_hyper, filter_parallel,         &
                                    apply_dirichlet_condition_in=apply_dirichlet)

    ! Set system size and degrees of freedom
    system_size = a_mat%ng
    n_dof = system_size

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
    integer, allocatable :: recv_counts(:), recv_disp(:)
    integer :: n_rhs, n_rhs_f, i_rhs
    integer :: in_local, in_global, index_n, id_master_in_world
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
    this%node_list%n_values = n_rhs_f+n_rhs
    do i=1, this%node_list%n_nodes
      call init_node(this%node_list%node(i), n_rhs_f+n_rhs)
    enddo

    this%rhs_vec%nrhs = (n_rhs + n_rhs_f)
    this%rhs_vec%n = this%system_size

    if (this%n_dof .ne. this%rhs_vec%n) then
      write(*,*) 'FATAL : n_tor*this%n_dof .ne. this%rhs_vec%n'
      write(*,*) this%n_dof,  this%rhs_vec%n
    endif

    if (associated(this%rhs_vec%val)) then
      deallocate(this%rhs_vec%val); this%rhs_vec%val => Null()
    endif
    allocate(this%rhs_vec%val(this%rhs_vec%n * this%rhs_vec%nrhs * n_tor))
    this%rhs_vec%val = 0.d0

    do i_rhs=1,n_rhs
      i_start =  this%rhs_vec%n * (i_rhs-1)
      do i_elm=1,this%element_list%n_elements
        do i=1,n_vertex_max
          inode = this%element_list%element(i_elm)%vertex(i)
          do j=1,n_degrees
            do in=1, n_tor
              index_large_i = n_tor * (this%node_list%node(inode)%index(j)-1) + in + i_start ! base index in the main matrix + rhs index
              this%rhs_vec%val(index_large_i) = this%rhs_vec%val(index_large_i) + this%rhs(j, i, i_elm, in, i_rhs)
            enddo
          enddo
        enddo
      enddo
    enddo

    do i_rhs=1,n_rhs_f
      ! Fill projection function part
      i_start =  this%rhs_vec%n * (n_rhs + i_rhs - 1)
      do i_elm=1,this%element_list%n_elements 
        do i=1,n_vertex_max
          inode = this%element_list%element(i_elm)%vertex(i)
          do j=1,n_degrees
            do in=1, n_tor
            index_large_i = n_tor * (this%node_list%node(inode)%index(j)-1) + in + i_start ! base index in the main matrix + rhs index
            this%rhs_vec%val(index_large_i) = this%rhs_vec%val(index_large_i) + this%rhs_f(j, i, i_elm, in, i_rhs)
            enddo
          enddo
        enddo
      enddo
    enddo

    if (this%my_id .eq. 0) then
      call MPI_Reduce(MPI_IN_PLACE,    this%rhs_vec%val,this%rhs_vec%n*this%rhs_vec%nrhs, MPI_REAL8, MPI_SUM, 0, this%mpi_comm_world, ierr)
    else
      call MPI_Reduce(this%rhs_vec%val,this%rhs_vec%val,this%rhs_vec%n*this%rhs_vec%nrhs, MPI_REAL8, MPI_SUM, 0, this%mpi_comm_world, ierr)
    endif

    if (allocated(this%rhs)) this%rhs = 0.d0
    if (allocated(this%rhs_f)) this%rhs_f = 0.d0

    call ieee_get_halting_mode(IEEE_USUAL, halt)
    call ieee_set_halting_mode(IEEE_USUAL, [.false., .false., .false.])

    call solve_sparse_projection_system(this%a_mat, this%rhs_vec, this%solver)

    call ieee_set_halting_mode(IEEE_USUAL, halt)

    call MPI_BARRIER(this%mpi_comm_world, ierr)

    ! Write the solution to the node_list
    if (this%my_id .eq. 0) then

      do i_var=1,n_rhs+n_rhs_f
    
        found_nan = .false.
        
        do i=1,this%node_list%n_nodes

          do k=1,n_degrees
        
            index = this%node_list%node(i)%index(k)
            
            do i_tor=1,n_tor
              this%node_list%node(i)%deltas(i_tor, k, i_var) = this%rhs_vec%val(n_tor*(index - 1) + i_tor + this%rhs_vec%n * (i_var - 1)) - this%node_list%node(i)%values(i_tor, k, i_var)
              this%node_list%node(i)%values(i_tor, k, i_var) = this%rhs_vec%val(n_tor*(index - 1) + i_tor + this%rhs_vec%n * (i_var - 1))
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

  subroutine assemble_projection_matrix(node_list, element_list,                             &
                              this_mpi_comm_world, a_mat,                                    &
                              filter, filter_hyper, filter_parallel,                         &
                              apply_dirichlet_condition_in)
  use data_structure, only: type_element, type_element_list, type_node, type_node_list, make_deep_copy_node, dealloc_node
  use mod_sparse_matrix, only: type_SP_MATRIX
  use mod_parameters, only: n_tor, n_vertex_max, n_degrees, n_coord_tor, n_plane
  use phys_module, only : F0, TWOPI, fix_axis_nodes, mode, n_tor_fft_thresh
  use iso_c_binding

  use mpi_mod
  use basis_at_gaussian
  use mod_basisfunctions

  implicit none

#ifdef USE_FFTW
  include 'fftw3.f03'
#endif

  type (type_node_list), intent(in)    :: node_list
  type (type_element_list), intent(in) :: element_list
  integer, intent(in)                  :: this_mpi_comm_world
  real*8, intent(in)                   :: filter
  real*8, intent(in)                   :: filter_hyper
  real*8, intent(in)                   :: filter_parallel
  logical, intent(in), optional        :: apply_dirichlet_condition_in

  type (type_element)      :: element
  type (type_node)         :: nodes(n_vertex_max)

  real*8, allocatable      :: ELM(:,:)
  real*8     :: wgauss2(n_gauss)
  real*8, dimension(n_plane,n_gauss,n_gauss) :: x_g,   x_s,   x_t,   x_ss,   x_tt,   x_st
  real*8, dimension(n_plane,n_gauss,n_gauss) :: y_g,   y_s,   y_t,   y_ss,   y_tt,   y_st
  real*8, dimension(n_plane,n_gauss,n_gauss) :: psi_g, psi_s, psi_t, psi_ss, psi_tt, psi_st
  real*8     :: v, v_s, v_t, v_ss, v_st, v_tt, v_x, v_y, v_xx, v_yy, v_p
  real*8     :: p, p_s, p_t, p_ss, p_st, p_tt, p_x, p_y, p_xx, p_yy, p_p
  real*8     :: wst, area, volume, xjac, xjac_x, xjac_y, psi_x, psi_y
  real*8     :: Bgrad_p, Bgrad_v_star, BB2
  integer    :: i, j, k, l, m, in, im, ilarge, index_large_i, index_large_k, inode, knode
  integer    :: nz_AA, n_AA, nz_bnd, i_elm, index_ij, index_kl, im_index, in_index, index1
  integer    :: ms, mt, mp, my_id, ierr
  logical    :: apply_dirichlet_condition
  logical    :: do_facto
  logical    :: use_fft
  real*8, allocatable :: ELM_term1(:,:,:), ELM_term2(:,:,:)
  real*8     :: in_fft(1:n_plane)
  complex*16 :: out_fft(1:n_plane)

  integer*8 :: fftw_plan

  integer :: ik, i_loc, j_loc

  type (type_SP_MATRIX) :: a_mat  !< Projection matrix

  a_mat%comm = MPI_COMM_WORLD

  call MPI_COMM_RANK(this_mpi_comm_world,  my_id  , ierr)

  nz_AA = element_list%n_elements * (n_vertex_max * n_degrees * n_tor)**2
  n_AA  = n_tor * maxval(node_list%node(1:node_list%n_nodes)%index(4))

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

  ! n_tor // n_coord_tor

  use_fft = (n_tor > n_tor_fft_thresh)
  ! Only perform the construction of the matrix on the host
  if (my_id .eq. 0) then

    allocate(ELM(n_tor*n_vertex_max*n_degrees,n_tor*n_vertex_max*n_degrees))

    ! Allocate space for elements

    allocate(  a_mat%val(nz_AA+nz_bnd),    a_mat%irn(nz_AA+nz_bnd),    a_mat%jcn(nz_AA+nz_bnd))

    a_mat%irn = 0
    a_mat%jcn = 0
    a_mat%val = 0.d0

    ! Copy wgauss into wgauss2 to get around gfortran not recognizing it as a shared
    ! thing https://groups.google.com/forum/#!topic/comp.lang.fortran/VKhoAm8m9KE
    wgauss2 = wgauss

    write(*,*) "Starting matrix construction with: ", element_list%n_elements, " elements"
    if (use_fft) then
        write(*,*) "Using FFT-accelerated path."
#ifdef USE_FFTW
        write(*,*) " with fftw library"
#else
        write(*,*) " with custom implementation"
#endif
    else
        write(*,*) "Using direct summation path."
    end if

#ifdef USE_FFTW
    if (use_fft) call dfftw_plan_dft_r2c_1d(fftw_plan, n_plane, in_fft, out_fft, FFTW_PATIENT)
#endif

    if (apply_dirichlet_condition) write(*,*) 'applying Dirichlet conditions'
    !$omp parallel do default(none) &
    !$omp shared(element_list, node_list,                                                 &
    !$omp        H, H_s, H_t, H_ss, H_st, H_tt, Hz, Hz_p, HZ_coord, a_mat, wgauss2, mode, &
    !$omp        filter, filter_hyper, filter_parallel, F0, fix_axis_nodes, use_fft, fftw_plan) &
    !$omp private(ELM, i_elm, element, i, j, k, l, ms, mt, in, im, mp,                    &
    !$omp         x_g, x_s, x_t, x_ss, x_st, x_tt,                                        &
    !$omp         y_g, y_s, y_t, y_ss, y_st, y_tt,                                        &
    !$omp         psi_g, psi_s, psi_t, psi_ss, psi_tt, psi_st,                            &
    !$omp         v, v_s, v_t, v_ss, v_st, v_tt, v_x, v_y, v_xx, v_yy, v_p,               &
    !$omp         p, p_s, p_t, p_ss, p_st, p_tt, p_x, p_y, p_xx, p_yy, p_p,               &
    !$omp         wst, xjac, xjac_x, xjac_y, psi_x, psi_y, BB2, Bgrad_p, Bgrad_v_star,    &
    !$omp         index_ij, index_kl, ilarge, in_index, im_index,                         &
    !$omp         inode, index_large_i, knode, index_large_k,                             &
    !$omp         ELM_term1, ELM_term2, in_fft, out_fft, ik, i_loc, j_loc)                &
    !$omp firstprivate(nodes)                                                             & 
    !$omp schedule(static) 
    do i_elm=1,element_list%n_elements
      
      ELM = 0.d0
      
      if (use_fft) then
        allocate(ELM_term1(n_plane, n_vertex_max*n_degrees, n_vertex_max*n_degrees))
        ELM_term1 = 0.d0
        if (filter_parallel .gt. 0.d0) then
          allocate(ELM_term2(n_plane, n_vertex_max*n_degrees, n_vertex_max*n_degrees))
          ELM_term2 = 0.d0
        endif
      end if

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
              do mp=1,n_plane
                do in=1,n_coord_tor
                  x_g(mp,ms,mt)  = x_g(mp,ms,mt)  + nodes(i)%x(in,j,1) * element%size(i,j) * H(i,j,ms,mt)    * HZ_coord(in,mp)
                  x_s(mp,ms,mt)  = x_s(mp,ms,mt)  + nodes(i)%x(in,j,1) * element%size(i,j) * H_s(i,j,ms,mt)  * HZ_coord(in,mp)
                  x_t(mp,ms,mt)  = x_t(mp,ms,mt)  + nodes(i)%x(in,j,1) * element%size(i,j) * H_t(i,j,ms,mt)  * HZ_coord(in,mp)

                  x_ss(mp,ms,mt) = x_ss(mp,ms,mt) + nodes(i)%x(in,j,1) * element%size(i,j) * H_ss(i,j,ms,mt) * HZ_coord(in,mp)
                  x_st(mp,ms,mt) = x_st(mp,ms,mt) + nodes(i)%x(in,j,1) * element%size(i,j) * H_st(i,j,ms,mt) * HZ_coord(in,mp)
                  x_tt(mp,ms,mt) = x_tt(mp,ms,mt) + nodes(i)%x(in,j,1) * element%size(i,j) * H_tt(i,j,ms,mt) * HZ_coord(in,mp)

                  y_g(mp,ms,mt)  = y_g(mp,ms,mt)  + nodes(i)%x(in,j,2) * element%size(i,j) * H(i,j,ms,mt)    * HZ_coord(in,mp)
                  y_s(mp,ms,mt)  = y_s(mp,ms,mt)  + nodes(i)%x(in,j,2) * element%size(i,j) * H_s(i,j,ms,mt)  * HZ_coord(in,mp)
                  y_t(mp,ms,mt)  = y_t(mp,ms,mt)  + nodes(i)%x(in,j,2) * element%size(i,j) * H_t(i,j,ms,mt)  * HZ_coord(in,mp)

                  y_ss(mp,ms,mt) = y_ss(mp,ms,mt) + nodes(i)%x(in,j,2) * element%size(i,j) * H_ss(i,j,ms,mt) * HZ_coord(in,mp)
                  y_st(mp,ms,mt) = y_st(mp,ms,mt) + nodes(i)%x(in,j,2) * element%size(i,j) * H_st(i,j,ms,mt) * HZ_coord(in,mp)
                  y_tt(mp,ms,mt) = y_tt(mp,ms,mt) + nodes(i)%x(in,j,2) * element%size(i,j) * H_tt(i,j,ms,mt) * HZ_coord(in,mp)
                enddo

                do im=1,n_tor
                  psi_g(mp,ms,mt)  = psi_g(mp,ms,mt)  + nodes(i)%values(im,j,1) * element%size(i,j) * H(i,j,ms,mt)   * HZ(im,mp)
                  psi_s(mp,ms,mt)  = psi_s(mp,ms,mt)  + nodes(i)%values(im,j,1) * element%size(i,j) * H_s(i,j,ms,mt) * HZ(im,mp)
                  psi_t(mp,ms,mt)  = psi_t(mp,ms,mt)  + nodes(i)%values(im,j,1) * element%size(i,j) * H_t(i,j,ms,mt) * HZ(im,mp)

                  psi_ss(mp,ms,mt) = psi_ss(mp,ms,mt) + nodes(i)%values(im,j,1) * element%size(i,j) * H_ss(i,j,ms,mt) * HZ(im,mp)
                  psi_st(mp,ms,mt) = psi_st(mp,ms,mt) + nodes(i)%values(im,j,1) * element%size(i,j) * H_st(i,j,ms,mt) * HZ(im,mp)
                  psi_tt(mp,ms,mt) = psi_tt(mp,ms,mt) + nodes(i)%values(im,j,1) * element%size(i,j) * H_tt(i,j,ms,mt) * HZ(im,mp)
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo

      do ms=1, n_gauss
        do mt=1, n_gauss

          wst = wgauss2(ms)*wgauss2(mt)

          do mp = 1, n_plane


            xjac    = x_s(mp,ms,mt)*y_t(mp,ms,mt)  - x_t(mp,ms,mt)*y_s(mp,ms,mt)
          
            xjac_x  = (x_ss(mp,ms,mt)*y_t(mp,ms,mt)**2 - y_ss(mp,ms,mt)*x_t(mp,ms,mt)*y_t(mp,ms,mt) - 2.d0*x_st(mp,ms,mt)*y_s(mp,ms,mt)*y_t(mp,ms,mt) &
                  + y_st(mp,ms,mt)*(x_s(mp,ms,mt)*y_t(mp,ms,mt) + x_t(mp,ms,mt)*y_s(mp,ms,mt))                                                      &
                  + x_tt(mp,ms,mt)*y_s(mp,ms,mt)**2 - y_tt(mp,ms,mt)*x_s(mp,ms,mt)*y_s(mp,ms,mt)) / xjac

            xjac_y  = (y_tt(mp,ms,mt)*x_s(mp,ms,mt)**2 - x_tt(mp,ms,mt)*y_s(mp,ms,mt)*x_s(mp,ms,mt) - 2.d0*y_st(mp,ms,mt)*x_t(mp,ms,mt)*x_s(mp,ms,mt) &
                  + x_st(mp,ms,mt)*(y_t(mp,ms,mt)*x_s(mp,ms,mt) + y_s(mp,ms,mt)*x_t(mp,ms,mt))                                                      &
                  + y_ss(mp,ms,mt)*x_t(mp,ms,mt)**2 - x_ss(mp,ms,mt)*y_t(mp,ms,mt)*x_t(mp,ms,mt)) / xjac


            psi_x = (  y_t(mp,ms,mt) * psi_s(mp,ms,mt) - y_s(mp,ms,mt) * psi_t(mp,ms,mt)) / xjac
            psi_y = (- x_t(mp,ms,mt) * psi_s(mp,ms,mt) + x_s(mp,ms,mt) * psi_t(mp,ms,mt)) / xjac

            BB2 = 1.d0
            if (filter_parallel .gt. 0.d0) BB2 = (F0*F0 + psi_x * psi_x + psi_y * psi_y )/x_g(mp,ms,mt)**2

            do i=1,n_vertex_max
              do j=1,n_degrees
                i_loc = n_degrees*(i-1) + j
                v_s  = H_s(i,j,ms,mt)  * element%size(i,j)
                v_t  = H_t(i,j,ms,mt)  * element%size(i,j)
                v_p  = H(i,j,ms,mt)    * element%size(i,j)
                v    = v_p
                v_ss = H_ss(i,j,ms,mt) * element%size(i,j)
                v_tt = H_tt(i,j,ms,mt) * element%size(i,j)
                v_st = H_st(i,j,ms,mt) * element%size(i,j)
                v_x = (  y_t(mp,ms,mt) * v_s - y_s(mp,ms,mt) * v_t) / xjac
                v_y = (- x_t(mp,ms,mt) * v_s + x_s(mp,ms,mt) * v_t) / xjac

                v_xx = (v_ss * y_t(mp,ms,mt)**2 - 2.d0*v_st * y_s(mp,ms,mt)*y_t(mp,ms,mt) + v_tt * y_s(mp,ms,mt)**2  &
                    + v_s * (y_st(mp,ms,mt)*y_t(mp,ms,mt) - y_tt(mp,ms,mt)*y_s(mp,ms,mt) )                          &
                    + v_t * (y_st(mp,ms,mt)*y_s(mp,ms,mt) - y_ss(mp,ms,mt)*y_t(mp,ms,mt) ) )  / xjac**2             &
                    - xjac_x * (v_s * y_t(mp,ms,mt) - v_t * y_s(mp,ms,mt)) / xjac**2

                v_yy = (v_ss * x_t(mp,ms,mt)**2 - 2.d0*v_st * x_s(mp,ms,mt)*x_t(mp,ms,mt) + v_tt * x_s(mp,ms,mt)**2  &
                    + v_s * (x_st(mp,ms,mt)*x_t(mp,ms,mt) - x_tt(mp,ms,mt)*x_s(mp,ms,mt) )                          &
                    + v_t * (x_st(mp,ms,mt)*x_s(mp,ms,mt) - x_ss(mp,ms,mt)*x_t(mp,ms,mt) ) )     / xjac**2          &
                    - xjac_y * (- v_s * x_t(mp,ms,mt) + v_t * x_s(mp,ms,mt) ) / xjac**2


                Bgrad_v_star = 0.d0
                if (filter_parallel .gt. 0.d0) Bgrad_v_star = ( F0 / x_g(mp,ms,mt) * v_p  +  v_x  * psi_y - v_y  * psi_x ) / x_g(mp,ms,mt)


                do k=1,n_vertex_max
                  do l=1,n_degrees
                    j_loc = n_degrees*(k-1) + l
                    p_s = h_s(k,l,ms,mt)   * element%size(k,l)
                    p_t = h_t(k,l,ms,mt)   * element%size(k,l)
                    p_p = h(k,l,ms,mt)     * element%size(k,l)
                    p   = p_p
                    p_ss = h_ss(k,l,ms,mt) * element%size(k,l)
                    p_tt = h_tt(k,l,ms,mt) * element%size(k,l)
                    p_st = h_st(k,l,ms,mt) * element%size(k,l)
                    p_x = (  y_t(mp,ms,mt) * p_s - y_s(mp,ms,mt) * p_t) / xjac
                    p_y = (- x_t(mp,ms,mt) * p_s + x_s(mp,ms,mt) * p_t) / xjac

                    p_xx = (p_ss * y_t(mp,ms,mt)**2 - 2.d0*p_st * y_s(mp,ms,mt)*y_t(mp,ms,mt) + p_tt * y_s(mp,ms,mt)**2  &
                        + p_s * (y_st(mp,ms,mt)*y_t(mp,ms,mt) - y_tt(mp,ms,mt)*y_s(mp,ms,mt) )                          &
                        + p_t * (y_st(mp,ms,mt)*y_s(mp,ms,mt) - y_ss(mp,ms,mt)*y_t(mp,ms,mt) ) )  / xjac**2             &
                        - xjac_x * (p_s * y_t(mp,ms,mt) - p_t * y_s(mp,ms,mt)) / xjac**2

                    p_yy = (p_ss * x_t(mp,ms,mt)**2 - 2.d0*p_st * x_s(mp,ms,mt)*x_t(mp,ms,mt) + p_tt * x_s(mp,ms,mt)**2  &
                        + p_s * (x_st(mp,ms,mt)*x_t(mp,ms,mt) - x_tt(mp,ms,mt)*x_s(mp,ms,mt) )                          &
                        + p_t * (x_st(mp,ms,mt)*x_s(mp,ms,mt) - x_ss(mp,ms,mt)*x_t(mp,ms,mt) ) )     / xjac**2          &
                        - xjac_y * (- p_s * x_t(mp,ms,mt) + p_t * x_s(mp,ms,mt) ) / xjac**2

                    Bgrad_p = 0.d0
                    if (filter_parallel .gt. 0.d0) Bgrad_p = ( F0 / x_g(mp,ms,mt) * p_p +  p_x * psi_y - p_y * psi_x ) / x_g(mp,ms,mt)

                    if (.not. use_fft) then
                      do im = 1, n_tor
                        do in = 1, n_tor
                          index_ij = n_tor*n_degrees*(i-1) + n_tor * (j-1) + im
                          index_kl = n_tor*n_degrees*(k-1) + n_tor * (l-1) + in
                          ELM(index_ij,index_kl) = ELM(index_ij,index_kl) &
                              + (p   * HZ(in,mp)) * (v   * HZ(im,mp)) * xjac * x_g(mp,ms,mt) * wst &
                              + filter          * ((p_x * HZ(in,mp)) * (v_x * HZ(im,mp)) + (p_y * HZ(in,mp)) * (v_y * HZ(im,mp))) * xjac * x_g(mp,ms,mt) * wst &
                              + filter_hyper    * ((v_xx* HZ(im,mp)) + (v_x* HZ(im,mp))/x_g(mp,ms,mt) + (v_yy* HZ(im,mp)))* &
                                                  ((p_xx* HZ(in,mp)) + (p_x* HZ(in,mp))/x_g(mp,ms,mt) + (p_yy* HZ(in,mp))) * xjac * x_g(mp,ms,mt) * wst &
                              + filter_parallel * (Bgrad_v_star * HZ_p(im,mp)) * (Bgrad_p * HZ_p(in,mp)) / BB2 * xjac * x_g(mp,ms,mt) * wst
                        enddo
                      enddo
                    else
                      ELM_term1(mp, i_loc, j_loc) = ELM_term1(mp, i_loc, j_loc) &
                                  + p * v * xjac * x_g(mp,ms,mt) * wst &
                                  + filter          * (p_x * v_x + p_y * v_y) * xjac * x_g(mp,ms,mt) * wst &
                                  + filter_hyper    * (v_xx + v_x/x_g(mp,ms,mt) + v_yy)*(p_xx + p_x/x_g(mp,ms,mt) + p_yy) * xjac * x_g(mp,ms,mt) * wst
                      if (filter_parallel .gt. 0.d0) &
                        ELM_term2(mp, i_loc, j_loc) = ELM_term2(mp, i_loc, j_loc) &
                                    + filter_parallel * Bgrad_v_star * Bgrad_p / BB2 * xjac * x_g(mp,ms,mt) * wst
                    endif
                  enddo
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
    
      if (use_fft) then
        ! FFT Assembly Block
        do i_loc = 1, n_vertex_max*n_degrees
          do j_loc = 1, n_vertex_max*n_degrees
            ! --- Process ELM_term1 ---
            in_fft(1:n_plane) = ELM_term1(1:n_plane, i_loc, j_loc)
#ifdef USE_FFTW
            call dfftw_execute_dft_r2c(fftw_plan, in_fft, out_fft)
#else
            call my_fft(in_fft, out_fft, n_plane)
#endif
            do m = 1, (n_tor+1)/2
              do k = 1, (n_tor+1)/2
                im_index = max(2*(m-1),1)
                in_index = max(2*(k-1),1)
                index_ij = n_tor*(i_loc-1) + im_index
                index_kl = n_tor*(j_loc-1) + in_index
                
                l = abs(k-m)
                ELM(index_ij,   index_kl)   = ELM(index_ij,   index_kl)   + real(out_fft(l+1))
                if (m > 1) ELM(index_ij+1, index_kl)   = ELM(index_ij+1, index_kl)   + imag(out_fft(l+1))*sign(1.d0,real(k-m))
                if (k > 1) ELM(index_ij,   index_kl+1) = ELM(index_ij,   index_kl+1) - imag(out_fft(l+1))*sign(1.d0,real(k-m))
                if (m > 1 .and. k > 1) ELM(index_ij+1, index_kl+1) = ELM(index_ij+1, index_kl+1) + real(out_fft(l+1))
                
                l = k+m-2
                if (l >= 0 .and. l < n_plane/2) then
                  ELM(index_ij,   index_kl)   = ELM(index_ij,   index_kl)   + real(out_fft(l+1))
                  if (m > 1) ELM(index_ij+1, index_kl)   = ELM(index_ij+1, index_kl)   - imag(out_fft(l+1))
                  if (k > 1) ELM(index_ij,   index_kl+1) = ELM(index_ij,   index_kl+1) - imag(out_fft(l+1))
                  if (m > 1 .and. k > 1) ELM(index_ij+1, index_kl+1) = ELM(index_ij+1, index_kl+1) - real(out_fft(l+1))
                endif
              enddo
            enddo

            ! --- Process ELM_term2 ---
            if (filter_parallel .gt. 0.d0) then
              in_fft(1:n_plane) = ELM_term2(1:n_plane, i_loc, j_loc)
#ifdef USE_FFTW
              call dfftw_execute_dft_r2c(fftw_plan, in_fft, out_fft)
#else
              call my_fft(in_fft, out_fft, n_plane)
#endif
              do m = 1, (n_tor+1)/2 
                do k = 1, (n_tor+1)/2
                  im_index = max(2*(m-1),1)
                  in_index = max(2*(k-1),1)
                  index_ij = n_tor*(i_loc-1) + im_index
                  index_kl = n_tor*(j_loc-1) + in_index
                  
                  l = abs(k-m)
                  ELM(index_ij,   index_kl)   = ELM(index_ij,   index_kl)   + real(out_fft(l+1)) * mode(in_index) * mode(im_index)
                  if (m > 1) ELM(index_ij+1, index_kl)   = ELM(index_ij+1, index_kl)   + imag(out_fft(l+1))*sign(1.d0,real(k-m)) * mode(in_index) * mode(im_index)
                  if (k > 1) ELM(index_ij,   index_kl+1) = ELM(index_ij,   index_kl+1) - imag(out_fft(l+1))*sign(1.d0,real(k-m)) * mode(in_index) * mode(im_index)
                  if (m > 1 .and. k > 1) ELM(index_ij+1, index_kl+1) = ELM(index_ij+1, index_kl+1) + real(out_fft(l+1)) * mode(in_index) * mode(im_index)
                  
                  l = k+m-2
                  if (l >= 0 .and. l < n_plane/2) then
                    ELM(index_ij,   index_kl)   = ELM(index_ij,   index_kl)   - real(out_fft(l+1)) * mode(in_index) * mode(im_index)
                    if (m > 1) ELM(index_ij+1, index_kl)   = ELM(index_ij+1, index_kl)   + imag(out_fft(l+1)) * mode(in_index) * mode(im_index)
                    if (k > 1) ELM(index_ij,   index_kl+1) = ELM(index_ij,   index_kl+1) + imag(out_fft(l+1)) * mode(in_index) * mode(im_index)
                    if (m > 1 .and. k > 1) ELM(index_ij+1, index_kl+1) = ELM(index_ij+1, index_kl+1) + real(out_fft(l+1)) * mode(in_index) * mode(im_index)
                  endif
                enddo
              enddo

            endif
          enddo
        enddo
        ELM = 0.5d0 * ELM
        deallocate(ELM_term1)
        if (filter_parallel .gt. 0.d0) deallocate(ELM_term2)
      endif

      ! Save contribution of this element in MUMPS format
      do i=1,n_vertex_max

        inode = element_list%element(i_elm)%vertex(i)
      
        do j=1,n_degrees

          do im =1, n_tor
        
            index_ij = n_tor*n_degrees*(i-1) + n_tor * (j-1) + im   ! index in the ELM matrix

            index_large_i = n_tor*(node_list%node(inode)%index(j)-1) + im   ! base index in the main matrix

            do k=1,n_vertex_max
          
              knode = element_list%element(i_elm)%vertex(k)
            
              do l=1,n_degrees

                do in =1, n_tor
            
                  index_kl = n_tor*n_degrees*(k-1) + n_tor * (l-1) + in   ! index in the ELM matrix

                  index_large_k = n_tor*(node_list%node(knode)%index(l)-1) + in   ! base index in the main matrix

                  ! Explicitly calculate the index
                  ilarge = in + n_tor*(l-1) + n_tor*(k-1)*n_degrees             &

                        + (im-1)    *   n_vertex_max *  n_degrees * n_tor       &
                        
                        + (j-1)     *   n_vertex_max *  n_degrees * n_tor **2   &
                        
                        + (i-1)     *   n_vertex_max * (n_degrees * n_tor)**2   &
                        
                        + (i_elm-1) * ((n_vertex_max *  n_degrees * n_tor)**2 )

                  a_mat%irn(ilarge)     = index_large_i
                  a_mat%jcn(ilarge)     = index_large_k

                  if( fix_axis_nodes .and.  (node_list%node(inode)%axis_node .and. (j .eq. 3 .or. j .eq. 4)) &
                    .and. (index_large_i .eq. index_large_k) ) then
                      a_mat%val(ilarge)   = 1.d12
                  else
                      a_mat%val(ilarge)     = ELM(index_ij,index_kl) * TWOPI / real(n_plane,8)
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
    ilarge = nz_AA

    if (apply_dirichlet_condition) then

      do i=1,node_list%n_nodes
        
        if ((node_list%node(i)%boundary .eq. 2) .or. (node_list%node(i)%boundary .eq. 3) .or. &
            (node_list%node(i)%boundary .eq. 5) .or. (node_list%node(i)%boundary .eq. 9)) then

          do j=1,3,2             ! order
            do k=1,2             ! variables

              index1 = node_list%node(i)%index(j)

              ilarge = ilarge + 1

              a_mat%irn(ilarge)     = n_tor*(index1-1) + k
              a_mat%jcn(ilarge)     = n_tor*(index1-1) + k
              a_mat%val(ilarge)     = 1.d12
            enddo
          enddo

        elseif ((node_list%node(i)%boundary .eq. 1) .or. (node_list%node(i)%boundary .eq. 3) .or. &
                (node_list%node(i)%boundary .eq. 4) .or. (node_list%node(i)%boundary .eq. 9)) then

          do j=1,2               ! order
            do k=1,2             ! variables

              index1 = node_list%node(i)%index(j)

              ilarge = ilarge + 1

              a_mat%irn(ilarge)     = n_tor*(index1-1) + k
              a_mat%jcn(ilarge)     = n_tor*(index1-1) + k
              a_mat%val(ilarge)     = 1.d12
            enddo
          enddo
      
        endif
      enddo

    endif

    nz_AA = ilarge

    a_mat%ng  = n_AA
    a_mat%nnz = nz_AA

#ifdef USE_FFTW
    if (use_fft) call dfftw_destroy_plan(fftw_plan)
#endif

    end if
  call MPI_Bcast(a_mat%ng, 1, MPI_INTEGER, 0, this_mpi_comm_world, ierr)
  call MPI_Bcast(a_mat%nnz, 1, MPI_INTEGER, 0, this_mpi_comm_world, ierr)

  end subroutine assemble_projection_matrix


subroutine my_fft(in_fft,out_fft,n)

  implicit none

  real*8     :: in_fft(*)
  complex*16 :: out_fft(*)

  integer    :: i, n
  real*8     :: tmp_fft(2*n+2)

  tmp_fft(1:n) = in_fft(1:n)

  call RFT2(tmp_fft,n,1)

  do i=1,n
    out_fft(i) = cmplx(tmp_fft(2*i-1),tmp_fft(2*i))
  enddo

  return
end subroutine my_fft

end module mod_coupled_projection

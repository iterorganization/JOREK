!> Contains functions to calculate the fields created by the plasma current alone and by external currents
!!
!!  * The routines B_plasma and psi_plasma calculate the plasma field given a set of (R,Z) points
!!  * Routines to predict the best coil currents for given fixed-bnd equilibrium are also avaialable
!!  * Functions to calculate elliptic integrals are accessible as well
module mod_plasma_response
  
  use constants
  use mod_parameters
  use data_structure
  use gauss
  use basis_at_gaussian
  use tr_module
  use phys_module
  use mod_vtk
  
  implicit none
  
  private
  
  public ::  Greens_functions,  psi_plasma,   B_plasma,  comelp
  public ::  find_Icoils, find_Icoils_JET
  public ::  t_pol_coil, t_coil, psi_coils
  public ::  read_coils, construct_test_coil, destruct_coils, log_coils, log_coil
  public ::  plasma_fields_at_xyz
  
  ! --- Derived datatypes
  type :: t_pol_coil
    integer             :: n_fila
    real*8, allocatable :: R_fila(:)
    real*8, allocatable :: Z_fila(:)
    real*8, allocatable :: weight(:)
  end type t_pol_coil
  
  type :: t_coil
    character(len=80)      :: name
    integer                :: coil_type
    integer                :: n_turns
    type(t_pol_coil) :: pol_coil
  end type t_coil
  
  integer, parameter :: C_POL_COIL  = 1       !< Value of coil_type indicating a poloidal field coil
  logical            :: debug=.true.        !< Output debugging information (fort.xxx files)
  logical            :: debug_norm=.false.   !< Output normalized magnetic field vectors
  logical            :: verbose=.false.      !< Output verbose information

  contains
  
    
  !!-------------------------------------------------------------------
  !> Calculates Green's fuctions G(r,r') for B_R, B_Z and \psi
  !!
  !! The Green's functions are used to calculate the fields by
  !! integrating current densities over areas
  !!
  !!      \psi(r) = \int G_psi(r,r') J(r') dA'
  !!-------------------------------------------------------------------
  pure subroutine Greens_functions(R, Z, R_p, Z_p, G_BR, G_BZ, G_psi)
  
    implicit none
    
    ! --- Routine parameters
    real*8,   intent(in)    :: R      !< R
    real*8,   intent(in)    :: Z      !< Z 
    real*8,   intent(in)    :: R_p    !< R'
    real*8,   intent(in)    :: Z_p    !< Z'
    real*8,   intent(inout) :: G_BR   !< Green's function for BR
    real*8,   intent(inout) :: G_BZ   !< Green's function for BZ
    real*8,   intent(inout) :: G_psi  !< Green's function for psi
    
    ! --- Local variables
    real*8   :: rho2, kk, Kellip_kk, Eellip_kk   
    
    ! --- Reference : Simple Analytic Expressions for the Magnetic Field of a Circular Current Loop, NASA
    rho2  =  (R_p+R)**2.d0 + (Z_p-Z)**2.d0      
    kk    =  sqrt( 4.d0*R_p*R / rho2 ) 
 
    call comelp( kk, Kellip_kk, Eellip_kk)  !--- calculate elliptic functions
    
    G_BR  = (0.5d0/PI) / sqrt( rho2 ) * ( Z-Z_p ) / R                      &
          * ( (R_p**2 + R**2 + (Z_p-Z)**2) / ((R_p-R)**2.d0 + (Z_p-Z)**2.d0) * Eellip_kk - Kellip_kk )
    
    G_BZ  = (0.5d0/PI) / sqrt( rho2 )                                      &
          * ( (R_p**2 - R**2 - (Z_p-Z)**2) / ((R_p-R)**2.d0 + (Z_p-Z)**2.d0) * Eellip_kk + Kellip_kk )
    
    G_psi = (0.5d0/PI) * sqrt(R_p*R)/kk * ( (2.d0-kk**2.d0)*Kellip_kk - 2.d0*Eellip_kk )    
  end subroutine Greens_functions
  
  
  
  
         
  !------------------------------------------------------------------
  !> Calculates psi_plasma at given (R,Z) points 
  !------------------------------------------------------------------
  subroutine psi_plasma(node_list,element_list,R0, Z0, psi_p)

    implicit none

    type (type_node_list),    intent(in) :: node_list
    type (type_element_list), intent(in) :: element_list
    type (type_element)      :: element
    type (type_node)         :: nodes(n_vertex_max)
    
    real*8, intent(in)    :: R0(:), Z0(:)
    real*8, intent(inout) :: psi_p(:)

    real*8     :: x_g(n_gauss,n_gauss),        x_s(n_gauss,n_gauss),        x_t(n_gauss,n_gauss)
    real*8     :: y_g(n_gauss,n_gauss),        y_s(n_gauss,n_gauss),        y_t(n_gauss,n_gauss)
    real*8     :: eq_g(n_plane,n_gauss,n_gauss)
    
    ! --- local variables    
    integer    :: i, j, ms, mt, iv, inode, ife, mp, in
    integer    :: ierr, n_cpu, my_id, ife_delta, ife_min, ife_max, omp_nthreads, omp_tid
    real*8     :: zj0, R, Z, wst, xjac, delta_phi
    real*8     :: G_BR, G_BZ, G_psi
    integer    :: n_points
    
    n_points = size(R0,1)
    
    psi_p     = 0.d0
    delta_phi = 2.d0 * PI / float(n_plane) / float(n_period)
        
    !--- Go through all the elements
    do ife = 1, element_list%n_elements
    
      element = element_list%element(ife)

      do iv = 1, n_vertex_max
        inode     = element%vertex(iv)
        nodes(iv) = node_list%node(inode)
      enddo
      
      x_g(:,:)    = 0.d0; x_s(:,:)    = 0.d0; x_t(:,:)    = 0.d0;
      y_g(:,:)    = 0.d0; y_s(:,:)    = 0.d0; y_t(:,:)    = 0.d0;
      
      !--- Calculate R,Z and derivatives at gausstian points
      do i=1,n_vertex_max
        do j=1,n_degrees

          do ms=1, n_gauss
            do mt=1, n_gauss

              x_g(ms,mt) = x_g(ms,mt) + nodes(i)%x(1,j,1) * element%size(i,j) * H(i,j,ms,mt)
              y_g(ms,mt) = y_g(ms,mt) + nodes(i)%x(1,j,2) * element%size(i,j) * H(i,j,ms,mt)

              x_s(ms,mt) = x_s(ms,mt) + nodes(i)%x(1,j,1) * element%size(i,j) * H_s(i,j,ms,mt)
              x_t(ms,mt) = x_t(ms,mt) + nodes(i)%x(1,j,1) * element%size(i,j) * H_t(i,j,ms,mt)
              y_s(ms,mt) = y_s(ms,mt) + nodes(i)%x(1,j,2) * element%size(i,j) * H_s(i,j,ms,mt)
              y_t(ms,mt) = y_t(ms,mt) + nodes(i)%x(1,j,2) * element%size(i,j) * H_t(i,j,ms,mt)

            enddo
          enddo
        enddo
      enddo
      
      eq_g(:,:,:) = 0.d0
            
      !--- Calculate the current at gaussian points
      do i=1,n_vertex_max
        do j=1,n_degrees

          do mp=1,n_plane
            do ms=1, n_gauss
              do mt=1, n_gauss
                do in=1,n_tor
                  eq_g(mp,ms,mt) = eq_g(mp,ms,mt) + nodes(i)%values(in,j,3) * element%size(i,j) * H(i,j,ms,mt)  * HZ(in,mp)
                enddo !---ntor
      	      enddo !---gauss
            enddo !---gauss
          enddo !---planes

        enddo !---order
      enddo !---vertex
      
      !---Do gaussian and toroidal planes integration
      do i=1, n_points
        do mp=1,n_plane
          do ms=1, n_gauss
            do mt=1, n_gauss

              wst  = wgauss(ms)*wgauss(mt)
              xjac = x_s(ms,mt)*y_t(ms,mt) - x_t(ms,mt)*y_s(ms,mt)
              R    = x_g(ms,mt)
              Z    = y_g(ms,mt)
              
              !--- Calculate Green's function            
              call Greens_functions(R0(i), Z0(i), R, Z, G_BR, G_BZ, G_psi)
              zj0  = eq_g(mp,ms,mt)
            
              !--- psi = \int Greens_funct * J_phi * dA       see (4.66 Computational Methods in P.Physics, Jardin)
              psi_p(i) = psi_p(i) - zj0 / R * G_psi *xjac * wst * delta_phi * n_period/ (2.d0 * PI)
            
            enddo
          enddo
        enddo
      enddo
      
    
    enddo !---elements

  end subroutine psi_plasma
 
 



 
  !-------------------------------------------------------------------------------------
  !> Calculates the magnetic field produced by plasma currents at arbitrary x,y,z points 
  !-------------------------------------------------------------------------------------
  subroutine plasma_fields_at_xyz(my_id, node_list,element_list, x,y,z, bx, by, bz)

    !$ use omp_lib
    use mpi_mod

    implicit none

    integer,                  intent(in) :: my_id
    type (type_node_list),    intent(in) :: node_list
    type (type_element_list), intent(in) :: element_list
    real*8,  intent(in)                  :: x(:), y(:), z(:)     ! Points where fields are calculated
    real*8,  intent(inout)               :: bx(:), by(:), bz(:)

    ! --- local variables    
    type (type_element)      :: element
    type (type_node)         :: nodes(n_vertex_max)
    
    real*8     :: x_g(n_gauss,n_gauss),        x_s(n_gauss,n_gauss),        x_t(n_gauss,n_gauss)
    real*8     :: y_g(n_gauss,n_gauss),        y_s(n_gauss,n_gauss),        y_t(n_gauss,n_gauss)
    real*8     :: eq_g(n_plane,n_gauss,n_gauss)
    
    integer    :: i, j, ms, mt, iv, inode, ife, mp, in
    integer    :: ierr, n_cpu, ife_delta, ife_min, ife_max, omp_nthreads, omp_tid
    real*8     :: zj0, R, xp,yp,zp, dd, wst, xjac, delta_phi, phi
    real*8     :: d_vec(3), J_vec(3), cross(3), dB(3)
    real*8     :: wgauss_copy(n_gauss)
    integer    :: n_points

    real*8, allocatable :: bx_tmp(:), by_tmp(:), bz_tmp(:)

    ! --- MPI initialization
    call MPI_COMM_SIZE(MPI_COMM_WORLD, n_cpu, ierr) ! number of MPI procs
    n_cpu = max(n_cpu,1)

    ife_delta = ceiling(float(element_list%n_elements) / n_cpu)
    ife_min   =      my_id     * ife_delta + 1
    ife_max   = min((my_id +1) * ife_delta, element_list%n_elements)

   
    n_points = size(x,1)
    allocate(bx_tmp(n_points), by_tmp(n_points), bz_tmp(n_points))
    
    bx     = 0.d0;  by     = 0.d0;  bz     = 0.d0;
    bx_tmp = 0.d0;  by_tmp = 0.d0;  bz_tmp = 0.d0;

    delta_phi = 2.d0 * PI / float(n_plane) 
 
    wgauss_copy = wgauss

    ! --- OpenMP parallelization of element loop
    !$omp parallel default(none)                                                           &
    !$omp   shared(my_id,element_list,node_list, H, H_s, H_t, HZ, ife_min, ife_max,        &
    !$omp          delta_phi, n_points, x, y, z, bx_tmp, by_tmp, bz_tmp, wgauss_copy)      &
    !$omp   private(ife,iv,inode,element,nodes,i,j, in, mp, ms, mt,                        &
    !$omp           x_g, y_g, x_s, y_s, x_t, y_t, xjac, eq_g, zj0, R, xp, yp, zp, dd, phi, &
    !$omp           d_vec, J_vec, cross, dB, wst, omp_nthreads,omp_tid)
    
#ifdef OPENMP
    omp_nthreads = omp_get_num_threads()
    omp_tid      = omp_get_thread_num()
#else
    omp_nthreads = 1
    omp_tid      = 0
#endif
    
    !$omp do reduction(+:bx_tmp, by_tmp, bz_tmp)     
       
    !--- Go through all the elements
    do ife = ife_min, ife_max
    
      element = element_list%element(ife)

      do iv = 1, n_vertex_max
        inode     = element%vertex(iv)
        nodes(iv) = node_list%node(inode)
      enddo
      
      x_g(:,:)    = 0.d0; x_s(:,:)    = 0.d0; x_t(:,:)    = 0.d0;
      y_g(:,:)    = 0.d0; y_s(:,:)    = 0.d0; y_t(:,:)    = 0.d0;
      
      !--- Calculate R,Z and derivatives at gausstian points
      do i=1,n_vertex_max
        do j=1,n_degrees

          do ms=1, n_gauss
            do mt=1, n_gauss

              x_g(ms,mt) = x_g(ms,mt) + nodes(i)%x(1,j,1) * element%size(i,j) * H(i,j,ms,mt)
              y_g(ms,mt) = y_g(ms,mt) + nodes(i)%x(1,j,2) * element%size(i,j) * H(i,j,ms,mt)

              x_s(ms,mt) = x_s(ms,mt) + nodes(i)%x(1,j,1) * element%size(i,j) * H_s(i,j,ms,mt)
              x_t(ms,mt) = x_t(ms,mt) + nodes(i)%x(1,j,1) * element%size(i,j) * H_t(i,j,ms,mt)
              y_s(ms,mt) = y_s(ms,mt) + nodes(i)%x(1,j,2) * element%size(i,j) * H_s(i,j,ms,mt)
              y_t(ms,mt) = y_t(ms,mt) + nodes(i)%x(1,j,2) * element%size(i,j) * H_t(i,j,ms,mt)

            enddo
          enddo
        enddo
      enddo
      
      eq_g(:,:,:) = 0.d0
            
      !--- Calculate the current at gaussian points
      do i=1,n_vertex_max
        do j=1,n_degrees

          do mp=1,n_plane
            do ms=1, n_gauss
              do mt=1, n_gauss
                do in=1,n_tor
                  eq_g(mp,ms,mt) = eq_g(mp,ms,mt) + nodes(i)%values(in,j,3) * element%size(i,j) * H(i,j,ms,mt)  * HZ(in,mp)
                enddo !---ntor
      	      enddo !---gauss
            enddo !---gauss
          enddo !---planes

        enddo !---order
      enddo !---vertex
      
      !---Do gaussian and toroidal planes integration
      do ms=1, n_gauss
        do mt=1, n_gauss

          wst  = wgauss_copy(ms)*wgauss_copy(mt)
          xjac = x_s(ms,mt)*y_t(ms,mt) - x_t(ms,mt)*y_s(ms,mt)
          R    = x_g(ms,mt)
          zp   = y_g(ms,mt)

          do mp=1,n_plane

            phi   =  float(mp-1) * delta_phi
            xp    =  R * Cos(phi)
            yp    = -R * Sin(phi)
                                
            zj0   =  eq_g(mp,ms,mt)
            J_vec =  zj0/R*(/ Sin(phi), Cos(phi), 0.d0 /)       

            ! --- Go over the given points
            do i=1, n_points

              d_vec(:)  = (/ xp-x(i), yp-y(i), zp-z(i) /)

              dd        = max(sqrt( sum( d_vec(:)**2.d0 ) ) , 1.d-9 )
    
              cross     = (/  d_vec(2)*J_vec(3) - d_vec(3)*J_vec(2),  &
                              d_vec(3)*J_vec(1) - d_vec(1)*J_vec(3),  &
                              d_vec(1)*J_vec(2) - d_vec(2)*J_vec(1) /)
    
              dB(:)     =  cross(:) / (dd**3.d0) / (4.d0*PI) * wst * xjac * R * delta_phi
    
              bx_tmp(i) = bx_tmp(i) + dB(1)
              by_tmp(i) = by_tmp(i) + dB(2)
              bz_tmp(i) = bz_tmp(i) + dB(3)
            enddo

          enddo
        enddo
      enddo
      
    
    enddo !---elements
    !$omp end do
    !$omp end parallel


    call MPI_AllReduce(bx_tmp,bx,n_points,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    call MPI_AllReduce(by_tmp,by,n_points,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    call MPI_AllReduce(bz_tmp,bz,n_points,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

    deallocate(bx_tmp, by_tmp, bz_tmp)   

  end subroutine plasma_fields_at_xyz
  
  
  
  
 
 
  !------------------------------------------------------------------
  !> Calculates B_plasma at given (R,Z) points 
  !------------------------------------------------------------------
  subroutine B_plasma(node_list,element_list,R0, Z0, B_p)

    implicit none

    type (type_node_list),    intent(in) :: node_list
    type (type_element_list), intent(in) :: element_list
    type (type_element)      :: element
    type (type_node)         :: nodes(n_vertex_max)
    
    real*8, intent(inout)    :: R0(:), Z0(:)
    real*8, intent(inout) :: B_p(:,:)

    real*8     :: x_g(n_gauss,n_gauss),        x_s(n_gauss,n_gauss),        x_t(n_gauss,n_gauss)
    real*8     :: y_g(n_gauss,n_gauss),        y_s(n_gauss,n_gauss),        y_t(n_gauss,n_gauss)
    real*8     :: eq_g(n_plane,n_gauss,n_gauss)
    
    ! --- local variables    
    integer    :: i, j, ms, mt, iv, inode, ife, mp, in
    integer    :: ierr, n_cpu, my_id, ife_delta, ife_min, ife_max, omp_nthreads, omp_tid
    real*8     :: zj0, R, Z, wst, xjac, delta_phi
    real*8     :: G_BR, G_BZ, G_psi
    integer    :: n_points
    
    n_points = size(R0,1)
    
    B_p       = 0.d0
    delta_phi = 2.d0 * PI / float(n_plane) / float(n_period)
        
    !--- Go through all the elements
    do ife = 1, element_list%n_elements
    
      element = element_list%element(ife)

      do iv = 1, n_vertex_max
        inode     = element%vertex(iv)
        nodes(iv) = node_list%node(inode)
      enddo
      
      x_g(:,:)    = 0.d0; x_s(:,:)    = 0.d0; x_t(:,:)    = 0.d0;
      y_g(:,:)    = 0.d0; y_s(:,:)    = 0.d0; y_t(:,:)    = 0.d0;
      
      !--- Calculate R,Z and derivatives at gausstian points
      do i=1,n_vertex_max
        do j=1,n_degrees

          do ms=1, n_gauss
            do mt=1, n_gauss

              x_g(ms,mt) = x_g(ms,mt) + nodes(i)%x(1,j,1) * element%size(i,j) * H(i,j,ms,mt)
              y_g(ms,mt) = y_g(ms,mt) + nodes(i)%x(1,j,2) * element%size(i,j) * H(i,j,ms,mt)

              x_s(ms,mt) = x_s(ms,mt) + nodes(i)%x(1,j,1) * element%size(i,j) * H_s(i,j,ms,mt)
              x_t(ms,mt) = x_t(ms,mt) + nodes(i)%x(1,j,1) * element%size(i,j) * H_t(i,j,ms,mt)
              y_s(ms,mt) = y_s(ms,mt) + nodes(i)%x(1,j,2) * element%size(i,j) * H_s(i,j,ms,mt)
              y_t(ms,mt) = y_t(ms,mt) + nodes(i)%x(1,j,2) * element%size(i,j) * H_t(i,j,ms,mt)

            enddo
          enddo
        enddo
      enddo
      
      eq_g(:,:,:) = 0.d0
            
      !--- Calculate the current at gaussian points
      do i=1,n_vertex_max
        do j=1,n_degrees

          do mp=1,n_plane
            do ms=1, n_gauss
              do mt=1, n_gauss
                do in=1,n_tor
                  eq_g(mp,ms,mt) = eq_g(mp,ms,mt) + nodes(i)%values(in,j,3) * element%size(i,j) * H(i,j,ms,mt)  * HZ(in,mp)
                enddo !---ntor
      	      enddo !---gauss
            enddo !---gauss
          enddo !---planes

        enddo !---order
      enddo !---vertex
      
      !---Do gaussian and toroidal planes integration
      do i=1, n_points
        do mp=1,n_plane
          do ms=1, n_gauss
            do mt=1, n_gauss

              wst  = wgauss(ms)*wgauss(mt)
              xjac = x_s(ms,mt)*y_t(ms,mt) - x_t(ms,mt)*y_s(ms,mt)
              R    = x_g(ms,mt)
              Z    = y_g(ms,mt)
              
              !--- Calculate Green's functions (Simple Analytic Expressions for the Magnetic Field of a Circular Current Loop, NASA)      
              call Greens_functions(R0(i), Z0(i), R, Z, G_BR, G_BZ, G_psi)  
                      
              zj0  = eq_g(mp,ms,mt)
            
              !B       =           !j_phi   !green   !dA
              B_p(i,1) = B_p(i,1) + zj0 / R * G_BR * xjac * wst * delta_phi * n_period/ (2.d0 * PI)
              B_p(i,2) = B_p(i,2) + zj0 / R * G_BZ * xjac * wst * delta_phi * n_period/ (2.d0 * PI)     
            
            enddo
          enddo
        enddo
      enddo
      
    
    enddo !---elements

  end subroutine B_plasma
  
 

  
  
    
  !> Find the best coil currents for a given fixed boundary equilibrium (using Btan)
  !---------------------------------------------------------------------
  subroutine find_Icoils(node_list,element_list,bnd_node_list,bnd_elm_list, coils)
        
    use vacuum
    use vacuum_response
    use mod_interp
    use mod_basisfunctions
    
    implicit none
    
    type (type_node_list),       intent(in) :: node_list
    type (type_element_list),    intent(in) :: element_list
    type(type_bnd_node_list),    intent(in) :: bnd_node_list         !< List of boundary elements
    type(type_bnd_element_list), intent(in) :: bnd_elm_list          !< List of grid elements
    type(t_coil),                intent(in) :: coils(:)
 
    ! --- Local variables

    type(type_bnd_element) :: bndelem_m
    integer  :: l_starwall, l_tor
    integer  :: m_bndelem, m_pt, m_elm, mv1
    integer  :: i_vertex, i_dof, i_node, i_node_bnd, i_resp, i_resp_old, i_resp_0
    real*8   :: i_size, basfunc_i
    real*8   :: H1(2,n_degrees_1d), H1_s(2,n_degrees_1d), H1_ss(2,n_degrees_1d)
    real*8   :: P, P_s, P_t, P_st, P_ss, P_tt
    real*8   :: R, R_s, R_t, R_st, R_ss, R_tt, Z, Z_s, Z_t, Z_st, Z_ss, Z_tt
    real*8   :: s_pt, t_pt, s_or_t ! s and t values at current point
    real*8   :: xjac               ! 2D Jacobian
    real*8   :: B_pol(2)           ! Poloidal magnetic field
    real*8   :: e_par(2)           ! Vector tangential to interface
    real*8   :: P_R, P_Z           ! dPsi/dR, dPsi/dZ
    real*8   :: R1, R2, Z1, Z2
    logical  :: s_const            ! Is the bound. elem. an s=const side of the 2D element?
  
    real*8,  allocatable :: R_vec(:), Z_vec(:), coeff(:,:), B_all(:,:), B_p(:,:), B_ext(:,:), Psitot(:)
    real*8,  allocatable :: Btan_ext(:), v_tan(:,:), weights(:), psi_c(:), psi_c_min(:), psi_p(:)
    real*8,  allocatable :: A_mat_min(:,:), RHS_min(:), RHS_min_per_turn(:), pf_current_Aturn(:)
    integer, allocatable :: ipiv(:)
    integer              :: n_points, n_points_elm, numb_coils
    integer              :: i, pt, count, i_c, j_c, info
    real*8               :: Btan_c, Bmax, Bc(2)
    
    l_starwall=1;  l_tor=1;

    write(*,*) ' '
    write(*,*) '************************************************************'
    write(*,*) '******** Calculate best I_coils from fixed bnd *************'
    write(*,*) '************************************************************'
    write(*,*) ' '
        
    !---- Calculate coefficients for each coil and point  **coeff(point,coil) **
    !------------------------------------------------------------------------
    
    numb_coils   = size(coils,1)
    write(*,*) ' Coils number = ', numb_coils
    
    n_points_elm = 6
    n_points     = n_points_elm * bnd_elm_list%n_bnd_elements  
    
    allocate( R_vec(n_points), Z_vec(n_points), coeff(n_points, numb_coils) )
    allocate( B_all(n_points,2), B_ext(n_points,2), B_p(n_points,2) , Psitot(n_points))
    allocate( v_tan(n_points,2), Btan_ext(n_points), weights(n_points) )    
    R_vec = 0.d0 ; Z_vec = 0.d0; coeff = 0.d0; B_all=0.d0; B_ext=0.d0; B_p=0.d0;
    v_tan = 0.d0 ; weights=0.d0;
    
    write(*,*) ' n_bnd_elm    = ',  bnd_elm_list%n_bnd_elements
    write(*,*) ' n_points     = ',  n_points
    
    count = 0
  
    ! --- For every boundary element, do...
    L_MB: do m_bndelem = 1, bnd_elm_list%n_bnd_elements
    
      bndelem_m = bnd_elm_list%bnd_element(m_bndelem)
      m_elm     = bnd_elm_list%bnd_element(m_bndelem)%element
      mv1       = bnd_elm_list%bnd_element(m_bndelem)%side
    
      R1 = node_list%node(bndelem_m%vertex(1))%x(1,1,1)
      Z1 = node_list%node(bndelem_m%vertex(1))%x(1,1,2)
      R2 = node_list%node(bndelem_m%vertex(2))%x(1,1,1)
      Z2 = node_list%node(bndelem_m%vertex(2))%x(1,1,2)
    
      ! --- For several points in the boundary element, do...
      L_MP: do m_pt = 1, n_points_elm
    
        count = count + 1
        
        ! --- Determine 1D basis function (and derivatives) at current point
        s_or_t = float(m_pt-1)/float(n_points_elm-1)
    
        call basisfunctions1(s_or_t, H1, H1_s, H1_ss)
    
        ! --- Which s and t values correspond to the current point and is the
        !     boundary element an s=const or t=const side of the 2D element?
        select case (mv1)
        case (1)
          s_pt = s_or_t;  t_pt = 0.d0;    s_const = .false.
        case (2)
          s_pt = 1.d0;    t_pt = s_or_t;  s_const = .true.
        case (3)
          s_pt = s_or_t;  t_pt = 1.d0;    s_const = .false.
        case (4)
          s_pt = 0.d0;    t_pt = s_or_t;  s_const = .true.
        end select
    
        ! --- Determine coordinate values (plus derivatives)
        call interp_RZ(node_list, element_list, m_elm, s_pt, t_pt, R, R_s, R_t, Z, Z_s, Z_t)
    
        ! --- 2D Jacobian
        xjac = R_s * Z_t - R_t * Z_s
    
        ! --- Tangential vector to the interface
        if ( s_const ) then
          e_par = (/ R_t, Z_t /) / sqrt( R_t**2 + Z_t**2 ) * (R_t * (R2-R1) + Z_t * (Z2-Z1))/abs(R_t * (R2-R1) + Z_t * (Z2-Z1))
        else
          e_par = (/ R_s, Z_s /) / sqrt( R_s**2 + Z_s**2 ) * (R_s * (R2-R1) + Z_s * (Z2-Z1))/abs(R_s * (R2-R1) + Z_s * (Z2-Z1))
        end if
    
        ! --- Psi value (plus derivatives) at current point (l_tor mode)
        call interp(node_list, element_list, m_elm, 1, l_tor, s_pt, t_pt, P, P_s, P_t, P_st, P_ss, P_tt)
        Psitot(count) = P
        ! --- Poloidal magnetic field at current point
        P_R   = (   P_s * Z_t - P_t * Z_s ) / xjac ! dPsi/dR
        P_Z   = ( - P_s * R_t + P_t * R_s ) / xjac ! dPsi/dZ
        B_pol = (/ P_Z, -P_R /) / R
        
        do i_c = 1, numb_coils        
          call B_coil_unit(coils(i_c), R, Z, Bc) 
          coeff(count, i_c) = Bc(1)*e_par(1) +  Bc(2)*e_par(2)  
        enddo  
        
        R_vec(count)   = R;     Z_vec(count) = Z;       B_all(count,:) = B_pol(:);   v_tan(count,:) = e_par(:); 
    
      end do L_MP
    
    end do L_MB
    
    
    !------------------------- Calculate B external from the fixed boundary
    !-----------------------------------------------------------------
    write(*,*) ' '
    write(*,*) 'Calculating external field contribution of fixed bnd ...'
   
    call B_plasma(node_list,element_list, R_vec, Z_vec, B_p)
    
    B_ext(:,:)  = B_all(:,:) - B_p(:,:) 
    Btan_ext(:) = B_ext(:,1)*v_tan(:,1) +  B_ext(:,2)*v_tan(:,2)  
    allocate( psi_c(n_points), psi_c_min(n_points),psi_p(n_points)) 
    !--- Check coils initial guess
    open(25,file='B_ext.txt',status="replace", position="append", action="write")
    do i = 1, n_points
      write(25,'(5ES14.6)') R_vec(i), Z_vec(i), Btan_ext(i), sum(coeff(i,:) * pf_coils(:)%current)
    enddo
    close(25)
    
    Bmax = maxval(abs(Btan_ext))
    !--- Calculate weights for each point
    do i = 1, n_points
      if ( abs(Btan_ext(i)) < 0.1*Bmax) then 
        weights(i) = 10.d0
      else 
        weights(i) = Bmax / abs(Btan_ext(i))
      endif
    enddo
   
    !------------------------- Solve minimization problem -------------------
    !------------------------------------------------------------------------
    ! Here we minimize  \sum_i ( \B_ext - \B_coils  )_i**2
    ! calculate RHS and A matrix
    write(*,*) ' '
    write(*,*) 'Solve minimization problem'
    
    allocate(RHS_min(numb_coils), A_mat_min(numb_coils,numb_coils), RHS_min_per_turn(numb_coils) )
    allocate(pf_current_Aturn(numb_coils))
    allocate(ipiv(numb_coils))
    RHS_min = 0.d0;     A_mat_min = 0.d0
    
    do i_c = 1, numb_coils
      do pt = 1, n_points
        RHS_min(i_c) = RHS_min(i_c) + Btan_ext(pt) * coeff(pt, i_c) * weights(pt)              
      enddo
      
      do j_c = 1, numb_coils
        do pt = 1, n_points  
          A_mat_min(i_c, j_c) = A_mat_min(i_c, j_c) + coeff(pt, i_c) * coeff(pt, j_c) * weights(pt)     
        enddo
      enddo
    enddo

   call dgesv( numb_coils, 1, A_mat_min, numb_coils, ipiv, RHS_min, numb_coils, info )
   write(*,*) 'info = ', info
   do i = 1, numb_coils
     RHS_min_per_turn(i) = RHS_min(i)/coils(i)%n_turns
     pf_current_Aturn(i) = pf_coils(i)%current*coils(i)%n_turns
   end do

   write(*,*) ' '
   write(*,*) ' Found total coil currents (stored in Icoils_found.txt) '
   write(*,*) ' Positive currents follow the JOREK +phi direction      '
   write(*,*) ' '

   open(26,file='Icoils_found.txt',status="replace", position="append", action="write")
   open(25,file='Psi_coils.txt',status="replace", position="append", action="write")

   call  psi_coils(coils,R_vec,Z_vec,psi_c, pf_current_Aturn)
   call  psi_coils(coils,R_vec,Z_vec,psi_c_min,RHS_min)
   call  psi_plasma(node_list,element_list,R_vec,Z_vec,psi_p)

   write(25,*) 'R    Z    psi_coil    psi_coil_min    psi_plasma    psi_tot'
    do i = 1, n_points
      write(25,'(6ES14.6)') R_vec(i), Z_vec(i), psi_c(i), psi_c_min(i), psi_p(i), Psitot(i)
    enddo
   write(26,*) '#current Aturns      current[A]/turn'
   do i=1, numb_coils
     write(26,'(2ES18.10)')   RHS_min(i),  RHS_min_per_turn(i)
     write(*, '(2ES18.10)')   RHS_min(i),  RHS_min_per_turn(i)
   enddo
   write(26,'(A,99ES14.6)') '# pf_coils%current = ', RHS_min_per_turn
   close(26)
   close(25)
   !----- Compare given and calculated currents
   write(*,*) ' '
   write(*,*) ' Initial coil currents, Calculated currents, Relative differences '
   write(*,*) ' Assuming pf_coils%fcurrent is given in A/turn '
   do i=1, numb_coils
     write(*,'(3ES14.6)') pf_coils(i)%current,  RHS_min_per_turn(i), (pf_coils(i)%current - RHS_min(i)/coils(i)%n_turns) / (pf_coils(i)%current + 0.1)
   enddo
   deallocate(RHS_min_per_turn, RHS_min, A_mat_min, ipiv)
   deallocate(psi_c, psi_c_min, psi_p, pf_current_Aturn)
   deallocate( R_vec, Z_vec, coeff)
   deallocate( B_all, B_ext, B_p , Psitot)
   deallocate( v_tan, Btan_ext, weights)    
  
  end subroutine find_Icoils 
  
  
  
  
  
  
  !!-------------------------------------------------------------------------------------------
  !> Find the best coil currents for a given fixed boundary equilibrium (using Btan)  for JET
  !! This is a separate routine as the 20 JET coils are constrained to 10 circuits and therefore
  !! 10 degrees of freedom. TO DO: A generalization for circuits and coils dof should be done
  !!-------------------------------------------------------------------------------------------
  subroutine find_Icoils_JET(node_list,element_list,bnd_node_list,bnd_elm_list, coils)
    
    use vacuum
    use vacuum_response
    use mod_interp
    use mod_basisfunctions
    
    type (type_node_list),       intent(in) :: node_list
    type (type_element_list),    intent(in) :: element_list
    type(type_bnd_node_list),    intent(in) :: bnd_node_list         !< List of boundary elements
    type(type_bnd_element_list), intent(in) :: bnd_elm_list          !< List of grid elements
    type(t_coil),                intent(in) :: coils(:)

 
    ! --- Local variables

    type(type_bnd_element) :: bndelem_m
    integer  :: l_starwall, l_tor
    integer  :: m_bndelem, m_pt, m_elm, mv1
    integer  :: i_vertex, i_dof, i_node, i_node_bnd, i_resp, i_resp_old, i_resp_0
    real*8   :: i_size, basfunc_i
    real*8   :: H1(2,n_degrees_1d), H1_s(2,n_degrees_1d), H1_ss(2,n_degrees_1d)
    real*8   :: P, P_s, P_t, P_st, P_ss, P_tt
    real*8   :: R, R_s, R_t, Z, Z_s, Z_t
    real*8   :: s_pt, t_pt, s_or_t ! s and t values at current point
    real*8   :: xjac               ! 2D Jacobian
    real*8   :: B_pol(2)           ! Poloidal magnetic field
    real*8   :: e_par(2)           ! Vector tangential to interface
    real*8   :: P_R, P_Z           ! dPsi/dR, dPsi/dZ
    real*8   :: R1, R2, Z1, Z2
    logical  :: s_const            ! Is the bound. elem. an s=const side of the 2D element?
  
    real*8,  allocatable :: R_vec(:), Z_vec(:), coeff(:,:), B_all(:,:), B_p(:,:), B_ext(:,:)
    real*8,  allocatable :: Btan_ext(:), v_tan(:,:), weights(:)
    real*8,  allocatable :: A_mat_min(:,:), RHS_min(:)
    integer, allocatable :: ipiv(:), alpha(:,:)
    integer              :: n_points, n_points_elm, numb_coils, n_eq
    integer              :: i, pt, count, i_c, j_c, info, i_eq_i, i_eq_j
    real*8               :: Btan_c, Bmax, cont_ip, cont_jp, Bc(2)
  
    l_starwall=1;  l_tor=1;
  
    write(*,*) ' '
    write(*,*) '************************************************************'
    write(*,*) '****** Calculate best JET I_coils from fixed bnd ***********'
    write(*,*) '************************************************************'
    write(*,*) ' '
        
    !---- Calculate coefficients for each coil and point  **coeff(point,coil) **
    !------------------------------------------------------------------------
    
    numb_coils   = size(coils, 1) 
    
    if (numb_coils /= 20) then
      write(*,*) 'For using this feature with JET you must use the 20 coils generated by write_jet_pfsystems'
      stop
    endif
    write(*,*) ' Coils number = ', numb_coils
    
    n_points_elm = 6
    n_points     = n_points_elm * bnd_elm_list%n_bnd_elements  
    
    allocate( R_vec(n_points), Z_vec(n_points), coeff(n_points, numb_coils) )
    allocate( B_all(n_points,2), B_ext(n_points,2), B_p(n_points,2) )
    allocate( v_tan(n_points,2), Btan_ext(n_points), weights(n_points) )    
    R_vec = 0.d0 ; Z_vec = 0.d0; coeff = 0.d0; B_all=0.d0; B_ext=0.d0; B_p=0.d0;
    v_tan = 0.d0 ; weights=0.d0;
    
    write(*,*) ' n_bnd_elm    = ',  bnd_elm_list%n_bnd_elements
    write(*,*) ' n_points     = ',  n_points
    
    count = 0
  
      ! --- For every boundary element, do...
    L_MB: do m_bndelem = 1, bnd_elm_list%n_bnd_elements
    
      bndelem_m = bnd_elm_list%bnd_element(m_bndelem)
      m_elm     = bnd_elm_list%bnd_element(m_bndelem)%element
      mv1       = bnd_elm_list%bnd_element(m_bndelem)%side
    
      R1 = node_list%node(bndelem_m%vertex(1))%x(1,1,1)
      Z1 = node_list%node(bndelem_m%vertex(1))%x(1,1,2)
      R2 = node_list%node(bndelem_m%vertex(2))%x(1,1,1)
      Z2 = node_list%node(bndelem_m%vertex(2))%x(1,1,2)
    
      ! --- For several points in the boundary element, do...
      L_MP: do m_pt = 1, n_points_elm
    
        count = count + 1
        
        ! --- Determine 1D basis function (and derivatives) at current point
        s_or_t = float(m_pt-1)/float(n_points_elm-1)
    
        call basisfunctions1(s_or_t, H1, H1_s, H1_ss)
    
        ! --- Which s and t values correspond to the current point and is the
        !     boundary element an s=const or t=const side of the 2D element?
        select case (mv1)
        case (1)
          s_pt = s_or_t;  t_pt = 0.d0;    s_const = .false.
        case (2)
          s_pt = 1.d0;    t_pt = s_or_t;  s_const = .true.
        case (3)
          s_pt = s_or_t;  t_pt = 1.d0;    s_const = .false.
        case (4)
          s_pt = 0.d0;    t_pt = s_or_t;  s_const = .true.
        end select
    
        ! --- Determine coordinate values (plus derivatives)
        call interp_RZ(node_list, element_list, m_elm, s_pt, t_pt, R, R_s, R_t, Z, Z_s, Z_t)
    
        ! --- 2D Jacobian
        xjac = R_s * Z_t - R_t * Z_s
    
        ! --- Tangential vector to the interface
        if ( s_const ) then
          e_par = (/ R_t, Z_t /) / sqrt( R_t**2 + Z_t**2 ) * (R_t * (R2-R1) + Z_t * (Z2-Z1))/abs(R_t * (R2-R1) + Z_t * (Z2-Z1))
        else
          e_par = (/ R_s, Z_s /) / sqrt( R_s**2 + Z_s**2 ) * (R_s * (R2-R1) + Z_s * (Z2-Z1))/abs(R_s * (R2-R1) + Z_s * (Z2-Z1))
        end if
    
        ! --- Psi value (plus derivatives) at current point (l_tor mode)
        call interp(node_list, element_list, m_elm, 1, l_tor, s_pt, t_pt, P, P_s, P_t, P_st, P_ss, P_tt)
    
        ! --- Poloidal magnetic field at current point
        P_R   = (   P_s * Z_t - P_t * Z_s ) / xjac ! dPsi/dR
        P_Z   = ( - P_s * R_t + P_t * R_s ) / xjac ! dPsi/dZ
        B_pol = (/ P_Z, -P_R /) / R
        
        do i_c = 1, numb_coils        
          call B_coil_unit(coils(i_c), R, Z, Bc) 
          coeff(count, i_c) = Bc(1)*e_par(1) +  Bc(2)*e_par(2)  
        enddo  
      
        R_vec(count)   = R;     Z_vec(count) = Z;       B_all(count,:) = B_pol(:);   v_tan(count,:) = e_par(:); 
    
      end do L_MP
    
    end do L_MB
    
    
    !------------------------- Calculate B external from the fixed boundary
    !-----------------------------------------------------------------
    write(*,*) ' '
    write(*,*) 'Calculating external field contribution of fixed bnd ...'
   
    call B_plasma(node_list,element_list, R_vec, Z_vec, B_p)
    
    B_ext(:,:)  = B_all(:,:) - B_p(:,:) 
    Btan_ext(:) = B_ext(:,1)*v_tan(:,1) +  B_ext(:,2)*v_tan(:,2)  
    
    Bmax = maxval(abs(Btan_ext))
    !--- Calculate weights for each point
    do i = 1, n_points
      if ( abs(Btan_ext(i)) < 0.1*Bmax) then 
        weights(i) = 10.d0
      else 
        weights(i) = Bmax / abs(Btan_ext(i))
      endif
    enddo
   
    !------------------------- Solve minimization problem -------------------
    !------------------------------------------------------------------------
    ! Here we minimize  \sum_i ( \B_ext - \B_coils  )_i**2
    ! calculate RHS and A matrix
    write(*,*) ' '
    write(*,*) 'Solve minimization problem'
    
    n_eq = 10
    
    allocate(RHS_min(n_eq), A_mat_min(n_eq,n_eq) )
    allocate(ipiv(n_eq))
    allocate(alpha(n_eq,numb_coils))
    RHS_min = 0.d0;     A_mat_min = 0.d0
    
    !--- Define relations between circuits and coils, multiply by number of turns
    alpha = 0.d0   !first index = circuit,  second = coil
    alpha(1, 1)   = 710.d0
    alpha(2, 2)   = 426.d0
    alpha(3, 3)   =  -8.d0
    alpha(3, 4)   = -20.d0
    alpha(3, 5)   =  -8.d0
    alpha(3, 6)   = -20.d0
    alpha(4, 7)   =  -8.d0
    alpha(4, 8)   =   8.d0
    alpha(3, 9)   =  30.d0
    alpha(3,10)   =  30.d0
    alpha(4,11)   = -20.d0
    alpha(4,12)   =  20.d0
    alpha(1,13)   =   2.d0
    alpha(1,14)   =   2.d0
    alpha(5,15)   =  61.d0
    alpha(5,16)   =  61.d0
    alpha(6,15)   = -61.d0
    alpha(6,16)   =  61.d0
    alpha(7,17)   =  15.99d0
    alpha(8,18)   =  15.d0
    alpha(9,19)   =  15.d0
    alpha(10,20)  =  21.d0
    
    !--- Calculate B_j 
    do i_eq_j = 1, n_eq
      do pt = 1, n_points
        do i_c = 1, numb_coils        
          RHS_min(i_eq_j) = RHS_min(i_eq_j) + Btan_ext(pt) * coeff(pt, i_c) * weights(pt) * alpha(i_eq_j, i_c)             
        enddo
      enddo
    enddo
    
    !--- Calculate A_ij
    do i_eq_i = 1, n_eq          ! loops for A matrix main indeces
      do i_eq_j = 1, n_eq   
        
        do pt = 1, n_points 
          
          cont_ip = 0.d0; cont_jp = 0.d0
          
          do i_c = 1, numb_coils           
            cont_ip = cont_ip + alpha(i_eq_i, i_c) * coeff(pt, i_c) 
            cont_jp = cont_jp + alpha(i_eq_j, i_c) * coeff(pt, i_c)  
          enddo
          
          A_mat_min(i_eq_i, i_eq_j) = A_mat_min(i_eq_i, i_eq_j) + cont_ip * cont_jp * weights(pt)
        
        enddo
        
      enddo
    enddo

   ! --- Solve linear system of equations for minimization 
   call dgesv( n_eq, 1, A_mat_min, n_eq, ipiv, RHS_min, n_eq, info )
   write(*,*) 'info = ', info
   
   ! --- Write out found currents
   write(*,*) ' '
   write(*,*) ' Found circuit currents (stored in Icircuits_found.txt) '

   open(25,file='Icircuits_found.txt',status="replace", position="append", action="write")
   do i=1, n_eq
     write(25,'(1ES18.10)') RHS_min(i)
     write(*, '(1ES18.10)') RHS_min(i)
   enddo
   close(25)
      
   write(*,*) ' '
   write(*,*) ' Found total coil currents (stored in Icoils_found.txt) '
   write(*,*) ' Positive currents follow the JOREK +phi direction      '
   write(*,*) ' '

   open(26,file='Icoils_found.txt',status="replace", position="append", action="write")
   do i=1, numb_coils
     write(26,'(1ES18.10)')sum( alpha(:,i) * RHS_min(:) )
     write(*, '(1ES18.10)')sum( alpha(:,i) * RHS_min(:) )
   enddo
   close(26)
  
   !--- Check coils initial guess
   open(27,file='B_ext_new.txt',status="replace", position="append", action="write")
   do i = 1, n_points
     write(27,'(4ES14.6)') R_vec(i), Z_vec(i), Btan_ext(i),  sum(coeff(i,:) * pf_coils(:)%current) 
   enddo
   close(27)

  end subroutine find_Icoils_JET
  
  
  
  
  
  
  
  ! --- Routine to calculate the elliptic integrals
  pure subroutine comelp ( hk, ck, ce )
  
    !*****************************************************************************80
    !
    !! COMELP computes complete elliptic integrals K(k) and E(k).
    !
    !  Licensing:
    !
    !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
    !    they give permission to incorporate this routine into a user program 
    !    provided that the copyright is acknowledged.
    !
    !  Modified:
    !
    !    07 July 2012
    !
    !  Author:
    !
    !    Shanjie Zhang, Jianming Jin
    !
    !  Reference:
    !
    !    Shanjie Zhang, Jianming Jin,
    !    Computation of Special Functions,
    !    Wiley, 1996,
    !    ISBN: 0-471-11963-6,
    !    LC: QA351.C45.
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) HK, the modulus.  0 <= HK <= 1.
    !
    !    Output, real ( kind = 8 ) CK, CE, the values of K(HK) and E(HK).
    !
      implicit none
    
      real*8, intent(in)    :: hk
      real*8, intent(inout) :: ck
      real*8, intent(inout) :: ce
      
      real ( kind = 8 ) ae
      real ( kind = 8 ) ak
      real ( kind = 8 ) be
      real ( kind = 8 ) bk
      real ( kind = 8 ) pk
    
      pk = 1.0D+00 - hk * hk
    
      if ( hk == 1.0D+00 ) then
    
        ck = 1.0D+300
        ce = 1.0D+00
    
      else
    
        ak = ((( &
            0.01451196212D+00   * pk &
          + 0.03742563713D+00 ) * pk &
          + 0.03590092383D+00 ) * pk &
          + 0.09666344259D+00 ) * pk &
          + 1.38629436112D+00
    
        bk = ((( &
            0.00441787012D+00   * pk &
          + 0.03328355346D+00 ) * pk &
          + 0.06880248576D+00 ) * pk &
          + 0.12498593597D+00 ) * pk &
          + 0.5D+00
    
        ck = ak - bk * log ( pk )
    
        ae = ((( &
            0.01736506451D+00   * pk &
          + 0.04757383546D+00 ) * pk &
          + 0.0626060122D+00  ) * pk &
          + 0.44325141463D+00 ) * pk &
          + 1.0D+00
    
        be = ((( &
            0.00526449639D+00   * pk &
          + 0.04069697526D+00 ) * pk &
          + 0.09200180037D+00 ) * pk &
          + 0.2499836831D+00  ) * pk
    
        ce = ae - be * log ( pk )
    
      end if
    
      return
  end  subroutine comelp






! --- Routines related to coils

  subroutine read_coils(coils, coil_file)
    implicit none
    
    ! --- Routine parameters
    type(t_coil), allocatable, intent(inout) :: coils(:)
    character(len=*),          intent(in)    :: coil_file
    ! --- Local variables
    integer :: file_handle = 21
    integer :: err, ncoils, i_c, i_f
    real*8, allocatable :: n_turns(:)
    real*8  :: n_turns_fila
    character(len=512) :: comment_line
    
    call destruct_coils(coils)
    open(file_handle, file=trim(coil_file), status='old', action='read', iostat=err)
    if ( err /= 0 ) then
      write(*,*) 'ERROR: Cannot open coil file "'//trim(coil_file)//'".'
      stop
    end if
    
    read(file_handle,'(a)') comment_line
    read(file_handle,*) ncoils
    
    allocate( coils(ncoils), n_turns(ncoils) )
    
    read(file_handle,'(a)') comment_line
    do i_c = 1, ncoils
      read(file_handle,*) coils(i_c)%pol_coil%n_fila, coils(i_c)%name, n_turns(i_c)
      
      coils(i_c)%coil_type = C_POL_COIL
      coils(i_c)%name = trim(adjustl(coils(i_c)%name))
      coils(i_c)%n_turns  =  n_turns(i_c)
      allocate( coils(i_c)%pol_coil%R_fila(coils(i_c)%pol_coil%n_fila) )
      allocate( coils(i_c)%pol_coil%Z_fila(coils(i_c)%pol_coil%n_fila) )
      allocate( coils(i_c)%pol_coil%weight(coils(i_c)%pol_coil%n_fila) )
    end do
    
    do i_c = 1, ncoils
      read(file_handle,'(a)') comment_line
      read(file_handle,'(a)') comment_line
      read(file_handle,'(a)') comment_line
      do i_f = 1, coils(i_c)%pol_coil%n_fila
        read(file_handle,*) coils(i_c)%pol_coil%R_fila(i_f), coils(i_c)%pol_coil%Z_fila(i_f),      &
          n_turns_fila
        
        coils(i_c)%pol_coil%weight(i_f) = n_turns_fila / n_turns(i_c)
      end do
    end do
    
    close(file_handle)
    
  end subroutine read_coils
  
  
  
  subroutine construct_test_coil(coils)
    implicit none
    
    ! --- Routine parameters
    type(t_coil), allocatable, intent(inout) :: coils(:)
    
    call destruct_coils(coils)
    allocate( coils(1) )
    coils(1)%name = 'Axisymmetric Test-Coil'
    coils(1)%coil_type = C_POL_COIL
    coils(1)%pol_coil%n_fila = 1
    allocate( coils(1)%pol_coil%R_fila(1) )
    allocate( coils(1)%pol_coil%Z_fila(1) )
    allocate( coils(1)%pol_coil%weight(1) )
    coils(1)%pol_coil%R_fila(:) = (/ 9.2 /)
    coils(1)%pol_coil%Z_fila(:) = (/ 1.0 /)
    coils(1)%pol_coil%weight(:) = (/ 1.0 /)
  end subroutine construct_test_coil
  
  
  
  subroutine destruct_coils(coils)
    implicit none
    
    ! --- Routine parameters
    type(t_coil), allocatable, intent(inout) :: coils(:)
    ! --- Local variables
    integer :: i_c
    
    if ( allocated(coils) ) then
      do i_c = 1, size(coils,1)
        if ( allocated(coils(i_c)%pol_coil%R_fila) ) deallocate(coils(i_c)%pol_coil%R_fila)
        if ( allocated(coils(i_c)%pol_coil%Z_fila) ) deallocate(coils(i_c)%pol_coil%Z_fila)
        if ( allocated(coils(i_c)%pol_coil%weight) ) deallocate(coils(i_c)%pol_coil%weight)
      end do
      deallocate(coils)
    end if
  end subroutine destruct_coils
  
  
  
  subroutine log_coils(coils,verbose)
    implicit none
    
    ! --- Routine parameters
    type(t_coil), allocatable, intent(in) :: coils(:)
    logical,                   intent(in) :: verbose
    ! --- Local variables
    integer :: i
    
    write(*,*)
    write(*,*) '+-- COILS --------------------------------------------------------'
    if ( allocated(coils) ) then
      write(*,'(1x,a,i4)') '| Number of coils:', size(coils,1)
      do i = 1, size(coils,1)
        call log_coil(coils(i),verbose)
      end do
    else
      write(*,*) '| [not allocated]'
    end if
    write(*,*) '+-----------------------------------------------------------------'
    write(*,*)
  end subroutine log_coils
  
  
  
  subroutine log_coil(coil,verbose)
    implicit none
    
    ! --- Routine parameters
    type(t_coil), intent(in) :: coil
    logical,      intent(in) :: verbose
    ! --- Local variables
    integer :: i
    
    write(*,'(1x,5a)') '| Coil "', trim(coil%name), '" (', trim(coil_type_name(coil%coil_type)), ')'
    if ( verbose .or. debug ) then
      if ( coil%coil_type == C_POL_COIL ) then
        write(*,*) '|   Filament         R               Z            weight'
        do i = 1, coil%pol_coil%n_fila
          write(*,'(1x,a,i10,3es16.7)') '|', i, coil%pol_coil%R_fila(i), coil%pol_coil%Z_fila(i),   &
            coil%pol_coil%weight(i)
          if ( debug ) write(37,*) coil%pol_coil%R_fila(i), coil%pol_coil%Z_fila(i)
        end do
        if ( debug ) write(37,*)
      else
        write(*,*) '|   ERROR: ILLEGAL COIL_TYPE!'
      end if
    end if
    
  end subroutine log_coil
 




  character(len=128) function coil_type_name(coil_type)
    implicit none
    
    ! --- Routine parameters
    integer, intent(in) :: coil_type
    
    if ( coil_type == C_POL_COIL ) then
      coil_type_name = 'poloidal field coil'
    else
      coil_type_name = 'ILLEGAL COIL TYPE'
    end if
    
  end function coil_type_name
 




  !------------------------------------------------------------------
  !> Calculates psi_coils at given (R,Z) points 
  !------------------------------------------------------------------
  subroutine psi_coils(coils, R0, Z0, psi_c, current_in)
    
    use vacuum, only : pf_coils

    implicit none

    type(t_coil),  intent(in) :: coils(:)

    real*8, intent(in)        :: R0(:), Z0(:)
    real*8, intent(inout)     :: psi_c(:)
    real*8, optional,intent(in)    :: current_in(:)
    
    ! --- local variables    
    integer    :: i_p, i_c, i_f 
    real*8     :: R_f, Z_f, I_coil
    real*8     :: G_BR, G_BZ, G_psi
    integer    :: n_points, n_coils
    
    n_points = size(R0,1)
    n_coils  = size(coils,1)    

    psi_c    = 0.d0
        
    do i_p = 1, n_points   ! --- loop over given RZ points 
      do i_c =1, n_coils     ! --- loop over coils

        I_coil = pf_coils(i_c)%current
        if (present(current_in)) then
          I_coil = current_in(i_c)
        end if
   
        do i_f = 1, coils(i_c)%pol_coil%n_fila ! (Loop over coil filaments)

          ! --- Position of filament point
          R_f   = coils(i_c)%pol_coil%R_fila(i_f)
          Z_f   = coils(i_c)%pol_coil%Z_fila(i_f)

           !--- Calculate Green's function            
          call Greens_functions(R0(i_p), Z0(i_p), R_f, Z_f, G_BR, G_BZ, G_psi)
            
          !--- psi = \int Greens_funct * I   see (4.66 Computational Methods in P.Physics, Jardin)
          psi_c(i_p) = psi_c(i_p) + G_psi * I_coil * coils(i_c)%pol_coil%weight(i_f) * mu_zero

        enddo
      enddo
    enddo
      
  end subroutine psi_coils
 
 
 

  !------------------------------------------------------------------
  !> Calculates B-field of a single coil with unitary current
  !------------------------------------------------------------------
  subroutine B_coil_unit(coil, R0, Z0, B)

    implicit none

    type(t_coil),  intent(in) :: coil

    real*8, intent(in)        :: R0, Z0
    real*8, intent(inout)     :: B(2)
    
    ! --- local variables    
    integer    :: i_p, i_c, i_f 
    real*8     :: R_f, Z_f
    real*8     :: G_BR, G_BZ, G_psi
    
    B = 0.d0     

    do i_f = 1, coil%pol_coil%n_fila ! (Loop over coil filaments)

      ! --- Position of filament point
      R_f   = coil%pol_coil%R_fila(i_f)
      Z_f   = coil%pol_coil%Z_fila(i_f)

      !--- Calculate Green's function            
      call Greens_functions(R0, Z0, R_f, Z_f, G_BR, G_BZ, G_psi)
            
      !--- B = \int Greens_funct * I   see (4.66 Computational Methods in P.Physics, Jardin)
      B(1)  = B(1)  -  G_BR *  coil%pol_coil%weight(i_f) * mu_zero
      B(2)  = B(2)  -  G_Bz *  coil%pol_coil%weight(i_f) * mu_zero

    enddo
      
  end subroutine B_coil_unit




   
end module mod_plasma_response

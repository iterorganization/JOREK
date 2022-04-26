!> This module contains some helper functions for particle projection
!> These are a few reference functions f_*, a function to project one of these,
!> functions to generate a square, polar and flux-aligned grid,
!> functions to calculate the projection RHS and a function to calculate the mean
!> and RMS error of a projected function.
module projection_helpers
use mod_project_particles
use data_structure
use mod_particle_types
implicit none

interface
  subroutine grid_flux_surface(xpoint,xcase,node_list,element_list,surface_list,n_flux,n_tht,xr1,sig1,xr2,sig2,refinement)
    use data_structure
    ! --- Routine parameters
    logical,                  intent(in)    :: xpoint, refinement
    type (type_node_list),    intent(inout) :: node_list
    type (type_element_list), intent(inout) :: element_list
    type (type_surface_list), intent(inout) :: surface_list
    integer,                  intent(in)    :: xcase
    integer,                  intent(in)    :: n_flux
    integer,                  intent(in)    :: n_tht
    real*8,                   intent(in)    :: xr1, xr2
    real*8,                   intent(in)    :: sig1, sig2
  end subroutine grid_flux_surface
end interface

contains

!> Functions to project
function f_0(R, Z)
  real*8, intent(in) :: R, Z
  real*8 :: f_0
  f_0 = 0.d0
end function f_0
function f_1(R, Z)
  real*8, intent(in) :: R, Z
  real*8 :: f_1
  f_1 = 1.d0
end function f_1
function f_R(R, Z)
  real*8, intent(in) :: R, Z
  real*8 :: f_R
  f_R = R
end function f_R
function f_RZ(R, Z)
  real*8, intent(in) :: R, Z
  real*8 :: f_RZ
  f_RZ = R*Z
end function f_RZ
function f_R4(R, Z)
  real*8, intent(in) :: R, Z
  real*8 :: f_R4
  f_R4 = R**4
end function f_R4
!> minor radius^2
function f_a2(R, Z)
  use phys_module, only: R_geo, Z_geo
  real*8, intent(in) :: R, Z
  real*8 :: f_a2
  f_a2 = norm2([(R-R_geo),Z-Z_geo])**1.8d0
end function f_a2


!> Create a simple square grid with n nodes in each dimension
subroutine default_square_grid(node_list, element_list, n)
  type(type_node_list), intent(out) :: node_list
  type(type_element_list), intent(out) :: element_list
  integer, intent(in) :: n !< number of nodes in each dimension
  real*8, parameter :: R_geo = 1.d0, Z_geo = 0.d0, amin = 0.5d0
  node_list%n_nodes = 0
  element_list%n_elements = 0
  call grid_bezier_square(n, n, R_geo-amin,R_geo+amin, Z_geo-amin, Z_geo+amin, .true., node_list, element_list)
  write(*,*) ' completed default_square_grid'
end subroutine

!> Create a simple polar grid with npol nodes in the poloidal direction, 30 radial
!> volume = 2 pi^2 R a^2
subroutine default_polar_grid(node_list, element_list, npol)
  use phys_module
  type(type_node_list), intent(out) :: node_list
  type(type_element_list), intent(out) :: element_list
  integer, intent(in) :: npol !< number of nodes in each dimension
  call preset_parameters()
  fbnd(1) = 2.d0
  fbnd(2:4) = 0.d0
  fpsi = 0.d0
  mf = 0
  n_radial = 30
  R_geo = 1.5
  Z_geo = 0.0
  amin = 1.0

  node_list%n_nodes = 0
  element_list%n_elements = 0
  call grid_polar_bezier(R_geo, Z_geo, amin, 0.d0, 0.d0, fbnd, fpsi, mf, n_radial, npol,    &
    node_list, element_list)
end subroutine default_polar_grid

subroutine default_flux_grid_31(node_list,element_list)
  type(type_node_list), intent(out) :: node_list
  type(type_element_list), intent(out) :: element_list
  type(type_node_list), allocatable, save :: cache_node_list
  type(type_element_list), allocatable, save :: cache_element_list
  logical, save :: saved = .false.
  if (.not. saved) then
    allocate(cache_node_list,cache_element_list)
    call default_flux_grid(cache_node_list, cache_element_list, 31)
    saved = .true.
  end if
  node_list = cache_node_list
  element_list = cache_element_list
end subroutine default_flux_grid_31

subroutine default_flux_grid_32(node_list,element_list)
  type(type_node_list), intent(out) :: node_list
  type(type_element_list), intent(out) :: element_list
  type(type_node_list), allocatable, save :: cache_node_list
  type(type_element_list), allocatable, save :: cache_element_list
  logical, save :: saved = .false.
  if (.not. saved) then
    allocate(cache_node_list,cache_element_list)
    call default_flux_grid(cache_node_list, cache_element_list, 32)
    saved = .true.
  end if
  node_list = cache_node_list
  element_list = cache_element_list
end subroutine default_flux_grid_32


!> Create a simple flux aligned grid with npol nodes in the poloidal direction, 40 radial
!> by calculating equilibrium and creating flux aligned grid (like in jorek2_main)
subroutine default_flux_grid(node_list, element_list, npol)
  use phys_module
  use mod_boundary
  use mpi_mod
  use mod_export_restart
  type(type_node_list), intent(out) :: node_list
  type(type_element_list), intent(out) :: element_list
  integer, intent(in) :: npol !< number of nodes in each dimension
  type (type_surface_list) :: surface_list
  type(type_bnd_node_list) :: bnd_node_list
  type(type_bnd_element_list) :: bnd_elm_list

  call tr_resetfile()
  call det_modes()
  ! Start with a polar grid
  call preset_parameters()
  ! Copy input file for simple case
  F0 = 3.0
  fbnd(1) = 2.d0
  fbnd(2:4) = 0.d0
  fpsi = 0.d0
  mf = 0
  n_radial = 30
  n_pol = 32
  R_geo = 1.5
  Z_geo = 0.0
  amin = 1.0
  ellip  = 1.
  tria_u = 0.
  tria_l = 0.
  quad_u = -0.
  quad_l = -0.
  xpoint = .false.

  rho_0 = 1.d0
  rho_1 = 1.d-1
  rho_coef(1)  =  -1.0d0
  rho_coef(2)  =  0.0
  rho_coef(3)  =  0.0
  rho_coef(4)  =  1.d0
  rho_coef(5)  =  5.d0
  T_0 = 2.d-3
  T_1 = 1.d-8
  T_coef(1) = -0.8d0
  T_coef(2) = +0.0d0
  T_coef(3) = 0.d0
  T_coef(4) = 1.d0
  T_coef(5) = 5.d0
  FF_0 = 2.d0
  FF_1 = 0.
  FF_coef(1)  = -1.0
  FF_coef(2)  = 0.d0
  FF_coef(3)  = 0.d0
  FF_coef(4)  = 0.03d0
  FF_coef(5)  = 5.0d0
  FF_coef(6)  = 1.d0
  FF_coef(7)  = 10.d0
  FF_coef(8)  = 1.d0
  D_par  = 0.d0
  D_perp = 1.d-5

  ZK_par  = 1.d0
  ZK_perp = 1.d-5
  
  eta   = 1.d-5
  visco = 1.d-6
  visco_par = 1.d-6

  eta_T_dependent = .false.
  visco_T_dependent = .false.

  heatsource     = 0.d0
  particlesource = 0.d0
  rst_hdf5 = 1


  node_list%n_nodes = 0
  element_list%n_elements = 0
  bnd_elm_list%n_bnd_elements  = 0

  use_mumps_eq  = .true.
  use_pastix_eq = .false.

  call define_boundary()

  call grid_polar_bezier(R_geo, Z_geo, amin, 0.d0, 0.d0, fbnd, fpsi, mf, n_radial, n_pol,    &
    node_list, element_list)

  call boundary_from_grid(node_list, element_list, bnd_node_list, bnd_elm_list, .false.)

  call initialise_mumps(MPI_COMM_WORLD)

  call equilibrium(0,node_list,element_list,bnd_node_list,bnd_elm_list,.false., 0, .true.) 

  ! Set parameters for making flux grid
  call grid_flux_surface(.false., 0, node_list, element_list, &
      surface_list, n_flux=30, n_tht=npol, xr1=xr1, sig1=sig1, xr2=xr2, sig2=sig2, refinement=.false.)
end subroutine default_flux_grid



!> Project a function onto the JOREK elements
subroutine project_f(node_list, element_list, f, filter, filter_hyper, integral)
  use mpi_mod
  type(type_node_list), intent(inout)    :: node_list
  type(type_element_list), intent(inout) :: element_list
  real*8, external                       :: f
  real*8, optional, intent(in)           :: filter, filter_hyper !< smoothing and hyper-smoothing
  real*8, optional, intent(out)          :: integral !< The integral of the projected function, from the weights
  real*8, dimension(:), allocatable      :: this_integral_weights
  type(DMUMPS_STRUC)                     :: p
  integer :: i, k, index, i_tor_local, n_tor_local, mpi_comm_n, mpi_comm_master, ierr
  real*8  :: my_filter, my_filter_hyper, area, volume

  call MPI_Comm_dup(MPI_COMM_WORLD, mpi_comm_n, ierr)
  call MPI_Comm_dup(MPI_COMM_WORLD, mpi_comm_master, ierr)

  i_tor_local     = 1
  n_tor_local     = 1

  my_filter       = 0.d0
  my_filter_hyper = 0.d0

  if (present(filter))       my_filter       = filter 
  if (present(filter_hyper)) my_filter_hyper = filter_hyper 
  
!  if (present(integral)) then  
  if (i_tor_local .eq. 1) then  
    call prepare_mumps_par_n0(node_list, element_list, n_tor_local, i_tor_local, mpi_comm_world, mpi_comm_n, mpi_comm_master, &
                              p,  area, volume, filter=my_filter, filter_hyper=my_filter_hyper, filter_parallel=0.d0, integral_weights=this_integral_weights)
  else
    call prepare_mumps_par(node_list, element_list, n_tor_local, i_tor_local, mpi_comm_world, mpi_comm_n, mpi_comm_master, &
                           p, filter=my_filter, filter_hyper=my_filter_hyper, filter_parallel=0.d0)
  endif

  ! Project manually
  p%JOB = 3
  p%icntl(21) = 0 ! solution is available only on host
  p%icntl(4)  = 1 ! print only errors

  ! Setup RHS by integrating manually
  allocate(p%rhs(p%n))
  call calc_rhs_f(node_list,element_list,f,p%rhs)

  call DMUMPS(p)

  if (present(integral)) integral = dot_product(p%rhs, this_integral_weights)

  do i=1,node_list%n_nodes
    do k=1,n_degrees
      index = 2*(node_list%node(i)%index(k)-1) + 1
      node_list%node(i)%values(1,k,1) = p%rhs(index)
    enddo
  enddo

  p%JOB=-2
  call DMUMPS(p)
end subroutine project_f

!> Calculate the right-hand side of a distribution f (toroidally symmetric)
subroutine calc_rhs_f(node_list,element_list,f,rhs)
  use basis_at_gaussian
  use phys_module, only : TWOPI
  type(type_node_list),    intent(inout) :: node_list
  type(type_element_list), intent(inout) :: element_list
  real*8, external   :: f
  real*8, dimension(:), intent(inout) :: rhs
  integer            :: i, j, k, m, index, i_elm, inode, ms, mt
  real*8, dimension(n_gauss,n_gauss) :: x_g, y_g, x_s, x_t, y_s, y_t
  real*8             :: wst, xjac, v
  type(type_node)    :: nodes(4)
  type(type_element) :: element

  rhs = 0.d0

  do i_elm=1,element_list%n_elements

    element = element_list%element(i_elm)
  
    do m=1,n_vertex_max
      nodes(m) = node_list%node(element%vertex(m))
    enddo

    ! Set up gauss points in this element
    x_g = 0.d0; x_s = 0.d0; x_t = 0.d0; y_g = 0.d0; y_s = 0.d0; y_t = 0.d0
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

    ! Perform gauss integration of RHS
    do ms=1, n_gauss
      do mt=1, n_gauss
    
        wst  = wgauss(ms)*wgauss(mt)
        xjac = x_s(ms,mt)*y_t(ms,mt) - x_t(ms,mt)*y_s(ms,mt)

        do i=1,n_vertex_max
          do j=1,n_degrees

            index = 2*(nodes(i)%index(j)-1) + 1

            v   = h(i,j,ms,mt)  * element%size(i,j)
            rhs(index) = rhs(index) + f(x_g(ms,mt), y_g(ms,mt)) * v * xjac * x_g(ms,mt) * wst * TWOPI
          enddo
        enddo
      enddo
    enddo
  enddo
end subroutine calc_rhs_f



subroutine elements_mean_rms(node_list, element_list, f, mean, rms, volume)
  use basis_at_gaussian
  use constants, only: TWOPI
  use mod_interp, only: interp
  type(type_node_list), intent(in) :: node_list
  type(type_element_list), intent(in) :: element_list
  real*8, external :: f
  real*8, intent(out) :: mean, rms
  real*8, intent(out), optional :: volume

  integer :: i_elm, m, i, j, ms, mt
  type(type_element) :: element
  type(type_node) :: nodes(4)
  real*8, dimension(n_gauss,n_gauss) :: x_g, x_s, x_t, y_g, y_s, y_t
  real*8 :: my_ref, wst, my_volume, xjac
  real*8 :: P, P_s, P_t, P_st, P_ss, P_tt

  call initialise_basis

  my_volume = 0.d0
  mean = 0.d0
  rms = 0.d0

  do i_elm=1,element_list%n_elements

    element = element_list%element(i_elm)
    do m=1,n_vertex_max
      nodes(m) = node_list%node(element%vertex(m))
    enddo

    ! Set up gauss points in this element
    x_g = 0.d0; x_s = 0.d0; x_t = 0.d0; y_g = 0.d0; y_s = 0.d0; y_t = 0.d0
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

    ! Perform gauss integration of LHS
    do ms=1, n_gauss
      do mt=1, n_gauss

        wst       = wgauss(ms)*wgauss(mt)
        xjac      =  x_s(ms,mt)*y_t(ms,mt) - x_t(ms,mt)*y_s(ms,mt)
        my_volume = my_volume + TWOPI * x_g(ms,mt) * xjac * wst

        ! calculate contribution to integral of this point
        call interp(node_list, element_list, i_elm, 1, 1, Xgauss(ms), Xgauss(mt), P, P_s, P_t, P_st, P_ss, P_tt)
      
        rms  = rms  + (P-f(x_g(ms,mt),y_g(ms,mt)))**2 * xjac * TWOPI * x_g(ms,mt) * wst
        mean = mean + P * xjac * TWOPI * x_g(ms,mt) * wst
        
        if (xjac .lt. 0) write(*,*) i_elm, ms, mt, xjac, x_g(ms,mt), wst
        
      enddo
    enddo
  enddo

  rms = sqrt(rms / my_volume)
  mean = mean/my_volume
  if (present(volume)) volume = my_volume
end subroutine elements_mean_rms
end module projection_helpers

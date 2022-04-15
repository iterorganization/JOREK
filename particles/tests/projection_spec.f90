!> This module contains some testcases for projections, ensuring
!> that the projection matrix, RHS and MUMPS work for these cases.
!>
!> It contains tests of projecting zero, 1, x, xy, x^4 onto a square, flux-aligned or circular grid
module projection_spec
use mod_project_particles
use data_structure
use mod_particle_types
use projection_helpers
use fruit
implicit none

logical, parameter :: write_proj_output = .true. !< Set to true to write restart files with the projected density
#ifdef PROJECTION_EXTRATEST
  logical, parameter :: EXTRATEST = .true. !< Set to .true. to do flux-aligned projection tests
#else
  logical, parameter :: EXTRATEST = .false. !< Set to .true. to do flux-aligned projection tests
#endif

contains

!> Actions to perform before any of these tests
subroutine setup_projection_spec
  use basis_at_gaussian, only: initialise_basis
  call initialise_basis !< Calculate the basis functions at the gaussian points
end subroutine setup_projection_spec


!> Project zero onto a square grid
subroutine test_project_0_square_10_10
  use nodes_elements
  use mpi_mod
  write(*,*) '  test_project_0_square_10_10'
  call default_square_grid(node_list, element_list, 10)
  call project_f_with_assert_and_write(node_list, element_list, f_0, 0.d0, 0.d0, '0_square_10_10')
  write(*,*) ' completed test_project_0_square_10_10'
end subroutine test_project_0_square_10_10

!> Project one onto a square grid
subroutine test_project_1_square_10_10
  use nodes_elements
  integer :: index(1:10*10*4), i, j

  call default_square_grid(node_list, element_list, 10)

  ! include a grid_bezier_square test here
  ! verify that no node shares the same index and all indices are used exactly once
  index = 0
  do i=1,node_list%n_nodes
    do j=1,n_degrees
      index(node_list%node(i)%index(j)) = index(node_list%node(i)%index(j)) + 1
    end do
  end do
  call assert_equals(size(index,1), count(index .gt. 0), 'all indices must be used')
  call assert_equals(0, count(index .gt. 1), 'no duplicate indices in this grid')
  call project_f_with_assert_and_write(node_list, element_list, f_1, 1.d0, 0.d0, '1_square_10_10')
  write(*,*) ' completed test_project_1_square_10_10'
end subroutine test_project_1_square_10_10

!> Project one onto two polar grids. One with an even number of elements in the poloidal direction
!> and one with an odd number of elements. For the polar grid this should not matter much.
subroutine test_project_1_polar_30_32
  use phys_module
  use nodes_elements
  call default_polar_grid(node_list, element_list, 32)
  call project_f_with_assert_and_write(node_list, element_list, f_1, 1.d0, 0.d0, '1_polar_30_32')
  write(*,*) ' completed test_project_1_polar_30_31'
end subroutine test_project_1_polar_30_32

subroutine test_project_1_polar_30_31
  use phys_module
  use nodes_elements
  call default_polar_grid(node_list, element_list, 31)
  call project_f_with_assert_and_write(node_list, element_list, f_1, 1.d0, 0.d0, '1_polar_30_31')
  write(*,*) ' completed test_project_1_polar_30_31'
end subroutine test_project_1_polar_30_31

!> Test projection onto a flux-aligned grid with odd number of nodes
subroutine test_project_1_flux_40_31
  use phys_module
  use nodes_elements
  if (.not. EXTRATEST) return
  call default_flux_grid_31(node_list, element_list)
  call project_f_with_assert_and_write(node_list, element_list, f_1, 1.d0, 0.d0, '1_flux_40_31')
  write(*,*) ' completed test_project_1_flux_40_31'
end subroutine test_project_1_flux_40_31
!> Test projection onto a flux-aligned grid with even number of nodes
!> The flux-aligned grid has some tricks to have continuity on axis, which only work for an even
!> number of elements.
subroutine test_project_1_flux_40_32
  use phys_module
  use nodes_elements
  if (.not. EXTRATEST) return
  call default_flux_grid_32(node_list, element_list)
  call project_f_with_assert_and_write(node_list, element_list, f_1, 1.d0, 0.d0, '1_flux_40_32')
  write(*,*) ' completed test_project_1_flux_40_32'
end subroutine test_project_1_flux_40_32


!> Project R onto a square grid
subroutine test_project_R_square_10_10
  use nodes_elements
  call default_square_grid(node_list, element_list, 10)
  call project_f_with_assert_and_write(node_list, element_list, f_R, 26.d0/24.d0, 0.d0, 'R_square_10_10')
  write(*,*) ' completed test_project_R_square_10_10'
end subroutine test_project_R_square_10_10

!> Project RZ onto a square grid
subroutine test_project_RZ_square_10_10
  use nodes_elements
  call default_square_grid(node_list, element_list, 10)
  call project_f_with_assert_and_write(node_list, element_list, f_RZ, 0.d0, 0.d0, 'RZ_square_10_10')
  write(*,*) ' completed test_project_RZ_square_10_10'
end subroutine test_project_RZ_square_10_10

!> Project R^4 onto a square grid
subroutine test_project_R4_square_10_10
  use nodes_elements
  call default_square_grid(node_list, element_list, 10)
  call project_f_with_assert_and_write(node_list, element_list, f_R4, 1.89583333333333d0, 0.d0, 'R4_square_10_10', rms_tol=3d-6)
  write(*,*) ' completed test_project_R4_square_10_10'
end subroutine test_project_R4_square_10_10

!> Project R^4 onto a few square grids and verify convergence with n
subroutine test_project_R4_square_n
  use nodes_elements
  integer, parameter :: n_max = 40
  integer            :: i
  character(len=13)  :: s
  real*8             :: tol
  do i=3,n_max,4
    write(s,'(A,i0.3)') "R4_square_", i
    call default_square_grid(node_list, element_list, i)
    tol = 0.06d0/(real(i)**4)
    call project_f_with_assert_and_write(node_list, element_list, f_R4, 1.89583333333333d0, 0.d0, s, rms_tol=tol)
  end do
  write(*,*) ' completed test_project_R4_square_n'
end subroutine test_project_R4_square_n



!> Peaked gaussian at R=1.0 with width a.
!> On a grid from 0.5 to 1.5 in R, -0.5 to 0.5 in Z the integral is given by
!> \[
!> a^2 2\pi^2 \erf\left(\frac{0.5}{a}\right) \erf\left(\frac{0.5}{a}\right)
!> \]
!> for a=5e-3 this yields 2pi*0.0000785398 (where the volume is 2pi)
!> the RMS value is the integral of (f_peak - <mean>)^2 over the volume
!> and is given by 0.0000397457
function f_peak(R, Z)
  real*8, intent(in) :: R, Z
  real*8 :: f_peak
  real*8, parameter :: a = 5d-3
  f_peak = exp(-((R-1.d0)**2+Z**2)/a**2)
end function f_peak

!> Project a peaked function onto a square grid
subroutine test_project_peak_square_20_20
  use nodes_elements
  real*8 :: ref = 0.785398d-4
  real*8 :: rms = 0.397457d-4
  call default_square_grid(node_list, element_list, 20)
  call project_f_with_assert_and_write(node_list, element_list, f_peak, ref, rms,&
      'peak_square_20_20', mean_tol=ref, rms_tol=3d-6, filter=1d-4, hyper_filter=4d-8)
  ! the mean-tol is enormous since we do not reproduce this well for very peaked distributions, with ngauss=4
  ! if ngauss=8 we get to within 20%, but the RMS is much larger
end subroutine test_project_peak_square_20_20




!> Helper function:
!> Project a function f onto grid in node_list and element_list
!> and test for mean and RMS value. Optionally write to file for visual inspection.
subroutine project_f_with_assert_and_write(node_list, element_list, f, mean, RMS, name, mean_tol, rms_tol, filter, hyper_filter)
  use mpi_mod
  type(type_node_list), intent(inout)    :: node_list
  type(type_element_list), intent(inout) :: element_list
  real*8, external   :: f
  real*8, intent(in) :: mean
  real*8, intent(in) :: rms
  real*8, intent(in), optional :: mean_tol, rms_tol
  real*8, intent(in), optional :: filter, hyper_filter !< smoothing and hyper-smoothing factor
  real*8  :: my_filter, my_hyper_filter
  character(len=*), intent(in) :: name
  real*8  :: m, e, my_rms_tol, my_mean_tol
  integer :: ierr

  my_filter       = 0.d0
  my_hyper_filter = 0.d0

  if (present(filter))       my_filter       = filter
  if (present(hyper_filter)) my_hyper_filter = hyper_filter
  
  call project_f(node_list, element_list, f, filter, hyper_filter)
  
  ! test rms
  my_mean_tol = 1d-12
  my_rms_tol  = 1d-12
  if (present(mean_tol)) my_mean_tol = mean_tol
  if (present(rms_tol))  my_rms_tol  = rms_tol

  call elements_mean_rms(node_list, element_list, f, m, e)
  
  call assert_equals(mean, m, my_mean_tol, 'mean value M')
  call assert_equals(RMS, e, my_rms_tol, 'rms error ok')
  
  write(*,*) 'project_f : mean = ',mean,m
  write(*,*) 'project_f : RMS  = ',RMS,e

  if (write_proj_output) then
    call write_particle_distribution_to_h5(node_list, element_list, filename=name//'.h5', n_fields=1, time=0.d0)
  end if

end subroutine project_f_with_assert_and_write



!> Test the exact form of the projection matrix for a simple square grid.
!> Reference integrals calculated with Mathematica.
subroutine test_projection_matrix_square_2_2
  use mpi_mod
  use basis_at_gaussian
  use nodes_elements
  use phys_module, only : TWOPI
  type(DMUMPS_STRUC)  :: p
  integer             :: i, j

  integer, parameter  :: n_R = 2, n_Z = 2
  real*8, parameter   :: R_geo = 1.d0, Z_geo = 0.d0, amin = 0.5d0

  real*8, parameter   :: ref(16) = [0.1008163, 0.0159184,  0.0142177,   0.0022449, & ! index 1
                                    0.0477551,-0.0110544,  0.00673469, -0.00155896,& ! index 2
                                    0.034898,  0.0055102, -0.00840136, -0.00132653,& ! index 3 (but node 4, because index is switched with matrix order)
                                    0.0165306,-0.00382653,-0.00397959,  0.000921202] ! index 4 (but node 3)
  real*8, parameter   :: tol = 1d-6
  real*8              :: area, volume
  real*8, allocatable :: integral_weights(:)
  integer :: mpi_comm_n, mpi_comm_master, i_tor_local, n_tor_local, ierr

  call MPI_Comm_dup(MPI_COMM_WORLD, mpi_comm_n, ierr)
  call MPI_Comm_dup(MPI_COMM_WORLD, mpi_comm_master, ierr)

  i_tor_local     = 1
  n_tor_local     = 1

  node_list%n_nodes       = 0
  element_list%n_elements = 0
  
  call grid_bezier_square(n_R, n_Z, R_geo-amin,R_geo+amin, Z_geo-amin, Z_geo+amin, .true., node_list, element_list)
 
  if (i_tor_local .eq. 1) then 
  
    call prepare_mumps_par_n0(node_list, element_list, n_tor_local, i_tor_local, mpi_comm_world, mpi_comm_n, mpi_comm_master, &
                              p,  area, volume, filter=0.d0, filter_hyper=0.d0, filter_parallel=0.d0, integral_weights=integral_weights )

  else

    call prepare_mumps_par(node_list, element_list, n_tor_local, i_tor_local, mpi_comm_world, mpi_comm_n, mpi_comm_master, &
                           p, filter=0.d0, filter_hyper=0.d0, filter_parallel=0.d0)
  endif

  do i=1,size(p%irn), 2
    if (p%irn(i) == 1) call assert_equals(ref((p%jcn(i)-1)/2+1)*TWOPI, p%A(i), tol, 'matrix element must match reference')
  end do
  
  write(*,*) 'completed test_projection_matrix_square_2_2'

end subroutine test_projection_matrix_square_2_2



!> Test the construction of the projection matrix with and without openmp for a
!> simple grid. This is quite a slow test so it is in the EXTRATEST suite
subroutine test_omp_projection_matrix_construction
  use basis_at_gaussian
  !$use omp_lib
  use nodes_elements
  type(DMUMPS_STRUC) :: p_seq, p_par

  integer, parameter :: n_R = 10, n_Z = 10
  real*8, parameter  :: R_geo = 1.d0, Z_geo = 0.d0, amin = 0.5d0

  integer            :: i, j, n_threads
  real*8, allocatable, dimension(:,:) :: A_par, A_seq
  integer            :: mpi_comm_n, mpi_comm_master, i_tor_local, n_tor_local, ierr
  character(len=11)  :: s
  n_threads = 1

  if (.not. EXTRATEST) return

  call MPI_Comm_dup(MPI_COMM_WORLD, mpi_comm_n, ierr)
  call MPI_Comm_dup(MPI_COMM_WORLD, mpi_comm_master, ierr)

  i_tor_local     = 1
  n_tor_local     = 1

  node_list%n_nodes = 0
  element_list%n_elements = 0
  call grid_bezier_square(n_R, n_Z, R_geo-amin,R_geo+amin, Z_geo-amin, Z_geo+amin, .true., node_list, element_list)

  ! Get the openmp number of threads
  !$omp parallel
    !$n_threads = omp_get_num_threads()
  !$omp end parallel
  call prepare_mumps_par(node_list, element_list, n_tor_local, i_tor_local, mpi_comm_world, mpi_comm_n, mpi_comm_master, &
                         p_par, filter=0d0, filter_hyper=0.d0, filter_parallel=0.d0, skip_factorisation=.true.)

  !$ call omp_set_num_threads(1)
  call prepare_mumps_par(node_list, element_list, n_tor_local, i_tor_local, mpi_comm_world, mpi_comm_n, mpi_comm_master, &
                         p_seq, filter=0d0, filter_hyper=0.d0, filter_parallel=0.d0, skip_factorisation=.true.)
  !$ call omp_set_num_threads(n_threads)

  ! Check that p_par and p_seq contain the same matrix
  ! a bit slow perhaps. Generate one check per matrix element
  allocate(A_seq(minval(p_par%irn):maxval(p_par%irn),minval(p_par%jcn):maxval(p_par%jcn)))
  A_seq = 0.d0
  allocate(A_par, source=A_seq)
  A_par = 0.d0

  do i=1,size(p_par%A)
    A_par(p_par%irn(i),p_par%jcn(i)) = A_par(p_par%irn(i),p_par%jcn(i)) + p_par%A(i)
    A_seq(p_seq%irn(i),p_seq%jcn(i)) = A_seq(p_seq%irn(i),p_seq%jcn(i)) + p_seq%A(i)
  end do
  
  do i=minval(p_par%irn),maxval(p_par%irn)
    do j=minval(p_par%jcn),maxval(p_par%jcn)
      write(s,"(i5,A1,i5)") i,j
      if (A_par(i,j) .ne. 0.d0 .or. A_seq(i,j) .ne. 0.d0) call assert_equals(A_seq(i,j), A_par(i,j), s)
    end do
  end do
end subroutine test_omp_projection_matrix_construction
end module projection_spec

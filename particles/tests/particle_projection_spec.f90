!> This module contains some testcases for projecting particles.
!> First we test the projection and scaling with n.
!> Then we look into the effect of the filter parameter s.
module particle_projection_spec
use mod_project_particles
use data_structure
use mod_particle_types
use mod_pcg32_rng
use mod_sobseq_rng
use projection_spec, only: elements_mean_rms, f_1, default_square_grid, default_polar_grid, &
    default_flux_grid_32, default_flux_grid_31, calc_rhs_f
use constants, only: TWOPI
use fruit
implicit none

logical, parameter :: write_proj_output = .false. !< Set to true to write restart files with the projected density
logical, parameter :: EXTRATEST = .false. !< Set to .true. to do flux-aligned grid projection tests

contains

subroutine setup_particle_projection_spec
  use basis_at_gaussian
  call initialise_basis
end subroutine setup_particle_projection_spec

!> Project 10^3-10^5 particles generated with pcg onto square grid
subroutine test_square_10_10_pcg
  use nodes_elements
  call default_square_grid(node_list, element_list, 10)
  call project_n_square_10_10(node_list, element_list, [1000,10000,100000], pcg32_rng(), mean_tol=3d-8*[1,1,1], rms_tol=23d0/sqrt([1d3,1d4,1d5]))
  write(*,*) ' completed test_square_10_10_pcg'
end subroutine test_square_10_10_pcg

!> Project 10^3-10^5 particles generated with sobseq onto square grid
subroutine test_square_10_10_sob
  use nodes_elements
  call default_square_grid(node_list, element_list, 10)
  call project_n_square_10_10(node_list, element_list, [1000,10000,100000], sobseq_rng(), mean_tol=3d-8*[1,1,1], rms_tol=405d0/[1d3,1d4,1d5])
  write(*,*) ' completed test_square_10_10_sob'
end subroutine test_square_10_10_sob


!> Project 10^3-10^5 particles generated with sobseq onto polar grid
subroutine test_polar_30_22_sob
  use nodes_elements
  call default_polar_grid(node_list, element_list, 22)
  call project_n_polar_30_npol(node_list, element_list, [1000,10000,100000], sobseq_rng(), mean_tol=3d-5*[1,1,1], rms_tol=5d4/[1d3,1d4,1d5])
  write(*,*) ' completed test_polar_30_22_sob'
end subroutine test_polar_30_22_sob
!> Project 10^3-10^5 particles generated with sobseq onto polar grid
subroutine test_polar_30_21_sob
  use nodes_elements
  call default_polar_grid(node_list, element_list, 21)
  call project_n_polar_30_npol(node_list, element_list, [1000,10000,100000], sobseq_rng(), mean_tol=3d-5*[1,1,1], rms_tol=5d4/[1d3,1d4,1d5])
  write(*,*) ' completed test_polar_30_22_sob'
end subroutine test_polar_30_21_sob

!> Project 10^3-10^5 particles generated with sobseq onto flux grid (odd)
subroutine test_flux_40_31_pcg
  use phys_module
  use nodes_elements
  real*8 :: w
  if (.not. EXTRATEST) return
  call default_flux_grid_31(node_list, element_list)
  w=TWOPI**2*R_geo*amin**2/2.d0
  call project_n(node_list, element_list, [1000,10000,100000], pcg32_rng(), 'flux_40_31', volume=w, mean_tol=2d-5*[1,1,1], rms_tol=45d0/sqrt([1d3,1d4,1d5]))
  write(*,*) ' completed test_flux_40_31_pcg'
end subroutine test_flux_40_31_pcg
!> Project 10^3-10^5 particles generated with sobseq onto flux grid (even)
subroutine test_flux_40_32_pcg
  use phys_module
  use nodes_elements
  real*8 :: w
  if (.not. EXTRATEST) return
  call default_flux_grid_32(node_list, element_list)
  w=TWOPI**2*R_geo*amin**2/2.d0
  call project_n(node_list, element_list, [1000,10000,100000], pcg32_rng(), 'flux_40_32', volume=w, mean_tol=2d-5*[1,1,1], rms_tol=45d0/sqrt([1d3,1d4,1d5]))
  write(*,*) ' completed test_polar_40_32_pcg'
end subroutine test_flux_40_32_pcg


!> Create a RHS by integrating f_1 and with monte carlo methods and check that they are close
!> This guards against errors in node indices etc
subroutine rhs_convergence_square_10_10(node_list, element_list, n, rng)
  use mod_initialise_particles
  use mod_particle_sim
  use phys_module, only : TWOPI

  type(type_node_list),    intent(inout) :: node_list
  type(type_element_list), intent(inout) :: element_list
  integer, intent(in) :: n
  class(type_rng), intent(in) :: rng

  integer :: n_AA, ielm_out, ifail, i, j
  real*8  :: R_out, Z_out, s_out, t_out, tol
  real*8, allocatable, dimension(:) :: rhs_f
  character*8        :: n_s
  type(projection)   :: p
  type(particle_sim) :: sim
  integer :: i_elm, inode, index_ij, index_large_i
  real*8, dimension(:), allocatable :: my_rhs

  write(n_s, '(i8)') n

  allocate(p%node_list,    source=node_list)
  allocate(p%element_list, source=element_list)

  p%node_list    = node_list
  p%element_list = element_list

  p%n_tor_local = 1
  p%i_tor_local = 1

  allocate(p%f(1))
  p%f(1)%group = 1
  p%f(1)       = proj_f(proj_one, group = 1)

  allocate(sim%groups(1))

  n_AA = maxval(node_list%node(1:node_list%n_nodes)%index(4))
  allocate(rhs_f(2*n_AA), my_rhs(2*n_AA))

  my_rhs = 0.d0
  allocate(particle_fieldline::sim%groups(1)%particles(n))

  ! to prevent omp trouble
  call find_RZ(node_list,element_list,2.d0,1.d0,R_out,Z_out,ielm_out,s_out,t_out,ifail)

  call initialise_particles(sim%groups(1)%particles, node_list, element_list, rng)
  do i=1,n
    sim%groups(1)%particles(i)%weight = real(TWOPI,8)/real(n,8)
  end do

  call sample_rhs(p, sim)
  
  deallocate(sim%groups(1)%particles)

  ! Convert RHS per element to 1D
  do i_elm=1,element_list%n_elements
    do i=1,n_vertex_max
      inode = element_list%element(i_elm)%vertex(i)
      do j=1,n_order+1

        index_large_i = 2*(node_list%node(inode)%index(j)-1)+1  ! base index in the main matrix

        my_rhs(index_large_i) = my_rhs(index_large_i) + p%rhs_f(j,i,i_elm,1,1) 
      enddo
    enddo
  enddo

  call calc_rhs_f(node_list,element_list,f_1,rhs_f)  ! rhs_f has size 2*n_aa

  call assert_false(isnan(sum(rhs_f)), 'sum integrated rhs is not nan[n='//trim(adjustl(n_s))//']')
  call assert_false(isnan(sum(my_rhs)),'sum MC rhs is not nan[n='//trim(adjustl(n_s))//']')

  select type (rng)
  type is (pcg32_rng)
    tol = 7d0/sqrt(real(n)) * TWOPI
  type is (sobseq_rng)
    tol = 65d0/real(n)      * TWOPI
  end select

  call assert_equals(0.d0,    sum(abs(rhs_f-my_rhs)), tol, 'sum abs integrated - MC rhs [n='//trim(adjustl(n_s))//']')
  call assert_equals(0.d0, maxval(abs(rhs_f-my_rhs)), tol, 'max abs integrated - MC rhs [n='//trim(adjustl(n_s))//']')
end subroutine rhs_convergence_square_10_10

!> Test convergence of RHS for n particles
subroutine test_rhs_square_10_10_pcg
  use nodes_elements
  call default_square_grid(node_list, element_list, 10)
  call rhs_convergence_square_10_10(node_list, element_list,   1000, pcg32_rng())
  call rhs_convergence_square_10_10(node_list, element_list,  10000, pcg32_rng())
  call rhs_convergence_square_10_10(node_list, element_list, 100000, pcg32_rng())
  call rhs_convergence_square_10_10(node_list, element_list,1000000, pcg32_rng())
  write(*,*) ' completed test_rhs_square_10_10_pcg'
end subroutine test_rhs_square_10_10_pcg
!> Test convergence of RHS for n particles
subroutine test_rhs_square_10_10_sob
  use nodes_elements
  call default_square_grid(node_list, element_list, 10)
  call rhs_convergence_square_10_10(node_list, element_list,   1000, sobseq_rng())
  call rhs_convergence_square_10_10(node_list, element_list,  10000, sobseq_rng())
  call rhs_convergence_square_10_10(node_list, element_list, 100000, sobseq_rng())
  call rhs_convergence_square_10_10(node_list, element_list,1000000, sobseq_rng())
  write(*,*) ' completed test_rhs_square_10_10_sob'
end subroutine test_rhs_square_10_10_sob


!*************************
! filter parameter test
!*************************

!> Test convergence of RHS for 10000 particles with varying filter factor
subroutine test_polar_30_22_10000_sob_smoothing
  use phys_module
  use nodes_elements
  real*8 :: w, s, x
  integer :: i
  character(len=8) :: ss
  call default_polar_grid(node_list, element_list, 22)
  w=TWOPI**2*R_geo*amin**2/2.d0
  do i=1,7
    x = real(i-8)
    s = 10d0**x
    write(ss,'(g8.1)') s
    write(*,*) ' smoothing : s =',s
    call project_n(node_list, element_list, [10000], sobseq_rng(), 'polar_30_22_s'//ss, volume=w, smoothing=s, &
        ! Calculate error from fit (set rms_tol to 0 to get errors) * 1.2
        rms_tol=[10.d0**(-0.0738*x**2 - 0.972*x - 3.71)*1.2], mean_tol=[2d-5])
  end do
  
  write(*,*) ' completed test_polar_30_22_10000_sob_smoothing'
end subroutine test_polar_30_22_10000_sob_smoothing



!> Helper function to project n particles onto a 10x10 square grid
subroutine project_n_square_10_10(node_list, element_list, n, rng, rms_tol, mean_tol)
  use mod_rng
  class(type_rng), intent(in)            :: rng
  integer, intent(in), dimension(:)      :: n
  type(type_node_list), intent(inout)    :: node_list
  type(type_element_list), intent(inout) :: element_list
  real*8, intent(in), dimension(:)       :: rms_tol, mean_tol

  call project_n(node_list, element_list, n, rng, 'square_10_10', volume=TWOPI, rms_tol=rms_tol, mean_tol=mean_tol)
 
  write(*,*) ' completed project_n_square_10_10'
end subroutine project_n_square_10_10

!> Helper function to project n particles onto a 30_npol polar grid
subroutine project_n_polar_30_npol(node_list, element_list, n, rng, rms_tol, mean_tol)
  use phys_module
  use mod_rng
  class(type_rng), intent(in)            :: rng
  integer, intent(in), dimension(:)      :: n
  type(type_node_list), intent(inout)    :: node_list
  type(type_element_list), intent(inout) :: element_list
  real*8, intent(in), dimension(:)       :: rms_tol, mean_tol

  call project_n(node_list, element_list, n, rng, 'polar_30_32', volume=TWOPI**2*R_geo*amin**2/2.d0, rms_tol=rms_tol, mean_tol=mean_tol)

  write(*,*) ' completed : project_n_polar_30_npol'
end subroutine project_n_polar_30_npol

!> Helper function to project n particles onto a grid in node_list, element_list with optional smoothing
!> Creates a projection type behind the scenes and uses that.
!> We also need to create a particle-sim here
subroutine project_n(node_list, element_list, n, rng, name, volume, smoothing, rms_tol, mean_tol)
  use mod_initialise_particles
  use mod_particle_sim
  type(type_node_list), intent(inout)    :: node_list
  type(type_element_list), intent(inout) :: element_list
  integer, intent(in), dimension(:)      :: n
  class(type_rng), intent(in)            :: rng
  character(len=*), intent(in)           :: name
  real*8, intent(in)                     :: volume !< total volume
  real*8, intent(in), dimension(:)       :: rms_tol, mean_tol
  real*8, optional, intent(in)           :: smoothing

  type(projection) :: p
  type(particle_sim) :: sim

  integer :: i, j, ifail, ielm_out
  real*8 :: R_out, Z_out, s, s_out, t_out
  real*8 :: m, e !< mean, rms error
  character*8 :: n_s, tol_s, ss

  s = 0.d0
  ss = ''
  if (present(smoothing)) then
    s = smoothing
    write(ss,'(g8.1)') s
  end if
  p = new_projection(node_list, element_list, filter_n0=s,  f=[proj_f(proj_one, group = 1)])

  p%n_tor_local = 1
  p%i_tor_local = 1
 
  allocate(p%node_list,    source = node_list)
  allocate(p%element_list, source = element_list)

  allocate(sim%groups(1))

  do j=1,size(n)
    
    write(n_s, '(i8)') n(j)

    allocate(particle_kinetic::sim%groups(1)%particles(n(j)))
    
    call find_RZ(node_list,element_list,2.d0,1.d0,R_out,Z_out,ielm_out,s_out,t_out,ifail)
    
    call initialise_particles(sim%groups(1)%particles, node_list, element_list, rng)
    
    do i=1,n(j)
      sim%groups(1)%particles(i)%weight = volume/real(n(j))
    end do
    
    call p%do(sim) ! now results are in node_list
    
    deallocate(sim%groups(1)%particles)

    ! test rms
    call elements_mean_rms(p%node_list, p%element_list, f_1, m, e)
    write(tol_s, '(g8.1)') mean_tol(j)

    call assert_equals(1.d0, m, mean_tol(j), 'mean value 1 [n='//trim(adjustl(n_s))//' -> tol='//trim(tol_s)//']'//trim(ss))
    write(tol_s, '(g8.1)') rms_tol(j)

    call assert_equals(0.d0, e, rms_tol(j), 'rms value 0 [n='//trim(adjustl(n_s))//' -> tol='//trim(tol_s)//']'//trim(ss))

    if (write_proj_output) then
      call write_particle_distribution_to_h5(node_list, element_list, &
        filename='part_'//trim(adjustl(n_s))//'_'//trim(name)//'.h5', n_fields=1, time=0.d0)
    end if
  end do
  call p%close_mumps()
end subroutine project_n
end module particle_projection_spec

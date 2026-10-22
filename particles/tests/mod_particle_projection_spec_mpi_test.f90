module mod_particle_projection_spec_mpi_test
use data_structure
use fruit
implicit none
private
public :: run_fruit_particle_projection_spec_mpi
!> Variables --------------------------------------
!> Set to true to wrtie restart files with projected density
logical,parameter :: write_projection_output=.false.
logical,parameter :: impose_dirichlet=.false.
logical,parameter :: rng_n_streams_round_off=.true.
integer,parameter :: message_len=100
integer,parameter :: filename_len=100
integer,parameter :: master_rank=0
real*8,parameter  :: R_particle_in=2.d0
real*8,parameter  :: Z_particle_in=1.d0
real*8,parameter  :: test_time=0.d0
integer,dimension(3),parameter :: n_particles=(/10000,100000,1000000/)
integer,dimension(1),parameter :: nx=(/10/)
integer,dimension(1),parameter :: ny=(/10/)
integer,dimension(2),parameter :: nrad=(/30,40/)
integer,dimension(4),parameter :: npol=(/22,21,31,32/)
type(type_node_list),pointer    :: test_nodes
type(type_element_list),pointer :: test_elements
integer :: rank_loc,n_tasks_loc,ifail_loc
contains
!> Fruit basket -----------------------------------
subroutine run_fruit_particle_projection_spec_mpi(rank,n_tasks,ifail)
  implicit none
  integer,intent(inout) :: ifail
  integer,intent(in)    :: rank,n_tasks
  write(*,'(/A)') "  ... setting-up: particle projection spec mpi"
  call setup(rank,n_tasks,ifail)
  write(*,'(/A)') "  ... running: particle projection spec mpi"
  call run_test_case(test_particle_projection_square_10_10_pcg32,&
  'test_particle_projection_square_10_10_pcg32')
  call run_test_case(test_particle_projection_square_10_10_sobseq,&
  'test_particle_projection_square_10_10_sobseq')
  call run_test_case(test_particle_projection_polar_30_22_sobseq,&
  'test_particle_projection_polar_30_22_sobseq')
  call run_test_case(test_particle_projection_polar_30_21_sobseq,&
  'test_particle_projection_polar_30_21_sobseq')
  call run_test_case(test_particle_projection_flux_40_31_pcg32,&
  'test_particle_projection_flux_40_31_pcg32')
  call run_test_case(test_particle_projection_flux_40_32_pcg32,&
  'test_particle_projection_flux_40_32_pcg32')
  call run_test_case(test_particle_projection_polar_30_22_10000_sob_smoothing,&
  'test_particle_projection_polar_30_22_10000_sob_smoothing')
  call run_test_case(test_rhs_square_10_10_pcg32,&
  'test_rhs_square_10_10_pcg32')
  call run_test_case(test_rhs_square_10_10_sobseq,&
  'test_rhs_square_10_10_sobseq')
  write(*,'(/A)') "  ... tearing-down: particle projection spec mpi"
  call teardown(rank,n_tasks,ifail)
end subroutine run_fruit_particle_projection_spec_mpi

!> Set-up and tear-down ---------------------------
subroutine setup(rank,n_tasks,ifail)
  implicit none
  integer,intent(inout) :: ifail
  integer,intent(in)    :: rank,n_tasks
  rank_loc=rank; n_tasks_loc=n_tasks; ifail_loc=ifail;
  allocate(test_nodes); allocate(test_elements);
  call init_node_list(test_nodes, n_nodes_max, test_nodes%n_dof, n_var)
end subroutine setup

subroutine teardown(rank,n_tasks,ifail)
  implicit none
  integer,intent(inout) :: ifail
  integer,intent(in)    :: rank,n_tasks
  rank_loc = -1; n_tasks_loc = 0; ifail=ifail_loc;
  deallocate(test_nodes); deallocate(test_elements);
end subroutine teardown

!> Tests ------------------------------------------
!> Project 10^3-10^5 particles generated with pcg ont square grid
subroutine test_particle_projection_square_10_10_pcg32
  use constants,                         only: TWOPI
  use mod_pcg32_rng,                     only: pcg32_rng
  use mod_rhs_projections,               only: proj_one
  use mod_projection_helpers_test_tools, only: f_1,default_square_grid
  implicit none
  real*8,parameter              :: expect_mean=1.d0
  real*8,parameter              :: expect_rms=0.d0
  real*8,parameter              :: tol_rms_const=2.3d1
  real*8,parameter              :: volume=TWOPI
  real*8,dimension(3),parameter :: tol_mean=[3.d-8,3.d-8,3.d-8]
  real*8,dimension(3)           :: tol_rms
  character(len=message_len)    :: message
  character(len=filename_len)   :: filename
  tol_rms = tol_rms_const/sqrt(real(n_particles(1:3)*n_tasks_loc,kind=8))
  write(message,'(A,I0,A,I0,A,I0,A)') 'Error particle projection square nx: ',&
  nx(1),' ny: ',ny(1),' pcg32 rank: ',rank_loc,':'
  write(filename,'(A,I0,A,I0,A,I0,A)') '_test_projection_square_rank',&
  rank_loc,'_nx',nx(1),'_ny',ny(1),'_pcg32'
  call default_square_grid(rank_loc,n_tasks_loc,nx(1),nx(1),&
  test_nodes,test_elements,ifail_loc)
  call project_n(rank_loc,master_rank,n_tasks_loc,test_nodes,test_elements,&
  proj_one,f_1,n_particles(1:3),pcg32_rng(),volume,expect_mean,expect_rms,&
  tol_mean,tol_rms,trim(adjustl(message)),trim(adjustl(filename)),&
  ifail_loc,apply_dirichlet_in=impose_dirichlet,&
  write_particle_in=write_projection_output)
end subroutine test_particle_projection_square_10_10_pcg32

!> Project 10^3-10^5 particles generated with sobseq ont square grid
subroutine test_particle_projection_square_10_10_sobseq
  use constants,                         only: TWOPI
  use mod_sobseq_rng,                    only: sobseq_rng
  use mod_rhs_projections,               only: proj_one
  use mod_projection_helpers_test_tools, only: f_1,default_square_grid
  implicit none
  real*8,parameter              :: expect_mean=1.d0
  real*8,parameter              :: expect_rms=0.d0
  real*8,parameter              :: volume=TWOPI
  real*8,dimension(3),parameter :: tol_mean=[3.d-8,3.d-8,3.d-8]
  real*8,dimension(3),parameter :: tol_rms=[2.5d3/real(n_particles(1),kind=8),&
                                   2.5d3/real(n_particles(2),kind=8),&
                                   2.5d3/real(n_particles(3),kind=8)]
  character(len=message_len)    :: message
  character(len=filename_len)   :: filename
  write(message,'(A,I0,A,I0,A,I0,A)') 'Error particle projection square nx: ',&
  nx(1),' ny: ',ny(1),' sobseq rank: ',rank_loc,':'
  write(filename,'(A,I0,A,I0,A,I0,A)') '_test_projection_square_rank',&
  rank_loc,'_nx',nx(1),'_ny',ny(1),'_sobseq'
  call default_square_grid(rank_loc,n_tasks_loc,nx(1),nx(1),&
  test_nodes,test_elements,ifail_loc)
  call project_n(rank_loc,master_rank,n_tasks_loc,test_nodes,test_elements,&
  proj_one,f_1,n_particles(1:3),sobseq_rng(),volume,expect_mean,expect_rms,&
  tol_mean,tol_rms/real(n_tasks_loc,kind=8),trim(adjustl(message)),&
  trim(adjustl(filename)),ifail_loc,apply_dirichlet_in=impose_dirichlet,&
  write_particle_in=write_projection_output)
end subroutine test_particle_projection_square_10_10_sobseq

!> Project 10^3-10^5 particles generated with sobseq ont even polar grid
subroutine test_particle_projection_polar_30_22_sobseq
  use constants,                         only: TWOPI
  use phys_module,                       only: R_geo,amin
  use mod_sobseq_rng,                    only: sobseq_rng
  use mod_rhs_projections,               only: proj_one
  use mod_projection_helpers_test_tools, only: f_1,default_polar_grid
  implicit none
  real*8,parameter              :: expect_mean=1.d0
  real*8,parameter              :: expect_rms=0.d0
  real*8,dimension(3),parameter :: tol_mean=[3.d-5,3.d-5,3.d-5]
  real*8,dimension(3),parameter :: tol_rms=[5d4/real(n_particles(1),kind=8),&
                                   5d4/real(n_particles(2),kind=8),&
                                   5d4/real(n_particles(3),kind=8)]
  real*8                        :: volume
  character(len=message_len)    :: message
  character(len=filename_len)   :: filename
  write(message,'(A,I0,A,I0,A,I0,A)') 'Error particle projection polar nrad: ',&
  nrad(1),' npol: ',npol(1),' sobseq rank: ',rank_loc,':'
  write(filename,'(A,I0,A,I0,A,I0,A)') '_test_projection_polar_rank',&
  rank_loc,'_nrad',nrad(1),'_npol',npol(1),'_sobseq'
  call default_polar_grid(rank_loc,n_tasks_loc,npol(1),nrad(1),&
  test_nodes,test_elements,ifail_loc)
  volume=5d-1*R_geo*((TWOPI*amin)**2)!< compute the volume
  call project_n(rank_loc,master_rank,n_tasks_loc,test_nodes,test_elements,&
  proj_one,f_1,n_particles(1:3),sobseq_rng(),volume,expect_mean,expect_rms,&
  tol_mean,tol_rms/real(n_tasks_loc,kind=8),trim(adjustl(message)),&
  trim(adjustl(filename)),ifail_loc,apply_dirichlet_in=impose_dirichlet,&
  write_particle_in=write_projection_output)
end subroutine test_particle_projection_polar_30_22_sobseq

!> Project 10^3-10^5 particles generated with sobseq ont odd polar grid
subroutine test_particle_projection_polar_30_21_sobseq
  use constants,                         only: TWOPI
  use phys_module,                       only: R_geo,amin
  use mod_sobseq_rng,                    only: sobseq_rng
  use mod_rhs_projections,               only: proj_one
  use mod_projection_helpers_test_tools, only: f_1,default_polar_grid
  implicit none
  real*8,parameter              :: expect_mean=1.d0
  real*8,parameter              :: expect_rms=0.d0
  real*8,dimension(3),parameter :: tol_mean=[3.d-5,3.d-5,3.d-5]
  real*8,dimension(3),parameter :: tol_rms=[5d4/real(n_particles(1),kind=8),&
                                   5d4/real(n_particles(2),kind=8),&
                                   5d4/real(n_particles(3),kind=8)]
  real*8                        :: volume
  character(len=message_len)    :: message
  character(len=filename_len)   :: filename
  write(message,'(A,I0,A,I0,A,I0,A)') 'Error particle projection polar nrad: ',&
  nrad(1),' npol: ',npol(2),' sobseq rank: ',rank_loc,':'
  write(filename,'(A,I0,A,I0,A,I0,A)') '_test_projection_polar_rank',&
  rank_loc,'_nrad',nrad(1),'_npol',npol(2),'_sobseq'
  call default_polar_grid(rank_loc,n_tasks_loc,npol(2),nrad(1),&
  test_nodes,test_elements,ifail_loc)
  volume=5d-1*R_geo*((TWOPI*amin)**2)!< compute the volume
  call project_n(rank_loc,master_rank,n_tasks_loc,test_nodes,test_elements,&
  proj_one,f_1,n_particles(1:3),sobseq_rng(),volume,expect_mean,expect_rms,tol_mean,&
  tol_rms/real(n_tasks_loc,kind=8),trim(adjustl(message)),trim(adjustl(filename)),&
  ifail_loc,apply_dirichlet_in=impose_dirichlet,write_particle_in=write_projection_output)
end subroutine test_particle_projection_polar_30_21_sobseq

!> Project 10^3-10^5 particles generated with pcg32 ont odd flux grid
subroutine test_particle_projection_flux_40_31_pcg32
  use constants,                         only: TWOPI
  use phys_module,                       only: R_geo,amin
  use mod_pcg32_rng,                     only: pcg32_rng
  use mod_rhs_projections,               only: proj_one
  use mod_projection_helpers_test_tools, only: f_1,default_flux_grid
  implicit none
  real*8,parameter              :: expect_mean=1.d0
  real*8,parameter              :: expect_rms=0.d0
  real*8,parameter              :: tol_rms_const=4.5d1
  real*8,dimension(3),parameter :: tol_mean=[7.d-5,7.d-5,7.d-5]
  real*8                        :: volume
  real*8,dimension(3)           :: tol_rms
  character(len=message_len)    :: message
  character(len=filename_len)   :: filename
  tol_rms = tol_rms_const/sqrt(real(n_tasks_loc*n_particles(1:3),kind=8))
  write(message,'(A,I0,A,I0,A,I0,A)') 'Error particle projection flux nrad: ',&
  nrad(2),' npol: ',npol(3),' pcg32 rank: ',rank_loc,':'
  write(filename,'(A,I0,A,I0,A,I0,A)') '_test_projection_flux_rank',&
  rank_loc,'_nrad',nrad(2),'_npol',npol(3),'_pcg32'
  call default_flux_grid(rank_loc,n_tasks_loc,npol(3),nrad(2),&
  test_nodes,test_elements,ifail_loc)
  volume=5d-1*R_geo*((TWOPI*amin)**2)!< compute the volume
  call project_n(rank_loc,master_rank,n_tasks_loc,test_nodes,test_elements,&
  proj_one,f_1,n_particles(1:3),pcg32_rng(),volume,expect_mean,expect_rms,tol_mean,&
  tol_rms,trim(adjustl(message)),trim(adjustl(filename)),ifail_loc,&
  apply_dirichlet_in=impose_dirichlet,write_particle_in=write_projection_output)
end subroutine test_particle_projection_flux_40_31_pcg32

!> Project 10^3-10^5 particles generated with pcg32 ont even flux grid
subroutine test_particle_projection_flux_40_32_pcg32
  use constants,                         only: TWOPI
  use phys_module,                       only: R_geo,amin
  use mod_pcg32_rng,                     only: pcg32_rng
  use mod_rhs_projections,               only: proj_one
  use mod_projection_helpers_test_tools, only: f_1,default_flux_grid
  implicit none
  real*8,parameter              :: expect_mean=1.d0
  real*8,parameter              :: expect_rms=0.d0
  real*8,parameter              :: tol_rms_const=4.5d1
  real*8,dimension(3),parameter :: tol_mean=[7.d-5,7.d-5,7.d-5]
  real*8                        :: volume
  real*8,dimension(3)           :: tol_rms
  character(len=message_len)    :: message
  character(len=filename_len)   :: filename
  tol_rms = tol_rms_const/sqrt(real(n_tasks_loc*n_particles(1:3),kind=8))
  write(message,'(A,I0,A,I0,A,I0,A)') 'Error particle projection flux nrad: ',&
  nrad(2),' npol: ',npol(4),' pcg32 rank: ',rank_loc,':'
  write(filename,'(A,I0,A,I0,A,I0,A)') '_test_projection_flux_rank',&
  rank_loc,'_nrad',nrad(2),'_npol',npol(4),'_pcg32'
  call default_flux_grid(rank_loc,n_tasks_loc,npol(4),nrad(2),&
  test_nodes,test_elements,ifail_loc)
  volume=5d-1*R_geo*((TWOPI*amin)**2)!< compute the volume
  call project_n(rank_loc,master_rank,n_tasks_loc,test_nodes,test_elements,&
  proj_one,f_1,n_particles(1:3),pcg32_rng(),volume,expect_mean,expect_rms,tol_mean,&
  tol_rms,trim(adjustl(message)),trim(adjustl(filename)),ifail_loc,&
  apply_dirichlet_in=impose_dirichlet,write_particle_in=write_projection_output)
end subroutine test_particle_projection_flux_40_32_pcg32

!> Test convergence of RHS for 10000 particles with varying filter factor
subroutine test_particle_projection_polar_30_22_10000_sob_smoothing
  use constants,                         only: TWOPI
  use phys_module,                       only: R_geo,amin
  use mod_sobseq_rng,                    only: sobseq_rng
  use mod_rhs_projections,               only: proj_one
  use mod_projection_helpers_test_tools, only: f_1,default_polar_grid
  implicit none
  integer,parameter              :: n_trials=7
  real*8,parameter               :: expect_mean=1.d0
  real*8,parameter               :: expect_rms=0.d0
  real*8,dimension(1),parameter  :: tol_mean=(/2.d-5/)
  integer              :: ii
  real*8               :: weight,x,smoothing
  real*8,dimension(1)  :: tol_rms 
  character(len=message_len)    :: message
  character(len=filename_len)   :: filename
  call default_polar_grid(rank_loc,n_tasks_loc,npol(1),nrad(1),&
  test_nodes,test_elements,ifail_loc)
  weight = 5d-1*R_geo*((amin*TWOPI)**2)
  do ii=1,n_trials
    x = real(ii-n_trials+1,kind=8); smoothing =1d1**(ii-n_trials+1);
    write(message,'(A,I0,A,I0,A,I0,A)') 'Error particle projection polar smoothing nrad: ',&
    nrad(1),' npol: ',npol(1),' pcg32 rank: ',rank_loc,':'
    write(filename,'(A,I0,A,I0,A,I0,A)') '_test_projection_smoothing_rank',&
    rank_loc,'_nrad',nrad(1),'_npol',npol(1),'_pcg32'
    tol_rms = (/(10.d0**(-0.0738*x**2 - 0.972*x - 3.71))*1.75/)
    write(*,*) 'tol: ',tol_rms
    call project_n(rank_loc,master_rank,n_tasks_loc,test_nodes,test_elements,&
    proj_one,f_1,[n_particles(2)],sobseq_rng(),weight,expect_mean,expect_rms,tol_mean,&
    tol_rms,trim(adjustl(message)),trim(adjustl(filename)),ifail_loc,smoothing_in=smoothing,&
    apply_dirichlet_in=impose_dirichlet,write_particle_in=write_projection_output)   
  enddo
end subroutine test_particle_projection_polar_30_22_10000_sob_smoothing

!> Test convergence of RHS for n_particles on square grid PCG32
subroutine test_rhs_square_10_10_pcg32
  use constants,                         only: TWOPI
  use mod_pcg32_rng,                     only: pcg32_rng
  use mod_rhs_projections,               only: proj_one
  use mod_projection_helpers_test_tools, only: f_1,default_square_grid 
  implicit none
  real*8,parameter            :: expect_mean=0.d0
  integer                     :: ii
  real*8                      :: tol
  character(len=message_len)  :: message
  write(message,'(A,I0,A,I0,A,I0,A)') 'Error RHS convergence square nx: ',&
  nx(1),' ny: ',ny(1),' pcg32 rank: ',rank_loc,':'
  call default_square_grid(rank_loc,n_tasks_loc,nx(1),nx(1),&
  test_nodes,test_elements,ifail_loc)
  do ii=1,size(n_particles)
    tol = TWOPI*(8d0/sqrt(real(n_particles(ii),kind=8)))
    call rhs_convergence(rank_loc,n_tasks_loc,test_nodes,&
    test_elements,n_particles(ii),expect_mean,tol,f_1,proj_one,&
    pcg32_rng(),message,ifail_loc,apply_dirichlet_in=impose_dirichlet)
  enddo
end subroutine test_rhs_square_10_10_pcg32

!> Test convergence of RHS for n_particles on square grid sobseq
subroutine test_rhs_square_10_10_sobseq
  use constants,                         only: TWOPI
  use mod_sobseq_rng,                    only: sobseq_rng
  use mod_rhs_projections,               only: proj_one
  use mod_projection_helpers_test_tools, only: f_1,default_square_grid 
  implicit none
  real*8,parameter            :: expect_mean=0.d0
  integer                     :: ii
  real*8                      :: tol
  character(len=message_len)  :: message
  write(message,'(A,I0,A,I0,A,I0,A)') 'Error RHS convergence square nx: ',&
  nx(1),' ny: ',ny(1),' sobseq rank: ',rank_loc,':'
  call default_square_grid(rank_loc,n_tasks_loc,nx(1),nx(1),&
  test_nodes,test_elements,ifail_loc)
  do ii=1,size(n_particles)
    tol = TWOPI*(67d0/sqrt(real(n_particles(ii),kind=8)))
    call rhs_convergence(rank_loc,n_tasks_loc,test_nodes,&
    test_elements,n_particles(ii),expect_mean,tol,f_1,proj_one,&
    sobseq_rng(),message,ifail_loc,apply_dirichlet_in=impose_dirichlet)
  enddo
end subroutine test_rhs_square_10_10_sobseq

!> Tool -------------------------------------------
!> Helper function to project n particles onto a grid in node_list,
!> element list with optional smoothing. It creates a projection
!> type behind the scenes and uses that. We also need to create
!> a particle-sim here.
subroutine project_n(rank,master,n_tasks,node_list,element_list,proj_f_proj,&
f_proj,n,rng,volume,mean_expect,rms_expect,mean_tol,rms_tol,message,&
fname,ifail,smoothing_in,n_tor_local_in,i_tor_local_in,&
apply_dirichlet_in,write_particle_in,n_fields_write_in)
  use mpi_mod
  use data_structure
  use mod_rng,                           only: type_rng
  use initialisers_base,                 only: initialise_particles
  use mod_particle_sim,                  only: particle_sim
  use mod_project_particles,             only: projection,new_projection, write_particle_distribution_to_h5
  use mod_rhs_projections,               only: proj_f
  use mod_particle_types,                only: particle_kinetic
  use mod_projection_helpers_test_tools, only: elements_mean_rms
  implicit none
  type(type_node_list),intent(inout)    :: node_list
  type(type_element_list),intent(inout) :: element_list
  class(type_rng),intent(in)            :: rng
  integer,intent(inout)                 :: ifail
  integer,intent(in)                    :: rank,master,n_tasks
  integer,dimension(:),intent(in)       :: n
  real*8,intent(in)                     :: volume,mean_expect,rms_expect
  real*8,dimension(:),intent(in)        :: rms_tol,mean_tol
  character(len=*),intent(in)           :: message,fname
  integer,optional,intent(in)           :: n_tor_local_in,i_tor_local_in
  integer,optional,intent(in)           :: n_fields_write_in
  real*8,intent(in),optional            :: smoothing_in
  logical,intent(in),optional           :: apply_dirichlet_in,write_particle_in
  real*8,external                       :: proj_f_proj,f_proj
  type(projection)                      :: project
  type(particle_sim)                    :: sim
  integer                               :: ii,jj,kk,ielm_out
  integer                               :: i_tor_local,n_tor_local,n_fields_write
  real*8                                :: smoothing 
  real*8                                :: R_out,Z_out,s_out,t_out
  real*8                                :: mean,rms_error
  character*8                           :: number_particles,tol_s
  character*8                           :: smooth_string,group_string
  logical                               :: apply_dirichlet,write_particle

  !> initialisations
  smoothing = 0.d0; smooth_string = '';
  if(present(smoothing_in)) then
    smoothing = smoothing_in; write(smooth_string,'(g8.1)') smoothing_in;
  endif
  write_particle =.false.; 
  if(present(write_particle_in)) write_particle = write_particle_in;
  apply_dirichlet = .true.
  if(present(apply_dirichlet_in)) apply_dirichlet = apply_dirichlet_in
  n_fields_write = 1; if(present(n_fields_write_in)) n_fields_write = n_fields_write_in;
  n_tor_local = 1;    if(present(n_tor_local_in)) n_tor_local = n_tor_local_in;
  i_tor_local = 1;    if(present(i_tor_local_in)) i_tor_local = i_tor_local_in;
  project = new_projection(node_list,element_list,filter_n0=smoothing,&
  f=[proj_f(proj_f_proj,group=1)],do_dirichlet=apply_dirichlet)
  project%n_tor_local = n_tor_local; project%i_tor_local = i_tor_local;
  sim%my_id=rank; sim%n_mpi=n_tasks; allocate(sim%groups(1));
  !> fill the particle structure
  do kk=1,size(n)
    write(number_particles,'(I8)') n(kk)
    write(group_string,'(I8)') kk
    allocate(particle_kinetic::sim%groups(1)%particles(n(kk)))
    !> to prevent omp trouble (!?)
    call find_RZ(node_list,element_list,R_particle_in,Z_particle_in,&
    R_out,Z_out,ielm_out,s_out,t_out,ifail)
    call initialise_particles(sim%groups(1)%particles,node_list,element_list,&
    rng, rng_n_streams_round_off_in=rng_n_streams_round_off)
    do ii=1,n(kk)
      sim%groups(1)%particles(ii)%weight = volume/real(n_tasks*n(kk),kind=8)
    enddo
    call project%do(sim) !< project particles, results in node_list
    deallocate(sim%groups(1)%particles) !< cleanup
    !> perform checks on the mean and rms
    call elements_mean_rms(project%node_list,project%element_list,&
    f_proj,mean,rms_error)
    write(tol_s,'(g8.1)') mean_tol(kk)
    call assert_equals(mean_expect,mean,mean_tol(kk),message//&
    ' mean value not matched for n='//trim(adjustl(number_particles))//&
    ', tol: '//trim(tol_s)//' smoothing: '//trim(smooth_string))
    write(tol_s,'(g8.1)') rms_tol(kk)
    call assert_equals(rms_expect,rms_error,rms_tol(kk),message//&
    ' RMS value not matched for n='//trim(adjustl(number_particles))//&
    ', tol: '//trim(tol_s)//' smoothing: '//trim(smooth_string))
    if(write_particle) then
      call write_particle_distribution_to_h5(project%node_list,project%element_list,&
      filename='part_'//trim(adjustl(number_particles))//'_group_'//&
      trim(adjustl(group_string))//trim(fname)//'.h5',n_fields=n_fields_write,time=test_time)
    endif
  enddo
  deallocate(sim%groups); call project%close_projection()
end subroutine project_n

!> Create RHS by integrating f and with monte carlo methods and
!> check that they are close. This guards against errors in
!> node indices etc
subroutine rhs_convergence(rank,n_tasks,node_list,element_list,&
n_particles,mean_expect,tol,funct,f_proj,rng,message,ifail,&
n_tor_local_in,i_tor_local_in,smoothing_in,apply_dirichlet_in)
  use mpi_mod
  use constants,                         only: TWOPI
  use data_structure,                    only: type_node_list,type_element_list
  use mod_rng,                           only: type_rng
  use initialisers_base,                 only: initialise_particles
  use mod_particle_sim,                  only: particle_sim
  use mod_project_particles,             only: projection, new_projection, sample_rhs
  use mod_rhs_projections,               only: proj_f
  use mod_particle_types,                only: particle_fieldline
  use mod_projection_helpers_test_tools, only: calc_rhs_f
  implicit none
  type(type_node_list),intent(inout)    :: node_list
  type(type_element_list),intent(inout) :: element_list 
  integer,intent(inout)                 :: ifail
  class(type_rng),intent(in)            :: rng
  integer,intent(in)                    :: n_particles,rank,n_tasks
  real*8,intent(in)                     :: tol,mean_expect
  character(len=*),intent(in)           :: message
  real*8,external                       :: funct,f_proj
  integer,intent(in),optional           :: n_tor_local_in,i_tor_local_in
  real*8,intent(in),optional            :: smoothing_in
  logical,intent(in),optional           :: apply_dirichlet_in
  type(projection)                :: project
  type(particle_sim)              :: sim
  integer                         :: ii,jj,n_AA,ielm_out,i_elm
  integer                         :: inode,index_ij,index_large_i
  integer                         :: n_tor_local,i_tor_local
  real*8                          :: R_out,Z_out,s_out,t_out,smoothing
  real*8,dimension(:),allocatable :: rhs_f,my_rhs
  character*8                     :: n_particles_string
  logical                         :: apply_dirichlet
  !> initialisation
  n_tor_local=1; if(present(n_tor_local_in)) n_tor_local=n_tor_local_in;
  i_tor_local=1; if(present(i_tor_local_in)) i_tor_local=i_tor_local_in;
  smoothing=0.d0; if(present(smoothing_in)) smoothing=smoothing_in;
  apply_dirichlet=.true.; 
  if(present(apply_dirichlet_in)) apply_dirichlet=apply_dirichlet_in;
  write(n_particles_string,'(I8)') n_particles
  n_AA = maxval(node_list%node(1:node_list%n_nodes)%index(4)); 
  allocate(rhs_f(2*n_AA)); rhs_f=0.d0; allocate(my_rhs(2*n_AA)); my_rhs=0.d0;
  allocate(project%node_list,project%element_list)
  project%node_list=node_list; project%element_list=element_list;
  allocate(project%f(1)); project%f(1)%group=1; project%f(1)=proj_f(f_proj,group=1);
  project%n_tor_local = n_tor_local; project%i_tor_local = i_tor_local;
  !> initialise particles
  sim%my_id=rank; sim%n_mpi=n_tasks; allocate(sim%groups(1)); 
  allocate(particle_fieldline::sim%groups(1)%particles(n_particles));
  !> to prevent omp trouble (!?)
  call find_RZ(node_list,element_list,R_particle_in,Z_particle_in,&
  R_out,Z_out,ielm_out,s_out,t_out,ifail)
  call initialise_particles(sim%groups(1)%particles,node_list,element_list,&
  rng,rng_n_streams_round_off_in=rng_n_streams_round_off)
  do ii=1,n_particles
    sim%groups(1)%particles(ii)%weight = TWOPI/real(n_particles,kind=8)
  enddo
  call sample_rhs(project,sim) !< comput the rhs using the projection method
  !> cleanup
  deallocate(sim%groups(1)%particles); deallocate(sim%groups);
  !> convert RHS per element to 1D
  do i_elm=1,element_list%n_elements
    do ii=1,n_vertex_max
      inode = element_list%element(i_elm)%vertex(ii)
      do jj=1,n_degrees
        index_large_i = 2*(node_list%node(inode)%index(jj)-1)+1
        my_rhs(index_large_i) = my_rhs(index_large_i) + project%rhs_f(jj,ii,i_elm,1,1)
      enddo
    enddo
  enddo
  call calc_rhs_f(node_list,element_list,funct,rhs_f) !< rhs size 2*n_AA
  !> checks
  call assert_false(isnan(sum(rhs_f)),message//' sum integrated rhs is nan for n='//&
  trim(adjustl(n_particles_string)))
  call assert_false(isnan(sum(my_rhs)),message//' sum Monte Carlo rhs is nan for n='//&
  trim(adjustl(n_particles_string)))
  call assert_equals(mean_expect,sum(abs(rhs_f-my_rhs)),tol,message//' |integrated-MC rhs|_1 n='//&
  trim(adjustl(n_particles_string)))
  call assert_equals(mean_expect,maxval(abs(rhs_f-my_rhs)),tol,message//' |integrated-MC rhs|_inf n='//&
  trim(adjustl(n_particles_string)))
  call MPI_Barrier(MPI_COMM_WORLD,ifail)
  deallocate(rhs_f); deallocate(my_rhs);! call project%close_projection() !< cleanup
end subroutine rhs_convergence

!> ------------------------------------------------
end module mod_particle_projection_spec_mpi_test

!> Testing the coupling of the projections of particles to JOREK
module domain
  real*8 :: psi_start, psi_end, psi_n_start, psi_n_end
end module

program particle_diagno

use particle_tracer
use mod_pusher_tools, only : particle_position_to_gc
use mod_particle_diagnostics
use mpi
use mod_atomic_elements
use mod_particle_io
use mod_event
use mod_project_particles
use mod_particle_loop
use mod_gc_variational
use nodes_elements
use mod_jorek_timestepping
use mod_random_seed
use mod_interp, only: mode_moivre, interp_RZ, interp_0
use mod_basisfunctions
use basis_at_gaussian
use phys_module, only: F0, tstep, nstep, nout, restart, rho_0, rho_1, rho_coef
use phys_module, only: CENTRAL_MASS, CENTRAL_DENSITY, xcase, xpoint, index_now, index_start
use phys_module, only: n_particles, nstep_particles, nsubstep_particles, tstep_particles
use phys_module, only: filter_perp, filter_hyper, filter_par, filter_perp_n0, filter_hyper_n0, filter_par_n0
use phys_module, only: xtime, energies, mode, restart_particles
use constants,   only: MU_ZERO, MASS_PROTON, ATOMIC_MASS_UNIT, K_BOLTZ, EL_CHG
use mod_export_restart
use live_data

use mod_edge_domain
use mod_edge_elements
use data_structure, only: type_bnd_element_list, type_bnd_node_list 
use equil_info
use mod_boundary, only: boundary_from_grid
use domain

!$ use omp_lib

implicit none

type(event)                                       :: fieldreader, partreader, partwriter
type(count_action)                                :: counter
type(projection), target                          :: jorek_feedback, project_profiles
type(jorek_timestep_action), target               :: jorek_stepper
type(type_edge_domain), allocatable, dimension(:) :: edge_domains
type(edge_elements)                               :: D_edge
type(write_particle_diagnostics)                  :: diag
!type(type_bnd_element_list) :: bnd_elm_list !< List of boundary elements
!type(type_bnd_node_list)    :: bnd_node_list !< List of boundary nodes.

real*8    :: target_time, tstep_keep
real*8    :: physical_particles, weight
real*8    :: oldtime, step_rest_time, particle_step_time, particle_start_time, diag_time
real*8    :: rho_norm, t_norm, v_norm, E_norm, M_norm, zn_norm, Tev_norm, tstep_si, timesteps
real*8    :: E(3), B(3), psi, U, B_norm
real*8    :: psi_n, dpsi_n, ss, psi_axis, psi_bnd
integer   :: i, ifail
!$ real*8 :: w0, w1, mmm(3)

psi_n_start = 0.d0          ! limit domain to psi_start:psi_end (in normalised psi)
psi_n_end   = 1.d0 

call sim%initialize(num_groups=2)

deallocate(sim%groups)
allocate(sim%groups(0))

partreader = event(read_action(filename='part_restart.h5'))
call with(sim, partreader)

fieldreader = event(read_jorek_fields_interp_linear(basename='jorek', i=-1))
call with(sim, fieldreader)

call det_modes()

zn_norm   = CENTRAL_DENSITY * 1.d20                              ! (number) density normalisation
rho_norm  = CENTRAL_MASS * MASS_PROTON * zn_norm                 ! rho_SI = rho_norm * rho
t_norm    = sqrt((MU_ZERO * rho_norm))                           ! t_SI   = t_norm * t_jorek
Tev_norm  = 1.d0 / (2.d0 * EL_CHG * MU_ZERO * zn_norm)           ! T_ev   = Tev_norm * T (factor 2 for electron temperature)

if (sim%my_id .eq. 0) call boundary_from_grid(sim%fields%node_list, sim%fields%element_list, bnd_node_list, bnd_elm_list, .false.)

call broadcast_boundary(sim%my_id, bnd_elm_list, bnd_node_list)
call update_equil_state(sim%my_id, sim%fields%node_list, sim%fields%element_list, bnd_elm_list, xpoint, xcase)

psi_axis  = ES%psi_axis
psi_bnd   = ES%psi_bnd
psi_start = psi_axis + (psi_bnd - psi_axis) * psi_n_start
psi_end   = psi_axis + (psi_bnd - psi_axis) * psi_n_end

if (sim%my_id .eq. 0) write(*,'(A,6e14.6)') 'PSI_AXIS, PSI_BND : ',psi_axis, psi_bnd, psi_start, psi_end, psi_n_start, psi_n_end

call with(sim, counter)

select type (p_gc => sim%groups(1)%particles)
type is (particle_gc_vpar)

#ifdef __GFORTRAN__
  !$omp parallel do default(shared) & 
#else
  !$omp parallel do default(none)   &
  !$omp shared(sim, p_gc)           &
#endif
  !$omp private(i, E, B, psi, U)    &
  !$omp schedule(dynamic,10)  
  do i=1, size(p_gc,1)    
    call sim%fields%calc_EBpsiU(sim%time, p_gc(i)%i_elm, p_gc(i)%st, p_gc(i)%x(3), E, B, psi, U)
    p_gc(i)%B_norm = norm2(B)  
  enddo
end select

! SHOULD MAKE PROFILES USING THE MULTIPLE KINETIC POINTS PER GYROCENTRE

project_profiles = new_projection(sim%fields%node_list, sim%fields%element_list, &
                      filter    = filter_perp,    filter_hyper    = filter_hyper,    filter_parallel    = filter_par, &
                      filter_n0 = filter_perp_n0, filter_hyper_n0 = filter_hyper_n0, filter_parallel_n0 = filter_par_n0, &
                      f=[proj_f(proj_one, group = 1), proj_f(proj_Pressure, group = 1)], fractional_digits = 9,  &
                      do_zonal = .false., calc_integrals=.true., to_vtk=.true., to_h5=.false., basename='profiles', nsub=5)

call with(sim, project_profiles)

call sim%finalize


contains

subroutine itg_energy(node_list,element_list,psi_start,psi_end,W_kin)
!---------------------------------------------------------------
!
!---------------------------------------------------------------
use data_structure
use gauss
use basis_at_gaussian
use phys_module

implicit none

type (type_node_list)    :: node_list
type (type_element_list) :: element_list
type (type_element)      :: element
type (type_node)         :: nodes(n_vertex_max)

real*8     :: x_g(n_gauss,n_gauss),  x_s(n_gauss,n_gauss),  x_t(n_gauss,n_gauss)
real*8     :: y_s(n_gauss,n_gauss),  y_t(n_gauss,n_gauss)
real*8     :: eq_s(n_gauss,n_gauss), eq_t(n_gauss,n_gauss), psi_eq(n_gauss,n_gauss)
integer    :: i, j, k, in, ms, mt, iv, inode, ife, n_elements, i_tor
real*8     :: W_kin(n_tor), xjac, BigR, wst
real*8     :: u0_x, u0_y, psi_start, psi_end, psi_norm, psi_axis, psi_bnd
integer, parameter :: ivar_psi=1, ivar_u=2, ivar_rho=5

W_kin = 0.d0

#ifdef __GFORTRAN__
   !$omp parallel do default(shared) & ! workaround for Error: �__vtab_mod_pcg32_rng_Pcg32_rng� not specified in enclosing �parallel�
#else
   !$omp parallel do default(none) &
#endif
   !$omp schedule(dynamic,10)                                                         &
   !$omp shared(element_list, node_list, psi_start, psi_end, H, H_s, H_t)             &
   !$omp private(ife, element, nodes, iv, inode, x_g, x_s, x_t, y_s, y_t, eq_s, eq_t, &
   !$omp         i, j, ms, mt, in, wst, bigR, xjac, u0_x, u0_y, psi_eq )              &
   !$omp reduction(+:W_kin)
do ife =1,  element_list%n_elements

  element = element_list%element(ife)

  do iv = 1, n_vertex_max
    inode     = element%vertex(iv)
    nodes(iv) = node_list%node(inode)
  enddo

  x_g(:,:)  = 0.d0; x_s(:,:) = 0.d0; x_t(:,:) = 0.d0;
  y_s(:,:)  = 0.d0; y_t(:,:) = 0.d0;
  eq_s(:,:) = 0.d0; eq_t(:,:) = 0.d0; psi_eq(:,:) = 0.d0

  do i=1,n_vertex_max
    do j=1,n_order+1
      do ms=1, n_gauss
        do mt=1, n_gauss

          x_g(ms,mt) = x_g(ms,mt) + nodes(i)%x(1,j,1) * element%size(i,j) * H(i,j,ms,mt)
          x_s(ms,mt) = x_s(ms,mt) + nodes(i)%x(1,j,1) * element%size(i,j) * H_s(i,j,ms,mt)
          x_t(ms,mt) = x_t(ms,mt) + nodes(i)%x(1,j,1) * element%size(i,j) * H_t(i,j,ms,mt)
          y_s(ms,mt) = y_s(ms,mt) + nodes(i)%x(1,j,2) * element%size(i,j) * H_s(i,j,ms,mt)
          y_t(ms,mt) = y_t(ms,mt) + nodes(i)%x(1,j,2) * element%size(i,j) * H_t(i,j,ms,mt)

          psi_eq(ms,mt) = psi_eq(ms,mt) + nodes(i)%values(1,j,ivar_psi) * element%size(i,j) * H(i,j,ms,mt)
          
        enddo
      enddo
    enddo
  enddo

  do in=2,n_tor

    eq_s(:,:) = 0.d0; eq_t(:,:) = 0.d0

    do ms=1, n_gauss
      do mt=1, n_gauss

        do i=1,n_vertex_max
          do j=1,n_order+1

            eq_s(ms,mt)  = eq_s(ms,mt)  + nodes(i)%values(in,j,ivar_u) * element%size(i,j) * H_s(i,j,ms,mt)
            eq_t(ms,mt)  = eq_t(ms,mt)  + nodes(i)%values(in,j,ivar_u) * element%size(i,j) * H_t(i,j,ms,mt)
        
          enddo
        enddo

      enddo
    enddo

    do ms=1, n_gauss
      do mt=1, n_gauss
          
        if (.not. ((psi_eq(ms,mt) .gt. psi_start) .and. (psi_eq(ms,mt) .lt. psi_end))) cycle

        wst = wgauss(ms)*wgauss(mt)

        xjac = x_s(ms,mt)*y_t(ms,mt) - x_t(ms,mt)*y_s(ms,mt)
        BigR = x_g(ms,mt)

        u0_x  = (   y_t(ms,mt) * eq_s(ms,mt) - y_s(ms,mt) * eq_t(ms,mt) ) / xjac
        u0_y  = ( - x_t(ms,mt) * eq_s(ms,mt) + x_s(ms,mt) * eq_t(ms,mt) ) / xjac

        W_kin(in) = W_kin(in) + (u0_x*u0_x   + u0_y*u0_y) * BigR**3 * xjac * wst        

      enddo
    enddo

  enddo  ! n_tor loop

enddo    ! elements loop
!$omp end parallel do

do in=1,n_tor
  if (mode(in) .ne. 0) then
    W_kin(in) = 0.5d0 * W_kin(in)
  endif
enddo

return
end

end program particle_diagno

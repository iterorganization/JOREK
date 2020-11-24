module mod_particle_loop
use phys_module
use constants,   only: MU_ZERO, MASS_PROTON
use mpi
use particle_tracer
!$ use omp_lib

implicit none
private 
public loop_particle_kinetic

contains

subroutine loop_particle_kinetic(sim, jorek_feedback, rng, timesteps, n_steps, particle_start_time)

  use mod_project_particles
  use mod_random_seed
  use mod_interp, only: mode_moivre
  use mod_basisfunctions
  use mod_particle_sputtering, only: sample_fluid_particle_energy
  use mod_particle_types, only: copy_particle_kinetic_leapfrog

implicit none

class(particle_sim), target, intent(inout)                :: sim
type(projection), target, intent(inout)                   :: jorek_feedback
type(count_action)                                        :: counter
type(pcg32_rng), dimension(:), allocatable, intent(inout) :: rng
type(particle_kinetic_leapfrog)                           :: particle_tmp

real*8, intent(in)     :: timesteps, particle_start_time 

real*8,allocatable :: feedback_rhs(:,:,:,:,:)

real*8    :: n_norm, rho_norm, t_norm, v_norm, E_norm, M_norm
real*8    :: t, E(3), B(3), psi, U, n_e, T_e, rz_old(2), st_old(2)
real*8    :: v_temp(3), T_eV, K_eV, v_kin_temp, B_norm(3), v
real*8    :: R_g, Z_g, R_s, R_t, Z_s, Z_t, xjac, HZ(n_tor), HH(4,4), HH_s(4,4), HH_t(4,4)
real*8    :: ion_rate, ion_source, ion_prob, ion_ran(1), cx_ran(7), cx_source, cx_energy 
real*8    :: n_lost_ion, n_lost_ion_all
real*8    :: cx_prob, CX_rate
real*8    :: kinetic_energy, ion_energy
real*8    :: particle_source, velocity_par_source, energy_source
real*8    :: density_tot, density_in, density_out,  pressure, pressure_in, pressure_out
real*8    :: mom_par_tot, mom_par_in, mom_par_out, kin_par_tot, kin_par_out, kin_par_in
real*8    :: particles_remaining, momentum_remaining, energy_remaining
!$ real*8 :: w0, w1, mmm(3)

integer, intent(in)   :: n_steps
integer   :: i, j, k, l, m, i_elm_old, i_elm
integer   :: seed, i_rng, n_stream, ierr, nthreads
integer   :: i_tor, index_lm, i_elm_temp
integer   :: n_particles,n_reflect,ifail

!$ w0 = omp_get_wtime()

! step_rest_time is the amount of time that needs to be covered in addition to target_time - sim%time,
! i.e. (sim%time - step_rest_time(i)) is the 'real' particle time before the next loop

if (use_ccs .and. use_pcs) write(*,*) 'ERROR: pcs and ccs should not be used simultaneously'

n_norm   = CENTRAL_DENSITY * 1.d20                              ! (number) density normalisation
rho_norm = CENTRAL_MASS * MASS_PROTON * n_norm                  ! rho_SI = rho_norm * rho
t_norm   = sqrt((MU_ZERO * rho_norm))                           ! t_SI   = t_norm * t_jorek
v_norm   = 1.d0 / t_norm                                        ! V_SI   = v_norm * v_jorek
E_norm   = 1.5d0 / MU_ZERO                                      ! E_SI   = E_norm * E_jorek
M_norm   = rho_norm * v_norm                                    ! momentum normalisation

n_lost_ion = 0.d0
n_lost_ion_all = 0.d0

jorek_feedback%rhs_gather_time = jorek_feedback%rhs_gather_time + n_steps * timesteps

allocate(feedback_rhs,source=jorek_feedback%rhs)

jorek_feedback%rhs = 0.d0
feedback_rhs       = 0.d0

call with(sim, counter)

select type (particles => sim%groups(1)%particles)
type is (particle_kinetic_leapfrog)
if (sim%my_id .eq. 0) then
  write(*,*) "use_ionisation:         ", use_ionisation
  write(*,*) "use_cx:                 ", use_cx
  write(*,*) "use_sputtering:         ", use_sputtering
  write(*,*) "use_ncs:                ", use_ncs
  write(*,*) "use_ccs:                ", use_ccs
  write(*,*) "use_pcs:                ", use_pcs
  write(*,*) "timesteps:              ", timesteps 
  write(*,*) "---------------------------------"
  write(*,*) "Normalization:"
  write(*,'(A,e14.6)') ' N_norm   : ',N_norm
  write(*,'(A,e14.6)') ' rho_norm : ',rho_norm
  write(*,'(A,e14.6)') ' t_norm   : ',t_norm
  write(*,'(A,e14.6)') ' V_norm   : ',v_norm
  write(*,'(A,e14.6)') ' M_norm   : ',M_norm
  write(*,'(A,e14.6)') ' E_norm   : ',E_norm
endif
#ifdef __GFORTRAN__
      !$omp parallel do default(shared) & ! workaround for Error: ‘__vtab_mod_pcg32_rng_Pcg32_rng’ not specified in enclosing ‘parallel’
#else
      !$omp parallel do default(none) &
#endif
   !$omp schedule(dynamic,10) &
   !$omp shared(sim, particles, n_steps, timesteps, rng, particle_start_time, &
   !$omp rho_norm, t_norm, v_norm, E_norm, M_norm, N_norm, &
   !$omp use_cx, use_ionisation, use_sputtering, use_ncs, use_ccs, use_pcs, &
   !$omp jorek_feedback, CENTRAL_DENSITY, CENTRAL_MASS) &
   !$omp private(particle_tmp, i_rng, i,j,k,l,m, t, E, B, psi, U, rz_old, st_old, &
   !$omp i_elm_old, i_elm, n_e, T_e, ion_rate, ion_prob, ion_ran, ion_source, ion_energy, kinetic_energy,& 
   !$omp R_g, R_s, R_t, Z_g, Z_s, Z_t, xjac, HH, HH_s, HH_t, HZ, index_lm, &
   !$omp ifail, CX_rate, CX_prob, CX_source, CX_energy, v, &
   !$omp particle_source, velocity_par_source, energy_source, v_temp, K_eV, T_eV, cx_ran) &
   !$omp reduction(+:feedback_rhs)

   do j=1,size(particles,1)

      call copy_particle_kinetic_leapfrog(particles(j),particle_tmp)
!      i_rng = 1
    !$ i_rng = omp_get_thread_num()+1

      do k=1,n_steps

        if (particle_tmp%i_elm .le. 0) exit

        t = particle_start_time + (k-1)*timesteps

        call sim%fields%calc_EBpsiU(t, particle_tmp%i_elm, particle_tmp%st, particle_tmp%x(3), E, B, psi, U)

        rz_old    = particle_tmp%x(1:2)
        st_old    = particle_tmp%st
        i_elm_old = particle_tmp%i_elm
        
        ion_source = 0.d0
        ion_energy = 0.d0


       !< Ionisation block

       if (use_ionisation) then

          call sim%fields%calc_NeTe(t, particle_tmp%i_elm, particle_tmp%st, particle_tmp%x(3), n_e, T_e)

       
          call sim%groups(1)%ad%SCD%interp(int(particle_tmp%q), log10(n_e), log10(T_e), ion_rate) ! [m^3/s]
        
          ion_prob = 1.d0 - exp(-ion_rate * n_e * timesteps) ! [0] poisson point process, exponential 

          ! If the weight is to small throw away the particle with the probability, else reduce weight with ionising probability
          ion_source = 0.d0

          if (particle_tmp%weight .le. 1.0d7) then

            call rng(i_rng)%next(ion_ran)

            if (ion_ran(1) .le. ion_prob) then
              particle_tmp%i_elm  = 0
              ion_source = particle_tmp%weight
            else
              ion_source = 0.d0
            endif

          else 
            ion_source = particle_tmp%weight * ion_prob
            particle_tmp%weight = particle_tmp%weight * (1.d0 - ion_prob)
          endif 

          kinetic_energy = dot_product(particle_tmp%v,particle_tmp%v) *sim%groups(1)%mass * ATOMIC_MASS_UNIT /2.d0

          ion_energy     = kinetic_energy !- binding_energy

        endif ! use_ionisation

        !< Charge Exchange
        !< It is assumed that we will have a exchange between hydrogen isotopes

        v_temp    = particle_tmp%v
        cx_source = 0.d0
        cx_energy = 0.d0

        if (use_cx) then
  
          call sim%groups(1)%ad%CCD%interp(int(particle_tmp%q+1), log10(n_e), log10(T_e), CX_rate) ! [m^3/s]

          CX_prob = 1.d0 - exp(-CX_rate * n_e * timesteps)

          call rng(i_rng)%next(cx_ran)
 
          if (cx_ran(1) .le. CX_prob) then

            ! sample boltzman, randomize velocity
            T_eV = T_e * K_BOLTZ / EL_CHG

            ! sample from main plasma (should this not be a shifted Maxwellian?)

            call sample_fluid_particle_energy(T_eV, cx_ran(2:4), 1, K_eV) ! K_eV in eV. 

!THIS IS WRONG: use sample distorted maxwellian (or box-Mueller transform)
            v_temp    = sqrt(2.d0* K_eV *EL_CHG/(sim%groups(1)%mass * ATOMIC_MASS_UNIT)) * cx_ran(5:7)

            CX_source = particle_tmp%weight

            CX_energy   = 0.5d0 * sim%groups(1)%mass * ATOMIC_MASS_UNIT *  (dot_product(particle_tmp%v,particle_tmp%v) - dot_product(v_temp,v_temp))

          endif

        endif        ! use charge exchange

        energy_source = 0.d0
        particle_source = 0.d0 
        velocity_par_source = 0.d0

        if (use_ionisation .or. use_cx) then
          energy_source       = ion_source * ion_energy + cx_source * cx_energy
          
          particle_source     = ion_source * sim%groups(1)%mass * ATOMIC_MASS_UNIT

          velocity_par_source = ion_source * dot_product(B, particle_tmp%v)          * sim%groups(1)%mass * ATOMIC_MASS_UNIT &
               
               + cx_source  * dot_product(B, particle_tmp%v - v_temp) * sim%groups(1)%mass * ATOMIC_MASS_UNIT 
        
        endif        ! particle sources                               

        particle_tmp%v = v_temp 
       
        ! Push the particle and determine it's new location.
      
        if (particle_tmp%i_elm .gt. 0) then

          call boris_push_cylindrical(particle_tmp, sim%groups(1)%mass, E, B, timesteps)

          call find_RZ_nearby(sim%fields%node_list, sim%fields%element_list, rz_old(1), rz_old(2), st_old(1), st_old(2), i_elm_old, &
                              particle_tmp%x(1), particle_tmp%x(2), particle_tmp%st(1), particle_tmp%st(2), particle_tmp%i_elm, ifail)
        end if


        ! Calculate the projection of the ion source in real-time

        call basisfunctions(particle_tmp%st(1), particle_tmp%st(2), HH, HH_s, HH_t)
        call mode_moivre(particle_tmp%x(3), HZ)
        
        HZ(1) = HZ(1)*0.5d0 ! int cos^2(nx) from 0 to 2pi = pi for n > 0
        HZ(:) = HZ(:)/PI ! int 1 from 0 to 2pi = 2pi

        if (use_ncs) then
           do l=1,n_vertex_max
              do m=1,n_order+1

                 index_lm = (l-1)*(n_order+1) + m
                 
                 v   = HH(l,m) * sim%fields%element_list%element(i_elm_old)%size(l,m) 
                 
                 do i_tor=1,n_tor
                   feedback_rhs(m,l,i_elm,i_tor,1) = feedback_rhs(m,l,i_elm,i_tor,1) + HZ(i_tor) * v * particle_source     * t_norm / rho_norm
                   feedback_rhs(m,l,i_elm,i_tor,2) = feedback_rhs(m,l,i_elm,i_tor,2) + HZ(i_tor) * v * energy_source       * t_norm / E_norm
                   feedback_rhs(m,l,i_elm,i_tor,3) = feedback_rhs(m,l,i_elm,i_tor,3) + HZ(i_tor) * v * velocity_par_source * t_norm / m_norm
                 enddo
              enddo   !< order
           enddo      !< vertex
        end if        !< use_ncs


      end do ! steps

      if (use_ncs) then
         feedback_rhs(:,:,:,:,1:3) = feedback_rhs(:,:,:,:,1:3)/jorek_feedback%rhs_gather_time
      end if

      call copy_particle_kinetic_leapfrog(particle_tmp, particles(j))

      i_elm = particle_tmp%i_elm


      !< RHS for the ccs in jorek_feedback should be only filled in once, not accumulated over the particle stepping time like for the sources in the ncs.
      if (i_elm.gt.0) then

        call basisfunctions(particle_tmp%st(1), particle_tmp%st(2), HH, HH_s, HH_t)
        call mode_moivre(particle_tmp%x(3), HZ)
        
        HZ(1) = HZ(1)*0.5d0 ! int cos^2(nx) from 0 to 2pi = pi for n > 0
        HZ(:) = HZ(:)/PI    ! int 1 from 0 to 2pi = 2pi
       
        do l=1,n_vertex_max
          do m=1,n_order+1
             index_lm = (l-1)*(n_order+1) + m
             v   = HH(l,m) * sim%fields%element_list%element(i_elm)%size(l,m) 

             !< Ccs cannot be used in the same time as ncs for now
             if (use_pcs) then
                do i_tor=1,n_tor
                   feedback_rhs(m,l,i_elm,i_tor,4) = feedback_rhs(m,l,i_elm,i_tor,4) &
                          + HZ(i_tor) * v * dot_product(particle_tmp%v,particle_tmp%v) * sim%groups(1)%mass * mass_proton * particle_tmp%weight * mu_zero / 3.d0

                   feedback_rhs(m,l,i_elm,i_tor,5) = feedback_rhs(m,l,i_elm,i_tor,5) &
                          + HZ(i_tor) * v * el_chg * particle_tmp%q * particle_tmp%weight * particle_tmp%v(3) * mu_zero
                enddo
             !< Pcs cannot be used in the same time as ncs for now, and should never be used together with ccs
             elseif (use_ccs) then
                do i_tor=1,n_tor
                   feedback_rhs(m,l,i_elm,i_tor,4) = feedback_rhs(m,l,i_elm,i_tor,4) &
                          + HZ(i_tor) * v * el_chg * particle_tmp%q * particle_tmp%weight * sqrt(mu_zero / rho_norm)   

                   feedback_rhs(m,l,i_elm,i_tor,5) = feedback_rhs(m,l,i_elm,i_tor,5) &
                           + HZ(i_tor) * v * el_chg * particle_tmp%q * particle_tmp%weight * particle_tmp%v(1) * mu_zero

                   feedback_rhs(m,l,i_elm,i_tor,6) = feedback_rhs(m,l,i_elm,i_tor,6) &
                           + HZ(i_tor) * v * el_chg * particle_tmp%q * particle_tmp%weight * particle_tmp%v(2) * mu_zero

                   feedback_rhs(m,l,i_elm,i_tor,7) = feedback_rhs(m,l,i_elm,i_tor,7) &
                           + HZ(i_tor) * v * el_chg * particle_tmp%q * particle_tmp%weight * particle_tmp%v(3) * mu_zero
                enddo
             endif  !< use_pcs/ccs

          enddo     !< order
        enddo       !< vertex
      endif         !< ielm
    end do   ! particles
    !$omp end parallel do

    if (sim%my_id .eq. 0) write(*,*) "End of the particle loop"
  end select

  jorek_feedback%rhs = feedback_rhs

  deallocate(feedback_rhs)

  ! if (use_ionisation) then
  !   call MPI_REDUCE(n_lost_ion, n_lost_ion_all, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  !   if (sim%my_id .eq. 0) write(*,*) " Lost particles at t due to ionisation: ", sim%time, n_lost_ion_all
  ! end if 

  !$ w1 = omp_get_wtime()
  !$ mmm = mpi_minmeanmax(w1-w0)
  !$ if (sim%my_id .eq. 0) write(*,"(f10.7,A,3f9.4,A)") sim%time, " Particle stepping complete in ", mmm, "s"

  ! call Integrals_3D(sim%my_id, sim%fields%node_list, sim%fields%element_list, density_tot, density_in, density_out, &
  !                   pressure, pressure_in, pressure_out, kin_par_tot, kin_par_in, kin_par_out, mom_par_tot, mom_par_in, mom_par_out)

  particles_remaining = 0.d0
  momentum_remaining  = 0.d0
  energy_remaining    = 0.d0

!   select type (particles => sim%groups(1)%particles)
!   type is (particle_kinetic_leapfrog)

!     !$omp parallel do default(none) &
!     !$omp reduction(+:particles_remaining, momentum_remaining, energy_remaining) &
!     !$omp shared(sim, particles) &
!     !$omp private(j, E, B, psi, U, B_norm)
!     do j=1,size(particles,1)

!       if (particles(j)%i_elm .le. 0) cycle

!       call sim%fields%calc_EBpsiU(sim%time , particles(j)%i_elm, particles(j)%st, particles(j)%x(3), E, B, psi, U)
!       B_norm = B/norm2(B)

!       particles_remaining = particles_remaining + particles(j)%weight
!       momentum_remaining  = momentum_remaining  + particles(j)%weight * dot_product(B_norm,particles(j)%v) *sim%groups(1)%mass * ATOMIC_MASS_UNIT
!       energy_remaining    = energy_remaining    + particles(j)%weight * dot_product(particles(j)%v,particles(j)%v) *sim%groups(1)%mass * ATOMIC_MASS_UNIT /2.d0
! !      energy_remaining    = energy_remaining    + particles(j)%weight * 2.18d-15

!     enddo
!     !omp end parallel do

!     write(*,'(A,3e16.8)') 'REMAINING (START) : ',particles_remaining, momentum_remaining, energy_remaining

!   end select

  ! write(*,'(A,126e16.8)') ' TOTAL : ',sim%time,density_tot+particles_remaining/1.d20, density_tot, particles_remaining/1.d20, &
  !                                   mom_par_tot+momentum_remaining, mom_par_tot, momentum_remaining, &
  !                                   pressure+kin_par_tot+energy_remaining, pressure, energy_remaining, kin_par_tot
end subroutine

end module mod_particle_loop

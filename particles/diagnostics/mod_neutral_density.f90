!> contains a routine to get the neutral density to the aux node list, without having to call particle_evolution_loop
module mod_neutral_density
    use mod_project_particles, only: projection
    use mod_interp, only: mode_moivre
    use mod_basisfunctions
    use particle_tracer
    use phys_module, only: n_part_groups
    !$ use omp_lib
  
    implicit none
    
    integer, parameter :: i_neutral_n=1

    private
    public :: get_neutral_density
contains

!> first calculates the neutral density from all particles with ncs, then projects it
subroutine get_neutral_density(sim, neutral_density_proj)
  implicit none

  class(particle_sim), target, intent(inout)                :: sim
  type(projection),    target, intent(inout)                :: neutral_density_proj
  
  real*8,allocatable :: feedback_rhs(:,:,:,:,:)
  type (particle_group),         pointer :: part_group
  
  type(particle_kinetic_leapfrog) :: particle_tmp
  real*8    :: HZ(n_tor), HH(4,4), HH_s(4,4), HH_t(4,4)
  integer   :: j, l, m, i_tor, group_num


  neutral_density_proj%rhs = 0.d0

  allocate(feedback_rhs,source=neutral_density_proj%rhs)
  feedback_rhs       = 0.d0
  
  do group_num=1, n_part_groups
    part_group => sim%groups(group_num)
  
    if(part_group%coupling_scheme /= "ncs") cycle

#ifdef __GFORTRAN__
    !$omp parallel do default(shared) & ! workaround for Error: «__vtab_mod_pcg32_rng_Pcg32_rng» not specified in enclosing «parallel»
#else
    !$omp parallel do default(none)   &
#endif
    !$omp schedule(runtime)           &
    !$omp shared(sim, group_num)      &
    !$omp private(particle_tmp, j, l, m, i_tor, HH, HH_s, HH_t, HZ)  &
    !$omp reduction(+:feedback_rhs)
      do j=1,size(sim%groups(group_num)%particles,1)
        particle_tmp = sim%groups(group_num)%particles(j)
        
        if (particle_tmp%i_elm .le. 0) cycle
    
        !> Calculate the projection of the ion source in real-time
        call basisfunctions(particle_tmp%st(1), particle_tmp%st(2), HH, HH_s, HH_t)
        call mode_moivre(particle_tmp%x(3), HZ)
              
        do l=1,n_vertex_max
          do m=1,n_order+1
            do i_tor=1,n_tor
              feedback_rhs(m,l,particle_tmp%i_elm,i_tor,i_neutral_n) = feedback_rhs(m,l,particle_tmp%i_elm,i_tor,i_neutral_n) + HZ(i_tor) * HH(l,m) * sim%fields%element_list%element(particle_tmp%i_elm)%size(l,m) *particle_tmp%weight
            enddo
          enddo
        enddo
      end do   ! particles
      !$omp end parallel do

      neutral_density_proj%rhs(:,:,:,:,i_neutral_n) = neutral_density_proj%rhs(:,:,:,:,i_neutral_n) + feedback_rhs(:,:,:,:,i_neutral_n)   !< neutral density 
  enddo

  call neutral_density_proj%do(sim)
    
end subroutine get_neutral_density

end module mod_neutral_density

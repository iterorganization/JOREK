module mod_sources



implicit none



contains



!> Determine the heat and particle sources at a given position.
subroutine sources(xpoint2, xcase2, Z, Z_xpoint, psi, psi_axis, psi_bnd, particle_source, heat_source)

use phys_module

implicit none

! --- Routine parameters.
logical, intent(in)   :: xpoint2
integer, intent(in)   :: xcase2
real*8,  intent(in)   :: Z
real*8,  intent(in)   :: Z_xpoint(2)
real*8,  intent(in)   :: psi
real*8,  intent(in)   :: psi_axis
real*8,  intent(in)   :: psi_bnd
real*8,  intent(out)  :: particle_source
real*8,  intent(out)  :: heat_source

! --- Local variables
real*8 :: psi_n

psi_n = (psi - psi_axis) / (psi_bnd - psi_axis)

if (xpoint2) then
  if ((Z .lt. Z_xpoint(1)) .and. (psi_n .lt. 1.d0) ) then
    psi_n = 2.d0 - psi_n
  endif
endif

particle_source = particlesource * (0.5d0 - 0.5d0*tanh((psi_n - particlesource_psin)/particlesource_sig)) &
    + edgeparticlesource * (0.5d0 + 0.5d0*tanh((psi_n - edgeparticlesource_psin)/edgeparticlesource_sig))
heat_source     = heatsource     * (0.5d0 - 0.5d0*tanh((psi_n - heatsource_psin    )/heatsource_sig    ))

return
end subroutine sources



end module mod_sources

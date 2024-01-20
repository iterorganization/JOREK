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

particle_source = particlesource * (0.5d0 - 0.5d0*tanh((psi_n - particlesource_psin)/particlesource_sig)) &
     + edgeparticlesource * (0.5d0 + 0.5d0*tanh((psi_n - edgeparticlesource_psin)/edgeparticlesource_sig))&
     + particlesource_gauss * exp(-(psi_n - particlesource_gauss_psin)**2/(particlesource_gauss_sig**2))

heat_source = 0.d0

return
end subroutine sources



end module mod_sources

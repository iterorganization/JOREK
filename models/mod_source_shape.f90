module mod_source_shape

  use constants
  implicit none

  contains

  !> Determine source shape
  pure function source_shape(R,Z,phi,ns_R,ns_Z,ns_phi,ns_radius,ns_deltaphi,&
            psi,ns_psi,ns_grad_psi,ns_delta_minor_rad)
#if _OPENMP >= 201511
    !$omp declare simd
#endif
    implicit none
    
    real*8, intent(in)  :: R, Z, phi                 ! position where the source is calculated
    real*8, intent(in)  :: ns_R, ns_Z, ns_phi        ! position of the shard
    real*8, intent(in)  :: ns_radius, ns_deltaphi    ! extent of the ablation cloud

    real*8, intent(in)  :: psi,ns_psi,ns_grad_psi,ns_delta_minor_rad ! additional variables for a poloidally elongated source

    real*8              :: source_shape
    real*8              :: radius, ns_pol_shape
    real*8              :: dphi,   ns_tor_shape
    real*8              :: dminrad,ns_minrad_shape

   ! A gaussian shape is chosen poloidally
    radius = sqrt((R-ns_R)**2 + (Z-ns_Z)**2)
    ns_pol_shape = exp(-(radius/ns_radius)**2.d0)  

    ! A gaussian shape is chosen toroidally
    dphi = abs(phi - ns_phi)
    if (dphi .gt. PI) dphi = 2*PI - dphi  
    ns_tor_shape = exp(-(dphi/ns_deltaphi)**2.d0)

    source_shape = ns_pol_shape * ns_tor_shape

    ! Optionally, a gaussian shape is chosen in the minor radius direction
    if (ns_delta_minor_rad .gt. 0.) then
       dminrad = abs(psi-ns_psi)/ns_grad_psi
       ns_minrad_shape = exp(-(dminrad/ns_delta_minor_rad)**2.d0)
       source_shape = source_shape * ns_minrad_shape
    end if

  end function source_shape

end module mod_source_shape

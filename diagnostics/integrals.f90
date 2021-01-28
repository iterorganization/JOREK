!> Determines some integrals over the JOREK computational domain to determine the total current etc.
subroutine integrals(node_list, element_list, R_axis, Z_axis, psi_axis, R_xpoint, Z_xpoint, psi_xpoint, psi_limit, &
  aminor, Bgeo, current, beta_p, beta_t, beta_n, density, density_in, density_out, pressure,       &
  pressure_in, pressure_out, heat_src_in, heat_src_out, part_src_in, part_src_out)
use constants
use mod_parameters
use data_structure
use gauss
use basis_at_gaussian
use phys_module
use domains
use corr_neg
#if (JOREK_MODEL == 500 || JOREK_MODEL == 555)
  use mod_neutral_source, only: neutral_source, total_n_particles, total_n_particles_inj, total_n_particles_inj_all
#endif
#if (JOREK_MODEL == 501 || JOREK_MODEL == 502)
  use mod_injection_source, only: inj_source, total_n_particles, total_n_particles_inj, total_n_particles_inj_all
  use mod_impurity, only: radiation_function, radiation_function_linear
#endif


implicit none

! --- Routine parameters
type(type_node_list),    intent(in)    :: node_list
type(type_element_list), intent(in)    :: element_list
real*8,                  intent(in)    :: R_axis
real*8,                  intent(in)    :: Z_axis
real*8,                  intent(in)    :: psi_axis
real*8,                  intent(in)    :: R_xpoint(2)
real*8,                  intent(in)    :: Z_xpoint(2)
real*8,                  intent(in)    :: psi_xpoint(2)
real*8,                  intent(in)    :: psi_limit
real*8,                  intent(in)    :: aminor
real*8,                  intent(out)   :: Bgeo
real*8,                  intent(out)   :: current
real*8,                  intent(out)   :: beta_p
real*8,                  intent(out)   :: beta_t
real*8,                  intent(out)   :: beta_n
real*8,                  intent(out)   :: density
real*8,                  intent(out)   :: density_in
real*8,                  intent(out)   :: density_out
real*8,                  intent(out)   :: pressure
real*8,                  intent(out)   :: pressure_in
real*8,                  intent(out)   :: pressure_out
real*8,                  intent(out)   :: heat_src_in
real*8,                  intent(out)   :: heat_src_out
real*8,                  intent(out)   :: part_src_in
real*8,                  intent(out)   :: part_src_out

! --- Local variables
type (type_element)      :: element
type (type_node)         :: nodes(n_vertex_max)
real*8     :: x_g(n_gauss,n_gauss), x_s(n_gauss,n_gauss), x_t(n_gauss,n_gauss)
real*8     :: y_g(n_gauss,n_gauss), y_s(n_gauss,n_gauss), y_t(n_gauss,n_gauss)
real*8     :: eq_g(0:n_var,n_gauss,n_gauss), eq_s(0:n_var,n_gauss,n_gauss), eq_t(0:n_var,n_gauss,n_gauss)
integer    :: i, j, k, in, ms, mt, iv, inode, ife, n_elements
real*8     :: xjac, BigR, wst, P_int, C_intern, ZJ_0, PS_0, Volume, Area
real*8     :: rho_00, T_00, Ti_00, Te_00, current_in, current_out, rhon_00 
real*8     :: C_hel, P_hel, D_int, D_ext, P_ext, C_ext
real*8     :: part_src, heat_src, heat_src_i, heat_src_e

real*8     :: r0_corr, T0_corr, rn0_corr
! Additional diagnostic variables for impurity model
! See https://www.jorek.eu/wiki/doku.php?id=model500_501_555 for details
#if (JOREK_MODEL == 501)
real*8  :: local_radiation, local_E_ion, total_radiation, total_E_ion
real*8  :: local_radiation_phi(n_plane), total_radiation_phi(n_plane)
! Atomic physics coefficients:
!   -Mass ratio between main ions and impurites (m_i/m_imp)
real*8  :: m_i_over_m_imp
!   -Mean impurity ionization state
real*8  :: Z_imp, T0_Zimp, alpha_Zimp, Z_eff, eta_coef, ne_JOREK
!   -Coefficients related to Z_imp
real*8  :: alpha_imp, beta_imp
!   -Corrected plasma temperature and density for radiation calculation
real*8  :: Te_corr_eV, ne_SI, Te_eV
!   -Temporary variable for charge state distribution
real*8, allocatable :: P_imp(:)
real*8     :: E_ion, Lrad, E_ion_bg
integer*8  :: ion_i, ion_k, i_phi
#endif

write(*,*) '***************************************'
write(*,*) '* Integrals                           *'
write(*,*) '***************************************'

density  = 0.d0
pressure = 0.d0
D_int    = 0.d0
P_int    = 0.d0
C_intern = 0.d0
D_ext    = 0.d0
P_ext    = 0.d0
C_ext    = 0.d0
P_hel    = 0.d0
C_hel    = 0.d0
Volume   = 0.d0
Area     = 0.d0
heat_src = 0.d0
part_src = 0.d0
heat_src_in  = 0.d0
heat_src_out = 0.d0
part_src_in  = 0.d0
part_src_out = 0.d0

Bgeo = F0 / R_geo

do ife =1, element_list%n_elements

  element = element_list%element(ife)

  do iv = 1, n_vertex_max
    inode     = element%vertex(iv)
    nodes(iv) = node_list%node(inode)
  enddo

  x_g(:,:)    = 0.d0; x_s(:,:)    = 0.d0; x_t(:,:)    = 0.d0;
  y_g(:,:)    = 0.d0; y_s(:,:)    = 0.d0; y_t(:,:)    = 0.d0;
  eq_g(:,:,:) = 0.d0; eq_s(:,:,:) = 0.d0; eq_t(:,:,:) = 0.d0;

  do i=1,n_vertex_max
    do j=1,n_order+1
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

  eq_g(:,:,:) = 0.d0; eq_s(:,:,:) = 0.d0; eq_t(:,:,:) = 0.d0;

  do i=1,n_vertex_max
    do j=1,n_order+1
      do ms=1, n_gauss
        do mt=1, n_gauss

          do k=1,n_var
            eq_g(k,ms,mt)  = eq_g(k,ms,mt)  + nodes(i)%values(1,j,k) * element%size(i,j) * H(i,j,ms,mt)
            eq_s(k,ms,mt)  = eq_s(k,ms,mt)  + nodes(i)%values(1,j,k) * element%size(i,j) * H_s(i,j,ms,mt)
            eq_t(k,ms,mt)  = eq_t(k,ms,mt)  + nodes(i)%values(1,j,k) * element%size(i,j) * H_t(i,j,ms,mt)
          enddo

        enddo
      enddo
    enddo
  enddo
!--------------------------------------------------- sum over the Gaussian integration points

  do ms=1, n_gauss

    do mt=1, n_gauss

      wst = wgauss(ms)*wgauss(mt)

      xjac = x_s(ms,mt)*y_t(ms,mt) - x_t(ms,mt)*y_s(ms,mt)
      BigR = x_g(ms,mt)

      rho_00 = eq_g(var_rho,ms,mt)
      if ( (jorek_model .eq. 400) .or. (jorek_model .eq. 502) .or. (jorek_model .eq. 711) ) then
        Ti_00 = eq_g(var_Ti,ms,mt)
        Te_00 = eq_g(var_Te,ms,mt)
        T_00  = Ti_00 + Te_00
      else
        T_00  = eq_g(var_T,ms,mt)
      endif
      ZJ_0  = eq_g(var_zj,ms,mt)
      PS_0  = eq_g(var_psi,ms,mt) 
      
#if (JOREK_MODEL == 500) || (JOREK_MODEL == 501) || (JOREK_MODEL == 555)
      rhon_00 = eq_g(var_rhon,ms,mt)
      rn0_corr = corr_neg_dens1(rhon_00)
#endif

#if (JOREK_MODEL == 501)
      r0_corr = corr_neg_dens1(rho_00)
      if (T_min > T_1) then
        T0_corr = corr_neg_temp(T_00,(/5.d-1,5.d-1/),2.*T_min) ! For use in eta(T), visco(T), ...
      else
        T0_corr = corr_neg_temp(T_00,(/5.d-1,5.d-1/),2.*T_1) ! For use in eta(T), visco(T), ...
      end if

      !-------------------------------------------
      ! Atomic physics parameters for Impurities
      !-------------------------------------------

      select case ( trim(imp_type) )
        case('D2')
          m_i_over_m_imp = central_mass/2.  ! Deuterium mass = 2 u
        case('Ar')
          m_i_over_m_imp = central_mass/40. ! Argon mass = 40 u
        case('Ne')
          m_i_over_m_imp = central_mass/20. ! Neon mass = 20 u
        case default
          write(*,*) '!! Gas type "', trim(imp_type), '" unknown (in mod_injection_source.f90) !!'
          write(*,*) '=> We assume the gas is D2.'
          m_i_over_m_imp = central_mass/2.
      end select

      ! Te in eV:
      Te_corr_eV = T0_corr/(2.d0*EL_CHG*MU_ZERO*central_density*1.d20)
      Te_eV = T_00/(2.d0*EL_CHG*MU_ZERO*central_density*1.d20)
   
      if (allocated(imp_adas(1)%ionisation_energy)) then
   
        if (allocated(P_imp)) deallocate(P_imp)
 
        allocate(P_imp(0:imp_adas(1)%n_Z))
   
        call imp_cor(1)%interp_linear(density=20.,temperature=log10(Te_corr_eV*EL_CHG/K_BOLTZ),&
                                      p_out=P_imp,z_avg=Z_imp)
   
        ! Calculate the ionization potential energy and its derivative wrt. temperature
        E_ion     = 0.
        E_ion_bg  = 13.6
        do ion_i=1, imp_adas(1)%n_Z
          do ion_k=1, ion_i
            E_ion     = E_ion + P_imp(ion_i)*imp_adas(1)%ionisation_energy(ion_k)
          end do
        end do
        ! Convert from eV to SI unit
        E_ion     = E_ion * EL_CHG
        E_ion_bg  = E_ion_bg * EL_CHG
      else
        call imp_cor(1)%interp_linear(density=20.,temperature=log10(Te_corr_eV*EL_CHG/K_BOLTZ),z_avg=Z_imp)
        E_ion     = 0.
        E_ion_bg  = 0.
      end if

      alpha_imp    = 0.5*m_i_over_m_imp*(Z_imp+1.) - 1.
      beta_imp     = m_i_over_m_imp*Z_imp - 1.
      ne_SI        = (r0_corr + beta_imp * rn0_corr) * 1.d20 * central_density !electron density (SI)
      ne_JOREK     = r0_corr + beta_imp * rn0_corr ! Electron density in JOREK unit
      ne_JOREK     = corr_neg_dens(ne_JOREK,(/1.d-1,1.d-1/),1.d-3) ! Correction for negative electron density
                                                               ! Too small rho_1 will cause a problem

#endif

#if (JOREK_MODEL == 501)
      density  = density  + ((rho_00-rhon_00) + rhon_00*m_i_over_m_imp) * xjac * 2.d0 * PI * BigR * wst
#else
      density  = density  + rho_00       * xjac * 2.d0 * PI * BigR * wst
#endif
#if (JOREK_MODEL == 501)
      pressure = pressure + (rho_00 + alpha_imp*rhon_00) * T_00 * xjac * 2.d0 * PI * BigR * wst
#else
      pressure = pressure + rho_00 * T_00 * xjac * 2.d0 * PI * BigR * wst
#endif

      if ( in_plasma(node_list,element_list,x_g(ms,mt),y_g(ms,mt),eq_g(1,ms,mt),xpoint,&
        xcase,R_xpoint,Z_xpoint,psi_xpoint,psi_limit,R_axis,Z_axis,psi_axis) ) then
        
#if ( (JOREK_MODEL == 400) || (JOREK_MODEL == 401) || (JOREK_MODEL == 502) || (JOREK_MODEL == 711) )
        call sources(xpoint, xcase, y_g(ms,mt), Z_xpoint, eq_g(var_psi,ms,mt), psi_axis, &
          psi_limit, part_src, heat_src_i, heat_src_e)
        heat_src = heat_src_i + heat_src_e
#else
        call sources(xpoint, xcase, y_g(ms,mt), Z_xpoint, eq_g(var_psi,ms,mt), psi_axis, &
          psi_limit, part_src, heat_src)
#endif
        
#if (JOREK_MODEL == 501)
        ! --- 3D integrals
        D_int = D_int + ((rho_00-rhon_00) + rhon_00*m_i_over_m_imp) * xjac * 2.d0 * PI * BigR * wst
        P_int = P_int + (rho_00+alpha_imp*rhon_00) * T_00 * xjac * 2.d0 * PI * BigR * wst
        ! --- 2D integrals
        P_hel = P_hel + (rho_00+alpha_imp*rhon_00) * T_00 * xjac * wst
#else
        ! --- 3D integrals
        D_int = D_int + rho_00       * xjac * 2.d0 * PI * BigR * wst
        P_int = P_int + rho_00 * T_00 * xjac * 2.d0 * PI * BigR * wst
        ! --- 2D integrals
        P_hel = P_hel + rho_00 * T_00 * xjac * wst
#endif
        C_intern = C_intern + ZJ_0 /BigR  * xjac * 2.d0 * PI * BigR * wst
        Volume = Volume + 2.d0 * PI * BigR * xjac * wst
        heat_src_in = heat_src_in + 2.d0 * PI * BigR * xjac * wst * heat_src
        part_src_in = part_src_in + 2.d0 * PI * BigR * xjac * wst * part_src
        
        ! --- 2D integrals
        C_hel = C_hel + ZJ_0 /BigR  * xjac * wst
        Area   = Area   + xjac * wst
        
      else

#if (JOREK_MODEL == 501)
        D_ext = D_ext + ((rho_00-rhon_00) + rhon_00*m_i_over_m_imp) * xjac * 2.d0 * PI * BigR * wst
        P_ext = P_ext + (rho_00+alpha_imp*rhon_00) * T_00 * xjac * 2.d0 * PI * BigR * wst
#else
        D_ext = D_ext + rho_00       * xjac * 2.d0 * PI * BigR * wst
        P_ext = P_ext + rho_00 * T_00 * xjac * 2.d0 * PI * BigR * wst
#endif
        C_ext = C_ext + ZJ_0 /BigR  * xjac * 2.d0 * PI * BigR * wst
        
        heat_src_out = heat_src_out + 2.d0 * PI * BigR * xjac * wst * heat_src
        part_src_out = part_src_out + 2.d0 * PI * BigR * xjac * wst * part_src
      endif
      
    enddo
  enddo
enddo

density_in   = D_int
density_out  = D_ext
pressure_in  = P_int
pressure_out = P_ext
current_in   = -C_intern
current_out  = -C_ext

current = -C_hel / MU_zero
beta_p  = 8.d0 * PI * P_hel / (C_hel**2 )
beta_t  = 2.d0 * P_hel / Bgeo**2 / (Area)
beta_n  = 100.d0 * (4.*PI/10.) * beta_t / (MU_zero * abs(current) /  (aminor * Bgeo))

write(*,'(A,f16.7)')    ' psi_limit        : ',psi_limit
write(*,'(A,f16.7,A)')  ' current          : ',current/1.e6,' MA'
write(*,'(A,f16.7)')    ' beta_p           : ',beta_p
write(*,'(A,f16.7)')    ' beta_t           : ',beta_t
write(*,'(A,f16.7,A)')  ' beta_n           : ',beta_n,' [%]'
write(*,'(A,f16.7,A)')  ' Area             : ',area,' m^2'
write(*,'(A,f16.7,A)')  ' Volume           : ',volume,' m^3'
write(*,'(A,2es18.7,A)') ' Heat.src (in/out): ', heat_src_in, heat_src_out, '/ sqrt((mu_0)^3 rho_0) W'
write(*,'(A,2es18.7,A)') ' Part.src (in/out): ', part_src_in, part_src_out

write(*,'(A,5f10.5)') ' density  (total/in/out)  : ',density,  density_in,  density_out 
write(*,'(A,5f10.5)') ' pressure (total/in/out)  : ',pressure, pressure_in, pressure_out 
write(*,'(A,5f10.5)') ' current  (in/out)        : ',current_in, current_out 

return
end subroutine integrals

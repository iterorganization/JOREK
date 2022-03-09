subroutine Integrals_3D(my_id, node_list, element_list, density_tot, density_in, density_out, pressure, pressure_in, pressure_out, &
                        kin_par_tot, kin_par_in, kin_par_out, mom_par_tot, mom_par_in,mom_par_out)
!---------------------------------------------------------------
!
!---------------------------------------------------------------
use constants
use data_structure
use gauss
use basis_at_gaussian
use phys_module
use pellet_module
use mpi_mod
use domains
use corr_neg
use equil_info, only : get_psi_n, ES
!$ use omp_lib
#ifdef WITH_Neutrals
  use mod_neutral_source
#endif
#ifdef WITH_Impurities
  use mod_impurity
#endif
use mod_sources

implicit none

real*8, intent(out) :: density_tot, density_in, density_out,  pressure, pressure_in, pressure_out
real*8, intent(out) :: kin_par_tot, kin_par_in, kin_par_out, mom_par_tot, mom_par_in, mom_par_out

type (type_node_list)    :: node_list
type (type_element_list) :: element_list
type (type_element)      :: element
type (type_node)         :: nodes(n_vertex_max)

real*8  :: x_g(n_gauss,n_gauss),        x_s(n_gauss,n_gauss),        x_t(n_gauss,n_gauss)
real*8  :: y_g(n_gauss,n_gauss),        y_s(n_gauss,n_gauss),        y_t(n_gauss,n_gauss)
real*8  :: eq_g(n_plane,0:n_var,n_gauss,n_gauss), eq_s(n_plane,0:n_var,n_gauss,n_gauss)
real*8  :: eq_t(n_plane,0:n_var,n_gauss,n_gauss), eq_p(n_plane,0:n_var,n_gauss,n_gauss)
real*8  :: wgauss_copy(n_gauss)

real*8  :: particle_source, heat_source, heat_source_i, heat_source_e, xt, t_norm, rho_norm, rotation_source
real*8  :: eq_zne(n_gauss,n_gauss), eq_zTe(n_gauss,n_gauss)
real*8  :: dn_dpsi,dn_dz,dn_dpsi2,dn_dz2,dn_dpsi_dz,dn_dpsi3,dn_dpsi_dz2, dn_dpsi2_dz
real*8  :: dT_dpsi,dT_dz,dT_dpsi2,dT_dz2,dT_dpsi_dz,dT_dpsi3,dT_dpsi_dz2, dT_dpsi2_dz

integer :: i, j, k, in, ms, mt, mp, iv, inode, ife, n_elements, ifail
integer :: ierr, n_cpu, my_id, ife_delta, ife_min, ife_max, omp_nthreads, omp_tid
real*8  :: beta_p, beta_n, beta_t, aminor
real*8  :: xjac, BigR, wst, P_int, C_intern, zj0, ps0, r0, T0, T0e, Vol, Volume, Area, Bgeo, psi_limit
real*8  :: r0_corr, T0_corr

real*8  :: current_in, current_out, D_int, D_ext, P_ext, C_ext, P_max, delta_phi, phi, P_tot, D_tot
real*8  :: VP_int, VP_ext,TVP_int, TVP_ext, VK_int, VK_ext, vpar0, BB2, VP_tot, TVP_tot, VK_tot
real*8  :: kin_perp_in, kin_perp_out, kin_perp_tot
real*8  :: VM_int, VM_ext, VM_tot, mag_in, mag_out, mag_tot, J2_int, J2_ext, J2_tot, ohm_in, ohm_tot, ohm_out
real*8  :: H_int, H_ext, S_int, S_ext, heating_in, heating_out, source_in, source_out
real*8  :: dTdx, dTdy, drhodx, drhody, dPdx, dPdy, dpsidx, dpsidy, dudx, dudy
real*8  :: grad_psi, grad_P, grad_P_psi, gradP_psi_max, gradP_max
real*8  :: source_volume, source_pellet, eta_T_ohm
real*8  :: local_pellet_particles, local_plasma_particles, local_pellet_volume
real*8  :: local_n_particles_inj, local_n_particles, source_neutral, rn0, rho_bar

integer    :: spi_i


call MPI_COMM_SIZE(MPI_COMM_WORLD, n_cpu, ierr) ! number of MPI procs

n_cpu = max(n_cpu,1)

if (my_id .eq. 0) then
  write(*,*) '***************************************'
  write(*,*) '* Integrals  (3D)                     *'
  write(*,*) '***************************************'
  write(*,*) ' n_plane : ',n_plane
  write(*,*) ' n_cpu   : ',n_cpu
endif

wgauss_copy = wgauss

density_tot  = 0.d0
pressure = 0.d0
D_int    = 0.d0
P_int    = 0.d0
C_intern = 0.d0
H_int    = 0.d0
S_int    = 0.d0
VP_int   = 0.d0
VK_int   = 0.d0
VM_int   = 0.d0
J2_int   = 0.d0
TVP_int  = 0.d0
D_ext    = 0.d0
P_ext    = 0.d0
C_ext    = 0.d0
H_ext    = 0.d0
S_ext    = 0.d0
VP_ext   = 0.d0
VK_ext   = 0.d0
VM_ext   = 0.d0
J2_ext   = 0.d0
TVP_ext  = 0.d0
Vol      = 0.d0
P_tot    = 0.d0
D_tot    = 0.d0
VP_tot   = 0.d0
VK_tot   = 0.d0
VM_tot   = 0.d0
J2_tot   = 0.d0
TVP_tot  = 0.d0

local_pellet_particles = 0.d0
local_plasma_particles = 0.d0
local_pellet_volume    = 0.d0

#if (defined WITH_Neutrals) && (!defined WITH_Impurities)
local_n_particles_inj = 0.d0
local_n_particles     = 0.d0
#endif

Bgeo = F0 / R_geo

delta_phi = 2.d0 * PI / float(n_plane) / float(n_period)

P_max         = 0.d0
gradP_max     = 0.d0
gradP_psi_max = 0.d0

psi_limit = ES%psi_bnd

ife_delta = ceiling(float(element_list%n_elements) / n_cpu)
ife_min   =      my_id     * ife_delta + 1
ife_max   = min((my_id +1) * ife_delta, element_list%n_elements)

!$omp parallel default(none)                                                                   &
!$omp   shared(element_list,node_list, H, H_s, H_t, HZ, HZ_p, ife_min, ife_max, xpoint, xcase, &
!$omp          ES, my_id, use_pellet, psi_limit, delta_phi,                                    &
!$omp          D_tot, D_int, D_Ext, P_tot, P_int, P_ext, Vol, C_intern, C_ext, VP_ext, VP_int, &
!$omp          VK_ext, VK_int, VK_tot, VM_ext, VM_int, VM_tot, J2_tot, J2_ext, J2_int,         &
!$omp          TVP_int, TVP_ext, TVP_tot,                                                      &
!$omp          H_int, H_ext, S_int, S_ext, F0, VP_tot, eta_ohmic, eta_T_dependent,             &
!$omp          T_max_eta_ohm, T_0,                                                             &
!$omp          pellet_amplitude,pellet_R,pellet_Z,pellet_psi,pellet_phi,                       &
!$omp          pellet_radius, pellet_delta_psi, pellet_sig, pellet_length, pellet_ellipse, pellet_theta,  &
!$omp          central_density, pellet_particles,pellet_density, pellet_volume,                &
!$omp          local_pellet_particles, local_plasma_particles, local_pellet_volume,            &
#if (defined WITH_Neutrals) && (!defined WITH_Impurities)
!$omp          local_n_particles_inj, local_n_particles, ns_amplitude, ns_R, ns_Z,             &
!$omp          ns_phi, ns_radius, ns_deltaphi, ns_tor_norm, spi_tor_rot,                       &
!$omp          t_now, A_Dmv, K_Dmv, V_Dmv, P_Dmv, t_ns, L_tube, JET_MGI,ASDEX_MGI,             &
!$omp          central_mass, pellets, tor_frequency,                                           &
!$omp          ns_radius_ratio, ns_radius_min, spi_shard_file,                                 &
#endif
!$omp          wgauss_copy)                                                                    &
!$omp   private(ife,iv,inode,element,nodes,i,j, k,in, mp, ms, mt, spi_i,                       &
!$omp           x_g, y_g, x_s, y_s, x_t, y_t, xjac, eq_g, eq_s, eq_t, eq_p,                    &
!$omp           wst, BigR, r0, T0, T0e, zj0, ps0, dTdx, dTdy, drhodx, drhody, dpsidx, dpsidy, dudx, dudy,  &
!$omp           dpdx, dpdy, grad_P, grad_psi, grad_P_psi,gradP_max, gradP_psi_max, phi,        &
!$omp           P_max, source_pellet, source_volume, eq_zne, eq_zTe, vpar0, BB2, eta_T_ohm,    &
!$omp           heat_source, heat_source_i, heat_source_e, particle_source, rotation_source,   &
!$omp           dn_dpsi,dn_dz,dn_dpsi2,dn_dz2,dn_dpsi_dz,dn_dpsi3,dn_dpsi_dz2, dn_dpsi2_dz,    &
!$omp           dT_dpsi,dT_dz,dT_dpsi2,dT_dz2,dT_dpsi_dz,dT_dpsi3,dT_dpsi_dz2, dT_dpsi2_dz,    &
!$omp           r0_corr, T0_corr,                                                              &
#if (defined WITH_Neutrals) && (!defined WITH_Impurities)
!$omp           rn0, source_neutral,                                                           &
#endif
!$omp           omp_nthreads,omp_tid)


#ifdef OPENMP
omp_nthreads = omp_get_num_threads()
omp_tid      = omp_get_thread_num()
#else
omp_nthreads = 1
omp_tid      = 0
#endif

!$omp do reduction(+:local_pellet_particles, local_plasma_particles, local_pellet_volume,     &
#ifdef WITH_Neutrals
!$omp                local_n_particles_inj,  local_n_particles,                               &
#endif
!$omp                D_int, D_ext, P_int, H_int, S_int, H_ext, S_ext, P_ext, C_intern, C_ext, &
!$omp                TVP_int, TVP_ext, TVP_tot, VP_int, VP_ext, VP_tot, VK_tot, VK_int, VK_ext, VM_ext,                  &
!$omp                VM_int, VM_tot, Vol, P_tot, D_tot,J2_tot, J2_int, J2_ext)

do ife = ife_min, ife_max

  element = element_list%element(ife)

  do iv = 1, n_vertex_max
    inode     = element%vertex(iv)
    nodes(iv) = node_list%node(inode)
  enddo

  x_g(:,:)    = 0.d0; x_s(:,:)    = 0.d0; x_t(:,:)    = 0.d0;
  y_g(:,:)    = 0.d0; y_s(:,:)    = 0.d0; y_t(:,:)    = 0.d0;

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

  eq_g(:,:,:,:) = 0.d0; eq_s(:,:,:,:) = 0.d0; eq_t(:,:,:,:) = 0.d0; eq_p(:,:,:,:) = 0.d0;

  do i=1,n_vertex_max
    do j=1,n_order+1

      do mp=1,n_plane
        do ms=1, n_gauss
          do mt=1, n_gauss

            do k=1,n_var
              do in=1,n_tor
                eq_g(mp,k,ms,mt) = eq_g(mp,k,ms,mt) + nodes(i)%values(in,j,k) * element%size(i,j) * H(i,j,ms,mt)  * HZ(in,mp)
                eq_s(mp,k,ms,mt) = eq_s(mp,k,ms,mt) + nodes(i)%values(in,j,k) * element%size(i,j) * H_s(i,j,ms,mt)* HZ(in,mp)
                eq_t(mp,k,ms,mt) = eq_t(mp,k,ms,mt) + nodes(i)%values(in,j,k) * element%size(i,j) * H_t(i,j,ms,mt)* HZ(in,mp)
                eq_p(mp,k,ms,mt) = eq_p(mp,k,ms,mt) + nodes(i)%values(in,j,k) * element%size(i,j) * H(i,j,ms,mt)  * HZ_p(in,mp)
              enddo
            enddo

	  enddo
        enddo
      enddo

    enddo
  enddo

  do ms=1, n_gauss
    do mt=1, n_gauss

      call density(xpoint, xcase, y_g(ms,mt), ES%Z_xpoint, eq_g(1,1,ms,mt),ES%psi_axis,ES%psi_bnd,eq_zne(ms,mt), &
                   dn_dpsi,dn_dz,dn_dpsi2,dn_dz2,dn_dpsi_dz,dn_dpsi3,dn_dpsi_dz2, dn_dpsi2_dz)

#ifdef WITH_TiTe
      call temperature_e(xpoint, xcase, y_g(ms,mt), ES%Z_xpoint, eq_g(1,1,ms,mt),ES%psi_axis,ES%psi_bnd,eq_zTe(ms,mt), &
                       dT_dpsi,dT_dz,dT_dpsi2,dT_dz2,dT_dpsi_dz,dT_dpsi3,dT_dpsi_dz2, dT_dpsi2_dz)
#else
      call temperature(xpoint, xcase, y_g(ms,mt), ES%Z_xpoint, eq_g(1,1,ms,mt),ES%psi_axis,ES%psi_bnd,eq_zTe(ms,mt), &
                       dT_dpsi,dT_dz,dT_dpsi2,dT_dz2,dT_dpsi_dz,dT_dpsi3,dT_dpsi_dz2, dT_dpsi2_dz)
#endif

    enddo
  enddo

#ifndef WITH_TiTe
  eq_zTe = eq_zTe / 2.d0	! electron temperature
#endif
!--------------------------------------------------- sum over the Gaussian integration points

  do mp=1,n_plane

    phi       = 2.d0*PI*float(mp-1)/float(n_plane) / float(n_period)

    do ms=1, n_gauss

      do mt=1, n_gauss

        wst = wgauss_copy(ms)*wgauss_copy(mt)

        xjac = x_s(ms,mt)*y_t(ms,mt) - x_t(ms,mt)*y_s(ms,mt)
        BigR = x_g(ms,mt)

        r0      = eq_g(mp,var_rho,ms,mt)
        r0_corr = corr_neg_dens1(r0)
#ifdef WITH_TiTe
        T0      = eq_g(mp,var_Ti,ms,mt)
        T0_corr = corr_neg_temp1(T0)
        T0e     = corr_neg_temp1(eq_g(mp,var_Te,ms,mt))
#else
        T0      = eq_g(mp,var_T,ms,mt)
        T0_corr = corr_neg_temp1(T0)
        T0e     = eq_g(mp,var_T,ms,mt) /2.d0
#endif
        zj0    = eq_g(mp,var_zj,ms,mt)
        ps0    = eq_g(mp,var_psi,ms,mt)

#ifdef WITH_Vpar
        vpar0 = eq_g(mp,var_Vpar,ms,mt)
#else
        vpar0 = 0.d0
#endif

#if (defined WITH_Neutrals) && (!defined WITH_Impurities)
        rn0    = eq_g(mp,var_rhon,ms,mt)
#endif

        ! --- Eta for ohmic heating
        if ( eta_T_dependent .and. T0_corr <= T_max_eta_ohm ) then
          eta_T_ohm     = eta_ohmic   * (T0_corr/T_0)**(-1.5d0)
        else if ( eta_T_dependent .and. T0_corr > T_max_eta_ohm ) then
          eta_T_ohm     = eta_ohmic   * (T_max_eta_ohm/T_0)**(-1.5d0)
        else
          eta_T_ohm     = eta_ohmic
        end if

        dTdx   = (   y_t(ms,mt) * eq_s(mp,6,ms,mt) - y_s(ms,mt) * eq_t(mp,6,ms,mt) ) / xjac
        dTdy   = ( - x_t(ms,mt) * eq_s(mp,6,ms,mt) + x_s(ms,mt) * eq_t(mp,6,ms,mt) ) / xjac
        drhodx = (   y_t(ms,mt) * eq_s(mp,5,ms,mt) - y_s(ms,mt) * eq_t(mp,5,ms,mt) ) / xjac
        drhody = ( - x_t(ms,mt) * eq_s(mp,5,ms,mt) + x_s(ms,mt) * eq_t(mp,5,ms,mt) ) / xjac

        dpsidx = (   y_t(ms,mt) * eq_s(mp,1,ms,mt) - y_s(ms,mt) * eq_t(mp,1,ms,mt) ) / xjac
        dpsidy = ( - x_t(ms,mt) * eq_s(mp,1,ms,mt) + x_s(ms,mt) * eq_t(mp,1,ms,mt) ) / xjac

        dudx = (   y_t(ms,mt) * eq_s(mp,2,ms,mt) - y_s(ms,mt) * eq_t(mp,2,ms,mt) ) / xjac
        dudy = ( - x_t(ms,mt) * eq_s(mp,2,ms,mt) + x_s(ms,mt) * eq_t(mp,2,ms,mt) ) / xjac

        BB2 = (F0*F0 + dpsidx*dpsidx + dpsidy*dpsidy) / BigR**2

        dPdx = r0 * dTdx + T0 * drhodx
        dPdy = r0 * dTdy + T0 * drhody

        grad_P   = sqrt(dPdx**2   + dPdy**2)
        grad_psi = sqrt(dpsidx**2 + dpsidy**2)

        grad_P_psi = (dPdx * dpsidx + dPdy * dpsidy)/grad_psi

#ifdef WITH_TiTe
        call sources(xpoint, xcase, y_g(ms,mt), ES%Z_xpoint, ps0, ES%psi_axis, psi_limit, &
                     particle_source,heat_source_i,heat_source_e)
        heat_source = heat_source_i + heat_source_e
#else
        call sources(xpoint, xcase, y_g(ms,mt), ES%Z_xpoint, ps0, ES%psi_axis, psi_limit, &
                     particle_source,heat_source)
#endif

        P_tot  = P_tot  + r0 * T0 * xjac * BigR * wst * delta_phi
        D_tot  = D_tot  + r0      * xjac * BigR * wst * delta_phi
        VP_tot = VP_tot + r0 * vpar0**2 * BB2       * xjac * BigR * wst * delta_phi
        TVP_tot= TVP_tot+ r0 * vpar0    * sqrt(BB2) * xjac * BigR * wst * delta_phi
        VK_tot = VK_tot + r0 * (dudx**2 + dudy**2) * BigR**2 * xjac * BigR * wst * delta_phi
        VM_tot = VM_tot + (dpsidx**2+dpsidy**2)/BigR**2 * xjac * BigR * wst * delta_phi
        J2_tot = J2_tot + eta_T_ohm * (ZJ0/BigR)**2 * xjac * BigR * wst * delta_phi

        P_max = max(P_max,r0 * T0)

        gradP_max     = max(gradP_max,grad_P)
        gradP_psi_max = max(gradP_psi_max,grad_P_psi)

        if (use_pellet) then

          call pellet_source2(pellet_amplitude,pellet_R,pellet_Z,pellet_psi,pellet_phi, &
                              pellet_radius, pellet_delta_psi, pellet_sig, pellet_length, pellet_ellipse, pellet_theta, &
                              x_g(ms,mt),y_g(ms,mt), ps0, phi, r0_corr, T0_corr/2.d0, central_density, &
                              pellet_particles, pellet_density, pellet_volume, source_pellet, source_volume)

          local_pellet_particles = local_pellet_particles + source_pellet * bigR * xjac * wst * delta_phi
          local_plasma_particles = local_plasma_particles + r0            * bigR * xjac * wst * delta_phi
          local_pellet_volume    = local_pellet_volume    + source_volume * bigR * xjac * wst * delta_phi

        endif

#if (defined WITH_Neutrals) && (!defined WITH_Impurities)
        !--- Calculate the neutral injection rate and the number of neutrals in the plasma

        source_neutral = 0.d0

        call total_neutral_source(x_g(ms,mt),y_g(ms,mt),phi,ps0,source_neutral)

        local_n_particles_inj = local_n_particles_inj + 0.5d0 * central_density * 1.d20 * source_neutral * bigR * xjac * wst * delta_phi / sqrt(MU_ZERO*central_mass*MASS_PROTON*central_density*1.d20)
        local_n_particles     = local_n_particles     + central_density * 1.d20 * rn0 * bigR * xjac * wst * delta_phi

#endif

        if (in_plasma(node_list,element_list,x_g(ms,mt),y_g(ms,mt),ps0,xpoint,xcase,ES%R_xpoint,ES%Z_xpoint,ES%psi_xpoint,psi_limit,ES%R_axis, ES%Z_axis,ES%psi_axis)) then

          D_int = D_int + r0        * xjac * BigR * wst * delta_phi
          P_int = P_int + r0 * T0   * xjac * BigR * wst * delta_phi
          C_intern = C_intern - zj0 /BigR * xjac *        wst * delta_phi    ! 2D integral
          Vol   = Vol   +             xjac * BigR * wst * delta_phi
          H_int = H_int + heat_source     * xjac * BigR * wst * delta_phi
          S_int = S_int + particle_source * xjac * BigR * wst * delta_phi
          VP_int = VP_int + r0 * vpar0**2 * BB2       * xjac * BigR * wst * delta_phi          
          TVP_int= TVP_int+ r0 * vpar0    * sqrt(BB2) * xjac * BigR * wst * delta_phi
          VK_int = VK_int + r0 * (dudx**2 + dudy**2) * BigR**2 * xjac * BigR * wst * delta_phi
          VM_int = VM_int + (dpsidx**2+dpsidy**2)/BigR**2 * xjac * BigR * wst * delta_phi
          J2_int = J2_int + eta_T_ohm * (ZJ0/BigR)**2 * xjac * BigR * wst * delta_phi

        else

          D_ext = D_ext + r0         * xjac * BigR * wst * delta_phi
          P_ext = P_ext + r0   * T0  * xjac * BigR * wst * delta_phi
          C_ext = C_ext - zj0 / BigR * xjac *        wst * delta_phi  ! 2D integral
          H_ext = H_ext + heat_source     * xjac * BigR * wst * delta_phi
          S_ext = S_ext + particle_source * xjac * BigR * wst * delta_phi
          VP_ext = VP_ext + r0 * vpar0**2 * BB2       * xjac * BigR * wst * delta_phi
          TVP_ext= TVP_ext+ r0 * vpar0    * sqrt(BB2) * xjac * BigR * wst * delta_phi
          VK_ext = VK_ext + r0 * (dudx**2 + dudy**2) * BigR**2 * xjac * BigR * wst * delta_phi
          VM_ext = VM_ext + (dpsidx**2+dpsidy**2)/BigR**2 * xjac * BigR * wst * delta_phi
          J2_ext = J2_ext + eta_T_ohm * (ZJ0/BigR)**2 * xjac * BigR * wst * delta_phi

        endif

      enddo
    enddo
  enddo

enddo
!$omp end do
!$omp end parallel

call MPI_AllReduce(D_int,density_in,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
call MPI_AllReduce(D_ext,density_out,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
call MPI_AllReduce(P_int,pressure_in,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
call MPI_AllReduce(P_ext,pressure_out,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
call MPI_AllReduce(C_intern,current_in,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
call MPI_AllReduce(C_ext,current_out,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
call MPI_AllReduce(Vol,Volume,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
call MPI_AllReduce(D_tot,density_tot,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
call MPI_AllReduce(P_tot,pressure,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
call MPI_AllReduce(H_ext,heating_out,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
call MPI_AllReduce(H_int,heating_in,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
call MPI_AllReduce(S_ext,source_out,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
call MPI_AllReduce(S_int,source_in,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

call MPI_AllReduce(VK_int,kin_perp_in,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
call MPI_AllReduce(VK_ext,kin_perp_out,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
call MPI_AllReduce(VK_tot,kin_perp_tot,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
call MPI_AllReduce(VM_int,mag_in,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
call MPI_AllReduce(VM_ext,mag_out,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
call MPI_AllReduce(VM_tot,mag_tot,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
call MPI_AllReduce(J2_int,ohm_in,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
call MPI_AllReduce(J2_ext,ohm_out,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
call MPI_AllReduce(J2_tot,ohm_tot,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

call MPI_AllReduce(VP_int,kin_par_in,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
call MPI_AllReduce(VP_ext,kin_par_out,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
call MPI_AllReduce(VP_tot,kin_par_tot,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
call MPI_AllReduce(TVP_int,mom_par_in,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
call MPI_AllReduce(TVP_ext,mom_par_out,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
call MPI_AllReduce(TVP_tot,mom_par_tot,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

if (use_pellet) then
  call MPI_AllReduce(local_pellet_particles,total_pellet_particles,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  call MPI_AllReduce(local_plasma_particles,total_plasma_particles,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  call MPI_AllReduce(local_pellet_volume,total_pellet_volume,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
endif

#if (defined WITH_Neutrals) && (!defined WITH_Impurities)
  call MPI_AllReduce(local_n_particles_inj, total_n_particles_inj,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  call MPI_AllReduce(local_n_particles, total_n_particles,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
#endif

rho_norm = central_density*1.d20 * central_mass * 1.67d-27
t_norm   = sqrt(MU_zero*rho_norm)

current_in  = n_period * current_in  / MU_zero / (2.d0 * PI)
current_out = n_period * current_out / MU_zero / (2.d0 * PI)
density_tot = n_period * density_tot * central_density
density_in  = n_period * density_in  * central_density
density_out = n_period * density_out * central_density
pressure    = n_period * pressure    / MU_zero * 1.5d0
pressure_in = n_period * pressure_in / MU_zero * 1.5d0
pressure_out= n_period * pressure_out/ MU_zero * 1.5d0
kin_par_tot = n_period * kin_par_tot / MU_zero * 0.5d0
kin_par_in  = n_period * kin_par_in  / MU_zero * 0.5d0
kin_par_out = n_period * kin_par_out / MU_zero * 0.5d0
mom_par_tot = n_period * mom_par_tot * rho_norm / sqrt(MU_zero * rho_norm)
mom_par_in  = n_period * mom_par_in  * rho_norm / sqrt(MU_zero * rho_norm)
mom_par_out = n_period * mom_par_out * rho_norm / sqrt(MU_zero * rho_norm)
kin_perp_tot= n_period * kin_perp_tot / MU_zero * 0.5d0
kin_perp_in = n_period * kin_perp_in  / MU_zero * 0.5d0
kin_perp_out= n_period * kin_perp_out / MU_zero * 0.5d0
mag_tot     = n_period * mag_tot      / MU_zero * 0.5d0
mag_in      = n_period * mag_in       / MU_zero * 0.5d0
mag_out     = n_period * mag_out     / MU_zero * 0.5d0
ohm_tot     = n_period * ohm_tot     / MU_zero / t_norm
ohm_in      = n_period * ohm_in      / MU_zero / t_norm
ohm_out     = n_period * ohm_out     / MU_zero / t_norm
heating_out = n_period * heating_out / MU_zero / t_norm * 1.5d0
heating_in  = n_period * heating_in  / MU_zero / t_norm * 1.5d0
source_out  = n_period * source_out  * central_density / t_norm
source_in   = n_period * source_in   * central_density / t_norm

if (my_id .eq. 0) then

  if (index_start .gt. 0) then
    xt = xtime(index_start)
  else
    xt = 0.d0
  endif

  write(*,'(A,3e14.6,A)') ' Time : ',xt,xt*t_norm,t_norm, ' [s]'
  write(*,'(A,4e14.6)')   ' Integrals_3D, PELLET : ',pellet_volume, total_pellet_volume, total_pellet_particles, total_plasma_particles
  write(*,'(A,2e14.6,A)') ' Volume   : ',xt,volume,' [m^3]'
  write(*,'(A,4e14.6,A)') ' density  (total/in/out) : ',xt,density_tot,  density_in,  density_out,' [10^20/m^3]'
  write(*,'(A,4e14.6,A)') ' pressure (total/in/out) : ',xt,pressure/1.d6, pressure_in/1.d6, pressure_out/1.d6,' [MJ]'
  write(*,'(A,4e14.6,A)') ' kinetic parallel (total/in/out) : ',xt,kin_par_tot/1.d6, kin_par_in/1.d6, kin_par_out/1.d6,' [MJ]'
  write(*,'(A,4e14.6,A)') ' parallel momentum (total/in/out) : ',xt,mom_par_tot, mom_par_in, mom_par_out,' [kg m/s]'


  write(*,'(A,4e14.6,A)') ' kinetic perp (total/in/out) : ',xt,kin_perp_tot/1.d6, kin_perp_in/1.d6, kin_perp_out/1.d6,' [MJ]'
  write(*,'(A,4e14.6,A)') ' magnetic (total/in/out) : ',xt,mag_tot/1.d6, mag_in/1.d6, mag_out/1.d6,' [MJ]'
  write(*,'(A,3e14.6,A)') ' current  (in/out)       : ',xt,current_in/1.d6, current_out/1.d6, ' [MA]'
  write(*,'(A,3e14.6,A)') ' heating  (in/out)       : ',xt,heating_in/1d6, heating_out/1.d6 ,' [MW]'
  write(*,'(A,3e14.6,A)') ' source   (in/out)       : ',xt,source_in, source_out,' [10^20/m^3/s]'
  write(*,'(A,4e14.6,A)') ' Ohmic    (in/out)       : ',xt,Ohm_tot/1.d6, Ohm_in/1.d6, Ohm_out/1.d6,' [MW]'

#ifndef fullmhd
  write(*,'(A,2e14.6)') ' li(3)    : ',xt, 2.d0 * mag_in/0.5  /(current_in**2 * R_geo * MU_zero)
  write(*,'(A,2e14.6)') ' betap(1) : ',xt, 4.d0 * pressure_in/1.5d0/(R_geo * current_in**2 * MU_zero)
#endif

  write(*,'(A120)') 'sum ,time ,density_tot, pressure, Wkin_par, Wkin_perp, Wmag, Ohm, heating, source'

  write(*,'(A,20e14.6)') 'sum ',xt,density_tot,pressure/1.d6,kin_par_tot/1.d6,kin_perp_tot/1.d6,mag_tot/1.d6, &
                                 Ohm_tot/1.d6,heating_in/1d6+heating_out/1.d6 ,source_in+source_out

#if (defined WITH_Neutrals) && (!defined WITH_Impurities)
  write(*,'(A,4e14.6)')   ' Integrals_3D, MGI : ', total_n_particles_inj, total_n_particles
#endif

endif

return

end

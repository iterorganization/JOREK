module mod_integrate_recomb

implicit none
contains

#include "corr_neg_include.f90"
subroutine integrate_recombination(my_id ,n_mpi, rec_rate_local, rec_v_R, rec_v_Z, rec_v_phi,volume_check, energy_neutrals, energy_radiation)
!---------------------------------------------------------------
!
!---------------------------------------------------------------
use data_structure
use gauss
use basis_at_gaussian
use phys_module
use nodes_elements
use mpi
use constants, only : TWOPI
! use mod_ionisation_recombination, only : rec_rate_global
! use mod_ionisation_recombination, only : rec_rate_local, rec_rate_global, rec_mom_local,rec_energy_local, rec_v_R, rec_v_Z, rec_v_phi
use mod_atomic_coeff_deuterium, only: rec_rate_to_kinetic 
!$ use omp_lib

implicit none

type (type_element)      :: element
type (type_node)         :: nodes(n_vertex_max)

integer    :: n_mpi, my_id
real*8, dimension(:,:), allocatable, intent(out) :: rec_rate_local , rec_v_R, rec_v_Z, rec_v_phi 
real*8, dimension(:,:), allocatable, optional, intent(out) :: volume_check, energy_neutrals, energy_radiation    

real*8     :: x_g(n_gauss,n_gauss),        x_s(n_gauss,n_gauss),        x_t(n_gauss,n_gauss)
real*8     :: y_g(n_gauss,n_gauss),        y_s(n_gauss,n_gauss),        y_t(n_gauss,n_gauss)
real*8, dimension(n_plane,n_var,n_gauss,n_gauss) :: eq_g, eq_s, eq_t

! real*8     :: eq_g(n_var,n_gauss,n_gauss), eq_s(n_var,n_gauss,n_gauss), eq_t(n_var,n_gauss,n_gauss)
real*8     :: density_eq(n_gauss,n_gauss), eq_g1(n_var,n_gauss,n_gauss), Fprofile(n_gauss,n_gauss)

integer    :: i, j, k, in, ms, mt,mp, iv, inode, ife,ielm, n_elements !,n_mpi
real*8     :: xjac, BigR, wst,delta_phi
real*8     :: ps0_x, ps0_y, u0_x, u0_y, r0, T0,vpar0
real*8     :: T0_corr, r0_corr

! rec_rate parameters
!real*8, intent(in)    :: Te0, ne0                        ! Electron temperature in JOREK units
real*8  :: Sion_T , dSion_dT           ! Normalized ionization coefficient and its temperature derivative
real*8  :: Srec_T , dSrec_dT           ! Normalized recombination coefficient and its temperature derivative
real*8  :: LradDcont_T, dLradDcont_dT 
real*8  :: LradDcont_corr, dLradDcont_dT_corr

!real*8     :: Sum_rec(n_gauss,n_gauss)
integer    :: missing, loc_rec_elms
integer, dimension(:), allocatable :: local_rec_elements
!define mpi local elements

! local_rec_elements tells every MPI process how many elements he gets.
allocate(local_rec_elements(n_mpi))
n_elements = element_list%n_elements
missing = mod(n_elements,n_mpi) !< tells us if n_elements can be evenly shared over all MPI processes

!> give all MPI processes the same amount of elements.
!> distribute the "missing" elements over the first MPI processes until we have all elements covered.
do i=1,n_mpi 
  local_rec_elements(i) = floor(n_elements/real(n_mpi,8))
  if (missing .gt. 0) then
      missing = missing - 1
      local_rec_elements(i) = local_rec_elements(i) +1  
  endif
enddo !n_mpi
!write(*,*) "length local rec list", local_rec_elements(my_id+1)


!if not allocated, allocate rec_variables of size (n_elements)
if(.not. allocated(rec_rate_local)) then
  allocate(rec_rate_local(local_rec_elements(my_id+1), n_plane  )) !< local_rec_elements bcause it's local
  allocate(rec_v_R(local_rec_elements(my_id+1), n_plane    ))
  allocate(rec_v_Z(local_rec_elements(my_id+1), n_plane    ))
  allocate(rec_v_phi(local_rec_elements(my_id+1), n_plane  ))

  allocate(volume_check(local_rec_elements(my_id+1), n_plane))
  allocate(energy_neutrals(local_rec_elements(my_id+1), n_plane))
  allocate(energy_radiation(local_rec_elements(my_id+1), n_plane))
endif

!> can now be done local?
rec_rate_local(:,:)   = 0.d0
!> momentum density in R,Z, phi direction. rho_rec * v
rec_v_R(:,:)          = 0.d0
rec_v_Z(:,:)          = 0.d0
rec_v_phi(:,:)        = 0.d0
!> volume check. 
volume_check(:,:)     = 0.d0
energy_neutrals(:,:)  = 0.d0
energy_radiation(:,:) = 0.d0

delta_phi     = TWOPI / real(n_plane,8) / real(n_period,8)
!HZ_p,n_plane,n_gauss,n_order,n_vertex_max,TWOPI
!$omp parallel do default(none)                                           &
!$omp schedule(static, 100)                                               &
!$omp   shared(local_rec_elements,my_id,n_mpi, volume_check,              &
!$omp          energy_neutrals, energy_radiation, rec_rate_local,         &
!$omp          rec_v_R,rec_v_Z,rec_v_phi,element_list,node_list, H, H_s,  & 
!$omp          H_t, HZ, tstep,F0, delta_phi, gamma                        &
!$omp          )                                                          &
!$omp   private(ife,ielm,iv,i,j,k,ms,mt,mp,in,                            &
!$omp           inode,element,                                            &
!$omp           x_g, y_g, x_s, y_s, x_t, y_t, xjac, eq_g, eq_s, eq_t,     &
!$omp           wst, BigR, r0, T0,  ps0_x,ps0_y ,u0_x,u0_y,vpar0,         &
!$omp           r0_corr, T0_corr, Sion_T, dSion_dT, Srec_T, dSrec_dT,     &
!$omp           LradDcont_T, dLradDcont_dT, LradDcont_corr,               &
!$omp           dLradDcont_dT_corr)                                       &
!$omp   firstprivate(nodes) 
! The firstprivate is used to explicitly force making a copy of nodes where the values member is not allocated, 
! which is then allocated in the make_deep_copy subroutine. This is done because there is no guarantee when OMP 
! makes a private copy for each thread, whether a deep or shallow copy is made (if the values member of nodes 
! point to the same memory location, or if new memory is allocated for the nodes on each thread), and would be 
! problematic is the nodes from each thread modifies the same memory location.
!> loop over all local recombination elements
do ife =1,  local_rec_elements(my_id+1) !element_list%n_elements !n_local_rec_elms

  !> real element number
  !  Note: By intention, we split the element order over the MPI processes 
  ! e.g. if n_mpi =4 and n_elements = 18 ,then my_id = 0 gets ielm 1,5,9,13,17. my_id = 1 gets ielm 2,6,10,14,18
  ! my_id = 3 gets ielm 3,7,11,15. my_id = 4 gets ielm 4,8,12,16.
  ! because we give the first MPI processes the most work.
  ! The recombination rate is localized in the plasma. So neighboring elements have similar rates.
  ! So to share the load (for the high recombination areas) over the MPI processes, 
  ! we distribute elements close to each other over different MPI processes.
  !if you want to change this, you should probably also change local_rec_elements, and the use of ielm in do_1particle_recombination
  ielm    = (my_id+1) + n_mpi*(ife - 1) !< this splits nearby elements over different nodes
  element = element_list%element(ielm)

  do iv = 1, n_vertex_max
    inode     = element%vertex(iv)
    call make_deep_copy_node(node_list%node(inode), nodes(iv))
  enddo
  
!----------------------- from elt_matrix_fft  
!---------------------------------------------------- value of (x,y) and derivatives on Gaussian points
x_g  = 0.d0; x_s  = 0.d0; x_t  = 0.d0!; x_st  = 0.d0; x_ss  = 0.d0; x_tt  = 0.d0;
y_g  = 0.d0; y_s  = 0.d0; y_t  = 0.d0!; y_st  = 0.d0; y_ss  = 0.d0; y_tt  = 0.d0;
eq_g = 0.d0; eq_s = 0.d0; eq_t = 0.d0 !; eq_st = 0.d0; eq_ss = 0.d0; eq_tt = 0.d0; eq_p = 0.d0;
do i=1,n_vertex_max
  do j=1,n_order+1
    do ms=1, n_gauss
      do mt=1, n_gauss

        x_g(ms,mt)  = x_g(ms,mt)  + nodes(i)%x(1,j,1) * element%size(i,j) * H(i,j,ms,mt)
        x_s(ms,mt)  = x_s(ms,mt)  + nodes(i)%x(1,j,1) * element%size(i,j) * H_s(i,j,ms,mt)
        x_t(ms,mt)  = x_t(ms,mt)  + nodes(i)%x(1,j,1) * element%size(i,j) * H_t(i,j,ms,mt)

        y_g(ms,mt)  = y_g(ms,mt)  + nodes(i)%x(1,j,2) * element%size(i,j) * H(i,j,ms,mt)
        y_s(ms,mt)  = y_s(ms,mt)  + nodes(i)%x(1,j,2) * element%size(i,j) * H_s(i,j,ms,mt)
        y_t(ms,mt)  = y_t(ms,mt)  + nodes(i)%x(1,j,2) * element%size(i,j) * H_t(i,j,ms,mt)

      end do !ms=1, n_gauss
    end do ! mt=1, n_gauss

    do ms=1, n_gauss
      do mt=1, n_gauss
        do k=1,n_var ![1,2,5,6,7]
          do in=1,n_tor
            do mp=1,n_plane
              eq_g(mp,k,ms,mt) = eq_g(mp,k,ms,mt) + nodes(i)%values(in,j,k) * element%size(i,j) * H(i,j,ms,mt)  * HZ(in,mp)
              eq_s(mp,k,ms,mt) = eq_s(mp,k,ms,mt) + nodes(i)%values(in,j,k) * element%size(i,j) * H_s(i,j,ms,mt)* HZ(in,mp)
              eq_t(mp,k,ms,mt) = eq_t(mp,k,ms,mt) + nodes(i)%values(in,j,k) * element%size(i,j) * H_t(i,j,ms,mt)* HZ(in,mp)

            enddo

          enddo

        enddo

      enddo
    enddo
  enddo
enddo

!------------------------end from elt_matrix_fft 

!--------------------------------------------------- sum over the Gaussian integration points
  do mp=1, n_plane

    do ms=1, n_gauss

      do mt=1, n_gauss

        wst = wgauss(ms)*wgauss(mt)

        xjac = x_s(ms,mt)*y_t(ms,mt) - x_t(ms,mt)*y_s(ms,mt)
        BigR = x_g(ms,mt)

        ps0_x  = (   y_t(ms,mt) * eq_s(mp,1,ms,mt) - y_s(ms,mt) * eq_t(mp,1,ms,mt) ) / xjac
        ps0_y  = ( - x_t(ms,mt) * eq_s(mp,1,ms,mt) + x_s(ms,mt) * eq_t(mp,1,ms,mt) ) / xjac
        u0_x  = (   y_t(ms,mt) * eq_s(mp,2,ms,mt) - y_s(ms,mt) * eq_t(mp,2,ms,mt) ) / xjac
        u0_y  = ( - x_t(ms,mt) * eq_s(mp,2,ms,mt) + x_s(ms,mt) * eq_t(mp,2,ms,mt) ) / xjac

        r0    = eq_g(mp,5,ms,mt)
        r0_corr = corr_neg_dens1(r0)
		
#ifdef WITH_TiTe
        T0       = eq_g(mp,var_Te,ms,mt)
        T0_corr  = corr_neg_temp1(T0*2.d0)/2.d0
        vpar0    = eq_g(mp,var_vpar,ms,mt)
        call rec_rate_to_kinetic(r0, T0, Sion_T, dSion_dT, Srec_T, dSrec_dT, LradDcont_T, dLradDcont_dT, LradDcont_corr, dLradDcont_dT_corr)
#else
        T0       = eq_g(mp,var_T,ms,mt)
        T0_corr  = corr_neg_temp1(T0)
        vpar0    = eq_g(mp,var_vpar,ms,mt)
        call rec_rate_to_kinetic(r0, 0.5d0*T0, Sion_T, dSion_dT, Srec_T, dSrec_dT, LradDcont_T, dLradDcont_dT, LradDcont_corr, dLradDcont_dT_corr)

        ! --- Transform derivatives on Te to derivatives in total T	
        dSrec_dT           = dSrec_dT      / 2.d0	
        dLradDcont_dT      = dLradDcont_dT / 2.d0
        dLradDcont_dT_corr = dLradDcont_dT_corr / 2.d0
#endif
    
        !> neutral density gain in element due to recombination
        rec_rate_local(ife,mp)   = rec_rate_local(ife,mp)+ (Srec_T * r0_corr * r0_corr)                                      *BigR *xjac *tstep * delta_phi *wst ! rho_rec
        !> momentum density in R,Z, phi direction. rho_rec * v
        rec_v_R(ife,mp)          = rec_v_R(ife,mp)       + (Srec_T * r0_corr * r0_corr) * (-BigR*u0_y  + vpar0/BigR * ps0_y) *BigR *xjac *tstep * delta_phi *wst !rho_rec*v_R
        rec_v_Z(ife,mp)          = rec_v_Z(ife,mp)       + (Srec_T * r0_corr * r0_corr) * (+ BigR*u0_x - vpar0/BigR * ps0_x) *BigR *xjac *tstep * delta_phi *wst !rho_rec*v_Z
        rec_v_phi(ife,mp)        = rec_v_phi(ife,mp)     + (Srec_T * r0_corr * r0_corr) * (+ F0*vpar0/BigR)                  *BigR *xjac *tstep * delta_phi *wst !rho_rec*v_phi
        !> volume check. 
        volume_check(ife,mp)     = volume_check(ife,mp)  + (1.d0)                                                            *BigR *xjac        * delta_phi *wst
        energy_neutrals(ife,mp)  = energy_neutrals(ife,mp)+ (gamma-1.d0) * 0.5d0 *T0 * r0_corr * r0_corr  * Srec_T           *BigR *xjac *tstep * delta_phi *wst 
        energy_radiation(ife,mp) = energy_radiation(ife,mp)+ r0_corr * r0_corr  * LradDcont_corr                             *BigR *xjac *tstep * delta_phi *wst       

      enddo !mt
    enddo !ms
  enddo !mp !in, n_tor
  do iv = 1, n_vertex_max
    call dealloc_node(nodes(iv))
  enddo

enddo !ife
!$omp end parallel do

! --- Volume integrals above were over wedge from 0 to (2 * PI) / n_period
! --- We need to account for all n_period wedges
rec_rate_local   = n_period * rec_rate_local
rec_v_R          = n_period * rec_v_R
rec_v_Z          = n_period * rec_v_Z
rec_v_phi        = n_period * rec_v_phi
volume_check     = n_period * volume_check
energy_neutrals  = n_period * energy_neutrals
energy_radiation = n_period * energy_radiation

!return
end subroutine

end module mod_integrate_recomb

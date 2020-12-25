!> Particle pusher module with the modified Qin variational scheme.
!> This module contains routines for pushing particles in the RZPhi
!>  (cylindrical) coordinate system. 
module mod_gc_variational
  use mod_particle_types
  use constants, only: EL_CHG, ATOMIC_MASS_UNIT
  implicit none  
  public copy_particle_gc_vpar, copy_particle_gc_Qin, copy_particle_gc_Qin_to_vpar
  public initialise_gc_Qin, push_gc_Qin, push_gc_rk4
    
contains

subroutine copy_particle_gc_vpar(particle_in, particle_out)
type(particle_gc_vpar), intent(in)  :: particle_in
type(particle_gc_vpar)              :: particle_out
    
particle_out%i_elm = particle_in%i_elm
particle_out%mu    = particle_in%mu
particle_out%x     = particle_in%x
particle_out%vpar  = particle_in%vpar
particle_out%st    = particle_in%st
particle_out%q     = particle_in%q

return
end

subroutine copy_particle_gc_Qin(particle_in,particle_out)
class(particle_gc_Qin), intent(in)  :: particle_in
class(particle_gc_Qin), intent(out) :: particle_out
  
particle_out%i_elm = particle_in%i_elm
particle_out%mu    = particle_in%mu
particle_out%x     = particle_in%x
particle_out%vpar  = particle_in%vpar
particle_out%st    = particle_in%st
particle_out%q     = particle_in%q

particle_out%x_m      = particle_in%x_m
particle_out%vpar_m   = particle_in%vpar_m
particle_out%Astar_m  = particle_in%Astar_m 
particle_out%Astar_k  = particle_in%Astar_k 
particle_out%dAstar_k = particle_in%dAstar_k
particle_out%Bn_k     = particle_in%Bn_k
particle_out%dBn_k    = particle_in%dBn_k
particle_out%Bnorm_k  = particle_in%Bnorm_k

return
end

subroutine copy_particle_gc_Qin_to_vpar(particle_in,particle_out)
use mod_particle_types
class(particle_gc_Qin), intent(in)   :: particle_in
class(particle_gc_vpar), intent(out) :: particle_out
  
particle_out%i_elm = particle_in%i_elm
particle_out%mu    = particle_in%mu
particle_out%x     = particle_in%x
particle_out%vpar  = particle_in%vpar
particle_out%st    = particle_in%st
particle_out%q     = particle_in%q
  
return
end

subroutine copy_particle_gc_Vpar_to_Qin(particle_in,particle_out)
use mod_particle_types
class(particle_gc_Vpar), intent(in) :: particle_in
class(particle_gc_Qin), intent(out) :: particle_out
    
particle_out%i_elm = particle_in%i_elm
particle_out%mu    = particle_in%mu
particle_out%x     = particle_in%x
particle_out%vpar  = particle_in%vpar
particle_out%st    = particle_in%st
particle_out%q     = particle_in%q

return
end
   
subroutine initialise_gc_Qin(fields, particle_Qin, mass, timestep)
use mod_particle_types
use mod_fields, only: fields_base
class(fields_base)      :: fields
type(particle_gc_Qin)  :: particle_Qin
real*8, intent(in)     :: timestep ! [s]
real*8, intent(in)     :: mass     ! [amu]
type(particle_gc_Vpar) :: particle_Vpar

real*8 :: A_m(3), dA_m(3,3), B_m(3), dB_m(3,3), E_m(3), bn_m, dbn_m(3), bnorm_m(3), dbnorm_m(3,3)
real*8 :: A_k(3), dA_k(3,3), B_k(3), dB_k(3,3), E_k(3), bn_k, dbn_k(3), bnorm_k(3), dbnorm_k(3,3)
real*8 :: qom, time_0

time_0 = 0.d0
qom = particle_Qin%q * EL_CHG / (mass * ATOMIC_MASS_UNIT) 

call fields%calc_Qin(time_0, particle_Qin%i_elm, particle_Qin%st, particle_Qin%x(3), &
                     A_k, dA_k, B_k, dB_k, Bnorm_k, dBnorm_k, bn_k, dbn_k, E_k)
!call fields%calc_Qin_analytic(particle_Qin%x(1),particle_Qin%x(2), particle_Qin%x(3), A_k, dA_k, B_k, dB_k, Bnorm_k, dBnorm_k, bn_k, dbn_k, E_k)

particle_Qin%Astar_k  = qom *  A_k + particle_Qin%vpar *  bnorm_k
particle_Qin%dAstar_k = qom * dA_k + particle_Qin%vpar * dbnorm_k 
particle_Qin%Bn_k     = Bn_k   
particle_Qin%dBn_k    = dBn_k   
particle_Qin%Bnorm_k  = Bnorm_k 

call copy_particle_gc_Qin_to_Vpar(particle_Qin, particle_Vpar)

call push_gc_rk4(fields, particle_Vpar, mass, -timestep, 1)      ! should be only the very first call

particle_Qin%x_m    = particle_Vpar%x
particle_Qin%vpar_m = particle_Vpar%vpar

call fields%calc_Qin(time_0, particle_Vpar%i_elm, particle_Vpar%st, particle_Vpar%x(3), &
                     A_m, dA_m, B_m, dB_m, Bnorm_m, dBnorm_m, bn_m, dbn_m, E_m)
!call fields%calc_Qin_analytic(particle_Vpar%x(1), particle_Vpar%x(2), particle_Vpar%x(3), A_m, dA_m, B_m, dB_m, Bnorm_m, dBnorm_m, bn_m, dbn_m, E_m)

particle_Qin%Astar_m = qom * A_m + particle_Qin%vpar_m *  Bnorm_m
 
return
end


subroutine initialise_gc_Qin2(fields, particle_Qin, mass, timestep)
use mod_particle_types
use mod_fields, only: fields_base
class(fields_base)      :: fields
type(particle_gc_Qin)  :: particle_Qin
real*8, intent(in)     :: timestep ! [s]
real*8, intent(in)     :: mass     ! [amu]
type(particle_gc_Vpar) :: particle_Vpar

real*8 :: A_m(3), dA_m(3,3), B_m(3), dB_m(3,3), E_m(3), bn_m, dbn_m(3), bnorm_m(3), dbnorm_m(3,3)
real*8 :: A_k(3), dA_k(3,3), B_k(3), dB_k(3,3), E_k(3), bn_k, dbn_k(3), bnorm_k(3), dbnorm_k(3,3)
real*8 :: qom, time_0

time_0 = 0.d0
qom = particle_Qin%q * EL_CHG / (mass * ATOMIC_MASS_UNIT) 

call fields%calc_Qin(time_0, particle_Qin%i_elm, particle_Qin%st, particle_Qin%x(3), &
                     A_m, dA_m, B_m, dB_m, Bnorm_m, dBnorm_m, bn_m, dbn_m, E_m)

particle_Qin%Astar_m  =  qom * A_m + particle_Qin%vpar *  Bnorm_m

call copy_particle_gc_Qin_to_Vpar(particle_Qin, particle_Vpar)

call push_gc_rk4(fields, particle_Vpar, mass, timestep, 1)      ! should be only the very first call

particle_Qin%x_m    = particle_Qin%x
particle_Qin%vpar_m = particle_Qin%vpar
particle_Qin%x      = particle_Vpar%x
particle_Qin%vpar   = particle_Vpar%vpar
particle_Qin%i_elm  = particle_Vpar%i_elm
particle_Qin%st     = particle_Vpar%st

call fields%calc_Qin(time_0, particle_Qin%i_elm, particle_Qin%st, particle_Qin%x(3), &
                     A_k, dA_k, B_k, dB_k, Bnorm_k, dBnorm_k, bn_k, dbn_k, E_k)

particle_Qin%Astar_k  = qom *  A_k + particle_Qin%vpar *  bnorm_k
particle_Qin%dAstar_k = qom * dA_k + particle_Qin%vpar * dbnorm_k 
particle_Qin%Bn_k     = Bn_k   
particle_Qin%dBn_k    = dBn_k   
particle_Qin%Bnorm_k  = Bnorm_k 

return
end

subroutine push_gc_Qin(fields, particle_Qin, mass, timestep, n_steps)
use nodes_elements
use mod_find_rz_nearby
use mod_fields, only: fields_base
implicit none
class(fields_base)      :: fields
type(particle_gc_Qin)  :: particle_Qin
real*8, intent(in)     :: timestep ! [s]
real*8, intent(in)     :: mass     ! [amu]
integer, intent(in)    :: n_steps
real*8                 :: time_0   ! not yet implemented
real*8                 :: qom

real*8 :: x_k(3), vpar_k, Astar_k(3), st_k(2), x_p(3), st_p(2)
real*8 :: A_p(3), dA_p(3,3), B_p(3), dB_p(3,3), E_p(3), bn_p, dbn_p(3), bnorm_p(3), dbnorm_p(3,3)
real*8 :: Astar_p(3), dAstar_p(3,3)
real*8 :: newton_rhs(3), newton_matrix(3,3), residue, mp_k, wp_k, mp_p, wp_p, mp_m, wp_m
integer :: info, ifail, ipiv(3), iter, i, j, it, i_elm_k, i_elm_p
  
qom = particle_Qin%q * EL_CHG / (mass * ATOMIC_MASS_UNIT) 

do it=1, n_steps

  x_k     = particle_Qin%x
  vpar_k  = particle_Qin%vpar
  Astar_k = particle_Qin%Astar_k
  st_k    = particle_Qin%st
  i_elm_k = particle_Qin%i_elm
  
  !----------------------------- initial guess from linearised equations
  do j=1,3
    newton_rhs(j) = - 2.d0 * particle_Qin%bnorm_k(j) * (particle_Qin%vpar_m + vpar_k) &
                    + 2.d0 * timestep * particle_Qin%mu * particle_Qin%dBn_k(j)
    do i=1,3
      newton_matrix(i,j) = particle_Qin%dAstar_k(i,j) - particle_Qin%dAstar_k(j,i) &
                         - 2.d0 * particle_Qin%bnorm_k(j) *  particle_Qin%bnorm_k(i) / timestep
    enddo
  enddo
  call dgesv(3, 1, newton_matrix, 3, ipiv, newton_rhs, 3, info)

  particle_Qin%x    =   particle_Qin%x_m + newton_rhs
  particle_Qin%vpar = - particle_Qin%vpar_m - 2.d0 * particle_Qin%vpar &
                    + 2.d0 * dot_product(particle_Qin%bnorm_k, particle_Qin%x - particle_Qin%x_m) / timestep

  call find_RZ_nearby(node_list, element_list, x_k(1), x_k(2), st_k(1), st_k(2), i_elm_k, &
                      particle_Qin%x(1), particle_Qin%x(2), particle_Qin%st(1), particle_Qin%st(2), particle_Qin%i_elm, ifail)

  call fields%calc_Qin(time_0, particle_Qin%i_elm, particle_Qin%st, particle_Qin%x(3), &
                       A_p, dA_p, B_p, dB_p, Bnorm_p, dBnorm_p, bn_p, dbn_p, E_p)  
  !call fields%calc_Qin_analytic(particle_Qin%x(1), particle_Qin%x(2), particle_Qin%x(3), A_p, dA_p, B_p, dB_p, Bnorm_p, dBnorm_p, bn_p, dbn_p, E_p)  

  Astar_p  = qom *  A_p + particle_Qin%vpar *  bnorm_p
  dAstar_p = qom * dA_p + particle_Qin%vpar * dbnorm_p

  !---------------------------- Newton iterations
  do iter =1,5
 
     do i=1,3

      newton_rhs(i) = (Astar_p(i) - particle_Qin%Astar_m(i)) + 2.d0 * timestep * particle_Qin%mu * particle_Qin%dBn_k(i) &

                    - dot_product(particle_Qin%dAstar_k(:,i), particle_Qin%x - particle_Qin%x_m)
!     if (i .eq. 3) write(*,*) it, iter, norm2(newton_rhs)
     do j=1,3
        
        newton_matrix(j,i) = + particle_Qin%dAstar_k(i,j) - qom * dA_p(j,i)           &
                             
                             + dbnorm_p(j,i) * (particle_Qin%vpar_m + 2.d0 * vpar_k)  &
        
                             - 2.d0 * bnorm_p(j) * particle_Qin%bnorm_k(i) / timestep &

                             - 2.d0 * dbnorm_p(j,i) * dot_product(particle_Qin%bnorm_k, particle_Qin%x - particle_Qin%x_m) / timestep
      enddo
    enddo

    call dgesv(3, 1, newton_matrix, 3, ipiv, newton_rhs, 3, info)

    x_p     = particle_Qin%x
    st_p    = particle_Qin%st
    i_elm_p = particle_Qin%i_elm
    
    particle_Qin%x    = particle_Qin%x + newton_rhs
    particle_Qin%vpar = - particle_Qin%vpar_m - 2.d0 * vpar_k + 2.d0 * dot_product(particle_Qin%bnorm_k, particle_Qin%x - particle_Qin%x_m) / timestep
  
    call find_RZ_nearby(node_list, element_list, x_p(1), x_p(2), st_p(1), st_p(2), i_elm_p, &
                        particle_Qin%x(1), particle_Qin%x(2), particle_Qin%st(1), particle_Qin%st(2), particle_Qin%i_elm, ifail)

    call fields%calc_Qin(time_0, particle_Qin%i_elm, particle_Qin%st, particle_Qin%x(3), &
                         A_p, dA_p, B_p, dB_p, Bnorm_p, dBnorm_p, bn_p, dbn_p, E_p)  
    !call fields%calc_Qin_analytic(particle_Qin%x(1), particle_Qin%x(2), particle_Qin%x(3), A_p, dA_p, B_p, dB_p, Bnorm_p, dBnorm_p, bn_p, dbn_p, E_p)  
  
    Astar_p  = qom *  A_p + particle_Qin%vpar *  bnorm_p
    dAstar_p = qom * dA_p + particle_Qin%vpar * dbnorm_p

  enddo ! newton iterations


  do j=1,3    
    newton_rhs(j) = (Astar_p(j) - particle_Qin%Astar_m(j)) + 2.d0 * timestep * particle_Qin%mu * particle_Qin%dBn_k(j) &
                 - dot_product(particle_Qin%dAstar_k(:,j), particle_Qin%x - particle_Qin%x_m)
  enddo
  residue = norm2(newton_rhs)
!  if (residue .gt. 1.d-6) write(*,'(A,i4,e14.6)') ' Qin : residue : ',it, residue

!  mp_p = Astar_p(3)
!  wp_p = 0.25d0 * (vpar_k**2 + particle_Qin%vpar**2) + 0.5d0 * particle_Qin%mu * (particle_Qin%bn_k + bn_p)

  particle_Qin%x_m      = x_k
  particle_Qin%vpar_m   = vpar_k
  particle_Qin%Astar_m  = Astar_k 
  particle_Qin%Astar_k  = Astar_p 
  particle_Qin%dAstar_k = dAstar_p 
  particle_Qin%Bn_k     = Bn_p
  particle_Qin%dBn_k    = dBn_p
  particle_Qin%Bnorm_k  = Bnorm_p

enddo

return
end

subroutine push_gc_rk4(fields, particle_gc, mass, timestep, n_steps)
use nodes_elements
use mod_find_rz_nearby
use mod_fields, only: fields_base
implicit none
class(fields_base)     :: fields
type(particle_gc_vpar) :: particle_gc
real*8, intent(in)     :: timestep ! [s]
real*8, intent(in)     :: mass     ! [amu]
integer, intent(in)    :: n_steps
real*8 :: A_0(3), dA_0(3,3), B_0(3), dB_0(3,3), bn_0, dbn_0(3), Bnorm_0(3), dBnorm_0(3,3), E_0(3)
real*8 :: A_1(3), dA_1(3,3), B_1(3), dB_1(3,3), bn_1, dbn_1(3), Bnorm_1(3), dBnorm_1(3,3), E_1(3)
real*8 :: A_2(3), dA_2(3,3), B_2(3), dB_2(3,3), bn_2, dbn_2(3), Bnorm_2(3), dBnorm_2(3,3), E_2(3)
real*8 :: A_3(3), dA_3(3,3), B_3(3), dB_3(3,3), bn_3, dbn_3(3), Bnorm_3(3), dBnorm_3(3,3), E_3(3)
real*8 :: x_0(3), x_1(3), x_2(3), x_3(3), u_0, u_1, u_2, u_3, time_0, time_1, time_2, time_3
real*8 :: delta_x1(3), delta_x2(3), delta_x3(3), delta_x4(3)
real*8 :: delta_u1, delta_u2, delta_u3, delta_u4, qom, p_phi, energy
type(particle_gc_vpar) :: p_0, p_1, p_2, p_3
integer :: i, ifail

qom = particle_gc%q * EL_CHG / (mass * ATOMIC_MASS_UNIT) 

call copy_particle_gc_vpar(particle_gc,p_0)
call copy_particle_gc_vpar(particle_gc,p_1)
call copy_particle_gc_vpar(particle_gc,p_2)
call copy_particle_gc_vpar(particle_gc,p_3)

call fields%calc_RK4(time_0, p_0%i_elm, p_0%st, p_0%x(3), A_0, dA_0, B_0, dB_0, Bnorm_0, dBnorm_0, bn_0, dbn_0, E_0)
!call fields%calc_RK4_analytic(p_0%x(1), p_0%x(2), p_0%x(3), A_0, dA_0, B_0, dB_0, Bnorm_0, dBnorm_0, bn_0, dbn_0, E_0)

do i =1, n_steps
  
  call rk4_step(p_0%x, p_0%vpar, qom, p_0%mu, E_0, B_0, Bnorm_0, dBnorm_0, dBn_0, delta_x1, delta_u1)

  p_1%x    = p_0%x    + 0.5d0 * timestep * delta_x1
  p_1%vpar = p_0%vpar + 0.5d0 * timestep * delta_u1

  call find_RZ_nearby(node_list, element_list, p_0%x(1), p_0%x(2), p_0%st(1), p_0%st(2), p_0%i_elm, &
                                               p_1%x(1), p_1%x(2), p_1%st(1), p_1%st(2), p_1%i_elm, ifail)
    
  call fields%calc_RK4(time_1, p_1%i_elm, p_1%st, p_1%x(3), A_1, dA_1, B_1, dB_1, Bnorm_1, dBnorm_1, bn_1, dbn_1, E_1)
!  call fields%calc_RK4_analytic(p_1%x(1), p_1%x(2), p_1%x(3), A_1, dA_1, B_1, dB_1, Bnorm_1, dBnorm_1, bn_1, dbn_1, E_1)

  call rk4_step(p_1%x, p_1%vpar, qom, p_1%mu, E_1, B_1, Bnorm_1, dBnorm_1, dBn_1, delta_x2, delta_u2)

  p_2%x    = p_0%x    + 0.5d0 * timestep * delta_x2
  p_2%vpar = p_0%vpar + 0.5d0 * timestep * delta_u2

  call find_RZ_nearby(node_list, element_list, p_0%x(1), p_0%x(2), p_0%st(1), p_0%st(2), p_0%i_elm, &
                                               p_2%x(1), p_2%x(2), p_2%st(1), p_2%st(2), p_2%i_elm, ifail)
  
  call fields%calc_RK4(time_2, p_2%i_elm, p_2%st, p_2%x(3), A_2, dA_2, B_2, dB_2, Bnorm_2, dBnorm_2, bn_2, dbn_2, E_2)
!  call fields%calc_RK4_analytic(p_2%x(1), p_2%x(2), p_2%x(3), A_2, dA_2, B_2, dB_2, Bnorm_2, dBnorm_2, bn_2, dbn_2, E_2)
                           
  call rk4_step(p_2%x, p_2%vpar, qom, p_2%mu, E_2, B_2, Bnorm_2, dBnorm_2, dBn_2, delta_x3, delta_u3)

  p_3%x    = p_0%x    + timestep * delta_x3
  p_3%vpar = p_0%vpar + timestep * delta_u3

  call find_RZ_nearby(node_list, element_list, p_0%x(1), p_0%x(2), p_0%st(1), p_0%st(2), p_0%i_elm, &
                                               p_3%x(1), p_3%x(2), p_3%st(1), p_3%st(2), p_3%i_elm, ifail)
  
  call fields%calc_RK4(time_3, p_3%i_elm, p_3%st, p_3%x(3), A_3, dA_3, B_3, dB_3, Bnorm_3, dBnorm_3, bn_3, dbn_3, E_3)
!  call fields%calc_RK4_analytic(p_3%x(1), p_3%x(2), p_3%x(3), A_3, dA_3, B_3, dB_3, Bnorm_3, dBnorm_3, bn_3, dbn_3, E_3)
    
  call rk4_step(p_3%x, p_3%vpar, qom, p_3%mu, E_3, B_3, Bnorm_3, dBnorm_3, dBn_3, delta_x4, delta_u4)
                                
  p_0%x    = p_0%x    + timestep * (delta_x1 + 2.d0*delta_x2 + 2.d0*delta_x3 + delta_x4)/6.d0
  p_0%vpar = p_0%vpar + timestep * (delta_u1 + 2.d0*delta_u2 + 2.d0*delta_u3 + delta_u4)/6.d0

  call find_RZ_nearby(node_list, element_list, &
                      particle_gc%x(1), particle_gc%x(2), particle_gc%st(1), particle_gc%st(2), particle_gc%i_elm, &
                      p_0%x(1),  p_0%x(2),  p_0%st(1),  p_0%st(2),  p_0%i_elm, ifail)

enddo

call copy_particle_gc_vpar(p_0,particle_gc)

return
end

function cross(A,B) result(AcrossB)
implicit none
real*8 :: A(3), B(3), AcrossB(3)
    AcrossB(1) = A(2)*B(3) - A(3)*B(2)
    AcrossB(2) = A(3)*B(1) - A(1)*B(3)
    AcrossB(3) = A(1)*B(2) - A(2)*B(1)
  return
end

function rot_tmp(x,A,dA) result(rotA)
implicit none
real*8, intent(in) :: x(3), A(3), dA(3,3)
real*8             :: rotA(3)
  rotA(1) = dA(3,2) - dA(2,3) / x(1)
  rotA(2) = dA(1,3) - dA(3,1) - A(3) / x(1)
  rotA(3) = dA(2,1) - dA(1,2) 
return
end    

subroutine rk4_step(x, vpar, qom, zmu, E, B, Bnorm, dBnorm, dB, delta_x, delta_u)
implicit none
real*8 :: x(3), vpar, qom, zmu, A(3), dA(3,3), E(3), Bnorm(3), dBnorm(3,3), B(3), dB(3)
real*8 :: Bstar(3), Estar(3), Bpar_star, delta_x(3), delta_u
    
  Bstar     = B + vpar * rot_tmp(x,Bnorm,dBnorm) / qom
  Bpar_star = dot_product(Bstar,Bnorm)
  Estar     = E(3) - zmu * dB /qom
        
  delta_x = (Bstar * vpar  - cross(Bnorm, Estar)) / Bpar_star
  delta_u = dot_product(Bstar,Estar) * qom        / Bpar_star

  delta_x(3) = delta_x(3) / x(1)

return
end

end module
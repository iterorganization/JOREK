!> mod_particle_assert_equal contains variables and procedure
!> for checking equalities between particles
module mod_particle_assert_equal
use fruit
use mod_assert_equals_tools, only: assert_equals_rel_error  
use mod_particle_types, only: particle_base
implicit none

private
public :: assert_equal_particle, assert_equal_rel_error_particle, set_tol, reset_tol_real4, reset_tol_real8

!> Variables ------------------------------------------------------!
!real*4,parameter :: tol_real4=real(1.d-5,kind=4)
real*4 :: tol_real4=real(1.d-5,kind=4)
!real*8,parameter :: tol_real8=1.d-15
real*8 :: tol_real8=1.d-15

!> Interfaces -----------------------------------------------------

interface assert_equal_particle
  module procedure assert_equal_particle_single
  module procedure assert_equal_particle_list 
end interface assert_equal_particle

interface assert_equal_rel_error_particle
  module procedure assert_equal_rel_error_particle_single
  module procedure assert_equal_rel_error_particle_list 
end interface assert_equal_rel_error_particle

interface set_tol
   module procedure set_tol_real4
   module procedure set_tol_real8
end interface set_tol
   
contains

!> Procedures -----------------------------------------------------

  !> Set the 64 bit tolerance
  subroutine set_tol_real8(tol)
    implicit none
    real*8, intent(in) :: tol
    tol_real8 = tol
  end subroutine set_tol_real8

  !> Set the 32 bit tolerance
  subroutine set_tol_real4(tol)
    implicit none
    real*4, intent(in) :: tol
    tol_real4 = tol
  end subroutine set_tol_real4
 
  !> reset the 64 bit tolerance
  subroutine reset_tol_real8()
    implicit none
    tol_real8 = 1.d-15
  end subroutine reset_tol_real8
  
  !> reset the 32 bit tolerance
  subroutine reset_tol_real4()
    implicit none
    tol_real4=real(1.d-5,kind=4)
  end subroutine reset_tol_real4
  
!> compare two particle lists using absolute error
subroutine assert_equal_particle_list(n_particles,particle_list_1,particle_list_2)
  implicit none
  class(particle_base),dimension(n_particles),intent(in) :: particle_list_1
  class(particle_base),dimension(n_particles),intent(in) :: particle_list_2
  integer,intent(in) :: n_particles
  integer :: ii
#ifndef CUDA_KERNELS
!  !$omp parallel do default(private) shared(n_particles,&
!  !$omp particle_list_1,particle_list_2)
#endif
  do ii=1,n_particles
    call assert_equal_particle_single(particle_list_1(ii),particle_list_2(ii))
 enddo
#ifndef CUDA_KERNELS
! !$omp end parallel do
#endif
end subroutine assert_equal_particle_list

!> compare two particle lists using relative error
subroutine assert_equal_rel_error_particle_list(n_particles,particle_list_1,particle_list_2)
  implicit none
  class(particle_base),dimension(n_particles),intent(in) :: particle_list_1
  class(particle_base),dimension(n_particles),intent(in) :: particle_list_2
  integer,intent(in) :: n_particles
  integer :: ii
#ifndef CUDA_KERNELS
!  !$omp parallel do default(private) shared(n_particles,&
!  !$omp particle_list_1,particle_list_2)
#endif
  do ii=1,n_particles
    call assert_equal_rel_error_particle_single(particle_list_1(ii),particle_list_2(ii))
 enddo
#ifndef CUDA_KERNELS
!  !$omp end parallel do
#endif
end subroutine assert_equal_rel_error_particle_list

!> compare two particles using absolute error
subroutine assert_equal_particle_single(particle_1,particle_2)
  implicit none
  class(particle_base),intent(in) :: particle_1,particle_2
  logical,dimension(8) :: lfails
  call assert_equal_particle_fieldline(particle_1,particle_2,lfails(1)) 
  call assert_equal_particle_gc(particle_1,particle_2,lfails(2))
  call assert_equal_particle_gc_vpar(particle_1,particle_2,lfails(3))
  call assert_equal_particle_gc_Qin(particle_1,particle_2,lfails(4))
  call assert_equal_particle_kinetic(particle_1,particle_2,lfails(5))
  call assert_equal_particle_kinetic_leapfrog(particle_1,particle_2,lfails(6))
  call assert_equal_particle_kinetic_relativistic(particle_1,particle_2,lfails(7))
  call assert_equal_particle_gc_relativistic(particle_1,particle_2,lfails(8))
  if(all(lfails))  call assert_equal_particle_base(particle_1,particle_2)
end subroutine assert_equal_particle_single

!> compare two particles using relative error
subroutine assert_equal_rel_error_particle_single(particle_1,particle_2)
  implicit none
  class(particle_base),intent(in) :: particle_1,particle_2
  logical,dimension(8) :: lfails
  logical, parameter   :: rel_error = .true.
  
  call assert_equal_particle_base(particle_1,particle_2,rel_error)
  call assert_equal_particle_fieldline(particle_1,particle_2,lfails(1),rel_error) 
  call assert_equal_particle_gc(particle_1,particle_2,lfails(2),rel_error)
  call assert_equal_particle_gc_vpar(particle_1,particle_2,lfails(3),rel_error)
  call assert_equal_particle_gc_Qin(particle_1,particle_2,lfails(4),rel_error)
  call assert_equal_particle_kinetic(particle_1,particle_2,lfails(5),rel_error)
  call assert_equal_particle_kinetic_leapfrog(particle_1,particle_2,lfails(6),rel_error)
  call assert_equal_particle_kinetic_relativistic(particle_1,particle_2,lfails(7),rel_error)
  call assert_equal_particle_gc_relativistic(particle_1,particle_2,lfails(8),rel_error)
!  if(all(lfails))  call assert_equal_particle_base(particle_1,particle_2,rel_error)
end subroutine assert_equal_rel_error_particle_single

!> compare particle_gc_relativistic
subroutine assert_equal_particle_gc_relativistic(particle_1,particle_2,lfail,rel_error)
  use mod_particle_types, only: particle_gc_relativistic
  implicit none
  class(particle_base),intent(in) :: particle_1,particle_2
  logical,intent(out)          :: lfail
  logical,optional,intent(in)  :: rel_error
  logical                      :: check_rel_error

  check_rel_error = .false.
  if(present(rel_error)) check_rel_error = rel_error
  lfail = .true.
  select type(p_1=>particle_1)
    type is (particle_gc_relativistic)
    select type (p_2=>particle_2)
      type is (particle_gc_relativistic)
         if(check_rel_error) then
            call assert_equals_rel_error(2,p_1%p,p_2%p,tol_real8,&   
                 "Error particle gc relativistic: momenta p mismatch!")
         else
            call assert_equals(p_1%p,p_2%p,2,tol_real8,&   
                 "Error particle gc relativistic: momenta p mismatch!")
         endif
         call assert_equals(int(p_1%q,kind=4),int(p_2%q,kind=4),&
              "Error particle gc relativistic: charge q mismatch!")
      lfail = .false.
    end select
  end select
end subroutine assert_equal_particle_gc_relativistic

!> compare particle_kinetic_relativistic
subroutine assert_equal_particle_kinetic_relativistic(particle_1,particle_2,lfail,rel_error)
  use mod_particle_types, only: particle_kinetic_relativistic
  implicit none
  class(particle_base),intent(in) :: particle_1,particle_2
  logical,intent(out)           :: lfail
  logical,optional, intent(in)  :: rel_error
  logical                       :: check_rel_error

  check_rel_error = .false.
  if(present(rel_error)) check_rel_error = rel_error
  lfail = .true.
  select type (p_1=>particle_1)
    type is (particle_kinetic_relativistic)
    select type (p_2=>particle_2)
      type is (particle_kinetic_relativistic)
      if(check_rel_error) then 
         call assert_equals_rel_error(3,p_1%p,p_2%p,tol_real8,&               
              "Error particle kinetic relativistic: momentum p mismatch!")
      else
         call assert_equals(p_1%p,p_2%p,3,tol_real8,&
              "Error particle kinetic relativistic: momentum p mismatch!")
      endif
      call assert_equals(int(p_1%q,kind=4),int(p_2%q,kind=4),&
           "Error particle kinetic relativistic: charge q mismatch!")
      lfail = .false.
    end select
  end select
end subroutine assert_equal_particle_kinetic_relativistic

!> compare particle_kinetic_leapfrog
subroutine assert_equal_particle_kinetic_leapfrog(particle_1,particle_2,lfail,rel_error)
  use mod_particle_types, only: particle_kinetic_leapfrog
  implicit none
  class(particle_base),intent(in) :: particle_1,particle_2
  logical,intent(out)           :: lfail
  logical,optional, intent(in)  :: rel_error
  logical                       :: check_rel_error

  check_rel_error = .false.
  if(present(rel_error)) check_rel_error = rel_error
  lfail = .true.
  select type (p_1=>particle_1)
    type is (particle_kinetic_leapfrog)
    select type (p_2=>particle_2)
      type is (particle_kinetic_leapfrog)
      if(check_rel_error) then 
         call assert_equals_rel_error(3,p_1%v,p_2%v,tol_real8,&
              "Error particle kinetic leapfrog: velocity v mismatch!")
      else
         call assert_equals(p_1%v,p_2%v,3,tol_real8,&
              "Error particle kinetic leapfrog: velocity v mismatch!")
      endif
      call assert_equals(int(p_1%q,kind=4),int(p_2%q,kind=4),&
           "Error particle kinetic leapfrog: charge q mismatch!")
      lfail = .false.
    end select
  end select
end subroutine assert_equal_particle_kinetic_leapfrog

!> compare particle_kinetic
subroutine assert_equal_particle_kinetic(particle_1,particle_2,lfail,rel_error)
  use mod_particle_types, only: particle_kinetic
  implicit none
  class(particle_base),intent(in) :: particle_1,particle_2
  logical,intent(out)           :: lfail
  logical,optional, intent(in)  :: rel_error
  logical                       :: check_rel_error

  check_rel_error = .false.
  if(present(rel_error)) check_rel_error = rel_error
  lfail = .true.
  select type (p_1=>particle_1)
    type is (particle_kinetic)
    select type (p_2=>particle_2)
      type is (particle_kinetic)
      if(check_rel_error) then 
         call assert_equals_rel_error(3,p_1%v,p_2%v,tol_real8,&
              "Error particle kinetic: velocity v mismatch!")
      else
         call assert_equals(p_1%v,p_2%v,3,tol_real8,&
              "Error particle kinetic: velocity v mismatch!")
      endif
      call assert_equals(int(p_1%q,kind=4),int(p_2%q,kind=4),&
           "Error particle kinetic: charge q mismatch!")
      lfail = .false.
    end select
  end select
end subroutine assert_equal_particle_kinetic

!> compare particle_gc_Qin
subroutine assert_equal_particle_gc_Qin(particle_1,particle_2,lfail,rel_error)
  use mod_particle_types, only: particle_gc_Qin
  implicit none
  class(particle_base),intent(in) :: particle_1,particle_2
  logical,intent(out)           :: lfail
  logical,optional, intent(in)  :: rel_error
  logical                       :: check_rel_error

  check_rel_error = .false.
  if(present(rel_error)) check_rel_error = rel_error
  lfail = .true.
  select type (p_1=>particle_1)
    type is (particle_gc_Qin)
      select type (p_2=>particle_2)
        type is (particle_gc_Qin)
        if(check_rel_error) then 
           call assert_equals_rel_error(p_1%vpar,p_2%vpar,tol_real8,&
                "Error particle gc Qin: velocity vpar mismatch!")
           call assert_equals_rel_error(p_1%mu,p_2%mu,tol_real8,&
                "Error particle gc Qin: moment mu mismatch!")
           call assert_equals_rel_error(3,p_1%x_m,p_2%x_m,tol_real8,&
                "Error particle gc Qin: x_m position mismatch!")
           call assert_equals_rel_error(p_1%vpar_m,p_2%vpar_m,tol_real8,&
                "Error particle gc Qin: vpar_m velocity mismatch!") 
           call assert_equals_rel_error(3,p_1%Astar_m,p_2%Astar_m,tol_real8,&
                "Error particle gc Qin: Astar_m potential mismatch!")
           call assert_equals_rel_error(3,p_1%Astar_k,p_2%Astar_k,tol_real8,&
                "Error particle gc Qin: Astar_k potential mismatch!")
           call assert_equals_rel_error(3,3,p_1%dAstar_k,p_2%dAstar_k,tol_real8,&
                "Error particle gc Qin: dAstar_k potential der. mismatch!")
           call assert_equals_rel_error(p_1%Bn_k,p_2%Bn_k,tol_real8,&
                "Error particle gc Qin: Bn_k intensity mismatch!")
           call assert_equals_rel_error(3,p_1%dBn_k,p_2%dBn_k,tol_real8,&
                "Error particle gc Qin: dBn_k intensity der. mismatch!")
           call assert_equals_rel_error(3,p_1%Bnorm_k,p_2%Bnorm_k,tol_real8,&
                "Error particle gc Qin: Bnorm_k field mismatch!")
           call assert_equals_rel_error(3,p_1%E_k,p_2%E_k,tol_real8,&
                "Error particle gc Qin: E_k field mismatch!")
           lfail = .false.
        else
           call assert_equals(p_1%vpar,p_2%vpar,tol_real8,&
                "Error particle gc Qin: velocity vpar mismatch!")
           call assert_equals(p_1%mu,p_2%mu,tol_real8,&
                "Error particle gc Qin: moment mu mismatch!")
           call assert_equals(p_1%x_m,p_2%x_m,3,tol_real8,&
                "Error particle gc Qin: x_m position mismatch!")
           call assert_equals(p_1%vpar_m,p_2%vpar_m,tol_real8,&
                "Error particle gc Qin: vpar_m velocity mismatch!") 
           call assert_equals(p_1%Astar_m,p_2%Astar_m,3,tol_real8,&
                "Error particle gc Qin: Astar_m potential mismatch!")
           call assert_equals(p_1%Astar_k,p_2%Astar_k,3,tol_real8,&
                "Error particle gc Qin: Astar_k potential mismatch!")
           call assert_equals(p_1%dAstar_k,p_2%dAstar_k,3,3,tol_real8,&
                "Error particle gc Qin: dAstar_k potential der. mismatch!")
           call assert_equals(p_1%Bn_k,p_2%Bn_k,tol_real8,&
                "Error particle gc Qin: Bn_k intensity mismatch!")
           call assert_equals(p_1%dBn_k,p_2%dBn_k,3,tol_real8,&
                "Error particle gc Qin: dBn_k intensity der. mismatch!")
           call assert_equals(p_1%Bnorm_k,p_2%Bnorm_k,3,tol_real8,&
                "Error particle gc Qin: Bnorm_k field mismatch!")
           call assert_equals(p_1%E_k,p_2%E_k,3,tol_real8,&
                "Error particle gc Qin: E_k field mismatch!")
        endif
        call assert_equals(int(p_1%q,kind=4),int(p_2%q,kind=4),&
             "Error particle gc Qin: charge q mismatch!")
     end select
  end select
end subroutine assert_equal_particle_gc_Qin

!> compare particle_gc_vpar
subroutine assert_equal_particle_gc_vpar(particle_1,particle_2,lfail,rel_error)
  use mod_particle_types, only: particle_gc_vpar
  implicit none
  class(particle_base),intent(in) :: particle_1,particle_2
  logical,intent(out)           :: lfail
  logical,optional, intent(in)  :: rel_error
  logical                       :: check_rel_error

  check_rel_error = .false.
  if(present(rel_error)) check_rel_error = rel_error
  lfail = .true.
  select type (p_1=>particle_1)
    type is (particle_gc_vpar)
    select type (p_2=>particle_2)
      type is (particle_gc_vpar)
      if(check_rel_error) then 
         call assert_equals_rel_error(p_1%vpar,p_2%vpar,tol_real8,&
              "Error particle gc vpar: parallel velocity vpar mismatch!")
         call assert_equals_rel_error(p_1%mu,p_2%mu,tol_real8,&
              "Error particle gc vpar: magnetic moment mu mistmatch!")
      else
         call assert_equals(p_1%vpar,p_2%vpar,tol_real8,&
              "Error particle gc vpar: parallel velocity vpar mismatch!")
         call assert_equals(p_1%mu,p_2%mu,tol_real8,&
              "Error particle gc vpar: magnetic moment mu mistmatch!")
      endif
      call assert_equals(int(p_1%q,kind=4),int(p_2%q,kind=4),&
           "Error particle gc vpar: charge q mistmatch!")
      lfail = .false.
    end select
  end select
end subroutine assert_equal_particle_gc_vpar

!> compare particle_gc
subroutine assert_equal_particle_gc(particle_1,particle_2,lfail,rel_error)
  use mod_particle_types, only: particle_gc
  implicit none
  class(particle_base),intent(in) :: particle_1,particle_2
  logical,intent(out) :: lfail
  logical,optional, intent(in)  :: rel_error
  logical                       :: check_rel_error

  check_rel_error = .false.
  if(present(rel_error)) check_rel_error = rel_error
  lfail = .true.
  select type (p_1=>particle_1)
    type is (particle_gc)
    select type (p_2=>particle_2)
      type is (particle_gc)
      if(check_rel_error) then 
         call assert_equals_rel_error(p_1%E,p_2%E,tol_real8,&
              "Error particle gc: energy E mistmatch!")
         call assert_equals_rel_error(p_1%mu,p_2%mu,tol_real8,&
              "Error particle gc: magnetic moment mu mistmatch!")
      else
         call assert_equals(p_1%E,p_2%E,tol_real8,&
           "Error particle gc: energy E mistmatch!")
         call assert_equals(p_1%mu,p_2%mu,tol_real8,&
              "Error particle gc: magnetic moment mu mistmatch!")
      endif
      call assert_equals(int(p_1%q,kind=4),int(p_2%q,kind=4),&
           "Error particle gc: charge q mistmatch!")
      lfail = .false.
    end select
  end select
end subroutine assert_equal_particle_gc

!> compare particle fieldlines
subroutine assert_equal_particle_fieldline(particle_1,particle_2,lfail,rel_error)
  use mod_particle_types, only: particle_fieldline
  implicit none
  class(particle_base),intent(in) :: particle_1,particle_2
  logical,intent(out)           :: lfail
  logical,optional, intent(in)  :: rel_error
  logical                       :: check_rel_error

  check_rel_error = .false.
  if(present(rel_error)) check_rel_error = rel_error
  lfail = .true.
  select type (p_1=>particle_1)
    type is (particle_fieldline)
    select type (p_2=>particle_2)
      type is (particle_fieldline)
      if(check_rel_error) then 
         !> B_hat is not stored in HDF5
         call assert_equals_rel_error(3,p_1%B_hat_prev,p_2%B_hat_prev,&
              tol_real8,"Error particle fieldline: B_hat_prev mistmatch!")
         call assert_equals_rel_error(p_1%v,p_2%v,tol_real8,&
              "Error particle fieldline: v mistmatch!")
      else
         !> B_hat is not stored in HDF5
         call assert_equals(p_1%B_hat_prev,p_2%B_hat_prev,3,&
              tol_real8,"Error particle fieldline: B_hat_prev mistmatch!")
         call assert_equals(p_1%v,p_2%v,tol_real8,&
              "Error particle fieldline: v mistmatch!")
      endif
      lfail = .false.
    end select
  end select
end subroutine assert_equal_particle_fieldline

!> compare particle base class
subroutine assert_equal_particle_base(particle_1,particle_2,rel_error)
  implicit none
  !> inputs
  class(particle_base),intent(in) :: particle_1,particle_2
  logical,optional, intent(in)  :: rel_error
  logical                       :: check_rel_error

  check_rel_error = .false.
  if(present(rel_error)) check_rel_error = rel_error
  !> check particle relative error
  if(check_rel_error) then 
     call assert_equals_rel_error(3,particle_1%x,particle_2%x,tol_real8,&
          "Error particle base: x position relative mismatch!")
     call assert_equals_rel_error(2,particle_1%st,particle_2%st,tol_real8,&
          "Error particle base: st position relative mistmatch!")
     call assert_equals_rel_error(particle_1%weight,particle_2%weight,tol_real8,&
          "Error particle base: weight relative mistmatch!")
  else   !> check particle absolute error
     call assert_equals(particle_1%x,particle_2%x,3,tol_real8,&
          "Error particle base: x position mismatch!")
     call assert_equals(particle_1%st,particle_2%st,2,tol_real8,&
          "Error particle base: st position mistmatch!")
     call assert_equals(particle_1%weight,particle_2%weight,tol_real8,&
          "Error particle base: weight mistmatch!")
  endif
  call assert_equals(particle_1%i_elm,particle_2%i_elm,&
       "Error particle base: i_elm element index mistmatch!")
  call assert_equals(particle_1%i_life,particle_2%i_life,&
       "Error particle base: i_life index mistmatch!")
  call assert_equals(particle_1%t_birth,particle_2%t_birth,tol_real4,&
       "Error particle base: t_birth time birth mistmatch!")
end subroutine assert_equal_particle_base

!>-----------------------------------------------------------------
end module mod_particle_assert_equal

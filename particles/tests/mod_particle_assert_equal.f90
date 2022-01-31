!> mod_particle_assert_equal contains variables and procedure
!> for checking equalities between particles
module mod_particle_assert_equal
use fruit
use mod_particle_types, only: particle_base
implicit none

private
public :: assert_equal_particle

!> Variables ------------------------------------------------------
real*4,parameter :: tol_real4=real(1.d-5,kind=4)
real*8,parameter :: tol_real8=1.d-15

!> Interfaces -----------------------------------------------------

interface assert_equal_particle
  module procedure assert_equal_particle_single
  module procedure assert_equal_particle_list 
end interface assert_equal_particle

contains

!> Procedures -----------------------------------------------------

!> compare two particle lists
subroutine assert_equal_particle_list(n_particles,particle_list_1,particle_list_2)
  implicit none
  class(particle_base),dimension(n_particles),intent(in) :: particle_list_1
  class(particle_base),dimension(n_particles),intent(in) :: particle_list_2
  integer,intent(in) :: n_particles
  integer :: ii
  !$omp parallel do default(private) shared(n_particles,&
  !$omp particle_list_1,particle_list_2)
  do ii=1,n_particles
    call assert_equal_particle_single(particle_list_1(ii),particle_list_2(ii))
  enddo
  !$omp end parallel do
end subroutine assert_equal_particle_list

!> compare two particles 
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

!> compare particle_gc_relativistic
subroutine assert_equal_particle_gc_relativistic(particle_1,particle_2,lfail)
  use mod_particle_types, only: particle_gc_relativistic
  implicit none
  class(particle_base),intent(in) :: particle_1,particle_2
  logical,intent(out) :: lfail
  lfail = .true.
  call assert_equal_particle_base(particle_1,particle_2)
  select type(p_1=>particle_1)
    type is (particle_gc_relativistic)
    select type (p_2=>particle_2)
      type is (particle_gc_relativistic)
      call assert_equals(p_1%p,p_2%p,2,tol_real8,&
      "Error particle gc relativistic: momenta p mismatch!")
      call assert_equals(int(p_1%q,kind=4),int(p_2%q,kind=4),&
      "Error particle gc relativistic: charge q mismatch!")
      lfail = .false.
    end select
  end select
end subroutine assert_equal_particle_gc_relativistic

!> compare particle_kinetic_relativistic
subroutine assert_equal_particle_kinetic_relativistic(particle_1,particle_2,lfail)
  use mod_particle_types, only: particle_kinetic_relativistic
  implicit none
  class(particle_base),intent(in) :: particle_1,particle_2
  logical,intent(out) :: lfail
  lfail = .true.
  call assert_equal_particle_base(particle_1,particle_2)
  select type (p_1=>particle_1)
    type is (particle_kinetic_relativistic)
    select type (p_2=>particle_2)
      type is (particle_kinetic_relativistic)
      call assert_equals(p_1%p,p_2%p,3,tol_real8,&
      "Error particle kinetic relativistic: momentum p mismatch!")
      call assert_equals(int(p_1%q,kind=4),int(p_2%q,kind=4),&
      "Error particle kinetic relativistic: charge q mismatch!")
      lfail = .false.
    end select
  end select
end subroutine assert_equal_particle_kinetic_relativistic

!> compare particle_kinetic_leapfrog
subroutine assert_equal_particle_kinetic_leapfrog(particle_1,particle_2,lfail)
  use mod_particle_types, only: particle_kinetic_leapfrog
  implicit none
  class(particle_base),intent(in) :: particle_1,particle_2
  logical,intent(out) :: lfail
  lfail = .true.
  call assert_equal_particle_base(particle_1,particle_2)
  select type (p_1=>particle_1)
    type is (particle_kinetic_leapfrog)
    select type (p_2=>particle_2)
      type is (particle_kinetic_leapfrog)
      call assert_equals(p_1%v,p_2%v,3,tol_real8,&
      "Error particle kinetic leapfrog: velocity v mismatch!")
      call assert_equals(int(p_1%q,kind=4),int(p_2%q,kind=4),&
      "Error particle kinetic leapfrog: charge q mismatch!")
      lfail = .false.
    end select
  end select
end subroutine assert_equal_particle_kinetic_leapfrog

!> compare particle_kinetic
subroutine assert_equal_particle_kinetic(particle_1,particle_2,lfail)
  use mod_particle_types, only: particle_kinetic
  implicit none
  class(particle_base),intent(in) :: particle_1,particle_2
  logical,intent(out) :: lfail
  lfail = .true.
  call assert_equal_particle_base(particle_1,particle_2)
  select type (p_1=>particle_1)
    type is (particle_kinetic)
    select type (p_2=>particle_2)
      type is (particle_kinetic)
      call assert_equals(p_1%v,p_2%v,3,tol_real8,&
      "Error particle kinetic: velocity v mismatch!")
      call assert_equals(int(p_1%q,kind=4),int(p_2%q,kind=4),&
      "Error particle kinetic: charge q mismatch!")
      lfail = .false.
    end select
  end select
end subroutine assert_equal_particle_kinetic

!> compare particle_gc_Qin
subroutine assert_equal_particle_gc_Qin(particle_1,particle_2,lfail)
  use mod_particle_types, only: particle_gc_Qin
  implicit none
  class(particle_base),intent(in) :: particle_1,particle_2
  logical,intent(out) :: lfail
  lfail = .true.
  call assert_equal_particle_base(particle_1,particle_2)
  select type (p_1=>particle_1)
    type is (particle_gc_Qin)
      select type (p_2=>particle_2)
        type is (particle_gc_Qin)
        call assert_equals(p_1%vpar,p_2%vpar,tol_real8,&
        "Error particle gc Qin: velocity vpar mismatch!")
        call assert_equals(p_1%mu,p_2%mu,tol_real8,&
        "Error particle gc Qin: moment mu mismatch!")
        call assert_equals(int(p_1%q,kind=4),int(p_2%q,kind=4),&
        "Error particle gc Qin: charge q mismatch!")
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
        lfail = .false.
      end select
  end select
end subroutine assert_equal_particle_gc_Qin

!> compare particle_gc_vpar
subroutine assert_equal_particle_gc_vpar(particle_1,particle_2,lfail)
  use mod_particle_types, only: particle_gc_vpar
  implicit none
  class(particle_base),intent(in) :: particle_1,particle_2
  logical,intent(out) :: lfail
  lfail = .true.
  call assert_equal_particle_base(particle_1,particle_2)
  select type (p_1=>particle_1)
    type is (particle_gc_vpar)
    select type (p_2=>particle_2)
      type is (particle_gc_vpar)
      call assert_equals(p_1%vpar,p_2%vpar,tol_real8,&
      "Error particle gc vpar: parallel velocity vpar mismatch!")
      call assert_equals(p_1%mu,p_2%mu,tol_real8,&
      "Error particle gc vpar: magnetic moment mu mistmatch!")
      call assert_equals(int(p_1%q,kind=4),int(p_2%q,kind=4),&
      "Error particle gc vpar: charge q mistmatch!")
      lfail = .false.
    end select
  end select
end subroutine assert_equal_particle_gc_vpar

!> compare particle_gc
subroutine assert_equal_particle_gc(particle_1,particle_2,lfail)
  use mod_particle_types, only: particle_gc
  implicit none
  class(particle_base),intent(in) :: particle_1,particle_2
  logical,intent(out) :: lfail
  lfail = .true.
  call assert_equal_particle_base(particle_1,particle_2)
  select type (p_1=>particle_1)
    type is (particle_gc)
    select type (p_2=>particle_2)
      type is (particle_gc)
      call assert_equals(p_1%E,p_2%E,tol_real8,&
      "Error particle gc: energy E mistmatch!")
      call assert_equals(p_1%mu,p_2%mu,tol_real8,&
      "Error particle gc: magnetic moment mu mistmatch!")
      call assert_equals(int(p_1%q,kind=4),int(p_2%q,kind=4),&
      "Error particle gc: charge q mistmatch!")
      lfail = .false.
    end select
  end select
end subroutine assert_equal_particle_gc

!> compare particle fieldlines
subroutine assert_equal_particle_fieldline(particle_1,particle_2,lfail)
  use mod_particle_types, only: particle_fieldline
  implicit none
  class(particle_base),intent(in) :: particle_1,particle_2
  logical,intent(out) :: lfail
  lfail = .true.
  call assert_equal_particle_base(particle_1,particle_2)
  select type (p_1=>particle_1)
    type is (particle_fieldline)
    select type (p_2=>particle_2)
      type is (particle_fieldline)
      !> B_hat is not stored in HDF5
      call assert_equals(p_1%B_hat_prev,p_2%B_hat_prev,3,&
      tol_real8,"Error particle fieldline: B_hat_prev mistmatch!")
      call assert_equals(p_1%v,p_2%v,tol_real8,&
      "Error particle fieldline: v mistmatch!")
      lfail = .false.
    end select
  end select
end subroutine assert_equal_particle_fieldline

!> compare particle base class
subroutine assert_equal_particle_base(particle_1,particle_2)
 implicit none
 !> inputs
 class(particle_base),intent(in) :: particle_1,particle_2
 !> check particle 
 call assert_equals(particle_1%x,particle_2%x,3,tol_real8,&
 "Error particle base: x position mistmatch!")
 call assert_equals(particle_1%st,particle_2%st,2,tol_real8,&
 "Error particle base: st position mistmatch!")
 call assert_equals(particle_1%weight,particle_2%weight,tol_real8,&
 "Error particle base: weight mistmatch!")
 call assert_equals(particle_1%i_elm,particle_2%i_elm,&
 "Error particle base: i_elm element index mistmatch!")
 call assert_equals(particle_1%i_life,particle_2%i_life,&
 "Error particle base: i_life index mistmatch!")
 call assert_equals(particle_1%t_birth,particle_2%t_birth,tol_real4,&
 "Error particle base: t_birth time birth mistmatch!")
end subroutine assert_equal_particle_base

!>-----------------------------------------------------------------
end module mod_particle_assert_equal

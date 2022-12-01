! This is a program to output the radiation power function and its components for a range of hard-coded parameters.
! The species and the ADAS data need to be specified below for the program to work properly.
program radiation_function_diagno
use mpi_mod
use mod_openadas
use mod_coronal
use constants
implicit none

real*8, parameter :: n_e = 1d19 ! m^-3
real*8, parameter :: T_min = 0.1, T_max = 50000 ! eV
integer, parameter :: n_t = 10000
integer :: i, j

integer :: i_t, ierr, provided
real*8 :: T_e, maxdensity
real*8 :: temp, GRCout
real*8, allocatable :: PLT(:), PRB(:)
type(adf11_all) :: ad
type(coronal) :: cor
character*120 :: directory

real*8,allocatable    :: P_imp(:), dP_imp_dT(:)
real*8                :: Z_imp, dZ_imp_dT, d2Z_imp_dT2, Prad, dPrad_dT

!call MPI_Init_thread(MPI_THREAD_MULTIPLE, provided, ierr)

directory = '${JOREK_ADAS_DIR}'
ad = read_adf11(0, '96_ne', directory)
cor = coronal(ad)

allocate(P_imp(0:ad%n_Z), dP_imp_dT(0:ad%n_Z))

!allocate(PLT(0:ad%n_Z))
!allocate(PRB(0:ad%n_Z))
!PLT = 0.d0
!PRB = 0.d0

maxdensity = ubound(cor%Prad,1)

open(1, file = 'T_Prad.dat', status = 'replace')  
write(1,*) 'Temperature (eV)'
write(1,*) 10**(cor%temperature)*K_BOLTZ/EL_CHG
!write(1,*) 'Density (m^-3)'
!write(1,*) 10**(cor%density)
write(1,*) 'Radiated power (Wm^3)' !Radiation per ion density per electron density
write(1,*) cor%Prad(maxdensity, :)/(10**(cor%density(maxdensity)))
close(1)


open(2, file = 'GRC_ACD.dat', status = 'replace')
do i=0,ad%n_Z
  write(2,*) 'Charge state ',i
  write(2,*) 'Real'
  do j=1,size(ad%ACD%temperature,1)
    write(2,*) ad%ACD%GRC(maxdensity,j,i), ad%ACD%temperature(j)
  end do
  write(2,*) 'Spline'
  do i_t=1,n_t
    temp = T_min * (T_max/T_min)**(real(i_t - 1)/real(n_t - 1))
    call ad%ACD%interp(i,maxdensity,log10(temp*EL_CHG/K_BOLTZ),GRCout)     !ad%PLT%GRCFspline
    write(2,*) GRCout, temp
  end do
end do
close(2)



open(3, file = 'GRC_SCD.dat', status = 'replace')
do i=0,ad%n_Z
  write(3,*) 'Charge state ',i
  write(3,*) 'Real'
  do j=1,30
    write(3,*) ad%SCD%GRC(maxdensity,j,i), ad%SCD%temperature(j) 
  end do
  write(3,*) 'Spline'
  do i_t=1,n_t
    temp = T_min * (T_max/T_min)**(real(i_t - 1)/real(n_t - 1))
    call ad%SCD%interp(i,maxdensity,log10(temp*EL_CHG/K_BOLTZ),GRCout)     !ad%PLT%GRCFspline
    write(3,*) GRCout, temp
  end do
end do
close(3)

open(4, file = 'GRC_PRB.dat', status = 'replace')
do i=0,ad%n_Z
  write(4,*) 'Charge state ',i
  write(4,*) 'Real'
  do j=1,30
    write(4,*) ad%PRB%GRC(maxdensity,j,i), ad%PRB%temperature(j) 
  end do
  write(4,*) 'Spline'
  do i_t=1,n_t
    temp = T_min * (T_max/T_min)**(real(i_t - 1)/real(n_t - 1))
    call ad%PRB%interp(i,maxdensity,log10(temp*EL_CHG/K_BOLTZ),GRCout)     !ad%PLT%GRCFspline
    write(4,*) GRCout, temp
  end do
end do
close(4)


open(5, file = 'GRC_PLT.dat', status = 'replace')
do i=0,ad%n_Z
  write(5,*) 'Charge state ',i
  write(5,*) 'Real'
  do j=1,30
    write(5,*) ad%PLT%GRC(maxdensity,j,i), ad%PLT%temperature(j) 
  end do
  write(5,*) 'Spline'
  do i_t=1,n_t
    temp = T_min * (T_max/T_min)**(real(i_t - 1)/real(n_t - 1))
    call ad%PLT%interp(i,maxdensity,log10(temp*EL_CHG/K_BOLTZ),GRCout)     !ad%PLT%GRCFspline
    write(5,*) GRCout, temp
  end do
end do
close(5)

open(2, file = 'dZ_imp_dT.dat', status = 'replace')
do i_t=1,n_t
  temp = T_min * (T_max/T_min)**(real(i_t - 1)/real(n_t - 1))
  call cor%interp_linear(density=20.,temperature=log10(temp*EL_CHG/K_BOLTZ),&
                         p_out=P_imp,p_Te_out=dP_imp_dT,z_avg=Z_imp,z_avg_Te=dZ_imp_dT,&
                         z_avg_TeTe=d2Z_imp_dT2)     !ad%PLT%GRCFspline
  write(2,*) temp, dZ_imp_dT
end do
close(2)

open(2, file = 'Z_imp.dat', status = 'replace')
do i_t=1,n_t
  temp = T_min * (T_max/T_min)**(real(i_t - 1)/real(n_t - 1))
  call cor%interp_linear(density=20.,temperature=log10(temp*EL_CHG/K_BOLTZ),&
                         p_out=P_imp,p_Te_out=dP_imp_dT,z_avg=Z_imp,z_avg_Te=dZ_imp_dT,&
                         z_avg_TeTe=d2Z_imp_dT2)     !ad%PLT%GRCFspline
  write(2,*) temp, Z_imp
end do
close(2)

! logspace in T
!do i_t=1, n_t
!  T_e = T_min * (T_max/T_min)**(real(i_t - 1)/real(n_t - 1))

!  write(*,*) T_e, cor%interp(0, log10(n_e),log10(T_e*EL_CHG/K_BOLTZ))*n_e
!end do

end program radiation_function_diagno

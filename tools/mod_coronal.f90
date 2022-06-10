!> module takes the OPEN-ADAS data to calculate the coronal equil-
!> ibrium temperature and radiation.
!> If you need time-dependent solutions of the corona matrix timestepping look
!> in the revision history for this file.
module mod_coronal
use mod_openadas
use mod_interp_splinear
use constants
implicit none
private
public coronal
public specific_coronal_equilibrium
public coronal_prad

!> Coronal equilibrium datatype
type coronal
  integer :: n_Z                          !< Atomic number
  real*8, allocatable :: density(:)       !< log10 density (m^-3)
  real*8, allocatable :: temperature(:)   !< log10 temperature (K)
  real*8, allocatable :: Z(:,:,:)         !< Charge state (e) density for specific densities, temperatures and charge states [i_n, i_T, i_q]
  real*8, allocatable :: Z_1T(:,:,:)      !< First temperature derivative of charge state (e) density for specific densities, temperatures and charge states [i_n, i_T, i_q]
  real*8, allocatable :: Prad(:,:)        !< log10 Radiated power per ion (W) for the above densities and temperatures [i_n, i_T]
  real*8, allocatable :: Prad_1T(:,:)     !< First temperature gradient of log10 Radiated power per ion (W) for the above densities and temperatures [i_n, i_T]
  real*8, allocatable :: Z_avg_CE(:,:)    !< The average charge of impurity for the above densities and temperatures [i_n, i_T]
  real*8, allocatable :: Z_avg_1T_CE(:,:) !< The first temperature gradient of the average charge for the above densities and temperatures [i_n, i_T]
  real*8, allocatable :: Z_avg_2T_CE(:,:) !< The second temperature gradient of the average charge for the above densities and temperatures [i_n, i_T]

  type(Fspline)       :: ZFspline         !< Spline functions for effective charge
  type(Fspline)       :: PradFspline      !< Spline functions for CE radiation function
  type(Fspline), allocatable :: PFspline(:)  !< Spline functions for each charge state
contains
  procedure :: interp => interpolate_coronal_spl
  procedure :: interp_gradients => interpolate_coronal_gradients
  procedure :: interp_linear => interpolate_coronal
end type coronal
interface coronal
  module procedure coronal_equilibrium
end interface coronal

contains

!> Radiated power in a specific coronal equilibrium configuration and temperature
function coronal_Prad(ad, density, temperature, fractions, neutral_density)
type (ADF11_all), intent(in)            :: ad              !< ADF11 datatype
real*8, intent(in)                      :: density         !< log10 density in m^-3
real*8, intent(in)                      :: temperature     !< log10 electron temperature in K
real*8, intent(in), dimension(0:ad%n_Z) :: fractions       !< Fractional charge states. Should sum to 1 but we do not check it!
real*8, intent(in), optional            :: neutral_density !< log10 neutral density in m^-3
real*8                                  :: coronal_Prad    !< Output power in W / atom


real*8, dimension(0:ad%n_Z) :: rad
real*8, dimension(0:ad%n_Z) :: rad_RC
real*8  :: density_n, rad_RB, rad_LT
integer :: iz

do iz=0,ad%n_Z
  call ad%PRB%interp(iz, density, temperature, rad_RB)
  call ad%PLT%interp(iz, density, temperature, rad_LT)
  call ad%PRC%interp(iz, density, temperature, rad_RC(iz))
  rad(iz) = rad_RB + rad_LT
enddo ! radiation emitted by atoms at level iz
! PRB and PLT should also be multiplied by n_e, and PRC with neutral density

! If the neutral Hydrogen density is present use it
if (present(neutral_density)) then
  density_n = neutral_density
else ! otherwise set it to some extremely low value
  density_n = -99.d0 ! this is log10 of density
endif

coronal_Prad = dot_product(fractions(0:ad%n_Z), rad*10.d0**density + rad_RC*10.d0**density_n)
end function coronal_Prad

!> The gradients of several quantities in a specific coronal equilibrium configuration and temperature
subroutine coronal_gradients(ad, cor, n_d, n_T)
type (ADF11_all), intent(in)            :: ad !< ADF11 datatype
type (coronal), intent(inout)           :: cor !< Coronal equilibrium datatypei
integer, intent(in)                     :: n_d, n_T

real*8, dimension(0:cor%n_Z) :: dp_dT
real*8  :: dPrad_dT, dZ_avg_dT, d2Z_avg_dT2
real*8  :: temperature_cor, density_cor
integer :: iz, m, k

do m = 1, n_d
  do k = 1, n_T

    density_cor     = cor%density(m)
    temperature_cor = cor%temperature(k)

    call cor%interp(density=density_cor,temperature=temperature_cor, p_Te_out=dp_dT,&
                    z_Te_out=dZ_avg_dT, z_TeTe_out=d2Z_avg_dT2, rad_Te_out=dPrad_dT)
    
    cor%Z_1T(m,k,:)       = dp_dT
    cor%Prad_1T(m,k)      = dPrad_dT
    cor%Z_avg_1T_CE(m,k)  = dZ_avg_dT
    cor%Z_avg_2T_CE(m,k)  = d2Z_avg_dT2

  enddo
enddo

end subroutine coronal_gradients


!> Calculate the coronal equilibrium values at specific values of density and temperature
function specific_coronal_equilibrium(ad, density, temperature) result(fractions)
type (ADF11_all), intent(in) :: ad !< ADF11 datatype
real*8, intent(in) :: density !< log10 density in m^-3
real*8, intent(in) :: temperature !< log10 temperature in K
real*8, dimension(0:ad%n_Z) :: fractions

integer :: iz
real*8 :: ion_rate, rec_rate

fractions(0) = 1.d0
do iz=1,ad%n_Z
  call ad%SCD%interp(iz-1, density, temperature, ion_rate) ! ionizing to level iz (0 is neutral)
  call ad%ACD%interp(iz,   density, temperature, rec_rate) ! recombining from iz+1
  fractions(iz) = fractions(iz-1) * ion_rate/rec_rate
end do
fractions = fractions/sum(fractions)
end function specific_coronal_equilibrium



!> Calculate the coronal equilibrium values at specific values of density and temperature
function coronal_equilibrium(ad) result(cor)
use constants
type (ADF11_all), intent(in) :: ad !< ADF11 datatype
type (coronal)               :: cor !< Coronal equilibrium datatype

real*8, dimension(0:ad%n_Z) :: p
integer :: n_d, n_T, iz, k, m
real*8 :: ion_rate, rec_rate
real*8, allocatable :: Z_eff(:,:)

cor%n_Z = ad%n_Z
n_d = 10
n_T = 1000

allocate(cor%density(n_d), cor%temperature(n_T), cor%Z(n_d,n_T,0:cor%n_Z), cor%Prad(n_d,n_T))
allocate(cor%Z_1T(n_d,n_T,0:cor%n_Z),cor%Prad_1T(n_d,n_T))
allocate(cor%Z_avg_CE(n_d,n_T),cor%Z_avg_1T_CE(n_d,n_T),cor%Z_avg_2T_CE(n_d,n_T))
allocate(Z_eff(n_d,n_T))
Z_eff = 0.0

allocate(cor%PFspline(0:cor%n_Z))
do iz=0,ad%n_Z
  call AllocFspline(cor%PFspline(iz),n_T,n_d)
end do
call AllocFspline(cor%ZFspline,n_T,n_d)
call AllocFspline(cor%PradFspline,n_T,n_d)

do m=1, n_d
  cor%density(m) = 18.d0 + real(m-1,8)/real(n_d-1,8) * (21.d0-18.d0) ! log10 [m^-3], linear between 18 and 21
end do
do k=1, n_T
  cor%temperature(k) = log10( 1.d0 + exp(log(4.d4)*float(k-1)/(float(n_T-1))) - 1.d0 ) + log10(EL_CHG) - log10(K_BOLTZ) ! in log10 [K], 1 to 40000 eV in logscale
end do

cor%ZFspline%xspline    = cor%temperature
cor%ZFspline%ylinear    = cor%density
cor%PradFspline%xspline = cor%temperature
cor%PradFspline%ylinear = cor%density

do iz=0,ad%n_Z
  cor%PFspline(iz)%xspline = cor%temperature
  cor%PFspline(iz)%ylinear = cor%density
end do

do m=1, n_d
  do k=1, n_T
    p(0) = 1.d0
    do iz=1,ad%n_Z
      ! The flux of particles from state iz to iz-1 is given by
      ! p(iz-1) * scd(iz-1) - p(iz) * acd(iz)
      call ad%SCD%interp(iz-1, cor%density(m), cor%temperature(k), ion_rate) ! ionizing from level iz-1 to iz
      call ad%ACD%interp(iz, cor%density(m), cor%temperature(k), rec_rate) ! recombining from iz to iz - 1
      p(iz) = p(iz-1) * ion_rate/rec_rate
    end do

    cor%Z(m,k,:) = p / sum(p)
    do iz=1,ad%n_Z
      Z_eff(m,k) = Z_eff(m,k) + cor%Z(m,k,iz) * real(iz,8)
    end do
    cor%Z_avg_CE(m,k) = Z_eff(m,k)
    cor%Prad(m,k) = coronal_Prad(ad, cor%density(m), cor%temperature(k), p/sum(p)) ! Do not set neutral density yet
  enddo

  cor%ZFspline%xspline = cor%temperature

  call spline(n_T,cor%ZFspline%xspline,Z_eff(m,:),0.d0,0.d0,2,&
                  cor%ZFspline%Aspline(m,:),cor%ZFspline%Bspline(m,:),&
                  cor%ZFspline%Cspline(m,:),cor%ZFspline%Dspline(m,:))

  call spline(n_T,cor%PradFspline%xspline,cor%Prad(m,:),0.d0,0.d0,2,&
                  cor%PradFspline%Aspline(m,:),cor%PradFspline%Bspline(m,:),&
                  cor%PradFspline%Cspline(m,:),cor%PradFspline%Dspline(m,:))
  do iz=0,ad%n_Z
    call spline(n_T,cor%PFspline(iz)%xspline,cor%Z(m,:,iz),0.d0,0.d0,2,&
                    cor%PFspline(iz)%Aspline(m,:),cor%PFspline(iz)%Bspline(m,:),&
                    cor%PFspline(iz)%Cspline(m,:),cor%PFspline(iz)%Dspline(m,:))
  end do
enddo

if (any(cor%Z .lt. 0.d0)) then
  write(*,*) "Z prob below zero", count(cor%Z .lt. 0.d0), minval(cor%Z), minloc(cor%Z)
end if

call ConstructFspline(cor%ZFspline,Z_eff)
call ConstructFspline(cor%PradFspline,cor%Prad)
do iz=0,ad%n_Z
  call ConstructFspline(cor%PFspline(iz),cor%Z(:,:,iz))
end do

call coronal_gradients(ad, cor, n_d, n_T)

end function coronal_equilibrium


!> Linear interpolation of coronal model charge at specific density and temperature

subroutine interpolate_coronal(cor, density, temperature, p_out, p_Te_out, &
                               z_avg, z_avg_Te, z_avg_TeTe, rad_out, rad_Te_out)
class(coronal), intent(in)      :: cor         !< Coronal equilibrium type
real*8, intent(in)              :: density     !< log10 density (m^-3)
real*8, intent(in)              :: temperature !< log10 temperature (K)
real*8, intent(out), optional, dimension(0:cor%n_Z) :: p_out, p_Te_out !< distribution of charge states (sum = 1)
real*8, intent(out), optional   :: z_avg, z_avg_Te, z_avg_TeTe !< Average charge according to coronal equilibrium and its derivatives
real*8, intent(out), optional   :: rad_out, rad_Te_out !< radiated power according to coronal equilibrium and its derivatives

real*8, dimension(0:cor%n_Z)    :: p, dp_dT    !< distribution of charge states (sum = 1)
real*8, dimension(0:cor%n_Z)    :: Z           !< The charge number at each charge state
integer                         :: iz

p = L2D2interp(cor%density,cor%temperature,cor%n_Z+1,cor%Z(:,:,:),density,temperature)
dp_dT = L2D2interp(cor%density,cor%temperature,cor%n_Z+1,cor%Z_1T(:,:,:),density,temperature)

if (abs(sum(p)-1.)>1.d-3) then
  write(*,*) "WARNING: Interpolation returns non-unity total fractional abundance, probably approaching ADAS parameter boundary!"
  write(*,*) "sum(p)=",sum(p),", log10(T_e(K))=",temperature,", log10(n_e(m^-3))=",density
endif

if (present(p_out))    p_out    = p/sum(p)
if (present(p_Te_out)) p_Te_out = dp_dT/sum(p)
if (present(z_avg)) then
  z_avg = L2Dinterp(cor%density,cor%temperature,cor%Z_avg_CE(:,:),density,temperature)
endif
if (present(z_avg_Te)) then
  z_avg_Te = L2Dinterp(cor%density,cor%temperature,cor%Z_avg_1T_CE(:,:),density,temperature)
endif
if (present(z_avg_TeTe)) then
  z_avg_TeTe = L2Dinterp(cor%density,cor%temperature,cor%Z_avg_2T_CE(:,:),density,temperature)
endif

!if (present(z_avg)) then
!  do iz=0,cor%n_Z
!    Z(iz) = real(iz,8)
!    if (p(iz)<0.d0) p(iz)=0.d0
!  enddo
!  z_avg = dot_product(p/sum(p),Z)
!endif

if (present(rad_out)) then
  rad_out = L2Dinterp(cor%density,cor%temperature,cor%Prad(:,:),density,temperature)
  !call SL2Dinterp(cor%PradFspline,temperature,density,fout=rad_out)
endif
if (present(rad_Te_out)) then
  rad_Te_out = L2Dinterp(cor%density,cor%temperature,cor%Prad_1T(:,:),density,temperature)
endif
end subroutine interpolate_coronal


!> Linear interpolation of coronal model charge at specific density and temperature.
!> Evaluate the gradients only
subroutine interpolate_coronal_gradients(cor, density, temperature, p_Te_out, p_Ne_out, z_avg_Te, z_avg_Ne)
class(coronal), intent(in)      :: cor !< Coronal equilibrium type
real*8, intent(in)              :: density !< log10 density (m^-3)
real*8, intent(in)              :: temperature !< log10 temperature (K)
real*8, intent(out), optional, dimension(0:cor%n_Z) :: p_Te_out, p_Ne_out !< gradient of distribution of charge states (sum = 1) to Te and Ne
real*8, intent(out), optional   :: z_avg_Te, z_avg_Ne !< effective charge gradient according to coronal equilibrium

real*8, dimension(0:cor%n_Z)    :: p !< distribution of charge states (sum = 1)
real*8, dimension(0:cor%n_Z)    :: p_Te, p_Ne !< gradient of distribution of charge states (sum = 1) to Te and Ne
real*8, dimension(0:cor%n_Z)    :: Z !< The charge number at each charge state
integer                         :: iz

p    = L2D2interp(cor%density,cor%temperature,cor%n_Z+1,cor%Z(:,:,:),density,temperature)
p_Te = L2D2interp_grad(cor%density,cor%temperature,cor%n_Z+1,cor%Z(:,:,:),density,temperature,1)
p_Ne = L2D2interp_grad(cor%density,cor%temperature,cor%n_Z+1,cor%Z(:,:,:),density,temperature,2)

do iz = 0, cor%n_z
  if (p(iz)<0.d0) p(iz)=0.d0
end do

! Converting log gradient to real gradient
p_Te = p_Te / (log(10.d0)*10.d0**temperature)
p_Ne = p_Ne / (log(10.d0)*10.d0**density)

if (abs(sum(p)-1.)>1.d-3) then
  write(*,*) "WARNING: Interpolation returns non-unity total fractional abundance, probably approaching ADAS parameter boundary!"
  write(*,*) "sum(p)=",sum(p),", log10(T_e(K))=",temperature,", log10(n_e(m^-3))=",density
endif

if (present(p_Te_out)) then
  p_Te_out = p_Te/sum(p)
endif

if (present(p_Ne_out)) then
  p_Ne_out = p_Ne/sum(p)
endif

if (present(z_avg_Te)) then
  do iz=0,cor%n_Z
    Z(iz) = real(iz,8)
  enddo
  z_avg_Te =  dot_product(p_Te/sum(p),Z)
endif

if (present(z_avg_Ne)) then
  do iz=0,cor%n_Z
    Z(iz) = real(iz,8)
  enddo
  z_avg_Ne =  dot_product(p_Ne/sum(p),Z)
endif

end subroutine interpolate_coronal_gradients

! Spline-linear interpolation of the coronal equilibrium for both value and
! gradients
subroutine interpolate_coronal_spl(cor, density, temperature, p_out, p_Te_out, p_Ne_out, z_out,&
                                         z_Te_out, z_TeTe_out, z_Ne_out, rad_out, rad_Te_out, rad_Ne_out)
class(coronal), intent(in)      :: cor !< Coronal equilibrium type
real*8, intent(in)              :: density !< log10 density (m^-3)
real*8, intent(in)              :: temperature !< log10 temperature (K)
real*8, intent(out), optional, dimension(0:cor%n_Z) :: p_out, p_Te_out, p_Ne_out !< gradient of distribution of charge states (sum = 1) to Te and Ne
real*8, intent(out), optional   :: z_out, z_Te_out, z_TeTe_out, z_Ne_out !< effective charge gradient according to coronal equilibrium
real*8, intent(out), optional   :: rad_out, rad_Te_out, rad_Ne_out !< density multiplied radiation function and gradients

real*8, dimension(0:cor%n_Z)    :: p !< distribution of charge states (sum = 1)
real*8, dimension(0:cor%n_Z)    :: p_Te, p_Ne !< gradient of distribution of charge states (sum = 1) to Te and Ne
real*8, dimension(0:cor%n_Z)    :: p_TeTe, Z_p
integer                         :: iz
logical                         :: Z_flag !If true, calculate Z_eff from distribution, otherwise interpolate 
real*8                          :: z, z_Te, z_Ne, z_TeTe ! local variables preparing for output
real*8                          :: rad, rad_Te, rad_Ne ! local variables preparing for output

Z_flag = .false.
if (present(p_out) .or. present(p_Te_out) .or. present(p_Ne_out)) then
  do iz = 0, cor%n_z
    call SL2Dinterp(cor%PFspline(iz),temperature,density,fout=p(iz),dfout_dx=p_Te(iz),&
                                                         dfout_dy=p_Ne(iz),d2fout_dx2=p_TeTe(iz))
    if (p(iz)<0.d0) p(iz)=0.d0
  end do

  ! Converting log gradient to real gradient
  p_Te = p_Te / (log(10.d0)*10.d0**temperature)
  p_Ne = p_Ne / (log(10.d0)*10.d0**density)
  p_TeTe = p_TeTe / (log(10.d0)**2 * 10.d0**(2.d0*temperature)) - p_Te/(10.d0**temperature)

  if (abs(sum(p)-1.)>1.d-3) then
    write(*,*) "WARNING: Interpolation returns non-unity total fractional abundance, probably approaching ADAS parameter boundary!"
    write(*,*) "sum(p)=",sum(p),", log10(T_e(K))=",temperature,", log10(n_e(m^-3))=",density
  endif

  if (present(p_out))    p_out    = p/sum(p)
  if (present(p_Te_out)) p_Te_out = p_Te/sum(p)
  if (present(p_Ne_out)) p_Ne_out = p_Ne/sum(p)

  if (temperature < log10(EL_CHG)-log10(K_BOLTZ) .or. temperature > log10(4.d4)+log10(EL_CHG)-log10(K_BOLTZ)) Z_flag = .true.

end if

if (present(z_out) .or. present(z_Te_out) .or. present(z_TeTe_out) .or. present(p_Ne_out)) then
  if (Z_flag) then
    do iz=0,cor%n_Z
      Z_p(iz) = real(iz,8)
    enddo
    z        =  dot_product(p/sum(p),Z_p)
    z_Te     =  dot_product(p_Te/sum(p),Z_p)
    z_Ne     =  dot_product(p_Ne/sum(p),Z_p)
    Z_TeTe   =  dot_product(p_TeTe/sum(p),Z_p)
  else
    call SL2Dinterp(cor%ZFspline,temperature,density,fout=z,dfout_dx=z_Te,dfout_dy=z_Ne,d2fout_dx2=z_TeTe)
  
    ! Converting log gradient to real gradient
    z_Te = z_Te / (log(10.d0)*10.d0**temperature)
    z_Ne = z_Ne / (log(10.d0)*10.d0**density)
    z_TeTe = (z_TeTe / (log(10.d0)**2 * 10.d0**(2.d0*temperature))) - z_Te/(10.d0**temperature)
  end if

  if (present(z_out))      z_out      = z
  if (present(z_Te_out))   z_Te_out   = z_Te
  if (present(z_Ne_out))   z_Ne_out   = z_Ne
  if (present(z_TeTe_out)) z_TeTe_out = z_TeTe
end if

if (present(rad_out) .or. present(rad_Te_out) .or. present(rad_Ne_out)) then
  call SL2Dinterp(cor%PradFspline,temperature,density,fout=rad,dfout_dx=rad_Te,dfout_dy=rad_Ne)

  ! Converting log gradient to real gradient
  rad_Te = rad_Te / (log(10.d0)*10.d0**temperature)
  rad_Ne = rad_Ne / (log(10.d0)*10.d0**density)

  if (present(rad_out))    rad_out    = rad
  if (present(rad_Te_out)) rad_Te_out = rad_Te
  if (present(rad_Ne_out)) rad_Ne_out = rad_Ne
end if

end subroutine interpolate_coronal_spl


end module mod_coronal

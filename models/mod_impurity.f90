module mod_impurity

  use constants

  contains 

  subroutine init_imp_adas(my_id)

    use phys_module
    use mod_openadas
    use mod_coronal
    use mod_parameters, only: with_impurities

    implicit none

    integer, intent(in) :: my_id
    integer             :: i

    character(len=512)  :: adas_suffix     !The suffix of adas data file to be read

    ! Temporary variable for charge state distribution
    integer             :: i_T, i_ion
    real*8, allocatable :: dP_imp_dT(:), P_imp(:)
    real*8              :: Z_imp

    if (allocated(imp_adas)) then
      deallocate(imp_adas)
    end if

    allocate (imp_adas(n_adas))  !< Dynamically allocate memeries for adas data

    if (allocated(imp_cor)) then
      deallocate(imp_cor)
    end if

    allocate (imp_cor(n_adas))  !< Dynamically allocate memeries for adas data
    if (nimp_bg(1) .gt. 0 .or. with_impurities) then
      do i=1, n_adas
        select case ( trim(imp_type(i)) )
          case('C')
            adas_suffix = '96_c'
          case('H')
            adas_suffix = '12_h'
          case('D')
            adas_suffix = '12_h' ! We are using the hydrogen data here, in want of deuterium data
          case('Ar')
            adas_suffix = '89_ar'
          case('Ne')
            adas_suffix = '96_ne'
          case('W')
            adas_suffix = '50_w'
          case('Be')
            adas_suffix = '96_be'
          case('N')
            adas_suffix = '96_n'	  
          case default
            write(*,*) "Unrecognized species, terminating."
            adas_suffix = 'none'
            deallocate(imp_cor)
            deallocate(imp_adas)
            stop
        end select

        imp_adas(i) = read_adf11(my_id, trim(adas_suffix),trim(adas_dir))
        imp_cor(i)  = coronal(imp_adas(i))

        
        ! This is to output a coronal equilibrium charge distribution as a
        ! function of temperature assuming constant density
        if (my_id == 0) call output_coronal(imp_cor(i),i)
      end do
    end if

    if (allocated(xtime_radiation)) call tr_deallocate(xtime_radiation,"xtime_radiation",CAT_GRID)
    if (nstep .gt. 0) call tr_allocate(xtime_radiation,1,nstep,"xtime_radiation")
    if (allocated(xtime_rad_power)) call tr_deallocate(xtime_rad_power,"xtime_rad_power",CAT_GRID)
    if (nstep .gt. 0) call tr_allocate(xtime_rad_power,1,nstep,"xtime_rad_power")
    if (allocated(xtime_E_ion)) call tr_deallocate(xtime_E_ion,"xtime_E_ion",CAT_GRID)
    if (nstep .gt. 0) call tr_allocate(xtime_E_ion,1,nstep,"xtime_E_ion")
    if (allocated(xtime_E_ion_power)) call tr_deallocate(xtime_E_ion_power,"xtime_E_ion_power",CAT_GRID)
    if (nstep .gt. 0) call tr_allocate(xtime_E_ion_power,1,nstep,"xtime_E_ion_power")

  end subroutine init_imp_adas

  subroutine radiation_function(ad,cor, density, temperature, Lrad, dLrad_dTe, dLrad_dNe)

    use phys_module
    use mod_openadas
    use mod_coronal
    use mod_interp_splinear

    implicit none

    type(adf11_all), intent(in) :: ad
    type(coronal), intent(in)   :: cor
    real*8, intent(in)          :: density !< log10 density in m^-3
    real*8, intent(in)          :: temperature !< log10 electron temperature in K

    real*8, intent(out)         :: Lrad ! value of radiation function
    real*8, intent(out), optional :: dLrad_dTe, dLrad_dNe ! derivatives of radiation functioni

    real*8                      :: rad!Local density multiplied radiation function
    real*8                      :: radRB, dradRB_dT, dradRB_dn, radLT, dradLT_dT, dradLT_dn
    real*8, dimension(0:ad%n_Z) :: rad_p, drad_dT, drad_dn
    real*8, dimension(0:cor%n_Z):: p          !< charge state distribution
    real*8, dimension(0:cor%n_Z):: p_Te, p_Ne !< gradient of distribution of charge states (sum = 1) to Te and Ne
    integer :: iz
    
    call cor%interp(density,temperature,rad_out=rad)
    Lrad = rad / (10.0**density) ! This is to recover the radiation coefficient
    if (present(dLrad_dTe) .or. present(dLrad_dNe)) then
      call cor%interp(density,temperature,p_out=p,p_Te_out=p_Te,p_Ne_out=p_Ne)
      do iz=0,ad%n_Z
        call ad%PRB%interp(iz,density,temperature,GRC_out=radRB,dGRC_dT_out=dradRB_dT,dGRC_dn_out=dradRB_dn)
        call ad%PLT%interp(iz,density,temperature,GRC_out=radLT,dGRC_dT_out=dradLT_dT,dGRC_dn_out=dradLT_dn)
        rad_p(iz)   = radRB + radLT
        drad_dT(iz) = dradRB_dT + dradLT_dT
        drad_dn(iz) = dradRB_dn + dradLT_dn
      enddo ! radiation emitted by atoms at level iz
      if (present(dLrad_dTe)) dLrad_dTe = dot_product(p_Te,rad_p) + dot_product(p,drad_dT)
      if (present(dLrad_dNe)) dLrad_dTe = dot_product(p_Ne,rad_p) + dot_product(p,drad_dn)
    end if

  end subroutine radiation_function

  subroutine radiation_function_linear(ad,cor, density, temperature, opt_ju, Lrad, dLrad_dT)

    use phys_module
    use mod_openadas
    use mod_coronal
    use mod_interp_splinear 

    implicit none

    type(adf11_all), intent(in) :: ad
    type(coronal), intent(in)   :: cor
    real*8, intent(in)          :: density !< log10 density in m^-3
    real*8, intent(in)          :: temperature !< log10 electron temperature in K
    logical, intent(in)         :: opt_ju !Convert outputs into jorek units if .true.

    real*8, intent(out)         :: Lrad ! value of radiation function
    real*8, intent(out), optional :: dLrad_dT ! derivative of radiation function wrt. T if with_TiTe=.f. and wrt. Te if with_TiTe=.t.

    real*8                      :: rad!Local density multiplied radiation function
    real*8                      :: radRB, radLT
    real*8, dimension(0:ad%n_Z) :: rad_p, drad_dT, dradRB_dT, dradLT_dT
    real*8, dimension(0:cor%n_Z):: p          !< charge state distribution
    real*8, dimension(0:cor%n_Z):: p_Te       !< gradient of distribution of charge states (sum = 1) to Te and Ne
    integer :: iz
    real*8                      :: coef_rad_imp ! coefficient to transform Lrad from SI to JU 
    
    call cor%interp_linear(density,temperature,rad_out=rad)
    Lrad = rad / (10.0**density) ! This is to recover the radiation coefficient
    if (present(dLrad_dT)) then
      call cor%interp_linear(density,temperature,p_out=p,p_Te_out=p_Te)
      dradRB_dT = ad%PRB%interp_grad_T(density,temperature) !Loglog gradient still!!!
      dradLT_dT = ad%PLT%interp_grad_T(density,temperature) !Loglog gradient still!!!
      do iz=0,ad%n_Z
        radRB     = ad%PRB%interp_linear(iz,density,temperature)
        radLT     = ad%PLT%interp_linear(iz,density,temperature)
        rad_p(iz)   = radRB + radLT
        drad_dT(iz) = dradRB_dT(iz) * radRB / (10.0**temperature) &
                      + dradLT_dT(iz) * radLT / (10.0**temperature) ! Convert to normal gradient
      enddo ! radiation emitted by atoms at level iz     
      if (present(dLrad_dT)) dLrad_dT = dot_product(p_Te,rad_p) + dot_product(p,drad_dT)
    end if

    !---------------------------------------------------------------------------                  
    ! --- Some post-processing to convert units, check NaNs, etc, before output
    !---------------------------------------------------------------------------  
    ! Normalization coefficient for radiation rate from SI units (W.m^3) to JOREK units:
    coef_rad_imp = (GAMMA-1.d0)*MU_ZERO**1.5d0*(central_mass*MASS_PROTON)**0.5d0 &
                       *(central_density*1.d20)**2.5d0

    if (opt_ju) then !Convert to JOREK units
      Lrad = Lrad * coef_rad_imp
      if (present(dLrad_dT)) then
        ! Convert gradient wrt. to T from 1/K into 1/eV
        dLrad_dT = dLrad_dT * coef_rad_imp *  EL_CHG / K_BOLTZ 
        ! ...and now from 1/eV into 1/(JOREK units)
        if (with_TiTe) then
          dLrad_dT = dLrad_dT / (EL_CHG*MU_ZERO*central_density*1.d20)
        else
          dLrad_dT = dLrad_dT / (2.d0*EL_CHG*MU_ZERO*central_density*1.d20)
        endif
      end if
    end if
    if (Lrad < 0.) then
      Lrad = 0.
      if (present(dLrad_dT)) dLrad_dT = 0.
    end if  
    if (Lrad/=Lrad) then
      write(*,*) "WARNING: Lrad ", Lrad
      stop
    end if
    if (present(dLrad_dT) .and. dLrad_dT/=dLrad_dT) then
      write(*,*) "WARNING: dLrad_dT ", dLrad_dT
      stop
    end if
  end subroutine radiation_function_linear

  !> This is to output a coronal equilibrium charge distribution as a
  !> function of temperature assuming constant density 10^20/m^3
  !> to a file charge_distribution.dat
  !> plot with gnuplot like
  !> 
  !> set logscale y
  !> p for [i=2:20] 'charge_distribution.dat' u 1:i t ''.(i-2) w l
  subroutine output_coronal(cor,i_imp)
  use phys_module, only: imp_type
  use mod_coronal
  class(coronal), intent(in)      :: cor !< Coronal equilibrium type
  integer, intent(in)             :: i_imp !< Coronal equilibrium type
  
  ! Temporary variable for charge state distribution
  integer             :: i_T, i_ion
  real*8, allocatable :: P_imp(:)
  real*8              :: Te_eV, Lrad
  real*8              :: Z_eff
  character(len=256)  :: filename
  
  
  write (filename,'(a,a)') trim(imp_type(i_imp)), "_charge_distribution.dat"
  open(20,file=trim(filename))
  
  write(20,'(4A22)',advance='no') 'temperature (log10(K))', 'charge states', 'summation', 'effective charge'
  write(20,'(A22)') 'radiation function'
  
  do i_T = 1, size(cor%temperature,1)
    Te_eV = cor%temperature(i_T)
  
    if (allocated(P_imp)) deallocate(P_imp)
  
    allocate(P_imp(0:cor%n_Z))
    call cor%interp(density=20.d0,temperature=Te_eV,p_out=P_imp,z_out=Z_eff,rad_out=Lrad)
    Lrad = Lrad / (1.d20) ! This is to recover the radiation coefficient
    write(20,'(f12.3)',advance='no') Te_eV
    do i_ion = 0, cor%n_Z
      write(20,'(f12.5)',advance='no') P_imp(i_ion)
    end do
    write(20,'(f12.5)',advance='no') sum(P_imp)
    write(20,'(f12.5)',advance='no') Z_eff
    write(20,'(e14.6)') Lrad
  end do
  close (20)

end subroutine output_coronal
end module mod_impurity

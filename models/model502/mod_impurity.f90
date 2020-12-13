module mod_impurity

  use constants

  contains 

  subroutine init_imp_adas(my_id)

    use phys_module
    use mod_openadas
    use mod_coronal

    implicit none

    integer, intent(in) :: my_id
    integer             :: err_alloc, i

    character(len=512)  :: adas_suffix     !The suffix of adas data file to be read

    ! Temporary variable for charge state distribution
    integer             :: i_T, i_ion
    real*8, allocatable :: dP_imp_dT(:), P_imp(:)
    real*8              :: Z_imp

    n_adas = 1 ! For now we only trace one species, in the future probably more 

    if (allocated(imp_adas)) then
      deallocate(imp_adas)
    end if

    allocate (imp_adas(n_adas),stat=err_alloc)  !< Dynamically allocate memeries for adas data

    if (err_alloc /= 0) then
      write(*,*) "Error when trying to dynamically allocate memeries for adas data.", my_id
      stop
    else
      if (allocated(imp_cor)) then
        deallocate(imp_cor)
      end if

      allocate (imp_cor(n_adas),stat=err_alloc)  !< Dynamically allocate memeries for adas data
      if (err_alloc /= 0) then
        write(*,*) "Error when trying to dynamically allocate memeries for CE vector.", my_id
        deallocate(imp_adas)
        stop
      else
        do i=1, n_adas
          select case ( trim(gas_type) )
            case('D2')
              write(*,*) "Deuterium adas calculation unsupported for now, terminating."
              adas_suffix = 'none'
              deallocate(imp_cor)
              deallocate(imp_adas)
              stop
            case('Ar')
              adas_suffix = '89_ar'
            case('Ne')
              adas_suffix = '96_ne'
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
          if (my_id == 0) call output_coronal(imp_cor(i))
        end do
      end if

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

  subroutine radiation_function_linear(ad,cor, density, temperature, Lrad, dLrad_dTe)

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
    real*8, intent(out), optional :: dLrad_dTe ! derivatives of radiation functioni

    real*8                      :: rad!Local density multiplied radiation function
    real*8                      :: radRB, radLT
    real*8, dimension(0:ad%n_Z) :: rad_p, drad_dT, dradRB_dT, dradLT_dT
    real*8, dimension(0:cor%n_Z):: p          !< charge state distribution
    real*8, dimension(0:cor%n_Z):: p_Te       !< gradient of distribution of charge states (sum = 1) to Te and Ne
    integer :: iz
    
    call cor%interp_linear(density,temperature,rad_out=rad)
    Lrad = rad / (10.0**density) ! This is to recover the radiation coefficient
    if (present(dLrad_dTe)) then
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
      if (present(dLrad_dTe)) dLrad_dTe = dot_product(p_Te,rad_p) + dot_product(p,drad_dT)
    end if

  end subroutine radiation_function_linear

end module mod_impurity

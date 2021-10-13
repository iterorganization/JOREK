!**********************************************************************
!* program to extract data from a JOREK2 restart file                 *
!*                                                                    *     
!* This diagnostic has to be run with MPI (mpirun/mpiexec/srun        *
!* depending on the system). For details, see:                        *
!* https://www.jorek.eu/wiki/doku.php?id=diagnostics#diagnostics      *
!**********************************************************************

program jorek2_diagno_spi
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
use data_structure
use phys_module
use basis_at_gaussian
use pellet_module
use mpi_mod
use mod_boundary, only: boundary_from_grid 
use mod_import_restart
use mod_impurity
use mod_integrals3D
use mod_expression, only: exprs_all_int, init_expr, t_expr_list

implicit none

type (type_node_list)        :: node_list
type (type_element_list)     :: element_list
type (type_bnd_node_list)    :: bnd_node_list
type (type_bnd_element_list) :: bnd_elm_list
type (t_expr_list)           :: expr_list
real*8, allocatable          :: res(:)
integer                      :: units


integer :: i, in, i_tor, i_spi
real*8  :: growth_kin, growth_mag,density,density_in,density_out,pressure,pressure_in,pressure_out
real*8  :: Rplot(2), Zplot(2)
real*8  :: psi_axis,R_axis,Z_axis,s_axis,t_axis
integer :: ifail, ierr, i_elm_axis
integer :: required, provided, StatInfo
real*8  :: spi_abl_rate_tot, spi_abl_tot
real*8  :: spi_abl_bg_rate_tot, spi_abl_bg_tot
real*8  :: spi_x, spi_y


write(*,*) '***************************************'
write(*,*) '* JOREK2_diagno                       *'
write(*,*) '***************************************'

call init_expr()
allocate(res(exprs_all_int%n_expr+1))
res = 0.d0


#ifdef FUNNELED
  required = MPI_THREAD_FUNNELED
#else
  required = MPI_THREAD_MULTIPLE
#endif
call MPI_Init_thread(required, provided, StatInfo)


call initialise_parameters(0, "__NO_FILENAME__")

do i_tor=1, n_tor
  mode(i_tor) = + int(i_tor / 2) * n_period
  write(*,*) ' toroidal mode numbers : ',i_tor,mode(i_tor)
enddo

call import_restart(node_list,element_list, 'jorek_restart', rst_format, ierr, .true.)

call initialise_basis()                              ! define the basis functions at the Gaussian points

open(20,file='energies.txt')

write(20,'(A,25(A11,i3.3))') '      i      time',('          M',n_period*((in-1)/2),in=1,n_tor,2), &
                                                 ('          K',n_period*((in-1)/2),in=1,n_tor,2)

do i=2,index_start

 !Growth_mag  = 0.5d0*log(abs(energies(n_tor,1,i)/energies(n_tor,1,i-1))) &
 !            / (xtime(i)-xtime(i-1))
 !Growth_kin  = 0.5d0*log(abs(energies(n_tor,2,i)/energies(n_tor,2,i-1))) &
 !            / (xtime(i)-xtime(i-1))

 !write(*,'(i7,f12.3,200e14.6)') i,xtime(i),energies(1:n_tor,:,i),growth_mag,growth_kin

 write(20,'(i7,f12.3,200e14.6)') i,xtime(i),energies(1,1,i),(energies(in,1,i)+energies(in+1,1,i),in=2,n_tor,2), &
                                            energies(1,2,i),(energies(in,2,i)+energies(in+1,2,i),in=2,n_tor,2)

enddo
close(20)

if (use_pellet) then

  open(20,file="pellet.txt")

  write(20,'(A,25(A11,i3.3))') '      i,      time,    pellet_R,    pellet_Z,   pellet_psi,  particles,   ablation'

  do i=1,index_start
    write(20,'(i7,f12.3,200e14.6)') i,xtime(i),xtime_pellet_R(i),xtime_pellet_Z(i),xtime_pellet_psi(i),xtime_pellet_particles(i), &
                                    xtime_phys_ablation(i)
  enddo
  close(20)

endif

open(20,file="rad_history.dat")

write(20,'(3A20)') 'time', 'tot_rad_power (MW)', 'total_radiation (MJ)'

do i=1,index_start
  write(20,'(i7,f12.3,2e14.6)') i,xtime(i), xtime_rad_power(i)/1.d6, xtime_radiation(i)/1.d6
enddo
close(20)

open(20,file="thermal_history.dat")

write(20,'(3A20)') 'time', 'e_th_energy (MJ)', 'i_th_energy (MJ)'

do i=1,index_start
  write(20,'(i7,f12.3,2e14.6)') i,xtime(i), thermal_e_tot_t(i)/1.d6, thermal_i_tot_t(i)/1.d6
enddo
close(20)

if (using_spi) then

  open(20,file="abl_history.dat")

  write(20,'(A11)') 'time', 'total_abl_rate', 'total_abl_number'

  do i=1,index_start
    spi_abl_rate_tot = 0.0
    spi_abl_tot = 0.0
    spi_abl_bg_rate_tot = 0.0
    spi_abl_bg_tot = 0.0
    do i_spi = 1, n_spi_tot
      spi_abl_rate_tot = spi_abl_rate_tot + xtime_spi_ablation_rate(i_spi,i)
      spi_abl_tot = spi_abl_tot + xtime_spi_ablation(i_spi,i)
      spi_abl_bg_rate_tot = spi_abl_bg_rate_tot + xtime_spi_ablation_bg_rate(i_spi,i)
      spi_abl_bg_tot = spi_abl_bg_tot + xtime_spi_ablation_bg(i_spi,i)
    end do
    write(20,'(i7,f12.3,4e14.6)') i,xtime(i), spi_abl_rate_tot, spi_abl_tot, spi_abl_bg_rate_tot, spi_abl_bg_tot
  enddo
  close(20)

  open(20,file="fragments_position.dat")

  do i_spi = 1, n_spi_tot
    write(20,'(i7,3f12.3,e14.6,f12.3)') i_spi, pellets(i_spi)%spi_R, pellets(i_spi)%spi_Z, pellets(i_spi)%spi_phi,&
                                               pellets(i_spi)%spi_radius, pellets(i_spi)%spi_species
  end do
  close(20)

  ! SPI fragments' position in (x,y,z) coordinate
  open(20,file="fragments_position_xyz.dat")

  do i_spi = 1, n_spi_tot
    spi_x  =   pellets(i_spi)%spi_R * cos(pellets(i_spi)%spi_phi)
    spi_y  = - pellets(i_spi)%spi_R * sin(pellets(i_spi)%spi_phi)
    write(20,'(i7,3e14.6,e14.6,e14.6)') i_spi, spi_x, spi_y, pellets(i_spi)%spi_Z,&
                                               pellets(i_spi)%spi_radius, pellets(i_spi)%spi_species
  end do
  close(20)


endif

#if (defined WITH_Neutrals) || (defined WITH_Impurities)
  ! --- Read ADAS data and generate coronal equilibrium is needed
  call init_imp_adas(0)
  if (output_prad_phi) then
    ! --- Determine boundary information from the grid
    call boundary_from_grid(node_list, element_list, bnd_node_list, bnd_elm_list, .false.)

    call int3d_new(0, node_list, element_list, bnd_node_list, bnd_elm_list, exprs_all_int, res, 1)
  endif
#endif
!if (use_pellet) then
!   pellet_volume = total_pellet_volume
!   call update_pellet(0,node_list,element_list)
!end if
!------------------lowshape3bis outside
!Rplot(1) = 3.0
!Rplot(2) = 3.676
!Zplot(1) = -2.066
!Zplot(2) = -1.9265

!------------------ lowshape3bis inside
!Rplot(1) = 2.3213
!Rplot(2) = 2.8511
!Zplot(1) = -1.9178
!Zplot(2) = -2.0339

!------------------ lowshape3bis midplane
!Rplot(1) = 1.88
!Rplot(2) = 4.88
!Zplot(1) = 0.09
!Zplot(2) = 0.09

!-----------------lowshape7,8 (midplane)
!Rplot(1) = 1.9
!Rplot(2) = 4.2
!Zplot(1) = 0.07
!Zplot(2) = 0.07

!call find_axis(0,node_list,element_list,psi_axis,R_axis,Z_axis,i_elm_axis,s_axis,t_axis,ifail)

!Rplot(1) = 1.0
!Rplot(2) = 3.5
!Zplot(1) = Z_axis
!Zplot(2) = Z_axis 

!call plot_profiles(node_list,element_list,Rplot,Zplot)


!call export_helena(node_list,element_list)

!----------------------------------------- plot profiles
!call begplt('profiles.ps')

!call plot_velocity_profile(node_list,element_list, 3.d0, 0.d0, 3.d0, 2.d0)

!call finplt

end program jorek2_diagno_spi

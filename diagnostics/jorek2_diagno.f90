!**********************************************************************
!* program to extract data from a JOREK2 restart file                 *
!*                                                                    *     
!* This diagnostic has to be run with MPI (mpirun/mpiexec/srun        *
!* depending on the system). For details, see:                        *
!* https://www.jorek.eu/wiki/doku.php?id=diagnostics#diagnostics      *
!**********************************************************************

program jorek2_diagno
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
use data_structure
use phys_module
use basis_at_gaussian
use pellet_module
use mpi_mod
use mod_import_restart
use mod_impurity

implicit none

type (type_node_list)    :: node_list
type (type_element_list) :: element_list
integer :: i, in, i_tor
real*8  :: growth_kin, growth_mag,density,density_in,density_out,pressure,pressure_in,pressure_out
real*8  :: Rplot(2), Zplot(2)
real*8  :: psi_axis,R_axis,Z_axis,s_axis,t_axis
integer :: ifail, my_id, ierr, i_elm_axis
integer :: required, provided, StatInfo


write(*,*) '***************************************'
write(*,*) '* JOREK2_diagno                       *'
write(*,*) '***************************************'

my_id=0

#ifdef FUNNELED
  required = MPI_THREAD_FUNNELED
#else
  required = MPI_THREAD_MULTIPLE
#endif
call MPI_Init_thread(required, provided, StatInfo)


call initialise_parameters(my_id, "__NO_FILENAME__")

do i_tor=1, n_tor
  mode(i_tor) = + int(i_tor / 2) * n_period
  write(*,*) ' toroidal mode numbers : ',i_tor,mode(i_tor)
enddo

call import_restart(node_list,element_list, 'jorek_restart', rst_format, ierr, .true.)

call initialise_basis                              ! define the basis functions at the Gaussian points

open(20,file='energies.txt')

write(20,'(A,25(A11,i3.3))') '      i      time',('          M',n_period*((in-1)/2),in=1,n_tor,2), &
                                                 ('          K',n_period*((in-1)/2),in=1,n_tor,2)

do i=2,index_start

 Growth_mag  = 0.5d0*log(abs(energies(n_tor,1,i)/energies(n_tor,1,i-1))) &
             / (xtime(i)-xtime(i-1))
 Growth_kin  = 0.5d0*log(abs(energies(n_tor,2,i)/energies(n_tor,2,i-1))) &
             / (xtime(i)-xtime(i-1))

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

#if (defined WITH_Neutrals) || (defined WITH_Impurities)
  ! --- Read ADAS data and generate coronal equilibrium if needed
  call init_imp_adas(my_id)
#endif


call Integrals_3D(my_id,node_list,element_list,density,density_in,density_out,pressure,pressure_in,pressure_out)

if (use_pellet) then
   pellet_volume = total_pellet_volume
   call update_pellet(my_id,node_list,element_list)
end if
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

!call find_axis(my_id,node_list,element_list,psi_axis,R_axis,Z_axis,i_elm_axis,s_axis,t_axis,ifail)

Rplot(1) = 1.0
Rplot(2) = 3.5
Zplot(1) = Z_axis
Zplot(2) = Z_axis 

call plot_profiles(node_list,element_list,Rplot,Zplot)


!call export_helena(node_list,element_list)

!----------------------------------------- plot profiles
!call begplt('profiles.ps')

!call plot_velocity_profile(node_list,element_list, 3.d0, 0.d0, 3.d0, 2.d0)

!call finplt

end program jorek2_diagno

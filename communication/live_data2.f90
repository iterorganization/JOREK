!> The module is a duplicate of live_data, in order to output additional data
!!
!! - The file <b>mac_vars2.dat</b> is created during the code run and filled with
!!   information about certain run parameters, energy timetraces, growth rates, etc.
!! - The input parameter phys_module::produce_live_data allows to switch the functionality of
!!   this module on (default) or off.
!! - The script extract_live_data.sh in the util/ folder allows to extract certain live data from
!!   the output file. Run it with option -h for usage information.
!! - The script plot_live_data.sh allows to plot live data, e.g., the energy
!!   time traces. Run it with option -h for usage information.
!!
module live_data2
  
#include "version.h"
  
  implicit none
  
  private
  public init_live_data2, write_live_data2, finalize_live_data2
  
  integer,           parameter :: LIVE_DATA_HANDLE2 = 66 !< File handle for live data file
  character(len=20), parameter :: LIVE_DATA_FILE2   = 'mac_vars2.dat' !< Live data file
  
  
  
  contains
  
  
  
  !> Open file, write out headers and some parameters.
  subroutine init_live_data2()
    
    use mod_parameters,    only: n_tor, n_plane, n_period, jorek_model, variable_names, n_var
    use phys_module,   only: produce_live_data, mode, mode_type, xpoint, xcase
    
    implicit none
    
    logical :: opened
    integer :: n, i
    
    if ( .not. produce_live_data ) return
    
    ! --- Check, that the file handle is not already in use.
    inquire(unit=LIVE_DATA_HANDLE2, opened=opened)
    if ( opened ) then
      write(*,*) 'WARNING: LIVE DATA CANNOT BE PRODUCED AS FILE HANDLE IS ALREADY IN USE!'
      produce_live_data = .false.
      return
    end if
    
    open(LIVE_DATA_HANDLE2, file=LIVE_DATA_FILE2, status='REPLACE', action='WRITE')
    
    ! --- Write some general information
    write(LIVE_DATA_HANDLE2,*)  '@rcs_version: ', RCS_VERSION
    write(LIVE_DATA_HANDLE2,'(A,I5)') '@jorek_model: ', jorek_model
    write(LIVE_DATA_HANDLE2,'(A,I5)') '@n_tor: ', n_tor
    write(LIVE_DATA_HANDLE2,'(A,I5)') '@n_plane: ', n_plane
    write(LIVE_DATA_HANDLE2,'(A,I5)') '@n_period: ', n_period
    write(LIVE_DATA_HANDLE2,'(A)') '@plottable: energies growth_rates times input_profiles'
    write(LIVE_DATA_HANDLE2,'(A,15(A11,1X))') '@variable_names: ', variable_names((/(i, i=1,n_var)/))
    
    ! --- Write file headers indicating what data is in the files.
    write(LIVE_DATA_HANDLE2,'(A,I5)') '@n_times: ', 1
    write(LIVE_DATA_HANDLE2,'(A)') '@times_xlabel: time step'
    write(LIVE_DATA_HANDLE2,'(A)') '@times_ylabel: normalized time'
    write(LIVE_DATA_HANDLE2,'(A)') '@times_logy: 0'
    write(LIVE_DATA_HANDLE2,'(A)') '@times: "step"     "time"'
    
    
    write(LIVE_DATA_HANDLE2,'(A,I5)') '@n_energies: ', 2*(n_tor+1)/2
    write(LIVE_DATA_HANDLE2,'(A)') '@energies_xlabel: normalized time'
    write(LIVE_DATA_HANDLE2,'(A)') '@energies_ylabel: normalized energy'
    write(LIVE_DATA_HANDLE2,'(A)') '@energies_logy: 1'
    write(LIVE_DATA_HANDLE2,'(A)',advance='no') '@energies: %"time"           '
    do n = 1, n_tor, 2
      write(LIVE_DATA_HANDLE2,'(A7,",",I2.2,A2,1x)',advance='no') '"A_{tem', mode(n), '}"'
    end do
    do n = 1, n_tor, 2
      write(LIVE_DATA_HANDLE2,'(A7,",",I2.2,A2,1x)',advance='no') '"A_{den', mode(n), '}"'
    end do
    write(LIVE_DATA_HANDLE2,*)
    
    write(LIVE_DATA_HANDLE2,'(A,I5)') '@n_growth_rates: ', 2*(n_tor+1)/2
    write(LIVE_DATA_HANDLE2,'(A)') '@growth_rates_xlabel: normalized time'
    write(LIVE_DATA_HANDLE2,'(A)') '@growth_rates_ylabel: normalized growth rate'
    write(LIVE_DATA_HANDLE2,'(A)') '@growth_rates_logy: 1'
    write(LIVE_DATA_HANDLE2,'(A)',advance='no') '@growth_rates: %"time"           '
    do n = 1, n_tor, 2
      write(LIVE_DATA_HANDLE2,'(A7,",",I2.2,A2,1x)',advance='no') '"G_{mag', mode(n), '}"'
    end do
    do n = 1, n_tor, 2
      write(LIVE_DATA_HANDLE2,'(A7,",",I2.2,A2,1x)',advance='no') '"G_{kin', mode(n), '}"'
    end do
    write(LIVE_DATA_HANDLE2,*)
    
    ! --- Call the model-specific part of the init_live_data routine
    call init_live_data_model(LIVE_DATA_HANDLE2) 
    
    close(LIVE_DATA_HANDLE2)
    
  end subroutine init_live_data2
  
  
  
  !> Write out data to text files during the code run.
  subroutine write_live_data2(index)
    
    use mod_parameters,  only: n_tor
    use phys_module, only: xtime, energies2, produce_live_data
    
    implicit none
    
    integer, intent(in) :: index !< Timestep index to write data for
    
    integer :: i, j
    real*8  :: e1, e2, growth_rate
    
    if ( .not. produce_live_data ) return
    
    open(LIVE_DATA_HANDLE2, file=LIVE_DATA_FILE2, status='OLD', position='APPEND', action='WRITE')
    
    ! --- Write data to the files.
    write(LIVE_DATA_HANDLE2,'(A,I6,1X,ES17.9)') '@times:', index, xtime(index)
    write(LIVE_DATA_HANDLE2,'(A,ES17.9)',advance='no') '@energies:', xtime(index)
    do j = 1, 2
      do i = 1, n_tor, 2
        write(LIVE_DATA_HANDLE2,'(ES17.9)',advance='no') sum(energies2(max(i-1,1):i,j,index))
      end do
    end do
    write(LIVE_DATA_HANDLE2,*)
    if ( index > 1 ) then
      write(LIVE_DATA_HANDLE2,'(A,ES17.9)',advance='no') '@growth_rates:', &
        (xtime(index)+xtime(index-1))/2.d0
      do j = 1, 2
        do i = 1, n_tor, 2
          e1 = sum(energies2(max(i-1,1):i,j,index))
          e2 = sum(energies2(max(i-1,1):i,j,index-1))
          growth_rate = 0.5d0 * ( log(e1) - log(e2) ) / (xtime(index)-xtime(index-1))
          write(LIVE_DATA_HANDLE2,'(ES17.9)',advance='no') growth_rate
        end do
      end do
    end if
    
    close(LIVE_DATA_HANDLE2)
    
  end subroutine write_live_data2
  
  
  
  !> Close file.
  subroutine finalize_live_data2()
    
    use phys_module, only: produce_live_data
    
    implicit none
    
    if ( .not. produce_live_data ) return
    
    ! -nothing to be done currently-
    
  end subroutine finalize_live_data2

end module live_data2

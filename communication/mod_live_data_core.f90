!> Simple dispatcher routine to run all appropriate live_data subroutines.
!> Note that this module does not care for MPI, i.e. be careful yourself
!> of running this only from task 0.
module mod_live_data_core
  use live_data
#ifdef JECCD
  use live_data2
  use live_data3
#ifdef JEC2DIAG
  use live_data4
#endif
#endif
  
  implicit none
  
contains
  
  !> Initialize various levels of live data reporting based on defines.
  subroutine init_live_data_all()
    call init_live_data()
#ifdef JECCD
    call init_live_data2()
    call init_live_data3()
#ifdef JEC2DIAG
    call init_live_data4()
#endif
#endif
  end subroutine init_live_data_all
  
  
  subroutine write_live_data_all(index_now)
    integer, intent(in) :: index_now
    call write_live_data(index_now)
#ifdef JECCD
    call write_live_data2(index_now)
    call write_live_data3(index_now)
#ifdef JEC2DIAG
    call write_live_data4(index_now)
#endif
#endif
  end subroutine write_live_data_all

  subroutine finalize_live_data_all()
    call finalize_live_data()
#ifdef JECCD
    call finalize_live_data2()
    call finalize_live_data3()
#ifdef JEC2DIAG
    call finalize_live_data4()
#endif
#endif
  end subroutine finalize_live_data_all
end module mod_live_data_core

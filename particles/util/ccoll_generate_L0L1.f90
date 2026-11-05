!< Simple program that generates a file "ccoll.data" that contains tabulated L0 and L1 values
!< that should be applicable for every fusion plasma.

program ccoll_generate_L0L1
  use mod_ccoll_relativistic
  
  implicit none

  type(ccoll_data) :: dat
  logical :: storage_file_on_disk
  character(20), parameter :: storage_file='ccolldata'

   real*8, parameter  :: uminxp  = -4.D0
   real*8, parameter  :: umaxxp  = 2.D0
   real*8, parameter  :: thminxp = -5.D0
   real*8, parameter  :: thmaxxp = -1.D0 
   integer, parameter :: nu      = 1001
   integer, parameter :: nth     = 801

   ! Initialize the look-up tables
  print*,''
  write(*,*) 'Initializing look-up tables...'
  inquire(file=storage_file,exist=storage_file_on_disk)
  if(storage_file_on_disk) then
     write(*,*) 'File already exists.'
  else
     dat = ccoll_compute_L0L1table(uminxp,umaxxp,thminxp,thmaxxp,nu,nth)
     call ccoll_write_L0L1table(dat,storage_file)
  end if
  print*,'Done!'

end program ccoll_generate_L0L1

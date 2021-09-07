!> module takes the OPEN-ADAS data to calculate the steady state (or time evolution) charge distribution
!> and average charge state as a function of temperature
module mod_openadas
use mod_interp_splinear
implicit none
private
public ADF11, ADF11_all, read_ADF11

!> Custom data structure containing relevant fields from ADF11 format files (unresolved case!)
type ADF11
  integer             :: n_Z !< Atomic number
  integer             :: izmin, izmax !< minimum and maximum value of z for which data is available
  real*8, allocatable :: density(:) !< log10 density (m^-3)
  real*8, allocatable :: temperature(:) !< log10 temperature (K)
  real*8, allocatable :: GRC(:,:,:) !< log10 of coefficient (parameters: d, T, z). Units:
  !< ACD, SCD: m3s-1 (for *CD ?) (converted from cm3s-1)
  !< PLT, PRB: Wm3 (for P* ?) (converted from Wcm3)
  !< z goes from 0 to n_Z inclusive
  type(Fspline), allocatable :: GRCFspline (:) !< spline functions for coefficients at each charge state
contains
  procedure :: interp => GRC_spl
  procedure :: interp_grad_T => dGRC_dT
  procedure :: interp_grad_n => dGRC_dn
  procedure :: interp_linear => GRC
end type ADF11

!> Compound datatype containing many type_ADF11
type ADF11_all
  integer     :: n_Z !< Atomic number
  type(ADF11) :: ACD !< Effective recombination coefficients
  type(ADF11) :: SCD !< Effective ionisation coefficients
  type(ADF11) :: CCD !< Charge exchange effective recombination coefficients
  type(ADF11) :: PLT !< Line power driven by excitation of dominant ions
  type(ADF11) :: PRB !< Continuum and line power driven by recombination and bremsstrahlung of dominant ions
  type(ADF11) :: PRC !< Line power due to charge transfer from thermal neutral hydrogen to dominant ions
  real*8, dimension(:), allocatable :: ionisation_energy !< energy in eV required to ionize to a level, indexed by the new charge state (i.e. 1 to 74 for W), no interpolation needed
end type ADF11_all
contains

!> Read ADF11 data files and import them into a type_ADF11
!> Tries to read ACD, SCD, CCD, PLT, PRB, PRC coefficients
!> if the files exist. Files of format acd$suffix.dat are read.
!> Suffix is usually of the form 50_w, 96_li
!> Try to also read the ionisation energy coefficients, but don't crash if they
!> are not present.
function read_adf11(my_id,suffix, directory) result(ad)
use constants
character(len=*), intent(in) :: suffix !< Usually year_atom (ex: 50_w, 96_li)
character(len=*), intent(in), optional :: directory
type(ADF11_all), target :: ad !< OpenAdas data type

type(ADF11), pointer :: a
integer :: i_ADF11
character*3, dimension(1:6), parameter :: ADF11_filenames = (/"acd", "scd", "ccd", "plt", "prb", "prc"/)
character*120 :: filename

integer, intent(in) :: my_id
integer :: i, ierr, n_d, n_T, k, q, i_n
logical :: file_exists, recombining

if (my_id .eq. 0) then
  write(*,'(A)') '*********************************'
  write(*,'(A)') '* Importing OpenAdas data       *'
  write(*,'(A,A,A)') '* open files ending in: ', suffix, '  *'
  if (present(directory)) write(*,'(A,A)') '* present in directory ', directory
  write(*,'(A)') '*********************************'
endif

do i_ADF11 = 1,size(ADF11_filenames,1)
  write(filename,"(A,A,A)") ADF11_filenames(i_ADF11), trim(suffix), '.dat'
  if (present(directory)) filename = trim(directory) // trim(filename)
  inquire(file=trim(filename), exist=file_exists)
  if (.not. file_exists) then
    write(*,*) "File not found for", filename
    cycle ! Skip this type of data
  end if

  if (my_id .eq. 0) write(*,"(A,A)",advance="no") "Reading data from ", trim(filename)
  open(10,file=trim(filename),action="read",status="old",iostat=ierr)
  if (ierr .ne. 0) then
    write(*,*) my_id, " failed with code ", ierr
    cycle
  endif

  ! Point a to right variable to read in
  select case (i_ADF11)
    case (1); a => ad%ACD; recombining=.true.
    case (2); a => ad%SCD; recombining=.false.
    case (3); a => ad%CCD; write(*,*) "Warning: CCD not implemented correctly yet"; recombining=.true.
    case (4); a => ad%PLT; recombining=.false.
    case (5); a => ad%PRB; recombining=.true.
    case (6); a => ad%PRC; write(*,*) "Warning: PRC not implemented correctly yet"; recombining=.true. ! see coronal model
  end select

  read(10,*)  a%n_z, n_d, n_T, a%izmin, a%izmax
  ad%n_z = a%n_z
  allocate(a%density(n_d), a%temperature(n_T), a%GRC(n_d,n_T,0:a%n_z))
  allocate(a%GRCFspline(0:a%n_z))
  
  read(10,*)
  read(10,*) a%density(:)
  ! Convert densities to log10 of m^-3 instead of cm^-3
  a%density = a%density + 6.d0
  read(10,*) a%temperature(:)
  ! Convert temperatures to log10 of K instead of eV
  ! From E eV = kB T
  a%temperature = a%temperature - log10(K_BOLTZ) + log10(EL_CHG) ! EL_CHG * 1 Volt actually

  a%GRC = -60.d0
  ! Note that recombination data is given as recombining FROM (Z=1 to Z=74)
  ! and ionisation data is given as ionising TO (Z=1 up to Z=74)
  ! renormalise this so the recombination data is given for the source z
  ! and the ionisation data is given for the source z also
  if (recombining) then
    do i = a%izmin, a%izmax ! 1 to n_z
      read(10,*)
      read(10,*) a%GRC(:,:,i)
    end do
  else
    do i = a%izmin-1, a%izmax-1 ! 0 to n-z - 1
      read(10,*)
      read(10,*) a%GRC(:,:,i)
    end do
  end if
  close(10)

  ! Convert GRC coefficients from cm to m
  a%GRC = a%GRC - 6.d0 ! because it is a logarithm. Conversion: /100.d0**3 (cm3s-1 => m3s-1)

  ! Allocate and construct the splines
  do i = 0, a%n_z
    call AllocFspline(a%GRCFspline(i),n_T,n_d)

    a%GRCFspline(i)%xspline = a%temperature
    a%GRCFspline(i)%ylinear = a%density

    call ConstructFspline(a%GRCFspline(i),a%GRC(:,:,i))
  end do   

  if (my_id .eq. 0) write(*,"(A)") " succeeded"
enddo

! Test if ACD and SCD were loaded at least
if (.not. (allocated(ad%ACD%density) .and. allocated(ad%SCD%density))) then
  write(*,*) my_id, "ACD and SCD not found, exiting"
  call exit(10)
else
  if (my_id .eq. 0) write(*,*) 'Done reading adas data for atomic number', ad%n_Z
endif


! Try to load the ionisation energies
! try 2 cases, first the full suffix and then the stripped suffix
! i.e. ion50_w.dat and ion_w.dat
do i=1,3,2 ! full, strip
  ! assume the suffix starts with 2 digits
  write(filename,"(A,A,A)") 'ion', trim(suffix(i:len(suffix))), '.dat'
  if (present(directory)) filename = trim(directory) // trim(filename)
  inquire(file=trim(filename), exist=file_exists)
  if (file_exists) then
    open(10,file=trim(filename),status="old",iostat=ierr, action="read")
    if (ierr .ne. 0) then
      write(*,*) my_id, " failed with code ", ierr
    else
      allocate(ad%ionisation_energy(1:ad%n_Z))
      do k=1,ad%n_Z
        read(10,*) q, ad%ionisation_energy(k)
        if (q + 1 .ne. k) then ! conversion from 0-based to 1-based indices for ionisation energies
          write(*,*) 'Mismatch in detected energy levels, ', q+1, k
          stop 1
        end if
      end do
      write(*,*) "Read ionisation energies from ", trim(filename)
      close(10)
      exit ! the loop, we have found a file
    endif
  else
    if (i .eq. 3) then
      write(*,*) "Cannot find ionisation data file ", trim(filename), "not loading ionisation energies"
    end if
  end if
end do
end function read_adf11



!> interpolation of log10 temperature gradient of GRC in density and temperature
!> Loglog gradient!
function dGRC_dT(a, density, temperature)
class(ADF11), intent(in) :: a           !< ADF11 datatype
real*8, intent(in)       :: density     !< log10 density in m^-3
real*8, intent(in)       :: temperature !< log10 temperature in K
real*8, dimension(0:a%n_Z) :: dGRC_dT !< Generalized Radiational Coefficient at this density and temperature
integer                  :: i_z     !< Index of charge state

! If GRC exists and we are looking for a Z that is nonzero
if (allocated(a%GRC)) then
  dGRC_dT = L2D2interp_grad(a%density,a%temperature,a%n_Z+1,a%GRC(:,:,0:a%n_Z),density,temperature,1)
!  do i_z = 0, a%n_z
!    call SL2Dinterp(a%GRCFspline(i_z),temperature,density,dfout_dx=dGRC_dT(i_z))
!  end do
else
  dGRC_dT = 0.d0
endif
end function dGRC_dT

!> interpolation of log10 density gradient of GRC in density and temperature.
!> Loglog gradient!
function dGRC_dn(a, density, temperature)
class(ADF11), intent(in) :: a           !< ADF11 datatype
real*8, intent(in)       :: density     !< log10 density in m^-3
real*8, intent(in)       :: temperature !< log10 temperature in K
real*8, dimension(0:a%n_Z) :: dGRC_dn !< Generalized Radiational Coefficient at this density and temperature
integer                  :: i_z     !< Index of charge state

! If GRC exists and we are looking for a Z that is nonzero
if (allocated(a%GRC)) then
  dGRC_dn = L2D2interp_grad(a%density,a%temperature,a%n_Z+1,a%GRC(:,:,0:a%n_Z),density,temperature,2)
!  do i_z = 0, a%n_z
!    call SL2Dinterp(a%GRCFspline(i_z),temperature,density,dfout_dy=dGRC_dn(i_z))
!  end do
else
  dGRC_dn = 0.d0
endif
end function dGRC_dn

!> interpolation of log10 values of GRC in density and temperature
function GRC(a, z, density, temperature)
class(ADF11), intent(in) :: a           !< ADF11 datatype
real*8, intent(in)            :: density     !< log10 density in m^-3
real*8, intent(in)            :: temperature !< log10 temperature in K
integer, intent(in)           :: z !< index in a%GRC(:,:,z) (is ionisation level or ionisation level - 1, 1:n_z)
real*8 :: GRC !< Generalized Radiational Coefficient at this density and temperature

! If GRC exists and we are looking for a Z that is nonzero
if (allocated(a%GRC) .and. z .le. ubound(a%GRC,3) .and. z .ge. lbound(a%GRC,3)) then
  GRC = 10.d0**L2Dinterp(a%density,a%temperature,a%GRC(:,:,z),density,temperature)
  !call SL2Dinterp(a%GRCFspline(z),temperature,density,fout=GRC)
  !GRC = 10.d0**GRC
else
  GRC = 0.d0
endif
end function GRC

!> interpolation of log10 values of GRC in density and temperature
subroutine GRC_spl(a, z, density, temperature, GRC_out, dGRC_dT_out, dGRC_dn_out)
class(ADF11), intent(in)      :: a           !< ADF11 datatype
real*8, intent(in)            :: density     !< log10 density in m^-3
real*8, intent(in)            :: temperature !< log10 temperature in K
integer, intent(in)           :: z !< index in a%GRC(:,:,z) (is ionisation level, 0:n_z)
real*8, intent(out), optional :: GRC_out !< Generalized Radiational Coefficient at this density and temperature
real*8, intent(out), optional :: dGRC_dT_out !< Temperature gradient of GRC
real*8, intent(out), optional :: dGRC_dn_out !< Density gradient of GRC
real*8                        :: GRC, dGRC_dT, dGRC_dn !< Local generalized Radiational Coefficient


! If GRC exists and we are looking for a Z that is nonzero
if (allocated(a%GRC) .and. z .le. ubound(a%GRC,3) .and. z .ge. lbound(a%GRC,3)) then
  call SL2Dinterp(a%GRCFspline(z),temperature,density,fout=GRC,dfout_dx=dGRC_dT,dfout_dy=dGRC_dn)
  GRC     = 10.d0**GRC ! Converting from log10
  dGRC_dT = dGRC_dT * GRC / (10.0**temperature) ! Converting from loglog, log10 terms cancel
  dGRC_dn = dGRC_dn * GRC / (10.0**density) ! Converting from loglog, log10 terms cancel
else
  GRC     = 0.d0
  dGRC_dT = 0.d0
  dGRC_dn = 0.d0
endif

if (present(GRC_out)) GRC_out = GRC
if (present(dGRC_dT_out)) dGRC_dT_out = dGRC_dT
if (present(dGRC_dn_out)) dGRC_dn_out = dGRC_dn

end subroutine GRC_spl


end module mod_openadas

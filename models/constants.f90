!> Module containing constants which are used in the code
module constants
  implicit none
  public
  
#include "version.h"

  ! @name Mathematical and physical constants
  real*8,  parameter :: PI            = 3.1415926535897932385d0  !< PI
  real*8,  parameter :: TWOPI         = 6.2831853071795864769d0  !< TWOPI
  real*8,  parameter :: MU_ZERO       = 4.d-7*PI                 !< Magnetic constant  [Vs/Am]
  real*8,  parameter :: EPS_ZERO      = 8.854187817d-12          !< Vacuum permittivity [F/m]
  real*8,  parameter :: EL_CHG        = 1.602176565d-19          !< Elementary charge  [C]
  real*8,  parameter :: K_BOLTZ       = 1.3806488d-23            !< Boltzmann constant [J/K]
  real*8,  parameter :: C_LIGHT       = 299792458.               !< Speed of light [m/s]
  real*8,  parameter :: MASS_PROTON   = 1.67262178d-27           !< proton mass [kg]
  real*8,  parameter :: ATOMIC_MASS_UNIT = 1.660539040d-27       !< standardized mass unit [kg]
  real*8,  parameter :: MASS_ELECTRON = 9.10938291d-31           !< electron mass [kg]
  real*8,  parameter :: SPEED_OF_LIGHT = 2.997924580105029d+8    !< speed of light in (m/s)
  real*8,  parameter :: MOLE_NUMBER   = 6.02214076d23            !< The Avogadro constant

  !> @name Constants which describe the domain of a certain position (used by function which_domain)
  integer, parameter :: DOMAIN_PLASMA         = 0    !< Plasma region
  integer, parameter :: DOMAIN_SOL            = 1    !< Scrape-off layer
  integer, parameter :: DOMAIN_OUTER_SOL      = 2    !< Outer scrape-off layer (double-null)
  integer, parameter :: DOMAIN_UPPER_PRIVATE  = 3    !< Upper private flux region
  integer, parameter :: DOMAIN_LOWER_PRIVATE  = 4    !< Lower private flux region

  !> @name Parameters which describe the X-point case
  integer, parameter :: LOWER_XPOINT          = 1
  integer, parameter :: UPPER_XPOINT          = 2
  integer, parameter :: DOUBLE_NULL           = 3
  integer, parameter :: SYMMETRIC_XPOINT      = 100  ! Used for grid construction purposes; do not use as value for xcase in the input file!

end module constants

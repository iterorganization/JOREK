subroutine neo_coef( xpoint2, xcase2, Z, Z_xpoint, psi_node ,psi_axis,psi_bnd, &
     amu_neo_node, aki_neo_node)
  use phys_module


  implicit none

  integer :: i
! all the following variables are defined in mod_phys_module.f90
!  integer    :: num_neo_len
!  real*8, allocatable   :: num_neo_psi(:)  ! Radial positions of profile points (PsiN values)
!  real*8, allocatable   :: num_aki_value(:)  !numerical aki profile (PsiN values)
!  real*8, allocatable   :: num_amu_value(:)!numerical amu profile (PsiN values)
!  real*8, allocatable :: aki_neo_prof(:),amu_neo_prof(:)! aki_neo, amu_neo -profiles



! --- Routine parameters
!integer,  intent(in)   :: num_neo_len
!character, intent(in)  :: neo_file
logical, intent(in)     :: xpoint2
integer, intent(in)     :: xcase2
real*8,  intent(in)     :: Z
real*8,  intent(in)     :: Z_xpoint(2)
real*8,  intent(in)     :: psi_node
real*8,  intent(in)     :: psi_axis
real*8,  intent(in)     :: psi_bnd
real*8,  intent(out)    :: amu_neo_node, aki_neo_node  ! neoclassical coeffs at the given position.

! --- local variables
real*8  :: psi_n, prof1, prof2, sigz
real*8  :: tanh2, atn_z, tanh2_u, atn_z_u
! for interpolating numerical profiles
integer :: left, right, mid
real*8  :: aux1, aux2
real*8  :: Z_star, Z_star_u
psi_n = (psi_node - psi_axis) / (psi_bnd - psi_axis)
psi_n = max( min(psi_n, 2.), 0. )
  
! --- Interpolate profiles to position psi_n by bisections.
left  = 1
right = num_neo_len
do
  if ( right == left + 1 ) exit
  mid = (left + right) / 2
  if ( num_neo_psi(mid) >= psi_n ) then
    right = mid
  else
    left = mid
  end if
end do
aux1 = (psi_n - num_neo_psi(left)) / (num_neo_psi(right) - num_neo_psi(left))
aux2 = (1. - aux1)
prof1 = num_aki_value(left)*aux2 + num_aki_value(right)*aux1
prof2 = num_amu_value(left)*aux2 + num_amu_value(right)*aux1


! --- Additional explicit dependence of the profile on Z to ensure that the profile is
!     approximately zero in the private flux region below the x-point.
if (xpoint2) then
 
  sigz = 0.1d0

  if (xcase2 .eq. LOWER_XPOINT) then
     atn_z_u = 1.d0
  else
    Z_star_u = (Z-Z_xpoint(2))/sigz
    Z_star_u = min( max( Z_star_u, -40.d0), 40.d0) ! avoid floating-point exceptions
    
    tanh2_u  = tanh(Z_star_u)

    atn_z_u  = (0.5d0 - 0.5d0*tanh2_u)
  endif
  if (xcase .eq. UPPER_XPOINT) then
     atn_z = 1.d0
  else
    Z_star = (Z_xpoint(1)-Z)/sigz
    Z_star = min( max( Z_star, -40.d0), 40.d0) ! avoid floating-point exceptions

    tanh2  = tanh(Z_star)
    
    atn_z  = (0.5d0 - 0.5d0*tanh2)
  endif

  aki_neo_node = prof1 * atn_z * atn_z_u
  amu_neo_node = prof2 * atn_z * atn_z_u
else
  aki_neo_node = prof1
  amu_neo_node = prof2
end if

! set ki=0 and  mu_neo=0 in the SOL and the private region
if ( (psi_n .gt. 1.d0 ) &
     .or. ( xpoint2 .and. (xcase2 .ne. UPPER_XPOINT) .and. (Z .lt. Z_xpoint(1)) ) &
     .or. ( xpoint2 .and. (xcase2 .ne. LOWER_XPOINT) .and. (Z .gt. Z_xpoint(2)) ) ) then
  aki_neo_node = 0.d0
  amu_neo_node = 0.d0
endif

end subroutine neo_coef

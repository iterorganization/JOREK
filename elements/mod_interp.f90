!> Interpolation in the JOREK elements for values, deltas or x
module mod_interp
use data_structure
use mod_basisfunctions
use mod_parameters, only: n_period, n_tor
implicit none
private
public :: interp !< interp a specific harmonic in finite elements
public :: interp_delta !< interp a specific harmonic in finite elements, of the deltas
public :: interp_0 !< interp variable only, no derivatives at a specific position in domain
public :: interp_0_delta !< interp variable only, no derivatives at a specific position in domain, of the deltas
public :: interp_RZ !< Interpolate space only
public :: interp_PRZ !< interp variable + pos at values or deltas
public :: sincosperiod_moivre, mode_moivre !< public for regtesting, used by interp_PRZ

interface interp_RZ
  module procedure interp_RZ_0, interp_RZ_1, interp_RZ_2
end interface interp_RZ

interface interp_PRZ
  module procedure interp_PRZ_0, interp_PRZ_1, interp_PRZ_2
end interface interp_PRZ

contains

!> This subroutine interpolates some variables at a specific position within one element at a given position (s,t)
pure subroutine interp_PRZ_0(node_list, element_list, i_elm, i_v, n_v, s, t, phi, P, R, Z, deltas)
type (type_node_list),    intent(in)  :: node_list
type (type_element_list), intent(in)  :: element_list
integer,                  intent(in)  :: i_elm
integer,                  intent(in)  :: n_v, i_v(n_v)
real*8,                   intent(in)  :: s, t, phi
real*8,                   intent(out) :: P(n_v)
real*8,                   intent(out) :: R, Z
logical, optional, intent(in)         :: deltas

! --- Local variables
real*8  :: H(4,4), HZ(n_tor), dHZ(n_tor)
integer :: kv, iv, kf, i
real*8  :: values(n_tor,n_order+1,n_v,n_vertex_max)
real*8  :: xR(n_order+1,n_vertex_max), xZ(n_order+1,n_vertex_max)
real*8  :: sizes(n_order+1), v, vp
logical :: my_deltas

call basisfunctions(s,t,H)
H = transpose(H)
call sincosperiod_moivre(phi, HZ, dHZ) ! dHZ unused

P = 0.d0

my_deltas = .false.
if (present(deltas)) then
  if (deltas) my_deltas = .true.
end if

! Preload values and premultiply with sizes(:,kv)
do kv = 1,n_vertex_max  ! 4 vertices
  iv = element_list%element(i_elm)%vertex(kv)
  sizes(:) = element_list%element(i_elm)%size(kv,:)

  if (my_deltas) then
    do i = 1, n_v
      do kf=1,n_order+1
        values(1:n_tor,kf,i,kv) = node_list%node(iv)%deltas(1:n_tor,kf,i_v(i)) * sizes(kf)
      end do
    end do
  else
    do i = 1, n_v
      do kf=1,n_order+1
        values(1:n_tor,kf,i,kv) = node_list%node(iv)%values(1:n_tor,kf,i_v(i)) * sizes(kf)
      end do
    end do
  end if
  xR(:,kv) = node_list%node(iv)%x(1,:,1) * sizes(:)
  xZ(:,kv) = node_list%node(iv)%x(1,:,2) * sizes(:)
end do

! together 7%
R   = sum(xR*H)
Z   = sum(xZ*H)

! 40% exec time
do kv = 1, n_vertex_max
  do i = 1, n_v
    do kf = 1, n_order+1
      v = dot_product(values(1:n_tor,kf,i,kv),HZ(1:n_tor))
      P(i)     = P(i)     + v * H(kf, kv)
    enddo
  enddo
enddo
end subroutine interp_PRZ_0

!> This subroutine interpolates some variables at a specific position within one element at a given position (s,t)
pure subroutine interp_PRZ_1(node_list, element_list, i_elm, i_v, n_v, s, t, phi, P, P_s, P_t, P_phi, R, R_s, R_t, Z, Z_s, Z_t, deltas)
type (type_node_list),    intent(in)  :: node_list
type (type_element_list), intent(in)  :: element_list
integer,                  intent(in)  :: i_elm
integer,                  intent(in)  :: n_v, i_v(n_v)
real*8,                   intent(in)  :: s, t, phi
real*8,                   intent(out) :: P(n_v), P_s(n_v), P_t(n_v), P_phi(n_v)
real*8,                   intent(out) :: R, R_s, R_t, Z, Z_s, Z_t
logical, optional, intent(in)         :: deltas

! --- Local variables
real*8  :: H(4,4), H_s(4,4), H_t(4,4), HZ(n_tor), dHZ(n_tor)
integer :: kv, iv, kf, i
real*8  :: values(n_tor,n_order+1,n_v,n_vertex_max)
real*8  :: xR(n_order+1,n_vertex_max), xZ(n_order+1,n_vertex_max)
real*8  :: sizes(n_order+1), v, vp
logical :: my_deltas

! 7% exec time
call basisfunctions_T(s,t,H,H_s,H_t)

P = 0.d0; P_s = 0.d0; P_t = 0.d0; P_phi = 0.d0

! 7% exec time
call sincosperiod_moivre(phi, HZ, dHZ)

my_deltas = .false.
if (present(deltas)) then
  if (deltas) my_deltas = .true.
end if

! 30% exec time
! Preload values and premultiply with sizes(:,kv)
do kv = 1,n_vertex_max  ! 4 vertices
  iv = element_list%element(i_elm)%vertex(kv)
  sizes(:) = element_list%element(i_elm)%size(kv,:)

  if (my_deltas) then
    do i = 1, n_v
      do kf=1,n_order+1
        values(1:n_tor,kf,i,kv) = node_list%node(iv)%deltas(1:n_tor,kf,i_v(i)) * sizes(kf)
      end do
    end do
  else
    do i = 1, n_v
      do kf=1,n_order+1
        values(1:n_tor,kf,i,kv) = node_list%node(iv)%values(1:n_tor,kf,i_v(i)) * sizes(kf)
      end do
    end do
  end if
  xR(:,kv) = node_list%node(iv)%x(1,:,1) * sizes(:)
  xZ(:,kv) = node_list%node(iv)%x(1,:,2) * sizes(:)
end do

! together 7%
R   = sum(xR*H)
R_s = sum(xR*H_s)
R_t = sum(xR*H_t)
Z   = sum(xZ*H)
Z_s = sum(xZ*H_s)
Z_t = sum(xZ*H_t)

! 40% exec time
do kv = 1, n_vertex_max
  do i = 1, n_v
    do kf = 1, n_order+1
      v = dot_product(values(1:n_tor,kf,i,kv),HZ(1:n_tor))
      P(i)     = P(i)     + v * H(kf, kv)
      P_s(i)   = P_s(i)   + v * H_s(kf, kv)
      P_t(i)   = P_t(i)   + v * H_t(kf, kv)
      vp = dot_product(values(1:n_tor,kf,i,kv),dHZ(1:n_tor))
      P_phi(i) = P_phi(i) + vp * H(kf, kv)
    enddo
  enddo
enddo
end subroutine interp_PRZ_1


!> This subroutine interpolates some variables at a specific position within one element at a given position (s,t)
pure subroutine interp_PRZ_2(node_list, element_list, i_elm, i_v, n_v, s, t, phi, &
        P, P_s, P_t, P_phi, P_st, P_ss, P_tt, P_sphi, P_tphi, P_phiphi, R, R_s, R_t, R_st, R_ss, R_tt, &
        Z, Z_s, Z_t, Z_st, Z_ss, Z_tt, deltas)
type (type_node_list),    intent(in)  :: node_list
type (type_element_list), intent(in)  :: element_list
integer,                  intent(in)  :: i_elm
integer,                  intent(in)  :: n_v, i_v(n_v)
real*8,                   intent(in)  :: s, t, phi
real*8, dimension(n_v),   intent(out) :: P, P_s, P_t, P_phi, P_st, P_ss, P_tt, P_sphi, P_tphi, P_phiphi
real*8,                   intent(out) :: R, R_s, R_t, R_st, R_ss, R_tt, Z, Z_s, Z_t, Z_st, Z_ss, Z_tt
logical, optional, intent(in)         :: deltas

! --- Local variables
real*8, dimension(n_order+1,n_vertex_max) :: H, H_s, H_t, H_st, H_ss, H_tt
real*8, dimension(n_tor) :: HZ(n_tor), dHZ(n_tor), ddHZ(n_tor)
integer :: kv, iv, kf, i
real*8  :: values(n_tor,n_order+1,n_v,n_vertex_max)
real*8  :: xR(n_order+1,n_vertex_max), xZ(n_order+1,n_vertex_max)
real*8  :: sizes(n_order+1), v, vp, vpp
logical :: my_deltas

call basisfunctions_T(s,t,H,H_s,H_t,H_st,H_ss,H_tt)

P = 0.d0; P_s = 0.d0; P_t = 0.d0; P_st = 0.d0; P_ss = 0.d0; P_tt = 0.d0
P_phi = 0.d0; P_sphi = 0.d0; P_tphi = 0.d0; P_phiphi = 0.d0

call sincosperiod_moivre(phi, HZ, dHZ)
do i=1,n_tor
  ddHZ(i) = HZ(i)*(n_period*(i/2))**2
end do

my_deltas = .false.
if (present(deltas)) then
  if (deltas) my_deltas = .true.
end if

! Preload values and premultiply with sizes(:,kv)
do kv = 1,n_vertex_max  ! 4 vertices
  iv = element_list%element(i_elm)%vertex(kv)
  sizes(:) = element_list%element(i_elm)%size(kv,:)

  if (my_deltas) then
    do i = 1, n_v
      do kf=1,n_order+1
        values(1:n_tor,kf,i,kv) = node_list%node(iv)%deltas(1:n_tor,kf,i_v(i)) * sizes(kf)
      end do
    end do
  else
    do i = 1, n_v
      do kf=1,n_order+1
        values(1:n_tor,kf,i,kv) = node_list%node(iv)%values(1:n_tor,kf,i_v(i)) * sizes(kf)
      end do
    end do
  end if
  xR(:,kv) = node_list%node(iv)%x(1,:,1) * sizes(:)
  xZ(:,kv) = node_list%node(iv)%x(1,:,2) * sizes(:)
end do

R    = sum(xR*H)
R_s  = sum(xR*H_s)
R_t  = sum(xR*H_t)
R_st = sum(xR*H_st)
R_ss = sum(xR*H_ss)
R_tt = sum(xR*H_tt)
Z    = sum(xZ*H)
Z_s  = sum(xZ*H_s)
Z_t  = sum(xZ*H_t)
Z_st = sum(xZ*H_st)
Z_ss = sum(xZ*H_ss)
Z_tt = sum(xZ*H_tt)

do kv = 1, n_vertex_max
  do i = 1, n_v
    do kf = 1, n_order+1
      v = dot_product(values(1:n_tor,kf,i,kv),HZ(1:n_tor))
      P(i)     = P(i)     + v * H(kf, kv)
      P_s(i)   = P_s(i)   + v * H_s(kf, kv)
      P_t(i)   = P_t(i)   + v * H_t(kf, kv)
      vp = dot_product(values(1:n_tor,kf,i,kv),dHZ(1:n_tor))
      P_phi(i) = P_phi(i) + vp * H(kf, kv)

      P_st(i)  = P_st(i)  + v * H_st(kf, kv)
      P_ss(i)  = P_ss(i)  + v * H_ss(kf, kv)
      P_tt(i)  = P_tt(i)  + v * H_tt(kf, kv)

      P_sphi(i)   = P_sphi(i)   + vp * H_s(kf, kv)
      P_tphi(i)   = P_tphi(i)   + vp * H_t(kf, kv)
      vpp = dot_product(values(1:n_tor,kf,i,kv),ddHZ(1:n_tor))
      P_phiphi(i) = P_phiphi(i) + vpp * H(kf, kv)
    enddo
  enddo
enddo
end subroutine interp_PRZ_2

! Apply De Moivre formula to calculate the series of sines.
! Assumes that mode is of the form [0 1 1 2 2 3 3 4 4] ([0 4 4 8 8 12 12])
! This is roughly 3-4 times faster in my tests than just calculating the sines
! and cosines (even when that is vectorized). Perhaps that changes for n_tor >> 10
! I tested n_tor = 17.
pure subroutine sincosperiod_moivre(phi,HZ,dHZ)
  integer, parameter :: n_mode = (n_tor-1)/2 ! number of modes excluding 0
  real*8, intent(in) :: phi
  real*8, intent(out) :: HZ(n_tor), dHZ(n_tor)

  integer :: i
  HZ(1) = 1.d0
  dHZ(1) = 0.d0
  if (n_mode .gt. 0) then
    HZ(2) = cos(n_period*phi)
    HZ(3) = sin(n_period*phi)
    dHZ(2) = HZ(3)*(-n_period)
    dHZ(3) = HZ(2)*(n_period)

    do i=2,n_mode

      call moivre(HZ(2),HZ(3), HZ(2*i-2),HZ(2*i-1), HZ(2*i),HZ(2*i+1))

      dHZ(2*i)   = HZ(2*i+1)*(-n_period*i)
      dHZ(2*i+1) = HZ(2*i)*(n_period*i)

    end do
  end if
end subroutine sincosperiod_moivre

pure subroutine moivre(ar,ai,br,bi,or,oi)
  real*8, intent(in) :: ar, ai !< real and imag part of e^(i x)
  real*8, intent(in) :: br, bi !< real and imag part of e^(i y)
  real*8, intent(out) :: or, oi !< real and imag part of e^(i (x+y))
  or = ar*br - ai*bi
  oi = ai*br + ar*bi
end subroutine moivre

! Apply De Moivre formula to calculate mode*phi
! Assumes that mode is of the form [0 1 1 2 2 3 3 4 4] ([0 4 4 8 8 12 12])
! This is roughly 3-4 times faster in my tests than just calculating the sines
! and cosines (even when that is vectorized).
pure subroutine mode_moivre(phi,HZ)
  integer, parameter :: n_mode = (n_tor-1)/2 ! number of modes excluding 0
  real*8, intent(in) :: phi
  real*8, intent(out) :: HZ(n_tor)

  integer :: i
  HZ(1) = 1.d0
  if (n_mode .gt. 0) then
  
    HZ(2) = cos(n_period*phi)
    HZ(3) = sin(n_period*phi)

    do i=2,n_mode
      call moivre(HZ(2),HZ(3), HZ(2*i-2),HZ(2*i-1), HZ(2*i),HZ(2*i+1))
    end do
  end if
end subroutine mode_moivre



!> subroutine calculates the interpolation within one element (i_elm) for a given position
!> (s,t) in the local coordinates
pure subroutine interp(node_list, element_list, i_elm, i_var, i_harm, s, t, P, P_s, P_t, P_st, P_ss, P_tt)
type (type_node_list),    intent(in)  :: node_list
type (type_element_list), intent(in)  :: element_list
integer,                  intent(in)  :: i_elm
integer,                  intent(in)  :: i_var
integer,                  intent(in)  :: i_harm
real*8,                   intent(in)  :: s
real*8,                   intent(in)  :: t
real*8,                   intent(out) :: P, P_s, P_t, P_st, P_ss, P_tt

! --- Local variables
real*8 :: G(4,4), G_s(4,4), G_t(4,4), G_st(4,4), G_ss(4,4), G_tt(4,4)
integer :: kv, iv, kf 

call basisfunctions(s,t,G, G_s, G_t, G_st, G_ss, G_tt)

P = 0.d0; P_s = 0.d0; P_t = 0.d0; P_st = 0.d0; P_ss = 0.d0; P_tt = 0.d0

do kv = 1,n_vertex_max  ! 4 vertices
  iv = element_list%element(i_elm)%vertex(kv)  ! the node number
  do kf = 1, n_order+1       ! 4 basis functions

#ifdef fullmhd
    if (i_var == 710) then
      P    = P    + node_list%node(iv)%Fprof_eq(kf) * element_list%element(i_elm)%size(kv,kf) * G(kv,kf)
      P_s  = P_s  + node_list%node(iv)%Fprof_eq(kf) * element_list%element(i_elm)%size(kv,kf) * G_s(kv,kf)
      P_t  = P_t  + node_list%node(iv)%Fprof_eq(kf) * element_list%element(i_elm)%size(kv,kf) * G_t(kv,kf)
      P_st = P_st + node_list%node(iv)%Fprof_eq(kf) * element_list%element(i_elm)%size(kv,kf) * G_st(kv,kf)
      P_ss = P_ss + node_list%node(iv)%Fprof_eq(kf) * element_list%element(i_elm)%size(kv,kf) * G_ss(kv,kf)
      P_tt = P_tt + node_list%node(iv)%Fprof_eq(kf) * element_list%element(i_elm)%size(kv,kf) * G_tt(kv,kf)
    elseif (i_var == 711) then
      P    = P    + node_list%node(iv)%psi_eq(kf) * element_list%element(i_elm)%size(kv,kf) * G(kv,kf)
      P_s  = P_s  + node_list%node(iv)%psi_eq(kf) * element_list%element(i_elm)%size(kv,kf) * G_s(kv,kf)
      P_t  = P_t  + node_list%node(iv)%psi_eq(kf) * element_list%element(i_elm)%size(kv,kf) * G_t(kv,kf)
      P_st = P_st + node_list%node(iv)%psi_eq(kf) * element_list%element(i_elm)%size(kv,kf) * G_st(kv,kf)
      P_ss = P_ss + node_list%node(iv)%psi_eq(kf) * element_list%element(i_elm)%size(kv,kf) * G_ss(kv,kf)
      P_tt = P_tt + node_list%node(iv)%psi_eq(kf) * element_list%element(i_elm)%size(kv,kf) * G_tt(kv,kf)
    else
#endif
      P    = P    + node_list%node(iv)%values(i_harm,kf,i_var) * element_list%element(i_elm)%size(kv,kf) * G(kv,kf)
      P_s  = P_s  + node_list%node(iv)%values(i_harm,kf,i_var) * element_list%element(i_elm)%size(kv,kf) * G_s(kv,kf)
      P_t  = P_t  + node_list%node(iv)%values(i_harm,kf,i_var) * element_list%element(i_elm)%size(kv,kf) * G_t(kv,kf)
      P_st = P_st + node_list%node(iv)%values(i_harm,kf,i_var) * element_list%element(i_elm)%size(kv,kf) * G_st(kv,kf)
      P_ss = P_ss + node_list%node(iv)%values(i_harm,kf,i_var) * element_list%element(i_elm)%size(kv,kf) * G_ss(kv,kf)
      P_tt = P_tt + node_list%node(iv)%values(i_harm,kf,i_var) * element_list%element(i_elm)%size(kv,kf) * G_tt(kv,kf)
#ifdef fullmhd
    endif
#endif
  end do
end do
end subroutine interp


!> subroutine calculates the interpolation within one element (i_elm) for a given position
!> (s,t) in the local coordinates, of the deltas instead of the values
pure subroutine interp_delta(node_list, element_list, i_elm, i_var, i_harm, s, t, P, P_s, P_t, P_st, P_ss, P_tt)
type (type_node_list),    intent(in)  :: node_list
type (type_element_list), intent(in)  :: element_list
integer,                  intent(in)  :: i_elm
integer,                  intent(in)  :: i_var
integer,                  intent(in)  :: i_harm
real*8,                   intent(in)  :: s
real*8,                   intent(in)  :: t
real*8,                   intent(out) :: P, P_s, P_t, P_st, P_ss, P_tt

! --- Local variables
real*8 :: G(4,4), G_s(4,4), G_t(4,4), G_st(4,4), G_ss(4,4), G_tt(4,4)
integer :: kv, iv, kf 

call basisfunctions(s,t,G, G_s, G_t, G_st, G_ss, G_tt)

P = 0.d0; P_s = 0.d0; P_t = 0.d0; P_st = 0.d0; P_ss = 0.d0; P_tt = 0.d0

do kv = 1,n_vertex_max  ! 4 vertices
  iv = element_list%element(i_elm)%vertex(kv)  ! the node number
  do kf = 1, n_order+1       ! 4 basis functions
    P    = P    + node_list%node(iv)%deltas(i_harm,kf,i_var) * element_list%element(i_elm)%size(kv,kf) * G(kv,kf)
    P_s  = P_s  + node_list%node(iv)%deltas(i_harm,kf,i_var) * element_list%element(i_elm)%size(kv,kf) * G_s(kv,kf)
    P_t  = P_t  + node_list%node(iv)%deltas(i_harm,kf,i_var) * element_list%element(i_elm)%size(kv,kf) * G_t(kv,kf)
    P_st = P_st + node_list%node(iv)%deltas(i_harm,kf,i_var) * element_list%element(i_elm)%size(kv,kf) * G_st(kv,kf)
    P_ss = P_ss + node_list%node(iv)%deltas(i_harm,kf,i_var) * element_list%element(i_elm)%size(kv,kf) * G_ss(kv,kf)
    P_tt = P_tt + node_list%node(iv)%deltas(i_harm,kf,i_var) * element_list%element(i_elm)%size(kv,kf) * G_tt(kv,kf)
  end do
end do
end subroutine interp_delta


!> This subroutine interpolates some variables at a specific position within one element at a given position (s,t)
pure subroutine interp_0(node_list, element_list, i_elm, i_v, n_v, s, t, phi, P)
type (type_node_list),    intent(in)  :: node_list
type (type_element_list), intent(in)  :: element_list
integer,                  intent(in)  :: i_elm
integer,                  intent(in)  :: n_v, i_v(n_v)
real*8,                   intent(in)  :: s, t, phi
real*8,                   intent(out) :: P(n_v)

real*8  :: H(4,4), ss, mode
integer :: kv, iv, kf, m, i, i_harm, i_tor

call basisfunctions(s,t,H)

P = 0.d0

do kv = 1,n_vertex_max  ! 4 vertices
  iv = element_list%element(i_elm)%vertex(kv)  ! the node number
  do kf = 1, n_order+1       ! 4 basis functions
    ss  = element_list%element(i_elm)%size(kv,kf)
    do i = 1, n_v
      P(i)    = P(i)   + node_list%node(iv)%values(1,kf,i_v(i)) * ss * H(kv,kf)
      do i_tor = 1, (n_tor-1)/2
        i_harm = 2*i_tor
        mode = i_tor * n_period
        P(i)    = P(i)   + node_list%node(iv)%values(i_harm,kf,i_v(i))   * ss * H(kv,kf)   * cos(mode*phi)
        P(i)    = P(i)   + node_list%node(iv)%values(i_harm+1,kf,i_v(i)) * ss * H(kv,kf)   * sin(mode*phi)
      end do
    end do
  end do
end do
end subroutine interp_0



!> This subroutine interpolates some variables at a specific position within one element at a given position (s,t)
!> of the deltas, not the values
pure subroutine interp_0_delta(node_list, element_list, i_elm, i_v, n_v, s, t, phi, P)
type (type_node_list),    intent(in)  :: node_list
type (type_element_list), intent(in)  :: element_list
integer,                  intent(in)  :: i_elm
integer,                  intent(in)  :: n_v, i_v(n_v)
real*8,                   intent(in)  :: s, t, phi
real*8,                   intent(out) :: P(n_v)

real*8  :: H(4,4), ss, mode
integer :: kv, iv, kf, m, i, i_harm, i_tor

call basisfunctions(s,t,H)

P = 0.d0

do kv = 1,n_vertex_max  ! 4 vertices
  iv = element_list%element(i_elm)%vertex(kv)  ! the node number
  do kf = 1, n_order+1       ! 4 basis functions
    ss  = element_list%element(i_elm)%size(kv,kf)
    do i = 1, n_v
      P(i)    = P(i)   + node_list%node(iv)%deltas(1,kf,i_v(i)) * ss * H(kv,kf)
      do i_tor = 1, (n_tor-1)/2
        i_harm = 2*i_tor
        mode = i_tor * n_period
        P(i)    = P(i)   + node_list%node(iv)%deltas(i_harm,kf,i_v(i))   * ss * H(kv,kf)   * cos(mode*phi)
        P(i)    = P(i)   + node_list%node(iv)%deltas(i_harm+1,kf,i_v(i)) * ss * H(kv,kf)   * sin(mode*phi)
      end do
    end do
  end do
end do
end subroutine interp_0_delta




!> Calculates the interpolation within one element (i_elm) for a given position (s,t) in local coordinates
pure subroutine interp_RZ_0(node_list,element_list,i_elm,s,t,R,Z)
type (type_node_list),    intent(in)  :: node_list
type (type_element_list), intent(in)  :: element_list
integer,                  intent(in)  :: i_elm
real*8,                   intent(in)  :: s,t
real*8,                   intent(out) :: R, Z

! --- Local variables
real*8  :: G(4,4)
real*8  :: xx1, xx2, ss
integer :: kv, iv, kf

call basisfunctions(s,t,G)

R = 0.d0; Z = 0.d0

do kv = 1,n_vertex_max  ! 4 vertices
  iv = element_list%element(i_elm)%vertex(kv)  ! the node number
  do kf = 1, n_order+1       ! 4 basis functions
    xx1 = node_list%node(iv)%x(1,kf,1)
    xx2 = node_list%node(iv)%x(1,kf,2)
    ss  = element_list%element(i_elm)%size(kv,kf)
    
    R    = R    + xx1 * ss * G(kv,kf)
    Z    = Z    + xx2 * ss * G(kv,kf)
  end do
end do

end subroutine interp_RZ_0



!> This subroutine interpolates space a specific position within one element at a given position (s,t)
pure subroutine interp_RZ_1(node_list, element_list, i_elm, s, t, R, R_s, R_t, Z, Z_s, Z_t)
type (type_node_list),    intent(in)  :: node_list
type (type_element_list), intent(in)  :: element_list
integer,                  intent(in)  :: i_elm
real*8,                   intent(in)  :: s, t
real*8,                   intent(out) :: R, R_s, R_t, Z, Z_s, Z_t

! --- Local variables
real*8  :: H(4,4), H_s(4,4), H_t(4,4)
integer :: kv, iv
real*8  :: xR(n_order+1,n_vertex_max), xZ(n_order+1,n_vertex_max)
real*8  :: sizes(n_order+1)

call basisfunctions_T(s,t,H,H_s,H_t)

! Preload values and premultiply with sizes(:,kv)
do kv = 1,n_vertex_max  ! 4 vertices
  iv = element_list%element(i_elm)%vertex(kv)
  sizes(:) = element_list%element(i_elm)%size(kv,:)
  xR(:,kv) = node_list%node(iv)%x(1,:,1) * sizes(:)
  xZ(:,kv) = node_list%node(iv)%x(1,:,2) * sizes(:)
end do

R   = sum(xR*H)
R_s = sum(xR*H_s)
R_t = sum(xR*H_t)
Z   = sum(xZ*H)
Z_s = sum(xZ*H_s)
Z_t = sum(xZ*H_t)

end subroutine interp_RZ_1



!> Calculates the interpolation within one element (i_elm) for a given position (s,t) in local coordinates
pure subroutine interp_RZ_2(node_list,element_list,i_elm,s,t,R,R_s,R_t,R_st,R_ss,R_tt,Z,Z_s,Z_t,Z_st,Z_ss,Z_tt)
type (type_node_list),    intent(in)  :: node_list
type (type_element_list), intent(in)  :: element_list
integer,                  intent(in)  :: i_elm
real*8,                   intent(in)  :: s,t
real*8,                   intent(out) :: R, R_s, R_t, R_st, R_ss, R_tt, Z, Z_s, Z_t, Z_st, Z_ss, Z_tt

! --- Local variables
real*8  :: H(4,4), H_s(4,4), H_t(4,4), H_st(4,4), H_ss(4,4), H_tt(4,4)
integer :: kv, iv
real*8  :: xR(n_order+1,n_vertex_max), xZ(n_order+1,n_vertex_max)
real*8  :: sizes(n_order+1)

call basisfunctions_T(s,t,H,H_s,H_t,H_st,H_ss,H_tt)

! Preload values and premultiply with sizes(:,kv)
do kv = 1,n_vertex_max  ! 4 vertices
  iv = element_list%element(i_elm)%vertex(kv)
  sizes(:) = element_list%element(i_elm)%size(kv,:)
  xR(:,kv) = node_list%node(iv)%x(1,:,1) * sizes(:)
  xZ(:,kv) = node_list%node(iv)%x(1,:,2) * sizes(:)
end do

R    = sum(xR*H)
R_s  = sum(xR*H_s)
R_t  = sum(xR*H_t)
R_st = sum(xR*H_st)
R_ss = sum(xR*H_ss)
R_tt = sum(xR*H_tt)
Z    = sum(xZ*H)
Z_s  = sum(xZ*H_s)
Z_t  = sum(xZ*H_t)
Z_st = sum(xZ*H_st)
Z_ss = sum(xZ*H_ss)
Z_tt = sum(xZ*H_tt)

end subroutine interp_RZ_2
end module mod_interp

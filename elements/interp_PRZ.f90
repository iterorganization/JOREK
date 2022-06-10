!> This subroutine interpolates some variables at a specific position within one element at a given position (s,t)
subroutine interp_PRZ(node_list, element_list, i_elm, i_v, n_v, s, t, phi, P, P_s, P_t, P_phi, R, R_s, R_t, Z, Z_s, Z_t)

use data_structure
use phys_module, only : mode
use mod_basisfunctions
use mod_parameters, only: n_degrees
implicit none

! --- Routine parameters
type (type_node_list),    intent(in)  :: node_list
type (type_element_list), intent(in)  :: element_list
integer,                  intent(in)  :: i_elm
integer,                  intent(in)  :: n_v, i_v(n_v)
real*8,                   intent(in)  :: s, t, phi
real*8,                   intent(out) :: P(n_v), P_s(n_v), P_t(n_v)
real*8,                   intent(out) :: R, R_s, R_t, Z, Z_s, Z_t
real*8,                   intent(out) :: P_phi(n_v)

! --- Local variables
real*8  :: H(4,n_degrees), H_s(4,n_degrees), H_t(4,n_degrees), xx1, xx2, ss
integer :: kv, iv, kf, m, i, i_harm, i_tor

call basisfunctions3(s,t,H,H_s,H_t)

P = 0.d0; P_s = 0.d0; P_t = 0.d0;
R = 0.d0; R_s = 0.d0; R_t = 0.d0;
Z = 0.d0; Z_s = 0.d0; Z_t = 0.d0;
P_phi = 0.d0

do kv = 1,n_vertex_max  ! 4 vertices

  iv = element_list%element(i_elm)%vertex(kv)  ! the node number

  do kf = 1, n_degrees       ! basis functions

    xx1 = node_list%node(iv)%x(1,kf,1)
    xx2 = node_list%node(iv)%x(1,kf,2)
    ss  = element_list%element(i_elm)%size(kv,kf)

    R    = R    + xx1 * ss * H(kv,kf)
    R_s  = R_s  + xx1 * ss * H_s(kv,kf)
    R_t  = R_t  + xx1 * ss * H_t(kv,kf)

    Z    = Z    + xx2 * ss * H(kv,kf)
    Z_s  = Z_s  + xx2 * ss * H_s(kv,kf)
    Z_t  = Z_t  + xx2 * ss * H_t(kv,kf)

    do i = 1, n_v

      P(i)    = P(i)   + node_list%node(iv)%values(1,kf,i_v(i)) * ss * H(kv,kf)
      P_s(i)  = P_s(i) + node_list%node(iv)%values(1,kf,i_v(i)) * ss * H_s(kv,kf)
      P_t(i)  = P_t(i) + node_list%node(iv)%values(1,kf,i_v(i)) * ss * H_t(kv,kf)

      do i_tor = 1, (n_tor-1)/2

        i_harm = 2*i_tor

        P(i)    = P(i)   + node_list%node(iv)%values(i_harm,kf,i_v(i))   * ss * H(kv,kf)   * cos(mode(i_harm)*phi)
        P_s(i)  = P_s(i) + node_list%node(iv)%values(i_harm,kf,i_v(i))   * ss * H_s(kv,kf) * cos(mode(i_harm)*phi)
        P_t(i)  = P_t(i) + node_list%node(iv)%values(i_harm,kf,i_v(i))   * ss * H_t(kv,kf) * cos(mode(i_harm)*phi)
        P_phi(i) = P_phi(i) + node_list%node(iv)%values(i_harm,kf,i_v(i)) &
              * ss * H(kv,kf) * sin(mode(i_harm)*phi) * (-mode(i_harm))

	P(i)    = P(i)   + node_list%node(iv)%values(i_harm+1,kf,i_v(i)) * ss * H(kv,kf)   * sin(mode(i_harm+1)*phi)
        P_s(i)  = P_s(i) + node_list%node(iv)%values(i_harm+1,kf,i_v(i)) * ss * H_s(kv,kf) * sin(mode(i_harm+1)*phi)
        P_t(i)  = P_t(i) + node_list%node(iv)%values(i_harm+1,kf,i_v(i)) * ss * H_t(kv,kf) * sin(mode(i_harm+1)*phi)
        P_phi(i) = P_phi(i) + node_list%node(iv)%values(i_harm+1,kf,i_v(i)) &
            * ss * H(kv,kf) * cos(mode(i_harm)*phi) * (mode(i_harm))

      enddo

    enddo

  end do

end do

return
end subroutine interp_PRZ

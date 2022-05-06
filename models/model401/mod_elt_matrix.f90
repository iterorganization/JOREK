! This module contains nothing (just a wrapper) but it is needed by
! construct_matrix for the other models.
! Can be removed once the other models have also combined element_matrix and
! element_matrix_fft.
module mod_elt_matrix
contains

  subroutine element_matrix(element, nodes, xpoint2, xcase2, R_axis, Z_axis, psi_axis, psi_bnd, R_xpoint, Z_xpoint, &
                            ELM, RHS, tid, i_tor_min, i_tor_max, aux_nodes)
  !--------------------------------------------------------------------------
  ! This is just a wrapper to the real routine since I combined both into one
  !--------------------------------------------------------------------------

    use data_structure
    use mod_elt_matrix_fft

    implicit none

    type (type_element)                      :: element
    type (type_node)                         :: nodes(n_vertex_max)
    type (type_node),optional                :: aux_nodes(n_vertex_max)

    integer                                  :: xcase2
    logical                                  :: xpoint2
    real*8                                   :: R_axis, Z_axis, psi_axis, psi_bnd, R_xpoint(2), Z_xpoint(2)
    real*8, dimension (:,:), allocatable     :: ELM
    real*8, dimension (:)  , allocatable     :: RHS
    integer, intent(in)                      :: tid, i_tor_min, i_tor_max

    call element_matrix_fft(element, nodes, xpoint2, xcase2, R_axis, Z_axis, psi_axis, psi_bnd, R_xpoint, Z_xpoint, ELM, RHS, tid,   &
                            thread_struct(tid)%ELM_p, thread_struct(tid)%ELM_n, thread_struct(tid)%ELM_k, thread_struct(tid)%ELM_kn, &
                            thread_struct(tid)%RHS_p, thread_struct(tid)%RHS_k, thread_struct(tid)%eq_g, thread_struct(tid)%eq_s,    &
                            thread_struct(tid)%eq_t, thread_struct(tid)%eq_p, thread_struct(tid)%eq_ss, thread_struct(tid)%eq_st,    &
                            thread_struct(tid)%eq_tt, thread_struct(tid)%delta_g, thread_struct(tid)%delta_s,                        & 
                            thread_struct(tid)%delta_t, i_tor_min, i_tor_max, aux_nodes, thread_struct(tid)%ELM_pnn)

    return

  end subroutine element_matrix

end module mod_elt_matrix


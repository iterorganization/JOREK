module mod_assembly
implicit none
!> This module provides routines to factorize code in boundary conditions from each models

contains
  !>
  !! Subroutine: boundary_conditions_add_one_entry
  !!
  !! Add one entry to the product and/or harminic matrix if index_node is local.
  !!
  !! @param index_node  row node index
  !! @param k           row var index
  !! @param in          row tor index
  !! @param index_node2 col node index
  !! @param k2          col var index
  !! @param in2         col tor index
  !! @param zbig        value
  !! @param solve_only  Do not add to harmonic matrix if .true.
  !! @param gmres       Do not add to product matrix if .false.
  !! @param use_murge   Use murge interface.
  !! @param use_murge_element Murge interface with elementary matrices.
  !! @param index_min   Minimal local element index
  !! @param index_max   Maximal local element index
  !!
  subroutine boundary_conditions_add_one_entry( &
       &   index_node,  k,  in,                 &
       &   index_node2, k2, in2,                &
       &   zbig, solve_only, gmres,             &
       &   index_min, index_max,                & 
       &   ijA_index, ijA_size, irn_jcn,        & 
       &   irn, jcn, A_mat, i_tor_min, i_tor_max)
    use mod_parameters
    use mod_locate_irn_jcn
    use mod_integer_types

    integer,               intent(in)                 :: k,  in
    integer,               intent(in)                 :: k2, in2
    integer,               intent(in)                 :: index_node
    integer,               intent(in)                 :: index_node2
    real*8,                intent(in)                 :: zbig
    logical,               intent(in)                 :: solve_only, gmres
    integer,               intent(in)                 :: index_min, index_max
    integer,               intent(in)                 :: i_tor_min, i_tor_max 
    integer(kind=int_all), intent(in),    allocatable :: ijA_index(:,:), ijA_size(:), irn_jcn(:,:) 
    integer(kind=int_all), intent(inout), allocatable :: irn(:), jcn(:) 
    real*8,                intent(inout), allocatable :: A_mat(:) 
    
    logical                                           :: is_local
    integer(kind=int_all)                             :: ija_position, ilarge_vp
    integer                                           :: n_tor_local

    if ( (k==0) .or. (k2==0) ) return ! ignore calls for model family extensions not in use (variable number zero)

    n_tor_local = i_tor_max - i_tor_min +1
    if ((index_node .ge. index_min) .and. (index_node .le. index_max)) then

       call locate_irn_jcn(index_node,index_node2,index_min,index_max,ijA_position,& 
                                    ijA_index, ijA_size, irn_jcn)
                             
       !-------- index dans A_mat
       ilarge_vp  = ijA_position  - 1 + ((k-1)*n_tor_local + in-i_tor_min ) * n_var*n_tor_local + (k2-1)*n_tor_local + in2&
                    -i_tor_min + 1 
                               
                             
       irn(ilarge_vp) =  n_tor_local * n_var * (index_node -1) + (k -1)*n_tor_local + in - i_tor_min + 1
       jcn(ilarge_vp) =  n_tor_local * n_var * (index_node2-1) + (k2-1)*n_tor_local + in2 - i_tor_min + 1
       A_mat(ilarge_vp)   = ZBIG
    endif
  end subroutine boundary_conditions_add_one_entry

  !>
  !! Subroutine: boundary_conditions_add_RHS
  !!
  !! Add one entry to the product and/or harminic matrix if index_node is local.
  !!
  !! @param index_node  row node index
  !! @param k           row var index
  !! @param in          row tor index
  !! @param index_node2 col node index
  !! @param k2          col var index
  !! @param in2         col tor index
  !! @param zbig        value
  !! @param solve_only  Do not add to harmonic matrix if .true.
  !! @param gmres       Do not add to product matrix if .false.
  !! @param index_min   Minimal local element index
  !! @param index_max   Maximal local element index
  !!
  subroutine boundary_conditions_add_RHS( &
       &   index_node, k, in,             &
       &   index_min, index_max,          &
       &   rhs_loc,  val,                 &
       &   i_tor_min, i_tor_max)
    use mod_parameters
    use mod_integer_types
    integer,               intent(in)    :: k,  in
    integer,               intent(in)    :: index_node
    integer,               intent(in)    :: index_min, index_max
    integer,               intent(in)    :: i_tor_min, i_tor_max 
    real*8,                intent(in)    :: val
    real*8,                intent(inOUT) :: rhs_loc(*) 
    integer                              :: n_tor_local 
    logical                              :: is_local

    if ( (k==0) ) return ! ignore calls for model family extensions not in use (variable number zero)

    n_tor_local = i_tor_max - i_tor_min +1
    if ((index_node .ge. index_min) .and. (index_node .le. index_max)) then
       RHS_loc(n_tor_local*n_var * (index_node-1) + (k-1)*n_tor_local + in - i_tor_min + 1) = val
    endif
  end subroutine boundary_conditions_add_RHS
end module mod_assembly

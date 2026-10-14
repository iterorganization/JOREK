!*******************************************************************************
!* Module: mod_boundary_ndotB                                                  *
!*                                                                             *
!* Store and retrieve n.B/|B| values per boundary node for stellarator SBC.    *
!* Uses MPI_Allreduce to share values across all ranks.                        *
!* Supports toroidal variation: stores n.B per (node, plane), then computes   *
!* the (nonlinear) v_par target in physical space and Fourier-decomposes      *
!* -- the target cannot be obtained by transforming n.B's own Fourier         *
!* coefficients, since tanh(.) does not commute with the DFT.                 *
!*******************************************************************************

module mod_boundary_ndotB

  implicit none
  private
  
  public :: init_boundary_ndotB
  public :: accumulate_ndotB_at_node_plane
  public :: finalize_boundary_ndotB
  public :: get_vpar_target_for_column
  public :: get_vpar_target_dT_for_column
  
  ! Per-node-plane storage (physical space, toroidal variation)
  real*8, allocatable, save :: ndotB_per_node_plane(:,:)        ! (n_nodes, n_plane)
  real*8, allocatable, save :: ndotB_count_per_node_plane(:,:)  ! (n_nodes, n_plane)
  
  ! Fourier coefficients of vpar_target (tanh smoothing applied in physical space first) and its T-derivative
  real*8, allocatable, save :: vpar_target_fourier_cos(:,:)     ! (n_nodes, n_tor)
  real*8, allocatable, save :: vpar_target_fourier_sin(:,:)     ! (n_nodes, n_tor)
  real*8, allocatable, save :: vpar_target_dT_fourier_cos(:,:)
  real*8, allocatable, save :: vpar_target_dT_fourier_sin(:,:)

contains

subroutine init_boundary_ndotB(n_nodes, n_plane_in, n_tor_in)
  implicit none
  integer, intent(in) :: n_nodes, n_plane_in, n_tor_in
  
  if (allocated(ndotB_per_node_plane))       deallocate(ndotB_per_node_plane)
  if (allocated(ndotB_count_per_node_plane)) deallocate(ndotB_count_per_node_plane)
  if (allocated(vpar_target_fourier_cos))    deallocate(vpar_target_fourier_cos)
  if (allocated(vpar_target_fourier_sin))    deallocate(vpar_target_fourier_sin)
  if (allocated(vpar_target_dT_fourier_cos)) deallocate(vpar_target_dT_fourier_cos)
  if (allocated(vpar_target_dT_fourier_sin)) deallocate(vpar_target_dT_fourier_sin)
  
  allocate(ndotB_per_node_plane(n_nodes, n_plane_in));       ndotB_per_node_plane       = 0.d0
  allocate(ndotB_count_per_node_plane(n_nodes, n_plane_in)); ndotB_count_per_node_plane = 0.d0
  allocate(vpar_target_fourier_cos(n_nodes, n_tor_in));      vpar_target_fourier_cos    = 0.d0
  allocate(vpar_target_fourier_sin(n_nodes, n_tor_in));      vpar_target_fourier_sin    = 0.d0
  allocate(vpar_target_dT_fourier_cos(n_nodes, n_tor_in));   vpar_target_dT_fourier_cos = 0.d0
  allocate(vpar_target_dT_fourier_sin(n_nodes, n_tor_in));   vpar_target_dT_fourier_sin = 0.d0
  
end subroutine init_boundary_ndotB


subroutine accumulate_ndotB_at_node_plane(inode, mp, ndotB_value)
  !---------------------------------------------------------------------------
  ! Accumulate ndotB at a specific node and toroidal plane. Called once per
  ! boundary Gauss point (ms==1 only, see mod_boundary_matrix_open.f90) during
  ! element assembly; summed here, divided by the hit-count in finalize().
  !---------------------------------------------------------------------------
  implicit none
  integer, intent(in) :: inode, mp
  real*8, intent(in) :: ndotB_value
  
  if (.not. allocated(ndotB_per_node_plane)) return
  if (inode < 1 .or. inode > ubound(ndotB_per_node_plane,1)) return
  if (mp    < 1 .or. mp    > ubound(ndotB_per_node_plane,2)) return
  
  ndotB_per_node_plane(inode, mp)       = ndotB_per_node_plane(inode, mp)       + ndotB_value
  ndotB_count_per_node_plane(inode, mp) = ndotB_count_per_node_plane(inode, mp) + 1.d0
  
end subroutine accumulate_ndotB_at_node_plane

subroutine finalize_boundary_ndotB()
  use mpi
  use nodes_elements
  use mod_parameters, only: n_period
  use corr_neg,       only: corr_neg_temp, dcorr_neg_temp_dT
  use constants, only: PI
  use phys_module,    only: vpar_sbc_alpha0, &
                            GAMMA, vpar_sbc_enable, vpar_sbc_T_floor
  use mod_model_settings, only: var_T, var_Vpar
  implicit none
  integer :: i, mp, in, ierr, my_id, n_nodes_loc, n_plane_loc, n_tor_loc, n_with_data
  real*8  :: phi, ndotB_val, cos_n, sin_n
  real*8  :: alpha0_rad
  real*8  :: vpar_target_val, dvpar_target_dT_val, dT_local_dT_DOF, cs, T_local, cs_for_deriv
  real*8  :: ndotB_min, ndotB_max, ndotB_sum
  real*8, allocatable :: ndotB_plane_global(:,:), count_plane_global(:,:)
  logical, save :: diag_printed = .false.
  logical, save :: warned_2T    = .false.
  
  if (.not. allocated(ndotB_per_node_plane)) return

  n_nodes_loc = ubound(ndotB_per_node_plane,1)
  n_plane_loc = ubound(ndotB_per_node_plane,2)
  n_tor_loc   = ubound(vpar_target_fourier_cos,2)
  
  call MPI_Comm_rank(MPI_COMM_WORLD, my_id, ierr)

  allocate(ndotB_plane_global(n_nodes_loc, n_plane_loc))
  allocate(count_plane_global(n_nodes_loc, n_plane_loc))

  call MPI_Allreduce(ndotB_per_node_plane, ndotB_plane_global, n_nodes_loc*n_plane_loc, &
                      MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_Allreduce(ndotB_count_per_node_plane, count_plane_global, n_nodes_loc*n_plane_loc, &
                      MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  
  alpha0_rad = vpar_sbc_alpha0 * PI / 180.d0

  vpar_target_fourier_cos    = 0.d0
  vpar_target_fourier_sin    = 0.d0
  vpar_target_dT_fourier_cos = 0.d0
  vpar_target_dT_fourier_sin = 0.d0

  n_with_data = 0
  ndotB_sum   = 0.d0
  ndotB_min   =  1.d10
  ndotB_max   = -1.d10

  do i = 1, n_nodes_loc

    ! Average the raw accumulation into physical-space n.B per plane,
    ! tracking global min/max/avg for the diagnostic print below.
    do mp = 1, n_plane_loc
      if (count_plane_global(i,mp) > 0.5d0) then
        ndotB_per_node_plane(i,mp) = ndotB_plane_global(i,mp) / count_plane_global(i,mp)
        n_with_data = n_with_data + 1
        ndotB_sum   = ndotB_sum + ndotB_per_node_plane(i,mp)
        ndotB_min   = min(ndotB_min, ndotB_per_node_plane(i,mp))
        ndotB_max   = max(ndotB_max, ndotB_per_node_plane(i,mp))
      else
        ndotB_per_node_plane(i,mp) = 0.d0
      endif
    enddo

    ! Local sound speed: always from the node's own T DOF, floor-corrected
    ! near zero. 2T not yet implemented, use floor value as placeholder
    if (var_T .gt. 0) then
      T_local         = corr_neg_temp(node_list%node(i)%values(1,1,var_T))
      dT_local_dT_DOF = dcorr_neg_temp_dT(node_list%node(i)%values(1,1,var_T))
    else
      T_local         = vpar_sbc_T_floor
      dT_local_dT_DOF = 0.d0
      if (vpar_sbc_enable .and. my_id == 0 .and. .not. warned_2T) then
        write(*,'(A)') "WARNING: vpar_sbc_enable=.true. in a 2T build (var_T<=0)."
        write(*,'(A)') "  v_par SBC sound speed is NOT coupled to Ti/Te yet -- using vpar_sbc_T_floor as a placeholder."
        warned_2T = .true.
      endif
    endif

    ! Compute sound speed
    cs           = sqrt(GAMMA * T_local)
    cs_for_deriv = sqrt(GAMMA * max(T_local, vpar_sbc_T_floor))  ! protects 1/cs in the Jacobian only

    do in = 1, n_tor_loc
      do mp = 1, n_plane_loc
        phi       = 2.d0 * PI * dble(mp-1) / dble(n_plane_loc * n_period)
        ndotB_val = ndotB_per_node_plane(i, mp)
        
        ! Compute vpar_target in physical space
        ! Smooth: vpar = cs * tanh(ndotB / sin(alpha0)); avoids Gibbs at sign flips
        vpar_target_val     = cs * tanh(ndotB_val / sin(alpha0_rad))

        ! T derivative with chain rule
        dvpar_target_dT_val = (GAMMA / (2.d0*cs_for_deriv)) * tanh(ndotB_val / sin(alpha0_rad)) &
                                * dT_local_dT_DOF
        
        ! Fourier coefficient for mode (in-1) [0-indexed internally]
        ! in=1 is n=0 mode (constant), in=2 is n=1 mode, etc.
        cos_n = cos(dble(in-1) * dble(n_period) * phi)
        sin_n = sin(dble(in-1) * dble(n_period) * phi)
        
        ! Fourier transform of vpar_target (for BC application - nonlinear!) and derivative
        vpar_target_fourier_cos(i,in)    = vpar_target_fourier_cos(i,in)    + vpar_target_val    * cos_n
        vpar_target_fourier_sin(i,in)    = vpar_target_fourier_sin(i,in)    + vpar_target_val    * sin_n
        vpar_target_dT_fourier_cos(i,in) = vpar_target_dT_fourier_cos(i,in) + dvpar_target_dT_val * cos_n
        vpar_target_dT_fourier_sin(i,in) = vpar_target_dT_fourier_sin(i,in) + dvpar_target_dT_val * sin_n
      enddo
      
      ! Normalize by number of planes (DFT normalization)
      vpar_target_fourier_cos(i,in)    = vpar_target_fourier_cos(i,in)    / dble(n_plane_loc)
      vpar_target_fourier_sin(i,in)    = vpar_target_fourier_sin(i,in)    / dble(n_plane_loc)
      vpar_target_dT_fourier_cos(i,in) = vpar_target_dT_fourier_cos(i,in) / dble(n_plane_loc)
      vpar_target_dT_fourier_sin(i,in) = vpar_target_dT_fourier_sin(i,in) / dble(n_plane_loc)
      
      ! For n=0 mode (in=1), the normalization is just the average
      ! For n>0 modes, multiply by 2 (standard DFT convention for real signals)
      if (in > 1) then
        vpar_target_fourier_cos(i,in)    = vpar_target_fourier_cos(i,in)    * 2.d0
        vpar_target_fourier_sin(i,in)    = vpar_target_fourier_sin(i,in)    * 2.d0
        vpar_target_dT_fourier_cos(i,in) = vpar_target_dT_fourier_cos(i,in) * 2.d0
        vpar_target_dT_fourier_sin(i,in) = vpar_target_dT_fourier_sin(i,in) * 2.d0
      endif
    enddo

    if (my_id == 0 .and. node_list%node(i)%boundary .ne. 0 .and. .not. diag_printed) then
      write(*,"(A)") "Mod_boundary_ndotB:"
      write(*,'(A,I6,A,E12.4,A,E12.4,A,E12.4,A,E12.4)') &
        " i=", i, " T_local=", T_local, " cs=", cs, &
        " vpar_target_fourier_cos(DC)=", vpar_target_fourier_cos(i,1), &
        " vpar_target_fourier_sin(DC)=", vpar_target_fourier_sin(i,1)
      if (var_Vpar .gt. 0) then
        write(*,'(A,E12.4,A,E12.4)') &
          "   vpar_target_dT_fourier_cos(DC)=", vpar_target_dT_fourier_cos(i,1), &
          " vpar_current(DC)=", node_list%node(i)%values(1,1,var_Vpar)
      endif
      if (T_local .lt. vpar_sbc_T_floor) then
        write(*,'(A,E12.4,A,E12.4)') "WARNING: node T_local=", T_local, &
          " below vpar_sbc_T_floor=", vpar_sbc_T_floor, " -- v_par SBC Jacobian capped for stability."
      endif
      diag_printed = .true.
    endif

  enddo

  deallocate(ndotB_plane_global, count_plane_global)

  if (my_id == 0 .and. n_with_data > 0) then
    write(*,'(A,I6,A,E12.4,A,E12.4,A,E12.4)') &
      " ndotB stats (MPI reduced, per node*plane): n_samples=", n_with_data, &
      " min=", ndotB_min, " max=", ndotB_max, " avg=", ndotB_sum/dble(n_with_data)
  endif
  
end subroutine finalize_boundary_ndotB


function get_vpar_target_for_column(inode, in) result(vpar_target)
  !---------------------------------------------------------------------------
  ! Return vpar_target for JOREK column index 'in' at boundary node.
  ! JOREK basis: in=1 DC, (2,3) first harmonic cos/-sin, (4,5) second, etc.
  !---------------------------------------------------------------------------
  implicit none
  integer, intent(in) :: inode, in
  real*8 :: vpar_target
  integer :: k_fourier

  if (.not. allocated(vpar_target_fourier_cos)) then
    vpar_target = 0.d0; return
  endif
  if (inode < 1 .or. inode > ubound(vpar_target_fourier_cos,1)) then
    vpar_target = 0.d0; return
  endif

  if (in .eq. 1) then
    vpar_target = vpar_target_fourier_cos(inode, 1)   ! DC; sin(0)=0 identically
    return
  endif

  k_fourier = in / 2 + 1
  if (k_fourier > ubound(vpar_target_fourier_cos,2)) then
    vpar_target = 0.d0; return
  endif

  if (mod(in, 2) .eq. 0) then
    vpar_target = vpar_target_fourier_cos(inode, k_fourier)
  else
    vpar_target = -vpar_target_fourier_sin(inode, k_fourier)   ! JOREK basis is -sin
  endif

end function get_vpar_target_for_column

function get_vpar_target_dT_for_column(inode, in) result(vpar_target_dT)
  !---------------------------------------------------------------------------
  ! Return target_T for derivatives, for JOREK column index 'in' at boundary node.
  implicit none
  integer, intent(in) :: inode, in
  real*8 :: vpar_target_dT
  integer :: k_fourier

  if (.not. allocated(vpar_target_dT_fourier_cos)) then
    vpar_target_dT = 0.d0; return
  endif
  if (inode < 1 .or. inode > ubound(vpar_target_dT_fourier_cos,1)) then
    vpar_target_dT = 0.d0; return
  endif

  if (in .eq. 1) then
    vpar_target_dT = vpar_target_dT_fourier_cos(inode, 1)
    return
  endif

  k_fourier = in / 2 + 1
  if (k_fourier > ubound(vpar_target_dT_fourier_cos,2)) then
    vpar_target_dT = 0.d0; return
  endif

  if (mod(in, 2) .eq. 0) then
    vpar_target_dT = vpar_target_dT_fourier_cos(inode, k_fourier)
  else
    vpar_target_dT = -vpar_target_dT_fourier_sin(inode, k_fourier)
  endif

end function get_vpar_target_dT_for_column

end module mod_boundary_ndotB
module mod_boundary_conditions
contains

!*******************************************************************************
!* Subroutine: boundary_condition                                              *
!*******************************************************************************
!*                                                                             *
!* Add boundary condition on the matrix.                                       *
!*                                                                             *
!* Parameters:                                                                 *
!*   my_id        - Identifier of the node in MPI_COMM_WORLD                   *
!*   node_list    - List of nodes                                              *
!*   element_list - List of all elements                                       *
!*   local_elms   - List of local elements                                     *
!*   n_local_elms - Number of local elements                                   *
!*   index_min    - Minimal index of local elements                            *
!*   index_max    - Maximal index of local elements (                          *
!*   xpoint2      -                                                            *
!*   xcase2       -                                                            *
!*   psi_axis     -                                                            *
!*   psi_bnd      -                                                            *
!*   Z_xpoint     -                                                            *
!*                                                                             *
!*******************************************************************************

subroutine boundary_conditions( my_id, node_list, element_list, bnd_node_list, local_elms,& 
                                n_local_elms, index_min, index_max, rhs_loc, xpoint2,     &
                                xcase2, R_axis, Z_axis, psi_axis, psi_bnd,                &
                                R_xpoint, Z_xpoint, psi_xpoint, a_mat)

use mod_assembly, only : boundary_conditions_add_one_entry, boundary_conditions_add_RHS
use data_structure
use vacuum, ONLY: is_freebound
use phys_module, only: F0, GAMMA, freeboundary, RMP_on, psi_RMP_cos, dpsi_RMP_cos_dR, dpsi_RMP_cos_dZ, &
       psi_RMP_sin, dpsi_RMP_sin_dR, dpsi_RMP_sin_dZ, t_now, RMP_growth_rate, RMP_ramp_up_time,            &
       RMP_start_time, tstep, RMP_har_cos, RMP_har_sin, T_min,                                             &
       mach_one_bnd_integral, Vpar_smoothing, vpar_smoothing_coef, no_mach1_bc, Mach1_openBC, Mach1_fix_B, &
       Number_RMP_harmonics, RMP_har_cos_spectrum,RMP_har_sin_spectrum, grid_to_wall, n_wall_blocks, keep_n0_const
use tr_module
use mpi_mod
use mod_basisfunctions
use mod_interp
use mod_integer_types
use mod_node_indices

implicit none

! --- Routine parameters
integer,                            intent(in)    :: my_id
type (type_node_list),              intent(in)    :: node_list
type (type_element_list),           intent(in)    :: element_list
type (type_bnd_node_list),          intent(in)    :: bnd_node_list
integer,                            intent(in)    :: local_elms(*)
integer,                            intent(in)    :: n_local_elms
integer,                            intent(in)    :: index_min
integer,                            intent(in)    :: index_max
logical,                            intent(in)    :: xpoint2
integer,                            intent(in)    :: xcase2
real*8,                             intent(in)    :: R_axis
real*8,                             intent(in)    :: Z_axis
real*8,                             intent(in)    :: psi_axis
real*8,                             intent(in)    :: psi_bnd
real*8,                             intent(in)    :: R_xpoint(2)
real*8,                             intent(in)    :: Z_xpoint(2)
real*8,                             intent(in)    :: psi_xpoint(2)
real*8,                             intent(inout) :: rhs_loc(*)
type(type_SP_MATRIX)                              :: a_mat

! Internal parameters
real*8  :: zbig, zbig_backup
real*8  :: R, R_b, R_s, R_t, R_bb, R_st, R_tt, R_mid, R_center
real*8  :: Z, Z_b, Z_s, Z_t, Z_bb, Z_st, Z_tt, Z_mid, Z_center
real*8  :: xjac, xjac_b
real*8  :: A3, A3_s, A3_t, A3_tt, A3_R, A3_Z, A3_p, A3_b, A3_bb, A32_b, A3_st, A3_Rb, A3_Zb
real*8  :: AR, AR_s, AR_t, AR_tt, AR_R, AR_Z, AR_p, AR_b, AR_bb, AR2_b, AR_st,        AR_Zb
real*8  :: AZ, AZ_s, AZ_t, AZ_tt, AZ_R, AZ_Z, AZ_p, AZ_b, AZ_bb, AZ2_b, AZ_st, AZ_Rb
real*8  :: uR, uR_b
real*8  :: uZ, uZ_b
real*8  :: up, up_b
real*8  :: T0, T0_b
real*8  :: normal(2), normal_direction(2), cs_direction
real*8  :: grad_s(2), grad_t(2), grad_b(2)
real*8  :: BR, BZ, Bp, BB2, Btot, Fprof, Fprof_b, B_dot_n
real*8  :: A3_R_A3b, A3_Z_A3b, BR_A3b, BZ_A3b, BB2_A3b, beta_A3b
real*8  :: A3_R_A3c, A3_Z_A3c, BR_A3c, BZ_A3c, BB2_A3c, beta_A3c
real*8  :: AR_R_ARb, AR_Z_ARb, Bp_ARb,         BB2_ARb, beta_ARb
real*8  :: AR_R_ARc, AR_Z_ARc, Bp_ARc,         BB2_ARc, beta_ARc
real*8  :: AZ_R_AZb, AZ_Z_AZb, Bp_AZb,         BB2_AZb, beta_AZb
real*8  :: AZ_R_AZc, AZ_Z_AZc, Bp_AZc,         BB2_AZc, beta_AZc
real*8  :: BR_b, BZ_b, Bp_b, BB2_b
real*8  :: A3_Rb_A3b, A3_Rb_A3c, A3_Rb_A3st
real*8  :: A3_Zb_A3b, A3_Zb_A3c, A3_Zb_A3st
real*8  :: AR_Zb_ARb, AR_Zb_ARc, AR_Zb_ARst
real*8  :: AZ_Rb_AZb, AZ_Rb_AZc, AZ_Rb_AZst
real*8  :: BR_b_A3b, BR_b_A3c, BR_b_A3st
real*8  :: BZ_b_A3b, BZ_b_A3c, BZ_b_A3st
real*8  :: Bp_b_ARb, Bp_b_ARc, Bp_b_ARst
real*8  :: Bp_b_AZb, Bp_b_AZc, Bp_b_AZst
real*8  :: BB2_b_A3b,BB2_b_A3c,BB2_b_A3st
real*8  :: BB2_b_ARb,BB2_b_ARc,BB2_b_ARst
real*8  :: BB2_b_AZb,BB2_b_AZc,BB2_b_AZst
real*8  :: beta_b_A3b, beta_b_A3c, beta_b_A3st
real*8  :: beta_b_ARb, beta_b_ARc, beta_b_ARst
real*8  :: beta_b_AZb, beta_b_AZc, beta_b_AZst 

real*8  :: BC_tmp
integer :: var_VVV(3)

real*8  :: Cs,   Cs_T,   Cs_b,   Cs_b_T,   Cs_b_Tb
real*8  :: beta, beta_T, beta_b, beta_b_T, beta_b_Tb

real*8  :: element_size_s, element_size_t, element_size_0, element_size_2
real*8  :: H1(2,n_degrees_1d), H1_s(2,n_degrees_1d), H1_ss(2,n_degrees_1d)

integer :: i, in, iv, iv2, iv3, inode, inode2, inode3, k
integer :: j, err, itest, i_mid, i_bnd, idir, iv_dir, iv_perp_dir, k_max
integer :: index_large_i, index_node, index_node2, index_node3, index_node4, ielm
integer(kind=int_all) :: ijA_position,ijA_position2
integer :: ilarge2, ilarge_vv, ilarge_vT, ilarge_vus, ilarge_vn
integer :: ilarge_vsvs, ilarge_vsTs, ilarge_vsT, ilarge_vut, ilarge_vtvt, ilarge_vtTt, ilarge_vtT
integer :: ilarge_vp, ilarge_vp2
integer :: ierr
logical :: apply_psi_BC, apply_current_BC, s_constant_boundary, t_constant_boundary, apply_cs, apply_dirichlet_all

real*8  :: factor, factor_b
real*8  :: c_1, c_2, c_3
real*8  :: bn, bn_b, bn_1, bn_2
real*8  :: dl, dl_b
real*8  :: hfact_b

integer :: n_rmp_harm, N_rmp_har_block_size
real*8, allocatable :: psi_RMP_cos1(:),dpsi_RMP_cos_dR1(:),dpsi_RMP_cos_dZ1(:)
real*8, allocatable :: psi_RMP_sin1(:),dpsi_RMP_sin_dR1(:),dpsi_RMP_sin_dZ1(:)
real*8  :: establish_RMP
real*8  :: delta_psi_rmp, delta_psi_rmp_dR, delta_psi_rmp_dZ, delta_psi_rmp_ds, delta_psi_rmp_dt, psi_test, sigmo_fonc

integer :: node_indices( (n_order+1)/2, (n_order+1)/2 ), index_tmp, kk, ll
logical, parameter :: include_2nd_derivatives = .false.


RMPspectrum: if (RMP_on .and. (n_tor .ge. 3)) then !*****
  
! for the moment it's done in a way that all RMP harmonics follow each other,i.e. n=2,n=3,n=4... 
! if you want for example n=2 and n=4 RMP you should consider n=2,3,4, but put zeros at the boundary in the input file for n=3 RMP
! example: ntor=13 and nperiod=1(so taking into account, toroidal numbers n=0,1,2....6) and  n=2 and n=3 are toroidal numbers of RMPs, 
! so Number_RMP_harmonics=2, RMP_har_cos_spectrum(1)=4,RMP_har_sin_spectrum(1)=5,RMP_har_cos_spectrum(2)=6,RMP_har_sin_spectrum(2)=7.  
  
  call tr_allocate(psi_RMP_cos1,1, bnd_node_list%n_bnd_nodes*Number_RMP_harmonics,"psi_RMP_cos1",CAT_UNKNOWN)
  call tr_allocate(dpsi_RMP_cos_dR1,1,bnd_node_list%n_bnd_nodes*Number_RMP_harmonics,"dpsi_RMP_cos_dR1",CAT_UNKNOWN)
  call tr_allocate(dpsi_RMP_cos_dZ1,1,bnd_node_list%n_bnd_nodes*Number_RMP_harmonics,"dpsi_RMP_cos_dZ1",CAT_UNKNOWN)
  call tr_allocate(psi_RMP_sin1,1,bnd_node_list%n_bnd_nodes*Number_RMP_harmonics,"psi_RMP_sin1",CAT_UNKNOWN)
  call tr_allocate(dpsi_RMP_sin_dR1,1,bnd_node_list%n_bnd_nodes*Number_RMP_harmonics,"dpsi_RMP_sin_dR1",CAT_UNKNOWN)
  call tr_allocate(dpsi_RMP_sin_dZ1,1,bnd_node_list%n_bnd_nodes*Number_RMP_harmonics,"dpsi_RMP_sin_dZ1",CAT_UNKNOWN)
  
  N_rmp_har_block_size=bnd_node_list%n_bnd_nodes
    
  psi_test =  node_list%node(bnd_node_list%bnd_node(1)%index_jorek)%values(RMP_har_cos_spectrum(1),1,1)
  
  ! if necessary, replace by:
  ! psi_test =  node_list%node(bnd_node_list%bnd_node(1)%index_jorek)%values(min(RMP_har_cos_spectrum(1), n_tor),1,1)
  write (*,*) 'psi_bnd at previous time step', psi_test
    
  if (abs(psi_test) .le. abs(psi_RMP_cos(1))) then
    sigmo_fonc = ( 1.d0 + exp(-RMP_growth_rate*( t_now - RMP_start_time - RMP_ramp_up_time/2.d0 )))**(-1) &
               - ( 1.d0 + exp(-RMP_growth_rate*( 0.d0 - RMP_ramp_up_time/2.d0 )))**(-1) 
    establish_RMP = (RMP_growth_rate*sigmo_fonc*(1-sigmo_fonc)+1.e-6)*tstep 
  else
    establish_RMP = 0.d0
  endif
  ! Other possibility (simpler) : if ( (t_now - RMP_start_time) .ge. 2.2*RMP_ramp_up_time/2.d0 ) then establish_RMP =0.0
  
  do j=1, bnd_node_list%n_bnd_nodes*Number_RMP_harmonics  
    psi_RMP_cos1(j)     = psi_RMP_cos(j)     * establish_RMP
    dpsi_RMP_cos_dR1(j) = dpsi_RMP_cos_dR(j) * establish_RMP
    dpsi_RMP_cos_dZ1(j) = dpsi_RMP_cos_dZ(j) * establish_RMP
    psi_RMP_sin1(j)     = psi_RMP_sin(j)     * establish_RMP
    dpsi_RMP_sin_dR1(j) = dpsi_RMP_sin_dR(j) * establish_RMP
    dpsi_RMP_sin_dZ1(j) = dpsi_RMP_sin_dZ(j) * establish_RMP
  end do

  if (my_id == 0) then
    write (*,*) 'psi_RMP_cos1(1) and derivatives after multiplication in boundary conditions'
    write (*,*) psi_RMP_cos1(1), dpsi_RMP_cos_dR1(1), dpsi_RMP_cos_dZ1(1)
    write (*,*) 'establish_RMP', establish_RMP
  endif

end if RMPspectrum

zbig        = 1.d12
zbig_backup = zbig

! --- calculate node_indices
call calculate_node_indices(node_indices)

do i=1, n_local_elms !=== do elements

  ielm = local_elms(i)

  i_bnd = 0

  do iv=1, n_vertex_max 
    inode = element_list%element(ielm)%vertex(iv)
    if (node_list%node(inode)%boundary .ne. 0) i_bnd = i_bnd + 1
  enddo
  
  if (i_bnd .lt. 2) cycle           

  R_mid = 0.d0; Z_mid = 0.d0; R_center = 0.d0; Z_center = 0.d0
  
  iv2 = 0; iv3 = 0

  do iv=1, n_vertex_max ! check vertices for being a boundary point

    inode = element_list%element(ielm)%vertex(iv)

    if (node_list%node(inode)%boundary .eq. 0) cycle 

    do idir=1, 2        ! check the two directions

      R_mid = node_list%node(inode)%x(1,1,1)
      Z_mid = node_list%node(inode)%x(1,1,2)

      if (idir .eq. 1) then
        iv2 = mod(iv  ,4) + 1
        iv3 = mod(iv+2,4) + 1
      else
        iv2 = mod(iv+2,4) + 1
        iv3 = mod(iv  ,4) + 1
      endif

      inode2 = element_list%element(ielm)%vertex(iv2)
      inode3 = element_list%element(ielm)%vertex(iv3)

      if (node_list%node(inode2)%boundary .eq. 0) cycle

      if ((iv*iv2 .eq. 2) .or. (iv*iv2 .eq. 12)) then
        s_constant_boundary = .false.
        t_constant_boundary = .true.
        iv_dir      = 2
        iv_perp_dir = 3
      elseif ((iv*iv2 .eq. 6) .or. (iv*iv2 .eq. 4)) then
        s_constant_boundary = .true.
        t_constant_boundary = .false.
        iv_dir      = 3
        iv_perp_dir = 2
      else
        write(*,*) 'THIS SHOULD NOT BE POSSIBLE'
      endif


      R_center = node_list%node(inode3)%x(1,1,1)
      Z_center = node_list%node(inode3)%x(1,1,2)

      normal_direction = (/R_mid - R_center, Z_mid - Z_center /) / norm2((/R_mid - R_center, Z_mid - Z_center /))

      apply_cs             = .false.
      apply_dirichlet_all  = .false.

      ! --- Dirichlet for normal nodes?

      if (                                              &
               (node_list%node(inode)%boundary .eq.  2) &
          .or. (node_list%node(inode)%boundary .eq. 12) &
      ) then
        apply_dirichlet_all = .true.
      endif

      ! --- Mach1 for normal nodes?
      if (                                              &
               (node_list%node(inode)%boundary .eq.  1) &
          .or. (node_list%node(inode)%boundary .eq.  4) &
          .or. (node_list%node(inode)%boundary .eq.  5) &
          .or. (node_list%node(inode)%boundary .eq. 11) &
          .or. (node_list%node(inode)%boundary .eq. 15) &
      ) then
        apply_cs = .true.
      endif

      ! --- Dirichlet and Mach-1 for corner nodes? UNCOMMENT ONLY ONE OF THEM!
      if (                                              &
               (node_list%node(inode)%boundary .eq.  3) &
          .or. (node_list%node(inode)%boundary .eq.  9) &
          .or. (node_list%node(inode)%boundary .eq. 19) &
          .or. (node_list%node(inode)%boundary .eq. 20) &
          .or. (node_list%node(inode)%boundary .eq. 21) &
      ) then
        ! --- UNCOMMENT ONLY ONE OF THEM!
        apply_dirichlet_all = .true.
        !apply_cs = .true.
      endif

      ! --- special user-requests to control Mach-1
      if (no_mach1_bc) apply_cs = .false.
      if (no_mach1_bc) apply_dirichlet_all = .true.
      if (Mach1_openBC) apply_cs = .false.


      do in=a_mat%i_tor_min, a_mat%i_tor_max  ! === do n_tor
      
        if (keep_n0_const  .and.  in .eq. 1 ) then
          zbig = 1.d15
        else
          zbig = zbig_backup
        endif

        ! --- start RMPs (for both directions)
        ! --- WARNING: THIS IS NOT YET ADAPTED TO FMHD!!! NEED PSI BUT ALSO AR AND AZ UPDATE
        if (RMP_on ) then

          do n_rmp_harm=1, Number_RMP_harmonics !=== do RMP harmonics

            if (((in.eq.RMP_har_cos_spectrum(n_rmp_harm)) .or. (in.eq.RMP_har_sin_spectrum(n_rmp_harm))) &
               .and. (.not. freeboundary)) then
              ! in .eq. RMP_har_cos corresponds to cos(n_perturbation)
              ! in .eq. RMP_har_sin corresponds to sin(n_perturbation)
                                 
              R   = node_list%node(inode)%x(1,1,1) 
              R_s = node_list%node(inode)%x(1,iv_dir,1) 
              Z   = node_list%node(inode)%x(1,1,2) 
              Z_s = node_list%node(inode)%x(1,iv_dir,2) 
                  
              if (in.eq.RMP_har_cos_spectrum(n_rmp_harm)) then
                delta_psi_rmp = psi_RMP_cos1(node_list%node(inode)%boundary_index +N_rmp_har_block_size*(n_rmp_harm-1))
                delta_psi_rmp_dR = dpsi_RMP_cos_dR1(node_list%node(inode)%boundary_index+N_rmp_har_block_size*(n_rmp_harm-1))
                delta_psi_rmp_dZ = dpsi_RMP_cos_dZ1(node_list%node(inode)%boundary_index+N_rmp_har_block_size*(n_rmp_harm-1))
              else 
                delta_psi_rmp = psi_RMP_sin1(node_list%node(inode)%boundary_index+N_rmp_har_block_size*(n_rmp_harm-1))
                delta_psi_rmp_dR = dpsi_RMP_sin_dR1(node_list%node(inode)%boundary_index+N_rmp_har_block_size*(n_rmp_harm-1))
                delta_psi_rmp_dZ = dpsi_RMP_sin_dZ1(node_list%node(inode)%boundary_index+N_rmp_har_block_size*(n_rmp_harm-1))
              endif
                  
              delta_psi_rmp_ds = delta_psi_rmp_dR * R_s + delta_psi_rmp_dZ * Z_s

              ! --- Apply RMP to variable value
              index_node = node_list%node(inode)%index(1) 
              call boundary_conditions_add_one_entry(                &
                     index_node, var_A3, in, index_node, var_A3, in, &
                     zbig, index_min, index_max, a_mat)
              call boundary_conditions_add_RHS(                      &
                     index_node, var_A3, in, index_min, index_max,   &
                     RHS_loc, ZBIG * delta_psi_rmp, a_mat%i_tor_min, a_mat%i_tor_max)
                  
              ! --- Apply RMP to variable derivative
              index_node2 = node_list%node(inode)%index(iv_dir)
              call boundary_conditions_add_one_entry(                 &
                     index_node2, var_A3, in, index_node2, var_A3, in,&
                     zbig, index_min, index_max, a_mat)
              call boundary_conditions_add_RHS(                       &
                     index_node2, var_A3, in, index_min, index_max,   &
                     RHS_loc, ZBIG * delta_psi_rmp_ds, a_mat%i_tor_min, a_mat%i_tor_max)

            endif !=== endif selection RMP harmonics
          enddo   !=== enddo RMP harmonics   
        endif     !=== endif RMP

 
        do k=1, n_var ! === do variables
                                                                                                 
          ! --- Decide when A3,AR,AZ need BCs
          apply_psi_BC = .false.
          if (     (k == var_A3) &
              .or. (k == var_AR) & ! will be needed eventually for RMPs
              .or. (k == var_AZ) & ! will be needed eventually for RMPs
             ) then                        
            if ( (RMP_on) .and. (in .lt. RMP_har_cos_spectrum(1))                    )   apply_psi_BC = .true.
            if ( (RMP_on) .and. (in .gt. RMP_har_sin_spectrum(Number_RMP_harmonics)) )   apply_psi_BC = .true.
            if ( (.not. RMP_on) .and. (in .ge. 2)              )                         apply_psi_BC = .true.
            if (in .eq. 1)                                                               apply_psi_BC = .true.
            if (is_freebound(in,k))                                                      apply_psi_BC = .false.                     
          endif
                
          ! --- Apply Dirichlet
          if (        apply_psi_BC      &
                 .or. ((k .eq. var_rho)  .and. apply_dirichlet_all)  &
                 .or. ((k .eq. var_T)    .and. apply_dirichlet_all)  &
                 .or. ((k .eq. var_uR)   .and. apply_dirichlet_all)  &
                 .or. ((k .eq. var_uZ)   .and. apply_dirichlet_all)  &
                 .or. ((k .eq. var_up)   .and. apply_dirichlet_all)  &
                 .or. ((k .eq. var_rhon) .and. apply_dirichlet_all)  &
              ) then

            ! --- Fix derivatives in one direction
            do kk = 1,(n_order+1)/2
              if ( (iv_dir .eq. 3) .and. (kk .gt. 1) ) cycle ! do only t-derivatives and node value
              do ll = 1,(n_order+1)/2
                if ( (iv_dir .eq. 2) .and. (ll .gt. 1) ) cycle ! do only s-derivatives and node value
                index_tmp = node_indices(kk,ll)
                index_node = node_list%node(inode)%index(index_tmp)
                call boundary_conditions_add_one_entry(                 &
                       index_node, k, in, index_node, k, in,            &
                       zbig, index_min, index_max, a_mat)
              enddo
            enddo

            ! --- Fix A3 also across boundary???
            if (k .eq. var_A3) then
              index_node = node_list%node(inode)%index(iv_perp_dir)
              call boundary_conditions_add_one_entry(                 &
                     index_node, k, in, index_node, k, in,            &
                     zbig, index_min, index_max, a_mat)
              index_node = node_list%node(inode)%index(4)
              call boundary_conditions_add_one_entry(                 &
                     index_node, k, in, index_node, k, in,            &
                     zbig, index_min, index_max, a_mat)
            endif

          endif ! --- apply Dirichlet
        enddo ! ---  variables

        ! --- Apply Mach1
        if ( (.not. mach_one_bnd_integral) .and. apply_cs) then

          ! --- Element sizes
          call basisfunctions1(0.d0, H1, H1_s, H1_ss)
          element_size_s = element_list%element(ielm)%size(iv,2) * H1_s(1,2)
          element_size_t = element_list%element(ielm)%size(iv,3) * H1_s(1,2)
          if ((iv .eq. 2) .or. (iv .eq. 3))  element_size_s = - element_size_s 
          if ((iv .eq. 3) .or. (iv .eq. 4))  element_size_t = - element_size_t 
          element_size_0 =   element_list%element(ielm)%size(iv, iv_dir) * H1_s(1,2) 
          element_size_2 =   element_list%element(ielm)%size(iv2,iv_dir) * H1_s(1,2) 

          if (t_constant_boundary .and. ((iv  .eq. 2) .or. (iv  .eq. 3))) element_size_0 = - element_size_0
          if (s_constant_boundary .and. ((iv  .eq. 3) .or. (iv  .eq. 4))) element_size_0 = - element_size_0
          if (t_constant_boundary .and. ((iv2 .eq. 2) .or. (iv2 .eq. 3))) element_size_2 = - element_size_2
          if (s_constant_boundary .and. ((iv2 .eq. 3) .or. (iv2 .eq. 4))) element_size_2 = - element_size_2

          ! --- Only use this with Vpar smoothing
          if (.not. vpar_smoothing) then
            element_size_0 = 1.d0
            element_size_s = 1.d0
            element_size_t = 1.d0
          endif
          
          ! --- Variables construction
          index_node  = node_list%node(inode)%index(1)             ! position of value
          index_node2 = node_list%node(inode)%index(iv_dir)        ! position of main deriative
          index_node3 = node_list%node(inode)%index(iv_perp_dir)   ! position of other deriative
          index_node4 = node_list%node(inode)%index(4)             ! position of cross deriative

          R        = node_list%node(inode)%x(1,1,1)
          R_b      = node_list%node(inode)%x(1,iv_dir,1) * element_size_0
          R_s      = node_list%node(inode)%x(1,2,1)      * element_size_s
          R_t      = node_list%node(inode)%x(1,3,1)      * element_size_t    
          R_st     = node_list%node(inode)%x(1,4,1)
          R_bb     = + element_list%element(ielm)%size(iv ,1)      * node_list%node(inode )%x(1,1,1)      * H1_ss(1,1)  &
                     + element_list%element(ielm)%size(iv ,iv_dir) * node_list%node(inode )%x(1,iv_dir,1) * H1_ss(1,2)  &
                     + element_list%element(ielm)%size(iv2,1)      * node_list%node(inode2)%x(1,1,1)      * H1_ss(2,1)  &
                     + element_list%element(ielm)%size(iv2,iv_dir) * node_list%node(inode2)%x(1,iv_dir,1) * H1_ss(2,2)  

          Z        = node_list%node(inode)%x(1,1,2)
          Z_b      = node_list%node(inode)%x(1,iv_dir,2) * element_size_0
          Z_s      = node_list%node(inode)%x(1,2,2)      * element_size_s
          Z_t      = node_list%node(inode)%x(1,3,2)      * element_size_t    
          Z_st     = node_list%node(inode)%x(1,4,2)
          Z_bb     = + element_list%element(ielm)%size(iv ,1)      * node_list%node(inode )%x(1,1,     2) * H1_ss(1,1)  &
                     + element_list%element(ielm)%size(iv ,iv_dir) * node_list%node(inode )%x(1,iv_dir,2) * H1_ss(1,2)  &
                     + element_list%element(ielm)%size(iv2,1)      * node_list%node(inode2)%x(1,1,     2) * H1_ss(2,1)  &
                     + element_list%element(ielm)%size(iv2,iv_dir) * node_list%node(inode2)%x(1,iv_dir,2) * H1_ss(2,2)  

          xjac     =  R_s*Z_t - R_t*Z_s

          T0       = max(node_list%node(inode)%values(1,1,var_T), T_min)
          T0_b     = node_list%node(inode)%values(1,iv_dir,var_T)    * element_size_0 

          uR       = node_list%node(inode)%values(1,1,var_uR)
          uR_b     = node_list%node(inode)%values(1,iv_dir,var_uR) * element_size_0 

          uZ       = node_list%node(inode)%values(1,1,var_uZ)
          uZ_b     = node_list%node(inode)%values(1,iv_dir,var_uZ) * element_size_0 

          up       = node_list%node(inode)%values(1,1,var_up)
          up_b     = node_list%node(inode)%values(1,iv_dir,var_up) * element_size_0 

          A3       = node_list%node(inode)%values(1,1,var_A3)
          A3_b     = node_list%node(inode)%values(1,iv_dir,var_A3)  * element_size_0 
          A3_s     = node_list%node(inode)%values(1,2,var_A3)       * element_size_s
          A3_t     = node_list%node(inode)%values(1,3,var_A3)       * element_size_t
          A3_st    = node_list%node(inode)%values(1,4,var_A3)       * element_size_s * element_size_t
          A3_bb    = element_list%element(ielm)%size(iv ,1)      * node_list%node(inode )%values(1,1,var_psi)      * H1_ss(1,1) &
                   + element_list%element(ielm)%size(iv ,iv_dir) * node_list%node(inode )%values(1,iv_dir,var_psi) * H1_ss(1,2) &
                   + element_list%element(ielm)%size(iv2,1)      * node_list%node(inode2)%values(1,1,var_psi)      * H1_ss(2,1) &
                   + element_list%element(ielm)%size(iv2,iv_dir) * node_list%node(inode2)%values(1,iv_dir,var_psi) * H1_ss(2,2)
          A3_R     = (   Z_t * A3_s - Z_s * A3_t ) / xjac
          A3_Z     = ( - R_t * A3_s + R_s * A3_t ) / xjac
          A3_p     = 0.d0 ! WARNING: ARE WE ASSUMING TOROIDAL AXISYMMETRY???
          A32_b    = node_list%node(inode2)%values(1,iv_dir,var_psi) * element_size_2 

          AR       = node_list%node(inode)%values(1,1,var_AR)
          AR_b     = node_list%node(inode)%values(1,iv_dir,var_AR)  * element_size_0 
          AR_s     = node_list%node(inode)%values(1,2,var_AR)       * element_size_s
          AR_t     = node_list%node(inode)%values(1,3,var_AR)       * element_size_t
          AR_st    = node_list%node(inode)%values(1,4,var_AR)       * element_size_s * element_size_t
          AR_bb    = element_list%element(ielm)%size(iv ,1)      * node_list%node(inode )%values(1,1,var_AR)      * H1_ss(1,1) &
                   + element_list%element(ielm)%size(iv ,iv_dir) * node_list%node(inode )%values(1,iv_dir,var_AR) * H1_ss(1,2) &
                   + element_list%element(ielm)%size(iv2,1)      * node_list%node(inode2)%values(1,1,var_AR)      * H1_ss(2,1) &
                   + element_list%element(ielm)%size(iv2,iv_dir) * node_list%node(inode2)%values(1,iv_dir,var_AR) * H1_ss(2,2)
          AR_R     = (   Z_t * AR_s - Z_s * AR_t ) / xjac
          AR_Z     = ( - R_t * AR_s + R_s * AR_t ) / xjac
          AR_p     = 0.d0 ! WARNING: ARE WE ASSUMING TOROIDAL AXISYMMETRY???

          AZ       = node_list%node(inode)%values(1,1,var_AZ)
          AZ_b     = node_list%node(inode)%values(1,iv_dir,var_AZ)  * element_size_0 
          AZ_s     = node_list%node(inode)%values(1,2,var_AZ)       * element_size_s
          AZ_t     = node_list%node(inode)%values(1,3,var_AZ)       * element_size_t
          AZ_st    = node_list%node(inode)%values(1,4,var_AZ)       * element_size_s * element_size_t
          AZ_bb    = element_list%element(ielm)%size(iv ,1)      * node_list%node(inode )%values(1,1,var_AZ)      * H1_ss(1,1) &
                   + element_list%element(ielm)%size(iv ,iv_dir) * node_list%node(inode )%values(1,iv_dir,var_AZ) * H1_ss(1,2) &
                   + element_list%element(ielm)%size(iv2,1)      * node_list%node(inode2)%values(1,1,var_AZ)      * H1_ss(2,1) &
                   + element_list%element(ielm)%size(iv2,iv_dir) * node_list%node(inode2)%values(1,iv_dir,var_AZ) * H1_ss(2,2)
          AZ_R     = (   Z_t * AZ_s - Z_s * AZ_t ) / xjac
          AZ_Z     = ( - R_t * AZ_s + R_s * AZ_t ) / xjac
          AZ_p     = 0.d0 ! WARNING: ARE WE ASSUMING TOROIDAL AXISYMMETRY???

          ! --- Magnetic field (ARE WE ASSUMING TOROIDAL AXISYMMETRY???)
          Fprof = node_list%node(inode)%Fprof_eq(1)
          Fprof_b = node_list%node(inode)%Fprof_eq(iv_dir)
          ! --- Use only the initial Bphi for numerical stability?
          if (Mach1_fix_B) then
            BR    = + A3_Z / R
            BZ    = - A3_R / R
            Bp    = Fprof / R
          else
            BR    = ( A3_Z - AZ_p )/ R
            BZ    = ( AR_p - A3_R )/ R
            Bp    = ( AZ_R - AR_Z )    + Fprof / R
          endif
          BB2   = BR**2 + BZ**2 + Bp**2
          Btot  = sqrt(BB2)

          ! --- Normal to boundary and direction
          grad_s = (/  Z_t,  -R_t /) / xjac
          grad_t = (/ -Z_s,   R_s /) / xjac
          if (s_constant_boundary) then
            grad_b = grad_s
          elseif (t_constant_boundary) then
            grad_b = grad_t
          endif
          normal     = dot_product(grad_b,normal_direction) * grad_b      ! outward pointing normal
          normal     = normal / norm2(normal)
          B_dot_n = (BR * normal(1) + BZ * normal(2))
          cs_direction = B_dot_n / abs(B_dot_n)

          ! --- Linearisation w.r.t. b-derivatives of A. Note, we also consider
          ! --- c-derivatives with index_node3 (ie. cross-derivative)
          ! --- A_R = (   Z_t * A_s - Z_s * A_t ) / xjac
          ! --- A_Z = ( - R_t * A_s + R_s * A_t ) / xjac
          if (s_constant_boundary) then
            ! --- A3                   ! --- AR                   ! --- AZ
            A3_R_A3b = - Z_s / xjac ;  AR_R_ARb = - Z_s / xjac ;  AZ_R_AZb = - Z_s / xjac
            A3_Z_A3b = + R_s / xjac ;  AR_Z_ARb = + R_s / xjac ;  AZ_Z_AZb = + R_s / xjac
            A3_R_A3c = + Z_t / xjac ;  AR_R_ARc = + Z_t / xjac ;  AZ_R_AZc = + Z_t / xjac
            A3_Z_A3c = - R_t / xjac ;  AR_Z_ARc = - R_t / xjac ;  AZ_Z_AZc = - R_t / xjac
          elseif (t_constant_boundary) then
            A3_R_A3b = + Z_t / xjac ;  AR_R_ARb = + Z_t / xjac ;  AZ_R_AZb = + Z_t / xjac
            A3_Z_A3b = - R_t / xjac ;  AR_Z_ARb = - R_t / xjac ;  AZ_Z_AZb = - R_t / xjac
            A3_R_A3c = - Z_s / xjac ;  AR_R_ARc = - Z_s / xjac ;  AZ_R_AZc = - Z_s / xjac
            A3_Z_A3c = + R_s / xjac ;  AR_Z_ARc = + R_s / xjac ;  AZ_Z_AZc = + R_s / xjac
          endif
          BR_A3b =   A3_Z_A3b / R
          BZ_A3b = - A3_R_A3b / R
          BR_A3c =   A3_Z_A3c / R
          BZ_A3c = - A3_R_A3c / R
          if (Mach1_fix_B) then
            Bp_ARb = 0.d0 ; Bp_ARc = 0.d0
            Bp_AZb = 0.d0 ; Bp_AZc = 0.d0
          else
            Bp_ARb = - AR_Z_ARb ; Bp_ARc = - AR_Z_ARc
            Bp_AZb = + AZ_R_AZb ; Bp_AZc = + AZ_Z_AZc
          endif
          BB2_A3b   = 2.0*(BR*BR_A3b + BZ*BZ_A3b) 
          BB2_A3c   = 2.0*(BR*BR_A3c + BZ*BZ_A3c)
          BB2_ARb   = 2.0*Bp*Bp_ARb
          BB2_ARc   = 2.0*Bp*Bp_ARc
          BB2_AZb   = 2.0*Bp*Bp_AZb
          BB2_AZc   = 2.0*Bp*Bp_AZc

          ! --- b-derivatives for psi (these are the most important terms for numerical stability!)
          if (s_constant_boundary) then
            xjac_b  =  R_st*Z_t + R_s*Z_bb - R_bb*Z_s - R_t*Z_st
            A3_Rb   = (   Z_t  * A3_st - Z_s  * A3_bb ) / xjac &
                     +(   Z_bb * A3_s  - Z_st * A3_t  ) / xjac &
                     -(   Z_t  * A3_s  - Z_s  * A3_t  ) / xjac**2 * xjac_b
            A3_Zb   = ( - R_t  * A3_st + R_s  * A3_bb ) / xjac &
                     +( - R_bb * A3_s  + R_st * A3_t  ) / xjac &
                     -( - R_t  * A3_s  + R_s  * A3_t  ) / xjac**2 * xjac_b
            AZ_Rb   = (   Z_t  * AZ_st - Z_s  * AZ_bb ) / xjac &
                     +(   Z_bb * AZ_s  - Z_st * AZ_t  ) / xjac &
                     -(   Z_t  * AZ_s  - Z_s  * AZ_t  ) / xjac**2 * xjac_b
            AR_Zb   = ( - R_t  * AR_st + R_s  * AR_bb ) / xjac &
                     +( - R_bb * AR_s  + R_st * AR_t  ) / xjac &
                     -( - R_t  * AR_s  + R_s  * AR_t  ) / xjac**2 * xjac_b
          elseif (t_constant_boundary) then
            xjac_b  =  R_bb*Z_t + R_s*Z_st - R_st*Z_s - R_t*Z_bb
            A3_Rb   = (   Z_t  * A3_bb - Z_s  * A3_st ) / xjac &
                     +(   Z_st * A3_s  - Z_bb * A3_t  ) / xjac &
                     -(   Z_t  * A3_s  - Z_s  * A3_t  ) / xjac**2 * xjac_b
            A3_Zb   = ( - R_t  * A3_bb + R_s  * A3_st ) / xjac &
                     +( - R_st * A3_s  + R_bb * A3_t  ) / xjac &
                     -( - R_t  * A3_s  + R_s  * A3_t  ) / xjac**2 * xjac_b
            AZ_Rb   = (   Z_t  * AZ_bb - Z_s  * AZ_st ) / xjac &
                     +(   Z_st * AZ_s  - Z_bb * AZ_t  ) / xjac &
                     -(   Z_t  * AZ_s  - Z_s  * AZ_t  ) / xjac**2 * xjac_b
            AR_Zb   = ( - R_t  * AR_bb + R_s  * AR_st ) / xjac &
                     +( - R_st * AR_s  + R_bb * AR_t  ) / xjac &
                     -( - R_t  * AR_s  + R_s  * AR_t  ) / xjac**2 * xjac_b
          endif
          BR_b  = + A3_Zb / R - A3_Z * R_b / R**2
          BZ_b  = - A3_Rb / R + A3_R * R_b / R**2
          if (Mach1_fix_B) then
            Bp_b  = + Fprof_b / R - Fprof * R_b / R**2
          else
            Bp_b  = ( AZ_Rb - AR_Zb )    + Fprof_b / R - Fprof * R_b / R**2
          endif
          BB2_b = 2.0 * (BR_b*BR + BZ_b*BZ + Bp_b*Bp)

          ! --- Linearisation of b-derivatives terms
          if (s_constant_boundary) then
            A3_Rb_A3b   = ( - Z_st ) / xjac &
                         -( - Z_s  ) / xjac**2 * xjac_b
            A3_Rb_A3c   = ( + Z_bb ) / xjac &
                         -( + Z_t  ) / xjac**2 * xjac_b
            A3_Rb_A3st  = ( + Z_t  ) / xjac 
            A3_Zb_A3b   = ( + R_st ) / xjac &
                         -( + R_s  ) / xjac**2 * xjac_b
            A3_Zb_A3c   = ( - R_bb ) / xjac &
                         -( - R_t  ) / xjac**2 * xjac_b
            A3_Zb_A3st  = ( - R_t  ) / xjac 
            AZ_Rb_AZb   = ( - Z_st ) / xjac &
                         -( - Z_s  ) / xjac**2 * xjac_b
            AZ_Rb_AZc   = (   Z_bb ) / xjac &
                         -(   Z_t  ) / xjac**2 * xjac_b
            AZ_Rb_AZst  = (   Z_t  ) / xjac 
            AR_Zb_ARb   = ( + R_st ) / xjac &
                         -( + R_s  ) / xjac**2 * xjac_b
            AR_Zb_ARc   = ( - R_bb ) / xjac &
                         -( - R_t  ) / xjac**2 * xjac_b
            AR_Zb_ARst  = ( - R_t  ) / xjac 
          elseif (t_constant_boundary) then
            A3_Rb_A3b   = (   Z_st ) / xjac &
                         -(   Z_t  ) / xjac**2 * xjac_b
            A3_Rb_A3c   = ( - Z_bb ) / xjac &
                         -( - Z_s  ) / xjac**2 * xjac_b
            A3_Rb_A3st  = ( - Z_s  ) / xjac 
            A3_Zb_A3b   = ( - R_st ) / xjac &
                         -( - R_t  ) / xjac**2 * xjac_b
            A3_Zb_A3c   = ( + R_bb ) / xjac &
                         -( + R_s  ) / xjac**2 * xjac_b
            A3_Zb_A3st  = ( + R_s  ) / xjac 
            AZ_Rb_AZb   = (   Z_st ) / xjac &
                         -(   Z_t  ) / xjac**2 * xjac_b
            AZ_Rb_AZc   = ( - Z_bb ) / xjac &
                         -( - Z_s  ) / xjac**2 * xjac_b
            AZ_Rb_AZst  = ( - Z_s  ) / xjac 
            AR_Zb_ARb   = ( - R_st ) / xjac &
                         -( - R_t  ) / xjac**2 * xjac_b
            AR_Zb_ARc   = ( + R_bb ) / xjac &
                         -( + R_s  ) / xjac**2 * xjac_b
            AR_Zb_ARst  = ( + R_s  ) / xjac 
          endif
          BR_b_A3b  = + A3_Zb_A3b / R - A3_Z_A3b * R_b / R**2
          BR_b_A3c  = + A3_Zb_A3c / R - A3_Z_A3c * R_b / R**2
          BR_b_A3st = + A3_Zb_A3st/ R
          BZ_b_A3b  = - A3_Rb_A3b / R + A3_R_A3b * R_b / R**2
          BZ_b_A3c  = - A3_Rb_A3c / R + A3_R_A3c * R_b / R**2
          BZ_b_A3st = - A3_Rb_A3st/ R
          if (Mach1_fix_B) then
            Bp_b_ARb  = 0.d0 ; Bp_b_AZb  = 0.d0
            Bp_b_ARc  = 0.d0 ; Bp_b_AZc  = 0.d0
            Bp_b_ARst = 0.d0 ; Bp_b_AZst = 0.d0
          else
            Bp_b_ARb  = ( - AR_Zb_ARb ) ; Bp_b_AZb  = ( AZ_Rb_AZb )
            Bp_b_ARc  = ( - AR_Zb_ARc ) ; Bp_b_AZc  = ( AZ_Rb_AZc )
            Bp_b_ARst = ( - AR_Zb_ARst) ; Bp_b_AZst = ( AZ_Rb_AZst)
          endif
          BB2_b_A3b  = 2.0 * (BR_b_A3b *BR + BR_b*BR_A3b  + BZ_b_A3b *BZ + BZ_b*BZ_A3b )
          BB2_b_A3c  = 2.0 * (BR_b_A3c *BR + BR_b*BR_A3c  + BZ_b_A3c *BZ + BZ_b*BZ_A3c )
          BB2_b_A3st = 2.0 * (BR_b_A3st*BR                + BZ_b_A3st*BZ               )
          BB2_b_ARb  = 2.0 * (Bp_b_ARb *Bp + Bp_b*Bp_ARb )  ;  BB2_b_AZb  = 2.0 * (Bp_b_AZb *Bp + Bp_b*Bp_AZb)
          BB2_b_ARc  = 2.0 * (Bp_b_ARc *Bp + Bp_b*Bp_ARc )  ;  BB2_b_AZc  = 2.0 * (Bp_b_AZc *Bp + Bp_b*Bp_AZc)
          BB2_b_ARst = 2.0 * (Bp_b_ARst*Bp               )  ;  BB2_b_AZst = 2.0 * (Bp_b_AZst*Bp              )

          ! --- Important: we assume that B is ~constant.
          Cs = (gamma * (T0))**0.5
          beta = Cs * cs_direction / sqrt(BB2)

          ! --- b-derivatives and c-derivatives of A
          beta_A3b = - 0.5 * Cs * cs_direction * BB2**(-1.5) * BB2_A3b
          beta_A3c = - 0.5 * Cs * cs_direction * BB2**(-1.5) * BB2_A3c
          beta_ARb = - 0.5 * Cs * cs_direction * BB2**(-1.5) * BB2_ARb
          beta_ARc = - 0.5 * Cs * cs_direction * BB2**(-1.5) * BB2_ARc
          beta_AZb = - 0.5 * Cs * cs_direction * BB2**(-1.5) * BB2_AZb
          beta_AZc = - 0.5 * Cs * cs_direction * BB2**(-1.5) * BB2_AZc

          ! --- T derivative
          Cs_T   = 0.5 * gamma * (gamma * (T0))**(-0.5)
          beta_T = Cs_T * cs_direction / sqrt(BB2)

          ! --- Vector along element side
          Cs_b = 0.5 * gamma * (T0_b) * (gamma * (T0))**(-0.5)
          beta_b = Cs_b * cs_direction / sqrt(BB2) - 0.5 * Cs * cs_direction * BB2_b * BB2**(-1.5)

          ! --- Vector along element side, T derivative
          Cs_b_T  = - 0.25 * gamma**2 * (T0_b) * (gamma * (T0))**(-1.5)
          Cs_b_Tb = 0.5 * gamma * (gamma * (T0))**(-0.5)
          beta_b_T  = Cs_b_T  * cs_direction / sqrt(BB2) - 0.5 * Cs_T * cs_direction * BB2_b * BB2**(-1.5)
          beta_b_Tb = Cs_b_Tb * cs_direction / sqrt(BB2)

          ! --- Vector along element side, A3 derivatives
          beta_b_A3b  = - 0.5 * Cs_b * cs_direction * BB2_A3b   * BB2**(-1.5) &
                        - 0.5 * Cs   * cs_direction * BB2_b_A3b * BB2**(-1.5) &
                        + 0.75* Cs   * cs_direction * BB2_b     * BB2**(-2.5) * BB2_A3b
          beta_b_A3c  = - 0.5 * Cs_b * cs_direction * BB2_A3c   * BB2**(-1.5) &
                        - 0.5 * Cs   * cs_direction * BB2_b_A3c * BB2**(-1.5) &
                        + 0.75* Cs   * cs_direction * BB2_b     * BB2**(-2.5) * BB2_A3c
          beta_b_A3st = - 0.5 * Cs   * cs_direction * BB2_b_A3st* BB2**(-1.5) 

          ! --- Vector along element side, AR derivatives
          beta_b_ARb  = - 0.5 * Cs_b * cs_direction * BB2_ARb   * BB2**(-1.5) &
                        - 0.5 * Cs   * cs_direction * BB2_b_ARb * BB2**(-1.5) &
                        + 0.75* Cs   * cs_direction * BB2_b     * BB2**(-2.5) * BB2_ARb
          beta_b_ARc  = - 0.5 * Cs_b * cs_direction * BB2_ARc   * BB2**(-1.5) &
                        - 0.5 * Cs   * cs_direction * BB2_b_ARc * BB2**(-1.5) &
                        + 0.75* Cs   * cs_direction * BB2_b     * BB2**(-2.5) * BB2_ARc
          beta_b_ARst = - 0.5 * Cs   * cs_direction * BB2_b_ARst* BB2**(-1.5) 

          ! --- Vector along element side, AZ derivatives
          beta_b_AZb  = - 0.5 * Cs_b * cs_direction * BB2_AZb   * BB2**(-1.5) &
                        - 0.5 * Cs   * cs_direction * BB2_b_AZb * BB2**(-1.5) &
                        + 0.75* Cs   * cs_direction * BB2_b     * BB2**(-2.5) * BB2_AZb
          beta_b_AZc  = - 0.5 * Cs_b * cs_direction * BB2_AZc   * BB2**(-1.5) &
                        - 0.5 * Cs   * cs_direction * BB2_b_AZc * BB2**(-1.5) &
                        + 0.75* Cs   * cs_direction * BB2_b     * BB2**(-2.5) * BB2_AZc
          beta_b_AZst = - 0.5 * Cs   * cs_direction * BB2_b_AZst* BB2**(-1.5) 

          ! --- Curve along boundary
          dl       = sqrt(R_b**2 + Z_b**2)
          dl_b     = (R_b*R_bb + Z_b*Z_bb) / dl

          ! --- B.dot.n
          bn     = dot_product( (/BR,BZ/), normal ) /  (R*Btot)
          bn_b   = 1.d0 / (Btot*dl*R) * (A3_bb - A3_b * dl_b /dl )

          bn_1 = + A3_b /(R*Btot*dl)
          bn_2 = + A32_b/(R*Btot*dl)

          if ((s_constant_boundary) .and. ((iv  .eq. 1) .or. (iv  .eq. 4))) bn_1 = - bn_1
          if ((t_constant_boundary) .and. ((iv  .eq. 3) .or. (iv  .eq. 4))) bn_1 = - bn_1
          if ((s_constant_boundary) .and. ((iv2 .eq. 1) .or. (iv2 .eq. 4))) bn_2 = - bn_2
          if ((t_constant_boundary) .and. ((iv2 .eq. 3) .or. (iv2 .eq. 4))) bn_2 = - bn_2

          ! --- Smoothing parameters
          c_1 = vpar_smoothing_coef(1); c_2 = vpar_smoothing_coef(2); c_3 = vpar_smoothing_coef(3)
          if ((vpar_smoothing) .and. (bn_1*bn_2 .lt. 0.d0)) then
            if (c_2 .gt. 0d0) then
              factor    = 0.25d0 * ( 1.d0 + tanh( (abs(bn) - c_1) / c_2 ) )**2 - c_3
              factor_b  = 0.5d0  * ( 1.d0 + tanh( (abs(bn) - c_1) / c_2 ) )           & 
                        * (bn_b * bn/abs(bn) /c_2) /(cosh( (abs(bn) - c_1) / c_2 ) )**2
             else
               factor    = tanh(bn/c_1)
               factor_b  = bn_b /c_1 / cosh(bn/c_1)**2
               cs_direction = 1.d0                            
            endif                       
            Hfact_b   = factor * R_b / R  + factor_b
          else
            factor   = 1.d0
            factor_b = 0.d0
            Hfact_b   = 0.d0
          endif

          ! --- Set arrays for uR,uZ,uP to avoid duplications
          var_VVV = (/var_uR, var_uZ, var_up/)
          do k=1,3
 
            ! -------------------------------------
            ! --- Variable Value ------------------
            ! -------------------------------------

            ! --- RHS for variable value
            if (in .eq. 1) then
              if (var_VVV(k) == var_uR) BC_tmp = uR    - factor * beta       * BR
              if (var_VVV(k) == var_uZ) BC_tmp = uZ    - factor * beta       * BZ
              if (var_VVV(k) == var_up) BC_tmp = up    - factor * beta       * Bp
              call boundary_conditions_add_RHS( index_node, var_VVV(k), in, index_min, index_max, RHS_loc, &
                                                - Zbig * BC_tmp, a_mat%i_tor_min, a_mat%i_tor_max)
            else
              call boundary_conditions_add_RHS( index_node, var_VVV(k), in, index_min, index_max, RHS_loc, &
                                                0.d0, a_mat%i_tor_min, a_mat%i_tor_max)
            endif
            
            ! --- LHS for variable value d/dVVV
            call boundary_conditions_add_one_entry( index_node, var_VVV(k), in, index_node, var_VVV(k), in, &
                                                    zbig,                                                   &
                                                    index_min, index_max, a_mat)
            
            ! --- LHS for variable value d/dT
            if (var_VVV(k) == var_uR) BC_tmp =       - factor * beta_T    * BR
            if (var_VVV(k) == var_uZ) BC_tmp =       - factor * beta_T    * BZ
            if (var_VVV(k) == var_up) BC_tmp =       - factor * beta_T    * Bp
            call boundary_conditions_add_one_entry( index_node, var_VVV(k), in, index_node, var_T, in, &
                                                    zbig * BC_tmp,                                      &
                                                    index_min, index_max, a_mat)

            ! --- LHS for variable value d/dA3_b
            if (var_VVV(k) == var_uR) BC_tmp =       - factor * beta_A3b    * BR - factor * beta    * BR_A3b
            if (var_VVV(k) == var_uZ) BC_tmp =       - factor * beta_A3b    * BZ - factor * beta    * BZ_A3b
            if (var_VVV(k) == var_up) BC_tmp =       - factor * beta_A3b    * Bp
            call boundary_conditions_add_one_entry( index_node, var_VVV(k), in, index_node2, var_A3, in, &
                                                    zbig * BC_tmp * element_size_0,                      &
                                                    index_min, index_max, a_mat)

            ! --- LHS for variable value d/dA3_c
            if (var_VVV(k) == var_uR) BC_tmp =       - factor * beta_A3c    * BR - factor * beta    * BR_A3c
            if (var_VVV(k) == var_uZ) BC_tmp =       - factor * beta_A3c    * BZ - factor * beta    * BZ_A3c
            if (var_VVV(k) == var_up) BC_tmp =       - factor * beta_A3c    * Bp
            call boundary_conditions_add_one_entry( index_node, var_VVV(k), in, index_node3, var_A3, in, &
                                                    zbig * BC_tmp * element_size_0,                      &
                                                    index_min, index_max, a_mat)

            ! --- LHS for variable value d/dAR_b
            if (var_VVV(k) == var_uR) BC_tmp =       - factor * beta_ARb    * BR
            if (var_VVV(k) == var_uZ) BC_tmp =       - factor * beta_ARb    * BZ
            if (var_VVV(k) == var_up) BC_tmp =       - factor * beta_ARb    * Bp - factor * beta    * Bp_ARb
            call boundary_conditions_add_one_entry( index_node, var_VVV(k), in, index_node2, var_AR, in, &
                                                    zbig * BC_tmp * element_size_0,                      &
                                                    index_min, index_max, a_mat)

            ! --- LHS for variable value d/dAR_c
            if (var_VVV(k) == var_uR) BC_tmp =       - factor * beta_ARc    * BR
            if (var_VVV(k) == var_uZ) BC_tmp =       - factor * beta_ARc    * BZ
            if (var_VVV(k) == var_up) BC_tmp =       - factor * beta_ARc    * Bp - factor * beta    * Bp_ARc
            call boundary_conditions_add_one_entry( index_node, var_VVV(k), in, index_node3, var_AR, in, &
                                                    zbig * BC_tmp * element_size_0,                      &
                                                    index_min, index_max, a_mat)

            ! --- LHS for variable value d/dAZ_b
            if (var_VVV(k) == var_uR) BC_tmp =       - factor * beta_AZb    * BR
            if (var_VVV(k) == var_uZ) BC_tmp =       - factor * beta_AZb    * BZ
            if (var_VVV(k) == var_up) BC_tmp =       - factor * beta_AZb    * Bp - factor * beta    * Bp_AZb
            call boundary_conditions_add_one_entry( index_node, var_VVV(k), in, index_node2, var_AZ, in, &
                                                    zbig * BC_tmp * element_size_0,                      &
                                                    index_min, index_max, a_mat)

            ! --- LHS for variable value d/dAZ_c
            if (var_VVV(k) == var_uR) BC_tmp =       - factor * beta_AZc    * BR
            if (var_VVV(k) == var_uZ) BC_tmp =       - factor * beta_AZc    * BZ
            if (var_VVV(k) == var_up) BC_tmp =       - factor * beta_AZc    * Bp - factor * beta    * Bp_AZc
            call boundary_conditions_add_one_entry( index_node, var_VVV(k), in, index_node3, var_AZ, in, &
                                                    zbig * BC_tmp * element_size_0,                      &
                                                    index_min, index_max, a_mat)
                        
            ! -------------------------------------
            ! --- Variable Derivative -------------
            ! -------------------------------------

            ! --- RHS for variable derivative
            if (in .eq. 1) then
              if (var_VVV(k) == var_uR) BC_tmp = uR_b - factor * beta_b     * BR  - Hfact_b * beta       * BR  - factor * beta     * BR_b
              if (var_VVV(k) == var_uZ) BC_tmp = uZ_b - factor * beta_b     * BZ  - Hfact_b * beta       * BZ - factor * beta     * BZ_b
              if (var_VVV(k) == var_up) BC_tmp = up_b - factor * beta_b     * Bp  - Hfact_b * beta       * Bp - factor * beta     * Bp_b
              call boundary_conditions_add_RHS( index_node2, var_VVV(k), in, index_min, index_max, RHS_loc, &
                                                - Zbig * BC_tmp, a_mat%i_tor_min, a_mat%i_tor_max)
            else
              call boundary_conditions_add_RHS( index_node2, var_VVV(k), in, index_min, index_max, RHS_loc, &
                                                0.d0, a_mat%i_tor_min, a_mat%i_tor_max) 
            endif
            
            ! --- LHS for variable derivative d/dVVV_b
            call boundary_conditions_add_one_entry( index_node2, var_VVV(k), in, index_node2, var_VVV(k), in, &
                                                    zbig * element_size_0,                                    &
                                                    index_min, index_max, a_mat)
            
            ! --- LHS for variable derivative d/dT_b
            if (var_VVV(k) == var_uR) BC_tmp =       - factor * beta_b_Tb * BR
            if (var_VVV(k) == var_uZ) BC_tmp =       - factor * beta_b_Tb * BZ
            if (var_VVV(k) == var_up) BC_tmp =       - factor * beta_b_Tb * Bp
            call boundary_conditions_add_one_entry( index_node2, var_VVV(k), in, index_node2, var_T, in, &
                                                    zbig * BC_tmp * element_size_0,                       &
                                                    index_min, index_max, a_mat)

            ! --- LHS for variable derivative d/dT
            if (var_VVV(k) == var_uR) BC_tmp =       - factor * beta_b_T  * BR  - Hfact_b * beta_T    * BR 
            if (var_VVV(k) == var_uZ) BC_tmp =       - factor * beta_b_T  * BZ  - Hfact_b * beta_T    * BZ
            if (var_VVV(k) == var_up) BC_tmp =       - factor * beta_b_T  * Bp  - Hfact_b * beta_T    * Bp
            call boundary_conditions_add_one_entry( index_node2, var_VVV(k), in, index_node, var_T, in, &
                                                    zbig * BC_tmp * element_size_0,                      &
                                                    index_min, index_max, a_mat)

            ! --- LHS for variable derivative d/dA3b
            if (var_VVV(k) == var_uR) BC_tmp = - factor * beta_b_A3b * BR     - Hfact_b * beta_A3b   * BR     - factor * beta_A3b * BR_b     &
                                               - factor * beta_b     * BR_A3b - Hfact_b * beta       * BR_A3b - factor * beta     * BR_b_A3b 
            if (var_VVV(k) == var_uZ) BC_tmp = - factor * beta_b_A3b * BZ     - Hfact_b * beta_A3b   * BZ     - factor * beta_A3b * BZ_b     &
                                               - factor * beta_b     * BZ_A3b - Hfact_b * beta       * BZ_A3b - factor * beta     * BZ_b_A3b 
            if (var_VVV(k) == var_up) BC_tmp = - factor * beta_b_A3b * Bp     - Hfact_b * beta_A3b   * Bp     - factor * beta_A3b * Bp_b     
            call boundary_conditions_add_one_entry( index_node2, var_VVV(k), in, index_node2, var_A3, in,&
                                                    zbig * BC_tmp * element_size_0,                      &
                                                    index_min, index_max, a_mat)

            ! --- LHS for variable derivative d/dA3c
            if (var_VVV(k) == var_uR) BC_tmp = - factor * beta_b_A3c * BR     - Hfact_b * beta_A3c   * BR     - factor * beta_A3c * BR_b     &
                                               - factor * beta_b     * BR_A3c - Hfact_b * beta       * BR_A3c - factor * beta     * BR_b_A3c 
            if (var_VVV(k) == var_uZ) BC_tmp = - factor * beta_b_A3c * BZ     - Hfact_b * beta_A3c   * BZ     - factor * beta_A3c * BZ_b     &
                                               - factor * beta_b     * BZ_A3c - Hfact_b * beta       * BZ_A3c - factor * beta     * BZ_b_A3c 
            if (var_VVV(k) == var_up) BC_tmp = - factor * beta_b_A3c * Bp     - Hfact_b * beta_A3c   * Bp     - factor * beta_A3c * Bp_b     
            call boundary_conditions_add_one_entry( index_node2, var_VVV(k), in, index_node3, var_A3, in,&
                                                    zbig * BC_tmp * element_size_0,                      &
                                                    index_min, index_max, a_mat)

            ! --- LHS for variable derivative d/dA3st
            if (var_VVV(k) == var_uR) BC_tmp = - factor * beta_b_A3st* BR                                                                    &
                                                                                                              - factor * beta     * BR_b_A3st
            if (var_VVV(k) == var_uZ) BC_tmp = - factor * beta_b_A3st* BZ                                                                    &
                                                                                                              - factor * beta     * BZ_b_A3st
            if (var_VVV(k) == var_up) BC_tmp = - factor * beta_b_A3st* Bp 
            call boundary_conditions_add_one_entry( index_node2, var_VVV(k), in, index_node4, var_A3, in,&
                                                    zbig * BC_tmp * element_size_0,                      &
                                                    index_min, index_max, a_mat)

            ! --- LHS for variable derivative d/dARb
            if (var_VVV(k) == var_uR) BC_tmp = - factor * beta_b_ARb * BR     - Hfact_b * beta_ARb   * BR     - factor * beta_ARb * BR_b     
            if (var_VVV(k) == var_uZ) BC_tmp = - factor * beta_b_ARb * BZ     - Hfact_b * beta_ARb   * BZ     - factor * beta_ARb * BZ_b     
            if (var_VVV(k) == var_up) BC_tmp = - factor * beta_b_ARb * Bp     - Hfact_b * beta_ARb   * Bp     - factor * beta_ARb * Bp_b     &
                                               - factor * beta_b     * Bp_ARb - Hfact_b * beta       * Bp_ARb - factor * beta     * Bp_b_ARb 
            call boundary_conditions_add_one_entry( index_node2, var_VVV(k), in, index_node2, var_AR, in,&
                                                    zbig * BC_tmp * element_size_0,                      &
                                                    index_min, index_max, a_mat)

            ! --- LHS for variable derivative d/dARc
            if (var_VVV(k) == var_uR) BC_tmp = - factor * beta_b_ARc * BR     - Hfact_b * beta_ARc   * BR     - factor * beta_ARc * BR_b     
            if (var_VVV(k) == var_uZ) BC_tmp = - factor * beta_b_ARc * BZ     - Hfact_b * beta_ARc   * BZ     - factor * beta_ARc * BZ_b     
            if (var_VVV(k) == var_up) BC_tmp = - factor * beta_b_ARc * Bp     - Hfact_b * beta_ARc   * Bp     - factor * beta_ARc * Bp_b     &
                                               - factor * beta_b     * Bp_ARc - Hfact_b * beta       * Bp_ARc - factor * beta     * Bp_b_ARc 
            call boundary_conditions_add_one_entry( index_node2, var_VVV(k), in, index_node3, var_AR, in,&
                                                    zbig * BC_tmp * element_size_0,                      &
                                                    index_min, index_max, a_mat)

            ! --- LHS for variable derivative d/dARst
            if (var_VVV(k) == var_uR) BC_tmp = - factor * beta_b_ARst* BR 
            if (var_VVV(k) == var_uZ) BC_tmp = - factor * beta_b_ARst* BZ
            if (var_VVV(k) == var_up) BC_tmp = - factor * beta_b_ARst* Bp                                                                    &
                                                                                                              - factor * beta     * Bp_b_ARst
            call boundary_conditions_add_one_entry( index_node2, var_VVV(k), in, index_node4, var_AR, in,&
                                                    zbig * BC_tmp * element_size_0,                      &
                                                    index_min, index_max, a_mat)

            ! --- LHS for variable derivative d/dAZb
            if (var_VVV(k) == var_uR) BC_tmp = - factor * beta_b_AZb * BR     - Hfact_b * beta_AZb   * BR     - factor * beta_AZb * BR_b     
            if (var_VVV(k) == var_uZ) BC_tmp = - factor * beta_b_AZb * BZ     - Hfact_b * beta_AZb   * BZ     - factor * beta_AZb * BZ_b     
            if (var_VVV(k) == var_up) BC_tmp = - factor * beta_b_AZb * Bp     - Hfact_b * beta_AZb   * Bp     - factor * beta_AZb * Bp_b     &
                                               - factor * beta_b     * Bp_AZb - Hfact_b * beta       * Bp_AZb - factor * beta     * Bp_b_AZb 
            call boundary_conditions_add_one_entry( index_node2, var_VVV(k), in, index_node2, var_AZ, in,&
                                                    zbig * BC_tmp * element_size_0,                      &
                                                    index_min, index_max, a_mat)

            ! --- LHS for variable derivative d/dAZc
            if (var_VVV(k) == var_uR) BC_tmp = - factor * beta_b_AZc * BR     - Hfact_b * beta_AZc   * BR     - factor * beta_AZc * BR_b     
            if (var_VVV(k) == var_uZ) BC_tmp = - factor * beta_b_AZc * BZ     - Hfact_b * beta_AZc   * BZ     - factor * beta_AZc * BZ_b     
            if (var_VVV(k) == var_up) BC_tmp = - factor * beta_b_AZc * Bp     - Hfact_b * beta_AZc   * Bp     - factor * beta_AZc * Bp_b     &
                                               - factor * beta_b     * Bp_AZc - Hfact_b * beta       * Bp_AZc - factor * beta     * Bp_b_AZc 
            call boundary_conditions_add_one_entry( index_node2, var_VVV(k), in, index_node3, var_AZ, in,&
                                                    zbig * BC_tmp * element_size_0,                      &
                                                    index_min, index_max, a_mat)

            ! --- LHS for variable derivative d/dAZst
            if (var_VVV(k) == var_uR) BC_tmp = - factor * beta_b_AZst* BR 
            if (var_VVV(k) == var_uZ) BC_tmp = - factor * beta_b_AZst* BZ
            if (var_VVV(k) == var_up) BC_tmp = - factor * beta_b_AZst* Bp                                                                    &
                                                                                                              - factor * beta     * Bp_b_AZst
            call boundary_conditions_add_one_entry( index_node2, var_VVV(k), in, index_node4, var_AZ, in,&
                                                    zbig * BC_tmp * element_size_0,                      &
                                                    index_min, index_max, a_mat)

            ! --- Fix derivatives in one direction
            do kk = 1,(n_order+1)/2
              do ll = 1,(n_order+1)/2
                if ( (iv_dir .eq. 2) .and. (ll .gt. 1) ) cycle ! do only pure s derivatives, not cross _st
                if ( (iv_dir .eq. 3) .and. (kk .gt. 1) ) cycle ! do only pure t derivatives, not cross _st
                if ( (iv_dir .eq. 2) .and. (kk .lt. 3) ) cycle ! do only node value, 1st and 2nd derivatives, fix the rest
                if ( (iv_dir .eq. 3) .and. (ll .lt. 3) ) cycle ! do only node value, 1st and 2nd derivatives, fix the rest
                index_tmp = node_indices(kk,ll)
                index_node = node_list%node(inode)%index(index_tmp)
                call boundary_conditions_add_one_entry(                        &
                       index_node, var_VVV(k), in, index_node, var_VVV(k), in, &
                       zbig, index_min, index_max, a_mat)
              enddo
            enddo

          enddo !=== var_VVV (3 variables of uR, uZ, up)

        endif   !=== apply_cs
        
      enddo     !=== enddo loop n_tor
    
    enddo       !=== enddo directions (i_dir)

  enddo         !=== enddo vertex
 
enddo           !=== do elements

if (RMP_on) then
  if (allocated(psi_RMP_cos1))         call tr_deallocate(psi_RMP_cos1,"psi_RMP_cos1",CAT_UNKNOWN)
  if (allocated(dpsi_RMP_cos_dR1))     call tr_deallocate(dpsi_RMP_cos_dR1,"dpsi_RMP_cos_dR1",CAT_UNKNOWN)
  if (allocated(dpsi_RMP_cos_dZ1))     call tr_deallocate(dpsi_RMP_cos_dZ1,"dpsi_RMP_cos_dZ1",CAT_UNKNOWN)
  if (allocated(psi_RMP_sin1))         call tr_deallocate(psi_RMP_sin1,"psi_RMP_sin1",CAT_UNKNOWN)
  if (allocated(dpsi_RMP_sin_dR1))     call tr_deallocate(dpsi_RMP_sin_dR1,"dpsi_RMP_sin_dR1",CAT_UNKNOWN)
  if (allocated(dpsi_RMP_sin_dZ1))     call tr_deallocate(dpsi_RMP_sin_dZ1,"dpsi_RMP_sin_dZ1",CAT_UNKNOWN)
endif

return
end subroutine boundary_conditions 

end module mod_boundary_conditions

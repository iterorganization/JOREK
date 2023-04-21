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
       mach_one_bnd_integral, Vpar_smoothing, vpar_smoothing_coef, no_mach1_bc,                            &
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
real*8  :: zbig, zbig_backup,  Ti0, Te0, Vpar0, bigR
real*8  :: R_s, R_t, Z, Z_s, Z_t, R_tt, Z_tt, ps0, ps0_s, ps0_t, ps0_tt, ps0_x, ps0_y, direction, xjac
real*8  :: ps0_b, Ti0_b, Te0_b, u0_b, Vpar0_b, Vpar0_bb, Ti0_bb, Te0_bb, u0_bb, R_b, Z_b, R_bb, Z_bb, ps0_bb, grad_b(2)
real*8  :: Btot, grad_psi, u0_s, u0_t, u0_x, u0_y
real*8  :: element_size_s, element_size_t, element_size_0, element_size_3
real*8  :: H1(2,n_degrees_1d), H1_s(2,n_degrees_1d), H1_ss(2,n_degrees_1d)
integer :: i, in, iv, iv2, iv3, inode, inode2, inode3, k
integer :: index_large_i, index_node, index_node2, index_node3, ielm
integer(kind=int_all) :: ijA_position,ijA_position2
integer :: ilarge2, kv, kTi, kTe, ku, kn, ilarge_vv, ilarge_vT, ilarge_vus, ilarge_vn
integer :: ilarge_vsvs, ilarge_vsTs, ilarge_vsT, ilarge_vut, ilarge_vtvt, ilarge_vtTt, ilarge_vtT
integer :: ierr
logical :: apply_psi_BC, apply_current_BC, s_constant_boundary, t_constant_boundary, apply_cs, apply_dirichlet_1234, apply_dirichlet_all

real*8, allocatable :: psi_RMP_cos1(:),dpsi_RMP_cos_dR1(:),dpsi_RMP_cos_dZ1(:)
real*8, allocatable :: psi_RMP_sin1(:),dpsi_RMP_sin_dR1(:),dpsi_RMP_sin_dZ1(:)
real*8  :: Rnode, dRnode_ds, Znode, dZnode_ds, dRnode_dt, dZnode_dt, establish_RMP
real*8  :: delta_psi_rmp, delta_psi_rmp_dR, delta_psi_rmp_dZ, delta_psi_rmp_ds, delta_psi_rmp_dt, psi_test, sigmo_fonc
real*8  :: R_mid, Z_mid, R_center, Z_center, direction2, normal(2), normal_direction(2), grad_s(2), grad_t(2)
real*8  :: factor, factor_b, factor_bb, c_1, c_2, c_3, bn, dl, dl_b
real*8  :: cs0, cs0_T, cs0_TT, cs0_TTT
real*8  :: bn_b, bn_b_abs, hfact_b, hfact_bb, bn_1, bn_2, ps2_b, element_size_2
integer :: ilarge_vp, ilarge_vp2
integer :: kp, j, err, itest, i_mid, i_bnd, idir, iv_dir, iv_perp_dir, k_max
integer :: n_rmp_harm, N_rmp_har_block_size

real*8  :: R_out, Z_out, s_elm, t_elm, QR,QR_s,QR_t,QR_st,QR_ss,QR_tt,QZ,QZ_s,QZ_t,QZ_st,QZ_ss,QZ_tt
real*8  :: QPs0,QPs0_s,QPs0_t,QPs0_st,QPs0_ss,QPs0_tt
integer :: ifail, i_elm

real*8  ::   Mach1BC,   Mach1BC_v,   Mach1BC_T,   Mach1BC_u
real*8  ::  dMach1BC,  dMach1BC_v,  dMach1BC_Ti, dMach1BC_Te,  dMach1BC_Tb, dMach1BC_ubb
real*8  :: d2Mach1BC, d2Mach1BC_v, d2Mach1BC_T, d2Mach1BC_Tb, d2Mach1BC_Tbb

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
      apply_dirichlet_1234 = .true.
      apply_dirichlet_all  = .false.

      if (     (node_list%node(inode)%boundary .eq.  2) &
          .or. (node_list%node(inode)%boundary .eq.  3) &
          .or. (node_list%node(inode)%boundary .eq. 12) &
          .or. (node_list%node(inode)%boundary .eq. 20) &
          .or. (node_list%node(inode)%boundary .eq. 21)) &
      then
        apply_dirichlet_all = .true.
      endif

      if      ((node_list%node(inode)%boundary .eq.  1) &
          .or. (node_list%node(inode)%boundary .eq.  3) &
          .or. (node_list%node(inode)%boundary .eq.  4) &
          .or. (node_list%node(inode)%boundary .eq.  5) &
          .or. (node_list%node(inode)%boundary .eq.  9) &
          .or. (node_list%node(inode)%boundary .eq. 11) &
          .or. (node_list%node(inode)%boundary .eq. 15) &
          .or. (node_list%node(inode)%boundary .eq. 19)) &
      then
        apply_cs = .true.
      endif
      if (no_mach1_bc) apply_cs = .false.
      if (no_mach1_bc) apply_dirichlet_all = .true.


      do in=a_mat%i_tor_min, a_mat%i_tor_max  ! === do n_tor
      
        if (keep_n0_const  .and.  in .eq. 1 ) then
          zbig = 1.d15
        else
          zbig = zbig_backup
        endif

!================================= start RMPs (for both directions) ==================================
        if (RMP_on ) then

          do n_rmp_harm=1, Number_RMP_harmonics !=== do RMP harmonics

            if (((in.eq.RMP_har_cos_spectrum(n_rmp_harm)) .or. (in.eq.RMP_har_sin_spectrum(n_rmp_harm))) &
               .and. (.not. freeboundary)) then
                   ! in .eq. RMP_har_cos corresponds to cos(n_perturbation)
                   ! in .eq. RMP_har_sin corresponds to sin(n_perturbation)
                                 
              kp=var_psi    ! variable psi
              kv=var_psi    ! equation for psi
                   
              index_node = node_list%node(inode)%index(1)  !=== index in RHS (or matrix A not compressed)
                                          
              Rnode     = node_list%node(inode)%x(1,1,1) 
              dRnode_ds = node_list%node(inode)%x(1,iv_dir,1) 
              Znode     = node_list%node(inode)%x(1,1,2) 
              dZnode_ds = node_list%node(inode)%x(1,iv_dir,2) 
                  
              if (in.eq.RMP_har_cos_spectrum(n_rmp_harm)) then
                delta_psi_rmp = psi_RMP_cos1(node_list%node(inode)%boundary_index +N_rmp_har_block_size*(n_rmp_harm-1))
                delta_psi_rmp_dR = dpsi_RMP_cos_dR1(node_list%node(inode)%boundary_index+N_rmp_har_block_size*(n_rmp_harm-1))
                delta_psi_rmp_dZ = dpsi_RMP_cos_dZ1(node_list%node(inode)%boundary_index+N_rmp_har_block_size*(n_rmp_harm-1))
              else 
                delta_psi_rmp = psi_RMP_sin1(node_list%node(inode)%boundary_index+N_rmp_har_block_size*(n_rmp_harm-1))
                delta_psi_rmp_dR = dpsi_RMP_sin_dR1(node_list%node(inode)%boundary_index+N_rmp_har_block_size*(n_rmp_harm-1))
                delta_psi_rmp_dZ = dpsi_RMP_sin_dZ1(node_list%node(inode)%boundary_index+N_rmp_har_block_size*(n_rmp_harm-1))
              endif
                  
              delta_psi_rmp_ds = delta_psi_rmp_dR * dRnode_ds + delta_psi_rmp_dZ * dZnode_ds

              call boundary_conditions_add_one_entry(                &
                     index_node, kv, in, index_node, kp, in,         &
                     zbig, index_min, index_max, a_mat)

              call boundary_conditions_add_RHS(                      &
                     index_node, kv, in, index_min, index_max,       &
                     RHS_loc, ZBIG * delta_psi_rmp, a_mat%i_tor_min, a_mat%i_tor_max)
                  
              index_node2 = node_list%node(inode)%index(iv_dir)

              call boundary_conditions_add_one_entry(                 &
                     index_node2, kv, in, index_node2, kp, in,        &
                     zbig, index_min, index_max, a_mat)

              call boundary_conditions_add_RHS(                       &
                     index_node2, kv, in, index_min, index_max,       &
                     RHS_loc, ZBIG * delta_psi_rmp_ds, a_mat%i_tor_min, a_mat%i_tor_max)

            endif !=== endif selection RMP harmonics
        
          enddo   !=== enddo RMP harmonics   
        
        endif     !=== endif RMP

 
        do k=1, n_var ! === do variables
                                                                                                 
          !------------ Decide when Psi or Current need BCs --------------------------------------------------                      
          !----Psi
          apply_psi_BC = .false.
          if (k == var_psi) then                        
            if ( (RMP_on) .and. (in .lt. RMP_har_cos_spectrum(1))                    )   apply_psi_BC = .true.
            if ( (RMP_on) .and. (in .gt. RMP_har_sin_spectrum(Number_RMP_harmonics)) )   apply_psi_BC = .true.
            if ( (.not. RMP_on) .and. (in .ge. 2)              )                         apply_psi_BC = .true.
            if (in .eq. 1)                                                               apply_psi_BC = .true.
            if (is_freebound(in,k))                                                      apply_psi_BC = .false.                     
          endif
                
          !----Current
          apply_current_BC = .false.
          if (k == var_zj) then
            if ( .not. is_freebound(in,k) )   apply_current_BC = .true.
          endif
                
          if (        apply_psi_BC      &
                 .or. apply_current_BC  &
                 .or. ((k .eq. var_u)    .and. apply_dirichlet_1234) &
                 .or. ((k .eq. var_w)    .and. apply_dirichlet_1234) &
                 .or. ((k .eq. var_rho)  .and. apply_dirichlet_all)  &
                 .or. ((k .eq. var_T)    .and. apply_dirichlet_all)  &
                 .or. ((k .eq. var_Te)   .and. apply_dirichlet_all)  &
                 .or. ((k .eq. var_Ti)   .and. apply_dirichlet_all)  &
                 .or. ((k .eq. var_vpar) .and. apply_dirichlet_all)  &
                 .or. ((k .eq. var_rhon) .and. apply_dirichlet_all)  &
              ) then

!            if ((k.eq.7) .and. (node_list%node(inode)%boundary .eq. 3)) cycle  !=== better included for ITER extended wall

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


          endif

        enddo !=== variables

        if ((node_list%node(inode)%boundary .eq.  3) .and. (node_list%node(inode2)%boundary .eq.  2)) cycle

        ! --- Mach1 Boundary Conditions
        if ( (.not. mach_one_bnd_integral) .and. apply_cs) then

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
          
          if (n_order .ge. 5) then
            element_size_3 =  element_list%element(ielm)%size(iv, iv_dir+3) * H1_ss(1,2) 
            if (t_constant_boundary .and. ((iv  .eq. 2) .or. (iv  .eq. 3))) element_size_3 = - element_size_3
            if (s_constant_boundary .and. ((iv  .eq. 3) .or. (iv  .eq. 4))) element_size_3 = - element_size_3
          endif

          index_node    = node_list%node(inode)%index(1)             ! position of value
          index_node2   = node_list%node(inode)%index(iv_dir)        ! position of first deriative
          if (n_order .ge. 5) &
            index_node3 = node_list%node(inode)%index(iv_dir+3)      ! position of 2nd deriative

          ! --- Determine the direction of the BCs and apply smoothing factors if requested
          ps0       = node_list%node(inode)%values(1,1,var_psi)
          ps0_b     = node_list%node(inode)%values(1,iv_dir,var_psi)  * element_size_0 
          ps0_s     = node_list%node(inode)%values(1,2,var_psi)       * element_size_s
          ps0_t     = node_list%node(inode)%values(1,3,var_psi)       * element_size_t

          BigR      = node_list%node(inode)%x(1,1,1)
          R_b       = node_list%node(inode)%x(1,iv_dir,1) * element_size_0
          Z_b       = node_list%node(inode)%x(1,iv_dir,2) * element_size_0

          R_s       = node_list%node(inode)%x(1,2,1)      * element_size_s
          R_t       = node_list%node(inode)%x(1,3,1)      * element_size_t    
          Z_s       = node_list%node(inode)%x(1,2,2)      * element_size_s
          Z_t       = node_list%node(inode)%x(1,3,2)      * element_size_t    
          Z         = node_list%node(inode)%x(1,1,2)
          
          ps0_bb = element_list%element(ielm)%size(iv ,1)      * node_list%node(inode )%values(1,1,var_psi)      * H1_ss(1,1) &
                 + element_list%element(ielm)%size(iv ,iv_dir) * node_list%node(inode )%values(1,iv_dir,var_psi) * H1_ss(1,2) &
                 + element_list%element(ielm)%size(iv2,1)      * node_list%node(inode2)%values(1,1,var_psi)      * H1_ss(2,1) &
                 + element_list%element(ielm)%size(iv2,iv_dir) * node_list%node(inode2)%values(1,iv_dir,var_psi) * H1_ss(2,2)

          R_bb = + element_list%element(ielm)%size(iv ,1)      * node_list%node(inode )%x(1,1,1)      * H1_ss(1,1)  &
                 + element_list%element(ielm)%size(iv ,iv_dir) * node_list%node(inode )%x(1,iv_dir,1) * H1_ss(1,2)  &
                 + element_list%element(ielm)%size(iv2,1)      * node_list%node(inode2)%x(1,1,1)      * H1_ss(2,1)  &
                 + element_list%element(ielm)%size(iv2,iv_dir) * node_list%node(inode2)%x(1,iv_dir,1) * H1_ss(2,2)  

          Z_bb = + element_list%element(ielm)%size(iv ,1)      * node_list%node(inode )%x(1,1,     2) * H1_ss(1,1)  &
                 + element_list%element(ielm)%size(iv ,iv_dir) * node_list%node(inode )%x(1,iv_dir,2) * H1_ss(1,2)  &
                 + element_list%element(ielm)%size(iv2,1)      * node_list%node(inode2)%x(1,1,     2) * H1_ss(2,1)  &
                 + element_list%element(ielm)%size(iv2,iv_dir) * node_list%node(inode2)%x(1,iv_dir,2) * H1_ss(2,2)  

          ps2_b     = node_list%node(inode2)%values(1,iv_dir,var_psi) * element_size_2 
          
          xjac  =  R_s*Z_t - R_t*Z_s
          ps0_x = (   Z_t * ps0_s - Z_s * ps0_t ) / xjac
          ps0_y = ( - R_t * ps0_s + R_s * ps0_t ) / xjac

          grad_s = (/  Z_t,  -R_t /) / xjac
          grad_t = (/ -Z_s,   R_s /) / xjac

          if (s_constant_boundary) then
            grad_b = grad_s
          elseif (t_constant_boundary) then
            grad_b = grad_t
          endif

          normal     = dot_product(grad_b,normal_direction) * grad_b      ! outward pointing normal
          normal     = normal / norm2(normal)
          direction  = sign(1.d0,dot_product((/ps0_y,-ps0_x/),normal))
          
          grad_psi = sqrt(ps0_x**2 + ps0_y**2)
          Btot     = sqrt(F0**2 + ps0_x**2 + ps0_y**2) / BigR
          dl       = sqrt(R_b**2 + Z_b**2)
          dl_b     = (R_b*R_bb + Z_b*Z_bb) / dl

          bn     = dot_product( (/ps0_y,-ps0_x/), normal ) /  (BigR*Btot)  ! B�n/Btot
          bn_b   = 1.d0 / (Btot*dl*BigR) * (ps0_bb - ps0_b * dl_b /dl )

          bn_1 = + ps0_b/(BigR*Btot*dl)
          bn_2 = + ps2_b/(BigR*Btot*dl)

          if ((s_constant_boundary) .and. ((iv  .eq. 1) .or. (iv  .eq. 4))) bn_1 = - bn_1
          if ((t_constant_boundary) .and. ((iv  .eq. 3) .or. (iv  .eq. 4))) bn_1 = - bn_1
          if ((s_constant_boundary) .and. ((iv2 .eq. 1) .or. (iv2 .eq. 4))) bn_2 = - bn_2
          if ((t_constant_boundary) .and. ((iv2 .eq. 3) .or. (iv2 .eq. 4))) bn_2 = - bn_2

          ! --- Apply Smoothing?
          c_1 = vpar_smoothing_coef(1); c_2 = vpar_smoothing_coef(2); c_3 = vpar_smoothing_coef(3)
          if ((vpar_smoothing) .and. (bn_1*bn_2 .lt. 0.d0)) then
            if (c_2 .gt. 0d0) then
              factor    = 0.25d0 * ( 1.d0 + tanh( (abs(bn) - c_1) / c_2 ) )**2 - c_3
              factor_b  = 0.5d0  * ( 1.d0 + tanh( (abs(bn) - c_1) / c_2 ) )           & 
                        * (bn_b * bn/abs(bn) /c_2) /(cosh( (abs(bn) - c_1) / c_2 ) )**2
              factor_bb = 0.5d0 * (bn_b * bn/abs(bn) /c_2)**2 / (cosh( (abs(bn) - c_1) / c_2 ) )**4 &
                         - 1.0d0  * ( 1.d0 + tanh( (abs(bn) - c_1) / c_2 ) ) &
                        * (bn_b * bn/abs(bn) /c_2)**2 /(cosh( (abs(bn) - c_1) / c_2 ) )**3 * sinh( (abs(bn) - c_1) / c_2 )
             else
               factor    = tanh(bn/c_1)
               factor_b  = bn_b /c_1 / cosh(bn/c_1)**2
               factor_bb = -2.d0 * (bn_b/c_1)**2 / cosh(bn/c_1)**3 * sinh(bn/c_1)
               direction = 1.d0                            
            endif                       
          else
            factor    = 1.d0
            factor_b  = 0.d0
            factor_bb = 0.d0
          endif
          Hfact_b   = factor * R_b / BigR  + factor_b
          Hfact_bb  = factor * R_bb/ BigR - factor * R_b**2 / BigR**2  + factor_bb

          ! --- For the BC's the magnetic field is assumed constant, and we do
          ! --- not take derivatives of psi into account. Vpar, T and U are the
          ! --- only variables that will be used for the linearisation
          Ti0       = max(node_list%node(inode)%values(1,1,var_Ti), T_min)
          Ti0_b     = node_list%node(inode)%values(1,iv_dir,var_Ti)   * element_size_0 

          Te0       = max(node_list%node(inode)%values(1,1,var_Te), T_min)
          Te0_b     = node_list%node(inode)%values(1,iv_dir,var_Te)   * element_size_0 

          Vpar0     = node_list%node(inode)%values(1,1,var_vpar)
          Vpar0_b   = node_list%node(inode)%values(1,iv_dir,var_Vpar) * element_size_0 

          u0_b      = node_list%node(inode)%values(1,iv_dir,var_u)    * element_size_0 

          if (n_order .ge. 5) then
            Ti0_bb    = node_list%node(inode)%values(1,iv_dir+3,var_Ti)   * element_size_3 
            Te0_bb    = node_list%node(inode)%values(1,iv_dir+3,var_Te)   * element_size_3 
            Vpar0_bb  = node_list%node(inode)%values(1,iv_dir+3,var_Vpar) * element_size_3 
            u0_bb     = node_list%node(inode)%values(1,iv_dir+3,var_u)    * element_size_3 
          endif

          ! --- Mach1 BC's and derivatives
          cs0      =   sqrt(gamma*(Ti0+Te0))
          cs0_T    =   0.5d0  * gamma    / cs0
          cs0_TT   = - 0.25d0 * gamma**2 / cs0**3 
          cs0_TTT  = 3.d0/8.d0* gamma**3 / cs0**5 

          Mach1BC     = - Vpar0   + direction / Btot * factor  * cs0               + factor / Btot * BigR**2 * U0_b/ps0_b 
          Mach1BC_v   = - 1.0
          Mach1BC_T   =           + direction / Btot * factor  * cs0_T 
          Mach1BC_u   =                                                            + factor / Btot * BigR**2 * element_size_0/ps0_b 
          dMach1BC    = - Vpar0_b + direction / Btot * factor  * cs0_T * (Ti0_b+Te0_b)  &
                                  + direction / Btot * Hfact_b * cs0         
          dMach1BC_v  = - element_size_0
          dMach1BC_Ti =           + direction / Btot * factor  * cs0_TT* (Ti0_b+Te0_b)  &
                                  + direction / Btot * Hfact_b * cs0_T
          dMach1BC_Te =           + direction / Btot * factor  * cs0_TT* (Ti0_b+Te0_b)  &
                                  + direction / Btot * Hfact_b * cs0_T
          dMach1BC_Tb =           + direction / Btot * factor  * cs0_T * element_size_0


          if (n_order .ge. 5) then
            dMach1BC     = dMach1BC + factor / Btot * BigR**2 * U0_bb/ps0_b
            dMach1BC_ubb = + factor / Btot * BigR**2 * element_size_3/ps0_b
            d2Mach1BC    = - Vpar0_bb + direction / Btot * factor   * cs0_TT * (Ti0_b+Te0_b)**2   &
                                      + direction / Btot * factor   * cs0_T  * (Ti0_bb+Te0_bb)   !&
                                      !+ direction / Btot * Hfact_b  * cs0_T  * T0_b *2.0 !&
                                      !+ direction / Btot * Hfact_bb * cs0         
            d2Mach1BC_v  = - element_size_3
            d2Mach1BC_T  =            + direction / Btot * factor   * cs0_TTT * (Ti0_b+Te0_b)**2   &
                                      + direction / Btot * factor   * cs0_TT  * (Ti0_bb+Te0_bb)    &
                                      + direction / Btot * Hfact_b  * cs0_TT  * (Ti0_b+Te0_b) *2.0 &
                                      + direction / Btot * Hfact_bb * cs0_T 
            d2Mach1BC_Tb =            + direction / Btot * factor   * cs0_TT  * (Ti0_b+Te0_b) * 2.0 * element_size_0 &
                                      + direction / Btot * Hfact_b  * cs0_T                   * 2.0 * element_size_0 
            d2Mach1BC_Tbb=            + direction / Btot * factor   * cs0_T   * element_size_3 
          endif

          ! --- Apply Mach1
          ku = var_u
          kv = var_Vpar
          kTi = var_Ti
          kTe = var_Te

          ! --- Impose Mach1 on node values
          call boundary_conditions_add_one_entry(             &
               index_node, kv, in, index_node, kv, in,        &
               - zbig * Mach1BC_v,                            &
               index_min, index_max, a_mat)

          call boundary_conditions_add_one_entry(             &
               index_node, kv, in, index_node, kTi, in,       &
               - zbig * Mach1BC_T,                            &
               index_min, index_max, a_mat)

          call boundary_conditions_add_one_entry(             &
               index_node, kv, in, index_node, kTe, in,       &
               - zbig * Mach1BC_T,                            &
               index_min, index_max, a_mat)

          call boundary_conditions_add_one_entry(             &
               index_node,  kv, in, index_node2, ku, in,      &
               - zbig * Mach1BC_u,                            &
               index_min, index_max, a_mat)

          if (in .eq. 1) then
            call boundary_conditions_add_RHS(                        &
                   index_node, kv, in,index_min, index_max, RHS_loc, &
                   Zbig * Mach1BC,                                   &
                   a_mat%i_tor_min, a_mat%i_tor_max)
          else
            call boundary_conditions_add_RHS(                         &
                   index_node, kv, in, index_min, index_max, RHS_loc, &
                   0.d0,                                              &
                   a_mat%i_tor_min, a_mat%i_tor_max)
          endif
  
          ! --- Impose Mach1 on node derivatives
          call boundary_conditions_add_one_entry(               &
                 index_node2, kv, in, index_node2, kv, in,      &
                 - zbig * dMach1BC_v,                           &
                 index_min, index_max, a_mat)

          call boundary_conditions_add_one_entry(               &
                 index_node2, kv, in, index_node2, kTi, in,     &
                 - zbig * dMach1BC_Tb,                          &
                 index_min, index_max, a_mat)

          call boundary_conditions_add_one_entry(               &
                 index_node2, kv, in, index_node2, kTe, in,     &
                 - zbig * dMach1BC_Tb,                          &
                 index_min, index_max, a_mat)

          call boundary_conditions_add_one_entry(               &
                 index_node2, kv, in, index_node,  kTi, in,     &
                 - zbig * dMach1BC_Ti,                          & 
                 index_min, index_max, a_mat)

          call boundary_conditions_add_one_entry(               &
                 index_node2, kv, in, index_node,  kTe, in,     &
                 - zbig * dMach1BC_Te,                          & 
                 index_min, index_max, a_mat)

          if ( include_2nd_derivatives .and. (n_order .ge. 5) ) then
            call boundary_conditions_add_one_entry(               &
                   index_node2, kv, in, index_node3, ku, in,      &
                   - zbig * dMach1BC_ubb,                         &
                   index_min, index_max, a_mat)
          endif

          if (in .eq. 1) then
            call boundary_conditions_add_RHS(                          &
                   index_node2, kv, in, index_min, index_max, RHS_loc, &
                   Zbig * dMach1BC,                                    &
                   a_mat%i_tor_min, a_mat%i_tor_max)
          else
             call boundary_conditions_add_RHS(                         &
                   index_node2, kv, in, index_min, index_max, RHS_loc, &
                   0.d0,                                               &
                   a_mat%i_tor_min, a_mat%i_tor_max) 
          endif

          ! --- Impose Mach1 on node 2nd derivatives
!DOESNT WORK. I DONT KNOW WHY...
          if ( include_2nd_derivatives .and. (n_order .ge. 5) ) then
            call boundary_conditions_add_one_entry(               &
                   index_node3, kv, in, index_node3, kv,  in,     &
                   - zbig * d2Mach1BC_v,                          &
                   index_min, index_max, a_mat)
            call boundary_conditions_add_one_entry(               &
                   index_node3, kv, in, index_node , kTi, in,     &
                   - zbig * d2Mach1BC_T,                          &
                   index_min, index_max, a_mat)
            call boundary_conditions_add_one_entry(               &
                   index_node3, kv, in, index_node , kTe, in,     &
                   - zbig * d2Mach1BC_T,                          &
                   index_min, index_max, a_mat)
            call boundary_conditions_add_one_entry(               &
                   index_node3, kv, in, index_node2, kTi, in,     &
                   - zbig * d2Mach1BC_Tb,                         &
                   index_min, index_max, a_mat)
            call boundary_conditions_add_one_entry(               &
                   index_node3, kv, in, index_node2, kTe, in,     &
                   - zbig * d2Mach1BC_Tb,                         &
                   index_min, index_max, a_mat)
            call boundary_conditions_add_one_entry(               &
                   index_node3, kv, in, index_node3, kTi, in,     &
                   - zbig * d2Mach1BC_Tbb,                        &
                   index_min, index_max, a_mat)
            call boundary_conditions_add_one_entry(               &
                   index_node3, kv, in, index_node3, kTe, in,     &
                   - zbig * d2Mach1BC_Tbb,                        &
                   index_min, index_max, a_mat)
            if (in .eq. 1) then
              call boundary_conditions_add_RHS(                           &
                     index_node3, kv, in, index_min, index_max, RHS_loc,  &
                     Zbig * d2Mach1BC,                                    &
                     a_mat%i_tor_min, a_mat%i_tor_max)
            else
               call boundary_conditions_add_RHS(                         &
                     index_node3, kv, in, index_min, index_max, RHS_loc, &
                     0.d0,                                               &
                     a_mat%i_tor_min, a_mat%i_tor_max) 
            endif
          endif

          ! --- Fix derivatives in one direction
          k = var_Vpar
          do kk = 1,(n_order+1)/2
            do ll = 1,(n_order+1)/2
              if ( (iv_dir .eq. 2) .and. (ll .gt. 1) ) cycle ! do only pure s derivatives, not cross _st
              if ( (iv_dir .eq. 3) .and. (kk .gt. 1) ) cycle ! do only pure t derivatives, not cross _st
              if ( (iv_dir .eq. 2) .and. (kk .lt. 3) ) cycle ! do only node value, 1st and 2nd derivatives, fix the rest
              if ( (iv_dir .eq. 3) .and. (ll .lt. 3) ) cycle ! do only node value, 1st and 2nd derivatives, fix the rest
              index_tmp = node_indices(kk,ll)
              index_node = node_list%node(inode)%index(index_tmp)
              call boundary_conditions_add_one_entry(                 &
                     index_node, k, in, index_node, k, in,            &
                     zbig, index_min, index_max, a_mat)
            enddo
          enddo

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

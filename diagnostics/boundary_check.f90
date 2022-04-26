!> Compares \f$B_{tan}\f$ and \f$B_{tan,vacuum}\f$ at the interface.
!!
!! The relative difference is calculated according to:
!! \f$\frac{\int dA |B_{tan,plasma}-B_{tan,vacuum}|}{\int dA |B_{tan,vacuum}|}\f$.
subroutine boundary_check(my_id)
  
  use tr_module
  use data_structure,  only: type_bnd_element 
  use phys_module,     only: resistive_wall
  use nodes_elements,  only: node_list, bnd_node_list, element_list, bnd_elm_list
  use vacuum_response
  use mpi_mod
  use mod_interp
  use mod_basisfunctions
  use mod_parameters, only: n_degrees_1d
  
  implicit none
  
  integer,                      intent(in)    :: my_id

  ! --- Local variables
  integer, parameter     :: N_POINTS = 3 ! Number of evaluation points per element
  type(type_bnd_element) :: bndelem_m
  real*8, allocatable    :: B_par(:), B_par_v(:)
  real*8, allocatable    :: val_integral(:), err_integral(:)
  real*8, allocatable    :: psibnd_vec(:), dpsibnd_vec(:), psibnd_coils(:)
  integer  :: l_starwall, l_tor
  integer  :: m_bndelem, m_pt, m_elm, mv1
  integer  :: i_vertex, i_dof, i_node, i_node_bnd, i_resp, i_resp_old, i_resp_0
  real*8   :: i_size, basfunc_i
  real*8   :: H1(2,n_degrees_1d), H1_s(2,n_degrees_1d), H1_ss(2,n_degrees_1d)
  real*8   :: P, P_s, P_t, P_st, P_ss, P_tt
  real*8   :: R, R_s, R_t, Z, Z_s, Z_t
  real*8   :: s_pt, t_pt, s_or_t ! s and t values at current point
  real*8   :: xjac               ! 2D Jacobian
  real*8   :: B_pol(2)           ! Poloidal magnetic field
  real*8   :: e_par(2)           ! Vector tangential to interface
  real*8   :: P_R, P_Z           ! dPsi/dR, dPsi/dZ
  real*8   :: R1, R2, Z1, Z2, B_par_v_tmp
  logical  :: s_const            ! Is the bound. elem. an s=const side of the 2D element?
  integer  :: ierr,step          ! variables for parallel version 
  
  if(my_id == 0) then
    write(*,*) '************************************'
    write(*,*) '*    check boundary conditions     *'
    write(*,*) '************************************'
  endif

  call tr_allocate(val_integral,1,sr%n_tor,"val_integral",CAT_GRID)
  call tr_allocate(err_integral,1,sr%n_tor,"err_integral",CAT_GRID)
  call tr_allocate(psibnd_vec,1,n_dof_starwall,"psibnd_vec",CAT_GRID)
  call tr_allocate(psibnd_coils,1,n_dof_starwall,"psibnd_vec",CAT_GRID)
  call tr_allocate(dpsibnd_vec,1,n_dof_starwall,"dpsibnd_vec",CAT_GRID)
  call tr_allocate(B_par,1,sr%n_tor,"B_par",CAT_GRID)
  call tr_allocate(B_par_v,1,sr%n_tor,"B_par_v",CAT_GRID)
  val_integral(:) = 0.d0
  err_integral(:) = 0.d0
  B_par(:)        = 0.d0
  B_par_v(:)      = 0.d0
  i_resp_0        = 0 

  call det_psibnd_vec(bnd_node_list, node_list, psibnd_vec, dpsibnd_vec, psibnd_coils)

  ! --- For every boundary element, do...
  L_MB: do m_bndelem = 1, bnd_elm_list%n_bnd_elements

    bndelem_m = bnd_elm_list%bnd_element(m_bndelem)
    m_elm     = bnd_elm_list%bnd_element(m_bndelem)%element
    mv1       = bnd_elm_list%bnd_element(m_bndelem)%side

    R1 = node_list%node(bndelem_m%vertex(1))%x(1,1,1)
    Z1 = node_list%node(bndelem_m%vertex(1))%x(1,1,2)
    R2 = node_list%node(bndelem_m%vertex(2))%x(1,1,1)
    Z2 = node_list%node(bndelem_m%vertex(2))%x(1,1,2)

    ! --- For several points in the boundary element, do...
    L_MP: do m_pt = 1, N_POINTS

      B_par(:)   = 0.d0
      B_par_v(:) = 0.d0

      ! --- Determine 1D basis function (and derivatives) at current point
      s_or_t = float(m_pt-1)/float(N_POINTS-1)

      call basisfunctions1(s_or_t, H1, H1_s, H1_ss)

      ! --- Which s and t values correspond to the current point and is the
      !     boundary element an s=const or t=const side of the 2D element?
      select case (mv1)
      case (1)
        s_pt = s_or_t;  t_pt = 0.d0;    s_const = .false.
      case (2)
        s_pt = 1.d0;    t_pt = s_or_t;  s_const = .true.
      case (3)
        s_pt = s_or_t;  t_pt = 1.d0;    s_const = .false.
      case (4)
        s_pt = 0.d0;    t_pt = s_or_t;  s_const = .true.
      end select

      ! --- Determine coordinate values (plus derivatives)
      call interp_RZ(node_list, element_list, m_elm, s_pt, t_pt, R, R_s, R_t, Z, Z_s, Z_t)

      ! --- 2D Jacobian
      xjac = R_s * Z_t - R_t * Z_s

      ! --- Tangential vector to the interface
      if ( s_const ) then
        e_par = (/ R_t, Z_t /) / sqrt( R_t**2 + Z_t**2 ) * (R_t * (R2-R1) + Z_t * (Z2-Z1))/abs(R_t * (R2-R1) + Z_t * (Z2-Z1))
      else
        e_par = (/ R_s, Z_s /) / sqrt( R_s**2 + Z_s**2 ) * (R_s * (R2-R1) + Z_s * (Z2-Z1))/abs(R_s * (R2-R1) + Z_s * (Z2-Z1))
      end if

      step       = sr%a_ye%step ! for distributed matrices
      B_par_v(:) = 0.d0
      
      ! --- Select one STARWALL harmonic
!$omp parallel do                                                                                  &
!$omp   schedule(dynamic)                                                                          &
!$omp   default(none)                                                                              &
!$omp   shared(node_list, element_list, m_elm, s_pt, t_pt, sr, bndelem_m, bnd_node_list, H1,       &
!$omp     resistive_wall, starwall_equil_coils, my_id, psibnd_vec, psibnd_coils, wall_curr,        &
!$omp     bext_tan, I_coils, B_par_v, R, R_s, R_t, Z, Z_s, Z_t, &
!$omp     xjac, e_par, B_par, step)                                                                &
!$omp   private(l_starwall, i_vertex, i_dof, l_tor, P, P_s, P_t, P_st, P_ss, P_tt, P_R, P_Z, B_pol,&
!$omp     i_node, i_node_bnd, i_size, i_resp_old, i_resp, i_resp_0, basfunc_i, B_par_v_tmp)
      L_LS: do l_starwall = 1, sr%n_tor

        l_tor = sr%i_tor(l_starwall)

        ! --- Psi value (plus derivatives) at current point (l_tor mode)
        call interp(node_list, element_list, m_elm, 1, l_tor, s_pt, t_pt, P, P_s, P_t, P_st, P_ss, P_tt)

        ! --- Poloidal magnetic field at current point
        P_R   = (   P_s * Z_t - P_t * Z_s ) / xjac ! dPsi/dR
        P_Z   = ( - P_s * R_t + P_t * R_s ) / xjac ! dPsi/dZ
        B_pol = (/ P_Z, -P_R /) / R

        ! --- Tangential magnetic field B_{||} reconstructed from the plasma
!$omp critical
        B_par(l_starwall) = - sum( B_pol * e_par )
!$omp end critical

        ! --- Sum over boundary dofs at which response is calculated
        L_IV: do i_vertex = 1, 2 ! (loop over nodes in element m_bndelem)

          i_node      = bndelem_m%vertex(i_vertex)
          i_node_bnd  = bndelem_m%bnd_vertex(i_vertex)

          L_ID: do i_dof = 1, 2 ! (loop over node dofs)

            i_size      = bndelem_m%size(i_vertex,i_dof)

            i_resp_old  = response_index(i_node_bnd,l_starwall,i_dof)

            i_resp   = (bnd_node_list%bnd_node(i_node_bnd)%index_starwall(1) - 1)*sr%n_tor &
                     + bnd_node_list%bnd_node(i_node_bnd)%n_dof*(l_starwall-1) &
                     + bnd_node_list%bnd_node(i_node_bnd)%index_starwall(i_dof)-bnd_node_list%bnd_node(i_node_bnd)%index_starwall(1) + 1
                   
            i_resp_0 = response_index_eq(i_node_bnd,i_dof)

            ! --- Determine basis function
            basfunc_i = H1(i_vertex,i_dof) * i_size

            ! --- Determine B_{||,v} as prescribed by the vacuum.
            B_par_v_tmp = 0.d0
            if ( resistive_wall ) then
              if (  (l_tor == 1) .and. (.not. starwall_equil_coils)  )  then

                if (i_resp>=sr%a_ey%ind_start .AND. i_resp<=sr%a_ey%ind_end) then
                  B_par_v_tmp = basfunc_i * (                                                      &
                    + sum( sr%a_ee%loc_mat(i_resp-step*my_id,:) * (psibnd_vec(:)-psibnd_coils(:)) )&
                    + sum( sr%a_ey%loc_mat(i_resp-step*my_id,:) * wall_curr(:) )                   &
                    - sum( bext_tan(i_resp_0, :) * I_coils(:) ) )
                end if

              else

                if (i_resp>=sr%a_ey%ind_start .AND. i_resp<=sr%a_ey%ind_end) then
                  B_par_v_tmp = basfunc_i * (                                                      &
                    + sum( sr%a_ee%loc_mat(i_resp-step*my_id, :) * psibnd_vec(:) )                 &
                    + sum( sr%a_ey%loc_mat(i_resp-step*my_id, :) * wall_curr(:)  ) )
                end if
                 

              end if
            else ! if ( resistive_wall ) then

              if (  (l_tor == 1) .and. (.not. starwall_equil_coils)  )  then
                
                if (i_resp>=sr%a_ey%ind_start .AND. i_resp<=sr%a_ey%ind_end) then
                  B_par_v_tmp = basfunc_i * (                                                      &
                    sum( sr%a_id%loc_mat(i_resp-my_id*step,:) * (psibnd_vec(:) - psibnd_coils(:)) )&
                    - sum( bext_tan(i_resp_0, :) * I_coils(:) ) )
                end if

              else
                
                if (i_resp>=sr%a_ey%ind_start .AND. i_resp<=sr%a_ey%ind_end) then
                  B_par_v_tmp = basfunc_i *                                                        &
                    sum( sr%a_id%loc_mat(i_resp-step*my_id, :) * psibnd_vec(:) )
                end if
                
              end if
            end if
            
!$omp critical
            B_par_v(l_starwall) = B_par_v(l_starwall) + B_par_v_tmp
!$omp end critical

!            write(*,'(6i5,8e12.4)') l_starwall,i_vertex,i_dof,i_node,i_node_bnd,i_resp,i_size,basfunc_i
          end do L_ID
        end do L_IV

      end do L_LS
!$omp end parallel do
 
      call MPI_AllREDUCE(MPI_IN_PLACE,B_par_v, size(B_par_v), MPI_DOUBLE_PRECISION,MPI_SUM, MPI_COMM_WORLD,ierr)

      ! --- Debugging output
      if ( vacuum_debug .and. (my_id == 0) ) then
        write(88,'(20ES15.5)') (m_bndelem-1 + s_or_t)/REAL(bnd_elm_list%n_bnd_elements), B_par(:)
        write(89,'(20ES15.5)') (m_bndelem-1 + s_or_t)/REAL(bnd_elm_list%n_bnd_elements), B_par_v(:)
      end if

      ! --- Integration of B_par_v values and differences between B_par and B_par_v.
      val_integral(:) = val_integral(:) + abs( B_par_v(:) )
      err_integral(:) = err_integral(:) + abs( B_par(:) - B_par_v(:) )

    end do L_MP

  end do L_MB
  
  if ( (minval(abs(val_integral)) /= 0.d0) .and. (my_id == 0) ) then
    write(*,'(1x,A,99ES15.5)') 'Relative errors in harmonics:', err_integral(:) / val_integral(:), err_integral(:), val_integral(:)
  end if

  ! --- Debugging output
  if ( vacuum_debug .and. (my_id == 0) ) then
    write(88,*)
    write(88,*)
    write(89,*)
    write(89,*)
    if ( minval(abs(val_integral)) /= 0.d0 ) then ! (avoid division by zero in first timestep)
      write(87,'(20ES15.5)') err_integral(:) / val_integral(:)
    end if
  end if

  call tr_deallocate(psibnd_vec,"psibnd_vec",CAT_GRID)
  call tr_deallocate(dpsibnd_vec,"dpsibnd_vec",CAT_GRID)
  call tr_deallocate(B_par,"B_par",CAT_GRID)
  call tr_deallocate(B_par_v,"B_par_v",CAT_GRID)
  call tr_deallocate(val_integral,"val_integral",CAT_GRID)
  call tr_deallocate(err_integral,"err_integral",CAT_GRID)

end subroutine boundary_check

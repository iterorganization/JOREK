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
!*   index_max    - Maximal index of local elements                            *
!*   xpoint2      -                                                            *
!*   xcase2       -                                                            *
!*   psi_axis     -                                                            *
!*   psi_bnd      -                                                            *
!*   Z_xpoint     -                                                            *
!*   gmres        - boolean indicating if we are using GMRES method            *
!*   solve_only   - Indicate if we want to perform only solve                  *
!*                                                                             *
!*******************************************************************************
!* Subroutine: solve_Psi_boundary_eqn                                          *
!*******************************************************************************
!*                                                                             *
!* Solve the differential equation for Psi on boundary, ensuring that n.B = 0  *
!* (see eq 4.23 in N. Nikulsin's PhD thesis, doi:10.17617/2.3359934)           *
!* The solution is then used as an imhomogeneous Dirichlet b.c. for Psi        *
!*                                                                             *
!*******************************************************************************
!* Subroutine: setup_boundary_condition                                        *
!*******************************************************************************
!*                                                                             *
!* Projects the solution of the boundary Psi equation onto the JOREK finite    *
!*  element basis and writes the calculated dofs into the appropriate nodes.   *
!*                                                                             *
!*******************************************************************************
module mod_boundary_conditions
implicit none

real*8, dimension(:), allocatable :: Psi_b
real*8, dimension(:), allocatable :: J_b

private :: Zn, Zn_chi, Zm, Zm_tht, Zm_tht_2, HZn, Zjm, Zjm_tht

contains
  subroutine boundary_conditions( my_id, node_list, element_list, bnd_node_list, local_elms,          &
                                  n_local_elms, index_min, index_max, rhs_loc, xpoint2, xcase2,       & 
                                  R_axis, Z_axis, psi_axis, psi_bnd, R_xpoint, Z_xpoint, psi_xpoint,  &  
                                  gmres, solve_only, ijA_index, ijA_size, irn_jcn, irn, jcn,          & 
                                  A_mat, i_tor_min, i_tor_max )

    use data_structure
    use phys_module, only: F0, GAMMA, bc_natural_open
    use vacuum, only: is_freebound
    use mpi_mod
    use mod_locate_irn_jcn
    use mod_integer_types

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
    real*8,                             intent(inout) :: rhs_loc(*)
    logical,                            intent(in)    :: xpoint2
    integer,                            intent(in)    :: xcase2
    real*8,                             intent(in)    :: R_axis
    real*8,                             intent(in)    :: Z_axis
    real*8,                             intent(in)    :: psi_axis
    real*8,                             intent(in)    :: psi_bnd
    real*8,                             intent(in)    :: R_xpoint(2)
    real*8,                             intent(in)    :: Z_xpoint(2)
    real*8,                             intent(in)    :: psi_xpoint(2)
    logical,                            intent(in)    :: gmres
    logical,                            intent(in)    :: solve_only
    integer,                            intent(in)    :: i_tor_min, i_tor_max 
    integer(kind=int_all), allocatable, intent(in)    :: ijA_index(:,:), ijA_size(:), irn_jcn(:,:) 
    integer(kind=int_all), allocatable, intent(inout) :: irn(:), jcn(:) 
    real*8,                allocatable, intent(inout) :: A_mat(:) 

    ! Internal parameters
    real*8                :: zbig
    integer               :: i, in, iv, inode, k
    integer               :: ielm
    integer               :: index_node
    integer(kind=int_all) :: ijA_position
    integer               :: ilarge2
    integer               :: ierr, n_tor_local

    n_tor_local = i_tor_max - i_tor_min + 1
    zbig = 1.d12
       do i=1, n_local_elms

          ielm = local_elms(i)

          do iv=1, n_vertex_max

             inode = element_list%element(ielm)%vertex(iv)

             if (node_list%node(inode)%boundary .ne. 0) then

                do in=i_tor_min, i_tor_max 

                   do k=1, n_var
                      if (bc_natural_open .and. k .eq. var_zj) cycle

                      !------------------------------------ the open field lines (in case of x-point grid)
                      if ((node_list%node(inode)%boundary .eq. 1) .or. (node_list%node(inode)%boundary .eq. 3)) then

                         if ((k .eq. 1) .or. (k .eq. 2) .or. (k .eq. 3) .or. &
                              (k .eq. 4) .or. (k .eq. 5) .or. (k .eq. 6) ) then
 
                          if ( (.not. is_freebound(in,k)) ) then ! apply fixed boundary conditions where necessary

                            index_node = node_list%node(inode)%index(1)
                            if ((index_node .ge. index_min) .and. (index_node .le. index_max)) then

                               call locate_irn_jcn(index_node,index_node,index_min,index_max,ijA_position,ijA_index, ijA_size, irn_jcn)

                               ilarge2 = ijA_position - 1 + ((k-1)*n_tor_local + in-i_tor_min) * n_var*n_tor_local  & 
                                 +  (k-1)*n_tor_local + in - i_tor_min + 1

                               irn(ilarge2) = n_tor_local * n_var * (index_node-1) + (k-1)*n_tor_local + in - i_tor_min + 1
                               jcn(ilarge2) = n_tor_local * n_var * (index_node-1) + (k-1)*n_tor_local + in - i_tor_min + 1
                               A_mat(ilarge2)   = zbig

                            endif

                            index_node = node_list%node(inode)%index(2)

                            if ((index_node .ge. index_min) .and. (index_node .le. index_max)) then

                               call locate_irn_jcn(index_node,index_node,index_min,index_max,ijA_position,ijA_index, ijA_size, irn_jcn)

                               ilarge2 = ijA_position - 1 + ((k-1)*n_tor_local + in-i_tor_min) * n_var*n_tor_local   & 
                                 +  (k-1)*n_tor_local + in - i_tor_min + 1

                               irn(ilarge2) = n_tor_local * n_var * (index_node-1) + (k-1)*n_tor_local + in - i_tor_min + 1
                               jcn(ilarge2) = n_tor_local * n_var * (index_node-1) + (k-1)*n_tor_local + in - i_tor_min + 1
                               A_mat(ilarge2)    = zbig

                            endif

                          endif
                        endif
                      endif

                      !------------------------------------ wall aligned with fluxsurface (in case of x-point grid)
                      if ((node_list%node(inode)%boundary .eq. 2) .or. (node_list%node(inode)%boundary .eq. 3)) then

                         if ( (.not. is_freebound(in,k)) ) then ! apply fixed boundary conditions where necessary

                            index_node = node_list%node(inode)%index(1)

                            if ((index_node .ge. index_min) .and. (index_node .le. index_max)) then

                               call locate_irn_jcn(index_node,index_node,index_min,index_max,ijA_position,ijA_index, ijA_size, irn_jcn)

                               ilarge2 = ijA_position - 1 + ((k-1)*n_tor_local + in-i_tor_min) * n_var*n_tor_local   & 
                                 +  (k-1)*n_tor_local + in - i_tor_min + 1

                               irn(ilarge2) = n_tor_local * n_var * (index_node-1) + (k-1)*n_tor_local + in - i_tor_min + 1
                               jcn(ilarge2) = n_tor_local * n_var * (index_node-1) + (k-1)*n_tor_local + in - i_tor_min + 1
                               A_mat(ilarge2)   = zbig

                            endif

                            index_node = node_list%node(inode)%index(3)

                            if ((index_node .ge. index_min) .and. (index_node .le. index_max)) then

                               call locate_irn_jcn(index_node,index_node,index_min,index_max,ijA_position,ijA_index, ijA_size, irn_jcn)

                               ilarge2 = ijA_position - 1 + ((k-1)*n_tor_local + in-i_tor_min) * n_var*n_tor_local   & 
                                 +  (k-1)*n_tor_local + in - i_tor_min + 1

                               irn(ilarge2) = n_tor_local * n_var * (index_node-1) + (k-1)*n_tor_local + in - i_tor_min + 1
                               jcn(ilarge2) = n_tor_local * n_var * (index_node-1) + (k-1)*n_tor_local + in - i_tor_min + 1
                               A_mat(ilarge2)   = zbig
                            end if

                         endif

                      endif

                   enddo

                enddo
             endif
          enddo
       enddo

    return
  end subroutine boundary_conditions
  
  subroutine solve_J_boundary_eqn(node_list, element_list, boundary_list)
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Routine to calculate the value of the JOREK plasma current, zj, on the plasma boundary based on the GVEC magnetic field:
  !
  ! \int\int v*zj*J*dtheta*dphi = -\int\int v*zj_GVEC*J*dtheta*dphi
  !
  ! where v = Zj_m(theta)*Z_n(phi) is the test function, zj_GVEC = F0 / Bv^2 * (grad(Chi). curl(B_GVEC)),
  ! and J is the local Jacobian of the plasma boundary surface element.
  !
  ! Note that Zj_m is different to Z_m, used in solve_Psi_boundary_eqn, because the m=0 component needs to be included, when solving 
  ! for the current.
  !---------------------------------------------------------------------------------------------------------------------------------
    use constants, only: pi
    use mod_parameters, only: n_period, n_plane, n_order, n_tor, n_coord_tor
    use phys_module, only: F0, m_pol_bc
    use gauss
    use basis_at_gaussian
    use mod_interp
    use data_structure
    use mod_chi
    use mod_SolveMN
    implicit none

    type(type_node_list),        intent(in) :: node_list
    type(type_element_list), intent(in)     :: element_list
    type(type_bnd_element_list), intent(in) :: boundary_list

    real*8  :: phi, theta, element_size_ij, element_size_perp, BigR, grad_chi(3), j_gvec(3)
    integer :: N_tht, mp, ms, ielm, ibnd_elm, i, j, j2, n, m, nn, mm, in, ind1, ind2, info, m_pol_bc_jb

    type(type_bnd_element) :: bnd_element
    type(type_node)        :: node

    real*8  :: cross_deriv(3), dA, xjac

    real*8, dimension(:,:,:), allocatable :: x_g, x_s, x_p, y_g, y_s, y_p
    real*8                :: R,R_s,R_t,R_phi,R_st,R_ss,R_tt,R_sp,R_tp,R_pp
    real*8                :: Z,Z_s,Z_t,Z_p,Z_st,Z_ss,Z_tt,Z_sp,Z_tp,Z_pp
    real*8                :: BR, BR_R, BR_Z, BR_p
    real*8                :: BZ, BZ_R, BZ_Z, BZ_p
    real*8                :: BP, BP_R, BP_Z, BP_p
    real*8                :: BRg, BRg_s, BRg_t, BRg_st, BRg_ss, BRg_tt
    real*8                :: BZg, BZg_s, BZg_t, BZg_st, BZg_ss, BZg_tt
    real*8                :: BPg, BPg_s, BPg_t, BPg_st, BPg_ss, BPg_tt
    real*8                :: jz, JRg, JPg, JZg
    real*8, dimension(:,:),   allocatable :: Amat
    real*8, dimension(:),     allocatable :: RHS

    real*8, dimension(0:n_order-1,0:n_order-1,0:n_order-1) :: chi
    m_pol_bc_jb = m_pol_bc+1   ! include m=0 mode

    write(*,*) "***************************************"
    write(*,*) "*    Solving boundary J equation    *"
    write(*,*) "***************************************"
    write(*,'(A,I4,A,I4)') "Toroidal modes: ", n_tor
    write(*,'(A,I4,A,I4)') "Poloidal modes: ", m_pol_bc_jb

    N_tht = boundary_list%n_bnd_elements

    allocate(x_g(n_plane,n_gauss,N_tht)); allocate(x_s, x_p, mold=x_g)
    allocate(y_g(n_plane,n_gauss,N_tht)); allocate(y_s, y_p, mold=y_g)
    allocate(Amat(n_tor*m_pol_bc_jb,n_tor*m_pol_bc_jb))
    allocate(RHS(n_tor*m_pol_bc_jb))

    x_g  = 0.d0; x_s = 0.d0; x_p = 0.d0
    y_g  = 0.d0; y_s = 0.d0; y_p = 0.d0
    
    do ibnd_elm=1,N_tht
      bnd_element = boundary_list%bnd_element(ibnd_elm)

      do i=1,2
        node = node_list%node(bnd_element%vertex(i))

        do j=1,2
          j2 = bnd_element%direction(i,j)
          element_size_ij = bnd_element%size(i,j)

          do ms=1,n_gauss
            do mp=1,n_plane
              do in=1,n_coord_tor
                x_g(mp,ms,ibnd_elm)  = x_g(mp,ms,ibnd_elm)  + node%x(in,j2,1)*element_size_ij*H1(i,j,ms)   *HZ_coord(in,mp)
                x_s(mp,ms,ibnd_elm)  = x_s(mp,ms,ibnd_elm)  + node%x(in,j2,1)*element_size_ij*H1_s(i,j,ms) *HZ_coord(in,mp)
                x_p(mp,ms,ibnd_elm)  = x_p(mp,ms,ibnd_elm)  + node%x(in,j2,1)*element_size_ij*H1(i,j,ms)   *HZ_coord_p(in,mp)

                y_g(mp,ms,ibnd_elm)  = y_g(mp,ms,ibnd_elm)  + node%x(in,j2,2)*element_size_ij*H1(i,j,ms)   *HZ_coord(in,mp)
                y_s(mp,ms,ibnd_elm)  = y_s(mp,ms,ibnd_elm)  + node%x(in,j2,2)*element_size_ij*H1_s(i,j,ms) *HZ_coord(in,mp)
                y_p(mp,ms,ibnd_elm)  = y_p(mp,ms,ibnd_elm)  + node%x(in,j2,2)*element_size_ij*H1(i,j,ms)   *HZ_coord_p(in,mp)

              end do
            end do
          end do
        end do
      end do
    end do

    ! Integration over a single field period
    Amat = 0.d0; RHS = 0.d0
    open(21, file='j_boundary_points.dat', action='write', status='replace')
    do ibnd_elm=1,N_tht
      ielm = boundary_list%bnd_element(ibnd_elm)%element

      do ms=1,n_gauss
        theta = 2.d0*pi*(float(ibnd_elm-1) + xgauss(ms))/float(N_tht)

        do mp=1,n_plane
          BigR = x_g(mp,ms,ibnd_elm)
          phi = 2.d0*pi*float(mp-1)/float(n_plane*n_period)
          chi = get_chi(x_g(mp,ms,ibnd_elm),y_g(mp,ms,ibnd_elm),phi)
          grad_chi = (/ chi(1,0,0), chi(0,1,0), chi(0,0,1)/BigR /)

          ! Calculate local surface area jacobian from covariant components
          cross_deriv = (/(x_s(mp,ms,ibnd_elm)*sin(phi)*y_p(mp,ms,ibnd_elm)-(x_p(mp,ms,ibnd_elm)*sin(phi)+x_g(mp,ms,ibnd_elm)*cos(phi))*y_s(mp,ms,ibnd_elm)),   &
                        (-x_s(mp,ms,ibnd_elm)*y_p(mp,ms,ibnd_elm)*cos(phi)+(x_p(mp,ms,ibnd_elm)*cos(phi)-x_g(mp,ms,ibnd_elm)*sin(phi))*y_s(mp,ms,ibnd_elm)),  &
                        (x_s(mp,ms,ibnd_elm) * x_g(mp,ms,ibnd_elm))  /)
          dA = sqrt(cross_deriv(1)*cross_deriv(1) + cross_deriv(2)*cross_deriv(2) + cross_deriv(3)*cross_deriv(3)) 
          
          ! Calculate GVEC current from finite element representation of B_GVEC
          call interp_RZP(node_list,element_list,ielm,1.0,xgauss(ms),phi,R,R_s,R_t,R_phi,R_st,R_ss,R_tt,R_sp,R_tp,R_pp, &
                      Z,Z_s,Z_t,Z_p,Z_st,Z_ss,Z_tt,Z_sp,Z_tp,Z_pp)
          xjac  = R_s * Z_t - R_t * Z_s
          JRg = 0.d0; JZg = 0.d0; JPg = 0.d0
          do in=1, n_coord_tor
            call interp_gvec(node_list,element_list,ielm,1,1,in,1.d0,xgauss(ms),BRg,BRg_s,BRg_t,BRg_st,BRg_ss,BRg_tt)
            call interp_gvec(node_list,element_list,ielm,1,2,in,1.d0,xgauss(ms),BZg,BZg_s,BZg_t,BZg_st,BZg_ss,BZg_tt)
            call interp_gvec(node_list,element_list,ielm,1,3,in,1.d0,xgauss(ms),Bpg,Bpg_s,Bpg_t,Bpg_st,Bpg_ss,Bpg_tt)
            BR_R  = (   Z_t * BRg_s - Z_s * BRg_t )     / xjac * HZ_coord(in,mp)
            BR_Z  = ( - R_t * BRg_s + R_s * BRg_t )     / xjac * HZ_coord(in,mp)
            BR_p  = BRg * HZ_coord_p(in, mp) - R_phi * BR_R - Z_p * BR_Z
            BZ_R  = (   Z_t * BZg_s - Z_s * BZg_t )     / xjac * HZ_coord(in,mp)
            BZ_Z  = ( - R_t * BZg_s + R_s * BZg_t )     / xjac * HZ_coord(in,mp)
            BZ_p  = BZg * HZ_coord_p(in, mp) - R_phi * BZ_R - Z_p * BZ_Z
            BP_R  = (   Z_t * BPg_s - Z_s * BPg_t )     / xjac * HZ_coord(in,mp)
            BP_Z  = ( - R_t * BPg_s + R_s * BPg_t )     / xjac * HZ_coord(in,mp)
            BP_p  = BPg * HZ_coord_p(in, mp) - R_phi * BP_R - Z_p * BP_Z
            
            JRg = JRg + 1 / BigR * BZ_p - BP_Z
            JZg = JZg + 1 / BigR * (BPg + BigR * BP_R - BR_p)
            JPg = Jpg + BR_Z - BZ_R
          enddo
          j_gvec = (/ JRg, JZg, Jpg /)
          jz = F0 / (chi(1,0,0)**2 + chi(0,1,0)**2 + chi(0,0,1)**2/BigR**2) * dot_product(j_gvec,grad_chi)
          
          ! Output diagnostic of the boundary current
          write(21, *) theta, phi, jz

          ind1 = 1
          do n=1,n_tor
            do m=1,m_pol_bc_jb
              ! Add gvec current as source to RHS
              RHS(ind1) = RHS(ind1) + jz*HZn(n,phi)*Zjm(m,theta)*dA*wgauss(ms)
              
              ind2 = 1
              ! Loop over basis functions in the 2D Fourier decomposition of zj
              do nn=1,n_tor
                do mm=1,m_pol_bc_jb
                  Amat(ind1,ind2) = Amat(ind1,ind2) + HZn(n,phi)*Zjm(m,theta)*HZn(nn,phi)*Zjm(mm,theta)*dA*wgauss(ms)
          
                  ind2 = ind2 + 1
                end do
              end do
  
              ind1 = ind1 + 1
            end do
          end do
        end do
      end do
    end do
    close(21)

    ! Invert matrix and solve for J_b
    call SolveMN(Amat, info)
    if (info .ne. 0) then
      write(*,'(A,I3)') "ERROR: solve_J_boundary_eqn, matrix inversion failed: ", info
      stop
    end if
    if (allocated(J_b)) deallocate(J_b)
    allocate(J_b(n_tor*m_pol_bc_jb))
    J_b = matmul(Amat,RHS)

    ! Output coefficients as a diagnostic
    open(21,file='j_boundary_coeff.dat')
    do n=1,n_tor
      do m=1,m_pol_bc_jb
        write(21, "(2i6,1e16.8)") n, m, J_b((n-1)*m_pol_bc_jb + m)
      enddo
    enddo  
    close(21)
  end subroutine solve_J_boundary_eqn
  
  subroutine solve_Psi_boundary_eqn(node_list, boundary_list)
  !---------------------------------------------------------------------------------------------------------------------------------
  ! The Psi boundary equation is solved using the Fourier-Galerkin method here
  !
  ! \int\int v*dPsi/dtheta*dtheta*dchi = -\int\int v*J'grad(s).grad(chi)*dtheta*dchi
  !
  ! where v = Z_m(theta)*Z_n(chi) is the test function and Psi = \sum_{m,n} C_{m,n}*Z_m(theta)*Z_n(chi)
  !
  ! The theta integration is over [0,2*pi] and the chi integration is over [chi_0+2*pi*F_0/N_p, chi_0], where chi_0 is the value of
  !  chi at the initial chi=const surface and N_p is the number of periods. Note that the limits of the chi integration is from 
  !  high to low values, because chi was calculated using a lefthand coordinate system, while JOREK is righthanded. Chi is therefore
  !  decreasing along the integration path.
  !---------------------------------------------------------------------------------------------------------------------------------
    use constants, only: pi
    use mod_parameters, only: n_period, n_plane, n_order, n_tor, n_coord_tor
    use phys_module, only: F0, m_pol_bc
    use gauss
    use basis_at_gaussian
    use data_structure
    use mod_chi
    use mod_SolveMN
    implicit none

    type(type_node_list),        intent(in) :: node_list
    type(type_bnd_element_list), intent(in) :: boundary_list

    real*8  :: phi, theta, element_size_ij, element_size_perp, BigR, grad_chi(3), Jgrad_ps(3)
    integer :: N_tht, mp, ms, ielm, i, j, j2, n, m, nn, mm, in, ind1, ind2, info

    type(type_bnd_element) :: bnd_element
    type(type_node)        :: node

    real*8, dimension(:,:,:), allocatable :: x_g, x_s, x_p, y_g, y_s, y_p
    real*8, dimension(:,:),   allocatable :: Amat
    real*8, dimension(:),     allocatable :: RHS

    real*8, dimension(0:n_order-1,0:n_order-1,0:n_order-1) :: chi

    write(*,*) "***************************************"
    write(*,*) "*    Solving boundary Psi equation    *"
    write(*,*) "***************************************"
    write(*,'(A,I4,A,I4)') "Toroidal modes: ", n_tor
    write(*,'(A,I4,A,I4)') "Poloidal modes: ", m_pol_bc

    N_tht = boundary_list%n_bnd_elements

    allocate(x_g(n_plane,n_gauss,N_tht)); allocate(x_s, x_p, mold=x_g)
    allocate(y_g(n_plane,n_gauss,N_tht)); allocate(y_s, y_p, mold=y_g)
    allocate(Amat(n_tor*m_pol_bc,n_tor*m_pol_bc))
    allocate(RHS(n_tor*m_pol_bc))

    x_g  = 0.d0; x_s = 0.d0; x_p = 0.d0
    y_g  = 0.d0; y_s = 0.d0; y_p = 0.d0

    do ielm=1,N_tht
      bnd_element = boundary_list%bnd_element(ielm)

      do i=1,2
        node = node_list%node(bnd_element%vertex(i))

        do j=1,2
          j2 = bnd_element%direction(i,j)
          element_size_ij = bnd_element%size(i,j)

          do ms=1,n_gauss
            do mp=1,n_plane
              do in=1,n_coord_tor
                x_g(mp,ms,ielm)  = x_g(mp,ms,ielm)  + node%x(in,j2,1)*element_size_ij*H1(i,j,ms)   *HZ_coord(in,mp)
                x_s(mp,ms,ielm)  = x_s(mp,ms,ielm)  + node%x(in,j2,1)*element_size_ij*H1_s(i,j,ms) *HZ_coord(in,mp)
                x_p(mp,ms,ielm)  = x_p(mp,ms,ielm)  + node%x(in,j2,1)*element_size_ij*H1(i,j,ms)   *HZ_coord_p(in,mp)

                y_g(mp,ms,ielm)  = y_g(mp,ms,ielm)  + node%x(in,j2,2)*element_size_ij*H1(i,j,ms)   *HZ_coord(in,mp)
                y_s(mp,ms,ielm)  = y_s(mp,ms,ielm)  + node%x(in,j2,2)*element_size_ij*H1_s(i,j,ms) *HZ_coord(in,mp)
                y_p(mp,ms,ielm)  = y_p(mp,ms,ielm)  + node%x(in,j2,2)*element_size_ij*H1(i,j,ms)   *HZ_coord_p(in,mp)
              end do
            end do
          end do
        end do
      end do
    end do

    Amat = 0.d0; RHS = 0.d0

    ! Integration of the RHS: it is much easier to integrate over one period than over the manifold [0,2*pi]x[chi_0,chi_0+2*pi*F_0/N_p]
    ! Due to periodicity, integration over the manifold above is equivalent to integration over one period
    ! The Jacobian for dtheta*dchi -> dtheta*dphi is dchi/dphi
    do ielm=1,N_tht
      do ms=1,n_gauss
        theta = 2.d0*pi*(float(ielm-1) + xgauss(ms))/float(N_tht)

        do mp=1,n_plane
          BigR = x_g(mp,ms,ielm)
          phi = 2.d0*pi*float(mp-1)/float(n_plane*n_period)
          chi = get_chi(x_g(mp,ms,ielm),y_g(mp,ms,ielm),phi)
          grad_chi = (/ chi(1,0,0), chi(0,1,0), chi(0,0,1)/BigR /)
          ! -e_theta x e_phi = -J*grad(psi)
          Jgrad_ps = (/ -BigR*y_s(mp,ms,ielm), BigR*x_s(mp,ms,ielm), x_p(mp,ms,ielm)*y_s(mp,ms,ielm) - x_s(mp,ms,ielm)*y_p(mp,ms,ielm) /)
          Jgrad_ps = Jgrad_ps*N_tht/(2.d0*pi) ! Multiply by dt/dtheta since J = 1/(grad(psi).(grad(theta)xgrad(phi))) and (s,t,phi) CS used above
          ! Note that J'*dchi/dphi = J, where J' = 1/(grad(psi).(grad(theta)xgrad(chi))) and dchi/dphi is from the switch dtheta*dchi -> dtheta*dphi

          ind1 = 1
          do n=1,n_tor
            do m=1,m_pol_bc
              RHS(ind1) = RHS(ind1) + dot_product(Jgrad_ps,grad_chi)*Zn(n,chi(0,0,0))*Zm(m,theta)*wgauss(ms)
              ind1 = ind1 + 1
            end do
          end do
        end do
      end do
    end do
    RHS = RHS*(2.d0*pi)**2/float(n_period*n_plane*N_tht)

    ! The integrals in the matrix are done analytically since they only contain sines and cosines
    ! Here, it is better to integrate over the original manifold
    ind1 = 1
    ! Loop over test functions (rows in the matrix)
    do n=1,n_tor
      do m=1,m_pol_bc
        ind2 = 1
        ! Loop over basis functions in the 2D Fourier decomposition of Psi (columns in the matrix)
        do nn=1,n_tor
          do mm=1,m_pol_bc
            if (n .eq. nn) then ! The test and basis functions are orthogonal if their chi indices are different (no chi derivatives)
              if (mod(m,2) .eq. 1 .and. mm .eq. m + 1) then ! If theta t.f. is cos, theta b.f. must be sin w/ same mode number (derivative is cos)
                Amat(ind1,ind2) = -pi*(mm/2)
              else if (mod(m,2) .eq. 0 .and. mm .eq. m - 1) then ! If theta t.f. is sin, theta b.f. must be cos w/ same mode number (derivative is sin)
                Amat(ind1,ind2) =  pi*(m/2)
              end if

              ! Multiply by the result of integration over chi; additional factor of 2 if both chi t.f. and b.f. are 1
              if (n .eq. 1) then
                Amat(ind1,ind2) = Amat(ind1,ind2)*2.d0*pi*F0/n_period
              else
                Amat(ind1,ind2) = Amat(ind1,ind2)*pi*F0/n_period
              end if
            end if

            ind2 = ind2 + 1
          end do
        end do

        ind1 = ind1 + 1
      end do
    end do

    call SolveMN(Amat, info)
    if (info .ne. 0) then
      write(*,'(A,I3)') "ERROR: solve_Psi_boundary_eqn, matrix inversion failed: ", info
      stop
    end if

    if (allocated(Psi_b)) deallocate(Psi_b)
    allocate(Psi_b(n_tor*m_pol_bc))

    Psi_b = matmul(Amat,RHS)

    ! Output coefficients as a diagnostic
    open(21,file='psi_boundary_coeff.dat')
    do n=1,n_tor
      do m=1,m_pol_bc
        write(21, "(2i6,1e16.8)") n, m, Psi_b((n-1)*m_pol_bc + m)
      enddo
    enddo  
    close(21)
  end subroutine solve_Psi_boundary_eqn

  subroutine setup_boundary_condition(node_list, bnd_node_list, boundary_list)
    use constants, only: pi
    use mod_parameters, only: n_tor, n_order, n_plane, n_period, n_coord_tor, var_Psi
    use phys_module, only: F0, m_pol_bc
    use basis_at_gaussian
    use gauss
    use data_structure
    use mod_chi
    implicit none
    
    type(type_node_list),     intent(inout) :: node_list
    type(type_bnd_node_list), intent(in)    :: bnd_node_list
    type(type_bnd_element_list), intent(in) :: boundary_list
    
    type(type_bnd_element) :: bnd_element
    type(type_node)        :: node
    
    integer :: ibnd_elm, ielm, i, j, j2, ms, mp, im, in, n, m, ind, N_tht
    real*8  :: R, z, R_s, z_s, phi, theta, Psi, Psi_tht, Psi_chi, chi_tht, zj, zj_tht, element_size_ij
    
    ! Dofs for Psi and zj at a particular boundary node -- values of variable and its poloidal (tangential) derivative
    real*8, dimension(n_tor) :: Psi_dof, Psi_dof2
    real*8, dimension(n_tor) :: zj_dof, zj_dof2
    real*8, dimension(:,:,:), allocatable  :: zj_g, x_g, y_g
    
    real*8, dimension(0:n_order-1,0:n_order-1,0:n_order-1) :: chi
    
    N_tht = boundary_list%n_bnd_elements
    allocate(zj_g(n_plane,n_gauss,N_tht))
    allocate(x_g(n_plane,n_gauss,N_tht)); allocate(y_g, mold=x_g)
    zj_g = 0.0  

    do i=1,bnd_node_list%n_bnd_nodes
      in = bnd_node_list%bnd_node(i)%index_jorek
      j  = bnd_node_list%bnd_node(i)%direction(2)
      theta = 2.d0*pi*float(i-1)/float(bnd_node_list%n_bnd_nodes)
      Psi_dof = 0.d0; Psi_dof2 = 0.d0
      zj_dof = 0.d0; zj_dof2 = 0.d0
      zj = 0.d0;  zj_tht = 0.d0
      
      do mp=1,n_plane
        phi = 2.d0*pi*float(mp-1)/float(n_plane*n_period)

        R = 0.d0; R_s = 0.d0
        z = 0.d0; z_s = 0.d0
        Psi = 0.d0; Psi_tht = 0.d0; Psi_chi = 0.d0
        
        ! Find R, z and their poloidal derivatives at the current boundary node and poloidal plane
        do im=1,n_coord_tor
          R   = R   + node_list%node(in)%x(im,1,1)*HZ_coord(im,mp)
          R_s = R_s + node_list%node(in)%x(im,j,1)*HZ_coord(im,mp)*3.d0 ! Poloidal derivative at node is 3x the corresponding dof of that node
          z   = z   + node_list%node(in)%x(im,1,2)*HZ_coord(im,mp)
          z_s = z_s + node_list%node(in)%x(im,j,2)*HZ_coord(im,mp)*3.d0
        end do
        
        chi = get_chi(R,z,phi)
        chi_tht = (R_s*chi(1,0,0) + z_s*chi(0,1,0))*bnd_node_list%n_bnd_nodes/(2.d0*pi)
        
        ! Using the (theta,chi) 2D Fourier basis, calculate the values of Psi and its derivatives at the current point
        ind = 1
        do n=1,n_tor
          do m=1,m_pol_bc
            Psi = Psi + Psi_b(ind)*Zn(n,chi(0,0,0))*Zm(m,theta)
            Psi_tht = Psi_tht + Psi_b(ind)*Zn(n,chi(0,0,0))*Zm_tht(m,theta)
            Psi_chi = Psi_chi + Psi_b(ind)*Zn_chi(n,chi(0,0,0))*Zm(m,theta)
            ind = ind + 1
          end do
        end do
        
        ! Now integrate Psi and its theta derivative (in the (psi,theta,phi) CS) over phi
        do im=1,n_tor
          Psi_dof(im) = Psi_dof(im) + Psi*HZn(im,phi)
          Psi_dof2(im) = Psi_dof2(im) + (Psi_tht + Psi_chi*chi_tht)*HZn(im,phi)
        end do
      end do
      
      Psi_dof2 = Psi_dof2*2.d0*pi/bnd_node_list%n_bnd_nodes ! Multiply by dtheta/dt to get derivative wrt element local coordinate
      
      ! Combining two steps here (multiplying by delta_phi and dividing by integrals of 1, sin^2, cos^2 to get modes), several factors cancel
      Psi_dof(1) = Psi_dof(1)/n_plane
      Psi_dof2(1) = Psi_dof2(1)/(3.d0*n_plane) ! divide the derivative by 3 to get the Bezier dof
      if (n_tor .gt. 1) then
        Psi_dof(2:) = Psi_dof(2:)*2.d0/n_plane
        Psi_dof2(2:) = Psi_dof2(2:)*2.d0/(3.d0*n_plane)
      end if
      
      ! Multiply by F0 so that the internally stored Psi matches the Grad-Shafranov solution in the tokamak limit
      node_list%node(in)%values(:,1,var_Psi) = F0*Psi_dof
      node_list%node(in)%values(:,j,var_Psi) = F0*Psi_dof2

      ! Compute zj degrees of freedom - as these have been calculated in (psi, theta, phi) space, the contributions of the poloidal harmonics just needs to be summed
      ind = 1
      do im=1, n_tor
        zj = 0.0
        zj_tht = 0.0
        do m=1,m_pol_bc+1
          zj     = zj + J_b(ind)*Zjm(m,theta)
          zj_tht = zj_tht + J_b(ind)*Zjm_tht(m,theta)
          ind = ind + 1
        enddo
        node_list%node(in)%values(im,1,var_zj) = zj
        node_list%node(in)%values(im,j,var_zj) = zj_tht*2.d0*pi/bnd_node_list%n_bnd_nodes/3.d0
      enddo

    end do

    ! Get a diagnostic of the current on the boundary in JOREK notation
    do ielm=1,N_tht
      bnd_element = boundary_list%bnd_element(ielm)
      write(*,*) ielm
      do i=1,2
        node = node_list%node(bnd_element%vertex(i))
      
        do j=1,2
          j2 = bnd_element%direction(i,j)
          element_size_ij = bnd_element%size(i,j)
      
          do ms=1,n_gauss
            do mp=1,n_plane
              do in=1,n_tor
                zj_g(mp,ms,ielm)  = zj_g(mp,ms,ielm)  + node%values(in,j2,var_zj)*element_size_ij*H1(i,j,ms)   *HZ(in,mp)
              end do
            end do
          end do
        end do
      end do
    end do
    write(*,*) 'Writing output'
    open(21, file='j_boundary_points_JOREK.dat', action='write', status='replace')
    do ibnd_elm=1,N_tht
      
      do ms=1,n_gauss
        theta = 2.d0*pi*(float(ibnd_elm-1) + xgauss(ms))/float(N_tht)
    
        do mp=1,n_plane
          phi = 2.d0*pi*float(mp-1)/float(n_plane*n_period)

          write(21, *) theta, phi, zj_g(mp, ms, ibnd_elm)
        enddo
      enddo
    enddo
    close(21)

  end subroutine setup_boundary_condition

  pure real*8 function Zn(n, chi)
    use mod_parameters, only: n_period
    use phys_module, only: F0
    implicit none
    integer, intent(in) :: n
    real*8,  intent(in) :: chi

    if (n .eq. 1) then
      Zn = 1.d0
      return
    end if

    if (mod(n,2) .eq. 0) then
      Zn = cos(n_period*int(n/2)*chi/F0)
    else
      Zn = sin(n_period*int(n/2)*chi/F0)
    end if
  end function Zn

  pure real*8 function Zn_chi(n, chi)
    use mod_parameters, only: n_period
    use phys_module, only: F0
    implicit none
    integer, intent(in) :: n
    real*8,  intent(in) :: chi

    if (n .eq. 1) then
      Zn_chi = 0.d0
      return
    end if

    if (mod(n,2) .eq. 0) then
      Zn_chi = -n_period*int(n/2)*sin(n_period*int(n/2)*chi/F0)/F0
    else
      Zn_chi =  n_period*int(n/2)*cos(n_period*int(n/2)*chi/F0)/F0
    end if
  end function Zn_chi
  
  pure real*8 function Zjm(m, theta)
    implicit none
    integer, intent(in) :: m
    real*8,  intent(in) :: theta
    
    if (m .eq. 1) then
      Zjm = 1.d0
      return
    end if

    if (mod(m,2) .eq. 0) then
      Zjm = cos(int(m/2)*theta)
    else
      Zjm = sin(int(m/2)*theta)
    end if
  end function Zjm
  
  pure real*8 function Zjm_tht(m, theta)
    implicit none
    integer, intent(in) :: m
    real*8,  intent(in) :: theta
    
    if (m .eq. 1) then
      Zjm_tht = 0.d0
      return
    end if

    if (mod(m,2) .eq. 0) then
      Zjm_tht = -int(m/2)*sin(int(m/2)*theta)
    else
      Zjm_tht = int(m/2)*cos(int(m/2)*theta)
    end if
  end function Zjm_tht
  
  pure real*8 function Zm(m, theta)
    implicit none
    integer, intent(in) :: m
    real*8,  intent(in) :: theta

    if (mod(m,2) .eq. 1) then
      Zm = cos(int((m+1)/2)*theta)
    else
      Zm = sin(int((m+1)/2)*theta)
    end if
  end function Zm

  pure real*8 function Zm_tht(m, theta)
    implicit none
    integer, intent(in) :: m
    real*8,  intent(in) :: theta

    if (mod(m,2) .eq. 1) then
      Zm_tht = -int((m+1)/2)*sin(int((m+1)/2)*theta)
    else
      Zm_tht =  int((m+1)/2)*cos(int((m+1)/2)*theta)
    end if
  end function Zm_tht

  pure real*8 function Zm_tht_2(m, theta)
    implicit none
    integer, intent(in) :: m
    real*8,  intent(in) :: theta

    if (mod(m,2) .eq. 1) then
      Zm_tht_2 = -int((m+1)/2)**2*cos(int((m+1)/2)*theta)
    else
      Zm_tht_2 = -int((m+1)/2)**2*sin(int((m+1)/2)*theta)
    end if
  end function Zm_tht_2

  pure real*8 function HZn(n, phi)
    use mod_parameters, only: n_period
    implicit none
    integer, intent(in) :: n
    real*8,  intent(in) :: phi

    if (n .eq. 1) then
      HZn = 1.d0
      return
    end if

    if (mod(n,2) .eq. 0) then
      HZn = cos(n_period*int(n/2)*phi)
    else
      HZn = sin(n_period*int(n/2)*phi)
    end if
  end function HZn
end module mod_boundary_conditions

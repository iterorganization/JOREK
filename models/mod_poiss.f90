module mod_poiss
contains
subroutine Poisson(my_id,itype,node_list,element_list,bnd_node_list,bnd_elm_list,   &
                 ivar_in,ivar_out,i_harm, psi_axis,psi_bnd,xpoint,xcase,Z_xpoint, &
                 freeboundary_equil,refinement,iter)
!-------------------------------------------------------------------------------
! collect the element matrices into one large sparse matrix in coordinate format
!-------------------------------------------------------------------------------
use tr_module 
use data_structure
use phys_module, only: amix, amix_freeb, delta_psi_GS, newton_GS_freebnd, newton_GS_fixbnd, &
                       n_limiter, treat_axis, fix_axis_nodes, &
                       use_mumps_eq, use_pastix_eq, use_strumpack_eq
use equil_info,  only: ES
use vacuum_equilibrium, only: vacuum_equil
use mpi_mod
use mod_interp
use mod_basisfunctions
use mod_integer_types
use mod_node_indices

use mod_axis_treatment
#ifdef USE_PASTIX6
use mod_pastix
#endif
use mod_sparse_data, only: type_SP_SOLVER, mumps, pastix, strumpack
use mod_sparse, only: solve_sparse_system

implicit none

! --- Routine parameters
integer,                  intent(in)    :: my_id             ! MPI id
integer,                  intent(in)    :: itype             ! selects the physics model (GS, Laplace)
                                                             ! -1: GS_perturbation
                                                             ! -2: GS_inverse
                                                             ! +2: Poisson_inverse
                                                             ! +1 or +3: Poisson
                                                             ! 0: variable projection
type (type_node_list),    intent(inout) :: node_list
type (type_element_list), intent(inout) :: element_list
integer,                  intent(in)    :: ivar_in           ! index of the input variable
integer,                  intent(in)    :: ivar_out          ! index of the output variable
integer,                  intent(in)    :: i_harm            ! index of toroidal harmonic
integer,                  intent(in)    :: iter              ! the iteration number
integer,                  intent(in)    :: xcase              
logical,                  intent(in)    :: xpoint            
logical,                  intent(in)    :: freeboundary_equil
logical,                  intent(in)    :: refinement       

! --- Local variables
type (type_element)      :: element
type (type_node)         :: nodes(n_vertex_max)
type (type_element)      :: element_father
type (type_node)         :: nodes_father(n_vertex_max)
type (type_bnd_node_list)    :: bnd_node_list
type (type_bnd_element_list) :: bnd_elm_list

real*8   :: ELM(n_vertex_max*n_degrees,n_vertex_max*n_degrees), RHS(n_vertex_max*n_degrees)
real*8   :: ELM_axis(n_vertex_max*n_degrees,n_vertex_max*n_degrees),  ELM_bnd(n_vertex_max*n_degrees,n_vertex_max*n_degrees)
real*8   :: zbig, Z_xpoint(2), psi_axis, psi_bnd, psi_xpoint(2), R_xpoint(2), s_xpoint(2), t_xpoint(2)
real*8   :: R_axis, Z_axis, s_axis, t_axis
real*8   :: psi_axis_kl(n_vertex_max,n_degrees), psi_bnd_kl(n_vertex_max,n_degrees) 
real*8   :: amix_used
real*8   :: psi_lim, R_lim, Z_lim, R_out, Z_out, s_bnd, t_bnd, P_s,P_t,P_st,P_ss,P_tt
integer  :: i_elm_bnd, i_elm_axis, i_elm_xpoint(2), ifail
integer  :: n_AA, nz_AA, nz_AA_old, n_border, ilarge, ife, iv, i,j,k,l
integer  :: inode, index_large_i, knode, index_large_k, index_ij, index_kl, index, index_i
logical   :: newton_method_GS

real*8, dimension(4,n_degrees)   :: H, H_s, H_t, H_st
real*8, dimension(4,n_degrees)   :: G_axis, G_bnd, G_s, G_t, G_st, G_ss, G_tt
real*8                           :: lambda, mu
real*8                           :: Psi,dPsi_ds,dPsi_dt,d2Psi_dsdt
real*8                           :: dX_ds, dX_dt, dY_ds, dY_dt, d2X_dsdt, d2Y_dsdt, h_u, h_v, h_w
integer                          :: inode_father, Index_elm, i_father
integer, dimension(n_vertex_max) :: pr
integer, dimension(2)            :: parent
integer, dimension(n_vertex_max) :: node_out
integer:: nnz, ierr
integer*8 :: check_data
character*8 :: type
integer :: node_indices( (n_order+1)/2, (n_order+1)/2 )

type(type_SP_MATRIX) :: a_mat
type(type_RHS) :: rhs_vec, sol_vec
type(type_SP_SOLVER) :: solver
real*8 :: tmp

real*8 :: new_dofs(1:4), old_dofs(1:4)

if (my_id == 0) then
  write(*,*) '**************************************'
  write(*,*) '*            Poisson                 *'
  write(*,*) '**************************************'
  
  if (itype .eq. 0) then
    write(*,*)'*************************************'
    write(*,*)'*   Projection of Variable: ',ivar_out
    write(*,*)'*************************************'
  endif
  
  if (iter .le. 1) then
    write(*,*) ' i_type       : ',itype
    write(*,*) ' n_elements   : ',element_list%n_elements
    write(*,*) ' n_nodes      : ',node_list%n_nodes
    write(*,*) ' freeboundary_equil : ',freeboundary_equil
  endif

  newton_method_GS = newton_GS_fixbnd
  if (freeboundary_equil) newton_method_GS = newton_GS_freebnd
 
  if (newton_method_GS .and. (itype==-1)) then  
    nz_AA = 3 * element_list%n_elements * (n_vertex_max * n_degrees)**2  !factor 3 comes from axis and x-point contributions 
  else
    nz_AA = 1 * element_list%n_elements * (n_vertex_max * n_degrees)**2  
  endif

  call tr_debug_write("Deb_poisson",nz_AA)
  
  n_border = 0
  if (itype .ne. 0) then
    do i=1,node_list%n_nodes
      if(treat_axis)then
        ! --- Only one fixed for fixed-axis (only valid for G1-cases at the moment!!!)
        if (node_list%node(i)%axis_node    ) n_border = n_border+1
      else
        ! --- t-derivatives and cross derivatives are switched off on axis, so (n_order+1)/2 are not fixed
        if (node_list%node(i)%axis_node    ) n_border = n_border + n_degrees - (n_order+1)/2
      endif    
      ! --- on non-corner boundaries, only tangent derivatives are fixed, ie. (n_order+1)/2
      if (node_list%node(i)%boundary .eq. 1) n_border = n_border + (n_order+1)/2
      if (node_list%node(i)%boundary .eq. 2) n_border = n_border + (n_order+1)/2
      if (node_list%node(i)%boundary .eq. 4) n_border = n_border + (n_order+1)/2
      if (node_list%node(i)%boundary .eq. 5) n_border = n_border + (n_order+1)/2
      if (node_list%node(i)%boundary .eq.11) n_border = n_border + (n_order+1)/2
      if (node_list%node(i)%boundary .eq.12) n_border = n_border + (n_order+1)/2
      if (node_list%node(i)%boundary .eq.15) n_border = n_border + (n_order+1)/2
      ! --- on corner boundaries, derivatives in both are fixed, but not cross-derivatives (-1 is to avoid having value twice)
      if (node_list%node(i)%boundary .eq. 3) n_border = n_border + 2 * (n_order+1)/2 - 1
      if (node_list%node(i)%boundary .eq. 9) n_border = n_border + 2 * (n_order+1)/2 - 1
      if (node_list%node(i)%boundary .eq.19) n_border = n_border + 2 * (n_order+1)/2 - 1
      if (node_list%node(i)%boundary .eq.20) n_border = n_border + 2 * (n_order+1)/2 - 1
      if (node_list%node(i)%boundary .eq.21) n_border = n_border + 2 * (n_order+1)/2 - 1
    enddo
  endif
  
  if ((.not. freeboundary_equil) .or. (itype .ne. -1)) then
    nz_AA = nz_AA + n_border
  elseif  (freeboundary_equil .and. (itype .eq. -1)) then
    nz_AA = nz_AA + 128 * bnd_node_list%n_bnd_nodes**2
  endif
    
  n_AA = 0
  do inode = 1, node_list%n_nodes
      do k = 1, n_degrees 
        n_AA = max(n_AA,node_list%node(inode)%index(k))
      enddo
  enddo
  
  if (iter .le. 1) then
    write(*,*) ' number of unknowns      : ',n_AA, node_list%n_nodes * n_degrees
    write(*,*) ' number of boundary nodes: ',n_border
    write(*,*) ' nz_AA                   : ',nz_AA
  endif

  a_mat%ng  = n_AA
  a_mat%nnz  = nz_AA
  a_mat%comm = MPI_COMM_SELF

  if (associated(a_mat%irn)) call tr_deallocatep(a_mat%irn,"a_mat_eq",CAT_DMATRIX)
  if (associated(a_mat%jcn)) call tr_deallocatep(a_mat%jcn,"a_mat_eq",CAT_DMATRIX)
  if (associated(a_mat%val)) call tr_deallocatep(a_mat%val,"a_mat_eq",CAT_DMATRIX)

  call tr_allocatep(a_mat%val,int1,a_mat%nnz,"a_mat_eq",CAT_DMATRIX)
  call tr_allocatep(a_mat%irn,int1,a_mat%nnz,"a_mat_eq",CAT_DMATRIX)
  call tr_allocatep(a_mat%jcn,int1,a_mat%nnz,"a_mat_eq",CAT_DMATRIX)

  a_mat%irn(1:a_mat%nnz) = 0
  a_mat%jcn(1:a_mat%nnz) = 0
  a_mat%val(1:a_mat%nnz) = 0.d0

  rhs_vec%n = n_AA

  if (associated(rhs_vec%val)) call tr_deallocatep(rhs_vec%val,"rhs_eq",CAT_DMATRIX)
  call tr_allocatep(rhs_vec%val,int1,rhs_vec%n,"rhs_eq",CAT_DMATRIX)
  rhs_vec%val(1:rhs_vec%n) = 0.d0

  ilarge=0
  
  amix_used = amix
  
  if (itype .eq. -1) then     !--- if solving Grad-Shafranov
    
   
    !-- Find axis and x-poit matrix contributions (for Newton iterations)   
    call basisfunctions(ES%s_axis, ES%t_axis, G_axis, G_s, G_t, G_st, G_ss, G_tt)
    call basisfunctions(ES%s_bnd ,  ES%t_bnd,  G_bnd, G_s, G_t, G_st, G_ss, G_tt)
    do k=1,n_vertex_max  
      do l=1,n_degrees
        psi_axis_kl(k,l) =  G_axis(k,l) * element_list%element(ES%i_elm_axis)%size(k,l)  !--- matrix contributions of axis dofs
        psi_bnd_kl(k,l)  =  G_bnd(k,l)  * element_list%element(ES%i_elm_bnd)%size(k,l)   !--- matrix contributions of bnd point dofs
      enddo
    enddo

    ! --- ATTENTION: limiter free-boundary plasmas not working well yet with Newton method 
    if (.not. xpoint) psi_bnd_kl = 0.d0  
 
    if (.not. newton_method_GS) then
      psi_bnd_kl  = 0.d0
      psi_axis_kl = 0.d0
    endif
    
    if (freeboundary_equil) amix_used = amix_freeb
  
  endif   !--- end type -1 (GS equilibrium)
  
  do ife =1, element_list%n_elements
  
    element = element_list%element(ife)
    
    if (refinement) then                  ! no contribution from elements which have children
      if (element%n_sons .ne. 0) cycle
    endif
    
    if (refinement) then
      i_father= element_list%element(ife)%father
      if (i_father.ne. 0) then
        element_father = element_list%element(i_father)
      endif
  
      do iv = 1, n_vertex_max
  
        if (i_father.ne.0) then
          inode_father=element_father%vertex(iv)
          nodes_father(iv) = node_list%node(inode_father)
        endif
      enddo
    endif ! refinement
    
    do iv = 1, n_vertex_max
      inode     = element%vertex(iv)
      nodes(iv) = node_list%node(inode)
    enddo
  
    if (itype .eq. -1) then
      
      call element_matrix_GS_perturbation(xpoint,xcase,Z_xpoint,psi_axis,psi_bnd,element,nodes,ivar_in,ivar_out,i_harm,ELM,RHS, &
      psi_axis_kl, ELM_axis, psi_bnd_kl, ELM_bnd, newton_method_GS)
      
    elseif (itype .eq. -2) then
  
      call element_matrix_GS_inverse(xpoint,xcase,Z_xpoint,psi_axis,psi_bnd,element,nodes,ivar_in,ivar_out,i_harm,ELM,RHS)
  
    elseif (itype .eq. +2) then
  
      call element_matrix_Poisson_inverse(itype,element,nodes,ivar_in,ivar_out,i_harm,ELM,RHS)
  
    elseif (itype .eq. 0) then
  
      call element_matrix_projection(itype,element,nodes,ivar_in,ivar_out,i_harm,ELM,RHS)
  
    else
  
      call element_matrix_Poisson(itype,element,nodes,ivar_in,ivar_out,i_harm,ELM,RHS)
  
    endif  

    ! Transform basis functions for the axis nodes. This will solve for new degrees of freedom at the axis.    
    if( (treat_axis) .and. (nodes(1)%axis_node .or. nodes(2)%axis_node .or. nodes(3)%axis_node .or. nodes(4)%axis_node) ) then
      call transform_basis_for_axis_element_poisson(nodes, ELM, RHS, ivar_in, ivar_out, i_harm)
    endif
    
    if (refinement) then ! Processing  "constrained nodes"
      call Chgmt_node(ife,element,nodes,element_father,nodes_father,ELM,RHS,node_out) 
    else
      node_out = element%vertex
    endif
    
    do i=1,n_vertex_max
  
      inode = node_out(i)
  
      do j=1,n_degrees
  
        index_ij = (i-1)*n_degrees + j     ! index in the ELM matrix
  
        index_large_i = node_list%node(inode)%index(j)  ! base index in the main matrix
  
        rhs_vec%val(index_large_i) = rhs_vec%val(index_large_i) + RHS(index_ij)
  
        do k=1,n_vertex_max
  
          knode         =node_out(k)! element%vertex(k)
  
          do l=1,n_degrees
  
            index_kl = (k-1)*n_degrees + l
  
            index_large_k = node_list%node(knode)%index(l)  ! base index in the main matrix
  
            ilarge = ilarge +1
  
            a_mat%irn(ilarge) = index_large_i
            a_mat%jcn(ilarge) = index_large_k
            a_mat%val(ilarge) = ELM(index_ij,index_kl)
  
          enddo
        enddo

        if (newton_method_GS .and. (itype==-1)) then     ! newton method extra contributions
        
          ! --- Perturbed contribution of the magnetic axis
          do k=1,n_vertex_max
    
            knode = element_list%element(ES%i_elm_axis)%vertex(k)
    
            do l=1,n_degrees
    
              index_kl = (k-1)*n_degrees + l
    
              index_large_k = node_list%node(knode)%index(l)  ! base index in the main matrix
    
              ilarge = ilarge +1
    
              a_mat%irn(ilarge) = index_large_i
              a_mat%jcn(ilarge) = index_large_k
              a_mat%val(ilarge)   = ELM_axis(index_ij,index_kl)
    
            enddo
          enddo
          
          ! --- Perturbed contribution of the limiter/X-point
          do k=1,n_vertex_max
    
            knode = element_list%element(ES%i_elm_bnd)%vertex(k)
    
            do l=1,n_degrees
    
              index_kl = (k-1)*n_degrees + l
    
              index_large_k = node_list%node(knode)%index(l)  ! base index in the main matrix
    
              ilarge = ilarge +1
    
              a_mat%irn(ilarge) = index_large_i
              a_mat%jcn(ilarge) = index_large_k
              a_mat%val(ilarge)   = ELM_bnd(index_ij,index_kl)           
    
            enddo
          enddo

        endif ! newton method extra contributions
        
      enddo
    enddo
  
  enddo

  nz_AA_old = nz_AA
  nz_AA = ilarge
  a_mat%nnz = nz_AA
  
  zbig = 1.d10
  
end if ! my_id == 0

!----------------------- boundary conditions

if (freeboundary_equil .and. (itype .eq. -1)) then
  
  call vacuum_equil(my_id,node_list,bnd_node_list,bnd_elm_list,psi_axis,psi_bnd, a_mat, rhs_vec)
  
elseif (itype .ne. 0) then        ! apply fixed boundary conditions (not for variable projection)

  if (my_id == 0 ) then

    ! --- calculate node_indices
    call calculate_node_indices(node_indices)

    do i=1,node_list%n_nodes
  
      ! --- On axis, we fix the t-derivatives, plus all cross-derivatives
      if (node_list%node(i)%axis_node) then
      
        if (treat_axis) then ! For G1 elements only at the moment !
          ! penalize 4th DoF to enforce C0 continuity at the grid center        
          index_i = node_list%node(i)%index(4)  ! base index in the main matrix
          a_mat%irn(ilarge+1) = index_i
          a_mat%jcn(ilarge+1) = index_i
          a_mat%val(ilarge+1)   = zbig
          ilarge = ilarge + 1
        endif

        if(fix_axis_nodes)then
          do k = 1,(n_order+1)/2
            do l = 2,(n_order+1)/2 ! start t-index from 2 to keep only the pure s-derivatives
              index = node_indices(k,l)
              index_i = node_list%node(i)%index(index)  ! base index in the main matrix
              a_mat%irn(ilarge+1) = index_i
              a_mat%jcn(ilarge+1) = index_i
              a_mat%val(ilarge+1)   = zbig
              ilarge = ilarge + 1
            enddo
          enddo
        endif

      endif

      if (node_list%node(i)%boundary .ne. 0) then
  
        ! --- fix node value (index is always 1)
        index_i = node_list%node(i)%index(1)  ! base index in the main matrix
        a_mat%irn(ilarge+1) = index_i
        a_mat%jcn(ilarge+1) = index_i
        a_mat%val(ilarge+1)   = zbig
        ilarge = ilarge + 1
           
        if (     (node_list%node(i)%boundary .eq. 1) &
            .or. (node_list%node(i)%boundary .eq. 3) &
            .or. (node_list%node(i)%boundary .eq. 4) &
            .or. (node_list%node(i)%boundary .eq. 9) &
            .or. (node_list%node(i)%boundary .eq.11) &
            .or. (node_list%node(i)%boundary .eq.12) &
            .or. (node_list%node(i)%boundary .eq.19) &
            .or. (node_list%node(i)%boundary .eq.20) &
            .or. (node_list%node(i)%boundary .eq.21) &
        ) then
  
          ! --- Fix s-derivatives
          do k = 2,(n_order+1)/2 ! start from 2 because node value already fixed above
            l = 1 ! t-index = 1 to fix only s-derivatives
            index = node_indices(k,l)
            index_i = node_list%node(i)%index(index)  ! base index in the main matrix
            a_mat%irn(ilarge+1) = index_i
            a_mat%jcn(ilarge+1) = index_i
            a_mat%val(ilarge+1)   = zbig
            ilarge = ilarge + 1
          enddo

        endif
  
        if (     (node_list%node(i)%boundary .eq. 2) &
            .or. (node_list%node(i)%boundary .eq. 3) &
            .or. (node_list%node(i)%boundary .eq. 5) &
            .or. (node_list%node(i)%boundary .eq. 9) &
            .or. (node_list%node(i)%boundary .eq.15) &
            .or. (node_list%node(i)%boundary .eq.19) &
            .or. (node_list%node(i)%boundary .eq.20) &
            .or. (node_list%node(i)%boundary .eq.21) &
        ) then
  
          ! --- Fix t-derivatives
          k = 1 ! s-index = 1 to fix only t-derivatives
          do l = 2,(n_order+1)/2 ! start from 2 because node value already fixed above
            index = node_indices(k,l)
            index_i = node_list%node(i)%index(index)  ! base index in the main matrix
            a_mat%irn(ilarge+1) = index_i
            a_mat%jcn(ilarge+1) = index_i
            a_mat%val(ilarge+1)   = zbig
            ilarge = ilarge + 1
          enddo
      
        endif

      endif
    enddo
  
    nz_AA_old = nz_AA
    nz_AA     = ilarge
    a_mat%nnz = nz_AA
 
  end if ! my_id == 0
  
endif

if (my_id == 0) then

  solver%equilibrium = .true.
  solver%verbose = .false.
  if (use_strumpack_eq) then
    solver%library = strumpack
  elseif (use_mumps_eq) then
    solver%library = mumps
  elseif (use_pastix_eq) then
    solver%library = pastix
  endif
  call solve_sparse_system(a_mat, rhs_vec, rhs_vec, solver)
  call solver%finalize()
  
  call tr_debug_write("a_mat%ng",int(a_mat%ng))
  call tr_debug_write("a_mat%nnz",int(a_mat%nnz))
  
  do i=1,node_list%n_nodes
  
    if ((.not. refinement) .or. (refinement .and. (.not. node_list%node(i)%constrained)) ) then
 
      ! We need to transform the new dof to old ones on the axis.
      ! The respective RHS entries on the axis (which are shared)
      ! are updated during the transformation. At the end of the loop,
      ! we recover RHS entries so that they same can be used for all the axis nodes.            
      if(treat_axis .and. node_list%node(i)%axis_node)then
        do k=1,n_degrees
          index = node_list%node(i)%index(k)
          new_dofs(k) = rhs_vec%val(index)
        enddo
        call new_to_old_dofs_on_the_axis(node_list, i, new_dofs, old_dofs)
        do k=1,n_degrees
          index = node_list%node(i)%index(k)
          rhs_vec%val(index) = old_dofs(k)
        enddo
      endif
            
      do k=1,n_degrees
  
        index = node_list%node(i)%index(k)
  
        !--------------- for equation in perturbation form
        if (itype .eq. -1) then
          node_list%node(i)%deltas(i_harm,k,ivar_out) = rhs_vec%val(index)
          node_list%node(i)%values(i_harm,k,ivar_out) = node_list%node(i)%values(i_harm,k,ivar_out) &
                                                      + (1.d0 - amix_used) * rhs_vec%val(index)
        !--------------- Variable projection
        elseif (itype .eq. 0) then
          if (ivar_out .eq. 710) then
#ifdef fullmhd
            node_list%node(i)%Fprof_eq(k) = node_list%node(i)%Fprof_eq(k) + (1.d0 - amix_used) * rhs_vec%val(index)
#endif
          else
            node_list%node(i)%values(1,k,ivar_out) = node_list%node(i)%values(1,k,ivar_out) + (1.d0 - amix_used) * rhs_vec%val(index)
          endif
        !--------------- for equation on total flux
        else
          node_list%node(i)%deltas(i_harm,k,ivar_out) = node_list%node(i)%values(i_harm,k,ivar_out) - rhs_vec%val(index)
          node_list%node(i)%values(i_harm,k,ivar_out) = amix_used * node_list%node(i)%values(i_harm,k,ivar_out) &
                                                      + (1.d0 - amix_used) * rhs_vec%val(index)
        endif
        
      enddo    ! order

      ! recover RHS entries.
      if(treat_axis .and. node_list%node(i)%axis_node)then
        do k=1,n_degrees
          index = node_list%node(i)%index(k)
          rhs_vec%val(index) = new_dofs(k)
        enddo
      endif
      
    endif      ! refinement, constrained
  enddo        ! nodes
  
  !*************************************************************************
  ! Solutions at constrained nodes                                         *
  !*************************************************************************
  if (refinement) then
 
    do i = 1, node_list%n_nodes
  
      if (node_list%node(i)%constrained) then
  
        lambda = node_list%node(i)%ref_lambda
        mu     = node_list%node(i)%ref_mu
        index_elm = node_list%node(i)%parent_elem
        parent(1) = node_list%node(i)%parents(1)
        parent(2) = node_list%node(i)%parents(2)
  
        call basisfunctions(lambda, mu, H, H_s, H_t, H_st)
  
        do j = 1, n_vertex_max
          pr(j) = element_list%element(index_elm)%vertex(j)
        end do
  
        dx_ds = 0.
        dx_dt = 0.
        dy_ds = 0.
        dy_dt = 0.
        d2x_dsdt = 0.
        d2y_dsdt = 0.
  
        Psi = 0.
        dPsi_ds = 0.
        dPsi_dt = 0.
        d2Psi_dsdt = 0.
  
        do k = 1, n_vertex_max
  
          if ((pr(k)==parent(1)).or.(pr(k)==parent(2))) then
  
            do l = 1, n_degrees
    
              dx_ds = dx_ds + node_list%node(pr(k))%x(1,l,1) * H_s(k,l) 	&
              * element_list%element(index_elm)%size(k,l)
        dx_dt = dx_dt + node_list%node(pr(k))%x(1,l,1) * H_t(k,l) 	&
                    * element_list%element(index_elm)%size(k,l)
  
              dy_ds = dy_ds + node_list%node(pr(k))%x(1,l,2) * H_s(k,l) 	&
                    * element_list%element(index_elm)%size(k,l)
        dy_dt = dy_dt + node_list%node(pr(k))%x(1,l,2) * H_t(k,l) 	&
              * element_list%element(index_elm)%size(k,l)
  
              d2x_dsdt = d2x_dsdt + node_list%node(pr(k))%x(1,l,1) * H_st(k,l) 	&
                       * element_list%element(index_elm)%size(k,l)
              d2y_dsdt = d2y_dsdt + node_list%node(pr(k))%x(1,l,2) * H_st(k,l) 	&
                       * element_list%element(index_elm)%size(k,l)
  
              Psi = Psi  + node_list%node(pr(k))%values(i_harm,l,ivar_out)*H(k,l)	   &
                  * element_list%element(index_elm)%size(k,l)
  
              dPsi_ds = dPsi_ds  + node_list%node(pr(k))%values(i_harm,l,ivar_out)*H_s(k,l)	   &
                      * element_list%element(index_elm)%size(k,l)
  
              dPsi_dt = dPsi_dt  + node_list%node(pr(k))%values(i_harm,l,ivar_out)*H_t(k,l)	   &
                      * element_list%element(index_elm)%size(k,l)
  
              d2Psi_dsdt = d2Psi_dsdt  + node_list%node(pr(k))%values(i_harm,l,ivar_out)*H_st(k,l)	   &
                * element_list%element(index_elm)%size(k,l)
            enddo
          endif
        enddo
  
        h_u = 1
        h_v = 1
        h_w = h_u*h_v
  
        node_list%node(i)%values(i_harm,1,ivar_out) = Psi
        node_list%node(i)%values(i_harm,2,ivar_out) = (dPsi_ds) /(3.*h_u)
        node_list%node(i)%values(i_harm,3,ivar_out) = (dPsi_dt) /(3.*h_v)
        node_list%node(i)%values(i_harm,4,ivar_out) = (d2Psi_dsdt) /(9.*h_w)
  
      endif   ! constrained
    enddo     ! nodes
  endif       ! refinement

  if (associated(a_mat%irn)) call tr_deallocatep(a_mat%irn,"a_mat_eq",CAT_DMATRIX)
  if (associated(a_mat%jcn)) call tr_deallocatep(a_mat%jcn,"a_mat_eq",CAT_DMATRIX)
  if (associated(a_mat%val)) call tr_deallocatep(a_mat%val,"a_mat_eq",CAT_DMATRIX)
  if (associated(rhs_vec%val)) call tr_deallocatep(rhs_vec%val,"rhs_eq",CAT_DMATRIX)

end if ! my_id == 0
  
return
end subroutine poisson

end module mod_poiss

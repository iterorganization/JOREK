!> Functionality for improving numerical stability at the grid center.
!! - Switched on by the logical namelist parameter treat_axis.
!! - When used, the dofs of all axis nodes are replaced by the following:
!!   1: value
!!   2: d/dR
!!   3: d/dZ
!!   4: d^2/dRdZ (actually not used, but forced to zero)
!! - For further details, refer to the JOREK wiki at https://www.jorek.eu/wiki/doku.php?id=grid-axis
module mod_axis_treatment

contains

! Transform basis functions for the axis nodes. This will solve for new degrees of freedom at the axis.
! This is the 2D version of the transform, useful for solving the Grad-Shafranov equation in mod_poiss.
! See the next subroutine, transform_basis_for_axis_element, if multiple Fourier modes need to be treated.
subroutine transform_basis_for_axis_element_poisson(nodes, ELM, RHS, ivar_in, ivar_out, i_harm)
use data_structure
use phys_module
implicit none
type(type_node),       intent(inout) :: nodes(n_vertex_max)
real*8,          intent(inout) :: ELM(1:n_vertex_max*(n_order+1), 1:n_vertex_max*(n_order+1))
real*8,          intent(inout) :: RHS(1:n_vertex_max*(n_order+1))
integer ,        intent(in)    :: ivar_in, ivar_out, i_harm
! --- routine parameters
integer :: axis_vertex1, axis_vertex4, dof1, dof2, dof3, dof4
integer :: iv, io, jv, jo, index_iv_io, index_jv_jo
real*8  :: Ptrans(1:4, 1:4), Pmat(1:4, 1:4), rhs_i(1:4), elm_ij(1:4, 1:4)

axis_vertex1 = 1 ; axis_vertex4 = 4
dof1 = 1 ; dof2 = 2 ; dof3 = 3 ; dof4 = 4

do iv = 1, n_vertex_max

  ! initialize transpose of the transformation matrix to identity
  Ptrans = 0.d0
  do io = 1, n_order+1
    Ptrans(io, io) = 1.d0
  enddo

  ! determine transpose of the transformation matrix
  if(iv==axis_vertex1 .or. iv==axis_vertex4)then
    if(nodes(iv)%axis_dof == dof2) then      
      Ptrans(dof2,dof3) = nodes(iv)%x(1,dof3,1) ; Ptrans(dof2,dof4) = nodes(iv)%x(1,dof4,1)
      Ptrans(dof3,dof3) = nodes(iv)%x(1,dof3,2) ; Ptrans(dof3,dof4) = nodes(iv)%x(1,dof4,2)
      Ptrans(dof2,dof2) = 0.d0
      Ptrans(dof4,dof4) = 0.d0
    elseif(nodes(iv)%axis_dof == dof3)then
      Ptrans(dof2,dof2) = nodes(iv)%x(1,dof2,1) ; Ptrans(dof2,dof4) = nodes(iv)%x(1,dof4,1)
      Ptrans(dof3,dof2) = nodes(iv)%x(1,dof2,2) ; Ptrans(dof3,dof4) = nodes(iv)%x(1,dof4,2)
      Ptrans(dof3,dof3) = 0.d0
      Ptrans(dof4,dof4) = 0.d0
    else
      Ptrans = 0.d0
      do io = 1, n_order+1
        Ptrans(io, io) = 1.d0
      enddo
    endif
  endif

  ! Extract RHS associated with a vertex iv: rhs_iv
  do io = 1, n_order+1
    index_iv_io = (iv-1)*(n_order+1) + io
    rhs_i(io)   = RHS(index_iv_io)
  enddo
  ! transform test functions for rhs_v
  rhs_i = matmul( Ptrans, rhs_i)
  ! fill the updated entries in RHS vector
  do io = 1, n_order+1
    index_iv_io = (iv-1)*(n_order+1) + io    
    RHS(index_iv_io) = rhs_i(io)
  enddo

  do jv = 1, n_vertex_max

    ! initialize the transformation matrix to identity  
    Pmat = 0.d0
    do jo = 1, n_order+1
      Pmat(jo, jo) = 1.d0
    enddo

    ! determine the transformation matrix    
    if(jv==axis_vertex1 .or. jv==axis_vertex4)then
      if(nodes(jv)%axis_dof == dof2)then
        Pmat(dof3,dof2)  = nodes(jv)%x(1,dof3,1) ; Pmat(dof3,dof3) = nodes(jv)%x(1,dof3,2)
        Pmat(dof4,dof2)  = nodes(jv)%x(1,dof4,1) ; Pmat(dof4,dof3) = nodes(jv)%x(1,dof4,2)
        Pmat(dof2,dof2)  = 0.d0
        Pmat(dof4,dof4)  = 0.d0 
      elseif(nodes(jv)%axis_dof == dof3)then
        Pmat(dof2,dof2)  = nodes(jv)%x(1,dof2,1) ; Pmat(dof2,dof3) = nodes(jv)%x(1,dof2,2)
        Pmat(dof4,dof2)  = nodes(jv)%x(1,dof4,1) ; Pmat(dof4,dof3) = nodes(jv)%x(1,dof4,2)
        Pmat(dof3,dof3)  = 0.d0
        Pmat(dof4,dof4)  = 0.d0
      else
        Pmat = 0.d0
        do jo = 1, n_order+1
          Pmat(jo, jo) = 1.d0
        enddo
      endif
    endif

    ! Extract ELM associated with a vertex iv and jv: elm_iv_jv
    do io = 1, n_order+1
       index_iv_io = (iv-1)*(n_order+1) + io
       do jo = 1, n_order+1
          index_jv_jo = (jv-1)*(n_order+1) + jo
          elm_ij(io, jo) =  ELM(index_iv_io, index_jv_jo)
       enddo
    enddo   
    ! transform test and basis functions for elm_iv_jv 
    elm_ij = matmul(matmul(PTrans, elm_ij), Pmat)
    ! fill the updated entries in ELM matrix
    do io = 1, n_order+1
       index_iv_io = (iv-1)*(n_order+1) + io
       do jo = 1, n_order+1
          index_jv_jo = (jv-1)*(n_order+1) + jo
          ELM(index_iv_io, index_jv_jo) = elm_ij(io, jo)
       enddo
    enddo

  enddo ! loop over jv

enddo ! loop over iv

end subroutine transform_basis_for_axis_element_poisson

! Transform basis functions for the axis nodes. This will solve for new degrees of freedom at the axis.
! This subroutine handles 3D transforms, i.e. when multiple Fourier modes are present.
! Parameters:
!   nodes:        An array containing all the nodes of the element whose ELM and RHS are being transformed
!   ielm:         The index in the element list of the element whose ELM and RHS are being transformed
!   ELM:          The element matrix that is being transformed
!   RHS:          The RHS that is being transformed
!   i_v:          An array containing the indices of all variables involved in the transformation, can be a subset of the variables in the matrix
!   n_v:          The total number of variables in the matrix, must be >= size(i_v)
!   i_n:          An array containing the indices of all Fourier modes involved in the transformation, can be a subset of the modes in the matrix
!   n_tor_local:  The total number of Fourier modes in the matrix, must be >= size(i_n)
subroutine transform_basis_for_axis_element(nodes, ielm, ELM, RHS, i_v, n_v, i_n, n_tor_local)
use data_structure
implicit none
type(type_node),       intent(inout) :: nodes(n_vertex_max)
integer, dimension(:), intent(in) :: i_v, i_n
real*8, dimension(n_tor_local*n_vertex_max*(n_order+1)*n_v, n_tor_local*n_vertex_max*(n_order+1)*n_v), intent(inout) :: ELM
real*8, dimension(n_tor_local*n_vertex_max*(n_order+1)*n_v), intent(inout) :: RHS
integer, intent(in) :: ielm, n_v, n_tor_local

call transform_rhs_for_axis_element(nodes, ielm, RHS, i_v, n_v, i_n, n_tor_local)
call transform_elm_for_axis_element(nodes, ielm, ELM, i_v, n_v, i_n, n_tor_local)

end subroutine transform_basis_for_axis_element

! Transform the test and basis functions of the element matrix for the axis nodes, ELM <- P^T ELM P.
! Parameters as in transform_basis_for_axis_element.
subroutine transform_elm_for_axis_element(nodes, ielm, ELM, i_v, n_v, i_n, n_tor_local)
use data_structure
implicit none
type(type_node),       intent(in) :: nodes(n_vertex_max)
integer, dimension(:), intent(in) :: i_v, i_n
real*8, dimension(n_tor_local*n_vertex_max*(n_order+1)*n_v, n_tor_local*n_vertex_max*(n_order+1)*n_v), intent(inout) :: ELM
integer, intent(in) :: ielm, n_v, n_tor_local

! --- routine parameters
integer :: iv, io, jv, jo, index_iv_io, index_jv_jo
integer :: ivar, jvar, im, in
real*8  :: Pmat(1:4, 1:4, n_vertex_max), Ptrans(1:4, 1:4), elm_ij(1:4, 1:4)

call axis_element_transform_matrices(nodes, ielm, Pmat)

do iv = 1, n_vertex_max

  Ptrans = transpose(Pmat(:,:,iv))

  do jv = 1, n_vertex_max

    ! Extract ELM associated with a vertex iv and jv: elm_iv_jv
    do ivar = 1, size(i_v)
    do im   = 1, size(i_n)
    
    do jvar = 1, size(i_v)
    do in   = 1, size(i_n)
    
       do io = 1, n_order+1
          index_iv_io = n_tor_local*n_v*(n_order+1)*(iv-1) + n_tor_local*n_v*(io-1) + n_tor_local*(i_v(ivar)-1) + (i_n(im)-i_n(1)+1)
          do jo = 1, n_order+1
             index_jv_jo = n_tor_local*n_v*(n_order+1)*(jv-1) + n_tor_local*n_v*(jo-1) + n_tor_local*(i_v(jvar)-1) + (i_n(in)-i_n(1)+1)
             elm_ij(io, jo) =  ELM(index_iv_io, index_jv_jo)
          enddo
       enddo
       ! transform test and basis functions for elm_iv_jv
       elm_ij = matmul(matmul(PTrans, elm_ij), Pmat(:,:,jv))
       ! fill the updated entries in ELM matrix
       do io = 1, n_order+1
          index_iv_io = n_tor_local*n_v*(n_order+1)*(iv-1) + n_tor_local*n_v*(io-1) + n_tor_local*(i_v(ivar)-1) + (i_n(im)-i_n(1)+1)
          do jo = 1, n_order+1
             index_jv_jo = n_tor_local*n_v*(n_order+1)*(jv-1) + n_tor_local*n_v*(jo-1) + n_tor_local*(i_v(jvar)-1) + (i_n(in)-i_n(1)+1)
             ELM(index_iv_io, index_jv_jo) = elm_ij(io, jo)
          enddo
       enddo

    enddo ! loop over in
    enddo ! loop over jvar
       
    enddo ! loop over im
    enddo ! loop over ivar

  enddo ! loop over jv

enddo ! loop over iv

end subroutine transform_elm_for_axis_element

! Transform the test functions of the element RHS for the axis nodes, RHS <- P^T RHS.
! Parameters as in transform_basis_for_axis_element.
subroutine transform_rhs_for_axis_element(nodes, ielm, RHS, i_v, n_v, i_n, n_tor_local)
use data_structure
implicit none
type(type_node),       intent(in) :: nodes(n_vertex_max)
integer, dimension(:), intent(in) :: i_v, i_n
real*8, dimension(n_tor_local*n_vertex_max*(n_order+1)*n_v), intent(inout) :: RHS
integer, intent(in) :: ielm, n_v, n_tor_local

! --- routine parameters
integer :: iv, io, index_iv_io, ivar, im
real*8  :: Pmat(1:4, 1:4, n_vertex_max), Ptrans(1:4, 1:4), rhs_i(1:4)

call axis_element_transform_matrices(nodes, ielm, Pmat)

do iv = 1, n_vertex_max

  Ptrans = transpose(Pmat(:,:,iv))

  do ivar = 1, size(i_v)
  do im   = 1, size(i_n)
    ! Extract RHS associated with a vertex iv: rhs_iv
    do io = 1, n_order+1
      index_iv_io = n_tor_local*n_v*(n_order+1)*(iv-1) + n_tor_local*n_v*(io-1) + n_tor_local*(i_v(ivar)-1) + (i_n(im)-i_n(1)+1)
      rhs_i(io)   = RHS(index_iv_io)
    enddo
    ! transform test functions for rhs_v
    rhs_i = matmul( Ptrans, rhs_i)
    ! fill the updated entries in RHS vector
    do io = 1, n_order+1
      index_iv_io = n_tor_local*n_v*(n_order+1)*(iv-1) + n_tor_local*n_v*(io-1) + n_tor_local*(i_v(ivar)-1) + (i_n(im)-i_n(1)+1)
      RHS(index_iv_io) = rhs_i(io)
    enddo
  enddo ! loop over im
  enddo ! loop over ivar

enddo ! loop over iv

end subroutine transform_rhs_for_axis_element

! Transform the test functions of an element RHS stored as (dof, vertex, harmonic), as in the particle
! projections. Same transform as transform_rhs_for_axis_element, applied one harmonic at a time.
! This is needed because transform_rhs_for_axis_element expects the RHS to be stored as a flat vector,
! as in the fluid.
subroutine transform_projection_rhs_for_axis_element(node_list, element, ielm, RHS, n_tor_local)
use data_structure
implicit none
type(type_node_list), intent(in)    :: node_list
type(type_element),   intent(in)    :: element
integer,              intent(in)    :: ielm, n_tor_local
real*8,               intent(inout) :: RHS(n_degrees, n_vertex_max, n_tor_local)

type(type_node) :: nodes(n_vertex_max)
integer         :: iv, in

if (.not. any(node_list%node(element%vertex(1:n_vertex_max))%axis_node)) return

do iv = 1, n_vertex_max
  call make_deep_copy_node(node_list%node(element%vertex(iv)), nodes(iv))
enddo

do in = 1, n_tor_local
  call transform_rhs_for_axis_element(nodes, ielm, RHS(:,:,in), (/ 1 /), 1, (/ 1 /), 1)
enddo

do iv = 1, n_vertex_max
  call dealloc_node(nodes(iv))
enddo

end subroutine transform_projection_rhs_for_axis_element

! Transformation matrices P (old dofs = P new dofs) of the four vertices of element ielm.
! Only the axis vertices 1 and 4 are transformed, the others get the identity.
subroutine axis_element_transform_matrices(nodes, ielm, Pmat)
use data_structure
use phys_module
implicit none
type(type_node), intent(in)  :: nodes(n_vertex_max)
integer,         intent(in)  :: ielm
real*8,          intent(out) :: Pmat(1:4, 1:4, n_vertex_max)

integer :: axis_vertex1, axis_vertex4, iv, io
real*8  :: theta

axis_vertex1 = 1 ; axis_vertex4 = 4

do iv = 1, n_vertex_max
  theta = 0.d0
#if STELLARATOR_MODEL
  if (iv .le. 2) then
    theta = 2.d0*pi*float(mod(ielm-1,n_tht))/float(n_tht)
  else
    theta = 2.d0*pi*float(mod(ielm,n_tht))/float(n_tht)
  end if
#endif
  if (iv==axis_vertex1 .or. iv==axis_vertex4) then
    call axis_transform_matrix(nodes(iv), theta, Pmat(:,:,iv))
  else
    Pmat(:,:,iv) = 0.d0
    do io = 1, n_order+1
      Pmat(io, io, iv) = 1.d0
    enddo
  end if
enddo

end subroutine axis_element_transform_matrices

! Transformation matrix P (old dofs = P new dofs) of one node. It is the identity unless the node is
! an axis node (axis_dof 2 or 3). theta is the poloidal angle of the node, used only by the stellarator grid.
subroutine axis_transform_matrix(node, theta, Pmat)
use data_structure
use phys_module
implicit none
type(type_node), intent(in)  :: node
real*8,          intent(in)  :: theta
real*8,          intent(out) :: Pmat(1:4, 1:4)

integer :: dof1, dof2, dof3, dof4, jo

dof1 = 1 ; dof2 = 2 ; dof3 = 3 ; dof4 = 4

! initialize the transformation matrix to identity
Pmat = 0.d0
do jo = 1, n_order+1
  Pmat(jo, jo) = 1.d0
enddo

if (node%axis_dof == dof2) then
#if STELLARATOR_MODEL
  Pmat(dof3,dof2) = cos(theta)/(3.d0*float(n_flux-1)/bloating_factor); Pmat(dof3,dof3) = sin(theta)/(3.d0*float(n_flux-1)/bloating_factor)
  Pmat(dof4,dof2) = -2.d0*pi*sin(theta)/(9.d0*float(n_tht*(n_flux-1))/bloating_factor); Pmat(dof4,dof3) = 2.d0*pi*cos(theta)/(9.d0*float(n_tht*(n_flux-1))/bloating_factor)
#else
  Pmat(dof3,dof2) = node%x(1,dof3,1) ; Pmat(dof3,dof3) = node%x(1,dof3,2)
  Pmat(dof4,dof2) = node%x(1,dof4,1) ; Pmat(dof4,dof3) = node%x(1,dof4,2)
#endif
  Pmat(dof2,dof2) = 0.d0
  Pmat(dof4,dof4) = 0.d0
else if (node%axis_dof == dof3) then
#if STELLARATOR_MODEL
  Pmat(dof2,dof2) = cos(theta)/(3.d0*float(n_flux-1)/bloating_factor); Pmat(dof2,dof3) = sin(theta)/(3.d0*float(n_flux-1)/bloating_factor)
  Pmat(dof4,dof2) = -2.d0*pi*sin(theta)/(9.d0*float(n_tht*(n_flux-1))/bloating_factor); Pmat(dof4,dof3) = 2.d0*pi*cos(theta)/(9.d0*float(n_tht*(n_flux-1))/bloating_factor)
#else
  Pmat(dof2,dof2) = node%x(1,dof2,1) ; Pmat(dof2,dof3) = node%x(1,dof2,2)
  Pmat(dof4,dof2) = node%x(1,dof4,1) ; Pmat(dof4,dof3) = node%x(1,dof4,2)
#endif
  Pmat(dof3,dof3) = 0.d0
  Pmat(dof4,dof4) = 0.d0
end if

end subroutine axis_transform_matrix

! Since the axis treatment solves for new degrees of freedom, we need to
! transforms degrees of freedom to old ones.
subroutine new_to_old_dofs_on_the_axis(node_list, vg, new_dofs, old_dofs)
use data_structure
use phys_module
implicit none
type(type_node_list),    intent(inout) :: node_list
integer,   intent(in) :: vg
real*8,    intent(in) :: new_dofs(4)
real*8,    intent(out):: old_dofs(4)
real*8  :: Pmat(4,4), theta

theta = 2.d0*pi*float(mod(vg-1,n_tht))/float(n_tht)
call axis_transform_matrix(node_list%node(vg), theta, Pmat)

old_dofs = matmul(Pmat, new_dofs)

end subroutine new_to_old_dofs_on_the_axis

!> Add condition for the axis directly in the matrix.
!! - This aims to apply C0 continuity on the grid axis and is used only when treat_aixs=.t.
!! - The fourth dof is penalized (d^2/dRdZ) when the treat_axis option is used.
subroutine penalize_dof_on_axis(node_list, dof, element_list, local_elms, n_local_elms, index_min, index_max, a_mat)

  use mod_assembly, only : boundary_conditions_add_one_entry, boundary_conditions_add_RHS
  use data_structure
  use mod_locate_irn_jcn

  implicit none

  ! Subroutine parameters
  integer,                   intent(in)    :: local_elms(*)         !< List of local elements
  integer,                   intent(in)    :: n_local_elms          !< Number of local elements
  integer,                   intent(in)    :: index_min, index_max  !< Min/max index of local elements
  type (type_node_list),     intent(in)    :: node_list             !< List of nodes
  integer,                   intent(in)    :: dof                   !< which dof to penalize  
  type (type_element_list),  intent(in)    :: element_list          !< List of all elements
  type(type_SP_MATRIX)                     :: a_mat
  ! Internal parameters
  real*8  :: zbig
  integer :: i, in, iv, inode, k
  integer :: index_large_i, index_node, index_node2, ielm
  integer :: ilarge2, n_tor_local
  integer(kind=int_all):: ijA_position,ijA_position2

  n_tor_local = (a_mat%i_tor_max - a_mat%i_tor_min + 1)

  zbig = 1.d12
  do i=1, n_local_elms

    ielm = local_elms(i)

    do iv=1, n_vertex_max

      inode = element_list%element(ielm)%vertex(iv)

      if (node_list%node(inode)%axis_node) then

        do in=a_mat%i_tor_min, a_mat%i_tor_max
          do k=1, n_var

            index_node = node_list%node(inode)%index(dof)
            if ((index_node .ge. index_min) .and. (index_node .le. index_max)) then
              call locate_irn_jcn(index_node,index_node,index_min,index_max,ijA_position,a_mat)
              index_large_i = n_tor_local * n_var * (index_node - 1)
              ilarge2 = ijA_position - 1 + ((k-1)*n_tor_local + in-a_mat%i_tor_min) * n_var*n_tor_local &
                + (k-1)*n_tor_local + in - a_mat%i_tor_min + 1
              a_mat%irn(ilarge2) =  n_tor_local * n_var * (index_node-1) + (k-1)*n_tor_local + in - a_mat%i_tor_min + 1
              a_mat%jcn(ilarge2) =  n_tor_local * n_var * (index_node-1) + (k-1)*n_tor_local + in - a_mat%i_tor_min + 1
              a_mat%val(ilarge2)   = zbig
            end if

          enddo
        enddo

      endif

    enddo  ! n_vertex
  enddo ! n_elements

  return

end subroutine penalize_dof_on_axis

end module mod_axis_treatment

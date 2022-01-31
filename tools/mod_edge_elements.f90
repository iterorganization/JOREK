!> Module with edge elements which can be written to vtk
!>
!> We store these as linear 2D elements, with one direction in the poloidal plane
!> and one in the toroidal direction. Store M disconnected sets of elements.
!> At domain corners we do not require continuity.
!> Commonly there are 4 domain corners
!>
!> The edge element and node numbering is generated as follows:
!>
!> Get M lists of N_i JOREK elements to include.
!> They are identified by the index in the JOREK element list and the side of the element.
!> This is expected in a list of type(edge_domain).
!> Store the linear elements in edge_distribution as an array of integers i_elm of size N = sum_i N_i
!> and an array of integers i_side of size N.
!> The offset of every N_i in that array is stored in an array of dimension(M+1) with the last
!> element having value N (and the first 1, commonly).
!>
!> The numbering order of the subelements and nodes is as follows:
!> Sides first, then elements, then subdivisions, then toroidal subdivisions
!>
!> Pseudocode:
!> do i_list=1,M
!>   do i_edge_element=offset(i_list),offset(i_list+1)
!>     do i_sub=1,n_sub
!>       i_elm += n_plane * n_sub_toroidal
!>       i_node += n_plane * n_sub_toroidal
!>     end do
!>     if (i_edge_element .eq. offset(i_list+1)) i_node += n_plane * n_sub_toroidal
!>   end do
!> end do
!>
!>
!> Example for the case of 2 lists, 3 planes, 2 edge elements per list, n_sub=1 and n_sub_toroidal=2
!> The dots indicate wrap-around the toroidal direction phi. l is the length
!> in the poloidal direction along the edge elements.
!> Total number of elements = n_side * n_plane * n_sub * n_sub_toroidal * n_elements    = 2*3*2*2 = nel = 24
!> Total number of nodes = n_side * n_plane * n_sub * n_sub_toroidal * (n_elements + 1) = 2*3*2*3 = nnos = 36
!>
!>     13---14---15---16---17---18...13
!> el 2 |  7 |  8 |  9 | 10 | 11 | 12 |    side=1
!>      7----8----9---10---11---12....7
!> el 1 |  1 |  2 |  3 |  4 |  5 |  6 |    side=1
!>      1----2----3----4----5----6....1
!>     31---32---33---34---35---36...31
!> el 4 | 19 | 20 | 21 | 22 | 23 | 24 |    side=2
!>     25---26---27---28---29---30...25               ^ l
!> el 3 | 13 | 14 | 15 | 16 | 17 | 18 |    side=2     |
!>     19---20---21---22---23---24...19               +---> phi
!>
!>      |                             |
!>   phi=0                         phi=TWOPI/n_period
!> 
!>
!> The nsub subelements are chosen at uniform positions in the element-local coordinate s or t (depending on side)
!> We assume R, Z to be linear between elements, i.e. nsub should be relatively large compared to 1.
!> R, Z, phi per node is stored in the xyz array in vtk_grid.
!> The connectivity, i.e. which nodes belong in an element is stored in ien(1:4,1:nel) in vtk_grid.
!> The values are stored in edge_distribution_vtk%scalars(1:nnos,1:n_scalars)
!>
!> For the VTK output we repeat the elements and nodes in the toroidal direction if n_period .ne. 1.
module mod_edge_elements
use data_structure
use mod_particle_sim
use mod_fields
use mod_edge_domain

implicit none
private
public :: elm_coords
public :: edge_elements
public :: sample_edge_elements
public :: find_edge_element

!> 1D Finite elements on the plasma boundary on a connected patch
!> goes around in the toroidal direction (with n_period)
type :: edge_elements_patch
  integer, public :: nsub = 6 !< Number of subdivision per edge
  integer, public :: nsub_toroidal = 4 !< Number of subdivisions per plane (i.e. n_plane * nsub_toroidal points in toroidal dir)

  type(type_edge_domain) :: edge !< JOREK edge domain description

  real*4, allocatable  :: xyz(:,:) !< positions of points (3,nnos) | xyz is written in RZPhi coordinates
  real*8, allocatable  :: st(:,:) !< (2,nnos) | local coordinates of points
  integer,allocatable  :: i_elm_jorek_edge(:) !< Index of jorek element for each edge element (/node on/before an element)
  
  integer, allocatable :: ien(:,:) !< connectivity in vtk file (4,nel)
  real*4, dimension(:,:), allocatable :: scalars !< (nnos,n_scalars)
  character(len=36), dimension(:), allocatable :: scalar_names !< (n_scalars)
end type edge_elements_patch

type :: edge_elements
  type(edge_elements_patch), allocatable, dimension(:) :: patch !< n_domains different patches of edge elements
contains
  procedure :: prepare !< this%prepare(node_list, element_list, edge_domains, nsub, nsub_toroidal)
  procedure :: write !< this%write(filename)
  procedure :: sample => sample_edge_elements !< this%sample(i_scalar, n_samples, u, integral, xyz, st, ielm)
end type edge_elements

contains

!> Prepare the elements along the plasma edge. Writes the xyz positions of nodes
!> and the connectivity matrix to xyz and ien in this, as well as storage of the
!> elements mentioned in edge_domains. Combined into a single finite-element
!> surface, consisting of size(edge_domains,1) disconnected components.
subroutine prepare(this, node_list, element_list, edge_domains, nsub, nsub_toroidal)
  use mpi
  use mod_event
  use constants, only: TWOPI
  use mod_interp
  use mod_edge_domain 
  !$ use omp_lib
  class(edge_elements), intent(inout) :: this
  type(type_node_list), intent(in)    :: node_list
  type(type_element_list), intent(in) :: element_list
  type(type_edge_domain), intent(in), dimension(:) :: edge_domains
  integer, intent(in)                 :: nsub, nsub_toroidal
  integer :: n_jorek_elm, n_edge_elm, n_edge_node, n_domain
  integer :: i, j, k, i_domain, i_sub, i_sub_tor
  real*8 :: st(2), R, Z

  n_domain = size(edge_domains,1)
  allocate(this%patch(n_domain))
  
  do i_domain=1,n_domain
    this%patch(i_domain)%edge = edge_domains(i_domain)
    this%patch(i_domain)%nsub = nsub
    this%patch(i_domain)%nsub_toroidal = nsub_toroidal
    n_jorek_elm = size(edge_domains(i_domain)%i_elm,1)
    n_edge_elm = n_jorek_elm * nsub * n_plane * nsub_toroidal
    n_edge_node = (n_jorek_elm * nsub +1 )* n_plane * nsub_toroidal

    allocate(this%patch(i_domain)%xyz(3,n_edge_node), this%patch(i_domain)%ien(4,n_edge_elm))
    allocate(this%patch(i_domain)%st(2,n_edge_node))
    allocate(this%patch(i_domain)%i_elm_jorek_edge(n_edge_node))
    
    i = 1

    ! Node and element definition. We do not define the last row of nodes here, only the ones that map directly
    ! to an element.
    do j=0,n_jorek_elm-1
      do i_sub=1,nsub
        do i_sub_tor=1,n_plane*nsub_toroidal
          k = (i_sub-1)*(n_plane*nsub_toroidal) + j*nsub*n_plane*nsub_toroidal + i_sub_tor
          !< See [[coord_in_neighbour]] in mod_neighbours.f90 for explanation of sides
          st = elm_coords(edge_domains(i_domain)%i_side(i+j), real(i_sub-1,8)/real(nsub,8))
          this%patch(i_domain)%st(:,k) = st !< save st to be able to easily link the edge element node number to the JOREK coordinate
          this%patch(i_domain)%i_elm_jorek_edge(k) = edge_domains(i_domain)%i_elm(i+j)
          call interp_RZ(node_list, element_list, edge_domains(i_domain)%i_elm(i+j), st(1), st(2), R, Z)
            
          this%patch(i_domain)%xyz(:,k) = real([R, Z, TWOPI*real(i_sub_tor-1)/real(n_plane*nsub_toroidal*n_period)],4)
          ! ien contains the coordinates of the edge points in counter-clockwise order.
          ! Elements point outwards from the domain (since particle fluxes are positive)
          ! orientation of JOREK elements is counterclockwise in poloidal plane
          ! reverse 1-4 below to flip orientation

          ! If direction = -1 we flip the orientation of the element
          
          if (edge_domains(i_domain)%direction(i+j) .eq. 1) then
            this%patch(i_domain)%ien(1:4:1, k) = -1 + k + [&
                0, 1, n_plane*nsub_toroidal + 1, n_plane*nsub_toroidal]
            ! If it is the last element in the row we need to subtract n_plane*this%patch(i_domain)%nsub_toroidal from the rightmost 2 node indices
          else if (edge_domains(i_domain)%direction(i+j) .eq. -1) then
            this%patch(i_domain)%ien(4:1:-1, k) = -1 + k + [&
                0, 1, n_plane*nsub_toroidal + 1, n_plane*nsub_toroidal]          
          end if
          
          if (i_sub_tor .eq. n_plane*nsub_toroidal) then
            this%patch(i_domain)%ien(2:3, k) = this%patch(i_domain)%ien(2:3, k) - n_plane*nsub_toroidal
          end if
        end do
      end do
    end do
    ! Define now the boundary nodes on the 'top' of this domain
    do i_sub_tor=1,n_plane*nsub_toroidal
      k = nsub*(n_plane*nsub_toroidal) + (n_jorek_elm-1)*nsub*n_plane*nsub_toroidal + i_sub_tor
      !< See [[coord_in_neighbour]] in mod_neighbours.f90 for explanation of sides
      st = elm_coords(edge_domains(i_domain)%i_side(i+n_jorek_elm-1), real(i_sub-1,8)/real(nsub,8))
      this%patch(i_domain)%st(:,k) = st
      this%patch(i_domain)%i_elm_jorek_edge(k) = edge_domains(i_domain)%i_elm(i+n_jorek_elm-1)
      call interp_RZ(node_list, element_list, edge_domains(i_domain)%i_elm(i+n_jorek_elm-1), st(1), st(2), R, Z)
      this%patch(i_domain)%xyz(:,k) = real([R, Z, TWOPI*real(i_sub_tor-1)/real(n_plane*nsub_toroidal*n_period)],4)
    end do
  end do
end subroutine prepare

!> Calculate the coordinates a distance x along side i_side
pure function elm_coords(i_side, x) result(st)
  integer, intent(in) :: i_side
  real*8, intent(in) :: x
  real*8 :: st(2)
  select case (i_side)
  case (1)
    st = [x,0.d0]
  case (2)
    st = [1.d0,x]
  case (3)
    st = [1.d0-x,1.d0]
  case (4)
    st = [0.d0,1.d0-x]
  end select
end function elm_coords


!> Find the element containing a point. This is done by locating the side of the position
!> in the element provided and comparing s, t with nsub.
!>
!> This can lead to a failure mode if a particle is lost very close to the
!> of an element, as it might then be recognized to be on a different side,
!> and not included.
function find_edge_element(this, i_elm, s, t, phi) result(i_elm_edge)
  use mod_parameters, only: n_plane, n_period
  use constants, only: TWOPI
  class(edge_elements_patch), intent(in) :: this
  integer, intent(in)                    :: i_elm
  real*8, intent(in)                     :: s, t, phi
  integer                                :: i_elm_edge

  integer :: side, i_edge, n_edge, n_edge_elm, i
  real*8  :: s_edge
  ! See [[mod_boundary]] for an explanation
  if (s .gt. t) then
    if (1.d0 - s .gt. t) then
      side = 1
      s_edge = s
    else
      side = 2
      s_edge = t
    end if
  else
    if (1.d0 - s .gt. t) then
      side = 4
      s_edge = 1.d0-t
    else
      side = 3
      s_edge = 1.d0-s
    end if
  end if

  n_edge = size(this%edge%i_elm,1)
  n_edge_elm = n_edge * this%nsub * n_plane * this%nsub_toroidal
  ! Find the first place where i_elm matches and i_side matches
  i_edge = minloc([(i,i=1,n_edge)], mask=(this%edge%i_elm .eq. i_elm .and. this%edge%i_side .eq. side),dim=1)
  if (i_edge .eq. 0) then
    i_elm_edge = -999 ! not found
    return
  end if
  s_edge = s_edge*this%edge%direction(i_edge) ! flip direction if the element is cw
  ! i_edge is the number of the edge domain index. This is related to real element numbers by nsub, nsub_toroidal and n_plane
  if (s_edge .ge. 1.d0) s_edge = 0.9999999999d0 ! bad stuff if exactly 1
  i_elm_edge = 1+ (i_edge-1) * this%nsub * n_plane*this%nsub_toroidal + &
               floor(s_edge*real(this%nsub,8)) * n_plane*this%nsub_toroidal + & 
               floor(mod(phi, TWOPI/n_period) /TWOPI * real(n_period*n_plane*this%nsub_toroidal,8)) ! Modulus needed for unbounded phi
end function find_edge_element


!> Save an already-projected set to a vtk file with current parameters
!> loop over edge domains and write each of them
!>
!> First we need to get the total number of points and prepare the scalars and the total connectivity matrix
!>
!> Assumptions: nsub_toroidal and n_plane are the same for each patch
!>
!> Writes only what is present on MPI ID 0! You need to reduce this yourself to a single set.
subroutine write(this, filename)
  use mpi
  use mod_event
  use mod_coordinate_transforms
  use mod_vtk
  class(edge_elements), intent(in) :: this
  character(len=*), intent(in)     :: filename

  integer, parameter :: etype = 9 ! for vtk_quad
  integer :: n_scalars, nnos, nel, n_rows, n_cols
  integer :: n_scalars_i, nnos_i, nel_i
  real*4, allocatable :: vectors(:,:,:)
  character*36, allocatable :: vector_names(:)
  real*4, allocatable :: xyz(:,:)
  integer, allocatable :: ien(:,:)
  real*4, allocatable :: scalars(:,:)

  integer :: i, j, n, ierr, my_id

  call mpi_comm_rank(MPI_COMM_WORLD, my_id, ierr)

  nnos = 0
  nel = 0
  n_scalars = 0
  do i=1,size(this%patch,1)
    if (.not. allocated(this%patch(i)%scalars)) then
      write(*,*) "No scalars defined for edge element writing, skipping"
      cycle
    end if
    nnos = nnos + size(this%patch(i)%scalars,1)
    nel  = nel  + size(this%patch(i)%ien,2)
    n_scalars = max(n_scalars, size(this%patch(i)%scalars,2))
  end do

  allocate(scalars(nnos, n_scalars))
  scalars = 0
  nnos_i = 0
  do i=1,size(this%patch,1)
    ! create a single array of scalars
    n_scalars_i = size(this%patch(i)%scalars,2)
    scalars(nnos_i+1:size(this%patch(i)%scalars,1), 1:n_scalars_i) = this%patch(i)%scalars
    nnos_i = nnos_i + size(this%patch(i)%scalars,1)
  end do
    
  if (my_id .eq. 0) then
    ! write only on the host
    n_scalars = size(scalars,2)
    allocate(vectors(0,3,0)) !----------------------------- array of length 0 (0,3,0)
    allocate(vector_names(0))
    ! Hardcode scalar names of first patch...
    if (.not. allocated(this%patch(1)%scalar_names)) then
      write(*,*) 'Scalar names not allocated for mod_edge_elements/write, skipping'
      return
    end if

    ! VTK expects cartesian coordinates
    allocate(xyz(3,nnos))
    nnos = 0
    nnos_i = 0
    do i=1,size(this%patch,1)
      nnos_i = size(this%patch(i)%xyz,2)
      do j=nnos+1,nnos+nnos_i
        xyz(:,j) = cylindrical_to_cartesian(this%patch(i)%xyz(:,j))
      end do
      nnos = nnos+nnos_i
    end do
    
    
    ! If the constructed grid is fully toroidal connected the grid will something like the grid above
    ! However, if the grid is not fully toroidal, we don't want to connect the end points,
    ! So the last element of every row must be removed if n_period .ne. 1    
    
    
    ! remove connecting elements between endpoints if n_period != 1
    ! Every n_plane*n_sub_toroidal-th element disappears (those for which i_elm % n_plane*n_sub_toroidal == 0)
    ! we need to update the ien matrix with nnos_i to shift correctly
    if (n_period .ne. 1) then
      n_rows = nel/(n_plane*this%patch(1)%nsub_toroidal)
      n_cols = nel/n_rows
      allocate(ien(4,(n_rows)*(n_cols-1)))
      nel_i = 0
      nnos_i = 0
      do i=1,size(this%patch,1)
        n = size(this%patch(i)%ien,2)
        ! assume nsub_toroidal equal between patches
        !>     13---14---15---16---17---18...13
        !> el 2 |  7 |  8 |  9 | 10 | 11 | 12 |    side=1
        !>      7----8----9---10---11---12....7
        !> el 1 |  1 |  2 |  3 |  4 |  5 |  6 |    side=1
        !>      1----2----3----4----5----6....1
        ! select in this case all ien except 6 and 12
        ! so we select nel_i + 1  * (j-1) * nt : nel_i + (j)*nt - 1
        do j=1,n/n_cols ! row index
          ien(:,nel_i+1+(j-1)*(n_cols-1):nel_i + j*(n_cols-1)) = &
              this%patch(i)%ien(:,1+(j-1)*n_cols:j*n_cols-1) &
              + nnos_i ! shift upwards with number of nodes before us
        end do
        nel_i = nel_i + n*(n_rows-1)/(n_rows)
        nnos_i = nnos_i + size(this%patch(i)%xyz,2)
      end do
    else
      allocate(ien(4,nel))
      nel_i = 0
      nnos_i = 0
      do i=1,size(this%patch,1)
        n = size(this%patch(i)%ien,2)
        ien(:,nel_i+1:nel_i+n) = this%patch(i)%ien + nnos_i
        nel_i = nel_i + n
        nnos_i = nnos_i + size(this%patch(i)%xyz,2)
      end do
    endif

    call write_vtk(filename,xyz,&
      ien, etype,&
      this%patch(1)%scalar_names,scalars,&
      vector_names,vectors)
      
  end if
end subroutine write

!> Sample from the edge elements.
!> We need to first integrate in the toroidal direction for every poloidal point.
!> Then we integrate in the poloidal direction over every patch, to determine the
!> distribution between patches (and normalize all CDFs).
!> Finally we can sample with the following algorithm
!>
!> do i=1,n_samples
!>   select patch with u(1,i) in discrete array of CDFs
!>   sample linear with u(1,i) from CDF of selected patch
!>   sample linear with u(2,i) from CDF at selected point
!> end do
!>
!> The probability density on the surface is given by Gamma/integral
!> Calculate the integral as
!> \int_0^{2\pi}\int_l \Gamma(r,phi) r dr dphi
!> First we eliminate the phi direction leading to
!> \int_l 2\pi \tilde\Gamma(r) r dr
!> This we can then sample from.
!>
!> Note that the sampled positions are not exactly corresponding to the sampled
!> element-local coordinates!
subroutine sample_edge_elements(this, i_scalar, n_samples, u, integral, xyz, st, ielm)
  use constants, only: TWOPI
  use mod_sampling, only: sample_piecewise_linear
  class(edge_elements), intent(in) :: this
  integer, intent(in) :: i_scalar
  integer, intent(in) :: n_samples
  real*8, intent(in)  :: u(2,n_samples)
  real*8, intent(out) :: integral
  real*8, intent(out) :: xyz(3,n_samples)
  real*8, intent(out) :: st(2,n_samples)
  integer, intent(out) :: ielm(n_samples)


  type :: arr1
    real*8, allocatable :: a(:)
  end type arr1

  real*8, allocatable :: pdf_patch(:) !< per patch
  real*8, allocatable :: cdf_patch(:) !< per patch
  type(arr1), allocatable :: pdf_pol(:) !< poloidal direction, per patch

  real*8 :: dphi
  integer :: i, j, n, nphi, i_patch, n_patch, i_min, i_max
  real*8 :: u_tmp(2)
  real*8, allocatable :: dl(:), r(:), p_phi(:)
  real*8 :: i_r, f_min, f_max

  if (.not. allocated(this%patch)) then
    write(*,*) 'No patches allocated for sample_edge_elements, exiting'
    return
  end if
  n_patch = size(this%patch,1)
  allocate(pdf_patch(n_patch))
  allocate(cdf_patch(0:n_patch))
  cdf_patch(0) = 0.d0
  allocate(pdf_pol(n_patch))

  ! Prepare PDFs
  nphi = n_plane*this%patch(1)%nsub_toroidal
  dphi = TWOPI / (n_period*nphi)
  do i_patch=1,n_patch
    n = size(this%patch(i_patch)%st,2)/nphi
    allocate(pdf_pol(i_patch)%a(n))
    allocate(r(n))

    !$omp parallel do default(shared) private(j)
    do j=1,n
      ! The phi-integral is given by the sum over the elements (they are equal size)
      ! pdf_pol(i_patch)%a(j) = sum(this%patch(i_patch)%scalars(i_scalar,(j-1)*nphi+1:j*nphi))*dphi
      pdf_pol(i_patch)%a(j) = sum(this%patch(i_patch)%scalars((j-1)*nphi+1:j*nphi,i_scalar))*dphi
    end do
    !$omp end parallel do

    dl = norm2(this%patch(i_patch)%xyz(1:2,     1:nphi*(n-1):nphi) - &
               this%patch(i_patch)%xyz(1:2,nphi+1:nphi*n    :nphi), dim=1)
    r = this%patch(i_patch)%xyz(1,1:nphi*n:nphi)

    ! Adjust the pdf with the jacobian and length to turn it into an integral over 1-N
    pdf_pol(i_patch)%a = pdf_pol(i_patch)%a(1:n) * r(1:n) * [dl(1), ((dl(j)+dl(j+1))/2, j=1,n-2), dl(n-1)]


    ! Integrate this pdf over the patch
    ! Two contributions: rising basis function and descending
    ! integral for each element is given by 0.5*dl*a for both of these functions
    !pdf_patch(i_patch) = 0.5d0*(sum(pdf_pol(i_patch)%a(1:n-1)*dl*r(1:n-1)) + &
                                !sum(pdf_pol(i_patch)%a(2:n)  *dl*r(2:n)))
    pdf_patch(i_patch) = sum(pdf_pol(i_patch)%a(2:n-1)) + 0.5d0*(pdf_pol(i_patch)%a(1) + pdf_pol(i_patch)%a(n))
    cdf_patch(i_patch) = cdf_patch(i_patch-1) + pdf_patch(i_patch)
    deallocate(r)
  end do

  integral = cdf_patch(n_patch)
  if (integral .le. 1d-30) then
    ! create only lost particles
    ielm = 0
    xyz = 0
    st = 0
    return
  end if

  !$omp parallel do default(none) shared(n_samples, u, cdf_patch, integral, pdf_pol, xyz, st, ielm, nphi, dphi, this, i_scalar, &
  !$omp pdf_patch) &
  !$omp private(i, u_tmp, i_patch, i_r, i_min, i_max, f_min, f_max, p_phi)
  do i=1,n_samples
    u_tmp = u(:,i)
    i_patch = minloc(cdf_patch, mask=(u_tmp(1)*integral .le. cdf_patch),dim=1)-1 ! -1 is due to cdf_patch starting at 0
    
    ! We use a trick here. Instead of putting the `real` distance l, we use the node numbers
    ! instead. This way we can very easily interpolate later.
    ! rescale our random number to this range
    i_r = sample_piecewise_linear(size(pdf_pol(i_patch)%a,1), &
                                  [(real(j,8),j=1,size(pdf_pol(i_patch)%a,1))], &
                                  pdf_pol(i_patch)%a, (u_tmp(1)*integral-cdf_patch(i_patch-1))/pdf_patch(i_patch))

    ! floor(i_r) is the lower node, ceil(i_r) is the upper, and if they are the same we're exactly on a node
    i_min = floor(i_r)
    i_max = ceiling(i_r)
    f_min = 1-mod(i_r, 1.d0)
    f_max = 1-f_min

    ! Fill in the positions
    ielm(i) = this%patch(i_patch)%i_elm_jorek_edge(i_min*nphi)
    xyz(1:2,i) = this%patch(i_patch)%xyz(1:2,i_min*nphi)*f_min + &
                 this%patch(i_patch)%xyz(1:2,i_max*nphi)*f_max
    ! We need to be careful here, since st is defined in the element with the
    ! same number as the node! (but then do some nodes not have it defined?)
    ! So if the i_elm_jorek_edge is different between the two, then
    ! we need to flip the latter one, to get it from 0 to 1 or vice versa.
    ! we don't know if this is the s or t coordinate though...
    if (this%patch(i_patch)%i_elm_jorek_edge(i_min*nphi) .ne. &
        this%patch(i_patch)%i_elm_jorek_edge(i_max*nphi)) then
      associate(s0 => this%patch(i_patch)%st(:,i_min*nphi), &
                s1 => this%patch(i_patch)%st(:,i_max*nphi))
      if (s0(1) .ne. s1(1)) then
        st(:,i) = s0*f_min + &
                  [1.0-s1(1), s1(2)]*f_max
      elseif (s0(2) .ne. s1(2)) then
        st(:,i) = s0*f_min + &
                  [s1(1), 1.0-s1(2)]*f_max
      else
        ! This only happens when nsub=1, i.e. probably both
        ! s0(1) and s1(1) are 0, or both s0(2) and s1(2) are zero.
        ! since the other coordinate can be zero or one depending on the side
        ! we need to check the boundary information of the JOREK element to
        ! find out which side it was supposed to be
        ! we don't have this info here, so just print a warning and let people
        ! use nsub > 1 (which is really necessary for convergence anyway)
        write(*,*) "ERROR: cannot determine side of element in sample_edge_elements. Use nsub > 1" ! or save i_elm_jorek_side in this%patch
      end if
      end associate
    else
      st(:,i) = this%patch(i_patch)%st(:,i_min*nphi)*f_min + &
                this%patch(i_patch)%st(:,i_max*nphi)*f_max
    end if

    ! Sample the toroidal position phi
    allocate(p_phi(nphi+1))
    p_phi(1:nphi) = this%patch(i_patch)%scalars((i_min-1)*nphi+1:i_min*nphi,i_scalar)*f_min + &
                    this%patch(i_patch)%scalars((i_max-1)*nphi+1:i_max*nphi, i_scalar)*f_max
    p_phi(nphi+1) = p_phi(1) ! make it circular
    xyz(3,i) = sample_piecewise_linear(nphi+1, [(real(j-1,8)*dphi*n_period,j=1,nphi+1)], p_phi, u_tmp(2))
    deallocate(p_phi)

  end do
  !$omp end parallel do

end subroutine sample_edge_elements
end module mod_edge_elements

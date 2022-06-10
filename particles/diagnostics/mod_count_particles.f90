!> Module containing routines to project particles onto the JOREK finite elements
module mod_count_particles
use mod_io_actions
use data_structure
use mod_particle_sim
implicit none
private
public count_particles, write_particle_counts_to_vtk
public count_to_vtk

!> Action to project all particle distributions and save them to vtk
type, extends(io_action) :: count_to_vtk
  type(type_node_list), allocatable    :: node_list !< node lists to save particle projections in
  type(type_element_list), allocatable :: element_list
  integer :: nsub !< Number of subdivisions per element (1 = 1 subelement per element, 2 = 4 etc)
  real*4,allocatable, private  :: xyz (:,:) !< positions of points in vtk file
  integer,allocatable, private :: ien (:,:) !< connectivity in vtk file
contains
  procedure :: do => count_and_save_to_vtk
end type count_to_vtk
interface count_to_vtk
  module procedure new_count_to_vtk
end interface count_to_vtk

contains
!> Constructor for project_to_vtk
!> Be sure to use keyword arguments when initializing, to avoid confusion
function new_count_to_vtk(node_list, element_list, nsub, filename, basename, decimal_digits, fractional_digits) result(new)
  use mpi
  type(count_to_vtk) :: new
  type(type_node_list), intent(in)       :: node_list
  type(type_element_list), intent(in)    :: element_list
  integer, intent(in), optional          :: nsub !< number of subdivisions of the finite elements
  character(len=*), intent(in), optional :: filename
  character(len=*), intent(in), optional :: basename
  integer, intent(in), optional          :: decimal_digits
  integer, intent(in), optional          :: fractional_digits
  integer :: my_id, ierr
  allocate(new%node_list,    source=node_list)
  allocate(new%element_list, source=element_list)
  new%nsub = 4
  if (present(nsub)) new%nsub = nsub
  if (present(filename)) new%filename = filename
  if (present(basename)) then
    new%basename = basename
  else
    new%basename = 'count'
  end if
  if (present(decimal_digits)) new%decimal_digits = decimal_digits
  if (present(fractional_digits)) new%fractional_digits = fractional_digits
  new%extension ='.vtk'
  new%name = "CountToVtk"
  new%log = .true.

  ! Precalculate the node positions in the vtk file and the connectivity
  call MPI_COMM_RANK(MPI_COMM_WORLD, my_id, ierr)
  if (my_id .eq. 0) call prepare_write_particle_counts_to_vtk(new%node_list, new%element_list,new%nsub,new%xyz,new%ien)
end function new_count_to_vtk

!> Action for projecting all particles and writing output to vtk
subroutine count_and_save_to_vtk(this, sim, ev)
  use mpi
  use mod_event
  !$ use omp_lib
  class(count_to_vtk), intent(inout) :: this
  type(particle_sim), intent(inout)    :: sim
  type(event), intent(inout), optional :: ev
  integer :: i, my_id, ierr
  character(len=120) :: filename
  real*4, allocatable, dimension(:,:,:) :: counts !< element_id, i_tor, group_id (real because multiplied by toroidal basis function)

  ! Safety checks
  if (.not. allocated(sim%groups)) return

  allocate(counts(this%element_list%n_elements*this%nsub**2,n_tor,size(sim%groups)))
  counts = 0
  do i=1,size(sim%groups)
    if (.not. allocated(sim%groups(i)%particles)) cycle
    call count_particles(this%node_list, this%element_list, counts(:,:,i), &
        this%nsub, sim%groups(i)%particles)
    ! results are saved only on mpi process 0
  end do

  if (len_trim(this%filename) .eq. 0) then
    filename = this%get_filename(sim%time)
  else
    filename = this%filename
  end if

  call MPI_COMM_RANK(MPI_COMM_WORLD, my_id, ierr)
  if (my_id .eq. 0) then
    ! write only on the host
    call write_particle_counts_to_vtk(this%node_list, this%element_list, &
      trim(filename), this%nsub, counts, this%xyz, this%ien)
  end if
  deallocate(counts)
end subroutine count_and_save_to_vtk



!> Perform the actual projection of a set of particles on variable ivar_out in node_list.
subroutine count_particles(node_list, element_list, counts, nsub, particles)
use phys_module
use mpi
use mod_particle_types
!$ use omp_lib
implicit none

type (type_node_list), intent(inout) :: node_list !< A copy of the node list which will be used to save variables
type (type_element_list), intent(in) :: element_list
real*4, dimension(:,:), intent(out) :: counts
integer, intent(in) :: nsub
class (particle_base), intent(in), dimension(:)    :: particles
integer :: i_tor, m, i, i_s, i_t, ierr, my_id
real*4, dimension(:,:), allocatable :: sum_counts
real*4 :: v

call MPI_COMM_RANK(MPI_COMM_WORLD, my_id, ierr)

do i_tor=1, n_tor
  do m=1,size(particles,1)
    if (particles(m)%i_elm .le. 0) cycle
    ! number the elements from low t to high t and then from low s to high s
    !  __________ s=1
    ! |    |    |
    ! | 3  | 4  |
    ! |^^^^|^^^^|
    ! |_1__|_2__| s=0
    ! t=0     t=1
    i_s = ceiling(particles(m)%st(1)*real(nsub))
    i_t = ceiling(particles(m)%st(2)*real(nsub))
    i = (particles(m)%i_elm-1)*nsub**2 + (i_t-1)*nsub + i_s
    if (mode(i_tor) .gt. 1) then ! mode(1) = 0, mode(2) -> cos, mode(3) -> sin
      if (mod(i_tor,2) .eq. 0) then
        v = cos(mode(i_tor)*particles(m)%x(3))
      else
        v = sin(mode(i_tor)*particles(m)%x(3))
      endif
      v = v / PI ! int cos^2(nx) from 0 to 2pi = pi for n > 0
    else
      v = 1.d0 / TWOPI ! int 1 from 0 to 2pi = 2pi
    endif
    counts(i,i_tor) = counts(i,i_tor) + v/particles(m)%x(1) ! normalize with R to get a density
  enddo

enddo ! i_tor
allocate(sum_counts(size(counts,1),size(counts,1)), mold=counts)
call MPI_Reduce(counts,sum_counts,size(counts), &
    MPI_REAL4, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
if (my_id .eq. 0) then
  counts = sum_counts
end if
deallocate(sum_counts)
end subroutine count_particles


!> Calculate the structure of the vtk file without putting in any scalars
subroutine prepare_write_particle_counts_to_vtk(node_list,element_list,nsub,xyz,ien)
use data_structure
use mod_interp
implicit none

!> Input parameters
type(type_node_list), intent(in)    :: node_list
type(type_element_list), intent(in) :: element_list
integer, intent(in)                 :: nsub !< Number of subdivisions of each element
real*4,allocatable, intent(out)     :: xyz (:,:)
integer,allocatable, intent(out)    :: ien (:,:)

integer :: nnos, nnoel, nel, i, j, ielm, inode, k, j_sub, k_sub
real*8 :: s, t, R, R_s, R_t, Z, Z_s, Z_t

nnos = nsub**2*4*element_list%n_elements
allocate(xyz(3,nnos))

nnoel = 4
nel   = nsub*nsub*element_list%n_elements
allocate(ien(nnoel,nel))

inode   = 0
ielm    = 0
xyz     = 0
ien     = 0

! Create points for each element
do i=1,element_list%n_elements
  do j_sub=1,nsub ! number of subelement in s
    do k_sub=1,nsub ! number of subelement in t
      do j=1,2
        s = float(j_sub-2+j)/float(nsub)
        do k=1,2
          t = float(k_sub-2+k)/float(nsub)
          call interp_RZ(node_list,element_list,i,s,t,R,R_s,R_t,Z,Z_s,Z_t)
          inode = (i-1)*nsub*nsub*4 + (k_sub-1)*4*nsub + (j_sub-1)*4 + (j-1)*2+k
          xyz(1:3,inode) = real([R, Z, 0.d0], 4)
        enddo
      enddo
      ielm	  = ielm+1
      inode       = (i-1)*nsub*nsub*4 + (k_sub-1)*4*nsub + (j_sub-1)*4 + 1
      ien(1,ielm) = inode - 1 ! because vtk is 0-based
      ien(2,ielm) = inode + 1
      ien(3,ielm) = inode + 2
      ien(4,ielm) = inode
    enddo
  enddo
enddo  ! n_elements
end subroutine prepare_write_particle_counts_to_vtk


!> Helper subroutine to write a particle distribution, saved in variables 1:n,
!> to `filename`. `nsub` is the number of subdivisions to make per element.
!> Should be called only by process 1
subroutine write_particle_counts_to_vtk(node_list,element_list,filename,nsub,counts,xyz,ien)
use data_structure
use phys_module, only: mode
use mod_vtk
use mod_basisfunctions
use basis_at_gaussian
use gauss, only: Wgauss
!$ use omp_lib
implicit none

!> Input parameters
type(type_node_list), intent(in)    :: node_list
type(type_element_list), intent(in) :: element_list
character*(*), intent(in)           :: filename
integer, intent(in) :: nsub
real*4, intent(in)  :: counts(:,:,:) !< ielm_sub, i_tor, i_group
real*4, intent(in)  :: xyz(:,:)
integer, intent(in) :: ien(:,:)

integer :: nnos, i, j, k, l, m, inode, ivar, ielm, i_elm
real*4, allocatable :: scalars(:,:), vectors(:,:,:)
integer :: n_scalars, n_vectors = 0
character*36, allocatable :: vector_names(:), scalar_names(:)
type(type_element) :: element
type(type_node)    :: nodes(4)
real*8, dimension(n_gauss,n_gauss) :: x_s, x_t, y_s, y_t
integer :: ms, mt, index_ij
real*8 :: wst, xjac
real*4 :: area

integer, parameter :: etype = 9 ! for vtk_quad

n_scalars = n_tor * size(counts,3)
nnos = nsub**2*4*element_list%n_elements
allocate(scalars(nnos,n_scalars),vectors(nnos,3,n_vectors))
allocate(scalar_names(n_scalars),vector_names(n_vectors))
do i=1,size(counts,3)
  write(scalar_names(n_tor*(i-1)+1),'(A,i0.2)') "rho_", i
  do j=1,(n_tor-1)/2
    write(scalar_names(n_tor*(i-1)+j+1),"(A4,i0.2,A4,i0.2)") "rho_", i, "_cos", mode(2*j)
    write(scalar_names(n_tor*(i-1)+j+2),"(A4,i0.2,A4,i0.2)") "rho_", i, "_sin", mode(2*j+1)
  end do
end do

scalars = 0.0
vectors = 0.0

! Create points for each element
!$omp parallel do default(shared) shared(element_list,nsub,node_list,counts,scalars,h_s,h_t) & ! default shared because Wgauss is a constant
!$omp private(i,j,k,l,m,inode,ivar,ielm,wst,x_s,x_t,y_s,y_t,xjac,ms,mt,area,nodes,element)
do i_elm=1,element_list%n_elements
  area = 0.0
  element = element_list%element(i_elm)
  do m=1,n_vertex_max
    nodes(m) = node_list%node(element%vertex(m))
  enddo

  ! Set up gauss points in this element
  x_s = 0.d0; x_t = 0.d0; y_s = 0.d0; y_t = 0.d0
  do i=1,n_vertex_max
    do j=1,n_degrees
      do ms=1, n_gauss
        do mt=1, n_gauss
          x_s(ms,mt) = x_s(ms,mt) + nodes(i)%x(1,j,1) * element%size(i,j) * H_s(i,j,ms,mt)
          x_t(ms,mt) = x_t(ms,mt) + nodes(i)%x(1,j,1) * element%size(i,j) * H_t(i,j,ms,mt)

          y_s(ms,mt) = y_s(ms,mt) + nodes(i)%x(1,j,2) * element%size(i,j) * H_s(i,j,ms,mt)
          y_t(ms,mt) = y_t(ms,mt) + nodes(i)%x(1,j,2) * element%size(i,j) * H_t(i,j,ms,mt)
        enddo
      enddo
    enddo
  enddo

  ! Perform gauss integration of LHS
  do ms=1, n_gauss
    do mt=1, n_gauss
      wst = wgauss(ms)*wgauss(mt)
      xjac =  x_s(ms,mt)*y_t(ms,mt) - x_t(ms,mt)*y_s(ms,mt)
      area = area + real(xjac * wst, 4)
    enddo
  enddo

  do j=1,size(counts,3)
    do k=1,n_tor
      ivar = (j-1)*n_tor + k
      do l=1,nsub
        do m=1,nsub
          ! Calculate surface of element to convert to density
          ! Not super accurate since this is done on a per-element level
          ! It would be better to do this per subelement
          ielm  = (i_elm-1)*nsub*nsub+(l-1)*nsub+m
          inode = (ielm-1)*4 + 1
          scalars(inode,ivar)   = counts(ielm,k,j)/area
          scalars(inode+1,ivar) = counts(ielm,k,j)/area
          scalars(inode+2,ivar) = counts(ielm,k,j)/area
          scalars(inode+3,ivar) = counts(ielm,k,j)/area
        end do
      end do
    end do
  end do
end do ! n_elements
!$omp end parallel do

! ------------- Write to VTK
call write_vtk(filename,xyz,&
  ien, etype,&
  scalar_names,scalars,&
  vector_names,vectors)
deallocate(scalars,vectors,scalar_names,vector_names)
end subroutine write_particle_counts_to_vtk
end module mod_count_particles

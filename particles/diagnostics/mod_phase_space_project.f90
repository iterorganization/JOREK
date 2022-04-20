!> Module containing functions for phase-space diagnostics
!! Wiki page with more theoretical information and examples:
!!  https://www.jorek.eu/wiki/doku.php?id=particles_phase_space&s[]=phase&s[]=space
!!
!! Projecting particles on rectangular uniform grids for arbitrary grid coordinates and values
!! Using kernels such that correct integration is guaranteed (up to discrete effects)
!!  See https://en.wikipedia.org/wiki/Kernel_density_estimation for more information 
!! Usage can include diagnosing 4D distribution functions (example in particles/examples/tae_phase_space_project.f90
!!  with example analysis script in util/phase_space_diagnostic_example.py This'll plot power-exchange, integrate 
!!  power exchange to show kernel integration, show an interactive 4D distribution function (click the poloidal plot
!!  for local distribution function) and finally print the total amount of particles (obtained by 4D kernel integration).
!! Other example is using power exchange diagnostic in vpar, mu phase space.
!! Supports both event-based calling as the projection on the jorek grids (after having defined 
!!  the appropriate functions to calculate the grid quantities) and a manual projection (more
!!  useful inside an already existing particle loop)
!!
!! Interface:
!! call with(sim,phase_space_projection) will do the projection on the grids
!! call output_phase_project(this,ino,output_grids_in): outputting for HDF5
!! call project_single_particle_x(this,x_in,value_arr,proj_value): projecting single particle on the values arrays

module mod_phase_space_project
  use mod_io_actions
  use mod_particle_sim
  use mod_particle_types
  use mpi_mod
  use mod_project_particles
  use hdf5
  implicit none

  private
  public phase_space_projection, proj_ndim_f, proj_ndim_f_interface
  public new_phase_space_projection, output_phase_project, project_single_particle_x
  
!> A function that will provide ndim values back. This is preferrable to a list of functions that each
!! provide one value back because in many cases an expensive interpolation has to be done which can be reused for different coordinates
!! Examples inclue a vpar, mu plot which both use the B field at the particle location. Similar to proj_f_interface in mod_project_particles
  interface
    function proj_ndim_f_interface(ndim, sim,group, particle)
      import particle_sim, particle_base, phase_space_projection
      type(particle_sim),           intent(in) :: sim
      integer,                      intent(in) :: ndim,group
      class(particle_base),         intent(in) :: particle
      real*8                                   :: proj_ndim_f_interface(ndim)
    end function proj_ndim_f_interface
  end interface

  type :: proj_ndim_f
    procedure(proj_ndim_f_interface), nopass, pointer :: f
    integer                                           :: group ! group to apply function to
  end type proj_ndim_f

  interface proj_ndim_f
    module procedure new_proj_ndim_f !< The constructor for this type
  end interface proj_ndim_f
  
  



  !> Type definition. In principle should not be accesed outside of this module.
  type, extends(io_action) :: phase_space_projection

    ! As the amount of dimensions is not specified and the resolution is not necessarily the same in all directions
    ! one large 1D array will be used both for the n-dimensional grid and the meshgrid. The indices then can be calculated seperately.
    ! where possible (and useful) the conversion is precalculated (basically mixed-radix arrays )
    real*8,       dimension(:),   allocatable         :: grids,dx
    real*8,       dimension(:),   allocatable         :: values

    ! Dimensions and size of values array.
    integer                                           :: ndim, val_size
    ! Resolution of each dimensions, index in grids array where the grid starts for each dimension
    ! and increments of the indices of each dimension in the values array.
    integer,      dimension(:),   allocatable         :: res, previndex, mult_index

    ! Functions used to calculate the projected quantity and the grids/coordinates on which
    ! it is projected
    type(proj_f)                                      :: f_proj
    type(proj_ndim_f)                                  :: f_grids

    ! Quantities for defined particle shapes. 
    ! support(ndim): support in number of grid points in each dimension
    ! prevsupp(ndim): same as previndex, but for support array in calculating kernel
    ! multsupp(ndim): same as mult_index, but for support array in calculating kernel\
    ! bandwidths(ndim): shape of particle in each dimension
    ! totsupport: total amount of support points
    ! sumsupport: sum of support (again calculating kernel)
    ! my_id: MPI ID
    ! totsupp_to_nD(ndim,totsupport): precalculated matrix for converting total support
    !  iteration number to ndim values for ndim indices (for calculating kernel, 20% improvement)
    ! functions_present: describes whether functions are given for event-based calling
    integer,     dimension(:),    allocatable         :: support,prevsupp, multsupp
    real*8,      dimension(:),    allocatable         :: bandwidths
    integer                                           :: totsupport,sumsupport, my_id
    integer,     dimension(:,:),  allocatable         :: totsupp_to_nD
    logical                                           :: functions_present
    !character(len=80)                                 :: basename

  contains
    procedure :: do => project_phase_space
  end type phase_space_projection




contains



!> Constructor for the function that provides ndim values back.
function new_proj_ndim_f(f, group)
  type(proj_ndim_f)                                     :: new_proj_ndim_f
  procedure(proj_ndim_f_interface), pointer, intent(in) :: f
  integer, intent(in)                                   :: group
  new_proj_ndim_f%f => f
  new_proj_ndim_f%group = group
end function new_proj_ndim_f



!> Constructor of the phase-space projection.
!!   ndim             : amount of dimensions
!!   res(ndim)        : resolution (amount of points) of the grids in each dimension
!!   start(ndim)      : starting points of each grid
!!   end(ndim)        : end points of each grid
!!   bandwidths(ndim) : ndim bandwidths (shape extent) in each dimension (optional, but highly recommended. Otherwise only nearest-neighbour)
!!   f_proj           : function that calculates the projected quantity (defined in mod_particle_projection)
!!   f_grids          : function that calculates the positions in each dimension
!!
!! it is recommended to have at least support >= 3 to have integration (it warns) reasonably correct.
!! For many dimensions discretisation errors can accumulate.
function new_phase_space_projection(ndim,res,start,end,bandwidths,f_proj,f_grids,basename) result(new)
  type(phase_space_projection)                          :: new
  integer,                     intent(in)               :: ndim
  integer,                     intent(in)               :: res(ndim)
  real*8,                      intent(in)               :: start(ndim), end(ndim)
  real*8,                      intent(in), optional     :: bandwidths(ndim)
  type(proj_f),                intent(in), optional     :: f_proj
  type(proj_ndim_f),           intent(in), optional     :: f_grids
  character(len=*),            intent(in), optional      :: basename
  integer                                               :: it, jt, gridpoints, valuepoints,j,mult_tmp, my_id
  
  call MPI_COMM_RANK(MPI_COMM_WORLD,my_id,it)
  new%my_id = my_id
  
  ! Allocating the members of the phase-space projection. Have to allocate here due to varying ndim.
  allocate(new%res(ndim))
  allocate(new%previndex(ndim))
  allocate(new%mult_index(ndim))
  allocate(new%dx(ndim))
  if(present(basename))then
    new%basename = trim(basename)//"_"
  else
    new%basename = "proj_"
  endif
  
  if(my_id .eq. 0) then
    write(*,"(A,I1,A)") " PARTICLES: basename: Phase space projection with ",ndim, " dimensions."
    if(ndim>7) then
      ! I cannot foresee a situation where this would come up, but for completeness sake.
      write(*,*) "PARTICLES: ", ndim, "dimensions is not supported by some compilers. "
    endif
  endif

  
  new%ndim = ndim
  new%res = res
  if(present(f_proj) .and. present(f_grids)) then
    new%f_proj  = f_proj
    new%f_grids = f_grids
    if(my_id.eq. 0) write(*,*)"PARTICLES: Initialised projection functions"
    new%functions_present = .true.
  else
    if(my_id.eq. 0) then
      write(*,*)"PARTICLES: No projection function initialised. Disabling event-based calling..."
    endif
    new%functions_present = .false.
  endif

  ! Calculate the total values needed for each dimension
  gridpoints=sum(res)
  ! Calculate the total grid points (on a meshgrid)
  valuepoints=product(res)

  ! Array size of 1D array containing all values.
  new%val_size=valuepoints

  ! This is the array which contains the starting point of each grid.
  new%previndex=0

  ! Allocated on all mpi processes as this is not very large most of the time, reduction can be done later
  ! For many dimensions, this can be prohibitive. E.g. a 4D grids with (90,95,50,75) is already ~256MB
  allocate(new%grids(gridpoints))
  allocate(new%values(valuepoints))
  new%values = 0.d0
  new%grids  = 0.d0

  ! Pre-calculate dimensions multiplication factor. First dimension varies the slowest,
  ! last the quickest. Therefore, the index in the values array that is incremented by
  ! the first dimension increasing one is the multiplication of all other dimension
  ! and so on.
  mult_tmp=valuepoints
  do it=1,ndim
    mult_tmp=mult_tmp/res(it)
    new%mult_index(it)=mult_tmp
    if(it>1) then
      new%previndex(it)=sum(res(1:it-1)) ! Starting point of each of the grids in the grids array
    endif

    ! Fill the grids uniformly
    do j=1,res(it)
      new%grids(new%previndex(it)+j)=start(it)+(end(it)-start(it))/(res(it)-1)*(j-1)

    enddo ! Grids loop
    new%dx(it)=new%grids(new%previndex(it)+2)-new%grids(new%previndex(it)+1)
  enddo   ! Dimension loop

  allocate(new%bandwidths(new%ndim))
  allocate(new%support(new%ndim))
  allocate(new%prevsupp(new%ndim))
  allocate(new%multsupp(new%ndim))
  if (present(bandwidths)) then
    if(my_id .eq.0) write(*,*) "PARTICLES: Custom particle shape given."
    new%bandwidths = bandwidths
  else
    if(my_id .eq.0) write(*,*) "PARTICLES: Bandwidths not given, using default shape of 4 grid spacings."
    new%bandwidths = 4.d0*new%dx
  endif
  new%prevsupp=0
  do it=1,ndim
    new%support(it)=floor(new%bandwidths(it)/new%dx(it))+1!

    if (it > 1 ) then
      new%prevsupp(it)=sum(new%support(1:(it-1)))
    endif

    if(new%support(it) < 3 .and. my_id .eq. 0) then
      write(*,"(A,I1,A)") " PARTICLES: Support for dimension ", it, " is very low. Will not integrate or project correctly."
    endif
    if(my_id .eq. 0) then
      write(*,"(A,I1,A,E11.4,A,I4)") " PARTICLES: In dimension ",it," the particle width is ",new%bandwidths(it), " giving a support of",new%support(it)
    endif
  enddo
  new%totsupport = product(new%support)
  new%sumsupport = sum(new%support)
  mult_tmp = new%totsupport
  do it=1,ndim
    mult_tmp=mult_tmp/new%support(it)
    new%multsupp(it)=mult_tmp
  enddo

  allocate(new%totsupp_to_nD(new%ndim,new%totsupport))

  ! Precaculating conversion between iteration of total support to ndim indices. 
  do it=1,new%totsupport
    new%totsupp_to_nD(:, it) = support_to_nD(new,it)
  enddo

  if(my_id .eq. 0 ) write(*,*) "PARTICLES: Constructed phase space projection"

end function new_phase_space_projection

!> Subroutine for projecting the whole particle distribution
!! on the grids. Not suitable for averaging over single particle
!! motion. Usage can include projecting initialization quantities directly.
!! As the do member of the io-action phase_space_projection points to this
!! function, it is called if calling with(sim, phase_space_event).
subroutine project_phase_space(this, sim, ev)
  use mod_event
  class(phase_space_projection),  intent(inout)           :: this
  type(particle_sim),             intent(inout)           :: sim
  type(event),                    intent(inout), optional :: ev
  real*8                                                  :: x_part(this%ndim), proj_value
  real*8, dimension(:), allocatable                       :: val_tmp
  integer                                                 :: i, i_group_proj, i_group_grids,ierr


  allocate(val_tmp,source=this%values)
  val_tmp = 0.d0

  i_group_proj  = this%f_proj%group
  i_group_grids = this%f_proj%group
  if(i_group_proj .ne. i_group_grids) then
    write(*,*) "PARTICLES: Error: different groups for value projection and grid projection"
    call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
  endif
  if(this%functions_present) then
#ifdef __GFORTRAN__
    !$omp parallel do default(shared) &
#else
    !$omp parallel do default(none) &
    !$omp shared(sim,this,i_group_proj)    &
#endif
    !$omp private(x_part,proj_value,i) &
    !$omp schedule(dynamic,10) &
    !$omp reduction(+:val_tmp)
    do i=1,size(sim%groups(i_group_proj)%particles,1)
      proj_value = this%f_proj%f(sim,i_group_proj,sim%groups(i_group_proj)%particles(i))
      x_part     = this%f_grids%f(this%ndim,sim,i_group_proj,sim%groups(i_group_proj)%particles(i))
      call project_single_particle_x(this,x_part,val_tmp,proj_value)
    enddo
    !$omp end parallel do
    this%values = val_tmp
  else
    if(sim%my_id .eq. 0 ) write(*,*) "PARTICLES: Phase space projection functions not given, returning 0."
    this%values = 0.d0
  endif
end subroutine project_phase_space

!> Output to h5 file. Structure: /values for the meshgrid-evaluated values. /grids/grid_i for 1D grids /grids/mgrid_i for ndim meshgrids
!! (can be immediately plotted using these meshgrids)
subroutine output_phase_project(this,ino,output_grids_in)
  class(phase_space_projection), intent(inout)     :: this
  integer, intent(in)                              :: ino
  logical, intent(in), optional                    :: output_grids_in
  real*8, dimension(:), allocatable                :: val_output
  real*8, dimension(:,:), allocatable                :: grid_mesh
  character(len=1024)                              :: filename
  integer                                          :: my_id, ierr,i,index_arr_tmp(this%ndim)
  integer(HID_T)                                   :: file_id, group_id_grid,dspace,dset_id
  integer                                          :: ierrhdf5
  integer                                          :: it,j
  CHARACTER(len=8)                                 :: tmp_name
  integer                                          :: res_tmp(7),order_tmp(7),res_tmp2(this%ndim) ! Change this if more dimensions needed...
  logical                                          :: output_grids
  !  but this should never come up as particles live in 7D at most.

  if(present(output_grids_in)) then
    output_grids = output_grids_in
  else
    output_grids=.false.
  endif

  ! Output check for too large dimensional arrays for Fortran to handle.
  if (this%ndim > 7) then
    write(*,*) "Can't output for ndim > 7, exiting output"
    return
  endif

  ! For the reshaping of the arrays into ndim arrays. E.g. if 2D, it would look like [res1,res2,0,0,0,0,0]
  ! with HDF5 ignoring the 0 dimensions.
  res_tmp=0
  res_tmp(1:this%ndim)=(/ (this%res(this%ndim-I),I=0,this%ndim-1) /)

  call MPI_COMM_RANK(MPI_COMM_WORLD, my_id, ierr)

  ! Allocate only properly on root process
  if (my_id .eq. 0) then
    allocate(val_output, source=this%values)
    allocate(grid_mesh(size(this%values),this%ndim))
  else
    allocate(val_output(0))
    allocate(grid_mesh(0,0))
  endif

  ! Reduce output to the root process for output
  val_output=0.d0
  call MPI_Reduce(this%values,val_output,size(this%values), MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

  ! Output to HDF5
  if(my_id .eq. 0 ) then

    !Calculate meshgrids
    do it=1, size(this%values)
      index_arr_tmp = calc_reverse_index_phase_proj(this,it)
      do j=1, this%ndim
        grid_mesh(it,j) = this%grids(index_arr_tmp(j)+this%previndex(j))
      enddo
    enddo


    ! HDF5 file creation
    call h5open_f(ierrhdf5)
    write(filename,"(A,i5.5,A)") trim(this%basename) ,ino, ".h5"
    call H5Fcreate_f(filename,H5F_ACC_TRUNC_F, file_id, ierrhdf5)
    if(output_grids)then
      call h5gcreate_f(file_id, "grids", group_id_grid, ierrhdf5)
    endif


    res_tmp2=res_tmp(1:this%ndim)

    if(output_grids) then
      do it=1,this%ndim

        ! Output all the 1D grids in each dimensions under the /grids/ group
        write(tmp_name,fmt="(A,I1)") "grid_",it
        call h5screate_simple_f(1, [int(this%res(it),kind=HSIZE_T)], dspace, ierr)!, &
        call h5dcreate_f(group_id_grid, tmp_name, H5T_NATIVE_DOUBLE, dspace, &
                       dset_id, ierrhdf5)
        call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,this%grids((this%previndex(it)+1):(this%previndex(it)+this%res(it))),[int(this%res(it),kind=HSIZE_T)],ierr)
        call h5dclose_f(dset_id, ierr)
        call h5sclose_f(dspace,ierr)

        !Same but meshgrid
        call h5screate_simple_f(this%ndim, int(res_tmp2,kind=HSIZE_T),dspace,ierr )
        write(tmp_name,fmt="(A,I1)") "mgrid_",it
        call h5dcreate_f(group_id_grid, tmp_name, H5T_NATIVE_DOUBLE, dspace, &
                       dset_id, ierrhdf5)

        call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,RESHAPE(grid_mesh(:,it),res_tmp),int(res_tmp2,kind=HSIZE_T),ierr)
        call h5dclose_f(dset_id,ierr)
        call h5sclose_f(dspace,ierr)

      enddo ! dimension loop grids
    endif
    ! Output meshgrids of the value
    call h5screate_simple_f(this%ndim, int(res_tmp2,kind=HSIZE_T), dspace, ierr)
    call h5dcreate_f(file_id,"values", H5T_NATIVE_DOUBLE, dspace, dset_id, ierrhdf5)

    ! The trick to gaining arbitrary dimensional arrays into HDF5 is by using reshape immediately in the function
    ! as this saves us a maximum-dimension array allocation in the code.
    ! Beware: stack limits for reshape! set "ulimit -s unlimited"
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE,RESHAPE(val_output(1:size(val_output)),res_tmp),int(res_tmp2,kind=HSIZE_T),ierr)
    call h5dclose_f(dset_id,ierr)
    call h5sclose_f(dspace,ierr)
    if(output_grids) then
      call h5gclose_f(group_id_grid, ierr)
    endif
    call h5fclose_f(file_id,  ierrhdf5)
    call h5close_f(ierrhdf5)
    write(*,* ) "PARTICLES: Written Phase Space Projection to ", trim(filename)
  endif  ! mpi
  deallocate(val_output)
  deallocate(grid_mesh)
end subroutine output_phase_project

!> Subroutine for projecting one particle to the value_arr
!! this: phase space projection for the arrays needed to project
!! x_in: location of particle on the dimensions to project
!! value_arr: array to be projected on. Should be the same size as this%values
!! proj_value: projection value (i.e. if only particle weight, then just particle%weight)
subroutine project_single_particle_x(this,x_in,value_arr,proj_value)
  class(phase_space_projection), intent(in)        :: this
  real*8,                        intent(in)        :: x_in(this%ndim)
  real*8,                        intent(inout)     :: value_arr(this%val_size) ! Large!
  real*8,                        intent(in)        :: proj_value
  integer                                          :: i, j, main_ind, min_ind,min_ind_cont(this%ndim), mesh_tmp,ind_mesh_tmp, i_phase
  integer                                          :: index_tmp, indices_phase_grids_tmp(this%ndim)
  real*8                                           :: dx, x, minx, maxx, val_phase_grids_tmp, supp_weight(this%sumsupport), distance

  index_tmp=-1

  supp_weight=0.d0

  ! whole loop ~ <5% execution time
  do i=1,this%ndim
    dx = this%dx(i)
    minx=this%grids(this%previndex(i)+1)
    maxx=this%grids(this%previndex(i)+this%res(i))
    x = x_in(i)
    ! Not test if main x is in grid yet in order for particles outside the considered grid to
    ! contribute with their support (i.e. considering R=9.5-10.0, particle at 9.45 with 0.1
    ! bandwidth will still contribute slightly)

    ! Calculate minimum index (will be < 0 if particle outside of grid, but will be tested later)
    ! 1-indexing, should be 1 higher to indicate a true minimum index in the array. But in the support loop this is taken care of.
    min_ind = ceiling( (x-this%bandwidths(i)/2.d0-minx)/dx)
    min_ind_cont(i)=min_ind

    !Now, we add all the grid points in the particles shape into the grids
    do j=1,this%support(i)
      ! Only calculate weight if it makes sense
      if(min_ind+j>0 .and. min_ind+j < this%res(i)+1) then

        ! Distance w.r.t particle centre.
        distance=abs((x-this%grids(this%previndex(i)+min_ind+j))/this%bandwidths(i)*2.d0)

        if (distance > 1.d0) then
          supp_weight(j+this%prevsupp(i))=0.d0
        else
          ! This is the kernel, i.e. https://en.wikipedia.org/wiki/Kernel_(statistics)
          ! but we can only use finite support ones naturally (i.e. no Gaussian)
          ! Will be scaled later to ensure integration to 1.
          ! These are all valid kernels and there is no major difference in performance
          ! for diagnostics most likely
          !Uniform
          !supp_weight(j+this%prevsupp(i))=0.5d0
          !Linear
          !supp_weight(j+this%prevsupp(i))=1-abs(distance)
          !Epachnikov
          !supp_weight(j+this%prevsupp(i))=3.d0/4.d0*(1-distance**2.d0)
          !Quartic
          !supp_weight(j+this%prevsupp(i))=15.d0/16.d0*(1-distance**2)**2
          !Triweight
          supp_weight(j+this%prevsupp(i))=35.d0/32.d0*(1-distance**2)**3
        endif ! distance < 1.0
      endif ! min_ind+j > 0 && min_ind+j < this%res(i)+1
    enddo !support points
  enddo ! n_dim

  ! Now for every point in the total support calculate the total kernel and output that on the grids.
  do i=1,this%totsupport

    val_phase_grids_tmp=1.d0
    ! Total weight of the support point i is calculated by the product of the weights in all dimensions.
    do j=1,this%ndim
      val_phase_grids_tmp=val_phase_grids_tmp*supp_weight(this%totsupp_to_nD(j,i)+this%prevsupp(j))*1/this%bandwidths(j)*2.d0 !Bandwidth = 2*bandwidth_wiki !~13% execution time
    enddo
    ! We need to check that this point is in fact on the grids
    if (all(this%totsupp_to_nD(:,i)+min_ind_cont >0) .and. all(this%totsupp_to_nD(:,i)+min_ind_cont < this%res+1)) then ! 13% execution time
      index_tmp=calc_index_phase_proj(this, this%totsupp_to_nD(:,i)+min_ind_cont) ! ~ 5% execution time
      ! 43% execution time 
      value_arr(index_tmp)=value_arr(index_tmp)+1.d0*val_phase_grids_tmp*proj_value
    endif
  enddo ! support points
end subroutine project_single_particle_x


pure function support_to_nD(this, index_in) result(index_out)
  class(phase_space_projection),  intent(in)       :: this
  integer,                        intent(in)       :: index_in
  integer                                          :: index_out(this%ndim)
  integer                                          :: index_tmp,it
  index_out = 0
  index_tmp=index_in-1
  do it=1,this%ndim
    index_out(it) = index_tmp/this%multsupp(it)+1
    index_tmp = modulo(index_tmp, this%multsupp(it))
  enddo
end function support_to_nD
!> To go from 1D large meshgrid to ndim meshgrids function
function calc_index_phase_proj(this, index_arr) result(index_values)
  class(phase_space_projection), intent(in) :: this
  integer, intent(in)                       :: index_arr(this%ndim)
  integer                                   :: it
  integer                                   :: index_values,index_tmp

  index_tmp=1
  do it=1,this%ndim
    index_tmp=index_tmp+this%mult_index(it)*(index_arr(it)-1)
  enddo
  index_values=index_tmp
end function calc_index_phase_proj
!> To go from index in 1D array to ndim indices.
function calc_reverse_index_phase_proj(this, index_values) result(index_arr)
  class(phase_space_projection), intent(in) :: this
  integer                                   :: index_arr(this%ndim)
  integer                                   :: it, tmp,tmp2
  integer                                   :: index_values,index_tmp

  index_arr=0
  index_tmp=index_values-1
  do it=1,this%ndim
    tmp = floor(real(index_tmp)/real(this%mult_index(it))) ! First stride. E.g. 10 -> 2,2,2,2 -> 0,
    tmp2 = modulo(index_tmp, this%mult_index(it)) !Remainder after subtracting first stride.
    index_arr(it)=tmp+1
    index_tmp = tmp2

  enddo

end function calc_reverse_index_phase_proj

end module mod_phase_space_project

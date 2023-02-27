module mod_catalyst_adaptor
#ifdef USE_CATALYST

  use, intrinsic :: iso_c_binding
  use basis_at_gaussian
  use nodes_elements
  use mod_interp
  use mod_parameters, only: n_var, variable_names

  implicit none

  public :: catalyst_adaptor_initialise
  public :: catalyst_adaptor_execute
  public :: catalyst_adaptor_finalise
  public :: catalyst_adaptor

  private

  character(kind=c_char, len=65536), public:: catalyst_scripts !< Path to Catalyst scripts (multiple paths separated with ':')
  integer(c_int), public                   :: catalyst_nsub = 5 !< Number of subdivisions of each JOREK element
  integer                                  :: n_scalars !< The number of scalar variables passed to Catalyst
  integer                                  :: nnos !< Number of nodes in the Catalyst grid
  integer                                  :: nel !< Number of cells in the Catalyst grid
  integer                                  :: nnoel = 4 !< Number of points to define a cell
  real(c_float), dimension(:), allocatable :: coords_s !< local subdivision coordinate s
  real(c_float), dimension(:), allocatable :: coords_t !< local subdivision coordinate t
  integer                                  :: i_plane = 1 !< which plane to interpolate at

  interface

    subroutine catalyst_adaptor_initialise(a_catalyst_scripts) bind(C)
      use, intrinsic :: iso_c_binding
      character(kind=c_char), intent(in) :: a_catalyst_scripts 
    end subroutine catalyst_adaptor_initialise

    subroutine catalyst_adaptor_execute(a_step_index, a_time) bind(C)
      use, intrinsic :: iso_c_binding
      integer(c_int), intent(in) :: a_step_index
      real(c_double), intent(in) :: a_time
    end subroutine catalyst_adaptor_execute

    subroutine catalyst_adaptor_finalise() bind(C)
    end subroutine catalyst_adaptor_finalise

    ! empty function to get the dependency generator to work with catalyst_adaptor.cpp
    subroutine catalyst_adaptor() bind(C)
    end subroutine catalyst_adaptor

  end interface

  contains

    subroutine catalyst_get_params(a_nsub, a_n_elements, a_n_scalars) bind(C)
      use, intrinsic :: iso_c_binding
      integer(c_int), intent(out) :: a_nsub
      integer(c_int), intent(out) :: a_n_elements
      integer(c_int), intent(out) :: a_n_scalars

      a_nsub = catalyst_nsub
      a_n_elements = element_list%n_elements
      ! Just pass all variables as scalars for now (can change this later)
      n_scalars = n_var
      a_n_scalars = n_scalars
    end subroutine catalyst_get_params

    subroutine catalyst_get_scalar_name(a_scalar_name, a_iscalar) bind(C)
      use, intrinsic :: iso_c_binding
      character(kind=c_char), intent(out), dimension(36) :: a_scalar_name
      integer(c_int), intent(in) :: a_iscalar

      integer :: ichar
      character(kind=c_char,len=12) :: trimmed_name
      trimmed_name = trim(variable_names(a_iscalar)) // c_null_char
      do ichar=1,len(trimmed_name)
        a_scalar_name(ichar) = trimmed_name(ichar:ichar)
      enddo
    end subroutine catalyst_get_scalar_name

    subroutine catalyst_interp_grid(a_nnos, a_nel, a_coords_R, a_coords_Z, a_cell_points) bind(C)
      use, intrinsic :: iso_c_binding
      integer(c_int), intent(in) :: a_nnos
      integer(c_int), intent(in) :: a_nel
      real(c_float), intent(inout), dimension(a_nnos), target :: a_coords_R
      real(c_float), intent(inout), dimension(a_nnos), target :: a_coords_Z
      integer(c_int), intent(inout), dimension(nnoel * a_nel) :: a_cell_points

      integer :: i, j, ielm, inode, k
      real*8 :: s, t, R, Z

      nnos = a_nnos
      nel = a_nel
      inode = 0
      ielm = 0

      if (.not. allocated(coords_s)) allocate(coords_s(a_nnos))
      if (.not. allocated(coords_t)) allocate(coords_t(a_nnos))

      ! Create points for each element
      do i=1,element_list%n_elements
        do j=1,catalyst_nsub
          s = float(j-1)/float(catalyst_nsub-1)
          ! Create nsub^2 points per element at regularly spaced intervals
          do k=1,catalyst_nsub
            t = float(k-1)/float(catalyst_nsub-1)
            call interp_RZ(node_list,element_list,i,s,t,R,Z)
            inode = inode + 1
            a_coords_R(inode) = real(R, c_float)
            a_coords_Z(inode) = real(Z, c_float)
            coords_s(inode) = real(s, c_float)
            coords_t(inode) = real(t, c_float)
          enddo
        enddo

        ! Calculate connectivity of each subelement
        do j=1,catalyst_nsub-1
          do k=1,catalyst_nsub-1
            ielm = ielm + 1
            ! Hopefully Catalyst expects the same convention as VTK
            ! Conduit doesn't specify a convention
            a_cell_points(4*ielm-3) = inode - catalyst_nsub*catalyst_nsub + catalyst_nsub*(j-1) + k-1 ! 0 based indices as for VTK
            a_cell_points(4*ielm-2) = inode - catalyst_nsub*catalyst_nsub + catalyst_nsub*(j  ) + k-1
            a_cell_points(4*ielm-1) = inode - catalyst_nsub*catalyst_nsub + catalyst_nsub*(j  ) + k
            a_cell_points(4*ielm  ) = inode - catalyst_nsub*catalyst_nsub + catalyst_nsub*(j-1) + k
          enddo
        enddo
      enddo
    end subroutine catalyst_interp_grid

    subroutine catalyst_interp_scalar(a_scalar, a_iscalar) bind(C)
      use, intrinsic :: iso_c_binding
      real(c_float), intent(inout), dimension(nnos) :: a_scalar
      integer(c_int), intent(in)                    :: a_iscalar

      integer :: i !< JOREK element index
      integer :: inode, i_tor
      real*8 :: s, t, P, P_s, P_t, P_ss, P_st, P_tt

      i = 0
      a_scalar = 0.0
      do inode=1,nnos
        ! Increment the JOREK element index by 1 every nsub^2 nodes
        if (modulo(inode - 1, catalyst_nsub * catalyst_nsub) .eq. 0) then
          i = i + 1
        endif
        s = coords_s(inode)
        t = coords_t(inode)
        do i_tor=1,n_tor
          call interp(node_list,element_list,i,a_iscalar,i_tor,s,t,P,P_s,P_t,P_st,P_ss,P_tt)
          a_scalar(inode) = a_scalar(inode) + real(P * HZ(i_tor,i_plane), c_float)
        enddo ! i_tor
      enddo ! inode
      
    end subroutine catalyst_interp_scalar

#endif 
end module mod_catalyst_adaptor
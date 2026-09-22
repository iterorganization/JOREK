!> Module for creating an RTree of elements and searching in this tree
module mod_element_rtree
use iso_c_binding
use data_structure
use mpi
implicit none
public populate_element_rtree, nearby_elements, elements_containing_point, rtree_initialized, tree_slices

logical :: rtree_initialized = .false.
integer, parameter :: tree_slices = n_plane
integer, parameter :: n_between = 16

#if STELLARATOR_MODEL
integer, parameter :: ND = 3
#else
integer, parameter :: ND = 2
#endif

interface
  subroutine RZ_minmax(node_list, element_list, i_elm, Rmin, Rmax, Zmin, Zmax, i_plane_query_in)
    use data_structure, only: type_node_list, type_element_list
    implicit none
    type(type_node_list), intent(in)    :: node_list
    type(type_element_list), intent(in) :: element_list
    integer, intent(in)                 :: i_elm
    real*8 , intent(out)                :: Rmin, Rmax, Zmin, Zmax
    integer, intent(in), optional       :: i_plane_query_in
  end subroutine RZ_minmax
  subroutine RZP_minmax(node_list, element_list, i_elm, phi, Rmin, Rmax, Zmin, Zmax)
    use data_structure, only: type_node_list, type_element_list
    implicit none
    type(type_node_list), intent(in)    :: node_list
    type(type_element_list), intent(in) :: element_list
    integer, intent(in)                 :: i_elm
    real*8 , intent(in)                 :: phi
    real*8 , intent(out)                :: Rmin, Rmax, Zmin, Zmax
  end subroutine RZP_minmax
  !> Name is element_rtree to match filename `element_rtree.cpp`.
  !> `void PopulateTree(int n_elms, int *sizes, double **min, double **max)
  subroutine element_rtree(n_elms, sizes, min, max) bind(C, name="PopulateTree")
    import :: C_INT, C_PTR
    integer(C_INT), value, intent(in) :: n_elms
    integer(C_INT), intent(in), dimension(*) :: sizes
    type(C_PTR), intent(in), dimension(*) :: min
    type(C_PTR), intent(in), dimension(*) :: max
  end subroutine element_rtree
  !> `int NumElementsInRect(double minx, double miny, double maxx, double maxy)`
  function num_elements_in_rect(min, max) bind(C,name='NumElementsInRect')
    import C_DOUBLE, C_INT
    real(C_DOUBLE), intent(in), dimension(*)         :: min, max
    integer(C_INT) :: num_elements_in_rect
  end function num_elements_in_rect
  !> `int ElementsInRect(double minx, double miny, double maxx, double maxy, int *ielm)`
  function elements_in_rect(min, max, ielm) bind(C,name='ElementsInRect')
    import C_DOUBLE, C_INT
    real(C_DOUBLE), intent(in), dimension(*)          :: min, max
    integer(C_INT), intent(inout), dimension(*) :: ielm
    integer(C_INT) :: elements_in_rect
  end function elements_in_rect
end interface
type :: type_box_container
   real(C_DOUBLE), allocatable :: min_vals(:) ! Size: n_slices * 3
   real(C_DOUBLE), allocatable :: max_vals(:) ! Size: n_slices * 3
end type type_box_container

contains


subroutine populate_element_rtree(node_list, element_list, use_3D_rtree)
   use iso_c_binding
   use constants, only: PI
   use data_structure, only: type_node_list, type_element_list
   implicit none

  type(type_node_list),    intent(in) :: node_list
  type(type_element_list), intent(in) :: element_list
  logical, intent(in), optional       :: use_3D_rtree
  logical :: use_3D

! Determine whether to use 3D RTree based on argument or compilation flag
#if !STELLARATOR_MODEL
   use_3D = .false.
#else
   if (present(use_3D_rtree)) then
      use_3D = use_3D_rtree
   else
      use_3D = .true.
   end if
#endif

   if (use_3D) then
      call populate_element_rtree_3D(node_list, element_list)
   else
      call populate_element_rtree_2D(node_list, element_list)
   end if
  ! This is just a dispatcher to the appropriate version (2D or 3D) based on compilation flags.
end subroutine populate_element_rtree

!> Populate the RTree with the squares containing elements
!> x=R, y=Z
!>
!> Only one call to this routine can be made simultaneously! It uses a static
!> data structure under the hood.
subroutine populate_element_rtree_2D(node_list, element_list)
    use iso_c_binding
    use constants, only: PI
    use mpi
    
    type(type_node_list),    intent(in) :: node_list
    type(type_element_list), intent(in) :: element_list

    ! --- C Interface Variables ---
    integer(C_INT) :: n_elms_c
    integer(C_INT), allocatable, target :: sizes_c(:)
    type(C_PTR),    allocatable, target :: min_ptrs(:), max_ptrs(:)
    
    ! --- Data Container ---
    ! Must be target so C_LOC can point to its components
    type(type_box_container), allocatable, target :: element_data(:)

    ! --- Local Variables ---
    logical :: is_init
    integer :: i, n, my_id, ierr
    real*8  :: rmin, rmax, zmin, zmax

    n = element_list%n_elements
    n_elms_c = int(n, C_INT)

    ! 1. Allocate Top-Level Arrays
    allocate(sizes_c(n))
    allocate(min_ptrs(n))
    allocate(max_ptrs(n))
    allocate(element_data(n))
    
    call MPI_Initialized(is_init, ierr)

    if (is_init) then
      call MPI_COMM_RANK(MPI_COMM_WORLD, my_id, ierr)
    else
      my_id = 0
    endif

    if (my_id .eq. 0) then
      write(*,*) "Initializing 2D RTree (Axisymmetric)..."
    endif

    ! 2. Loop over elements
    do i = 1, n
       
       ! Calculate R/Z extents using standard JOREK routine (axisymmetric)
       call RZ_minmax(node_list, element_list, i, rmin, rmax, zmin, zmax)

       ! In 2D, we have exactly 1 box per element
       sizes_c(i) = 1

       allocate(element_data(i)%min_vals(ND))
       allocate(element_data(i)%max_vals(ND))

       ! 3. Define the Box
       ! R and Z come from the element geometry
       element_data(i)%min_vals(1) = real(rmin, C_DOUBLE)
       element_data(i)%min_vals(2) = real(zmin, C_DOUBLE)
       
       element_data(i)%max_vals(1) = real(rmax, C_DOUBLE)
       element_data(i)%max_vals(2) = real(zmax, C_DOUBLE)

#if STELLARATOR_MODEL
       ! Initialise phi component for stellarator geometry.
       element_data(i)%min_vals(3) = real(-PI, C_DOUBLE)
       element_data(i)%max_vals(3) = real(+PI, C_DOUBLE)
#endif

       ! 4. Store Pointers for C++
       min_ptrs(i) = C_LOC(element_data(i)%min_vals)
       max_ptrs(i) = C_LOC(element_data(i)%max_vals)

    end do

    ! 5. Call C++
    call element_rtree(n_elms_c, sizes_c, min_ptrs, max_ptrs)

    rtree_initialized = .true.

  end subroutine populate_element_rtree_2D

recursive subroutine get_element_critical_planes(node_list, element_list, i_elm, planes_out)
   use data_structure, only: type_node_list, type_element_list
   use mod_interp,     only: interp_RZP
   use constants,      only: PI
   implicit none

   ! --- Arguments ---
   type(type_node_list),    intent(in) :: node_list
   type(type_element_list), intent(in) :: element_list
   integer,                 intent(in) :: i_elm
   real*8, allocatable,    intent(out) :: planes_out(:)

   ! --- Local Parameters ---
   real*8, parameter :: tol_res  = 1d-15      ! Tolerance for gradient (Residual)
   real*8, parameter :: max_step = 0.2d0      ! Max Newton step size (Damping)
   integer, parameter :: max_iter = 10
   integer, parameter :: n_guesses_phi = n_coord_tor * n_coord_period * 2   ! Grid resolution for 1D scan
   integer, parameter :: n_guesses_st  = 5                                  ! Grid resolution for s/t in 2D scan
   integer, parameter :: n_guesses_phi_2d = 10                              ! Grid resolution for phi in 2D scan

   ! --- Local Variables ---
   integer :: i, j, k, n_found, iface, icorner, ig_phi, ig_s, ig_t
   real*8  :: period, phi, s, t, delta_sub, val_start, val_end
   real*8, allocatable :: temp_planes(:)

   ! Bracketing variables for 1D
   real*8  :: dphi_grid, phi_L, phi_R, phi_M, val_L, val_R, val_M, seed_guess, dst_grid, st_seed
   logical :: check_L, check_R

   ! Interp variables
   real*8 :: R, R_s, R_t, R_p, R_st, R_ss, R_tt, R_sp, R_tp, R_pp
   real*8 :: Z, Z_s, Z_t, Z_p, Z_st, Z_ss, Z_tt, Z_sp, Z_tp, Z_pp

   allocate(temp_planes(64))

   period = 2.0d0 * PI / real(n_coord_period, 8)
   n_found = 0

   ! Always add domain boundaries (start)
   call add_unique_plane(0.0d0)

   ! =========================================================================
   !  Loop over Target Functions: k=1 (R), k=2 (Z)
   ! =========================================================================
   do k = 1, 2
      ! Loop over Corners (s,t) = (0,0), (1,0), (0,1), (1,1)
      do icorner = 1, 4
         if (icorner == 1) then
            s = 0.0d0; t = 0.0d0
         else if (icorner == 2) then
            s = 1.0d0; t = 0.0d0
         else if (icorner == 3) then
            s = 0.0d0; t = 1.0d0
         else
            s = 1.0d0; t = 1.0d0
         end if

         ! Grid Search / Seeding for Newton
         ! We seed periodically to catch all critical points
         dphi_grid = period / real(n_guesses_phi, 8)
         do ig_phi = 0, n_guesses_phi - 1
            seed_guess = real(ig_phi, 8) * dphi_grid
            call solve_newton_1d(k, s, t, seed_guess)
         end do
      end do
      
      ! Loop over Faces (2D Newton) - refactored to reduce duplication
      dst_grid = 1.0d0 / real(n_guesses_st, 8)

      dphi_grid = period / real(n_guesses_phi_2d, 8)

      do iface = 1, 2
         do ig_s = 0, n_guesses_st - 1
            st_seed = (real(ig_s, 8) + 0.5d0) * dst_grid
            do ig_phi = 0, n_guesses_phi_2d - 1
               seed_guess = real(ig_phi, 8) * dphi_grid
               call solve_newton_2d(k, 1, real(iface - 1, 8), st_seed, seed_guess) ! s_fixed
               call solve_newton_2d(k, 2, real(iface - 1, 8), st_seed, seed_guess) ! t_fixed
            end do
         end do
      end do

      end do

   ! --- Final Sort and Pack ---
   call sort(temp_planes(1:n_found), n_found)

   ! Ensure space for one more
   if (n_found + 1 > size(temp_planes)) call resize_temp_planes()
   n_found = n_found + 1
   temp_planes(n_found) = period 

   allocate(planes_out(n_found + n_between * (n_found - 1)))

   ! Standard Intervals
   k = 0
   do i = 1, n_found - 1
      k = k + 1
      planes_out(k) = temp_planes(i)
      if (n_between > 0) then
         val_start = temp_planes(i)
         val_end   = temp_planes(i+1)
         delta_sub = val_end - val_start
         ! Only add intermediate planes if there is enough space
         if (delta_sub > 1d-10) then
            delta_sub = delta_sub / real(n_between + 1, 8)
            do j = 1, n_between
               k = k + 1
               planes_out(k) = val_start + real(j, 8) * delta_sub
            end do
         end if
      end if
   end do
   ! Add the last plane (period)
   planes_out(k+1) = temp_planes(n_found)

   ! shrink to size
   if (k + 1 < size(planes_out)) then
      planes_out = planes_out(1:k+1)
   end if

   contains

   ! =========================================================================
   !  Solver: 1D Bounded Newton
   !  Finds root of R_p or Z_p starting from p_seed
   ! =========================================================================
   subroutine solve_newton_1d(comp, s_fix, t_fix, p_in)
      integer, intent(in) :: comp ! 1=R, 2=Z
      real*8, intent(in) :: s_fix, t_fix, p_in
      
      real*8 :: p_curr, delta, d1, d2
      integer :: iter

      p_curr = p_in
      
      do iter = 1, max_iter
         call interp_RZP(node_list, element_list, i_elm, s_fix, t_fix, p_curr, &
                           R, R_s, R_t, R_p, R_st, R_ss, R_tt, R_sp, R_tp, R_pp, &
                           Z, Z_s, Z_t, Z_p, Z_st, Z_ss, Z_tt, Z_sp, Z_tp, Z_pp)
                           
         if (comp == 1) then
            d1 = R_p
            d2 = R_pp
         else
            d1 = Z_p
            d2 = Z_pp
         endif
         
         ! Check convergence (residual)
         if (abs(d1) < tol_res) exit
         
         ! Safety for division (d2 approx 0)
         if (abs(d2) < 1d-9) d2 = sign(1d-9, d2)
         
         delta = d1 / d2
         
         ! Damping / Limiting
         if (abs(delta) > max_step) delta = sign(max_step, delta)
         
         p_curr = p_curr - delta
      end do
      
      ! Add result to unique planes (handles normalization)
      call add_unique_plane(p_curr)
      
   end subroutine solve_newton_1d

   ! =========================================================================
   !  Solver: 2D Newton for Faces
   !  Solves { R_p=0, R_t=0 } or { Z_p=0, Z_t=0 } (if s fixed)
   !  Solves { R_p=0, R_s=0 } or { Z_p=0, Z_s=0 } (if t fixed)
   ! =========================================================================
   subroutine solve_newton_2d(comp, type, fix_val, var_seed, p_seed)
      integer, intent(in) :: comp ! 1=R, 2=Z
      integer, intent(in) :: type ! 1=s_fixed (vary t,p), 2=t_fixed (vary s,p)
      real*8, intent(in)  :: fix_val, var_seed, p_seed
      
      real*8 :: p_curr, var_curr
      real*8 :: J11, J12, J21, J22, det, invJ11, invJ12, invJ21, invJ22
      real*8 :: F1, F2, d_var, d_p
      real*8 :: s_loc, t_loc
      integer :: iter
      
      p_curr = p_seed
      var_curr = var_seed
      
      do iter = 1, max_iter
         if (type == 1) then ! s fixed
            s_loc = fix_val
            t_loc = var_curr
         else ! t fixed
            s_loc = var_curr
            t_loc = fix_val
         endif
         
         ! Enforce bounds on s/t [0,1]
         if (var_curr < 0.0d0) var_curr = 0.0d0
         if (var_curr > 1.0d0) var_curr = 1.0d0

         call interp_RZP(node_list, element_list, i_elm, s_loc, t_loc, p_curr, &
                           R, R_s, R_t, R_p, R_st, R_ss, R_tt, R_sp, R_tp, R_pp, &
                           Z, Z_s, Z_t, Z_p, Z_st, Z_ss, Z_tt, Z_sp, Z_tp, Z_pp)
         
         ! We want to find critical point wrt Phi AND the varying coordinate (s or t)
         ! So we solve Gradient = 0.
         ! Variables: x = [var, p]
         
         if (comp == 1) then ! R
            if (type == 1) then ! vary t, p
               F1 = R_t
               F2 = R_p
               ! Jacobian d(F1,F2)/d(t,p)
               J11 = R_tt
               J12 = R_tp
               J21 = R_tp
               J22 = R_pp
            else ! vary s, p
               F1 = R_s
               F2 = R_p
               ! Jacobian d(F1,F2)/d(s,p)
               J11 = R_ss
               J12 = R_sp
               J21 = R_sp
               J22 = R_pp
            endif
         else ! Z
            if (type == 1) then ! vary t, p
               F1 = Z_t
               F2 = Z_p
               ! Jacobian
               J11 = Z_tt
               J12 = Z_tp
               J21 = Z_tp
               J22 = Z_pp
            else ! vary s, p
               F1 = Z_s
               F2 = Z_p
               ! Jacobian
               J11 = Z_ss
               J12 = Z_sp
               J21 = Z_sp
               J22 = Z_pp
            endif
         endif
         
         ! Check Convergence
         if (sqrt(F1**2 + F2**2) < tol_res) exit
         
         ! Solve J * d = F
         det = J11*J22 - J12*J21
         if (abs(det) < 1d-12) det = sign(1d-12, det)
         
         invJ11 =  J22 / det
         invJ12 = -J12 / det
         invJ21 = -J21 / det
         invJ22 =  J11 / det
         
         d_var = invJ11 * F1 + invJ12 * F2
         d_p   = invJ21 * F1 + invJ22 * F2
         
         ! Damping
         if (abs(d_var) > max_step) d_var = sign(max_step, d_var)
         if (abs(d_p)   > max_step) d_p   = sign(max_step, d_p)
         
         var_curr = var_curr - d_var
         p_curr   = p_curr   - d_p
      end do
      
      ! Add result
      call add_unique_plane(p_curr)

   end subroutine solve_newton_2d

   ! =========================================================================
   !  Utility: Add Unique Plane
   ! =========================================================================
   subroutine add_unique_plane(p_in)
      real*8, intent(in) :: p_in
      real*8 :: p_norm
      integer :: i
      logical :: dup

      ! Normalize to [0, 2pi)
      p_norm = p_in - floor(p_in / period) * period
      
      dup = .false.
      do i = 1, n_found
         if (abs(temp_planes(i) - p_norm) < 1d-4) then
            dup = .true.; exit
         endif
         ! Check wraparound proximity
         if (abs(temp_planes(i) - (p_norm + period)) < 1d-4 .or. &
            abs(temp_planes(i) - (p_norm - period)) < 1d-4) then
            dup = .true.; exit
         endif
      end do
      
      if (.not. dup) then
         if (n_found + 1 > size(temp_planes)) call resize_temp_planes()
         n_found = n_found + 1
         temp_planes(n_found) = p_norm
      endif
   end subroutine add_unique_plane

   ! =========================================================================
   !  Utility: Resize Temp Planes Array
   ! =========================================================================
   subroutine resize_temp_planes()
      real*8, allocatable :: temp_arr(:)
      allocate(temp_arr(size(temp_planes) * 2))
      temp_arr(1:size(temp_planes)) = temp_planes
      call move_alloc(temp_arr, temp_planes)
   end subroutine resize_temp_planes

   ! =========================================================================
   !  Utility: Sort
   ! =========================================================================
   subroutine sort(arr, n)
      real*8, intent(inout) :: arr(:)
      integer, intent(in) :: n
      integer :: i, j
      real*8 :: temp
      do i = 1, n-1
         do j = i+1, n
            if (arr(j) < arr(i)) then
               temp = arr(i); arr(i) = arr(j); arr(j) = temp
            end if
         end do
      end do
   end subroutine sort

end subroutine get_element_critical_planes

!> Only one call to this routine can be made simultaneously! It uses a static
subroutine populate_element_rtree_3D(node_list, element_list)
  use iso_c_binding
  use constants, only: PI
  use mpi
  
  type(type_node_list),    intent(in) :: node_list
  type(type_element_list), intent(in) :: element_list

  ! Variables for the C Interface
  integer(C_INT) :: n_elms_c
  integer(C_INT), allocatable, target :: sizes_c(:)
  type(C_PTR),    allocatable, target :: min_ptrs(:), max_ptrs(:)
  
  ! Storage for the actual data (Jagged arrays)
  type(type_box_container), allocatable, target :: element_data(:)

  ! Local logic variables
  integer :: i, j, n, num_planes, num_boxes, my_id, ierr
  real*8, allocatable :: phi_planes(:) ! Array to hold critical angles for current element
  real*8 :: rmin, rmax, zmin, zmax, rmin_next, rmax_next, zmin_next, zmax_next
  real*8 :: phi_start, phi_end
  integer :: idx
  logical :: is_init

  n = element_list%n_elements
  n_elms_c = int(n, C_INT)

  ! Allocate arrays for the C interface
  allocate(sizes_c(n))
  allocate(min_ptrs(n))
  allocate(max_ptrs(n))
  
  ! Allocate the container that holds the jagged actual data
  allocate(element_data(n))

  call MPI_Initialized(is_init, ierr)

  if (is_init) then
    call MPI_COMM_RANK(MPI_COMM_WORLD, my_id, ierr)
  else
    my_id = 0
  endif

  if (my_id .eq. 0) then
     write(*,*) "Generating Adaptive RTree Slices..."
  endif

  ! Loop over all elements
  !$omp parallel do default(shared) &
  !$omp private(i, j, phi_planes, num_planes, num_boxes, rmin, rmax, zmin, zmax) &
  !$omp private(rmin_next, rmax_next, zmin_next, zmax_next, phi_start, phi_end, idx) &
  !$omp schedule(dynamic)
  do i = 1, n
      call get_element_critical_planes(node_list, element_list, i, phi_planes)
      ! For testing, use uniform planes
      
      num_planes = size(phi_planes)
      num_boxes  = num_planes - 1
      
      ! Store size for C interface
      sizes_c(i) = int(num_boxes, C_INT)
      
      ! Allocate memory for this specific element (ND = 3: R, Z, Phi)
      allocate(element_data(i)%min_vals(num_boxes * ND))
      allocate(element_data(i)%max_vals(num_boxes * ND))

      ! -----------------------------------------------------------------------
      ! STEP 2: Calculate Bounding Boxes between planes
      ! -----------------------------------------------------------------------
      
      ! Pre-calculate initial plane Min/Max R/Z
      call RZP_minmax(node_list, element_list, i, phi_planes(1), &
                                rmin, rmax, zmin, zmax)

      do j = 1, num_boxes
          phi_start = phi_planes(j)
          phi_end   = phi_planes(j+1)

          ! Calculate R/Z min/max at the NEXT plane
          call RZP_minmax(node_list, element_list, i, phi_planes(j+1), &
                                    rmin_next, rmax_next, zmin_next, zmax_next)

          idx = (j - 1) * ND 

          ! Min Coordinate (R, Z, Phi)
          ! Note: We take the min of the R/Z values at both ends of the slice
          ! (Assuming linear or monotonic interpolation between critical planes)
          element_data(i)%min_vals(idx + 1) = real(min(rmin, rmin_next), C_DOUBLE)
          element_data(i)%min_vals(idx + 2) = real(min(zmin, zmin_next), C_DOUBLE)
          element_data(i)%min_vals(idx + 3) = real(phi_start, C_DOUBLE)

          ! Max Coordinate (R, Z, Phi)
          element_data(i)%max_vals(idx + 1) = real(max(rmax, rmax_next), C_DOUBLE)
          element_data(i)%max_vals(idx + 2) = real(max(zmax, zmax_next), C_DOUBLE)
          element_data(i)%max_vals(idx + 3) = real(phi_end, C_DOUBLE)

          ! Move next to current for the next iteration
          rmin = rmin_next; rmax = rmax_next
          zmin = zmin_next; zmax = zmax_next
      end do

      ! Clean up the temporary planes array for this element
      if (allocated(phi_planes)) deallocate(phi_planes)

      ! -----------------------------------------------------------------------
      ! STEP 3: Assign C Pointers
      ! -----------------------------------------------------------------------
      ! We get the address of the first element of the allocated arrays
      min_ptrs(i) = C_LOC(element_data(i)%min_vals)
      max_ptrs(i) = C_LOC(element_data(i)%max_vals)

  end do
  !$omp end parallel do
  
   !   do i = 1, n
   !    num_boxes = sizes_c(i)
   !       write(*,*) 'boxes for element ', i, ' : ', num_boxes
   !       do j = 1, num_boxes
   !          idx = (j - 1) * ND 
   !          write(*,*) 'Box ', j, ': Rmin=', element_data(i)%min_vals(idx + 1), &
   !                      ' Rmax=', element_data(i)%max_vals(idx + 1), &
   !                      ' Zmin=', element_data(i)%min_vals(idx + 2), &
   !                      ' Zmax=', element_data(i)%max_vals(idx + 2), &
   !                      ' Phimin=', element_data(i)%min_vals(idx + 3), &
   !                      ' Phimax=', element_data(i)%max_vals(idx + 3)
   !       end do
   !    end do

  if (my_id .eq. 0) then
    write(*,*) "Populating C++ RTree..."
  endif
  
  ! Call the C function
  call element_rtree(n_elms_c, sizes_c, min_ptrs, max_ptrs)

  ! Cleanup (Assuming C++ copied the data, which ElementTree.Insert usually does)
  rtree_initialized = .true.

end subroutine populate_element_rtree_3D

!> Find probable neighbours of element i and return their indices.
!> This is done by taking the bounding box of an element and expanding
!> it slightly (10^-6 in absolute value) and returning all elements in
!> this box, except element i
subroutine nearby_elements(node_list, element_list, i_elm, i_nearby)
  use constants, only: PI
  type(type_node_list), intent(in)    :: node_list
  type(type_element_list), intent(in) :: element_list
  integer, intent(in)                 :: i_elm
  integer, dimension(:), allocatable, intent(out) :: i_nearby
  integer(C_int), dimension(:), allocatable       :: i_nearby_C

  real*8, dimension(n_vertex_max, 2) :: vertices
  real(C_DOUBLE), dimension(ND) :: min_bb, max_bb
  integer(C_INT) :: num_elements
  integer :: iv

  do iv=1,n_vertex_max
     vertices(iv,:) = get_vertex_pos_in_rtree_plane(node_list%node(element_list%element(i_elm)%vertex(iv))%x(1:n_coord_tor,1,1:2))
  enddo

  min_bb (1) = real(minval(vertices(:,1)) - 1d-6, kind=C_DOUBLE)
  max_bb (1) = real(maxval(vertices(:,1)) + 1d-6, kind=C_DOUBLE)

  min_bb (2) = real(minval(vertices(:,2)) - 1d-6, kind=C_DOUBLE)
  max_bb (2) = real(maxval(vertices(:,2)) + 1d-6, kind=C_DOUBLE)

#if STELLARATOR_MODEL
  min_bb (3) = real(0.d0, kind=C_DOUBLE)
  max_bb (3) = real(0.d0, kind=C_DOUBLE)
#endif

  num_elements = int(num_elements_in_rect(min_bb, max_bb))
  allocate(i_nearby(num_elements),i_nearby_C(num_elements))
  num_elements = int(elements_in_rect(min_bb, max_bb, i_nearby_C))
  i_nearby = int(i_nearby_C(1:num_elements))
end subroutine nearby_elements

!> Find elements that could probably contain this point.
subroutine elements_containing_point(R, Z, phi, i_elms)
  real*8, intent(in)                              :: R, Z, phi
  integer, dimension(:), allocatable, intent(out) :: i_elms

  real(C_DOUBLE), dimension(ND) :: min_bb, max_bb
  integer(C_INT) :: num_elements
  integer(C_int), dimension(:), allocatable :: i_nearby_C

  if (R .ne. R .or. Z .ne. Z .or. phi .ne. phi) then
    write(*,*) "Warning: NaN supplied for R or Z in elements_containing_point, returning 0 elements"
    allocate(i_elms(0))
    return
  end if

  if (.not. rtree_initialized) then
    write(*,*) "Warning: RTree not initialised, exiting"
    stop 11
  end if

  min_bb (1) = real(R, kind=C_DOUBLE)
  max_bb (1) = real(R, kind=C_DOUBLE)

  min_bb (2) = real(Z, kind=C_DOUBLE)
  max_bb (2) = real(Z, kind=C_DOUBLE)

#if STELLARATOR_MODEL
  min_bb (3) = real(phi, kind=C_DOUBLE)
  max_bb (3) = real(phi, kind=C_DOUBLE)
#endif

  ! This calls the search routine twice! To get around that either pass a large enough array alway
  ! or implement it as a mask/bitfield of the total number of elements.
  num_elements = int(num_elements_in_rect(min_bb, max_bb))
  allocate(i_nearby_C(num_elements), i_elms(num_elements))
  num_elements = int(elements_in_rect(min_bb, max_bb, i_nearby_C))
  i_elms = int(i_nearby_C(1:num_elements))
end subroutine elements_containing_point


pure function get_vertex_pos_in_rtree_plane(x) result(pos)
  use phys_module, only: i_plane_rtree
  use basis_at_gaussian, only: HZ_coord
  implicit none

  real*8, intent(in)    :: x(n_coord_tor, 2)
  real*8   :: pos(2)
  integer  :: i

  do i=1,2
    pos(i) = sum(x(1:n_coord_tor,i)*HZ_coord(1:n_coord_tor, i_plane_rtree))
  end do
end function get_vertex_pos_in_rtree_plane

!> calculates the basis functions at the Gaussian points
pure function get_arbitrary_HZ_coord(tree_slices_requested) result(arbitrary_HZ_coord)
  use constants, only: PI
  use phys_module, only: n_coord_tor, n_coord_period, mode_coord

  implicit none

  integer, intent(in) :: tree_slices_requested

  ! --- local variables
  integer :: k,i
  real*8  :: phi
  real*8 :: arbitrary_HZ_coord (n_coord_tor,tree_slices_requested)

  do k=1,tree_slices_requested

    phi = 2.d0*PI*float(k-1)/float(tree_slices_requested) / float(n_coord_period)

    arbitrary_HZ_coord(1,k)   = 1.d0

    do i=1,(n_coord_tor-1)/2
      arbitrary_HZ_coord(2*i,k)      = + cos(mode_coord(2*i)  *phi)
      arbitrary_HZ_coord(2*i+1,k)    = - sin(mode_coord(2*i+1)*phi)
    enddo
  enddo
end function get_arbitrary_HZ_coord


end module mod_element_rtree

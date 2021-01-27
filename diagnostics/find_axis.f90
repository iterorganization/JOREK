!> Locate the position of the magnetic axis.
subroutine find_axis(my_id, node_list, element_list, psi_axis, R_axis, Z_axis, i_elm_axis, s_axis, &
  t_axis, ifail)

use data_structure
use gauss
use basis_at_gaussian
use equil_info,  only: ES
use phys_module, only: R_geo, Z_geo, axis_srch_radius, R_axis_t, Z_axis_t, index_start  
use mod_interp

implicit none

interface
   subroutine mnewtax(node_list,element_list,i_elm, r, s, errx, errf, ifail)
     use data_structure
     type (type_node_list)    :: node_list
     type (type_element_list) :: element_list
     real*8    :: r, s, errf, errx
     integer   :: ifail, i_elm
   end subroutine mnewtax
end interface

! --- Routine parameters
integer,                 intent(in)  :: my_id        !< MPI proc number
type(type_node_list),    intent(in)  :: node_list    !< List of grid nodes
type(type_element_list), intent(in)  :: element_list !< List of grid elements
real*8,                  intent(out) :: psi_axis     !< Poloidal flux at axis
real*8,                  intent(out) :: R_axis       !< R-position of axis
real*8,                  intent(out) :: Z_axis       !< Z-position of axis
real*8,                  intent(out) :: s_axis       !< s-position of axis in Bezier element i_elm_axis
real*8,                  intent(out) :: t_axis       !< t-position of axis in Bezier element i_elm_axis
integer,                 intent(out) :: i_elm_axis   !< Bezier element, axis is located in
integer,                 intent(out) :: ifail        !< Error code

! --- Local variables
real*8  :: ps_x, ps_y, ps_s, ps_t, xjac
real*8  :: R, R_s, R_t, R_st, R_ss, R_tt, Z, Z_s, Z_t, Z_st, Z_ss, Z_tt, P, P_s, P_t, P_st, P_ss, P_tt
integer :: ij_axis(2), i, iv, ms, mt, kf, kv, i_tries, n_tries, i_elm_axis_init, min_indices(3)
real*8  :: x(2), s, t, xerr, ferr, s_axis_init, t_axis_init
real*8  :: R0, Z0, search_radius
logical :: found_axis, axis_in_rst_file
real*8,  allocatable :: grad_psi(:,:,:)
logical, allocatable :: include_pt(:,:,:)

if (my_id .eq. 0) then
  write(*,*) '*********************************'
  write(*,*) '*     find_axis                 *'
  write(*,*) '*********************************'
endif

n_tries = 500          ! --- number of attempts to find the axis 
found_axis = .false.

allocate(grad_psi(element_list%n_elements,4,4))            ! --- vector storing |grad_psi| at gaussian poitns
allocate(include_pt(element_list%n_elements,4,4))          ! --- vector storing if point should be considered or not
grad_psi    = 0.d0
include_pt  = .false.

ifail        = 1
ij_axis      = 1 
psi_axis     = 1.d20

! --- define geometrical limits to search for the axis
search_radius = axis_srch_radius
if (axis_srch_radius > 50.d0) search_radius = 0.2d0 * R_geo

! --- if restarting, check that the previous axis is in the restart file 
axis_in_rst_file = .false.
if( (index_start /= 0) .and. allocated(R_axis_t)) then
  if (R_axis_t(index_start) /= 0.d0) axis_in_rst_file = .true.  
endif

if (ES%initialized) then
  R0 = ES%R_axis
  Z0 = ES%Z_axis
else if ((.not. ES%initialized) .and. (index_start > 0) .and. axis_in_rst_file) then
  R0 = R_axis_t(index_start)
  Z0 = Z_axis_t(index_start)
else 
  R0 = R_geo
  Z0 = Z_geo
endif


! save |grad_psi| at gaussian points of all elements
do i=1,element_list%n_elements   ! --- loop over elements

  do ms = 1, 4           ! 4 Gaussian points
    do mt = 1, 4         ! 4 Gaussian points

      ps_s = 0.d0
      ps_t = 0.d0
      R_s  = 0.d0 
      Z_s  = 0.d0
      R_t  = 0.d0 
      Z_t  = 0.d0
      R    = 0.d0
      Z    = 0.d0

      do kf = 1, 4       ! 4 basis functions
        do kv = 1, 4     ! 4 vertices

          iv = element_list%element(i)%vertex(kv)

          ps_s = ps_s + node_list%node(iv)%values(1,kf,1) * element_list%element(i)%size(kv,kf) * H_s(kv,kf,ms,mt)
          ps_t = ps_t + node_list%node(iv)%values(1,kf,1) * element_list%element(i)%size(kv,kf) * H_t(kv,kf,ms,mt)

          R   = R   + node_list%node(iv)%x(1,kf,1) * element_list%element(i)%size(kv,kf) * H(kv,kf,ms,mt)
          Z   = Z   + node_list%node(iv)%x(1,kf,2) * element_list%element(i)%size(kv,kf) * H(kv,kf,ms,mt)

          R_s = R_s + node_list%node(iv)%x(1,kf,1) * element_list%element(i)%size(kv,kf) * H_s(kv,kf,ms,mt)
          Z_s = Z_s + node_list%node(iv)%x(1,kf,2) * element_list%element(i)%size(kv,kf) * H_s(kv,kf,ms,mt)
          R_t = R_t + node_list%node(iv)%x(1,kf,1) * element_list%element(i)%size(kv,kf) * H_t(kv,kf,ms,mt)
          Z_t = Z_t + node_list%node(iv)%x(1,kf,2) * element_list%element(i)%size(kv,kf) * H_t(kv,kf,ms,mt)

        enddo
      enddo

      xjac = R_s * Z_t - R_t * Z_s
      ps_x = (  ps_s * Z_t - ps_t * Z_s)/ xjac
      ps_y = (- ps_s * R_t + ps_t * R_s)/ xjac

      grad_psi(i,ms,mt) = sqrt(ps_x*ps_x + ps_y*ps_y)
      
      ! --- include points if they are within the search region
      include_pt(i,ms,mt) = sqrt( (R-R0)**2.d0 + (Z-Z0)**2.d0 ) < search_radius
          
    enddo
  enddo

enddo   ! --- end loop over elements


do i_tries=1,  n_tries  ! --- start attempts to find the axis

  ! --- min_indices = indices for gaussian point with min |grad_psi|,   (1) = element index, (2) = s-gaussian point index, (3) = t-gaussian point index
  min_indices(:) = minloc(grad_psi, mask=include_pt)  ! --- find minimum of |grad_psi|
  if (.not. any(include_pt)) min_indices = 0
  
  if ((min_indices(1) == 0) .and. (i_tries == 1)) then     ! --- if all elements are initially excluded, stop search and initialize values
    found_axis = .false.
    s_axis_init     = 0.d0
    t_axis_init     = 0.d0
    i_elm_axis_init = 1
    exit
  else if  (min_indices(1) == 0) then   ! --- if all elements have been excluded, exit search
    found_axis = .false.
    exit
  endif
  
  i_elm_axis = min_indices(1)    ! --- element with minimum |grad_psi|
  s          = Xgauss(min_indices(2)) 
  t          = Xgauss(min_indices(3))

  ! --- Find \grad_psi = 0 in i_elm_axis with Newton's method
  call mnewtax(node_list,element_list,i_elm_axis,s,t,xerr,ferr,ifail)

  if (ifail .ne. 0 ) then      ! --- if Newton's method failed, exclude element in next search
    include_pt(i_elm_axis,:,:) = .false.
  else
    found_axis = .true.
    s_axis     = s
    t_axis     = t
    exit
  endif
  
  if (i_tries == 1) then    ! --- save first attempt in case all the attempts fail
    s_axis_init     = s
    t_axis_init     = t
    i_elm_axis_init = i_elm_axis
  endif
  
enddo !--- end tries

if (.not. found_axis) then    ! --- if all the attempts to find axis failed, the axis is the initial solution
  s_axis     = s_axis_init
  t_axis     = t_axis_init
  i_elm_axis = i_elm_axis_init
endif

call interp(node_list,element_list,i_elm_axis,1,1,s_axis,t_axis,psi_axis,P_s,P_t,P_st,P_ss,P_tt)

call interp_RZ(node_list,element_list,i_elm_axis,s_axis,t_axis,R_axis,Z_axis)

if ((.not. found_axis) .and. (my_id .eq. 0) .and. (ES%initialized)) then
  write(*,*) '  WARNING: axis was not properly found after ', n_tries, ' attempts'
  write(*,*) '  Please check that you have set axis_srch_radius correctly in the input file.'
  write(*,*) '  If the expected magnetic axis lies outside the circle with radius'
  write(*,*) '  axis_srch_radius and center(R_geo,Z_geo) or previous axis (if defined),'
  write(*,*) '  you must increase axis_srch_radius. Now you are using:'
  write(*,*) '    axis_srch_radius = ',  search_radius
  write(*,*) '    R_geo, Z_geo     = ',  R_geo, Z_geo
endif
if (my_id .eq. 0) write(*,'(A,i6,4f14.8)') ' magnetic axis : ',i_elm_axis,R_axis,Z_axis,psi_axis

deallocate(include_pt, grad_psi)

return
end subroutine find_axis

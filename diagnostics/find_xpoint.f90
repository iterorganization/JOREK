!> Routine determines the position(s) of the xpoint(s).
subroutine find_xpoint(my_id,node_list,element_list,psi_xpoint,R_xpoint,Z_xpoint,i_elm_xpoint,s_xpoint,t_xpoint,xcase,ifail)

use data_structure
use gauss
use basis_at_gaussian
use phys_module, only: tokamak_device
use mod_interp

implicit none

! --- Routine parameters
integer,                  intent(in)    :: my_id
type (type_node_list),    intent(in)    :: node_list
type (type_element_list), intent(in)    :: element_list
real*8,                   intent(out)   :: psi_xpoint(2)
real*8,                   intent(out)   :: R_xpoint(2)
real*8,                   intent(out)   :: Z_xpoint(2)
integer,                  intent(out)   :: i_elm_xpoint(2)
real*8,                   intent(out)   :: s_xpoint(2)
real*8,                   intent(out)   :: t_xpoint(2)
integer,                  intent(in)    :: xcase
integer,                  intent(out)   :: ifail

! --- Local variables
real*8  :: ps_s, ps_t, ps_x, ps_y, xjac
real*8  :: R, R_s, R_t, Z, Z_s, Z_t, P, P_s, P_t, P_st, P_ss, P_tt
real*8  :: x(2), s, t, xerr, ferr, s_xp_init(2), t_xp_init(2)
integer :: ij_xpoint(2,2), i, iv, ms, mt, kf, kv, i_tries, n_tries
integer :: i_elm_xp_init(2), min_indices_lw(3), min_indices_up(3)
logical :: found_upper, found_lower
real*8,  allocatable :: grad_psi(:,:,:)
logical, allocatable :: include_pt_lw(:,:,:), include_pt_up(:,:,:)

if (my_id .eq. 0) then
  write(*,*) '*********************************'
  write(*,*) '*     find_xpoint               *'
  write(*,*) '*********************************'
endif

ifail   = 1
n_tries = 500

psi_xpoint = 0.
R_xpoint   = 0.;    Z_xpoint = 0.
s_xpoint   = 0.;    t_xpoint = 0.
i_elm_xpoint = 0

allocate(grad_psi      (element_list%n_elements,4,4))            ! --- vector storing |grad_psi| at gaussian poitns
allocate(include_pt_lw (element_list%n_elements,4,4))
allocate(include_pt_up (element_list%n_elements,4,4))
grad_psi    = 0.d0
include_pt_lw = .false.
include_pt_up = .false.

found_upper = .false. 
found_lower = .false.


do i=1,element_list%n_elements    ! --- loop over elements
  
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
    
      ! --- Look for the lower Xpoint
      if (xcase .ne. 2) then
        if (     ((tokamak_device(1:4) .ne. 'MAST') .and. (tokamak_device(1:7) .ne. 'COMPASS') .and. (Z .lt. -0.4d0)) &
            .or. ((tokamak_device(1:4) .eq. 'MAST') .and. (Z .lt. -0.4d0) .and. (R .gt. 0.45d0) .and. (R .lt. 1.d0))  &
            .or. ((tokamak_device(1:7) .eq. 'COMPASS') .and. (Z .lt. -0.2d0))) then
          include_pt_lw(i,ms,mt) = .true.        
        endif
      endif
      
      ! --- And for the upper Xpoint
      if (xcase .ne. 1) then
        if (     ((tokamak_device(1:4) .ne. 'MAST') .and. (Z .gt.  0.4d0)) &
            .or. ((tokamak_device(1:4) .eq. 'MAST') .and. (Z .gt.  0.4d0) .and. (R .gt. 0.45d0) .and. (R .lt. 1.d0)) ) then
          include_pt_up(i,ms,mt) = .true.
        endif
      endif

    enddo
  enddo

enddo    ! --- end loop over elements


if(xcase .ne. 2) then
  do i_tries=1,  n_tries  ! --- start attempts to find the lower x-point
    
    ! --- min_indices = indices for gaussian point with min |grad_psi|,   (1) = element index, (2) = s-gaussian point index, (3) = t-gaussian point index
    min_indices_lw(:) = minloc(grad_psi, mask=include_pt_lw)
    if (.not. any(include_pt_lw)) min_indices_lw = 0

    if ((min_indices_lw(1) == 0) .and. (i_tries == 1)) then     ! --- if all elements are initially excluded, stop search and initialize values
      found_lower      = .false.
      s_xp_init(1)     = 0.d0
      t_xp_init(1)     = 0.d0
      i_elm_xp_init(1) = 1
      exit
    else if  (min_indices_lw(1) == 0) then   ! --- if all elements have been excluded, exit search
      found_lower = .false.
      exit
    endif
    
    i_elm_xpoint(1) = min_indices_lw(1)    ! --- element with minimum |grad_psi|
    s = Xgauss(min_indices_lw(2)) 
    t = Xgauss(min_indices_lw(3))
    
    call mnewtax(node_list,element_list,i_elm_xpoint(1),s,t,xerr,ferr,ifail)
    if (ifail .ne. 0 ) then      ! --- if Newton's method failed, exclude element in next search
      include_pt_lw(i_elm_xpoint(1),:,:) = .false.
    else
      found_lower   = .true.
      s_xpoint(1)   = s
      t_xpoint(1)   = t
      exit
    endif
    if (i_tries == 1) then    ! --- save first attempt in case all the attempts fail
      s_xp_init(1)     = s
      t_xp_init(1)     = t
      i_elm_xp_init(1) = i_elm_xpoint(1)
    endif     
  enddo
endif

if(xcase .ne. 1) then

  do i_tries=1,  n_tries  ! --- start attempts to find the upper x-point

    ! --- min_indices = indices for gaussian point with min |grad_psi|,   (1) = element index, (2) = s-gaussian point index, (3) = t-gaussian point index
    min_indices_up(:) = minloc(grad_psi, mask=include_pt_up)
    if (.not. any(include_pt_up)) min_indices_up = 0
    
    if ((min_indices_up(1) == 0) .and. (i_tries == 1)) then     ! --- if all elements are initially excluded, stop search and initialize values
      found_upper      = .false.
      s_xp_init(2)     = 0.d0                             
      t_xp_init(2)     = 0.d0
      i_elm_xp_init(2) = 1
      exit
    else if  (min_indices_up(1) == 0) then   ! --- if all elements have been excluded, exit search
      found_upper     = .false.
      exit
    endif

    i_elm_xpoint(2) = min_indices_up(1)    ! --- element with minimum |grad_psi|
    s = Xgauss(min_indices_up(2)) 
    t = Xgauss(min_indices_up(3))
    
    call mnewtax(node_list,element_list,i_elm_xpoint(2),s,t,xerr,ferr,ifail)
    if (ifail .ne. 0 ) then       ! --- if Newton's method failed, exclude element in next search
      include_pt_up(i_elm_xpoint(2),:,:) = .false.
    else
      found_upper   = .true.
      s_xpoint(2)   = s
      t_xpoint(2)   = t
      exit
    endif 
    if (i_tries == 1) then    ! --- save first attempt in case all the attempts fail
      s_xp_init(2)     = s
      t_xp_init(2)     = t
      i_elm_xp_init(2) = i_elm_xpoint(2)
    endif
  enddo ! --- end attempts     

endif  
  


if(xcase .ne. 2) then
  if (.not. found_lower) then    ! --- if all the attempts failed, take the initial solution
    s_xpoint(1)     = s_xp_init(1)     
    t_xpoint(1)     = t_xp_init(1)     
    i_elm_xpoint(1) = i_elm_xp_init(1) 
  endif

  call interp(node_list,element_list,i_elm_xpoint(1),1,1,s_xpoint(1),t_xpoint(1),psi_xpoint(1),P_s,P_t,P_st,P_ss,P_tt)
  call interp_RZ(node_list,element_list,i_elm_xpoint(1),s_xpoint(1),t_xpoint(1),R_xpoint(1),R_s,R_t,Z_xpoint(1),Z_s,Z_t)

  xjac = R_s * Z_t - R_t * Z_s
  ps_x = (  P_s * Z_t - P_t * Z_s)/ xjac
  ps_y = (- P_s * R_t + P_t * R_s)/ xjac

  if (my_id .eq. 0) then
    write(*,'(A,i6,4f14.8)') ' Lower X-point : ',i_elm_xpoint(1),R_xpoint(1),Z_xpoint(1),psi_xpoint(1),sqrt(ps_x**2+ps_y**2)
  endif
  if ((.not. found_lower )) write(*,*) 'WARNING: lower X-point not properly found after ', n_tries, ' attempts'
endif

if(xcase .ne. 1) then 
  if (.not. found_upper) then    ! --- if all the attempts failed, take the initial solution
    s_xpoint(2)     = s_xp_init(2)     
    t_xpoint(2)     = t_xp_init(2)     
    i_elm_xpoint(2) = i_elm_xp_init(2) 
  endif
  
  call interp(node_list,element_list,i_elm_xpoint(2),1,1,s_xpoint(2),t_xpoint(2),psi_xpoint(2),P_s,P_t,P_st,P_ss,P_tt)
  call interp_RZ(node_list,element_list,i_elm_xpoint(2),s_xpoint(2),t_xpoint(2),R_xpoint(2),R_s,R_t,Z_xpoint(2),Z_s,Z_t)

  xjac = R_s * Z_t - R_t * Z_s
  ps_x = (  P_s * Z_t - P_t * Z_s)/ xjac
  ps_y = (- P_s * R_t + P_t * R_s)/ xjac

  if (my_id .eq. 0) then
    write(*,'(A,i6,4f14.8)') ' Upper X-point : ',i_elm_xpoint(2),R_xpoint(2),Z_xpoint(2),psi_xpoint(2),sqrt(ps_x**2+ps_y**2)
  endif  
  if ((.not. found_upper )) write(*,*) 'WARNING: upper X-point not properly found after ', n_tries, ' attempts'
endif

deallocate(include_pt_lw,include_pt_up, grad_psi)

return
end subroutine find_xpoint

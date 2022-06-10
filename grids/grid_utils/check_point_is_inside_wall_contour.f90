module check_point_is_inside_wall_contour

contains 

!> This routine checks whether a point is outside(=0), inside(=1)
subroutine check_point_is_inside_wall(R, Z, inside)

  use phys_module
  implicit none
  
  ! --- Routine parameters
  real*8,  intent(in)           :: R, Z
  integer, intent(inout)        :: inside
  
  call check_point_is_inside_contour(R, Z, n_limiter, R_limiter, Z_limiter, inside)
  
  return

end subroutine check_point_is_inside_wall
  
  

!> Same thing but for an arbitrary contour --- outside(=0), inside(=1)
! algorithm is originally from W. Randolph Franklin (I think)
subroutine check_point_is_inside_contour(R, Z, n_ctr, R_ctr, Z_ctr, inside)

  use phys_module
  implicit none
  
  ! --- Routine parameters
  real*8,  intent(in)           :: R, Z
  integer, intent(inout)        :: n_ctr
  real*8,  intent(in)           :: R_ctr(n_ctr), Z_ctr(n_ctr)
  integer, intent(inout)        :: inside
  
  ! --- Local variables
  integer       :: i, count
  real*8        :: R_int
  real*8        :: R_tmp1, Z_tmp1
  real*8        :: R_tmp2, Z_tmp2
  
  R_tmp1 = R_ctr(n_ctr)
  Z_tmp1 = Z_ctr(n_ctr)
  
  count = 0
  do i=1,n_ctr
    
    R_tmp2 = R_ctr(i)
    Z_tmp2 = Z_ctr(i)
    
    if ( (Z .ge. min(Z_tmp1,Z_tmp2)) .and. (Z .le. max(Z_tmp1,Z_tmp2)) .and. (R .lt. max(R_tmp1, R_tmp2)) ) then
      if (Z_tmp1 .ne. Z_tmp2) then
        R_int = (Z-Z_tmp1) * (R_tmp2-R_tmp1) / (Z_tmp2-Z_tmp1) + R_tmp1
        if ( (R .lt. R_int) .and. (Z .ne. Z_tmp2) ) count = count + 1
      endif
    endif
    
    R_tmp1 = R_tmp2
    Z_tmp1 = Z_tmp2
    
  enddo
  
  if (mod(count,2) .eq. 0) then
    inside = 0
  else
    inside = 1
  endif
  
  
  return

end subroutine check_point_is_inside_contour
  
  

end module check_point_is_inside_wall_contour

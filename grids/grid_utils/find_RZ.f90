subroutine find_RZ(node_list,element_list,R_find,Z_find,R_out,Z_out,ielm_out,s_out,t_out,ifail)
!-------------------------------------------------------------------------
!< Find all elements for which minmax is correct and run find_RZ_single on those.
!< Return the first result.
!-------------------------------------------------------------------------
use data_structure
use mod_element_rtree, only: elements_containing_point
implicit none

type (type_node_list), intent(in)    :: node_list
type (type_element_list), intent(in) :: element_list
real*8, intent(in)     :: R_find, Z_find
real*8, intent(out)    :: R_out,Z_out,s_out,t_out
integer, intent(inout) :: ielm_out
integer, intent(out)   :: ifail

integer :: k
integer, dimension(:), allocatable :: i_elms

ielm_out = 0
call elements_containing_point(R_find, Z_find, i_elms)

! then loop through all
do k=1,size(i_elms)
  call find_RZ_single(node_list,element_list,i_elms(k),R_find,Z_find,R_out,Z_out,ielm_out,s_out,t_out,ifail)
  if (ifail .eq. 0) exit
enddo

if (ielm_out .eq. 0) ifail = 99
if (ifail .eq. 999) ielm_out = 0 ! Otherwise testing ielm=0 on output does not
! work anymore (and we don't always check ifail)
end subroutine find_RZ

subroutine find_RZ_single(node_list,element_list,i_elm,R_find,Z_find,R_out,Z_out,ielm_out,s_out,t_out,ifail)
!-------------------------------------------------------------------------
!< solves two non-linear equations using Newtons method (from numerical recipes)
!< LU decomposition replaced by explicit solution of 2x2 matrix.
!<
!< finds the crossing of two coordinate lines given as a series of cubics in element
!< i_elm
!-------------------------------------------------------------------------
use data_structure
use mod_interp, only: interp_RZ
implicit none

type (type_node_list), intent(in)    :: node_list
type (type_element_list), intent(in) :: element_list
integer, intent(in)    :: i_elm
real*8, intent(in)     :: R_find, Z_find
real*8, intent(out)    :: R_out,Z_out,s_out,t_out
integer, intent(out)   :: ielm_out
integer, intent(out)   :: ifail

integer :: i, ntrial, istart
real*8  :: RRg1,dRRg1_dr,dRRg1_ds
real*8  :: ZZg1,dZZg1_dr,dZZg1_ds
real*8  :: tolx, tolf, errx, errf, temp, dis
real*8  :: x(2), FVEC(2), FJAC(2,2), p(2)

ntrial = 20
tolx = 1.d-8
tolf = 1.d-15

ielm_out = i_elm ! Since we only test a single element

do istart = 1,5

  if (istart .eq. 1) then
    x(1) = 0.5d0
    x(2) = 0.5d0
  elseif (istart .eq. 2) then
    x(1) = 0.75d0
    x(2) = 0.75d0
  elseif (istart .eq. 3) then
    x(1) = 0.75d0
    x(2) = 0.25d0
  elseif (istart .eq. 4) then
    x(1) = 0.25d0
    x(2) = 0.75d0
  elseif (istart .eq. 5) then
    x(1) = 0.25d0
    x(2) = 0.25d0
  endif

  ifail = 999

  do i=1,ntrial

    call interp_RZ(node_list,element_list,i_elm,x(1),x(2),RRg1,dRRg1_dr,dRRg1_ds, &
                                                    ZZg1,dZZg1_dr,dZZg1_ds)

    FVEC(1)   = RRg1 - R_find
    FVEC(2)   = ZZg1 - Z_find
    FJAC(1,1) = dRRg1_dr
    FJAC(1,2) = dRRg1_ds
    FJAC(2,1) = dZZg1_dr
    FJAC(2,2) = dZZg1_ds

    errf=abs(fvec(1))+abs(fvec(2))

!      write(*,'(A,i3,8e16.8)') ' newton   : ',i,errf,errx,x,RRg1,R_find,ZZg1,Z_find
!      write(*,'(A,i3,8e16.8)') ' newton   : ',i,dRRg1_dr,dRRg1_ds,dZZg1_dr,dZZg1_ds

    if (errf .le. tolf) then

      s_out     = x(1)
      t_out     = x(2)

      ielm_out  = i_elm
      R_out     = RRg1
      Z_out     = ZZg1

!        write(*,'(A,i3,4e16.8)') ' newton (1) : ',i,errf,errx,x

      ifail = 0
      return
    endif

    p = -fvec

    temp = p(1)
    dis  = fjac(2,2)*fjac(1,1)-fjac(1,2)*fjac(2,1)

    if (dis .ne. 0.d0) then
      p(1) = (fjac(2,2)*p(1)-fjac(1,2)*p(2))/dis
      p(2) = (fjac(1,1)*p(2)-fjac(2,1)*temp)/dis
    else
      exit
    endif

    errx=abs(p(1)) + abs(p(2))

    p = min(p,+0.25d0)
    p = max(p,-0.25d0)

    x = x + p

    x = max(x,+0.d0)
    x = min(x,+1.d0)

    if (errx .le. tolx) then

      s_out     = x(1)
      t_out     = x(2)

      ielm_out  = i_elm
      R_out     = RRg1
      Z_out     = ZZg1

!        write(*,'(A,i3,4e16.8)') ' newton (2) : ',i,errf,errx,x

      ifail = 0
      return
    endif

  enddo
enddo
end subroutine find_RZ_single

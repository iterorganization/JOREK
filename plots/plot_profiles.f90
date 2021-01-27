subroutine plot_profiles(node_list,element_list,Rplot,Zplot)
!****************************************************************************
!* plots 1D profile of the n=0 equilibrium flow                             *
!****************************************************************************
use tr_module 
use data_structure
use phys_module
use basis_at_gaussian
use mod_interp, only: interp

implicit none
type (type_node_list)    :: node_list
type (type_element_list) :: element_list

real*8  :: Rplot(2), Zplot(2), Rp, Zp, Rmin, Rmax, Zmin, Zmax, s_out, t_out, R_out, Z_out
real*8  :: P, P_s, P_t, P_st, P_ss, P_tt, U, U_s, U_t, U_st, U_ss, U_tt
real*8  :: R, R_s, R_t, R_st, R_ss, R_tt, Z, Z_s, Z_t, Z_st, Z_ss, Z_tt
real*8  :: ps0_x, ps0_y, u0_x, u0_y, xjac, ss
real*8, allocatable :: xp(:,:),yp(:,:)
integer :: i, k, i_elm, i_var, i_harm, ifail, iplot, nplot

Rmin =  1.d20
Rmax = -1.d20
Zmin =  1.d20
Zmax = -1.d20
do i=1,node_list%n_nodes
  Rmin = min(Rmin,node_list%node(i)%x(1,1,1))
  Rmax = max(Rmax,node_list%node(i)%x(1,1,1))
  Zmin = min(Zmin,node_list%node(i)%x(1,1,2))
  Zmax = max(Zmax,node_list%node(i)%x(1,1,2))
enddo

nplot  = 1001

call tr_allocate(xp,1,nplot,1,2,"xp",CAT_GRID)
call tr_allocate(yp,1,nplot,1,n_var,"yp",CAT_GRID)

!---------------------------------- plot perpendicular equilibrium flow

xp = 0.d0
yp = 0.d0

iplot = 0
do i=1,nplot

  ss = float(i-1)/float(nplot-1)

  Rp = Rplot(1) + ss * (Rplot(2) - Rplot(1))
  Zp = Zplot(1) + ss * (Zplot(2) - Zplot(1))

  call find_RZ(node_list,element_list,Rp,Zp,R_out,Z_out,i_elm,s_out,t_out,ifail)

  if (ifail .eq. 0) then

    iplot = iplot+1
    xp(iplot,1:2) = (/ R_out, Z_out /)

    do k=1,n_var

      do i_harm = 1, n_tor

        call interp(node_list,element_list,i_elm,k,i_harm,s_out,t_out,P,P_s,P_t,P_st,P_ss,P_tt)

        yp(iplot,k)   = yp(iplot,k) + P * HZ(i_harm,1)

      enddo

    enddo

  endif

enddo

open(21,file='profiles.txt')
do i = 1,iplot
   write(21,'(20e14.6)') xp(i,1:2), yp(i,1:n_var),yp(i,5)*yp(i,6)*yp(i,7)
enddo


call tr_deallocate(xp,"xp",CAT_GRID)
call tr_deallocate(yp,"yp",CAT_GRID)

return
end

subroutine plot_velocity_profile(node_list,element_list,Rp_start,Zp_start,Rp_end,Zp_end)
!****************************************************************************
!* plots 1D profile of the n=0 equilibrium flow                             *
!****************************************************************************
use tr_module 
use data_structure
use phys_module
use mod_interp

implicit none
type (type_node_list)    :: node_list
type (type_element_list) :: element_list

real*8  :: Rp_start, Zp_start, Rp_end, Zp_end
real*8  :: Rp, Zp, Rmin, Rmax, Zmin, Zmax, s_out, t_out, R_out, Z_out
real*8  :: P, P_s, P_t, P_st, P_ss, P_tt, U, U_s, U_t, U_st, U_ss, U_tt
real*8  :: R, R_s, R_t, Z, Z_s, Z_t
real*8  :: ps0_x, ps0_y, u0_x, u0_y, xjac
real*8, allocatable :: xp(:),yp(:),stmp(:),stmp2(:)
integer :: i, i_elm, i_var, i_harm, ifail, iplot, nplot, nplot2

write(*,'(A)') '***************************************'
write(*,'(A)') '*         velocity profiles           *'
write(*,'(A)') '***************************************'

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

i_var  = 2
i_harm = 1

call tr_allocate(xp,1,nplot,"xp",CAT_GRID)
call tr_allocate(yp,1,nplot,"yp",CAT_GRID)
call tr_allocate(stmp,1,nplot,"stmp",CAT_GRID)

xr1  = 9999.
sig1 = 9999.
xr2  = 99999.
sig2 = 99999.

call meshac2(nplot,stmp,xr1,xr2,sig1,sig2,0.6d0,1.0d0)  ! allows concentration of gridpoints


!---------------------------------- plot perpendicular equilibrium flow

iplot = 0

do i=1,nplot

  Rp = Rp_start + stmp(i) * (Rp_end - Rp_start)
  Zp = Zp_start + stmp(i) * (Zp_end - Zp_start)
  
  call find_RZ(node_list,element_list,Rp,Zp,R_out,Z_out,i_elm,s_out,t_out,ifail)

  if (ifail .eq. 0) then

    call interp(node_list,element_list,i_elm,1,i_harm,s_out,t_out,P,P_s,P_t,P_st,P_ss,P_tt)
    call interp(node_list,element_list,i_elm,i_var,i_harm,s_out,t_out,U,U_s,U_t,U_st,U_ss,U_tt)
    call interp_RZ(node_list,element_list,i_elm,s_out,t_out,R,R_s,R_t,Z,Z_s,Z_t)

    xjac  = R_s * Z_t - R_t * Z_s
    ps0_x = (   Z_t * P_s - Z_s * P_t ) / xjac
    ps0_y = ( - R_t * P_s + R_s * P_t ) / xjac
    u0_x  = (   Z_t * U_s - Z_s * U_t ) / xjac
    u0_y  = ( - R_t * U_s + R_s * U_t ) / xjac

    if ((xjac .lt. 1.d-6) .or. (abs(ps0_x*ps0_x+ps0_y*ps0_y) .lt. 1.d-6)) ifail = 999

    if (ifail .eq. 0) then
      
      iplot = iplot+1
      
      xp(iplot) = sqrt( (R_out-Rp_start)**2 + (Z_out-Zp_start)**2 )
      
      yp(iplot) = R * sqrt(u0_x * u0_x + u0_y * u0_y)
    
    endif

  endif

enddo

open(21,file='profilesflow')
write(21,'(A,4f10.5)') '     distance,      poloidal_velocity at :',Rp_start,Zp_start,Rp_end,Zp_end
do i = 1,iplot
   write(21,*) xp(i), yp(i)
enddo

if ( write_ps ) call lplot(1,1,1,xp,yp,iplot,1,'perp. flow',10,'distance',8,'flow',4)

call tr_deallocate(xp,"xp",CAT_GRID)
call tr_deallocate(yp,"yp",CAT_GRID)
call tr_deallocate(xp,"stmp",CAT_GRID)

return
end



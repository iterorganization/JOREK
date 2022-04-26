!> The routine finds fluxsurfaces by finding crossings with the edges of the elements
subroutine find_flux_surfaces(my_id,xpoint,xcase,node_list,element_list,surface_list)

use constants
use tr_module 
use data_structure
use grid_xpoint_data
use mod_interp
use phys_module, only:   SDN_threshold
use mod_newton_methods
use equil_info

implicit none

! --- Routine parameters
integer,                  intent(in)     :: my_id        !< MPI proc number
logical,                  intent(in)     :: xpoint
integer,                  intent(in)     :: xcase
type (type_node_list)   , intent(in)     :: node_list
type (type_element_list), intent(in)     :: element_list
type (type_surface_list), intent(inout)  :: surface_list

! --- Local variables
real*8  :: psimin, psimax, a0, a1, a2, a3
real*8  :: dpsi_dr(4),dpsi_ds(4)
real*8  :: p1, dp1, dp4, p4, p2, p3, r_psi(4), s_psi(4), tht(4)
real*8  :: s, s2, s3, r_tmp, s_tmp, psr_tmp, pss_tmp, ttmp, tt
real*8  :: psi_xpoint(2), r_av, s_av

real*8  :: RRg(4), ZZg(4)
real*8  :: distance, distance_max
integer :: kp1, k_keep
real*8  :: Rmin, Rmax, Zmin, Zmax
integer :: k_topleft, k_topright
integer :: k_botleft, k_botright
real*8  :: dpsi_dr_copy(4),dpsi_ds_copy(4)
real*8  :: r_psi_copy(4), s_psi_copy(4), tht_copy(4)
integer :: l, i_neigh, Xneigh, icount
integer :: i, j, k, ifound, iv, im, is, n1, n2, n3
integer :: ifail, itht(4), itmp
integer :: n_found
real*8  :: st_found(n_order)

if (my_id == 0) then
  write(*,*) '***********************************'
  write(*,*) '*   find_flux_surfaces            *'
  write(*,*) '***********************************'
endif
!write(*,*) ' n_psi : ',surface_list%n_psi
!write(*,*) ' values : ',surface_list%psi_values(1),surface_list%psi_values(surface_list%n_psi)

if (allocated(surface_list%flux_surfaces)) then
   call tr_unregister_mem(sizeof(surface_list%flux_surfaces),"surface_list%flux_surfaces")
   deallocate(surface_list%flux_surfaces)
end if

allocate(surface_list%flux_surfaces(1:surface_list%n_psi))
call tr_register_mem(sizeof(surface_list%flux_surfaces),"surface_list%flux_surfaces")

do j=1, surface_list%n_psi
  surface_list%flux_surfaces(j)%n_pieces = 0
  surface_list%flux_surfaces(j)%elm      = 0
  surface_list%flux_surfaces(j)%s        = 0
  surface_list%flux_surfaces(j)%t        = 0
enddo

if (xpoint) then
  psi_xpoint(1) = ES%psi_xpoint(1)
  psi_xpoint(2) = ES%psi_xpoint(2)
  ! if we have a symmetric double-null, force the single separatrix
  if (ES%active_xpoint .eq. SYMMETRIC_XPOINT) then
    psi_xpoint(1) = (ES%psi_xpoint(1)+ES%psi_xpoint(2))/2.d0
    psi_xpoint(2) = psi_xpoint(1)
  endif
endif


do i=1, element_list%n_elements
       
  call psi_minmax(node_list,element_list,i,psimin,psimax)
 
  do j=1, surface_list%n_psi

    ifound = 0

!    if ((surface_list%psi_values(j) .ge. psimin) .and. (surface_list%psi_values(j) .le. psimax)) then

      do iv=1, 4

        ! --- For bi-cubic elements
        if (n_order .eq. 3) then

          im = MOD(iv,4) + 1
          n1 = element_list%element(i)%vertex(iv)
          n2 = element_list%element(i)%vertex(im)
         
          if (node_list%node(n1)%axis_node .and. node_list%node(n2)%axis_node) cycle
         
          is = mod(iv+1,2) + 2
         
          p1  =  node_list%node(n1)%values(1,1,1)  * element_list%element(i)%size(iv,1)
          dp1 =  node_list%node(n1)%values(1,is,1) * element_list%element(i)%size(iv,is)
          p4  =  node_list%node(n2)%values(1,1,1)  * element_list%element(i)%size(im,1)
          dp4 =  node_list%node(n2)%values(1,is,1) * element_list%element(i)%size(im,is)
         
          p2  = p1 + dp1
          p3  = p4 + dp4
         
          a3 = -        p1 + 3.d0 * p2 - 3.d0 * p3 + p4
          a2 = + 3.d0 * p1 - 6.d0 * p2 + 3.d0 * p3
          a1 = - 3.d0 * p1 + 3.d0 * p2
          a0 =          p1                                       - surface_list%psi_values(j)
         
          call SOLVP3(a0,a1,a2,a3,s,s2,s3,ifail)

        ! --- For higher order, cubic root finder need to be replaced by Newton methods
        else
          im = MOD(iv,4) + 1
          n1 = element_list%element(i)%vertex(iv)
          n2 = element_list%element(i)%vertex(im)
          if (node_list%node(n1)%axis_node .and. node_list%node(n2)%axis_node) cycle
          call newton_1D_find_value_on_element_side(node_list,element_list,i, iv, var_psi, surface_list%psi_values(j), n_found, st_found)
          if (n_found .gt. 3) write(*,*)'Warning, found more than 3 intersections!',i,iv,n_found,surface_list%psi_values(j)
          ! SOLVP3 reverses the search on sides 3 and 4
          if (iv .ge. 3) then
            st_found(1:n_found) = 1.d0 - st_found(1:n_found)
          endif
          s  = 99.
          s2 = 999.
          s3 = 9999.
          if (n_found .ge. 1) s  = st_found(1)
          if (n_found .ge. 2) s2 = st_found(2)
          if (n_found .ge. 3) s3 = st_found(3)
        endif
          
        if ((s .ge. 0.d0) .and. (s .le. 1.d0)) then

          ifound = ifound + 1

!         write(*,*) ' first solution : ',s,iv,i

          call flux_surface_add_point(node_list,element_list,surface_list,s,i,iv,ifound,r_psi,s_psi,dpsi_dr,dpsi_ds)

        endif

        if ((s2 .ge. 0.d0) .and. (s2 .le. 1.d0)) then

          ifound = ifound + 1

!         write(*,*) ' second solution : ',s2,iv,i

          call flux_surface_add_point(node_list,element_list,surface_list,s2,i,iv,ifound,r_psi,s_psi,dpsi_dr,dpsi_ds)

          if (abs(s3) .le. 1.d0)         write(*,*) ' WARNING another solution : ',s3,i

        endif

        if ((s3 .ge. 0.d0) .and. (s3 .le. 1.d0)) then

          ifound = ifound + 1

!         write(*,*) ' third solution : ',s3,iv,i

          call flux_surface_add_point(node_list,element_list,surface_list,s3,i,iv,ifound,r_psi,s_psi,dpsi_dr,dpsi_ds)

        endif

      enddo ! end of 4 edges
  
      
      ! --- Normally this should not really happen, but it does, ie. one of the end points is on a corner...
      ! --- Or the surface is tangential to the edge of the element
      if (ifound .eq. 3) then
        do k=1,3
          call interp_RZ(node_list,element_list,i,r_psi(k),s_psi(k),RRg(k),ZZg(k))
          dpsi_dr_copy(k) = dpsi_dr(k)
          dpsi_ds_copy(k) = dpsi_ds(k)
          r_psi_copy(k)   = r_psi(k)  
          s_psi_copy(k)   = s_psi(k)  
          tht_copy(k)     = tht(k)    
        enddo
        distance_max = 0.d0
        do k=1,3
          kp1 = mod(k,3)+1
          distance = sqrt( (RRg(k)-RRg(kp1))**2 + (ZZg(k)-ZZg(kp1))**2 )
          if (distance .gt. distance_max) then
            distance_max = distance
            k_keep = k
          endif
        enddo
        k = k_keep
        kp1 = mod(k,3)+1
        dpsi_dr(1) = dpsi_dr_copy(k) ; dpsi_dr(2) = dpsi_dr_copy(kp1) 
        dpsi_ds(1) = dpsi_ds_copy(k) ; dpsi_ds(2) = dpsi_ds_copy(kp1) 
        r_psi(1)   = r_psi_copy(k)   ; r_psi(2)   = r_psi_copy(kp1)   
        s_psi(1)   = s_psi_copy(k)   ; s_psi(2)   = s_psi_copy(kp1)   
        tht(1)     = tht_copy(k)     ; tht(2)     = tht_copy(kp1)     
        ifound = 2
      endif
      
      if (ifound .eq. 2) then

        call flux_surface_add_line(node_list,element_list,surface_list,i,j,r_psi(1:2),s_psi(1:2),dpsi_dr(1:2),dpsi_ds(1:2))

      elseif (ifound .eq. 4) then
      
! complicated : 2 line pieces but which point belongs to which line piece?

        r_av = (r_psi(1)+r_psi(2)+r_psi(3)+r_psi(4))/4.d0
        s_av = (s_psi(1)+s_psi(2)+s_psi(3)+s_psi(4))/4.d0

        tht(1:4) = atan2(s_psi(1:4)-s_av,r_psi(1:4)-r_av)

        where (tht .lt. 0.d0) tht = tht + 2.d0*PI

        itht(1)= 1; itht(2) = 2; itht(3) = 3; itht(4) = 4

        if (tht(2) .lt. tht(1)) then
          itmp = itht(1); itht(1) = itht(2) ; itht(2) = itmp;
        endif
        if (tht(4) .lt. tht(3)) then
          itmp = itht(3); itht(3) = itht(4) ; itht(4) = itmp;
        endif
        if (tht(itht(3)) .lt. tht(itht(2))) then
          itmp = itht(2); itht(2) = itht(3) ; itht(3) = itmp;
        endif
        if (tht(itht(2)) .lt. tht(itht(1))) then
          itmp = itht(1); itht(1) = itht(2) ; itht(2) = itmp;
        endif
        if (tht(itht(4)) .lt. tht(itht(3))) then
          itmp = itht(3); itht(3) = itht(4) ; itht(4) = itmp;
        endif
        if (tht(itht(3)) .lt. tht(itht(2))) then
          itmp = itht(2); itht(2) = itht(3) ; itht(3) = itmp;
        endif
        

        if ((xpoint) .and. (     ((i .eq. ES%i_elm_xpoint(1)) .and. (xcase .ne. UPPER_XPOINT) .and. (surface_list%psi_values(j) .eq. psi_xpoint(1)) )  &
                            .or. ((i .eq. ES%i_elm_xpoint(2)) .and. (xcase .ne. LOWER_XPOINT) .and. (surface_list%psi_values(j) .eq. psi_xpoint(2)) )  ) ) then

          call flux_surface_add_line(node_list,element_list,surface_list,i,j,r_psi(itht(1:3:2)), &
                   s_psi(itht(1:3:2)),dpsi_dr(itht(1:3:2)),dpsi_ds(itht(1:3:2)))
          call flux_surface_add_line(node_list,element_list,surface_list,i,j,r_psi(itht(2:4:2)), &
                                s_psi(itht(2:4:2)),dpsi_dr(itht(2:4:2)),dpsi_ds(itht(2:4:2)))
          ! --- Because of the saddle, near the Xpoint, the derivatives of the spline can be very noisy
          ! --- Sometimes (usually) it leads to surface pieces going back on themselves, clearly incoherent
          ! --- The safest is to set the derivatives to zero
          ! --- The only problem is that we cannot ensure the two pieces will cross at the Xpoint
          ! --- However, with the derivatives taken from flux_surface_add_line, this is not ensured either, so...
          k = surface_list%flux_surfaces(j)%n_pieces-1
          surface_list%flux_surfaces(j)%s(2,k) = 0.d0
          surface_list%flux_surfaces(j)%s(4,k) = 0.d0
          surface_list%flux_surfaces(j)%t(2,k) = 0.d0
          surface_list%flux_surfaces(j)%t(4,k) = 0.d0
          k = surface_list%flux_surfaces(j)%n_pieces
          surface_list%flux_surfaces(j)%s(2,k) = 0.d0
          surface_list%flux_surfaces(j)%s(4,k) = 0.d0
          surface_list%flux_surfaces(j)%t(2,k) = 0.d0
          surface_list%flux_surfaces(j)%t(4,k) = 0.d0

        else

          ! This is a little tricky, we look if the element is a neighbour of one of the Xpoints
          Xneigh = 0
          do k=1,4
            i_neigh = element_list%element(i)%neighbours(k)
            if( (xcase .ne. UPPER_XPOINT) .and. (i_neigh .eq. ES%i_elm_xpoint(1)) ) then
              Xneigh = 1
              exit
            endif
            if( (xcase .ne. LOWER_XPOINT) .and. (i_neigh .eq. ES%i_elm_xpoint(2)) ) then
              Xneigh = 2
              exit
            endif
          enddo
          ! If it is a neighbour, then record all four intersections (also do that for cases where
          ! the element is ES%i_elm_xpoint, but the flux surface is not the LCFS)
          if( (Xneigh .gt. 0) &
            .or. ((i .eq. ES%i_elm_xpoint(1)) .and. (xcase .ne. UPPER_XPOINT)) & 
            .or. ((i .eq. ES%i_elm_xpoint(2)) .and. (xcase .ne. LOWER_XPOINT)) ) then
            do k=1,4
              call interp_RZ(node_list,element_list,i,r_psi(k),s_psi(k),RRg(k),ZZg(k))
            enddo
          endif
          ! Then, look if the element is above/below or right/left of ES%i_elm_xpoint, 
          ! and then reorder the points 1,2,3,4 so that 1,2 are always right/above Xpoint,
          ! and 3,4 are always left/below Xpoint
          if(Xneigh .gt. 0) then
            if( (maxval(RRg) .gt. ES%R_xpoint(Xneigh)) .and. (minval(RRg) .lt. ES%R_xpoint(Xneigh)) ) then
              icount = 0
              do k=1,4
                if(RRg(k) .gt. ES%R_xpoint(Xneigh)) then
                  icount = icount + 1
                  itht(icount) = k
                endif
              enddo
              do k=1,4
                if(RRg(k) .lt. ES%R_xpoint(Xneigh)) then
                  icount = icount + 1
                  itht(icount) = k
                endif
              enddo
            else
              icount = 0
              do k=1,4
                if(ZZg(k) .gt. ES%Z_xpoint(Xneigh)) then
                  icount = icount + 1
                  itht(icount) = k
                endif
              enddo
              do k=1,4
                if(ZZg(k) .lt. ES%Z_xpoint(Xneigh)) then
                  icount = icount + 1
                  itht(icount) = k
                endif
              enddo
            endif
          endif
          
          ! In the case where the element actually is ES%i_elm_xpoint, 
          ! but the flux surface is not the LCFS, we need to check if the line is right&left
          ! or above&below the Xpoint
          if( (Xneigh .eq. 0) &
            .and. (    ((i .eq. ES%i_elm_xpoint(1)) .and. (xcase .ne. UPPER_XPOINT)) &
                  .or. ((i .eq. ES%i_elm_xpoint(2)) .and. (xcase .ne. LOWER_XPOINT)) ) ) then
            if(i .eq. ES%i_elm_xpoint(1)) Xneigh = 1
            if(i .eq. ES%i_elm_xpoint(2)) Xneigh = 2
            if(surface_list%psi_values(j) .gt. psi_xpoint(Xneigh)) then
              icount = 0
              do k=1,4
                if(RRg(k) .gt. ES%R_xpoint(Xneigh)) then
                  icount = icount + 1
                  itht(icount) = k
                endif
              enddo
              do k=1,4
                if(RRg(k) .lt. ES%R_xpoint(Xneigh)) then
                  icount = icount + 1
                  itht(icount) = k
                endif
              enddo
            else
              icount = 0
              do k=1,4
                if(ZZg(k) .gt. ES%Z_xpoint(Xneigh)) then
                  icount = icount + 1
                  itht(icount) = k
                endif
              enddo
              do k=1,4
                if(ZZg(k) .lt. ES%Z_xpoint(Xneigh)) then
                  icount = icount + 1
                  itht(icount) = k
                endif
              enddo
            endif
          endif

          ! Then add the lines
          call flux_surface_add_line(node_list,element_list,surface_list,i,j,r_psi(itht(1:2)), &
                              s_psi(itht(1:2)),dpsi_dr(itht(1:2)),dpsi_ds(itht(1:2)))
          call flux_surface_add_line(node_list,element_list,surface_list,i,j,r_psi(itht(3:4)), &
                              s_psi(itht(3:4)),dpsi_dr(itht(3:4)),dpsi_ds(itht(3:4)))

        endif

      endif
     
!    endif

  enddo

enddo

return
end

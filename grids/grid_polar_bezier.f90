!> Defines a polar grid using Bezier finite elements (using the HELENA
!! cubic Hermite elements formulation)
subroutine grid_polar_bezier(Rgeo,Zgeo,amin,acentre,angle_start,fbnd,fpsi,mf,nr,np,node_list,element_list)

use constants
use tr_module
use mod_parameters
use data_structure
use mod_neighbours, only: update_neighbours
use phys_module, only: psi_axis_init, XR_r, SIG_r, XR_tht, SIG_tht, fix_axis_nodes, force_central_node

implicit none

! --- Routine parameters
real*8,                  intent(in)    :: Rgeo           !< R-position of geometric center
real*8,                  intent(in)    :: Zgeo           !< Z-position of geometric center
real*8,                  intent(in)    :: amin           !< minor radius
real*8,                  intent(in)    :: acentre        !< inner radius of grid
real*8,                  intent(in)    :: angle_start    !< poloidal angle of first element
real*8,                  intent(in)    :: fbnd(*)        !< Fourier series describing the radius as function
                                                         !!   of the poloidal angle
real*8,                  intent(in)    :: fpsi(*)        !< Fourier series of flux at the boundary
integer,                 intent(in)    :: mf             !< Number of Fourier modes in fbnd and fpsi
integer,                 intent(in)    :: nr             !< number of radial points, (nr-1) elements
integer,                 intent(in)    :: np             !< number of poloidal points, np elements
type(type_node_list),    intent(inout) :: node_list      !< list of nodes with grid information
type(type_element_list), intent(inout) :: element_list   !< list of elements with element information

! --- local variables
real*8 :: acentre2, radius2
real*8              :: si
real*8, allocatable :: RR(:,:),ZZ(:,:),PSI(:,:)
real*8              :: dt, ds, thtj, radius, rm, drm, drmt, drmtr, angle, psi_axis
real*8              :: delta_rm, delta_zm, delta_rp, delta_zp, dir_2, dir_3
integer             :: i, j, m, index, index0, node, k, iv, ivp, ivm, node_iv, node_ivp, node_ivm,i_sons
integer             :: n_element_start, n_node_start, n_index_start
real*8              :: abltg(3), dr_ds, dtht_dt
real*8, allocatable :: S1(:), S2(:), SP1(:), SP2(:), SP3(:), SP4(:)
real*8, allocatable :: T1(:), T2(:), TP1(:), TP2(:), TP3(:), TP4(:)
real*8, external    :: spwert
logical             :: skip_update_neighbours
logical             :: doing_polar_square


call tr_allocate(RR,1,4,1,nr*np,"RR",CAT_GRID)
call tr_allocate(ZZ,1,4,1,nr*np,"ZZ",CAT_GRID)
call tr_allocate(PSI,1,4,1,nr*np,"PSI",CAT_GRID)

dt = 2.d0*pi/real(np)
ds = 1.d0/real(nr-1)

n_element_start  = element_list%n_elements
n_node_start     = node_list%n_nodes
do i=n_node_start+1,n_nodes_max
  node_list%node(i)%x        = 0.d0
  node_list%node(i)%values   = 0.d0
  node_list%node(i)%index    = 0
  node_list%node(i)%boundary = 0
enddo

skip_update_neighbours = .false.
if ( n_element_start /= 0 ) skip_update_neighbours = .true. ! In such a case, the call to update_neighbours needs to be done in another routine like grid_bezier_square_polar

do i=n_element_start+1,n_elements_max
  element_list%element(i)%vertex     = 0
  element_list%element(i)%size       = 0.d0
  element_list%element(i)%neighbours = 0
enddo

n_index_start = 0
do i=1,n_node_start
  n_index_start = max(n_index_start,maxval(node_list%node(i)%index(:)))
enddo
doing_polar_square = .false.
if (n_index_start .gt. 0) doing_polar_square = .true.

write(*,*) '*************************************'
write(*,*) '*        grid_polar_bezier          *'
write(*,*) '*************************************'
write(*,*) ' existing number of elements          : ',n_element_start
write(*,*) ' existing number of nodes             : ',n_node_start
write(*,*) ' index_start                          : ',n_index_start

psi_axis = psi_axis_init ! Initial guess for Psi at the magnetic axis

call tr_allocate(S1,1,nr,"S1",CAT_GRID)
call tr_allocate(S2,1,nr,"S2",CAT_GRID)
call tr_allocate(SP1,1,nr,"SP1",CAT_GRID)
call tr_allocate(SP2,1,nr,"SP2",CAT_GRID)
call tr_allocate(SP3,1,nr,"SP3",CAT_GRID)
call tr_allocate(SP4,1,nr,"SP4",CAT_GRID)
S2 = 0.d0
call tr_allocate(T1,1,np+1,"T1",CAT_GRID)
call tr_allocate(T2,1,np+1,"T2",CAT_GRID)
call tr_allocate(TP1,1,np+1,"TP1",CAT_GRID)
call tr_allocate(TP2,1,np+1,"TP2",CAT_GRID)
call tr_allocate(TP3,1,np+1,"TP3",CAT_GRID)
call tr_allocate(TP4,1,np+1,"TP4",CAT_GRID)
T2 = 0.d0
do i=1,nr
  S1(i) = real(i-1)/real(nr-1)
enddo
do j=1,np+1
  T1(j) = real(j-1)/real(np)
enddo

call meshac2(nr,S2,XR_r(1),XR_r(2),SIG_r(1),SIG_r(2),0.6d0,1.0d0)
call spline(nr,S1,S2,0.d0,0.d0,2,SP1,SP2,SP3,SP4)

call meshac2(np+1,T2,XR_tht(1),XR_tht(2),SIG_tht(1),SIG_tht(2),0.6d0,1.0d0)
call spline(np+1,T1,T2,0.d0,0.d0,2,TP1,TP2,TP3,TP4)

acentre2 = acentre * 2.d0 / fbnd(1) !###

do i=1,nr

  si = spwert(nr,S1(i),SP1,SP2,SP3,SP4,S1,ABLTG)

  radius = ( acentre2 + (1.d0-acentre2) * si )
  radius2= si
  dr_ds  = (1.d0-acentre2) * abltg(1)
  
  do  j=1,np
  
    node   = np*(i-1) + j

    thtj     = angle_start + spwert(np+1,T1(j),TP1,TP2,TP3,TP4,T1,ABLTG) * 2.d0 * PI
    dtht_dt  = abltg(1)
    
    RR(1,node) = Rgeo + amin * radius * fbnd(1) * cos(thtj) / 2.d0
    RR(2,node) =        amin *          fbnd(1) * cos(thtj) / 2.d0
    RR(3,node) =      - amin * radius * fbnd(1) * sin(thtj) / 2.d0
    RR(4,node) =      - amin          * fbnd(1) * sin(thtj) / 2.d0
    ZZ(1,node) = Zgeo + amin * radius * fbnd(1) * sin(thtj) / 2.d0
    ZZ(2,node) =        amin *          fbnd(1) * sin(thtj) / 2.d0
    ZZ(3,node) =        amin * radius * fbnd(1) * cos(thtj) / 2.d0
    ZZ(4,node) =        amin *          fbnd(1) * cos(thtj) / 2.d0

    PSI(1,node) =        radius**8 * fpsi(1) / 2.d0 + psi_axis*(1.d0 - radius**2)
    PSI(2,node) =  8.d0* radius**7 * fpsi(1) / 2.d0 - psi_axis*2.d0* radius
    PSI(3,node) =  0.d0
    PSI(4,node) =  0.d0

!   write(*,'(A,2i6,12e16.8)') ' PSI(1,node) : ', node,1,PSI(1,node),radius,thtj,fpsi(1),fpsi(2)

    do m = 2, mf/2

      if (m .eq. 2) then
        rm   = radius2 * ( fbnd(2*M-1) * cos((M-1)*THTJ)           + fbnd(2*M) * sin((M-1)*THTJ) )
        drm  =           ( fbnd(2*M-1) * cos((M-1)*THTJ)           + fbnd(2*M) * sin((M-1)*THTJ) )
        drmt = radius2 * (-fbnd(2*M-1) * (M-1)*sin((M-1)*THTJ)     + fbnd(2*M) * (M-1)*cos((M-1)*THTJ))
        drmtr=           (-fbnd(2*M-1) * (M-1)*sin((M-1)*THTJ)     + fbnd(2*M) * (M-1)*cos((M-1)*THTJ))
      else
        rm   =      radius2**(M-1) * ( fbnd(2*M-1) * cos((M-1)*THTJ)       + fbnd(2*M) * sin((M-1)*THTJ) )
        drm  =(M-1)*radius2**(M-2) * ( fbnd(2*M-1) * cos((M-1)*THTJ)       + fbnd(2*M) * sin((M-1)*THTJ) )
        drmt =      radius2**(M-1) * (-fbnd(2*M-1) * (M-1)*sin((M-1)*THTJ) + fbnd(2*M) *(M-1)*cos((M-1)*THTJ))
        drmtr=(M-1)*radius2**(M-2) * (-fbnd(2*M-1) * (M-1)*sin((M-1)*THTJ) + fbnd(2*M) *(M-1)*cos((M-1)*THTJ))
      endif
      RR(1,node) = RR(1,node) + amin * rm  * cos(thtj)
      ZZ(1,node) = ZZ(1,node) + amin * rm  * sin(thtj)
      RR(2,node) = RR(2,node) + amin * drm * cos(thtj)
      ZZ(2,node) = ZZ(2,node) + amin * drm * sin(thtj)
      RR(3,node) = RR(3,node) - amin * rm  * sin(thtj) + amin * drmt  * cos(thtj)
      ZZ(3,node) = ZZ(3,node) + amin * rm  * cos(thtj) + amin * drmt  * sin(thtj)
      RR(4,node) = RR(4,node) - amin * drm * sin(thtj) + amin * drmtr * cos(thtj)
      ZZ(4,node) = ZZ(4,node) + amin * drm * cos(thtj) + amin * drmtr * sin(thtj) 

      PSI(1,node) = PSI(1,node) + radius**8              * (   fpsi(2*m-1) *            cos((m-1)*thtj)       &
                                                             + fpsi(2*m)   *            sin((m-1)*thtj) )
      PSI(2,node) = PSI(2,node) + 8.d0 * radius**7 * drm * (   fpsi(2*m-1) *            cos((m-1)*thtj)       &
                                                             + fpsi(2*m)   *            sin((m-1)*thtj) )
      PSI(3,node) = PSI(3,node) + radius**8              * ( - fpsi(2*m-1) * float(m-1)*sin((m-1)*thtj)       &
                                                             + fpsi(2*m)   * float(m-1)*cos((m-1)*thtj) )
      PSI(4,node) = PSI(4,node) + 8.d0 * radius**7 * drm * ( - fpsi(2*m-1) * float(m-1)*sin((m-1)*thtj)       &
                                                             + fpsi(2*m)   * float(m-1)*cos((m-1)*thtj) )

    enddo

    RR(2,node)  = RR(2,node)  * ds/2.d0 * dr_ds
    RR(3,node)  = RR(3,node)  * dt/2.d0 * dtht_dt
    RR(4,node)  = RR(4,node)  * ds/2.d0 * dt/2.d0 * dr_ds * dtht_dt
    ZZ(2,node)  = ZZ(2,node)  * ds/2.d0 * dr_ds
    ZZ(3,node)  = ZZ(3,node)  * dt/2.d0 * dtht_dt
    ZZ(4,node)  = ZZ(4,node)  * ds/2.d0 * dt/2.d0 * dr_ds * dtht_dt
    PSI(2,node) = PSI(2,node) * ds/2.d0 * dr_ds
    PSI(3,node) = PSI(3,node) * dt/2.d0 * dtht_dt
    PSI(4,node) = PSI(4,node) * ds/2.d0 * dt/2.d0 * dr_ds * dtht_dt

  enddo
enddo

element_list%n_elements  = n_element_start  + (nr-1)*np
node_list%n_nodes        = n_node_start     + nr*np

if ( node_list%n_nodes > n_nodes_max ) then
  write(*,*) 'ERROR in grid_polar_bezier: hard-coded parameter n_nodes_max is too small'
  stop
else if ( element_list%n_elements > n_elements_max ) then
  write(*,*) 'ERROR in grid_polar_bezier: hard-coded parameter n_elements_max is too small'
  stop
end if


do i=1,nr-1

 do j=1,np
   node  = np*(i-1) + j
   index = n_element_start + node
   element_list%element(index)%vertex(1) = n_node_start + (i-1)*np + j
   element_list%element(index)%vertex(4) = n_node_start + (i-1)*np + j + 1
   element_list%element(index)%vertex(3) = n_node_start + (i  )*np + j + 1
   element_list%element(index)%vertex(2) = n_node_start + (i  )*np + j

   if (j .eq. np) then          
                 element_list%element(Index)%vertex(4) = n_node_start + (i-1)*np + 1
                 element_list%element(Index)%vertex(3) = n_node_start + i*np    + 1
   endif    

           !Neighbours of the element (refinement procedure)

            if(i==1) then               
                 element_list%element(Index)%neighbours(4) = 0    
              else               
                 element_list%element(Index)%neighbours(4) = Index - np 
            end if       
            
            if(j==np) then
                 element_list%element(Index)%neighbours(3) = Index - np + 1  
              else
                 element_list%element(Index)%neighbours(3) = Index + 1       
            end if                  
            
            if(i==nr-1) then
                 element_list%element(Index)%neighbours(2) = 0   
              else   
                 element_list%element(Index)%neighbours(2) = Index + np 
            end if 
                
            if(j==1) then
                 element_list%element(Index)%neighbours(1) = Index + np -1 
              else   
                 element_list%element(Index)%neighbours(1) = Index -1       
            end if     
          
            ! Initialization of the genealogy  (refinement procedure)

            element_list%element(Index)%father = 0
            element_list%element(Index)%n_sons = 0
            do i_sons = 1, 4
                 element_list%element(Index)%sons(i_sons) = 0
            end do 
            
  enddo
enddo

 


!-------------------- translate cubic Hermite to Bezier parameters
!
!  type type_node                                      ! type definition of a node (i.e. a vertex)
!    real*8    :: x(n_order+1,ndim)                      ! x,y coordinates of points and additional nodal geometry
!    integer :: boundary                               ! = 1 for boundary nodes
!  endtype type_node                                   ! x(:,1) : position, x(:,2) : vector u, x(:,3) : vector v, x(4) : vector w
!
!  type type_element
!    integer :: vertex(n_vertex_max)
!    integer :: neighbours(n_vertex_max)
!    real*8    :: size(n_vertex_max,n_order+1)
!  endtype type_element
!-----------------------------------------------------------

do i=1,nr

 do j=1,np

   angle = 2.d0*PI * float(j-1)/float(np)

   index0 = np*(i-1) + j
   index  = n_node_start + np*(i-1) + j

   node_list%node(index)%X(1,1,1)        = RR(1,index0)
   node_list%node(index)%X(1,1,2)        = ZZ(1,index0)
   node_list%node(index)%values(1,1,1) = PSI(1,index0)

   node_list%node(index)%X(1,2,1)        = RR(2,index0)  * 2.d0/3.d0
   node_list%node(index)%X(1,2,2)        = ZZ(2,index0)  * 2.d0/3.d0
   node_list%node(index)%values(1,2,1) = PSI(2,index0) * 2.d0/3.d0

   node_list%node(index)%X(1,3,1)        = RR(3,index0)  * 2.d0/3.d0
   node_list%node(index)%X(1,3,2)        = ZZ(3,index0)  * 2.d0/3.d0
   node_list%node(index)%values(1,3,1) = PSI(3,index0) * 2.d0/3.d0

   node_list%node(index)%X(1,4,1)        = RR(4,index0)  * 4.d0/9.d0
   node_list%node(index)%X(1,4,2)        = ZZ(4,index0)  * 4.d0/9.d0
   node_list%node(index)%values(1,4,1) = PSI(4,index0) * 4.d0/9.d0

   node_list%node(index)%boundary = 0
   if (i .eq. nr) node_list%node(index)%boundary = 2

   node_list%node(index)%axis_node = .false.
   if ( fix_axis_nodes .and. (.not. doing_polar_square) .and. (i .eq. 1) ) node_list%node(index)%axis_node = .true.

   if (force_central_node .and. (.not. doing_polar_square) .and. (i.eq.1)) then

     node_list%node(index)%index(1) = 1

     if (j.eq.1) n_index_start = n_index_start + 1

     node_list%node(index)%index(2) = n_index_start + 1
     node_list%node(index)%index(3) = n_index_start + 2
     node_list%node(index)%index(4) = n_index_start + 3
     n_index_start = n_index_start + n_order

   else
     do k=1,n_order+1
       node_list%node(index)%index(k) = n_index_start + k
     enddo
     n_index_start = n_index_start + n_order+1
   endif
  
   node_list%node(index)%constrained=.false.
 enddo

enddo


do k=n_element_start+1 , element_list%n_elements   ! fill in the size of the elements

 do iv = 1, 4                    ! over 4 corners of an element

   ivp = mod(iv,4)   + 1         ! vertex with index one higher
   ivm = mod(iv+2,4) + 1         ! vertex with index one below

   node_iv  = element_list%element(k)%vertex(iv)
   node_ivp = element_list%element(k)%vertex(ivp)
   node_ivm = element_list%element(k)%vertex(ivm)

   if ((iv .eq. 1) .or. (iv .eq.3)) then

     delta_Rp = node_list%node(node_ivp)%X(1,1,1) - node_list%node(node_iv)%X(1,1,1)
     delta_Zp = node_list%node(node_ivp)%X(1,1,2) - node_list%node(node_iv)%X(1,1,2)
     dir_2    = delta_Rp * node_list%node(node_iv)%X(1,2,1) + delta_Zp * node_list%node(node_iv)%X(1,2,2)

     delta_Rm = node_list%node(node_ivm)%X(1,1,1) - node_list%node(node_iv)%X(1,1,1)
     delta_Zm = node_list%node(node_ivm)%X(1,1,2) - node_list%node(node_iv)%X(1,1,2)
     dir_3    = delta_Rm * node_list%node(node_iv)%X(1,3,1) + delta_Zm * node_list%node(node_iv)%X(1,3,2)

   else

     delta_Rp = node_list%node(node_ivp)%X(1,1,1) - node_list%node(node_iv)%X(1,1,1)
     delta_Zp = node_list%node(node_ivp)%X(1,1,2) - node_list%node(node_iv)%X(1,1,2)
     dir_3    = delta_Rp * node_list%node(node_iv)%X(1,3,1) + delta_Zp * node_list%node(node_iv)%X(1,3,2)

     delta_Rm = node_list%node(node_ivm)%X(1,1,1) - node_list%node(node_iv)%X(1,1,1)
     delta_Zm = node_list%node(node_ivm)%X(1,1,2) - node_list%node(node_iv)%X(1,1,2)
     dir_2    = delta_Rm * node_list%node(node_iv)%X(1,2,1) + delta_Zm * node_list%node(node_iv)%X(1,2,2)

   endif

   if (dir_2 .ne. 0.d0) then
     dir_2 = dir_2 / abs(dir_2)
   else
     dir_2 = 1.d0
   endif
   if (dir_3 .ne. 0.d0) then
     dir_3 = dir_3 / abs(dir_3)
   else
     dir_3 = -1.d0
     if (iv.eq.1) dir_3 = 1.d0              ! admittedly not very elegant
   endif

   element_list%element(k)%size(iv,1) = 1.d0
   element_list%element(k)%size(iv,2) = dir_2
   element_list%element(k)%size(iv,3) = dir_3
   element_list%element(k)%size(iv,4) = element_list%element(k)%size(iv,2) * element_list%element(k)%size(iv,3)
   if (fix_axis_nodes) then
      j = element_list%element(k)%vertex(iv)
      if (node_list%node(j)%axis_node) then
        element_list%element(k)%size(iv,3) = 0.d0
        element_list%element(k)%size(iv,4) = 0.d0
      endif
   endif

!   if ((RR(2,node_iv)**2 + ZZ(2,node_iv)**2) .eq. 0.) element_list%element(k)%size(iv,2) = dir_2
!   if ((RR(3,node_iv)**2 + ZZ(3,node_iv)**2) .eq. 0.) element_list%element(k)%size(iv,3) = dir_3
!   if ((RR(4,node_iv)**2 + ZZ(4,node_iv)**2) .eq. 0.) element_list%element(k)%size(iv,4) = dir_2 * dir_3

!    write(*,'(2i5,12e16.8)') k,iv,element_list%element(k)%size(iv,1:4)

 enddo

enddo

if ( .not. skip_update_neighbours ) call update_neighbours(node_list,element_list, force_rtree_initialize=.true.)
return
end subroutine grid_polar_bezier

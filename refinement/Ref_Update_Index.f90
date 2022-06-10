subroutine Ref_Update_Index( element_list,node_list)



 use mod_parameters
 use data_structure

implicit none

type (type_node_list)    :: node_list
type (type_element_list)  :: element_list

 
!integer, dimension(node_list%n_nodes) :: active_node
integer, dimension(4)                 ::son1,son2,son3,son4,son5
integer, dimension(5)                 :: ap
integer                               ::n_active_nodes
integer                               :: i,j,k,l,p,iv,ielm,inode,ov
 
   !n_active_nodes = 0
  
  !do i = 1, node_list%n_nodes
       
      
            !if((node_list%node(i)%constrained==.false.) )  then
                ! n_active_nodes = n_active_nodes + 1
	        ! active_node(i) = n_active_nodes	
           !do k=1,n_degrees
            !node_list%node(i)%index(k) =  n_degrees*(active_node(i)-1)+k 
           
           !enddo	
        ! else
	     ! active_node(i) = 0							                          
       !end if 
  
  !end do
!return
  !write(*,*) 'Noeuds actifs = ', n_active_nodes,node_list%n_nodes
 n_active_nodes = 0

 do ielm=1,element_list%n_elements

 if (element_list%element(ielm)%n_gen .eq.0)then
   
   
      inode=element_list%element(ielm)%vertex(1)
      n_active_nodes = n_active_nodes + 1
      !active_node(inode) = n_active_nodes
          do ov=1,n_degrees
            node_list%node(inode)%index(ov) =  n_degrees*(n_active_nodes-1)+ov
           enddo
    
  if(element_list%element(ielm)%n_sons .ne.0)then

   
  
   
   do iv=1,5
   
       ap(iv)=element_list%element(ielm)%contain_node(iv)
    if (ap(iv).ne.0) then
     if ((.not. node_list%node(ap(iv))%constrained) ) then
     n_active_nodes = n_active_nodes + 1
     !active_node(ap(iv)) = n_active_nodes
           do ov=1,n_degrees
            node_list%node(ap(iv))%index(ov) =  n_degrees*(n_active_nodes-1)+ov 
           enddo
     endif
    endif
   enddo
   
  do i=1,element_list%element(ielm)%n_sons
    son1(i)=element_list%element(ielm)%sons(i)

   
    do iv=1,5
   
       ap(iv)=element_list%element(son1(i))%contain_node(iv)
    if (ap(iv).ne.0) then
     if ((.not. node_list%node(ap(iv))%constrained) ) then
     n_active_nodes = n_active_nodes + 1
     !active_node(ap(iv)) = n_active_nodes
           do ov=1,n_degrees
            node_list%node(ap(iv))%index(ov) =  n_degrees*(n_active_nodes-1)+ov 
           enddo
     endif
    endif
   enddo
    
  do j=1,element_list%element(son1(i))%n_sons
    son2(j)=element_list%element(son1(i))%sons(j)
   
   do iv=1,5
   
       ap(iv)=element_list%element(son2(j))%contain_node(iv)
    if (ap(iv).ne.0) then
     if ((.not. node_list%node(ap(iv))%constrained) ) then
     n_active_nodes = n_active_nodes + 1
     !active_node(ap(iv)) = n_active_nodes
           do ov=1,n_degrees
            node_list%node(ap(iv))%index(ov) =  n_degrees*(n_active_nodes-1)+ov 
           enddo
     endif
    endif
   enddo

 
  do k=1,element_list%element(son2(j))%n_sons
    son3(k)=element_list%element(son2(j))%sons(k)
   
   do iv=1,5
   
       ap(iv)=element_list%element(son3(k))%contain_node(iv)
    if (ap(iv).ne.0) then
     if ((.not. node_list%node(ap(iv))%constrained) ) then
     n_active_nodes = n_active_nodes + 1
    ! active_node(ap(iv)) = n_active_nodes
           do ov=1,n_degrees
            node_list%node(ap(iv))%index(ov) =  n_degrees*(n_active_nodes-1)+ov 
           enddo
     endif
    endif
   enddo
  
     do l=1,element_list%element(son3(k))%n_sons
    son4(l)=element_list%element(son3(k))%sons(l)
   
   do iv=1,5
   
       ap(iv)=element_list%element(son4(l))%contain_node(iv)
    if (ap(iv).ne.0) then
     if ((.not. node_list%node(ap(iv))%constrained) ) then
     n_active_nodes = n_active_nodes + 1
     !active_node(ap(iv)) = n_active_nodes
           do ov=1,n_degrees
            node_list%node(ap(iv))%index(ov) =  n_degrees*(n_active_nodes-1)+ov 
           enddo
     endif
    endif
   enddo
  
   
  !************
  
     do p=1,element_list%element(son4(l))%n_sons
    son5(p)=element_list%element(son4(l))%sons(p)
   
   do iv=1,5
   
       ap(iv)=element_list%element(son5(p))%contain_node(iv)
    if (ap(iv).ne.0) then
     if ((.not. node_list%node(ap(iv))%constrained) ) then
     n_active_nodes = n_active_nodes + 1
    ! active_node(ap(iv)) = n_active_nodes
           do ov=1,n_degrees
            node_list%node(ap(iv))%index(ov) =  n_degrees*(n_active_nodes-1)+ov 
           enddo
     endif
    endif
   enddo
  
     enddo
    enddo
   enddo
  enddo
 enddo
   
  endif

    
 endif
           
 enddo
   
 do ielm=1,element_list%n_elements
  if (element_list%element(ielm)%n_gen .eq.0)then
   if(element_list%element(ielm)%neighbours(2).eq.0) then
      inode=element_list%element(ielm)%vertex(2)
      n_active_nodes = n_active_nodes + 1
      !active_node(inode) = n_active_nodes
           do ov=1,n_degrees
            node_list%node(inode)%index(ov) =  n_degrees*(n_active_nodes-1)+ov 
           enddo
    endif
  endif
  enddo
 write(*,*) 'Noeuds actifs = ', n_active_nodes,node_list%n_nodes

  return


end subroutine Ref_Update_Index

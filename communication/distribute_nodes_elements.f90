subroutine distribute_nodes_elements(my_id,m_cpu,n_cpu,node_list,element_list, direct_construction, &
                                    local_elms, n_local_elms, n_dof, index_min, index_max, restart, freeboundary)
!---------------------------------------------------------------------------------------------
! subroutine divides the nodes (not their individual dof) over n_cpu equal parts
!            builds local_elms, contain all elements with at least one node with 
!            one index between index_min and index_max
!---------------------------------------------------------------------------------------------
use data_structure
use mod_integer_types

implicit none

type (type_node_list)    :: node_list
type (type_element_list) :: element_list
type (type_surface_list) :: flux_list

logical, parameter :: DEBUG = .false.

integer               :: local_elms(*)
integer               :: my_id, n_cpu, m_cpu, n_local_elms, inode
integer(kind=int_all) :: n_dof
integer               :: index_total
integer               :: inext, i,j, k, iv,index1
integer               :: index_min(*), index_max(*)
logical               :: restart, freeboundary

logical :: elm_is_local, direct_construction
!integer, dimension(node_list%n_nodes) :: active_node
!integer                               :: n_active_nodes
integer :: mpi_distr_count, ib, l_index, ik

if (my_id .eq. 0) then 
  if (.not. direct_construction) then
    write(*,*) '************************************'
    write(*,*) '* distributing nodes global matrix *'
    write(*,*) '************************************'
  else
    write(*,*) '**************************************'
    write(*,*) '* distributing nodes harmonic matrix *'
    write(*,*) '**************************************'
  endif
endif
 

index_total = -1
do inode=1,node_list%n_nodes
  index_total = max(index_total,maxval(node_list%node(inode)%index))
enddo

index_min(1:n_cpu) = 0
index_max(1:n_cpu) = 0

!----------------------------- must really take into account the number of elements contributing to each node

if (.not. direct_construction) then ! global matrix construction
  
  index_min(1) = 1
  do i=1,n_cpu
    index_max(i) = (i * index_total) / n_cpu
  enddo
  do i=2,n_cpu
    index_min(i) = index_max(i-1) + 1
  enddo
  if (my_id .eq. n_cpu-1) index_max(my_id+1) = index_total
  if (DEBUG) write(*,'(A,3i6)') ' index_min,index_max:',my_id,index_min(my_id+1),index_max(my_id+1)
  
else ! harmonic matrix "direct" construction
  
  if (mod(my_id,m_cpu) .eq. 0) index_min(my_id+1) = 1
  do i=1,n_cpu
    index_max(i) = ((mod(i-1,m_cpu)+1) * index_total) / m_cpu
  enddo
  do i=2,n_cpu
    if (mod(i-1,m_cpu) .ne. 0) then
      index_min(i) = index_max(i-1) + 1
    endif
  enddo
  if (mod(my_id+1,m_cpu) .eq. 0) index_max(my_id+1) = index_total
  if (DEBUG) write(*,'(A,3i6)') ' index_min,index_max:',my_id,index_min(my_id+1),index_max(my_id+1)
  
end if 

n_dof = index_total * n_tor * n_var

!----------------------------------------------- find the elements that have a local node
inext = 0

do i = 1, element_list%n_elements

  ELM_is_local = .false.

  L_IV: do iv=1,n_vertex_max

    inode = element_list%element(i)%vertex(iv)

    do k=1, n_order+1

      if ( (node_list%node(inode)%index(k) .ge. index_min(my_id+1)) .and. &
           (node_list%node(inode)%index(k) .le. index_max(my_id+1)) ) then
        ELM_is_local = .true.
        exit L_IV
      endif
    enddo

    if(node_list%node(inode)%constrained) then
	     do j = 1, 2
		 index1 = node_list%node(inode)%parents(j)
		  do k=1, n_order+1

                       if ( (node_list%node(index1)%index(k) .ge. index_min(my_id+1)) .and. &
                        (node_list%node(index1)%index(k) .le. index_max(my_id+1)) ) then
                        ELM_is_local = .true.
                        exit L_IV
                      endif
                  enddo
	    end do     
	    	    
    end if

  end do L_IV
  if (ELM_is_local) then
     inext = inext + 1
     local_elms(inext) = i
  endif

enddo

n_local_ELMs = inext

if (DEBUG) write(*,'(i4,A,20i8)') my_id,' n_local_elms  : ',n_local_elms,element_list%n_elements

!--------------------- check distribution of boundary nodes 
mpi_distr_count=0
if (restart .and. freeboundary) then 
  do ib = 1, node_list%n_nodes	
    if ( node_list%node(ib)%boundary > 0 ) then
      do ik=1, n_order+1
        l_index = node_list%node(ib)%index(ik)
          if ((l_index .ge. index_min(my_id+1)) .and. (l_index .le. index_max(my_id+1))) then ! This MPI proc responsible?
            mpi_distr_count=mpi_distr_count+1
          end if
      end do
    end if
  end do
  if (DEBUG) write(*,*) 'task ', my_id, 'is responsible for ', mpi_distr_count, 'boundary nodes.'
  if (mpi_distr_count == 0) then 
    write(*,*) 'WARNING: boundary node indices seem to be unevenly distributed among the MPI tasks. To avoid this, use freeb_change_indices=.true. from the very beginning of your simulation."'
  end if 
end if
     


return
end

subroutine distribute_nodes_elements(my_id, m_cpu, index_size, node_list, element_list, direct_construction, &
                                    local_elms, n_local_elms, restart, freeboundary, a_mat)
!---------------------------------------------------------------------------------------------
! subroutine divides the nodes (not their individual dof) over index_size equal parts
!            builds local_elms, contain all elements with at least one node with 
!            one index between index_min and index_max
!---------------------------------------------------------------------------------------------
use data_structure
use mod_integer_types
use tr_module

implicit none

type (type_node_list)    :: node_list
type (type_element_list) :: element_list
type (type_surface_list) :: flux_list
type(type_SP_MATRIX)     :: a_mat

logical, parameter :: DEBUG = .false.

integer               :: local_elms(*)
integer, dimension(:), pointer :: index_min, index_max
integer               :: my_id, index_size, m_cpu, n_local_elms, inode
integer               :: index_total
integer               :: inext, i,j, k, iv,index1
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
    write(*,*) '*    distributing nodes PC matrix    *'
    write(*,*) '**************************************'
  endif
endif
 

index_total = -1
do inode=1,node_list%n_nodes
  index_total = max(index_total,maxval(node_list%node(inode)%index))
enddo
node_list%n_dof = index_total * n_tor * n_var

call tr_allocatep(index_min,1,index_size,"index_min",CAT_FEM)
call tr_allocatep(index_max,1,index_size,"index_max",CAT_FEM)

index_min(1:index_size) = 0
index_max(1:index_size) = 0


!----------------------------- must really take into account the number of elements contributing to each node

if (.not. direct_construction) then ! global matrix construction
  
  index_min(1) = 1
  do i=1,index_size
    index_max(i) = (i * index_total)/index_size
  enddo
  do i=2,index_size
    index_min(i) = index_max(i-1) + 1
  enddo
  if (my_id .eq. index_size-1) index_max(my_id+1) = index_total
  if (DEBUG) write(*,'(A,3i6)') ' index_min,index_max:', my_id, index_min(my_id+1), index_max(my_id+1)
  
else ! PC matrix "direct" construction
  
  if (mod(my_id,m_cpu) .eq. 0) index_min(my_id+1) = 1
  do i=1,index_size
    index_max(i) = ((mod(i-1,m_cpu)+1) * index_total)/m_cpu
  enddo
  do i=2,index_size
    if (mod(i-1,m_cpu) .ne. 0) then
      index_min(i) = index_max(i-1) + 1
    endif
  enddo
  if (mod(my_id+1,m_cpu) .eq. 0) index_max(my_id+1) = index_total
  if (DEBUG) write(*,'(A,3i6)') ' index_min,index_max:',my_id,index_min(my_id+1),index_max(my_id+1)
  
end if 

!----------------------------------------------- find the elements that have a local node
inext = 0

do i = 1, element_list%n_elements

  ELM_is_local = .false.

  L_IV: do iv=1,n_vertex_max

    inode = element_list%element(i)%vertex(iv)

    do k=1, n_degrees

      if ( (node_list%node(inode)%index(k) .ge. index_min(my_id+1)) .and. &
           (node_list%node(inode)%index(k) .le. index_max(my_id+1)) ) then
        ELM_is_local = .true.
        exit L_IV
      endif
    enddo

    if(node_list%node(inode)%constrained) then
	     do j = 1, 2
		 index1 = node_list%node(inode)%parents(j)
		  do k=1, n_degrees

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

if (DEBUG) write(*,'(i4,A,20i8)') my_id,' n_local_elms  : ',n_local_elms, element_list%n_elements

!--------------------- check distribution of boundary nodes 
mpi_distr_count=0
if (restart .and. freeboundary) then 
  do ib = 1, node_list%n_nodes	
    if ( node_list%node(ib)%boundary > 0 ) then
      do ik=1, n_degrees
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
     
a_mat%index_min => index_min
a_mat%index_max => index_max

a_mat%my_ind_min = a_mat%index_min(my_id + 1)
a_mat%my_ind_max = a_mat%index_max(my_id + 1)

return
end

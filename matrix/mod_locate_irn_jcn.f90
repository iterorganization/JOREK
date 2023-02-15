module mod_locate_irn_jcn
use mod_integer_types
implicit none
contains
subroutine locate_irn_jcn(index_node1,index_node2,index_min,index_max,ijA_position,a_mat)
use mod_integer_types
use data_structure, only: type_SP_MATRIX
!**************************************************************************
! subroutine finds the position in the global matrix of the index of      *
! node1 and node2 (this is the index per block)                           *
!                                                                         *
! search to be replaced by binary search                                  *
!**************************************************************************
integer               :: index_node1, index_node2, index_min, index_max, index1_local
integer(kind=int_all) :: ijA_position, i
logical               :: found_index
type(type_SP_MATRIX)  :: a_mat

found_index = .false.

index1_local = index_node1 - index_min + 1

!write(*,'(A,8i8)') ' LOCATE : ',index_node1,index_min,index_max,index1_local

do i=1,a_mat%ijA_size(index1_local)           ! replace by binary search?

  if (a_mat%irn_jcn(index1_local,i) .eq. index_node2) then
    ijA_position = a_mat%ijA_index(index1_local,i)
    found_index = .true.
    exit
  endif

enddo

if (.not.found_index) then

  write(*,*) ' FATAL locate_irn_jcn : index not found ',index_node1,index_node2

  do i=1,a_mat%ijA_size(index1_local)           ! replace by binary search?

    write(*,*) i, a_mat%irn_jcn(index1_local,i)
 stop
  enddo

endif

return
end subroutine locate_irn_jcn
end module mod_locate_irn_jcn

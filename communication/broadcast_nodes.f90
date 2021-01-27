subroutine Broadcast_nodes(my_id,node_list)
!----------------------------------------------------------
! subroutine to broadcast all the nodes in the node_list
!----------------------------------------------------------
use tr_module 
use data_structure
use mpi_mod
implicit none

type (type_node_list)    :: node_list

type (type_node)         :: anode
integer                  :: i, ierr, my_id, position, bufsize, IDBL_EXT, INT_EXT, ILOG_EXT
character, allocatable   :: buffer(:)

!  type type_node                                      ! type definition of a node (i.e. a vertex)
!    real*8     :: x(n_order+1,n_dim)                  ! x,y,z coordinates of points and additional nodal geometry
!    real*8     :: values(n_tor,n_order+1,n_var)
!    integer    :: index(n_order+1)                    ! the index in the main matrix
!    integer    :: boundary                            ! = 1, 2 or 3 for boundary nodes
!  endtype type_node                                   ! x(:,1) : position, x(:,2) : vector u, x(:,3) : vector v, x(4) : vector w

!  type type_node_list                                 ! type definition of a list of nodes
!    integer :: n_nodes                                ! the number of nodes in the list
!    type (type_node)     :: node(n_nodes_max)         ! an allocatable list of nodes
!  endtype type_node_list


call MPI_PACK_SIZE(1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,IDBL_EXT,ierr)
call MPI_PACK_SIZE(1,MPI_INTEGER,MPI_COMM_WORLD,INT_EXT,ierr)
call MPI_PACK_SIZE(1,MPI_LOGICAL,MPI_COMM_WORLD,ILOG_EXT,ierr)

call MPI_BCAST(node_list%n_nodes,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(node_list%n_dof,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

#ifdef fullmhd
bufsize = node_list%n_nodes * ((n_coord_tor*(n_order+1)*n_dim + 2*n_tor*(n_order+1)*n_var+2*(n_order+1)+2)*IDBL_EXT + (n_order+1 + 1+3+1)*INT_EXT + (2)*ILOG_EXT)
#elif altcs                          
bufsize = node_list%n_nodes * ((n_coord_tor*(n_order+1)*n_dim + 2*n_tor*(n_order+1)*n_var+2*(n_order+1)+2)*IDBL_EXT + (n_order+1 + 1+3+1)*INT_EXT + (2)*ILOG_EXT)
#else                                
bufsize = node_list%n_nodes * ((n_coord_tor*(n_order+1)*n_dim + 2*n_tor*(n_order+1)*n_var+2)*IDBL_EXT + (n_order+1 + 1+3+1)*INT_EXT + (2)*ILOG_EXT)
#endif


allocate(buffer(bufsize))
call tr_register_mem(bufsize,"bcastn_buffer")

if (my_id .eq. 0) then

  position = 0
  do i=1,node_list%n_nodes

    anode = node_list%node(i)

    call MPI_PACK(anode%x              ,n_coord_tor*(n_order+1)*n_dim      ,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
    call MPI_PACK(anode%values         ,n_tor*(n_order+1)*n_var,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
    call MPI_PACK(anode%deltas         ,n_tor*(n_order+1)*n_var,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
#ifdef fullmhd
    call MPI_PACK(anode%Fprof_eq       ,n_order+1,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
    call MPI_PACK(anode%psi_eq         ,n_order+1,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
#endif
    call MPI_PACK(anode%index          ,n_order+1,MPI_INTEGER,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
    call MPI_PACK(anode%boundary       ,1        ,MPI_INTEGER,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
    call MPI_PACK(anode%boundary_index ,1        ,MPI_INTEGER,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
    call MPI_PACK(anode%axis_node      ,1        ,MPI_LOGICAL,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
    call MPI_PACK(anode%constrained    ,1        ,MPI_LOGICAL,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
    call MPI_PACK(anode%parents(1:2)   ,2        ,MPI_INTEGER,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
    call MPI_PACK(anode%parent_elem    ,1        ,MPI_INTEGER,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
    call MPI_PACK(anode%ref_lambda     ,1        ,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
    call MPI_PACK(anode%ref_mu         ,1        ,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  
   
  enddo

endif

call MPI_BCAST(buffer,bufsize,MPI_PACKED,0,MPI_COMM_WORLD,ierr)

if (my_id .ne. 0) then

  position = 0
  do i=1,node_list%n_nodes

    call MPI_UNPACK(buffer,bufsize,position,anode%x              ,n_coord_tor*(n_order+1)*n_dim      ,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
    call MPI_UNPACK(buffer,bufsize,position,anode%values         ,n_tor*(n_order+1)*n_var,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
    call MPI_UNPACK(buffer,bufsize,position,anode%deltas         ,n_tor*(n_order+1)*n_var,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
#ifdef fullmhd
    call MPI_UNPACK(buffer,bufsize,position,anode%Fprof_eq       ,n_order+1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
    call MPI_UNPACK(buffer,bufsize,position,anode%psi_eq         ,n_order+1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
#endif
    call MPI_UNPACK(buffer,bufsize,position,anode%index          ,n_order+1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
    call MPI_UNPACK(buffer,bufsize,position,anode%boundary       ,1        ,MPI_INTEGER,MPI_COMM_WORLD,ierr)
    call MPI_UNPACK(buffer,bufsize,position,anode%boundary_index ,1        ,MPI_INTEGER,MPI_COMM_WORLD,ierr)
    call MPI_UNPACK(buffer,bufsize,position,anode%axis_node      ,1        ,MPI_LOGICAL,MPI_COMM_WORLD,ierr)
    call MPI_UNPACK(buffer,bufsize,position,anode%constrained    ,1        ,MPI_LOGICAL,MPI_COMM_WORLD,ierr)
    call MPI_UNPACK(buffer,bufsize,position,anode%parents(1:2)   ,2        ,MPI_INTEGER,MPI_COMM_WORLD,ierr)
    call MPI_UNPACK(buffer,bufsize,position,anode%parent_elem    ,1        ,MPI_INTEGER,MPI_COMM_WORLD,ierr)
    call MPI_UNPACK(buffer,bufsize,position,anode%ref_lambda     ,1        ,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
    call MPI_UNPACK(buffer,bufsize,position,anode%ref_mu         ,1        ,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
    node_list%node(i) = anode

  enddo

endif

call tr_unregister_mem(bufsize,"bcastn_buffer")
deallocate(buffer)

return
end subroutine Broadcast_nodes

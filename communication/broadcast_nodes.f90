subroutine Broadcast_nodes(my_id,node_list)
!----------------------------------------------------------
! subroutine to broadcast all the nodes in the node_list
!----------------------------------------------------------
use tr_module 
use data_structure
use mpi_mod
implicit none

integer, intent(in)                        :: my_id
type (type_node_list), intent(inout)       :: node_list

type (type_node)                           :: anode
integer                                    :: i, ierr, position, bufsize, IDBL_EXT, INT_EXT, ILOG_EXT
character, allocatable                     :: buffer(:)

call MPI_PACK_SIZE(1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,IDBL_EXT,ierr)
call MPI_PACK_SIZE(1,MPI_INTEGER,MPI_COMM_WORLD,INT_EXT,ierr)
call MPI_PACK_SIZE(1,MPI_LOGICAL,MPI_COMM_WORLD,ILOG_EXT,ierr)

call MPI_BCAST(node_list%n_values,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(node_list%n_nodes,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(node_list%n_dof,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

#ifdef STELLARATOR_MODEL
#ifndef USE_DOMM
! Model180 or Model183 with USE_EXT_FIELD: include j_field, b_field, b_vac_field
! For model183: r_tor_eq(n_degrees) + b_field(n_coord_tor*n_degrees*(n_dim+1)) + b_vac_field(n_coord_tor*n_degrees*(n_dim+1)) + j_source(n_tor*n_degrees)
! Buffer formula: n_coord_tor*n_degrees*(n_dim + 2*(n_dim+1)) + ... accounts for x(n_dim) + b_vac_field + ONE of (b_field or chi_correction)
! For model183 with USE_EXT_FIELD we now also pack b_field, so add +1*(n_dim+1) = 3*(n_dim+1) total
bufsize = node_list%n_nodes * ((n_coord_tor*n_degrees*(n_dim+3*(n_dim+1)+1) + 2*n_tor*n_degrees*node_list%n_values + n_tor*n_degrees + 2 + 2*n_degrees)*IDBL_EXT + (n_degrees + 1+3+1+1)*INT_EXT + (2)*ILOG_EXT)
#else
! USE_DOMM: model180 packs j_field+b_field+b_vac_field = 3*(n_dim+1), model183 only b_vac_field = 1*(n_dim+1)
! Use 3*(n_dim+1) to cover worst case (model180). Slight overallocation for model183 is harmless.
bufsize = node_list%n_nodes * ((n_coord_tor*n_degrees*(n_dim+3*(n_dim+1)) + 2*n_tor*n_degrees*node_list%n_values + n_tor*n_degrees + 2 + 2*n_degrees)*IDBL_EXT + (n_degrees + 1+3+1+1)*INT_EXT + (2)*ILOG_EXT)
#endif
#elif fullmhd
bufsize = node_list%n_nodes * ((n_coord_tor*n_degrees*n_dim + 2*n_tor*n_degrees*node_list%n_values+2*n_degrees+2)*IDBL_EXT + (n_degrees +1+3+1+1)*INT_EXT + (2)*ILOG_EXT)
#elif altcs                          
bufsize = node_list%n_nodes * ((n_coord_tor*n_degrees*n_dim + 2*n_tor*n_degrees*node_list%n_values+2*n_degrees+2)*IDBL_EXT + (n_degrees +1+3+1+1)*INT_EXT + (2)*ILOG_EXT)
#else                                
bufsize = node_list%n_nodes * ((n_coord_tor*n_degrees*n_dim + 2*n_tor*n_degrees*node_list%n_values+2)*IDBL_EXT + (n_degrees + 1+3+1+1)*INT_EXT + (2)*ILOG_EXT)
#endif

call init_node(anode, node_list%n_values)
allocate(buffer(bufsize))
call tr_register_mem(bufsize,"bcastn_buffer")

if (my_id .eq. 0) then

  position = 0

  do i=1,node_list%n_nodes

    call make_deep_copy_node(node_list%node(i), anode)

    call MPI_PACK(anode%x              ,n_coord_tor*n_degrees*n_dim      ,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
    call MPI_PACK(anode%values         ,n_tor*n_degrees*node_list%n_values,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
    call MPI_PACK(anode%deltas         ,n_tor*n_degrees*node_list%n_values,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
#ifdef STELLARATOR_MODEL
    call MPI_PACK(anode%r_tor_eq       ,n_degrees,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
#if JOREK_MODEL == 180
    call MPI_PACK(anode%pressure       ,n_degrees,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
    call MPI_PACK(anode%j_field        ,n_coord_tor*n_degrees*(n_dim+1),MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
    call MPI_PACK(anode%b_field        ,n_coord_tor*n_degrees*(n_dim+1),MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
#endif
#ifndef USE_DOMM
#ifdef USE_EXT_FIELD
    call MPI_PACK(anode%b_vac_field    ,n_coord_tor*n_degrees*(n_dim+1),MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
#else
    call MPI_PACK(anode%chi_correction ,n_coord_tor*n_degrees          ,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
#endif
#endif
    call MPI_PACK(anode%j_source       ,n_tor*n_degrees,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
#elif fullmhd
    call MPI_PACK(anode%Fprof_eq       ,n_degrees,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
    call MPI_PACK(anode%psi_eq         ,n_degrees,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
#endif
    call MPI_PACK(anode%index          ,n_degrees,MPI_INTEGER,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
    call MPI_PACK(anode%boundary       ,1        ,MPI_INTEGER,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
    call MPI_PACK(anode%boundary_index ,1        ,MPI_INTEGER,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
    call MPI_PACK(anode%axis_node      ,1        ,MPI_LOGICAL,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
    call MPI_PACK(anode%axis_dof       ,1        ,MPI_INTEGER,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
    call MPI_PACK(anode%constrained    ,1        ,MPI_LOGICAL,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
    call MPI_PACK(anode%parents(1:2)   ,2        ,MPI_INTEGER,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
    call MPI_PACK(anode%parent_elem    ,1        ,MPI_INTEGER,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
    call MPI_PACK(anode%ref_lambda     ,1        ,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
    call MPI_PACK(anode%ref_mu         ,1        ,MPI_DOUBLE_PRECISION,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
  
   
  enddo

  if (position > bufsize) then
    write(*,*) 'FATAL broadcast_nodes: MPI pack buffer overflow: position=', position, ', bufsize=', bufsize
    write(*,*) 'Bufsize formula is too small - update broadcast_nodes.f90'
    stop
  end if

endif

call MPI_BCAST(buffer,bufsize,MPI_PACKED,0,MPI_COMM_WORLD,ierr)

if (my_id .ne. 0) then

    if (.not. allocated(node_list%node)) call init_node_list(node_list, node_list%n_nodes, node_list%n_dof, node_list%n_values)

  position = 0
  do i=1,node_list%n_nodes

    call MPI_UNPACK(buffer,bufsize,position,anode%x              ,n_coord_tor*n_degrees*n_dim      ,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
    call MPI_UNPACK(buffer,bufsize,position,anode%values         ,n_tor*n_degrees*node_list%n_values,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
    call MPI_UNPACK(buffer,bufsize,position,anode%deltas         ,n_tor*n_degrees*node_list%n_values,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
#ifdef STELLARATOR_MODEL
    call MPI_UNPACK(buffer,bufsize,position,anode%r_tor_eq       ,n_degrees                        ,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
#if JOREK_MODEL == 180
    call MPI_UNPACK(buffer,bufsize,position,anode%pressure       ,n_degrees                        ,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
    call MPI_UNPACK(buffer,bufsize,position,anode%j_field        ,n_coord_tor*n_degrees*(n_dim+1)  ,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
    call MPI_UNPACK(buffer,bufsize,position,anode%b_field        ,n_coord_tor*n_degrees*(n_dim+1)  ,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
#endif
#ifndef USE_DOMM
#ifdef USE_EXT_FIELD
    call MPI_UNPACK(buffer,bufsize,position,anode%b_vac_field    ,n_coord_tor*n_degrees*(n_dim+1)  ,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
#else
    call MPI_UNPACK(buffer,bufsize,position,anode%chi_correction ,n_coord_tor*n_degrees            ,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
#endif
#endif
    call MPI_UNPACK(buffer,bufsize,position,anode%j_source       ,n_tor*n_degrees                  ,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
#elif fullmhd
    call MPI_UNPACK(buffer,bufsize,position,anode%Fprof_eq       ,n_degrees,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
    call MPI_UNPACK(buffer,bufsize,position,anode%psi_eq         ,n_degrees,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
#endif
    call MPI_UNPACK(buffer,bufsize,position,anode%index          ,n_degrees,MPI_INTEGER,MPI_COMM_WORLD,ierr)
    call MPI_UNPACK(buffer,bufsize,position,anode%boundary       ,1        ,MPI_INTEGER,MPI_COMM_WORLD,ierr)
    call MPI_UNPACK(buffer,bufsize,position,anode%boundary_index ,1        ,MPI_INTEGER,MPI_COMM_WORLD,ierr)
    call MPI_UNPACK(buffer,bufsize,position,anode%axis_node      ,1        ,MPI_LOGICAL,MPI_COMM_WORLD,ierr)
    call MPI_UNPACK(buffer,bufsize,position,anode%axis_dof       ,1        ,MPI_INTEGER,MPI_COMM_WORLD,ierr)
    call MPI_UNPACK(buffer,bufsize,position,anode%constrained    ,1        ,MPI_LOGICAL,MPI_COMM_WORLD,ierr)
    call MPI_UNPACK(buffer,bufsize,position,anode%parents(1:2)   ,2        ,MPI_INTEGER,MPI_COMM_WORLD,ierr)
    call MPI_UNPACK(buffer,bufsize,position,anode%parent_elem    ,1        ,MPI_INTEGER,MPI_COMM_WORLD,ierr)
    call MPI_UNPACK(buffer,bufsize,position,anode%ref_lambda     ,1        ,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
    call MPI_UNPACK(buffer,bufsize,position,anode%ref_mu         ,1        ,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)

    call make_deep_copy_node(anode, node_list%node(i))

  enddo

endif

call tr_unregister_mem(bufsize,"bcastn_buffer")
call dealloc_node(anode)
deallocate(buffer)

return
end subroutine Broadcast_nodes

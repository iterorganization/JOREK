module centralization_mod

implicit none

contains

  !> Collect a distributed matrix on a single MPI process (the master task of the MPI group)
  subroutine centralization_harmonic(my_id, my_id_n, n_cpu_n, MPI_COMM_N)

  use data_structure 
  use global_distributed_matrix
  use mumps_module  
  use mpi_mod
  use phys_module, only: centralize_harm_mat, use_strumpack
  use mod_integer_types
  
  implicit none
  
  ! --- Routine parameters
  integer, intent(in)                :: my_id, my_id_n, n_cpu_n, MPI_COMM_N
  
  ! --- Local variables
  integer(kind=int_all), allocatable :: nz_array(:), disp_array(:)
  integer(kind=int_all), parameter   :: Int1=1
  integer(kind=int_all)              :: nz_total, i
  integer                            :: i_cpu, ierr
  
 
  if (.not.centralize_harm_mat .and. use_strumpack) then

    mumps_par%nz  = nz_harm

    if (associated(mumps_par%A))   call tr_deallocatep(mumps_par%A,  "dh_mumps_par%A",  CAT_DMATRIX)
    if (associated(mumps_par%irn)) call tr_deallocatep(mumps_par%irn,"dh_mumps_par%irn",CAT_DMATRIX)
    if (associated(mumps_par%jcn)) call tr_deallocatep(mumps_par%jcn,"dh_mumps_par%jcn",CAT_DMATRIX)
    call tr_allocatep(mumps_par%A,  Int1,mumps_par%nz,"dh_mumps_par%A",  CAT_DMATRIX)
    call tr_allocatep(mumps_par%irn,Int1,mumps_par%nz,"dh_mumps_par%irn",CAT_DMATRIX)
    call tr_allocatep(mumps_par%jcn,Int1,mumps_par%nz,"dh_mumps_par%jcn",CAT_DMATRIX)
    
    do i = Int1, mumps_par%nz
      mumps_par%A(i)   = A_harm(i)
      mumps_par%irn(i)   = irn_harm(i)
      mumps_par%jcn(i)   = jcn_harm(i)
    enddo
    
  else

    allocate(nz_array  (n_cpu_n))
    allocate(disp_array(n_cpu_n))
      
    ! --- Determine matrix dimension n and number of nonzeros nz
    if (n_cpu_n > 1) then
      call MPI_GATHER(nz_harm, 1, MPI_INTEGER_ALL, nz_array, 1, MPI_INTEGER_ALL, 0, MPI_COMM_N, ierr) 
      disp_array = 0 
      nz_total   = 0
      if (my_id_n .eq. 0) then
        do i_cpu = 2, n_cpu_n  
          disp_array(i_cpu) = disp_array(i_cpu-1) + nz_array(i_cpu-1)
        enddo  
        nz_total = disp_array(n_cpu_n) + nz_array(n_cpu_n) 
      endif 
      mumps_par%nz  = nz_total   
    else
      mumps_par%nz  = nz_harm
    endif
  
    ! --- Allocate arrays for centralized matrix
    ! --- Centralize matrix (if it was not distributed, copy it into the right data structure)
    if (associated(mumps_par%A))   call tr_deallocatep(mumps_par%A,  "dh_mumps_par%A",  CAT_DMATRIX)
    if (associated(mumps_par%irn)) call tr_deallocatep(mumps_par%irn,"dh_mumps_par%irn",CAT_DMATRIX)
    if (associated(mumps_par%jcn)) call tr_deallocatep(mumps_par%jcn,"dh_mumps_par%jcn",CAT_DMATRIX)
    call tr_allocatep(mumps_par%A,  Int1,mumps_par%nz,"dh_mumps_par%A",  CAT_DMATRIX)
    call tr_allocatep(mumps_par%irn,Int1,mumps_par%nz,"dh_mumps_par%irn",CAT_DMATRIX)
    call tr_allocatep(mumps_par%jcn,Int1,mumps_par%nz,"dh_mumps_par%jcn",CAT_DMATRIX)

    if (n_cpu_n > 1) then
      call MPI_GATHERV(A_harm, nz_harm, MPI_DOUBLE_PRECISION, mumps_par%A, nz_array, disp_array,&
        MPI_DOUBLE_PRECISION, 0, MPI_COMM_N, ierr)
  
      call MPI_GATHERV(irn_harm, nz_harm, MPI_INTEGER_ALL, mumps_par%irn, nz_array, disp_array,&
        MPI_INTEGER_ALL, 0, MPI_COMM_N, ierr)
  
      call MPI_GATHERV(jcn_harm, nz_harm, MPI_INTEGER_ALL, mumps_par%jcn, nz_array, disp_array,&
        MPI_INTEGER_ALL, 0, MPI_COMM_N, ierr)
    else 

     do i = Int1, mumps_par%nz
       mumps_par%A(i)   = A_harm(i) 
       mumps_par%irn(i)   = irn_harm(i)
       mumps_par%jcn(i)   = jcn_harm(i)
     enddo

    endif

    if ( allocated(nz_array) )   deallocate(nz_array)
    if ( allocated(disp_array) ) deallocate(disp_array)
   
  endif ! centralize_harm_mat

  mumps_par%n   = ndof_harm
  
  if (associated(mumps_par%rhs)) call tr_deallocatep(mumps_par%rhs,"dh_mumps_par%rhs",CAT_DMATRIX)
  call tr_allocatep(mumps_par%rhs,Int1,mumps_par%n, "dh_mumps_par%rhs",CAT_DMATRIX)
  
  do i = 1, mumps_par%n
    mumps_par%rhs(i) = rhs_harm(i)
  enddo  

  if (allocated(A_harm)) call tr_deallocate(A_harm,"A_harm",CAT_DMATRIX)
  
  end subroutine centralization_harmonic

end module centralization_mod

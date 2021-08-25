module gmres_setup
  implicit none
contains
!> Setup MPI communicators for gmres solving in JOREK.
!> Example (2 harmonic, 6 mpi proc, number is my_id_*):
!>      MPI_COMM_WORLD MPI_COMM_N     MPI_COMM_TRANS MPI_COMM_MASTER
!> harm +-----------+  +-----------+  +-----------+  +-----------+
!>  1   | 3 | 4 | 5 |  | 0 | 1 | 2 |  | 1 | 1 | 1 |  |   |   |   |
!>      +-----------+  +-----------+  +-----------+  +-----------+
!>  0   | 0 | 1 | 2 |  | 0 | 1 | 2 |  | 0 | 0 | 0 |  | 0 | 1 | 2 |
!>      +-----------+  +-----------+  +-----------+  +-----------+
subroutine gmres_setup_jorek(my_id, n_cpu, i_tor, my_id_n, n_cpu_n, my_id_trans, n_cpu_trans, my_id_master,  &
                             MPI_COMM_N, MPI_COMM_TRANS, MPI_COMM_MASTER, MPI_GROUP_MASTER, MPI_GROUP_WORLD, &
                             my_family_id)
  use mod_parameters, only: n_tor
  use tr_module
  use mpi_mod
  implicit none
  integer, intent(in) :: my_id, n_cpu
  integer, allocatable, intent(inout) :: i_tor(:)
  integer, intent(out) :: my_id_n, n_cpu_n
  integer, intent(out) :: my_id_trans, n_cpu_trans
  integer, intent(out) :: my_id_master
  integer, intent(out) :: MPI_COMM_N !< group for each harmonic
  integer, intent(out) :: MPI_COMM_TRANS !< transversal groups (i.e. every 1st of MPI_COMM_N, every 2nd of MPI_COMM_N etc)
  integer, intent(out) :: MPI_COMM_MASTER !< Every first of MPI_COMM_N
  integer, intent(out) :: MPI_GROUP_MASTER !< subset of MPI_COMM_WORLD corresponding to MPI_COMM_MASTER
  integer, intent(out) :: MPI_GROUP_WORLD
  integer, intent(out) :: my_family_id
  integer :: N_masters, i_rank(n_tor), m_cpu
  integer :: i, ierr

  N_masters = (n_tor+1)/2
  if (MOD(n_cpu, N_masters) == 0) then
    M_cpu = n_cpu / (N_masters)
  else
    M_cpu = (n_cpu - MOD(n_cpu, N_masters))/N_masters +1
  end if

  if (allocated(i_tor)) call tr_deallocate(i_tor,"i_tor",CAT_UNKNOWN)
  call tr_allocate(i_tor,1,n_cpu,"i_tor",CAT_UNKNOWN)
  do i = 1, n_cpu 
    i_tor(i) =  MOD(i-1, M_cpu)+1
  end do
  call MPI_COMM_SPLIT(MPI_COMM_WORLD,i_tor(my_id+1),my_id,MPI_COMM_TRANS,ierr)

  do i=1,n_cpu
    i_tor(i) = ((i-1) - MOD(i-1, M_cpu))/ M_cpu  + 1
  enddo

  my_family_id = i_tor(my_id+1)

  call MPI_COMM_SPLIT(MPI_COMM_WORLD,i_tor(my_id+1),my_id,MPI_COMM_N,ierr)
  
  do i=1,N_masters
    i_rank(i) = (i-1) * M_cpu
  enddo

  call MPI_COMM_GROUP(MPI_COMM_WORLD,MPI_GROUP_WORLD,ierr)
  call MPI_GROUP_INCL(MPI_GROUP_WORLD,N_masters,i_rank,MPI_GROUP_MASTER,ierr)

  call MPI_COMM_CREATE(MPI_COMM_WORLD,MPI_GROUP_MASTER,MPI_COMM_MASTER,ierr)

  call MPI_COMM_RANK(MPI_COMM_N, my_id_n, ierr) ! id of this cpu in local comm
  call MPI_COMM_SIZE(MPI_COMM_N, n_cpu_n, ierr) ! number of local comm cpu
  call MPI_COMM_RANK(MPI_COMM_TRANS, my_id_trans, ierr) ! id of this proc in transverse comm
  call MPI_COMM_SIZE(MPI_COMM_TRANS, n_cpu_trans, ierr) ! num proc in transverse comm

  write(*,*) 'gmres_setup : m_cpu, n_cpu_n ',m_cpu, n_cpu_n

end subroutine gmres_setup_jorek
end module gmres_setup

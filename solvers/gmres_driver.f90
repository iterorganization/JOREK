!> Driver for the reverse communication GMRES routine from dPackgmres (CERFACS)
subroutine gmres_driver(my_id,my_id_n,i_tor,n_tor,MPI_COMM_N,MPI_COMM_MASTER,iter_gmres)

use tr_module 
use mumps_module
use global_distributed_matrix
use mpi_mod
use mod_clock
use phys_module,  only : gmres_4, gmres_tol, gmres_m, index_now
use mod_integer_types

implicit none
#include "r3_info.h"

interface 
   subroutine gmres_matrix_vector(size_x,x,size_y,y,my_id,my_id_n, i_tor, MPI_COMM_MASTER)      
     use mod_integer_types
     integer               :: i_tor(:), MPI_COMM_MASTER
     integer(kind=int_all) :: size_x,size_y
     real*8                :: x(size_x), y(size_y)
     integer               :: my_id, my_id_n
   end subroutine gmres_matrix_vector
end interface
integer               :: i_tor(:), j, my_id, my_id_n, my_id_master, MPI_COMM_N, MPI_COMM_MASTER
integer               :: revcom, iter_gmres, n_tor
integer               :: icntl(8), info(3)
integer(kind=int_all) :: i, m, colx, coly, colz, nbscal, lwork, n_dof, irc(5)

integer(kind=int_all), parameter   :: Int1=1

integer :: matvec, precondLeft, precondRight, dotProd, ierr
real*8  :: cntl(5), rinfo(2), sum, err, Bnorm, Xnorm
real*8, allocatable :: work(:), work_ndof(:), work_ndof2(:)
type(clcktype)           :: t0, t1
real*8                   :: tsecond
REAL*8, ALLOCATABLE      :: rhs_tmp(:)
real*8 ::ZERO, ONE
CHARACTER(LEN=128) :: fname

parameter (ZERO = 0.0d0, ONE = 1.0d0)
parameter (matvec=1, precondLeft=2, precondRight=3, dotProd=4)

!write(*,*) ' GMRES DRIVER : ',my_id,my_id_n
call r3_info_begin (r3_info_index_0, 'gmres_driver')  ! timing
call clck_time(t0)
call init_dgmres(icntl,cntl)

icntl(3) = 6            ! output unit
icntl(7) = iter_gmres   ! Maximum number of iterations
icntl(4) = 1            ! preconditioner (1) = left preconditioner
icntl(5) = 3            ! orthogonalization scheme
icntl(6) = 1            ! initial guess  (1) = user supplied guess
icntl(8) = 1            ! residual calculation strategy at restart

cntl(1) = gmres_tol         ! stopping tolerance
cntl(2) = 1.d0
cntl(3) = 1.d0
cntl(4) = gmres_4          ! 1.d0
cntl(5) = 1.d0

m = gmres_m 

n_dof = ndof_glob

lwork = m*m + m*(n_dof+5) + 6*n_dof + m + 1

call tr_allocate(work,Int1,lwork,"work",CAT_GMRES)
call tr_allocate(work_ndof,Int1,n_dof,"work_ndof",CAT_GMRES)
call tr_allocate(work_ndof2,Int1,n_dof,"work_ndof2",CAT_GMRES)

work(1:n_dof)         = deltas(1:n_dof)                     ! the initial guess
work(n_dof+1:2*n_dof) = RHS_glob(1:n_dof)                   ! the right hand side

work_ndof(Int1:n_dof) = work(Int1:n_dof)
work_ndof2(Int1:n_dof) = work(2*n_dof+Int1:3*n_dof)
call gmres_matrix_vector(n_dof,work_ndof,n_dof,work_ndof2,my_id,my_id_n, i_tor, MPI_COMM_MASTER)
work(Int1:n_dof) = work_ndof(Int1:n_dof)
work(2*n_dof+Int1:3*n_dof) = work_ndof2(Int1:n_dof)

if (my_id .eq. 0) then
  sum = 0.d0
  err = -1.d20
  Bnorm = 0.d0
  Xnorm = 0.d0
  do i=1,n_dof
    sum = sum + (work(2*n_dof+i)-work(n_dof+i))**2
    err = max(err,abs(work(2*n_dof+i)-work(n_dof+i)))
    Bnorm = Bnorm + RHS_glob(i)**2
    Xnorm = Xnorm + deltas(i)**2
  enddo

  write(*,'(A,4e16.8)') ' residu test before : ',sqrt(sum),err,sqrt(Bnorm),sqrt(Xnorm)
endif
if (my_id .eq. 0) then
   write(fname,'(A,I6.6,A1,I6.6)')"product_before",index_now
   call tr_vdump(fname, work(2*n_dof+1:3*n_dof), n_dof)
   write(fname,'(A,I6.6,A1,I6.6)')"X_before",index_now
   call tr_vdump(fname, RHS_glob, n_dof)
end if

!*****************************************
!** Reverse communication implementation
!*****************************************

10     call MPI_barrier(MPI_COMM_WORLD,ierr)

       if (my_id .eq. 0) then
         call drive_dgmres(n_dof,n_dof,m,lwork,work,irc,icntl,cntl,info,rinfo)
       endif

       call MPI_BCAST(irc,5,MPI_INTEGER_ALL,0,MPI_COMM_WORLD,ierr)
       revcom = irc(1)
       colx   = irc(2)
       coly   = irc(3)
       colz   = irc(4)
       nbscal = irc(5)

       if (revcom.eq.matvec) then                  ! perform the matrix vector product
                                                   ! work(colz) <-- A * work(colx)
        
         work_ndof(Int1:n_dof) = work(colx:colx+n_dof-Int1)
         work_ndof2(Int1:n_dof) = work(colz:colz+n_dof-Int1)
         call gmres_matrix_vector(n_dof,work_ndof,n_dof,work_ndof2,my_id,my_id_n, i_tor, MPI_COMM_MASTER)
         work(colx:colx+n_dof-Int1) = work_ndof(Int1:n_dof)
         work(colz:colz+n_dof-Int1) = work_ndof2(Int1:n_dof)
         goto 10

       else if (revcom.eq.precondLeft) then        ! perform the left preconditioning
                                                   ! work(colz) <-- M^{-1} * work(colx)
         call gmres_precondition(work(colx),work(colz),i_tor,my_id,my_id_n,MPI_COMM_MASTER,MPI_COMM_N)
         goto 10

       else if (revcom.eq.precondRight) then       ! perform the right preconditioning

         if (my_id .eq. 0) call dcopy(n_dof,work(colx),Int1,work(colz),Int1)

         goto 10

       else if (revcom.eq.dotProd) then            ! perform the scalar product
                                                   ! work(colz) <-- work(colx) work(coly)

         if (my_id .eq. 0) call dgemv('C',n_dof,nbscal,ONE, work(colx),n_dof,work(coly),Int1,ZERO,work(colz),Int1)
 
         goto 10

       endif

!******************************** end of GMRES reverse communication

if (my_id .eq. 0) write(*,*) 'check work : ',n_dof,maxval(abs(work(1:n_dof)))

if (my_id .eq. 0) deltas(1:n_dof) = work(1:n_dof)

work_ndof(Int1:n_dof) = deltas(Int1:n_dof)
work_ndof2(Int1:n_dof) = work(n_dof+Int1:2*n_dof)
call gmres_matrix_vector(n_dof,work_ndof,n_dof,work_ndof2,my_id,my_id_n, i_tor, MPI_COMM_MASTER)
deltas(Int1:n_dof) = work_ndof(Int1:n_dof)
work(n_dof+Int1:2*n_dof) = work_ndof2(Int1:n_dof)

if (my_id .eq. 0) then
  sum = 0.d0
  err = -1.d20
  Bnorm = 0.d0
  Xnorm = 0.d0
  do i=1,n_dof
    sum   = sum + (work(n_dof+i)-RHS_glob(i))**2
    err   = max(err,abs(work(n_dof+i)-RHS_glob(i)))
    Bnorm = Bnorm + RHS_glob(i)**2
    Xnorm = Xnorm + deltas(i)**2
  enddo
  write(*,'(A,4e16.8)') ' residu test after : ',sqrt(sum),err,sqrt(Bnorm),sqrt(Xnorm)
  call tr_locvnorms("gmres_rhs",rhs_glob,n_dof)
endif

iter_gmres = info(2) ! Actual number of iterations
call tr_debug_writei("nbiter_gmres", iter_gmres)
call MPI_BCAST(iter_gmres,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

call tr_deallocate(work,"work",CAT_GMRES)
call tr_deallocate(work_ndof,"work_ndof",CAT_GMRES)
call tr_deallocate(work_ndof2,"work_ndof2",CAT_GMRES)

call clck_time(t1)
call clck_ldiff(t0,t1,tsecond)
!if (my_id .eq. 0)  write(*,FMT_TIMING) my_id, '## Elapsed time gmres :', tsecond

call r3_info_end (r3_info_index_0)  ! timing
return
end

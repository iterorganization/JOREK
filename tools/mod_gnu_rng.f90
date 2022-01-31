! the mod_gnu_rng modules contains procedures
! helping the generation of random numbers
module mod_gnu_rng
implicit none

private
public :: gnu_rng_interval
public :: gnu_rng_array_norep
public :: set_seed_sys_time

! gnu_rng_interval: template interface for the
! gnu_rng_interval implementation
interface gnu_rng_interval
  module procedure gnu_rng_int_interval_single
  module procedure gnu_rng_int_interval_1d
  module procedure gnu_rng_int_interval_2d
  module procedure gnu_rng_real8_interval_single 
  module procedure gnu_rng_real8_interval_1d
  module procedure gnu_rng_real8_interval_2d
  module procedure gnu_rng_interval_1d_val_1d
end interface gnu_rng_interval

! gnu_rng_array_norep: template function for 
! procedures implementing the generation of
! random arrays without repetitions
interface gnu_rng_array_norep
  module procedure gnu_rng_int_array_1d_norep
end interface

interface set_seed_sys_time
  module procedure set_seed_sys_time_novar_1d
  module procedure set_seed_sys_time_1var_1d
  module procedure set_seed_sys_time_2vars_1d
end interface set_seed_sys_time

contains

! Methods for settings RNGs seeds ----------------------------
!> set a different seed based on the system time
!>   rn_size:  (integer) size of the random number vector
!>   interval: (integer)(2) interval for random integer generation
!> outputs
!>   seed_out: (integer*8)(optional) the set rng seed
subroutine set_seed_sys_time_novar_1d(interval)
  implicit none
  integer,dimension(2),intent(in) :: interval
  integer                         :: n_seeds,time
  integer,dimension(:),allocatable :: seed
  call random_seed(size=n_seeds)
  allocate(seed(n_seeds))
  call gnu_rng_interval(n_seeds,interval,seed)
  call system_clock(time)
  seed = seed*time
  call random_seed(put=seed)
  deallocate(seed)
end subroutine set_seed_sys_time_novar_1d

!> set a different seed based on the system time
!> for one variable (mpi rand or thread id)
!> inputs:
!>   rn_size:  (integer) size of the random number vector
!>   n_var:    (integer) number of variables
!>   interval: (integer)(2) interval for random integer generation
!>   var:      (integer) variable for setting the seed
!> outputs:
!>   n_seed_out: (integer) size of the random seed
!>   seed_out:  (integer)(n_seed_out)(optional) the set rng seed
subroutine set_seed_sys_time_1var_1d(interval,var)
  implicit none
  integer,intent(in)               :: var
  integer,dimension(2),intent(in)  :: interval
  integer                          :: n_seeds,time
  integer,dimension(:),allocatable :: seed
  call random_seed(size=n_seeds) 
  allocate(seed(n_seeds))
  call gnu_rng_interval(n_seeds,interval,seed)
  call system_clock(time)
  seed = sign(1,var)*(abs(var)+1)*seed*time
  call random_seed(put=seed)
  deallocate(seed)
end subroutine set_seed_sys_time_1var_1d

!> set a different seed based on the system time
!> for two variables (e.g. mpi rank and thread)
!> inputs:
!>   rn_size:  (integer) size of the random number vector
!>   n_var_1:  (integer) number of variables 1
!>   n_var_2   (integer) number of variables 2
!>   interval: (integer)(2) interval for random integer generation
!>   var1:     (integer) first variable setting rng
!>   var2:     (integer) second variable setting rng
!> outputs:
!>   n_seed_out: (integer)(optional) size of the random seed
!>   seed_out:   (integer)(n_seed_out)(optional) the set rng seed
subroutine set_seed_sys_time_2vars_1d(interval,var1,var2)
  implicit none
  integer,dimension(2),intent(in) :: interval
  integer,intent(in)              :: var1,var2
  integer                         :: n_seeds,time
  integer,dimension(:),allocatable :: seed
  call random_seed(size=n_seeds) 
  allocate(seed(n_seeds))
  call gnu_rng_interval(n_seeds,interval,seed)
  call system_clock(time)
  seed = sign(1,var2)*sign(1,var1)*(abs(var1)+1)*(abs(var2)+1)*seed*time
  call random_seed(put=seed)
  deallocate(seed)
end subroutine set_seed_sys_time_2vars_1d


! Uniform RNG ------------------------------------------------

! gnu_rng_int_interval_single computes a
! random integer from uniform distribution
! inputs:
!   interval: (2)(integer) minimum and maximum values
! outputs:
!   id: (index) index
subroutine gnu_rng_int_interval_single(interval,id)
  implicit none

  ! inputs
  integer,dimension(2),intent(in) :: interval
  ! outputs
  integer,intent(out) :: id
  ! variables
  real*8 :: rnd

   ! compute random integer 
   call random_number(rnd)
   id = floor(interval(1)+(interval(2)-interval(1)+1)*rnd)

end subroutine gnu_rng_int_interval_single

! gnu_rng_int_interval_1d computes an array of
! random integers from a uniform distribution
! inputs:
!   N:        (intrger) length of the array
!   interval: (integer)(2) minimum and maximum values
! outputs:
!   ids: (index)(N) array of random integers
subroutine gnu_rng_int_interval_1d(N,interval,ids)
  implicit none

  ! inputs
  integer,intent(in) :: N
  integer,dimension(2),intent(in) :: interval
  ! outputs
  integer,dimension(N),intent(out) :: ids
  ! variables
  real*8,dimension(N) :: rnds

  ! compute random integer
  call random_number(rnds)
  ids = floor(interval(1)+(interval(2)-interval(1)+1)*rnds)
end subroutine gnu_rng_int_interval_1d

! gnu_rng_int_interval_2d computes a 2D-array of
! random integers from a uniform distribution
! inputs:
!   N_rows:   (intrger) N# of array rows
!   N_cols:   (integer) N# of array columns
!   interval: (integer)(2) minimum and maximum values
! outputs:
!   ids: (index)(N_rows,N_cols) array of random integers
subroutine gnu_rng_int_interval_2d(N_rows,N_cols,interval,ids)
  implicit none

  ! inputs
  integer,intent(in) :: N_rows,N_cols
  integer,dimension(2),intent(in) :: interval
  ! outputs
  integer,dimension(N_rows,N_cols),intent(out) :: ids
  ! variables
  real*8,dimension(N_rows,N_cols) :: rnds

  ! compute random integer
  call random_number(rnds)
  ids = floor(interval(1)+(interval(2)-interval(1)+1)*rnds)
end subroutine gnu_rng_int_interval_2d

! gnu_rng_rng_interval_single computes a
! random double from uniform distribution
! inputs:
!   interval: (2)(real8) minimum and maximum values
! outputs:
!   rng_0d: (index) index
subroutine gnu_rng_real8_interval_single(interval,rng_0d)
  implicit none

  ! inputs
  real*8,dimension(2),intent(in) :: interval
  ! outputs
  real*8,intent(out) :: rng_0d

   ! compute random integer 
   call random_number(rng_0d)
   rng_0d = interval(1)+(interval(2)-interval(1))*rng_0d

end subroutine gnu_rng_real8_interval_single

! gnu_rng_interval_1d generates a 1D random number array
! using the gnu-fortran intrinsic function within a 
! predefined interval (uniform distribution)
! inputs:
!   N:            (int) length of the array
!   interval:     (2)(real8) interval within the value
!                            are generated
!   rng_array_1d: (N)(real8) 1D array of uniform random
!                            numbers within an interval
! outputs:
!   rng_array_1d: (N)(real8) 1D array of uniform random
!                            numbers within an interval
subroutine gnu_rng_real8_interval_1d(N,interval,rng_array_1d)
  implicit none

  ! inputs:
  integer,intent(in)                :: N
  real*8,dimension(2),intent(in)    :: interval
  ! inputs-outputs:
  real*8,dimension(N),intent(inout) :: rng_array_1d

  ! generate random number array
  call random_number(rng_array_1d)
  rng_array_1d = interval(1) + (interval(2)-&
  interval(1))*rng_array_1d

end subroutine gnu_rng_real8_interval_1d

! gnu_rng_interval_2d generates a 2D random number array
! using the gnu-fortran intrinsic function within a 
! predefined interval (uniform distribution)
! inputs:
!   N_rows:       (int) number of array rows
!   N_cols:       (int) number of array cols
!   interval:     (2)(real8) interval within the value
!                            are generated
!   rng_array_2d: (N_rows,N_cols)(real8) 2D array of 
!                            uniform random numbers 
!                            within an interval
! outputs:
!   rng_array_2d: (N_rows,N_cols)(real8) 2D array of 
!                            uniform random numbers 
!                            within an interval
subroutine gnu_rng_real8_interval_2d(N_rows,N_cols,interval,rng_array_2d)
  implicit none

  ! inputs:
  integer,intent(in)                            :: N_rows,N_cols
  real*8,dimension(2),intent(in)                :: interval
  ! inputs-outputs:
  real*8,dimension(N_rows,N_cols),intent(inout) :: rng_array_2d

  ! generate random number array
  call random_number(rng_array_2d)
  rng_array_2d = interval(1) + (interval(2)-&
  interval(1))*rng_array_2d

end subroutine gnu_rng_real8_interval_2d

! gnu_rng_interval_1d_val_1d generates a 1D random number array
! using the gnu-fortran intrinsic function within a 
! predefined set of intervals (uniform distribution)
! inputs:
!   N:               (integer) array length
!   interval_lowbnd: (N)(real8) interval lower bounds
!   interval_uppbnd: (N)(real8) interval upper bounds
!   rng_array_1d:    (N)(real8) 1D array of uniform random
!                               numbers within an intervals
! outputs:
!   rng_array_1d:    (N)(real8) 1D array of uniform random
!                               numbers within an intervals
subroutine gnu_rng_interval_1d_val_1d(N,interval_lowbnd,&
interval_uppbnd,rng_array_1d)
  implicit none

  ! inputs
  integer,intent(in)                :: N
  real*8,dimension(N),intent(in)    :: interval_lowbnd
  real*8,dimension(N),intent(in)    :: interval_uppbnd
  ! inputs-outputs:
  real*8,dimension(N),intent(inout) :: rng_array_1d

  ! generate random number array
  call random_number(rng_array_1d)
  rng_array_1d = interval_lowbnd + (interval_uppbnd-&
  interval_lowbnd)*rng_array_1d

end subroutine gnu_rng_interval_1d_val_1d

! Additional tools -------------------------------------------

! gnu_rng_int_array_1d_norep generates an array of random
! integer number without repetitions within an interval.
! inputs:
!   N:         (integer) number of elements
!   interval:  (2)(integer) minimum and maximum value
!   n_max:     (integer) maximum value
!   ierr:      (integer) error value
!   max_it_in: (integer)(optional) maximum number of iterations
! outputs:
!   ids:   (integer)(N) array of non repeated random integers
subroutine gnu_rng_int_array_1d_norep(N,interval,ids,ierr,max_it_in)
  implicit none

  ! inputs
  integer,intent(in) :: N
  integer,dimension(2),intent(in) :: interval
  integer,intent(in),optional :: max_it_in
  ! outputs
  integer,dimension(N),intent(out) :: ids
  ! inputs-outputs
  integer,intent(inout) :: ierr
  ! variables
  integer :: ii,jj,max_it
  real*8 :: rnd
  
  ! check for optional inputs
  max_it = 1000000
  if(present(max_it_in)) max_it = max_it_in 
  ! generate sequence of random numbers
  call gnu_rng_int_interval_1d(N,interval,ids)
  ! try to correct for repeated ids
  do jj=1,N
      ii=1
      do while((count(ids==ids(jj)).gt.1).and.(ii.le.max_it))
      call random_number(rnd)
      ids(jj) = floor(interval(1)+(interval(2)-interval(1)+1)*rnd)
      ii=ii+1
    enddo
    if(ii.gt.max_it) ierr=1
  enddo
end subroutine gnu_rng_int_array_1d_norep

end module mod_gnu_rng

!> This module contains testcases for the event system
module event_timestep_spec
use mod_event_timestep
use mod_event, only: TICK
use fruit
implicit none

real*8, parameter :: division_tol = 1d-6
contains

subroutine test_two_events_weighted_average
  real*8, dimension(1) :: pusher_timestep = [1.d0]
  real*8, dimension(2) :: event_start = [0.d0, 0.d0]
  real*8, dimension(2) :: event_step  = [0.8d0, 1.2d0]
  logical, dimension(2,1) :: constraints = .true.
  integer :: ierr
  write(*,*) ' starting  test_two_events_weighted_average'
  call fix_event_timestep(pusher_timestep, event_start, event_step, constraints, ierr)
  call assert_equals(0, ierr, "must run without error")
  call assert_equals(0.94669509594882717d0, pusher_timestep(1), TICK, "timestep should be set to weighted average")
  call assert_equals(0.94669509594882695d0, event_step(1), TICK, "timestep 1 should be set to weighted average")
  call assert_equals(0.94669509594882695d0, event_step(2), TICK, "timestep 2 should be set to weighted average")

  write(*,*) 'completed test_two_events_weighted_average'

end subroutine test_two_events_weighted_average

subroutine test_two_events_start_only
  real*8, dimension(1) :: pusher_timestep = [1.d0]
  real*8, dimension(2) :: event_start = [0.d0, 0.d0]
  real*8, dimension(2) :: event_step  = [huge(0.d0), huge(0.d0)]
  logical, dimension(2,1) :: constraints = .true.
  integer :: ierr
  call fix_event_timestep(pusher_timestep, event_start, event_step, constraints, ierr)
  call assert_equals(0, ierr, "must run without error")
  call assert_equals(1.d0, pusher_timestep(1), TICK, "timestep should be set to weighted average")
  call assert_equals(huge(0.d0), event_step(1), TICK, "timestep 1 should be set to weighted average")
  call assert_equals(huge(0.d0), event_step(2), TICK, "timestep 2 should be set to weighted average")

  write(*,*) 'completed test_two_events_start_only'

end subroutine test_two_events_start_only

subroutine test_two_events_start_only_two
  real*8, dimension(1) :: pusher_timestep = [1.d0]
  real*8, dimension(2) :: event_start = [0.d0, 1d0]
  real*8, dimension(2) :: event_step  = [huge(0.d0), huge(0.d0)]
  logical, dimension(2,1) :: constraints = .true.
  integer :: ierr
  call fix_event_timestep(pusher_timestep, event_start, event_step, constraints, ierr)
  call assert_equals(0, ierr, "must run without error")
  call assert_equals(1.d0, pusher_timestep(1), TICK, "timestep should be set to weighted average")
  call assert_equals(huge(0.d0), event_step(1), TICK, "timestep 1 should stay huge")
  call assert_equals(huge(0.d0), event_step(2), TICK, "timestep 2 should stay huge")
  call assert_equals(1d0, event_start(2), TICK, "start 2 should stay same")

  write(*,*) ' completed test_two_events_start_only_two'

end subroutine test_two_events_start_only_two

subroutine test_three_events_one_start
  real*8, dimension(1) :: pusher_timestep = [5d-9]
  real*8, dimension(3) :: event_start = [0.d0, 1.8082299789290283d-003, 0.d0]
  real*8, dimension(3) :: event_step  = [huge(0.d0), huge(0.d0), huge(0.d0)]
  logical, dimension(3,1) :: constraints = .true.
  integer :: ierr
  call fix_event_timestep(pusher_timestep, event_start, event_step, constraints, ierr)
  call assert_equals(0, ierr, "must run without error")
  call assert_equals(5d-9, pusher_timestep(1), TICK, "timestep should be set to weighted average")
  call assert_equals(0.d0, event_start(1), TICK, "don't alter zero start time")
  call assert_equals(1.8082299789290283d-003, event_start(2), TICK, "don't alter start time?")
  call assert_equals(0.d0, event_start(3), TICK, "don't alter zero start time")
  call assert_equals(huge(0.d0), event_step(1), TICK, "timestep 1 should be set to weighted average")
  call assert_equals(huge(0.d0), event_step(2), TICK, "timestep 2 should be set to weighted average")
  call assert_equals(huge(0.d0), event_step(3), TICK, "timestep 3 should be set to weighted average")

  write(*,*) ' completed  test_three_events_one_start'

end subroutine test_three_events_one_start


subroutine test_two_events_fit_exactly
  real*8, dimension(1) :: pusher_timestep = [1.d-3]
  real*8, dimension(2) :: event_start = [0.d0, 0.d0]
  real*8, dimension(2) :: event_step  = [0.8d0, 1.2d0]
  logical, dimension(2,1) :: constraints = .true.
  integer :: ierr
  call fix_event_timestep(pusher_timestep, event_start, event_step, constraints, ierr)
  call assert_equals(0, ierr, "must run without error")
  call assert_equals(1d-3, pusher_timestep(1), TICK, "timestep should be set to weighted average")
  call assert_equals(0.8d0, event_step(1), TICK, "timestep 1 should remain same")
  call assert_equals(1.2d0, event_step(2), TICK, "timestep 2 should remain same")

  write(*,*) ' completed test_two_events_fit_exactly'

end subroutine test_two_events_fit_exactly


subroutine test_two_events_fit_almost
  real*8, dimension(1) :: pusher_timestep = [1.1d-3]
  real*8, dimension(2) :: event_start = [0.d0, 0.d0]
  real*8, dimension(2) :: event_step  = [0.8d0, 1.2d0]
  logical, dimension(2,1) :: constraints = .true.
  integer :: ierr
  call fix_event_timestep(pusher_timestep, event_start, event_step, constraints, ierr)
  call assert_equals(0, ierr, "must run without error")
  call assert_equals(1.1001069111187256d-3, pusher_timestep(1), TICK, "timestep should be set to weighted average")
  call assert_equals(727d0, event_step(1)/pusher_timestep(1), division_tol, "timestep ratio should be integer")
  call assert_equals(1091d0, event_step(2)/pusher_timestep(1), division_tol, "timestep ratio should be integer")

  write(*,*) ' completed test_two_events_fit_almost'

end subroutine test_two_events_fit_almost


subroutine test_two_pushers_one_event
  real*8, dimension(2) :: pusher_timestep = [1.1d-3, 0.9d-3]
  real*8, dimension(1) :: event_start = [0.d0]
  real*8, dimension(1) :: event_step  = [1.0d0]
  logical, dimension(1,2) :: constraints = .true.
  integer :: ierr
  call fix_event_timestep(pusher_timestep, event_start, event_step, constraints, ierr)
  call assert_equals(0, ierr, "must run without error")
  call assert_equals(909d0, event_step(1)/pusher_timestep(1), division_tol, "timestep 1 should divide event")
  call assert_equals(1111d0, event_step(1)/pusher_timestep(2), division_tol, "timestep 2 should divide event")

  write(*,*) ' completed test_two_pushers_one_event'

end subroutine test_two_pushers_one_event


subroutine test_two_pushers_two_events
  real*8, dimension(2) :: pusher_timestep = [1.1d-3, 0.9d-3]
  real*8, dimension(2) :: event_start = [0.d0, 0.d0]
  real*8, dimension(2) :: event_step  = [1.0d0, 2.0d0]
  logical, dimension(2,2) :: constraints = .true.
  integer :: ierr
  call fix_event_timestep(pusher_timestep, event_start, event_step, constraints, ierr)
  call assert_equals(0, ierr, "must run without error")
  call assert_equals(909d0, event_step(1)/pusher_timestep(1), division_tol, "timestep 1 should divide event 1")
  call assert_equals(1111d0, event_step(1)/pusher_timestep(2), division_tol, "timestep 2 should divide event 1")
  call assert_equals(2*909d0, event_step(2)/pusher_timestep(1), division_tol, "timestep 1 should divide event 2")
  call assert_equals(2*1111d0, event_step(2)/pusher_timestep(2), division_tol, "timestep 2 should divide event 2")

  write(*,*) ' completed test_two_pushers_two_events'

end subroutine test_two_pushers_two_events

subroutine test_two_pushers_two_events_some_constraints
  real*8, dimension(2) :: pusher_timestep = [1.1d-3, 0.9d-3]
  real*8, dimension(2) :: event_start = [0.d0, 0.d0]
  real*8, dimension(2) :: event_step  = [1.0d0, 2.0d0]
  logical, dimension(2,2) :: constraints = .true.
  integer :: ierr
  constraints(1,2) = .false.
  call fix_event_timestep(pusher_timestep, event_start, event_step, constraints, ierr)
  call assert_equals(0, ierr, "must run without error")
  call assert_equals(909d0, event_step(1)/pusher_timestep(1), division_tol, "timestep 1 should divide event 1")
  call assert_equals(2*909d0, event_step(2)/pusher_timestep(1), division_tol, "timestep 1 should divide event 2")
  call assert_equals(2*1111d0, event_step(2)/pusher_timestep(2), division_tol, "timestep 2 should divide event 2")

  write(*,*) ' completed test_two_pushers_two_events_some_constraints'

end subroutine test_two_pushers_two_events_some_constraints

subroutine test_two_pushers_two_events_of_which_one_constrained
  real*8, dimension(2) :: pusher_timestep = [1.1d-3, 0.9d-3]
  real*8, dimension(2) :: event_start = [0.d0, 0.d0]
  real*8, dimension(2) :: event_step  = [1.0d0, 2.0d0]
  logical, dimension(2,2) :: constraints = .true.
  integer :: ierr
  call fix_event_timestep(pusher_timestep, event_start, event_step, constraints, ierr)
  call assert_equals(0, ierr, "must run without error")
  call assert_equals(2*909d0, event_step(2)/pusher_timestep(1), division_tol, "timestep 1 should divide event 2")
  call assert_equals(2*1111d0, event_step(2)/pusher_timestep(2), division_tol, "timestep 2 should divide event 2")

  write(*,*) ' completed test_two_pushers_two_events_of_which_one_constrained'

end subroutine test_two_pushers_two_events_of_which_one_constrained

!> Test when there are many many events at the same time
subroutine test_one_pusher_six_events
  real*8, dimension(1) :: pusher_timestep = [1.0d0]
  real*8, dimension(6) :: event_start = 0.d0
  real*8, dimension(6) :: event_step  = [1d0, 2d0, 3d0, 4d0, 5d0, 6d0]
  logical, dimension(6,1) :: constraints = .true.
  integer :: ierr
  call fix_event_timestep(pusher_timestep, event_start, event_step, constraints, ierr)
  call assert_equals(0, ierr, "must run without error")
  call assert_equals(1d0, event_step(1)/pusher_timestep(1), division_tol, "timestep 1 should divide event 1")
  call assert_equals(2d0, event_step(2)/pusher_timestep(1), division_tol, "timestep 1 should divide event 2")
  call assert_equals(3d0, event_step(3)/pusher_timestep(1), division_tol, "timestep 1 should divide event 3")
  call assert_equals(4d0, event_step(4)/pusher_timestep(1), division_tol, "timestep 1 should divide event 4")
  call assert_equals(5d0, event_step(5)/pusher_timestep(1), division_tol, "timestep 1 should divide event 5")
  call assert_equals(6d0, event_step(6)/pusher_timestep(1), division_tol, "timestep 1 should divide event 6")

  write(*,*) ' completed  test_one_pusher_six_events'

end subroutine test_one_pusher_six_events

!> Test when there are many many events at the same time, with a start time
subroutine test_one_pusher_six_events_with_start
  real*8, dimension(1) :: pusher_timestep = [1.0d0]
  real*8, dimension(6) :: event_start = [1d0, 2d0, 3d0, 4d0, 5d0, 6d0]
  real*8, dimension(6) :: event_step  = [1d0, 2d0, 3d0, 4d0, 5d0, 6d0]
  logical, dimension(6,1) :: constraints = .true.
  integer :: ierr
  call fix_event_timestep(pusher_timestep, event_start, event_step, constraints, ierr)
  call assert_equals(0, ierr, "must run without error")
  call assert_equals(1d0, event_step(1)/pusher_timestep(1), division_tol, "timestep 1 should divide event 1")
  call assert_equals(2d0, event_step(2)/pusher_timestep(1), division_tol, "timestep 1 should divide event 2")
  call assert_equals(3d0, event_step(3)/pusher_timestep(1), division_tol, "timestep 1 should divide event 3")
  call assert_equals(4d0, event_step(4)/pusher_timestep(1), division_tol, "timestep 1 should divide event 4")
  call assert_equals(5d0, event_step(5)/pusher_timestep(1), division_tol, "timestep 1 should divide event 5")
  call assert_equals(6d0, event_step(6)/pusher_timestep(1), division_tol, "timestep 1 should divide event 6")

  write(*,*) ' completed test_one_pusher_six_events_with_start'

end subroutine test_one_pusher_six_events_with_start

subroutine test_three_pushers_six_events_with_start
  real*8, dimension(3) :: pusher_timestep = [1d0, 5d-1, 1d-1]
  real*8, dimension(6) :: event_start = 0.d0 ![1d0, 2d0, 3d0, 4d0, 5d0, 6d0]
  real*8, dimension(6) :: event_step  = [1d0, 2d0, 3d0, 4d0, 5d0, 6d0]
  logical, dimension(6,3) :: constraints = .true.
  integer :: ierr
  call fix_event_timestep(pusher_timestep, event_start, event_step, constraints, ierr)
  call assert_equals(0, ierr, "must run without error")
  call assert_equals(1d0, event_step(1)/pusher_timestep(1), division_tol, "timestep 1 should divide event 1")
  call assert_equals(2d0, event_step(2)/pusher_timestep(1), division_tol, "timestep 1 should divide event 2")
  call assert_equals(3d0, event_step(3)/pusher_timestep(1), division_tol, "timestep 1 should divide event 3")
  call assert_equals(4d0, event_step(4)/pusher_timestep(1), division_tol, "timestep 1 should divide event 4")
  call assert_equals(5d0, event_step(5)/pusher_timestep(1), division_tol, "timestep 1 should divide event 5")
  call assert_equals(6d0, event_step(6)/pusher_timestep(1), division_tol, "timestep 1 should divide event 6")
  call assert_equals(2d0, event_step(1)/pusher_timestep(2), division_tol, "timestep 2 should divide event 1")
  call assert_equals(4d0, event_step(2)/pusher_timestep(2), division_tol, "timestep 2 should divide event 2")
  call assert_equals(6d0, event_step(3)/pusher_timestep(2), division_tol, "timestep 2 should divide event 3")
  call assert_equals(8d0, event_step(4)/pusher_timestep(2), division_tol, "timestep 2 should divide event 4")
  call assert_equals(1d1, event_step(5)/pusher_timestep(2), division_tol, "timestep 2 should divide event 5")
  call assert_equals(1.2d1, event_step(6)/pusher_timestep(2), division_tol, "timestep 2 should divide event 6")

  write(*,*) ' completed test_one_pusher_six_events_with_start'

end subroutine test_three_pushers_six_events_with_start

subroutine test_single_event_start
  real*8, dimension(1) :: pusher_timestep = [0.1d0]
  real*8, dimension(1) :: event_start = [1.d0]
  real*8, dimension(1) :: event_step = [huge(1.d0)]
  logical, dimension(1,1) :: constraints = .true.
  integer :: ierr
  call fix_event_timestep(pusher_timestep, event_start, event_step, constraints, ierr)
  call assert_equals(0, ierr, "must run without error")
  call assert_equals(1d1, event_start(1)/pusher_timestep(1), division_tol, "timestep 1 should divide event start 1")
  call assert_equals(huge(1.d0), event_step(1), TICK, "event huge step should not change")

  write(*,*) ' completed test_one_pusher_six_events_with_start'

end subroutine test_single_event_start

end module event_timestep_spec

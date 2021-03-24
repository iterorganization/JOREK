!> Module containing examples from fruit unit testing
module test_fruit_example
use fruit

!> You can have module-global variables and helper subroutines
logical :: is_initialized = .false.
contains

!> This runs before any test subroutines
subroutine setup_test_fruit_example
is_initialized = .true.
end subroutine setup_test_fruit_example

!> A demonstration of the different options available here
subroutine test_asserts
  call assert_true(.true., "Simple assert true example")
  call assert_true(is_initialized, "The initializer ran")
  call assert_equals(1, 1, "Equality assertion")
  call assert_equals(1.0, 1.2, 0.4, "Equality assertion with tolerance")
end subroutine test_asserts

!> This runs after any test subroutines
subroutine teardown_test_fruit_example
is_initialized = .false.
end subroutine teardown_test_fruit_example
end module test_fruit_example

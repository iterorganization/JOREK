!> Program to calculate the frictional slowdown of particles in a background plasma
!> without a temperature gradient or magnetic field.
!> Requires gnuplot to generate figures.
!> Figures are placed in media/tests/collisions/ *.png
program test_plot_collisions
  use mod_test_collisions
  use mod_collision_cases_homma2012
  use constants
  implicit none

  integer, parameter :: n_particles = 50000
  real*8, parameter  :: timesteps(1) = [2.03d-10]
  real*8, parameter  :: time_end = 7.034d-7
  real*8, parameter  :: time_out = 1d-8
  type(collision_case) :: my_case

  real*8 :: v_out(3,n_particles)
  real*8, dimension(n_cases) :: slope_t
  integer :: u, i, i_case

  do i_case=1,n_cases
    my_case = homma_2012_collision_case(i_case)
    call run_collision_test(my_case, n_particles, timesteps(1), nint(time_end/timesteps(1)), &
        nint(time_end/time_out), v_out)
    call rm('out'//case_numbers(i_case)//'.txt')
    open(newunit=u, file='out'//case_numbers(i_case)//'.txt', status='new')
    do i=1,nint(time_end/time_out)
      write(u,*) time_out*real(i)/my_case%slowdown_time, v_out(1,i), v_out(3,i)
    end do
    close(u)
    ! Normalize slope with slowdown time (x-axis), initial velocity (z-axis)
    slope_t(i_case) = f_theo(i_case)*1d-17 * my_case%slowdown_time /(183.84d0*ATOMIC_MASS_UNIT * my_case%v_0(3))
  end do

  ! Create a temporary file for gnuplot
  call rm('plot_collisions.gp')
  open(newunit=u, file='plot_collisions.gp', status='new')
  write(u,"(A)") '&
      set terminal png; &
      set key top right; &
      set xlabel "t/tau_s"'

  i_case = 1
  write(u,*) 'set ylabel "v_x/v_0"'
  write(u,*) 'set output "media/tests/collisions/'//case_numbers(i_case)//'_x.png"'
  write(u,*) 'plot "out'//case_numbers(i_case)//'.txt" u 1:2 w p t "Vx: Ref Cal."'

  write(u,*) 'set ylabel "v_z/v_0"'
  write(u,*) 'set output "media/tests/collisions/'//case_numbers(i_case)//'_z.png"'
  write(u,*) 'plot "out'//case_numbers(i_case)//'.txt" u 1:3 w p t "Vz: Ref Cal.", \'
  write(u,'(A,g11.4,A)') '     1+', slope_t(1), '*x w l t "Ref Theo."'

  i_case = 2!and 3
  write(u,*) 'set ylabel "v_z/v_0"'
  write(u,*) 'set output "media/tests/collisions/1_z.png"'
  write(u,*) 'plot "out1-1.txt" u 1:3 w p t "Vz: 1-1 Cal.", \'
  write(u,'(A,g11.4,A)') '     1+', slope_t(2), '*x w l t "Vz: 1-1 Theo.", \'
  write(u,*) '     "out1-2.txt" u 1:3 w p t "Vz: 1-2 Cal.", \'
  write(u,'(A,g11.4,A)') '     1+', slope_t(3), '*x w l t "Vz: 1-2 Theo.", \'
  write(u,*) '     "out0-0.txt" u 1:3 w p t "Vz: Ref Cal.", \'
  write(u,*) '     exp(-x) w l t "Vz: Ref Theo."'

  i_case = 4!and 5
  write(u,*) 'set ylabel "v_x/v_0"'
  write(u,*) 'set output "media/tests/collisions/2_x.png"'
  write(u,*) 'plot "out2-1.txt" u 1:2 w p t "Vx: 2-1 Cal.", \'
  write(u,'(A,g11.4,A)') '     ', slope_t(4), '*x w l t "Vx: 2-1 Theo.", \'
  write(u,*) '     "out2-2.txt" u 1:2 w p t "Vx: 2-2 Cal.", \'
  write(u,'(A,g11.4,A)') '     ', slope_t(5), '*x w l t "Vx: 2-2 Theo.", \'
  write(u,*) '     "out0-0.txt" u 1:2 w p t "Vx: Ref Cal."'

  i_case = 6!and 7
  write(u,*) 'set ylabel "v_z/v_0"'
  write(u,*) 'set output "media/tests/collisions/3_z.png"'
  write(u,*) 'plot "out3-1.txt" u 1:3 w p t "Vz: 3-1 Cal.", \'
  write(u,'(A,g11.4,A)') '     1+', slope_t(6), '*x w l t "Vz: 3-1 Theo.", \'
  write(u,*) '     "out3-2.txt" u 1:3 w p t "Vz: 3-2 Cal.", \'
  write(u,'(A,g11.4,A)') '     1+', slope_t(7), '*x w l t "Vz: 3-2 Theo."'

  i_case = 8!and 9
  write(u,*) 'set ylabel "v_z/v_0"'
  write(u,*) 'set output "media/tests/collisions/4_z.png"'
  write(u,*) 'plot "out4-1.txt" u 1:3 w p t "Vz: 4-1 Cal.", \'
  write(u,'(A,g11.4,A)') '     1+', slope_t(8), '*x w l t "Vz: 4-1 Theo.", \'
  write(u,*) '     "out4-2.txt" u 1:3 w p t "Vz: 4-2 Cal.", \'
  write(u,'(A,g11.4,A)') '     1+', slope_t(9), '*x w l t "Vz: 4-2 Theo.", \'
  write(u,*) '     "out1-2.txt" u 1:3 w p t "Vz: 1-2 Cal.", \'
  write(u,'(A,g11.4,A)') '     1+', slope_t(3), '*x w l notitle'

  close(u)
  call system('sleep 2s && gnuplot plot_collisions.gp')

  ! delete the created files again
  !call rm('plot_collisions.gp')
  do i_case=1,n_cases
    !call rm('out'//case_numbers(i_case)//'.txt')
  end do

contains

subroutine rm(file)
  character(len=*), intent(in) :: file
  integer :: u, stat
  open(newunit=u, iostat=stat, file=file, status='old')
  if (stat .eq. 0) close(u, status='delete')
end subroutine rm
end program test_plot_collisions

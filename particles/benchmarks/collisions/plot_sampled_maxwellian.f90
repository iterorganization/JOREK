!> Program to sample angle distribution for testing binary collision model
program plot_sampled_maxwellian
  use mod_collision_cases_homma2012
  use mod_collisions
  use mod_pcg32_rng
  use mod_random_seed
  use constants
  implicit none

  integer, parameter :: n_particles_batch = 500
  integer, parameter :: n_batches = 70000/n_particles_batch
  integer, parameter :: n_bins = 180
  real*8, parameter :: gradkT_z(4) = [0.d0, 1d1, 2d1, 3d1]*EL_CHG ! [J/m]
  type(collision_case) :: my_case
  integer :: i, j, u, ifail, i_batch
  type(pcg32_rng) :: rng
  real*8, dimension(:,:), allocatable :: rans
  real*8, dimension(:,:), allocatable :: v_p
  real*8, dimension(n_particles_batch) :: theta
  integer, dimension(n_bins) :: counts
  real :: binstart, binend, binmid, A, v_thermal
  character(len=2) :: i_s

  call rng%initialize(6, random_seed(), 0, 0, ifail)
  allocate(rans(6,n_particles_batch))

  !$omp parallel do default(private) num_threads(4) firstprivate(rng)
  do i=1,4
    my_case = homma_2012_collision_case(1)
    my_case%grad_kT = [0.d0, 0.d0, gradkT_z(i)]
    write(*,'(A,i3,3e16.8)') " grad_kT : ",i,my_case%grad_kT
    allocate(v_p(3,n_particles_batch))

    ! delete the created files so we won't get a runtime error if they exist
    call rm('plot_sampled_maxwellian.gp'//i_s)
    call rm('sampled_theta_histogram.txt'//i_s)

    counts = 0

    A = 75.d0/2.d0 * PI*sqrt(PI) * EPS_ZERO**2/my_case%coulomb_log * (my_case%kT)/(my_case%n_b*(real(my_case%q_b)*EL_CHG)**4)
    v_thermal = sqrt(my_case%kT/(my_case%m_b*ATOMIC_MASS_UNIT))

    do i_batch = 1,n_batches
      v_p = 0.d0
      do j=1,n_particles_batch
        call rng%next(rans(:,j))
      end do
      call sample_velocity_dist_unmagnetized(n_particles_batch, rans, my_case%coulomb_log, &
          my_case%kT, my_case%grad_kT, my_case%n_b, my_case%m_b, my_case%q_b, v_b=[0.d0,0.d0,0.d0], w_1=v_p)
      ! Calculate histogram
      theta = acos(v_p(3,:)/norm2(v_p,dim=1))
      do j=0,n_bins-1
        binstart = real(j)/real(n_bins)
        binend   = (real(j)+1.d0)/real(n_bins)
        counts(j+1) = counts(j+1) + count((theta/PI .ge. binstart) .and. (theta/PI .lt. binend) .and. &
            (abs(norm2(v_p,dim=1)/v_thermal-1.d0) .lt. 1d-2))
      end do
    end do

    open(newunit=u, file='sampled_theta_histogram.txt'//i_s, status='new')
    do j=0,n_bins-1
      binstart = real(j)/real(n_bins)
      binend   = (real(j)+1.d0)/real(n_bins)
      binmid   = (binstart+binend)*0.5d0
      write(u,*) binmid, real(counts(j+1))/real(n_particles_batch*n_batches), &
      0.02 * 1/(2*pi*sqrt(2*pi))* exp(-0.5d0) * (&
      !0.02 * 0.0385108d0 * ( &
        1.d0 + A*0.8d0*cos(binmid*PI)*norm2(my_case%grad_kT)) * sin(binmid*PI) * 2*PI * PI/real(n_bins)
      ! In this testcase coordinate systems 1 (lab frame) and 2 (grad_kT in z direction) are equal.
      ! Since the range in w is quite small we estimate the integral as (w2-w1)*f(w=1)*2pi
    end do
    close(u)

    ! Create a temporary file for gnuplot
    ! set terminal cairolatex, set output ... .tex
    open(newunit=u, file='plot_sampled_maxwellian.gp'//i_s, status='new')
    write(u,"(A,i0.2,A,i0.2,A)") '&
      set terminal png; &
      set key top right; &
      set output "sampled_theta_', nint(gradkT_z(i)/EL_CHG), '.png"; &
      set xlabel "$\\theta/\\pi$"; &
      set boxwidth; &
      plot "sampled_theta_histogram.txt'//i_s//'" u 1:2 w boxes t "$\\Nabla T = ', nint(gradkT_z(i)/EL_CHG), '$ eV/m $\\mathbf{e}_Z$ &
      &", "sampled_theta_histogram.txt'//i_s//'" u 1:3 w l t "Theo"'
    close(u)

    !$omp critical
    call system('gnuplot plot_sampled_maxwellian.gp'//i_s)
    !$omp end critical

    ! delete the created files again
    !call rm('plot_sampled_maxwellian.gp'//i_s)
    !call rm('sampled_theta_histogram.txt'//i_s)

    write(*,*) "Finished step ", i, "gradkT=", nint(gradKt_z(i)/EL_CHG)
    deallocate(v_p)
  end do
  !$omp end parallel do
contains

subroutine rm(file)
  character(len=*), intent(in) :: file
  integer :: u, stat
  open(newunit=u, iostat=stat, file=file, status='old')
  if (stat .eq. 0) close(u, status='delete')
end subroutine rm
end program plot_sampled_maxwellian

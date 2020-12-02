!> Module testing Coronal equilibrium for selected elements
!> W, Argon, Helium
module test_coronal_eq
use fruit
use mod_openadas
use mod_coronal
implicit none

character(len=5) :: sets(8) = ['89_b ', '96_c ', '96_li', '96_n ', '96_ne', '50_w ', '89_ar', '96_he']
type(ADF11_all)  :: adas(8)

contains

!> Download necessary files from the open-adas website
subroutine setup_test_coronal_eq
  integer :: i
  character(len=200) :: set
  set = ''
  do i=1,size(sets,1)
    set = trim(set) // ' ' // trim(sets(i))
  end do
  ! read all in one go
  call system('util/fetch_openadas.sh '//trim(set))
  ! Now we need to wait a bit for filesystem update
  call system('sleep 0.5')
  do i=1,size(sets,1)
    adas(i) = read_adf11(0,trim(sets(i)))
  end do
end subroutine setup_test_coronal_eq

!> For a single density test that <Z> is increasing with T
subroutine test_coronal_Z_increasing
  type(coronal) :: cor
  integer :: i, j, i_set
  real*8 :: Z_eff, Z_eff_old
  character(len=10) :: T_s

  do i_set = 1,size(sets,1)
    cor = coronal(adas(i_set))
    if (i_set .eq. 7) call output_coronal(cor)
    Z_eff_old = 0.d0
    do i=1,size(cor%temperature,1)
      Z_eff = sum(cor%Z(1, i, :) * [(j,j=0,cor%n_Z)])
      write(T_s,'(g10.3)') 10.d0**(cor%temperature(i))
      call assert_true(Z_eff .ge. Z_eff_old, '<Z> must be increasing with T ('//sets(i_set)//') T_e='//trim(T_s)//' K')
      Z_eff_old = Z_eff
    end do
  end do
end subroutine test_coronal_Z_increasing
end module test_coronal_eq

module mod_equations
  use mod_semianalytical
  use mod_parameters
  use data_structure, only: nbthreads
  implicit none

  type type_thread_eq
#ifdef DEBUG
    type(action), dimension(:), allocatable :: rhs1seq, rhs3seq
    type(action), dimension(:), allocatable :: amat11seq, amat13seq
    type(action), dimension(:), allocatable :: amat22seq
    type(action), dimension(:), allocatable :: amat33seq
    type(action), dimension(:), allocatable :: amat44seq
    type(action), dimension(:), allocatable :: amat55seq
    type(action), dimension(:), allocatable :: amat66seq
    type(action), dimension(:), allocatable :: aBv2seq, aBv2xseq, aBv2yseq, aBv2pseq
#endif

    ! Indices in eq array: variable index (see algexpr's below), R derivative order, z derivative order, phi derivative order,
    !   separation of terms with covariant phi derivatives in test function and unknown (FFT, not used in this model)
    real*8, dimension(:,:,:,:,:), allocatable :: eq
  end type type_thread_eq

  ! Values
  type(algexpr), parameter, private :: Psi0       = algexpr(basic=.true.,var=1)
  type(algexpr), parameter, private :: Phi0       = algexpr(basic=.true.,var=2)
  type(algexpr), parameter, private :: zj0        = algexpr(basic=.true.,var=3)
  type(algexpr), parameter, private :: w0         = algexpr(basic=.true.,var=4)
  type(algexpr), parameter, private :: rho0       = algexpr(basic=.true.,var=5)
  type(algexpr), parameter, private :: T0         = algexpr(basic=.true.,var=6)
  ! Test function
  type(algexpr), parameter, private :: v          = algexpr(basic=.true.,var=7)
  ! Unknowns
  type(algexpr), parameter, private :: Psi        = algexpr(basic=.true.,var=8)
  type(algexpr), parameter, private :: Phi        = algexpr(basic=.true.,var=8)
  type(algexpr), parameter, private :: zj         = algexpr(basic=.true.,var=8)
  type(algexpr), parameter, private :: w          = algexpr(basic=.true.,var=8)
  type(algexpr), parameter, private :: rho        = algexpr(basic=.true.,var=8)
  type(algexpr), parameter, private :: T          = algexpr(basic=.true.,var=8)
  ! Other quantities
  type(algexpr), parameter, private :: chi        = algexpr(basic=.true.,var=9)
  type(algexpr), parameter, private :: R          = algexpr(basic=.true.,var=10)
  ! Quantities imported from GVEC
  type(algexpr), parameter, private :: p0_gvec    = algexpr(basic=.true.,var=11)
  type(algexpr), parameter, private :: B0x_gvec   = algexpr(basic=.true.,var=12)
  type(algexpr), parameter, private :: B0y_gvec   = algexpr(basic=.true.,var=13)
  type(algexpr), parameter, private :: B0p_gvec   = algexpr(basic=.true.,var=14)
  ! Auxiliary variables (aux)
  type(algexpr), parameter, private :: Bv2        = algexpr(basic=.true.,var=15)

  type(algexpr), private :: rhs1, rhs2, rhs3, rhs4, rhs5, rhs6
  type(algexpr), private :: amat11, amat12, amat13, amat16
  type(algexpr), private :: amat21, amat22, amat23, amat24, amat25, amat26
  type(algexpr), private :: amat31, amat33
  type(algexpr), private :: amat42, amat44
  type(algexpr), private :: amat51, amat52, amat55
  type(algexpr), private :: amat61, amat62, amat63, amat65, amat66
  type(algexpr), private :: a_Bv2

  integer, parameter :: n_rhs = 2, n_amat = 7, n_aux = 4

  type(algexpr), private :: rhs2e, rhs3e, rhs4e, rhs5e, rhs6e
  type(algexpr), private :: amat22e, amat25e, amat26e
  type(algexpr), private :: amat31e
  type(algexpr), private :: amat52e, amat55e
  type(algexpr), private :: amat62e, amat65e, amat66e
  type(algexpr), private :: ea_Bv2x, ea_Bv2y, ea_Bv2p

  type(type_thread_eq), dimension(:), allocatable, target :: thread_eq

  contains

  subroutine init_equations()
    implicit none


    a_Bv2 = dx(chi)*dx(chi) + dy(chi)*dy(chi) + dp(chi)*dp(chi)/(R*R)

    rhs1 = (-Bv2)*inprod(v,Psi0)

    rhs3 = -dx(v)*(dy(chi)*B0p_gvec - dp(chi)*B0y_gvec/R) + dy(v)*(dx(chi)*B0p_gvec - dp(chi)*B0x_gvec/R) &
         - dp(v)*(dx(chi)*B0y_gvec - dy(chi)*B0x_gvec)/R

    amat11 = Bv2*inprod(v,Psi)
    amat13 = v*Bv2*zj

    amat22 = v*Phi

    amat33 = v*Bv2*zj

    amat44 = v*w

    amat55 = v*rho

    amat66 = v*T

    ea_Bv2x = Dexpand(deepcopy(dx(a_Bv2))); ea_Bv2y = Dexpand(deepcopy(dy(a_Bv2))); ea_Bv2p = Dexpand(deepcopy(dp(a_Bv2)))
  end subroutine init_equations

  subroutine init_eq_struct()
    use data_structure, only: nbthreads
    implicit none
    integer :: i

    if (.not. allocated(thread_eq)) then
      allocate(thread_eq(nbthreads))
      do i=1,nbthreads
        allocate(thread_eq(i)%eq(n_var+9,0:n_order-1,0:n_order-1,0:n_order-1,4))
#ifdef DEBUG
        allocate(thread_eq(i)%rhs1seq(countsubexprs(rhs1)))
        allocate(thread_eq(i)%rhs3seq(countsubexprs(rhs3)))
        allocate(thread_eq(i)%amat11seq(countsubexprs(amat11)))
        allocate(thread_eq(i)%amat13seq(countsubexprs(amat13)))
        allocate(thread_eq(i)%amat22seq(countsubexprs(amat22)))
        allocate(thread_eq(i)%amat33seq(countsubexprs(amat33)))
        allocate(thread_eq(i)%amat44seq(countsubexprs(amat44)))
        allocate(thread_eq(i)%amat55seq(countsubexprs(amat55)))
        allocate(thread_eq(i)%amat66seq(countsubexprs(amat66)))
        allocate(thread_eq(i)%aBv2seq(countsubexprs(a_Bv2)))
        allocate(thread_eq(i)%aBv2xseq(countsubexprs(ea_Bv2x)))
        allocate(thread_eq(i)%aBv2yseq(countsubexprs(ea_Bv2y)))
        allocate(thread_eq(i)%aBv2pseq(countsubexprs(ea_Bv2p)))
#endif
      end do
    end if
  end subroutine init_eq_struct

#ifdef DEBUG
  subroutine build_all_seq()
    use data_structure, only: nbthreads
    implicit none
    integer :: i

    do i=1,nbthreads
      call buildsequence(rhs1, thread_eq(i)%rhs1seq, thread_eq(i)%eq)
      call buildsequence(rhs3, thread_eq(i)%rhs3seq, thread_eq(i)%eq)

      call buildsequence(amat11, thread_eq(i)%amat11seq, thread_eq(i)%eq)
      call buildsequence(amat13, thread_eq(i)%amat13seq, thread_eq(i)%eq)

      call buildsequence(amat22, thread_eq(i)%amat22seq, thread_eq(i)%eq)

      call buildsequence(amat33, thread_eq(i)%amat33seq, thread_eq(i)%eq)

      call buildsequence(amat44, thread_eq(i)%amat44seq, thread_eq(i)%eq)

      call buildsequence(amat55, thread_eq(i)%amat55seq, thread_eq(i)%eq)

      call buildsequence(amat66, thread_eq(i)%amat66seq, thread_eq(i)%eq)

      call buildsequence(a_Bv2, thread_eq(i)%aBv2seq, thread_eq(i)%eq)
      call buildsequence(ea_Bv2x, thread_eq(i)%aBv2xseq, thread_eq(i)%eq)
      call buildsequence(ea_Bv2y, thread_eq(i)%aBv2yseq, thread_eq(i)%eq)
      call buildsequence(ea_Bv2p, thread_eq(i)%aBv2pseq, thread_eq(i)%eq)
    end do
  end subroutine build_all_seq
#endif

  subroutine get_rhs(rhs,varnames)
    implicit none
    type(algexpr), dimension(n_rhs), intent(out) :: rhs
    character(8),  dimension(n_rhs), intent(out) :: varnames

    rhs = (/ rhs1, rhs3 /)
    varnames = (/ "rhs_ij_1", "rhs_ij_3" /)
  end subroutine get_rhs

  subroutine get_amat(amat,varnames)
    implicit none
    type(algexpr), dimension(n_amat), intent(out) :: amat
    character(7),  dimension(n_amat), intent(out) :: varnames

    amat = (/ amat11,         amat13, &
                      amat22, &
                              amat33, &
                                      amat44, &
                                              amat55, &
                                                      amat66 /)
    varnames = (/ "amat_11",            "amat_13", &
                             "amat_22", &
                                        "amat_33", &
                                                   "amat_44", &
                                                              "amat_55", &
                                                                         "amat_66" /)
  end subroutine get_amat
  
  subroutine get_aux(aux,varnames)
    implicit none
    type(algexpr), dimension(n_aux), intent(out) :: aux
    character(14), dimension(n_aux), intent(out) :: varnames
    integer      :: i
    character(2) :: num

    aux = (/ a_Bv2, ea_Bv2x, ea_Bv2y, ea_Bv2p /)
    varnames = (/ "eq(15,0,0,0,:)", "eq(15,1,0,0,:)", "eq(15,0,1,0,:)", "eq(15,0,0,1,:)" /)
  end subroutine get_aux

  type(algexpr) function Bv_pbrack(a,b)
    implicit none
    type(algexpr), intent(in) :: a, b

    Bv_pbrack = ((dy(a)*dp(b) - dp(a)*dy(b))*dx(chi) + (dp(a)*dx(b) - dx(a)*dp(b))*dy(chi) + (dx(a)*dy(b) - dy(a)*dx(b))*dp(chi))/R
  end function Bv_pbrack

  type(algexpr) function Bv_parderiv(a)
    implicit none
    type(algexpr), intent(in) :: a

    Bv_parderiv = dx(a)*dx(chi) + dy(a)*dy(chi) + dp(a)*dp(chi)/(R*R)
  end function Bv_parderiv

  type(algexpr) function gradprod(a,b)
    implicit none
    type(algexpr), intent(in) :: a, b

    gradprod = dx(a)*dx(b) + dy(a)*dy(b) + dp(a)*dp(b)/(R*R)
  end function gradprod

  type(algexpr) function inprod(a,b)
    implicit none
    type(algexpr), intent(in) :: a, b

    inprod = gradprod(a,b) - Bv_parderiv(a)*Bv_parderiv(b)/Bv2
  end function inprod

  type(algexpr) function Lap(a)
    implicit none
    type(algexpr), intent(in) :: a

    Lap = dx(R*dx(a))/R + dy(dy(a)) + dp(dp(a))/(R*R)
  end function Lap

  type(algexpr) function pLap(a)
    implicit none
    type(algexpr), intent(in) :: a

    pLap = Lap(a) - Bv_parderiv(Bv_parderiv(a)/Bv2)
  end function pLap
end module mod_equations

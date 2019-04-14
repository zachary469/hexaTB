module fun_mod
use, intrinsic :: iso_fortran_env
use :: dtype_mod, only: dp
implicit none
public
real (kind=dp), parameter :: AA = 2.46_dp ! graphene lattice parameter
integer, parameter :: PARAMSIZE = 7 ! number of TB parameters in bilayer graphene model
integer, parameter :: N_SIZE = 4 ! size of bilayer graphene Hamiltonian
contains

! The following subroutine sets the parameters for bilayer graphene found in
! E. McCann and M. Koshino, Rep. Prog. Phys. 76, 056503 (2013)
  subroutine tb_params(params)
    real (kind=dp),intent(out) :: params(7)
    params(1)=0.0_dp ! U parameter
    params(2)=0.0_dp ! \delta_{AB} parameter
    params(3)=0.022_dp ! \Delta^{\prime} parameter
    params(4)=3.16_dp ! \gamma_0 parameter
    params(5)=0.381_dp ! \gamma_1 parameter
    params(6)=0.38_dp ! \gamma_3 parameter
    params(7)=0.14_dp ! \gamma_4 parameter
  end subroutine

! The following function computes Eq. (12) from
! E. McCann and M. Koshino, Rep. Prog. Phys. 76, 056503 (2013)
  function ff(aa, kk) result (ffres)
    implicit none
    complex (kind=dp) :: ffres ! function value
    real (kind=dp) :: phase ! temporary complex phase for evaluating terms in ff
    real (kind=dp), intent(in) :: aa ! lattice parameter input value
    real (kind=dp), dimension(:), intent(in) :: kk ! quasimomentum vector input value
! Below we assume that kk has at least a dimension of 2:
    phase = aa * kk(2) / sqrt(3.0_dp)
!    ffres = cmplx(cos(phase),sin(phase))
!    phase= -aa * kk(2) / (2.0_dp * sqrt(3.0_dp))
    ffres = cmplx(cos(phase),sin(phase)) + &
&( cmplx(cos(-phase/2.0_dp),sin(-phase/2.0_dp)) * cmplx(2.0_dp * cos(aa*kk(1)/2.0_dp),0.0_dp) )
  end function

! The following function computes Eq. (16) from
! E. McCann and M. Koshino, Rep. Prog. Phys. 76, 056503 (2013)
  function Ham(aa, kk, params, n_size) result (Hamres)
    implicit none
    integer, intent(in) :: n_size ! size of Hamiltonian
    complex (kind=dp) :: ffint ! local variable
    real (kind=dp), intent(in) :: aa ! lattice parameter input value
    real (kind=dp), dimension(:), intent(in) :: kk ! quasimomentum vector input value
    real (kind=dp), dimension(:), intent(in) :: params ! tight-binding parameters
    complex (kind=dp),allocatable :: Hamres(:,:) ! function value
    real (kind=dp) :: UU,delta_AB,delta_prime,gamma_0,gamma_1,gamma_3,gamma_4 ! TB model parameters
! Here we define the physically meaningful parameters of the TB model to make definition of
! the Hamiltonian clearer:
    UU=params(1)
    delta_AB=params(2)
    delta_prime=params(3)
    gamma_0=params(4)
    gamma_1=params(5)
    gamma_3=params(6)
    gamma_4=params(7)
    allocate(Hamres(n_size,n_size))
    Hamres=0.0_dp
    ffint=ff(aa,kk)
    Hamres(1,1)=(-UU+delta_AB)/2.0_dp
    Hamres(4,4)=-Hamres(1,1)
    Hamres(2,2)=delta_prime-(UU+delta_AB)/2.0_dp
    Hamres(3,3)=delta_prime+(UU+delta_AB)/2.0_dp
    Hamres(2,3)=gamma_1
    Hamres(3,2)=gamma_1
    Hamres(2,1)=-gamma_0*ffint
    Hamres(4,3)=-gamma_0*ffint
    Hamres(3,1)=gamma_4*ffint
    Hamres(4,2)=gamma_4*ffint
    Hamres(1,4)=-gamma_3*ffint
    Hamres(1,2)=conjg(Hamres(2,1))
    Hamres(3,4)=conjg(Hamres(4,3))
    Hamres(1,3)=conjg(Hamres(3,1))
    Hamres(2,4)=conjg(Hamres(4,2))
    Hamres(4,1)=conjg(Hamres(1,4))
    end function

end module fun_mod


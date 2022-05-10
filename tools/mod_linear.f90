!> This module contains procedures for computing the linear interpolation of
!> vectors and matrices
module mod_linear
  implicit none
  private
  public linear_interp_differentials, linear_interp_differentials_dt

contains

  !> This function performs a linear interpolation given values and
  !> their differentials
  !> inputs:
  !>   n:      (integer) number of values
  !>   y_new:  (real8)(n) values for interpolation at the end time
  !>   dy_new: (real8)(n) differentials at the end time
  !>   df:     (real8) normalised time step (t_new-t)/(t_new-t_old)
  !> outputs:
  !>   y: (real8)(n) interpolated values
#ifdef CUDA_KERNELS
  attributes(host,device) &
#endif
  pure function linear_interp_differentials(n,y_new,dy_new,df) result(y)
    !> declare input variables
    integer, intent(in)                    :: n
    real(kind=8), intent(in)               :: df
    real(kind=8), dimension(n), intent(in) :: y_new, dy_new
    !> declare output variables
    real(kind=8), dimension(n) :: y

    !> compute linear interpolation
    y = y_new - df*dy_new
    
  end function linear_interp_differentials
  
  !> This procedure computes the linear interpolation of first order derivatives
  !> inputs:
  !>   n:          (integer) number of derivatives
  !>   dy_new:     (real8)(n) differentials at the end time
  !>   inverse_dt: (real8)(n) inverse of the interval duration
  !> outputs:
  !>   dydt:       (real8)(n) linear interpolation of derivatives
#ifdef CUDA_KERNELS
  attributes(host,device) &
#endif
  pure function linear_interp_differentials_dt(n,dy_new,inverse_dt) &
    result(dydt)
    !> declare input variables
    integer, intent(in)                    :: n
    real(kind=8), dimension(n), intent(in) :: dy_new
    real(kind=8), intent(in)               :: inverse_dt
    !> declare output variables
    real(kind=8), dimension(n) :: dydt

    !> compute derivatives
    dydt = dy_new*inverse_dt
    
  end function linear_interp_differentials_dt
  
end module mod_linear

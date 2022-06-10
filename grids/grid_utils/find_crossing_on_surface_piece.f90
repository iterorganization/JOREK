subroutine find_crossing_on_surface_piece(node_list,element_list,surface,piece, R_c,Z_c, &
                                          R_out,Z_out, r_flux,s_flux, t_tht,ifail, gofast)
  !-------------------------------------------------------------------------
  ! solves two non-linear equations using Newtons method (from numerical recipes)
  ! LU decomposition replaced by explicit solution of 2x2 matrix.
  !
  ! finds the crossing of two coordinate lines given as a series of cubics
  !-------------------------------------------------------------------------
  use data_structure
  use mod_interp, only: interp_RZ
  use phys_module, only: R_geo, surface_cross_tol
  
  implicit none
  
  ! --- Routine parameters
  type (type_node_list),    intent(in)          :: node_list
  type (type_element_list), intent(in)          :: element_list
  type (type_surface),      intent(in)          :: surface
  integer,                  intent(in)          :: piece
  real*8,                   intent(inout)       :: R_c(4),   Z_c(4)
  real*8,                   intent(inout)       :: R_out,    Z_out
  real*8,                   intent(inout)       :: r_flux,s_flux,t_tht
  integer,                  intent(inout)       :: ifail
  logical,                  intent(in)          :: gofast

  ! --- Local variables
  integer :: i, ntrial, istart
  integer :: i_elm
  real*8  :: t_flux,dr_flux, ds_flux
  real*8  :: rr1, drr1, rr2, drr2, ss1, dss1, ss2, dss2
  real*8  :: RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss
  real*8  :: ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss
  real*8  :: RRg2,dRRg2_dr,dRRg2_ds,dRRg2_drs,dRRg2_drr,dRRg2_dss
  real*8  :: ZZg2,dZZg2_dr,dZZg2_ds,dZZg2_drs,dZZg2_drr,dZZg2_dss
  real*8  :: dRRg1_dt, dZZg1_dt, dRRg2_dt, dZZg2_dt
  real*8  :: RR_flux, dRR_flux, RR_tht, dRR_tht,  ZZ_flux, dZZ_flux, ZZ_tht, dZZ_tht
  real*8  :: x(2), FVEC(2), FJAC(2,2), p(2), x_previous(2)
  real*8  :: tolx, tolf, errx, errf, temp, dis, max_step, tol_far
  real*8  :: RR_mid, ZZ_mid
  real*8  :: RR_mid_surf, ZZ_mid_surf

  tol_far = 0.2 * R_geo

  x_previous = 0.d0
  
  ifail = 99
  
  if ((R_c(1) .eq. R_c(3)) .and. (Z_c(1) .eq. Z_c(3))) then
    ifail = 9
    return
  endif
  
  rr1  = surface%s(1,piece);   ss1  = surface%t(1,piece)
  drr1 = surface%s(2,piece);   dss1 = surface%t(2,piece)
  rr2  = surface%s(3,piece);   ss2  = surface%t(3,piece)
  drr2 = surface%s(4,piece);   dss2 = surface%t(4,piece)
  
  i_elm = surface%elm(piece)

  call interp_RZ(node_list,element_list,i_elm,rr1,ss1,RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss, &
                                                      ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss)
  call interp_RZ(node_list,element_list,i_elm,rr2,ss2,RRg2,dRRg2_dr,dRRg2_ds,dRRg2_drs,dRRg2_drr,dRRg2_dss, &
                                                      ZZg2,dZZg2_dr,dZZg2_ds,dZZg2_drs,dZZg2_drr,dZZg2_dss)
  dRRg1_dt = dRRg1_dr * drr1 + dRRg1_ds * dss1
  dZZg1_dt = dZZg1_dr * drr1 + dZZg1_ds * dss1
  dRRg2_dt = dRRg2_dr * drr2 + dRRg2_ds * dss2
  dZZg2_dt = dZZg2_dr * drr2 + dZZg2_ds * dss2
  
  ntrial = 20
  tolx = 1.d-6
  tolf = 1.d-12

  ! --- First check end points
  if ( sqrt( (RRg1-R_c(1))**2 + (ZZg1-Z_c(1))**2) .lt. tolx )  then
    t_flux = -1.d0
    t_tht  = -1.d0
  
    call CUB1D(rr1, drr1, rr2, drr2, t_flux, r_flux, dr_flux)
    call CUB1D(ss1, dss1, ss2, dss2, t_flux, s_flux, ds_flux)
    call CUB1D(RRg1, dRRg1_dt, RRg2, dRRg2_dt, t_flux, RR_flux, dRR_flux)
    call CUB1D(ZZg1, dZZg1_dt, ZZg2, dZZg2_dt, t_flux, ZZ_flux, dZZ_flux)
    call CUB1D(R_c(1), R_c(2), R_c(3), R_c(4), t_tht, RR_tht, dRR_tht)
    call CUB1D(Z_c(1), Z_c(2), Z_c(3), Z_c(4), t_tht, ZZ_tht, dZZ_tht)
  
    R_out     = 0.5d0*(RR_tht + RR_flux)
    Z_out     = 0.5d0*(ZZ_tht + ZZ_flux)

    ifail = 0
    return
  endif
  ! --- First check end points
  if ( sqrt( (RRg1-R_c(3))**2 + (ZZg1-Z_c(3))**2) .lt. tolx )  then
    t_flux = -1.d0
    t_tht  = +1.d0
  
    call CUB1D(rr1, drr1, rr2, drr2, t_flux, r_flux, dr_flux)
    call CUB1D(ss1, dss1, ss2, dss2, t_flux, s_flux, ds_flux)
    call CUB1D(RRg1, dRRg1_dt, RRg2, dRRg2_dt, t_flux, RR_flux, dRR_flux)
    call CUB1D(ZZg1, dZZg1_dt, ZZg2, dZZg2_dt, t_flux, ZZ_flux, dZZ_flux)
    call CUB1D(R_c(1), R_c(2), R_c(3), R_c(4), t_tht, RR_tht, dRR_tht)
    call CUB1D(Z_c(1), Z_c(2), Z_c(3), Z_c(4), t_tht, ZZ_tht, dZZ_tht)
  
    R_out     = 0.5d0*(RR_tht + RR_flux)
    Z_out     = 0.5d0*(ZZ_tht + ZZ_flux)

    ifail = 0
    return
  endif
  ! --- First check end points
  if ( sqrt( (RRg2-R_c(1))**2 + (ZZg2-Z_c(1))**2) .lt. tolx )  then
    t_flux = +1.d0
    t_tht  = -1.d0
  
    call CUB1D(rr1, drr1, rr2, drr2, t_flux, r_flux, dr_flux)
    call CUB1D(ss1, dss1, ss2, dss2, t_flux, s_flux, ds_flux)
    call CUB1D(RRg1, dRRg1_dt, RRg2, dRRg2_dt, t_flux, RR_flux, dRR_flux)
    call CUB1D(ZZg1, dZZg1_dt, ZZg2, dZZg2_dt, t_flux, ZZ_flux, dZZ_flux)
    call CUB1D(R_c(1), R_c(2), R_c(3), R_c(4), t_tht, RR_tht, dRR_tht)
    call CUB1D(Z_c(1), Z_c(2), Z_c(3), Z_c(4), t_tht, ZZ_tht, dZZ_tht)
  
    R_out     = 0.5d0*(RR_tht + RR_flux)
    Z_out     = 0.5d0*(ZZ_tht + ZZ_flux)

    ifail = 0
    return
  endif
  ! --- First check end points
  if ( sqrt( (RRg2-R_c(3))**2 + (ZZg2-Z_c(3))**2) .lt. tolx )  then
    t_flux = +1.d0
    t_tht  = +1.d0
  
    call CUB1D(rr1, drr1, rr2, drr2, t_flux, r_flux, dr_flux)
    call CUB1D(ss1, dss1, ss2, dss2, t_flux, s_flux, ds_flux)
    call CUB1D(RRg1, dRRg1_dt, RRg2, dRRg2_dt, t_flux, RR_flux, dRR_flux)
    call CUB1D(ZZg1, dZZg1_dt, ZZg2, dZZg2_dt, t_flux, ZZ_flux, dZZ_flux)
    call CUB1D(R_c(1), R_c(2), R_c(3), R_c(4), t_tht, RR_tht, dRR_tht)
    call CUB1D(Z_c(1), Z_c(2), Z_c(3), Z_c(4), t_tht, ZZ_tht, dZZ_tht)
  
    R_out     = 0.5d0*(RR_tht + RR_flux)
    Z_out     = 0.5d0*(ZZ_tht + ZZ_flux)

    ifail = 0
    return
  endif

  ! --- Then try by steps
  do istart = 1,5

    if (istart .eq. 1) then
      x(1) = 0.d0
      x(2) = 0.d0
    elseif (istart .eq. 2) then
      x(1) = -0.71d0
      x(2) = -0.71d0
    elseif (istart .eq. 3) then
      x(1) =  0.71d0
      x(2) = -0.71d0
    elseif (istart .eq. 4) then
      x(1) =  0.71d0
      x(2) =  0.71d0
    elseif (istart .eq. 5) then
      x(1) = -0.71d0
      x(2) =  0.71d0
    endif
  
    ifail = 999
    max_step = surface_cross_tol
  
    do i=1,ntrial
  
      t_flux = x(1)
      t_tht  = x(2)
  
      call CUB1D(RRg1, dRRg1_dt, RRg2, dRRg2_dt, t_flux, RR_flux, dRR_flux)
      call CUB1D(ZZg1, dZZg1_dt, ZZg2, dZZg2_dt, t_flux, ZZ_flux, dZZ_flux)

      call CUB1D(R_c(1), R_c(2), R_c(3), R_c(4), t_tht, RR_tht, dRR_tht)
      call CUB1D(Z_c(1), Z_c(2), Z_c(3), Z_c(4), t_tht, ZZ_tht, dZZ_tht)
  
      FVEC(1)   = RR_tht - RR_flux
      FVEC(2)   = ZZ_tht - ZZ_flux
      FJAC(1,1) = - dRR_flux
      FJAC(1,2) =   dRR_tht
      FJAC(2,1) = - dZZ_flux
      FJAC(2,2) =   dZZ_tht
  
      errf=abs(fvec(1))+abs(fvec(2))
  
      !if (i .eq. ntrial) write(*,'(A,i3,4e16.8)') ' newton   : ',i,errf,errx,x

      if (errf .le. tolf) then
  
        t_flux = x(1)
        t_tht  = x(2)
  
        call CUB1D(rr1, drr1, rr2, drr2, t_flux, r_flux, dr_flux)
        call CUB1D(ss1, dss1, ss2, dss2, t_flux, s_flux, ds_flux)
  
        R_out     = 0.5d0*(RR_tht + RR_flux)
        Z_out     = 0.5d0*(ZZ_tht + ZZ_flux)
  
        !write(*,'(A,i3,4e16.8)') ' newton (1) : ',i,errf,errx,x

        ifail = 0
        return
      endif
  
      p = -fvec
  
      temp = p(1)
      dis  =  fjac(2,2)*fjac(1,1)-fjac(1,2)*fjac(2,1)
      if (dis .eq. 0.d0) dis = 1.d-10
      p(1) = (fjac(2,2)*p(1)     -fjac(1,2)*p(2)     )/dis
      p(2) = (fjac(1,1)*p(2)     -fjac(2,1)*temp     )/dis
  
      errx=abs(p(1)) + abs(p(2))
  
      p = min(p,+0.25d0)
      p = max(p,-0.25d0)
  
      x = x + p
  
      ! Sometimes you need to look outside the grid...
      x = max(x,-max_step)
      x = min(x,+max_step)
      if(abs(x_previous(1)) .eq. max_step) max_step = max_step + 0.002d0
      if(abs(x_previous(2)) .eq. max_step) max_step = max_step + 0.002d0
      x_previous(1) = x(1)
      x_previous(2) = x(2)
  
  
      if (errx .le. tolx) then

        t_flux = x(1)
        t_tht  = x(2)
  
        call CUB1D(rr1, drr1, rr2, drr2, t_flux, r_flux, dr_flux)
        call CUB1D(ss1, dss1, ss2, dss2, t_flux, s_flux, ds_flux)
  
        R_out     = 0.5d0*(RR_tht + RR_flux)
        Z_out     = 0.5d0*(ZZ_tht + ZZ_flux)
  
        !write(*,'(A,i3,4e16.8)') ' newton (2) : ',i,errf,errx,x
  
        ifail = 0
        return
      endif
      
      if (gofast) then
        if ( (i .gt. 5) .and. (errx .gt. tol_far) ) exit
      endif
  
    enddo
  enddo
    
  !write(*,'(A,8e16.8)') ' crossing wrong exit ',x,errx,errf
  
  return

end


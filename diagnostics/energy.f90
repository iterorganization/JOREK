subroutine energy(W_mag,W_kin)
!---------------------------------------------------------------
!
!---------------------------------------------------------------
use data_structure
use gauss
use basis_at_gaussian
use phys_module
use nodes_elements

implicit none

type (type_element)      :: element
type (type_node)         :: nodes(n_vertex_max)

real*8     :: x_g(n_gauss,n_gauss),        x_s(n_gauss,n_gauss),        x_t(n_gauss,n_gauss)
real*8     :: y_g(n_gauss,n_gauss),        y_s(n_gauss,n_gauss),        y_t(n_gauss,n_gauss)
real*8     :: eq_g(n_var,n_gauss,n_gauss), eq_s(n_var,n_gauss,n_gauss), eq_t(n_var,n_gauss,n_gauss)
real*8     :: density_eq(n_gauss,n_gauss), eq_g1(n_var,n_gauss,n_gauss), Fprofile(n_gauss,n_gauss)
real*8     :: AR0_p(n_gauss,n_gauss), AZ0_p(n_gauss,n_gauss), AR0_Z, AZ0_R, A30_R, A30_Z, BR, BZ, Bp, factor

integer    :: i, j, k, in, ms, mt, iv, inode, ife, n_elements
real*8     :: W_kin(n_tor), W_mag(n_tor), xjac, BigR, wst
real*8     :: ps0_x, ps0_y, u0_x, u0_y

W_mag = 0.d0
W_kin = 0.d0

do ife =1,  element_list%n_elements

  element = element_list%element(ife)

  do iv = 1, n_vertex_max
    inode     = element%vertex(iv)
    nodes(iv) = node_list%node(inode)
  enddo

  x_g(:,:) = 0.d0;    x_s(:,:) = 0.d0;    x_t(:,:) = 0.d0;
  y_g(:,:) = 0.d0;    y_s(:,:) = 0.;      y_t(:,:) = 0.d0;
  eq_g(:,:,:) = 0.d0; eq_s(:,:,:) = 0.d0; eq_t(:,:,:) = 0.d0;
  Fprofile(:,:) = 0.d0

  do i=1,n_vertex_max
    do j=1,n_degrees
      do ms=1, n_gauss
        do mt=1, n_gauss

          x_g(ms,mt) = x_g(ms,mt) + nodes(i)%x(1,j,1) * element%size(i,j) * H(i,j,ms,mt)
          y_g(ms,mt) = y_g(ms,mt) + nodes(i)%x(1,j,2) * element%size(i,j) * H(i,j,ms,mt)

          x_s(ms,mt) = x_s(ms,mt) + nodes(i)%x(1,j,1) * element%size(i,j) * H_s(i,j,ms,mt)
          x_t(ms,mt) = x_t(ms,mt) + nodes(i)%x(1,j,1) * element%size(i,j) * H_t(i,j,ms,mt)
          y_s(ms,mt) = y_s(ms,mt) + nodes(i)%x(1,j,2) * element%size(i,j) * H_s(i,j,ms,mt)
          y_t(ms,mt) = y_t(ms,mt) + nodes(i)%x(1,j,2) * element%size(i,j) * H_t(i,j,ms,mt)

#ifdef fullmhd
          Fprofile(ms,mt) = Fprofile(ms,mt) + nodes(i)%Fprof_eq(j) * element%size(i,j) * H(i,j,ms,mt)  
#endif

        enddo
      enddo
    enddo
  enddo

  do in=1,n_tor

    eq_g(:,:,:) = 0.d0; eq_s(:,:,:) = 0.d0; eq_t(:,:,:) = 0.d0; AR0_p(:,:) = 0.d0; AZ0_p(:,:) = 0.d0

    do i=1,n_vertex_max
      do j=1,n_degrees
        do ms=1, n_gauss
          do mt=1, n_gauss

            do k=1,n_var
              eq_g(k,ms,mt)  = eq_g(k,ms,mt)  + nodes(i)%values(in,j,k) * element%size(i,j) * H(i,j,ms,mt)
              eq_s(k,ms,mt)  = eq_s(k,ms,mt)  + nodes(i)%values(in,j,k) * element%size(i,j) * H_s(i,j,ms,mt)
              eq_t(k,ms,mt)  = eq_t(k,ms,mt)  + nodes(i)%values(in,j,k) * element%size(i,j) * H_t(i,j,ms,mt)
            enddo
        
	    if (in .eq. 1) then
              density_eq(ms,mt) = abs(eq_g(5,ms,mt))
            endif

#ifdef fullmhd
            if (in .eq. 1) then
              do k = 1,n_var
                eq_g1(k,ms,mt) = eq_g(k,ms,mt) ! store the n=0 component
              enddo
            endif

            if ( mod(in,2) == 0 ) then  ! cosine
              ! in+1 is a sine, so d/dphi a cosine
              AR0_p(ms,mt) = AR0_p(ms,mt) + mode(in) * nodes(i)%values(in+1,j,var_AR) * element%size(i,j) * H(i,j,ms,mt)
              AZ0_p(ms,mt) = AZ0_p(ms,mt) + mode(in) * nodes(i)%values(in+1,j,var_AZ) * element%size(i,j) * H(i,j,ms,mt)
            elseif( mode(in) /= 0) then ! sine (for n=0 component AR0_p = AZ0_p = 0)
              ! in-1 is a cosine, so d/dphi a (-)sine
              AR0_p(ms,mt) = AR0_p(ms,mt) - mode(in) * nodes(i)%values(in-1,j,var_AR) * element%size(i,j) * H(i,j,ms,mt)
              AZ0_p(ms,mt) = AZ0_p(ms,mt) - mode(in) * nodes(i)%values(in-1,j,var_AZ) * element%size(i,j) * H(i,j,ms,mt)
            endif
#endif

          enddo
        enddo
      enddo
    enddo

!--------------------------------------------------- sum over the Gaussian integration points

    do ms=1, n_gauss

      do mt=1, n_gauss

        wst = wgauss(ms)*wgauss(mt)

        xjac = x_s(ms,mt)*y_t(ms,mt) - x_t(ms,mt)*y_s(ms,mt)
        BigR = x_g(ms,mt)

#ifdef fullmhd
        AR0_Z = ( - x_t(ms,mt) * eq_s(var_AR,ms,mt)  + x_s(ms,mt) * eq_t(var_AR,ms,mt) ) / xjac
        AZ0_R = (   y_t(ms,mt) * eq_s(var_AZ,ms,mt)  - y_s(ms,mt) * eq_t(var_AZ,ms,mt) ) / xjac
        A30_R = (   y_t(ms,mt) * eq_s(var_A3,ms,mt)  - y_s(ms,mt) * eq_t(var_A3,ms,mt) ) / xjac
        A30_Z = ( - x_t(ms,mt) * eq_s(var_A3,ms,mt)  + x_s(ms,mt) * eq_t(var_A3,ms,mt) ) / xjac


        BR     = ( A30_Z - AZ0_p(ms,mt) ) / BigR
        BZ     = ( AR0_p(ms,mt) - A30_R ) / BigR
        Bp     = ( AZ0_R - AR0_Z ) 
        if ( in == 1 ) Bp = Bp + Fprofile(ms,mt) / BigR

        factor = 2.d0 * PI * BigR * xjac * wst

        W_mag(in) = W_mag(in) + 0.5d0 * ( BR**2 + BZ**2 + Bp**2 ) * factor
        W_kin(in) = W_kin(in) + 0.5d0 * factor * eq_g1(var_rho,ms,mt)*( eq_g(var_uR,ms,mt)**2 + eq_g(var_uZ,ms,mt)**2 + eq_g(var_up,ms,mt)**2 )

        if (in /= 1) then ! a non-axisymmetric density harmonic adds to the kinetic energy when there is a nonzero n=0 velocity component
          W_kin(in) = W_kin(in) + 0.5d0 * factor * ( & 
        + eq_g(var_rho,ms,mt) *( eq_g1(var_uR,ms,mt)*eq_g(var_uR,ms,mt) + eq_g1(var_uZ,ms,mt)*eq_g(var_uZ,ms,mt) + eq_g1(var_up,ms,mt)*eq_g(var_up,ms,mt)  ) &
        + eq_g(var_rho,ms,mt) *( eq_g(var_uR,ms,mt) *eq_g1(var_uR,ms,mt)+ eq_g(var_uZ,ms,mt) *eq_g1(var_uZ,ms,mt)+ eq_g(var_up,ms,mt) *eq_g1(var_up,ms,mt) ) )
        endif

#else
        ps0_x = (   y_t(ms,mt) * eq_s(var_psi,ms,mt) - y_s(ms,mt) * eq_t(var_psi,ms,mt) ) / xjac
        ps0_y = ( - x_t(ms,mt) * eq_s(var_psi,ms,mt) + x_s(ms,mt) * eq_t(var_psi,ms,mt) ) / xjac
        u0_x  = (   y_t(ms,mt) * eq_s(var_u,  ms,mt) - y_s(ms,mt) * eq_t(var_u,ms,mt) ) / xjac
        u0_y  = ( - x_t(ms,mt) * eq_s(var_u,  ms,mt) + x_s(ms,mt) * eq_t(var_u,ms,mt) ) / xjac

        W_mag(in) = W_mag(in) +                     (ps0_x*ps0_x + ps0_y*ps0_y ) / BigR    * xjac * wst
        W_kin(in) = W_kin(in) + density_eq(ms,mt) * (u0_x*u0_x   + u0_y*u0_y)    * BigR**3 * xjac * wst
#endif

!        if (gamma /= 1.d0) then ! the internal energy density p/(gamma-1) may be absorbed in the kinetic energy
!                                  ( perhaps best is to output it as a separate energy in the future ? ) 
!          W_kin(in) = W_kin(in) + eq_g(var_rho,ms,mt)*eq_g(var_T,ms,mt) / ( gamma - 1.d0 ) * factor
!        endif


      enddo
    enddo
  enddo

enddo

do in=1,n_tor
  if (mode(in) .ne. 0) then
    W_mag(in) = 0.5d0 * W_mag(in)
    W_kin(in) = 0.5d0 * W_kin(in)
  endif
enddo
return
end

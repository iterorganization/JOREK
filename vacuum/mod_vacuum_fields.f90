!< Calculates fields created by wall and PF coil currents at arbitrary points
!< The points can be inside the plasma or in the vacuum, any xyz point but
!< DO NOT calculate the fields exactly at the STARWALL's coil/wall triangles
!< Otherwise singularities at those points may occur!
!< Additionally, wall forces are computed through the integral of the stress tensor
module mod_vacuum_fields

  use vacuum
  use vacuum_response, only: reconstruct_triangle_potentials, reconstruct_coil_potentials

  implicit none

  ! --- Scaled STARWALL wall triangle coordinates (n_tri, 3)
  real*8, allocatable :: xw_scaled(:,:), yw_scaled(:,:), zw_scaled(:,:)

   
  contains






  !< This routine calculates the total wall forces by integrating the force tensor on a closed surface
  !! outside the wall (doi:10.1088/1741-4326/aa8876)
  subroutine total_wall_forces(my_id, node_list, element_list, scale_fact, Fx, Fy, Fz)

    use constants
    use data_structure
    use mod_plasma_response,    only: plasma_fields_at_xyz    
    use phys_module,            only: F0

    implicit none

    ! --- External parameters
    integer,                   intent(in)    :: my_id
    type (type_node_list),     intent(in)    :: node_list
    type (type_element_list),  intent(in)    :: element_list   
    real*8,                    intent(in)    :: scale_fact
    real*8,                    intent(inout) :: Fx, Fy, Fz  ! --- The total force in SI units (cartesian components)

    ! --- Local parameters
    real*8               :: bx, by, bz
    real*8, allocatable  :: bx_c(:), by_c(:), bz_c(:)
    real*8, allocatable  :: bx_w(:), by_w(:), bz_w(:)
    real*8, allocatable  :: bx_p(:), by_p(:), bz_p(:)
    real*8, allocatable  ::  x_w(:),  y_w(:),  z_w(:)
    real*8               :: tri_area, B2, Bn, r1(3), r2(3), r3(3), r21(3), r32(3) 
    real*8               :: nx, ny, nz, r21_cross_r32(3), x_mid, y_mid, R_mid
    real*8               :: Bphi_x, Bphi_y
    integer              :: i

    ! --- Create a surface just outside the wall (made of triangles)
    call resize_starwall_wall(scale_fact)

    ! --- Obtain xyz points at the center of the triangles
    allocate( x_w(sr%ntri_w), y_w(sr%ntri_w), z_w(sr%ntri_w) )
    do i=1, sr%ntri_w
      x_w(i) = sum( xw_scaled(i,:) ) / 3.d0
      y_w(i) = sum( yw_scaled(i,:) ) / 3.d0
      z_w(i) = sum( zw_scaled(i,:) ) / 3.d0
    end do 

    ! --- Calculate the total magnetic field at the triangle centers
    allocate(bx_c(sr%ntri_w), by_c(sr%ntri_w), bz_c(sr%ntri_w))
    allocate(bx_w(sr%ntri_w), by_w(sr%ntri_w), bz_w(sr%ntri_w))
    allocate(bx_p(sr%ntri_w), by_p(sr%ntri_w), bz_p(sr%ntri_w))

    call coil_fields_at_xyz(my_id, x_w, y_w, z_w, bx_c, by_c, bz_c)
    call wall_fields_at_xyz(my_id, x_w, y_w, z_w, bx_w, by_w, bz_w)
    call plasma_fields_at_xyz(my_id, node_list,element_list, x_w, y_w, z_w, bx_p, by_p, bz_p)

    ! ---- Do integral of the force tensor over the triangle discretized surface
    Fx = 0.d0;    Fy = 0.d0;   Fz = 0.d0

    do i=1, sr%ntri_w

      r1(:)  = (/ xw_scaled(i,1), yw_scaled(i,1), zw_scaled(i,1) /)
      r2(:)  = (/ xw_scaled(i,2), yw_scaled(i,2), zw_scaled(i,2) /)
      r3(:)  = (/ xw_scaled(i,3), yw_scaled(i,3), zw_scaled(i,3) /)

      r21(:) = r1(:)-r2(:)
      r32(:) = r2(:)-r3(:)

      r21_cross_r32(:) = (/ r21(2)*r32(3) - r21(3)*r32(2), r21(3)*r32(1) - r21(1)*r32(3),          &
        r21(1)*r32(2) - r21(2)*r32(1) /)

      tri_area = sqrt(sum(r21_cross_r32**2.d0)) / 2.d0     

      nx = r21_cross_r32(1) / sqrt(sum(r21_cross_r32**2.d0)) 
      ny = r21_cross_r32(2) / sqrt(sum(r21_cross_r32**2.d0)) 
      nz = r21_cross_r32(3) / sqrt(sum(r21_cross_r32**2.d0)) 

      x_mid  =  sum( xw_scaled(i,:) ) / 3.d0
      y_mid  =  sum( yw_scaled(i,:) ) / 3.d0
      R_mid  = sqrt( x_mid**2.d0 + y_mid**2.d0 )

      Bphi_x = F0/R_mid * (  y_mid/R_mid ) ! Bphi * (-sin phi) 
      Bphi_y = F0/R_mid * ( -x_mid/R_mid ) ! Bphi * (-cos phi)

      bx = bx_p(i) + bx_c(i) + bx_w(i) + Bphi_x
      by = by_p(i) + by_c(i) + by_w(i) + Bphi_y
      bz = bz_p(i) + bz_c(i) + bz_w(i)

      Bn = bx*nx + by*ny + bz*nz
      B2 = bx**2.d0 + by**2.d0 + bz**2.d0 

      Fx = Fx + (Bn*bx - B2*0.5d0*nx) * tri_area / mu_zero
      Fy = Fy + (Bn*by - B2*0.5d0*ny) * tri_area / mu_zero
      Fz = Fz + (Bn*bz - B2*0.5d0*nz) * tri_area / mu_zero

    end do 

    ! --- Clean-up
    deallocate(bx_c, by_c, bz_c)
    deallocate(bx_w, by_w, bz_w)
    deallocate(bx_p, by_p, bz_p)
    deallocate( x_w,  y_w,  z_w)  
    deallocate(xw_scaled, yw_scaled, zw_scaled)

  end subroutine total_wall_forces





 
  !< This routine calculates the fields created by the STARWALL wall at given cartesian coordinates
  !< DO NOT calculate the fields exactly at the STARWALL's coil/wall triangles
  !< Otherwise singularities at those points may occur!
  subroutine wall_fields_at_xyz(my_id,x,y,z,bx,by,bz)

    use constants

    implicit none

    ! --- External parameters
    integer, intent(in)     :: my_id
    real*8,  intent(in)     :: x(:), y(:), z(:)     ! Points where fields are calculated
    real*8,  intent(inout)  :: bx(:), by(:), bz(:)

    ! --- Local parameters
    real*8,  allocatable    :: tripot_w(:)
    real*8,  allocatable    :: phi_w(:,:), x_w(:,:), y_w(:,:), z_w(:,:)
    real*8                  :: Iw_net_tor
    integer                 :: i, j, ipot
     
    call reconstruct_triangle_potentials(tripot_w, wall_curr, my_id, Iw_net_tor)

    if (sr%file_version < 5) then
      write(*,*) 'wall_fields_at_xyz not available for STARWALL reponse file version < 5'
      STOP
    endif

    allocate(phi_w(sr%ntri_w,3))
    allocate(x_w(sr%ntri_w,3), y_w(sr%ntri_w,3), z_w(sr%ntri_w,3))

    ! --- Total wall triangle potentials (including net currents)
    do i = 1, sr%ntri_w
      do j = 1, 3
        ipot       = sr%jpot_w(i,j)
        phi_w(i,j) = tripot_w(ipot) + Iw_net_tor*sr%phi0_w(i,j) 
        x_w(i,j)   = sr%xyzpot_w(ipot,1)
        y_w(i,j)   = sr%xyzpot_w(ipot,2)
        z_w(i,j)   = sr%xyzpot_w(ipot,3)
      end do
    end do

    phi_w   = phi_w / mu_zero  ! Wall current potentials in Amperes

    call triang_fields_at_xyz(my_id,x,y,z,x_w,y_w,z_w,phi_w,bx,by,bz)

    deallocate(phi_w, x_w, y_w, z_w)

  end subroutine wall_fields_at_xyz





  !< This routine calculates the fields created by the STARWALL coils at given cartesian coordinates
  !< DO NOT calculate the fields exactly at the STARWALL's coil/wall triangles
  !< Otherwise singularities at those points may occur!
  subroutine coil_fields_at_xyz(my_id,x,y,z,bx,by,bz, icoil)

    use constants

    implicit none

    ! --- External parameters
    integer,           intent(in)     :: my_id
    real*8,            intent(in)     :: x(:), y(:), z(:)     ! Points where fields are calculated
    real*8,            intent(inout)  :: bx(:), by(:), bz(:)
    integer, optional, intent(in)     :: icoil ! If present, gets fields only from coil number "icoil"

    ! --- Local parameters
    real*8,  allocatable    :: pot_c(:)
    real*8,  allocatable    :: phi_c(:,:)
    integer                 :: i, j, ipot, ntri_c, i_c
    integer                 :: i_tri_start, i_tri_end

    if (sr%ncoil < 1) then
      write(*,*) 'coil_fields_at_xyz needs STARWALL coils (ncoils < 1) detected'
      return
    endif
    
    call reconstruct_coil_potentials(pot_c, wall_curr, my_id)

    allocate(phi_c(sr%ntri_c,3))

    ! --- Distribute coil currents over the coil triangles
    i_tri_start = 1

    do i_c=1, sr%ncoil

      i_tri_end = i_tri_start + sr%jtri_c(i_c) - 1

      ! --- When icoil is given, set the other coil currents to 0
      if (present(icoil)) then
        if (i_c /= icoil) pot_c(i_c) = 0.d0
      endif

      do i = i_tri_start, i_tri_end
        do j = 1, 3
          phi_c(i,j) = sr%phi_coil(i,j) * pot_c(i_c)
        end do
      end do

      i_tri_start = i_tri_end + 1

    enddo

    phi_c   = phi_c / mu_zero  ! Wall current potentials in Amperes

    call triang_fields_at_xyz(my_id,x,y,z,sr%x_coil,sr%y_coil,sr%z_coil,phi_c,bx,by,bz)

    deallocate(phi_c)

  end subroutine coil_fields_at_xyz






  !< This routine calculates the fields produced by currents flowing on a set of triangles
  !! at given xyz points (taken from STARWALL)
  subroutine triang_fields_at_xyz(my_id,x,y,z,x_tri,y_tri,z_tri,phi_tri,bx,by,bz)

    use constants
    use mpi_mod
    !$ use omp_lib

    implicit none

    ! --- External parameters
    integer, intent(in)     :: my_id
    real*8,  intent(in)     :: x(:), y(:), z(:) ! Points where fields are calculated
    real*8,  intent(in)     :: x_tri(:,:), y_tri(:,:), z_tri(:,:), phi_tri(:,:)
    real*8,  intent(inout)  :: bx(:), by(:), bz(:)

    ! --- Local parameters
    integer :: ierr, n_cpu, k_delta, k_min, k_max
    integer :: i, j, k, np, ntri, omp_nthreads, omp_tid
    real*8  :: s1,s2,s3                                    &
              ,d221,d232,d213,al1,al2,al3                  &
              ,ata1,ata2,ata3,at                           &
              ,s21,s22,s23,dp1,dm1,dp2,dm2,dp3,dm3         &
              ,ap1,am1,ap2,am2,ap3,am3                     &
              ,h,ar1,ar2,ar3                               &
              ,x21,y21,z21,x32,y32,z32,x13,y13,z13,vx,vy,vz&
              ,tx1,ty1,tz1,tx2,ty2,tz2,tx3,ty3,tz3         &
              ,nx,ny,nz,pi41,area,d21,d32,d13,jx,jy,jz     &
              ,dep1,dep2,dep3,dem1,dem2,dem3
    real*8  :: Rcent, jphi, cosx, siny
    real*8  :: x1,y1,z1,x2,y2,z2,x3,y3,z3,sn
    real*8,  allocatable  :: bx_tmp(:), by_tmp(:), bz_tmp(:)

    np   = size(x,1)
    ntri = size(x_tri,1)

    pi41 = 0.125d0/asin(1.d0)

    allocate(bx_tmp(np), by_tmp(np), bz_tmp(np))

    bx     = 0.d0;  by     = 0.d0;  bz     = 0.d0;
    bx_tmp = 0.d0;  by_tmp = 0.d0;  bz_tmp = 0.d0;

    ! --- MPI initialization
    call MPI_COMM_SIZE(MPI_COMM_WORLD, n_cpu, ierr) ! number of MPI procs
    n_cpu = max(n_cpu,1)

    k_delta = ceiling(float(ntri) / n_cpu)
    k_min   =      my_id     * k_delta + 1
    k_max   = min((my_id +1) * k_delta, ntri)
    ! --- OpenMP parallelization of given points loop
    !$omp parallel default(none)                                                            &
    !$omp   shared(np,x_tri,y_tri,z_tri,x,y,z, k_min, k_max,pi41,phi_tri,                   &
    !$omp          bx_tmp, by_tmp, bz_tmp)                                                  &
    !$omp   private(i,x1,y1,z1,x2,y2,z2,x3,y3,z3,sn,h,s21,s22,s23,s1,s2,s3,al1,al2,al3,     &
    !$omp           ar1,ar2,ar3,dp1,dp2,dp3,dm1,dm2,dm3,ap1,ap2,ap3,dep1,dep2,dep3,         &
    !$omp           d21,d32,d13,tx2,ty2,tz2,tx3,ty3,tz3,d221, d232,d213,area,Rcent,cosx,    &
    !$omp           nx,ny,nz, jx,jy,jz, tx1,ty1,tz1,k,x21,y21,z21,x32,y32,z32,x13,y13,z13,  &
    !$omp           dem1,dem2,dem3,am1,am2,am3,ata1,ata2,ata3,at,vx,vy,vz,siny,jphi,        &
    !$omp           omp_nthreads,omp_tid)
    
#ifdef OPENMP
    omp_nthreads = omp_get_num_threads()
    omp_tid      = omp_get_thread_num()
#else
    omp_nthreads = 1
    omp_tid      = 0
#endif

    !$omp do reduction(+:bx_tmp, by_tmp, bz_tmp)     
 
    do k=k_min, k_max    ! --- integral over wall triangles

     !--- only use toroidal current, projection
      x21   = x_tri(k,2) - x_tri(k,1)
      y21   = y_tri(k,2) - y_tri(k,1)
      z21   = z_tri(k,2) - z_tri(k,1)
      x32   = x_tri(k,3) - x_tri(k,2)
      y32   = y_tri(k,3) - y_tri(k,2)
      z32   = z_tri(k,3) - z_tri(k,2)
      x13   = x_tri(k,1) - x_tri(k,3)
      y13   = y_tri(k,1) - y_tri(k,3)
      z13   = z_tri(k,1) - z_tri(k,3)
      d221  = x21**2+y21**2+z21**2
      d232  = x32**2+y32**2+z32**2
      d213  = x13**2+y13**2+z13**2
      d21   = sqrt(d221)
      d32   = sqrt(d232)
      d13   = sqrt(d213)
      nx    = -y21*z13 + z21*y13
      ny    = -z21*x13 + x21*z13
      nz    = -x21*y13 + y21*x13
      area  = 1./sqrt(nx*nx+ny*ny+nz*nz)

      ! --- The triangle current density
      jx = (x32*phi_tri(k,1)+x13*phi_tri(k,2)+x21*phi_tri(k,3))*area*pi41 
      jy = (y32*phi_tri(k,1)+y13*phi_tri(k,2)+y21*phi_tri(k,3))*area*pi41 
      jz = (z32*phi_tri(k,1)+z13*phi_tri(k,2)+z21*phi_tri(k,3))*area*pi41 

      ! --- Only use toroidal current (RMHD)
      Rcent = sqrt( (sum(x_tri(k,:))/3.d0)**2.d0 + (sum(y_tri(k,:))/3.d0)**2.d0 )

      cosx  =  sum(x_tri(k,:))/(3.d0*Rcent)
      siny  = -sum(y_tri(k,:))/(3.d0*Rcent)

      jphi  = -jx*siny - jy*cosx

      jx    = -jphi * siny
      jy    = -jphi * cosx
      jz    = 0.d0

      nx    = nx*area
      ny    = ny*area
      nz    = nz*area
      tx1   = (y32*nz-z32*ny)
      ty1   = (z32*nx-x32*nz)
      tz1   = (x32*ny-y32*nx)
      tx2   = (y13*nz-z13*ny)
      ty2   = (z13*nx-x13*nz)
      tz2   = (x13*ny-y13*nx)
      tx3   = (y21*nz-z21*ny)
      ty3   = (z21*nx-x21*nz)
      tz3   = (x21*ny-y21*nx)


      do i=1, np      ! --- go over given points

        x1    = x_tri(k,1) - x(i)
        y1    = y_tri(k,1) - y(i)
        z1    = z_tri(k,1) - z(i)
        x2    = x_tri(k,2) - x(i)
        y2    = y_tri(k,2) - y(i)
        z2    = z_tri(k,2) - z(i)
        x3    = x_tri(k,3) - x(i)
        y3    = y_tri(k,3) - y(i)
        z3    = z_tri(k,3) - z(i)
        sn    = nx*x1+ny*y1+nz*z1
        h     = abs(sn)
        s21   = x1**2+y1**2+z1**2
        s22   = x2**2+y2**2+z2**2
        s23   = x3**2+y3**2+z3**2
        s1    = sqrt(s21)
        s2    = sqrt(s22)
        s3    = sqrt(s23)
        al1   = alog((s2+s1+d21)/(s1+s2-d21))
        al2   = alog((s3+s2+d32)/(s3+s2-d32))
        al3   = alog((s1+s3+d13)/(s1+s3-d13))
        ar1   = x1*tx3+y1*ty3+z1*tz3
        ar2   = x2*tx1+y2*ty1+z2*tz1
        ar3   = x3*tx2+y3*ty2+z3*tz2
        dp1   = .5*(s22-s21+d221)
        dp2   = .5*(s23-s22+d232)
        dp3   = .5*(s21-s23+d213)
        dm1   = dp1-d221
        dm2   = dp2-d232
        dm3   = dp3-d213
        ap1   = ar1*dp1
        dep1  = ar1**2+h*d221*(h+s2)
        ap2   = ar2*dp2
        dep2  = ar2**2+h*d232*(h+s3)
        ap3   = ar3*dp3
        dep3  = ar3**2+h*d213*(h+s1)
        am1   = ar1*dm1
        dem1  = ar1**2+h*d221*(h+s1)
        am2   = ar2*dm2
        dem2  = ar2**2+h*d232*(h+s2)
        am3   = ar3*dm3
        dem3  = ar3**2+h*d213*(h+s3)
        ata1  = atan2(ap1*dem1-am1*dep1,dep1*dem1+ap1*am1)
        ata2  = atan2(ap2*dem2-am2*dep2,dep2*dem2+ap2*am2)
        ata3  = atan2(ap3*dem3-am3*dep3,dep3*dem3+ap3*am3)
        at    = sign(1.,sn)*(ata1+ata2+ata3)
        vx    = -nx*at + al1*tx3/d21+al2*tx1/d32+al3*tx2/d13
        vy    = -ny*at + al1*ty3/d21+al2*ty1/d32+al3*ty2/d13
        vz    = -nz*at + al1*tz3/d21+al2*tz1/d32+al3*tz2/d13

        bx_tmp(i) = bx_tmp(i) + vy*jz-vz*jy
        by_tmp(i) = by_tmp(i) + vz*jx-vx*jz
        bz_tmp(i) = bz_tmp(i) + vx*jy-vy*jx

      enddo
    enddo
    !$omp end do
    !$omp end parallel


    call MPI_AllReduce(bx_tmp,bx,np,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    call MPI_AllReduce(by_tmp,by,np,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    call MPI_AllReduce(bz_tmp,bz,np,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

    bx = -bx * mu_zero
    by = -by * mu_zero
    bz = -bz * mu_zero

    deallocate(bx_tmp, by_tmp, bz_tmp)   

  end subroutine triang_fields_at_xyz






  !< This routine re-calculates the STARWALL wall from wall Fourier harmonics and
  !! resizes it acording to a scale factor. This is useful to create a surface
  !! just outside the wall to perform wall forces integration. Note that using the
  !! same wall triangles would lead to singularities for the wall field calculation
  subroutine resize_starwall_wall(scale_fact)

    implicit none

    ! --- External parameters
    real*8, intent(in) :: scale_fact

    ! --- Local parameters
    real*8              :: pi2,fnu,alu,alv
    real*8              :: cm,cn,cov,siv,cou,siu,cop,sip,co,si
    real*8, allocatable ::  r_w(:),  x_w(:),  y_w(:),  z_w(:)
    real*8, allocatable :: rc_w(:), rs_w(:), zc_w(:), zs_w(:)
    integer             :: nwuv, i, j, kv, ku

    if (sr%file_version < 5) then
      write(*,*) 'STARWALL version not suitable for wall resizing'
      STOP
    endif
    if (sr%iwall /= 1) then
      write(*,*) 'Wall must be given in Fourier harmonics for resizing (iwall=1)'
      STOP
    endif

    nwuv  = sr%nwu*sr%nwv ! total number of wall nodes

    allocate (x_w(nwuv),y_w(nwuv),z_w(nwuv),r_w(nwuv))
    allocate (rc_w(sr%mn_w),rs_w(sr%mn_w),zc_w(sr%mn_w),zs_w(sr%mn_w))
    if (allocated(xw_scaled)) deallocate(xw_scaled)
    if (allocated(yw_scaled)) deallocate(yw_scaled)
    if (allocated(zw_scaled)) deallocate(zw_scaled)
    allocate (xw_scaled(sr%ntri_w,3),yw_scaled(sr%ntri_w,3),zw_scaled(sr%ntri_w,3))

    ! --- Resize the wall, wall construction copied from STARWALL
    pi2  =4.*asin(1.)
    fnu  = 1./float(sr%nwu) 
    alu  = pi2*fnu 
    alv  = pi2/float(sr%nwv)
    z_w = 0.
    r_w = 0.

    rc_w(1) =sr%rc_w(1)
    rs_w(1) =sr%rs_w(1)
    zc_w(1) =sr%zc_w(1)
    zs_w(1) =sr%zs_w(1)

    rc_w(2:sr%mn_w) =sr%rc_w(2:sr%mn_w)*scale_fact
    rs_w(2:sr%mn_w) =sr%rs_w(2:sr%mn_w)*scale_fact
    zc_w(2:sr%mn_w) =sr%zc_w(2:sr%mn_w)*scale_fact
    zs_w(2:sr%mn_w) =sr%zs_w(2:sr%mn_w)*scale_fact

    do  j =  1, sr%mn_w
      cm = sr%m_w(j)*pi2
      cn = sr%n_w_fourier(j)*pi2
      do  kv=1, sr%nwv
        cov = cos(alv*sr%n_w_fourier(j)*(kv-1))
        siv = sin(alv*sr%n_w_fourier(j)*(kv-1))
        do  ku=1,sr%nwu
          i = ku+sr%nwu*(kv-1)
          cou = cos(alu*sr%m_w(j)*(ku-1))
          siu = sin(alu*sr%m_w(j)*(ku-1))
    
          cop = cou*cov-siu*siv
          sip = siu*cov+cou*siv
    
          r_w(i) = r_w(i) + rs_w(j)*sip + rc_w(j)*cop
          z_w(i) = z_w(i) + zs_w(j)*sip + zc_w(j)*cop 
        end do
      end do
    end do
    
    do  kv = 1, sr%nwv
      co   = cos(alv*(kv-1)) 
      si   = sin(alv*(kv-1)) 
      do  ku = 1,sr%nwu
        i      = sr%nwu*(kv-1)+ku
        x_w(i)   =   co * r_w(i)
        y_w(i)   =   si * r_w(i)
      end do
    end do

    ! --- Reconstruct scaled wall triangle coordinates
    do i=1,2*nwuv
      xw_scaled(i,1) = x_w(sr%jpot_w(i,1))
      xw_scaled(i,2) = x_w(sr%jpot_w(i,2))
      xw_scaled(i,3) = x_w(sr%jpot_w(i,3))

      yw_scaled(i,1) = y_w(sr%jpot_w(i,1))
      yw_scaled(i,2) = y_w(sr%jpot_w(i,2))
      yw_scaled(i,3) = y_w(sr%jpot_w(i,3))

      zw_scaled(i,1) = z_w(sr%jpot_w(i,1))
      zw_scaled(i,2) = z_w(sr%jpot_w(i,2))
      zw_scaled(i,3) = z_w(sr%jpot_w(i,3))
    end do

    deallocate( x_w,  y_w,  z_w,  r_w)
    deallocate(rc_w, rs_w, zc_w, zs_w)

  end subroutine resize_starwall_wall


end module mod_vacuum_fields


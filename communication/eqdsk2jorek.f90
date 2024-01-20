      program eqdsk
!--------------------------------------------------------------------
! little program to construct an input file for jorek out of data
! in a eqdsk file
!                         Guido Huysmans,          date : 14-12-2010
!
! Some documentation can be found here: https://www.jorek.eu/wiki/doku.php?id=eqdsk2jorek.f90
!--------------------------------------------------------------------
implicit none

real*8,allocatable :: psi(:),p(:),f(:),q(:),rlim(:),zlim(:),rbnd(:), zbnd(:)
real*8,allocatable :: apsi(:),bpsi(:),cpsi(:),dpsi(:)
real*8,allocatable :: ap(:),bp(:),cp(:),dp(:)
real*8,allocatable :: af(:),bf(:),cf(:),df(:)
real*8,allocatable :: radius(:),theta(:),rad(:)
real*8,allocatable :: ar(:),br(:),cr(:),dr(:)
real*8,allocatable :: dpr(:),df2(:),dg(:),work(:),psirz(:,:)
real*8,allocatable :: xx(:),yy(:),zc(:), r_bnd(:), z_bnd(:), psi_bnd(:)
real*8,allocatable :: df2_ext(:),rho_ext(:),T_ext(:),psi_ext(:),p_ext(:)
real*8,allocatable :: tx(:),ty(:),c(:,:),wrk(:)
integer,allocatable :: iwrk(:)
real*8             :: angle, ellip, tria_u, tria_l, quad_u, quad_l, r0, z0, a0, PI
real*8             :: dummy(3), xdim,zdim,rzero,rgrid1,zmid,rmaxis,zmaxis,ssimag,ssibry,bcentr,a_minor
real*8             :: xip,xdum1,xdum2,xdum3,xdum4,xdum5
real*8             :: psi_sep, sig_sep, tanh1, zmu0, zn0, zmd, rho_bnd, T_bnd
real*8             :: xb ,xe, yb, ye, smth, fp, fout
real*8             :: B_scale, I_scale, R_scale, F_axis, factor, dfactor
real*8             :: ellip_in,tria_up_in,tria_low_in,quad_up_in,quad_low_in,r0_in,z0_in,a0_in
integer            :: mx,my,kx,ky,nxest,nyest,lwrk,kwrk,ier,iopt,nx,ny, i1, j1,iostatus, str_id,ierr
integer            :: nr, nz, n_psi, nbbs, limitr, i,j, nc, n_tht, n_sol, n_ext, ivtk, n_tht_in
character          :: AA*52, tokamak_name*50,boundary_type*100
character          :: buffer*80, lf*1, str1*12, str2*24, string_in*250,eqdsk_string_r_min*250

namelist /eqdsk2jorek_params/ tokamak_name,boundary_type,ellip_in,tria_up_in,&
                              tria_low_in,quad_up_in,quad_low_in,n_tht_in,r0_in,&
                              z0_in,a0_in,B_scale,I_scale,R_scale,smth,eqdsk_string_r_min

!----------------------------- read eqdsk file -----------
B_scale = 1.d0/1.d0  ! scaling factor for the vacuum toroidal field 
I_scale = 1.d0/1.d0  ! scaling factor for the toroidal current
R_scale = 1.d0/1.d0  ! scaling factor for the space coordinates 
!> Define the defaults JOREK boundaries for a specific tokamak: 'ITER' (default), 'JET', 
!> 'DIII-D' or define a boundary via user's inputs: 'USER_DEFINED' 
tokamak_name  = 'ITER'
!> Define which default boundary to be used (major radius, vertical position and 
!> minor radius are rescaled by the R_scale factor):
!> NOTE: the JOREK boundary is computed using:
!>   R = R_axis + r_minor*cos(theta+triangularity*sin(theta)+quadrangularity*sin(2*theta))
!>   Z = Z_axis + r_minor*ellipticity*sin(theta)
!> ITER: CLOSE_WALL_FIT, OUTSIDE_WALL (default)
!> JET: OUTSIDE_WALL, OUTSIDE_WALL_SHORT_LEG, CIRCULAR
!> DIII-D: OUTSIDE_WALL, NIMROD_M3DC1
!>
!> eqdsk_string_r_min: string of the EQDSK file identifying the plasma minor radius
!>   default value: 'MINOR RADIUS -> A [m]'
boundary_type = 'OUTSIDE_WALL'
eqdsk_string_r_min = 'MINOR RADIUS -> A [m]'
ellip_in    = 1.d0; tria_up_in  = 0.d0; tria_low_in = 0.d0;
quad_up_in  = 0.d0; quad_low_in = 0.d0; n_tht_in    = 259;
r0_in       = 3.d0; z0_in       = 0.d0; a0_in       = 1.d0;
smth = 1.d-6 ! Controls the tradeoff between closeness of fit and smoothness of fit. When too small, can lead to noise pick-up. When too large, can lead to inaccurate fit.
             ! May need hand tuning, based on a visual inspection of the output.
             ! For more details, see the documentation of regrid.f in libdierckx or the "Tips and Tricks" section of the Wiki page https://www.jorek.eu/wiki/doku.php?id=eqdsk2jorek.f90. 
write(*,*) ' EQDSK to JOREK2 '

! --- Read parameters from namelist file 'eqdsk2jorek.nml' if it exists
open(42, file='eqdsk2jorek.nml', action='read', status='old', iostat=ierr)
if ( ierr == 0 ) then
  write(*,*) 'Reading parameters from eqdsk2jorek.nml namelist.'
  read(42,eqdsk2jorek_params)
  close(42)
end if 

write(*,*) '   Tokamak = ', tokamak_name
write(*,*) '   Boundary type = ',boundary_type
write(*,*) '   EQDSK minor radius string = ',eqdsk_string_r_min

read(5,'(A52,2i4)') AA,nr,nz

write(*,*) AA
write(*,'(A,2i5)') ' nr, nz : ',nr,nz

read(5,'(5e16.9)') xdim,zdim,rzero,rgrid1,zmid
read(5,'(5e16.9)') rmaxis,zmaxis,ssimag,ssibry,bcentr
read(5,'(5e16.9)') xip,ssimag,xdum1,rmaxis,xdum2
read(5,'(5e16.9)') zmaxis,xdum3,ssibry,xdum4,xdum5

write(*,'(A,2f10.5,A)') ' xdim,  zdim : ',xdim,zdim,' m'
write(*,'(A,2f10.5,A)') ' rzero, zmid : ',rzero,zmid, ' m'
write(*,'(A,f10.5,A)')  ' xip         : ',xip/1e6,' MA'
write(*,'(A,f10.5,A)')  ' zmaxis      : ',zmaxis,' m'
write(*,'(A,f10.5,A)')  ' Bvac        : ',bcentr,' T'
write(*,'(A,3f10.5,A)') ' psi (axis,bnd) :',ssimag,ssibry,ssibry-ssimag,' Wb'

write(*,*) '   reading profiles'
  
n_psi=nr
allocate(f(n_psi),p(n_psi),df2(n_psi),dpr(n_psi),psirz(nr,nz),q(n_psi))

read(5,'(5e16.9)') (f(i),i=1,n_psi)
read(5,'(5e16.9)') (p(i),i=1,n_psi)
read(5,'(5e16.9)') (df2(i),i=1,n_psi)
read(5,'(5e16.9)') (dpr(i),i=1,n_psi)
read(5,'(5e16.9)') ((psirz(i,j),i=1,nr),j=1,nz)
read(5,'(5e16.9)') (q(i),i=1,n_psi)

write(*,*) '   reading limiter'

read(5,*)  nbbs,limitr
allocate(rbnd(nbbs),zbnd(nbbs))
read(5,'(5e16.9)') (rbnd(i),zbnd(i),i=1,nbbs)
allocate(rlim(limitr),zlim(limitr))
read(5,'(5e16.9)') (rlim(i),zlim(i),i=1,limitr)

! ------------------- find and read the plasma minor radius
iostatus = 0; str_id = 0;string_in = ''; 
do while(.true.)
  read(5,'(A)',IOSTAT=iostatus) string_in
  str_id = index(trim(string_in),trim(eqdsk_string_r_min));
  if(str_id.ne.0 .or. iostatus.ne.0) exit
enddo
string_in = trim(string_in(str_id+len(trim(eqdsk_string_r_min))+1:len(string_in)))
read(string_in,fmt=*) a_minor
if(iostatus.ne.0) write(*,*) '!!!! WARNING SOMETHING WENT WRONG READING THE PLASMA MINOR RADIUS !!!!' 
write(*,*),"a_minor:  ",a_minor,' m'

write(*,*) ' done reading'

!=============== scaling of equilibrium with B
bcentr = bcentr * B_scale
p      = p      * B_scale**2
dpr    = dpr    * B_scale
F      = F      * B_scale
dF2    = dF2    * B_scale 
psirz  = psirz  * B_scale 
xip    = xip    * B_scale

!=============== scaling of equilibrium with space dimension
xip    = xip    * R_scale
p      = p
dpr    = dpr    / R_scale**2
F      = F      * R_scale
dF2    = dF2
psirz  = psirz  * R_scale**2

rgrid1 = rgrid1 * R_scale
rzero  = rzero  * R_scale
xdim   = xdim   * R_scale
zmid   = zmid   * R_scale
zdim   = Zdim   * R_scale

!=============== scaling of equilibrium with plasma current
xip    = xip    * I_scale
p      = p      * I_scale**2
dpr    = dpr    * I_scale
dF2   = dF2     * I_scale
psirz  = psirz  * I_scale

F_axis = f(1)

do i=1, nr
  factor  = sqrt(1.d0 + F_axis**2/F(i)**2 * (1.d0/I_scale**2 - 1.d0))
  dfactor = -1.d0/(factor) * (1.d0/I_scale**2 - 1.d0) * F_axis**2/F(i)**4 * dF2(i)
  q(i)   = factor * q(i)
  F(i)   = factor * I_scale * F(i)
enddo

if ((R_scale .ne. 1.d0) .or. (B_scale .ne. 1.d0) .or. (I_scale .ne. 1.d0)) then
  write(*,'(A)')               '********************************************************'
  write(*,'(A)')               '  Equilibrium scaling : '
  write(*,'(A,f8.4,A,f8.4,A)') '    R_scale : ',R_scale,'    major radius   :',rzero,' [m]'
  write(*,'(A,f8.4,A,f8.4,A)') '    B_scale : ',B_scale,'    vacuum field   :',bcentr,' [T]'
  write(*,'(A,f8.4,A,f8.4,A)') '    I_scale : ',I_scale,'    plasma current :',xip/1d6,' [MA]'
  write(*,'(A)')               '********************************************************'
endif

allocate(xx(nr),yy(nz),psi(n_psi))
do i=1,n_psi
  psi(i) = real(i-1)/real(n_psi-1)
enddo     
do i=1,nr
  xx(i) = rgrid1 + xdim*real(i-1)/real(nr-1)
enddo     
do i=1,nz
  yy(i) = zmid + zdim*(real(i-1)/real(nz-1)-0.5)
enddo     

if (tokamak_name == 'ITER') then

  !--------------------close fit to ITER wall
  if (boundary_type == 'CLOSE_WALL_FIT') then
    ellip  = 2.0
    tria_u = 0.55
    tria_l = 0.65
    quad_u = -0.1
    quad_l = 0.15
    n_tht  = 257
    r0     = 6.2  * R_scale
    z0     = 0.1  * R_scale
    a0     = 2.25 * R_scale

  !-------------------- contour outside ITER wall
  else if(boundary_type == 'OUTSIDE_WALL') then
    ellip  = 2.1
    tria_u = 0.58
    tria_l = 0.65
    quad_u = -0.12
    quad_l = -0.
    n_tht  = 257
    r0     = 6.2   * R_scale
    z0     = -0.05 * R_scale
    a0     = 2.34  * R_scale
  else
    write(*,*) 'JOREK boundary not or wrongly specified, stopping' 
    stop
  endif

else if (tokamak_name == 'JET') then
  
  !-------------------- contour outside JET wall
  ! blue contour in https://www.jorek.eu/wiki/doku.php?id=eqdsk2jorek.f90
  if(boundary_type == 'OUTSIDE_WALL') then
    ellip  = 1.85
    tria_u = 0.4
    tria_l = 0.4
    quad_u = -0.2
    quad_l = -0.2
    n_tht  = 257
    r0     = 2.9  * R_scale
    z0     = 0.1  * R_scale
    a0     = 1.08 * R_scale

  !-------------------- contour to avoid too long divertor legs
  else if(boundary_type == 'OUTSIDE_WALL_SHORT_LEG') then
    ellip  = 1.7
    tria_u = 0.4
    tria_l = 0.4
    quad_u = -0.4
    quad_l = -0.2
    n_tht  = 257
    r0     = 2.85 * R_scale
    z0     = 0.15 * R_scale
    a0     = 1.1  * R_scale

  !-------------------- try circular plasmas
  else if(boundary_type == 'CIRCULAR') then
    ellip  = 1.
    tria_u = 0.
    tria_l = 0.
    quad_u = 0.
    quad_l = 0.
    n_tht  = 257
    r0     = rmaxis * R_scale
    z0     = zmaxis * R_scale
    a0     = a_minor * R_scale
  else
    write(*,*) 'JOREK boundary not or wrongly specified, stopping' 
    stop
  end if

else if (tokamak_name == 'DIII-D') then

  !-------------------- contour outside DIII-D wall
  if(boundary_type == 'OUTSIDE_WALL') then
    ellip  = 1.85
    tria_u = 0.4
    tria_l = 0.4
    quad_u = -0.2
    quad_l = -0.2
    n_tht  = 257
    r0     = 1.7 * R_scale
    z0     = 0.  * R_scale
    a0     = 0.7 * R_scale
  
  !-------------------- Atomic physics JOREK/NIMROD/M3D-C1 benchmark case (paper by B. Lyons)
  else if(boundary_type == 'NIMROD_M3DC1') then
    ellip  = 1.35/0.7
    tria_u = 0.3
    tria_l = 0.3
    quad_u = 0.
    quad_l = 0.
    n_tht   = 257
    r0     = 1.7 * R_scale
    z0     = 0.  * R_scale
    a0     = 0.7 * R_scale
  else
    write(*,*) 'JOREK boundary not or wrongly specified, stopping' 
    stop
  end if

else if(tokamak_name == 'USER_DEFINED') then
    ellip  = ellip_in
    tria_u = tria_up_in
    tria_l = tria_low_in
    quad_u = quad_up_in
    quad_l = quad_low_in
    n_tht  = n_tht_in
    r0     = r0_in * R_scale
    z0     = z0_in * R_scale
    a0     = a0_in * R_scale

else

  write(*,*) 'Tokamak name not or wrongly specified, stopping'
  stop

end if  

write(*,*) 'JOREK boundary parameters'
write(*,*) 'ellipticity: ',ellip
write(*,*) 'upper and lower triangularity: ',tria_u,' ',tria_l
write(*,*) 'upper and lower quadrangularity: ',quad_u,' ',quad_l
write(*,*) 'axis position, R = ',r0,'[m] Z =',z0,'[m]'
write(*,*) 'minor radius, a = ',a0,'[m]'
write(*,*) 'N# poloidal angles: ',n_tht
  
PI = 2.d0 * asin(1.d0)

allocate(r_bnd(n_tht),z_bnd(n_tht),psi_bnd(n_tht))

!--------------------------------- interpolate flux using Dierckx spline routine
iopt= 0 
mx = nr
my = nz
xb = xx(1)
xe = xx(nr)
yb = yy(1)
ye = yy(nz)
kx = 3
ky = 3
nxest = 3*nr/4 ! Upper bound for the number of knots used for the splines. We set it a bit smaller than nr to test the quality of the fit.
nyest = 3*nz/4
lwrk  = 4+nxest*(my+2*kx+5)+nyest*(2*ky+5)+mx*(kx+1)+my*(ky+1)+my+nxest
kwrk  = 3+mx+my+nxest+nyest
write(*,*) ' Interpolation smoothing parameter = ',smth

allocate(tx(nxest),ty(nyest),c(nxest,nyest),wrk(lwrk),iwrk(kwrk))

call regrid(iopt,mx,xx,my,yy,transpose(psirz),xb,xe,yb,ye,kx,ky,smth,nxest,nyest,nx,tx,ny,ty,c,fp,wrk,lwrk,iwrk,kwrk,ier)

if (ier > 0) then
  write(*,*) '!!!!! WARNING: Problem with the Dierckx spline interpolation !!!!!'
  write(*,*) '!!!!! You may need to tune smth and/or nxest and nyest.      !!!!!'
end if 
write(*,*) ' Dierckx ier   : ',ier
write(*,*) ' Dierckx fp    : ',fp
write(*,*) ' Dierckx nx,ny : ',nx,ny
if (ier > 0) then
  write(*,*) '!!!!! Exiting                                                !!!!!'
  stop
end if 

lwrk = mx*(kx+1)+my*(ky+1)
kwrk = mx+my
deallocate(wrk,iwrk)
allocate(wrk(lwrk),iwrk(kwrk))

!do i=1,nr
!  call bispev(tx,nx,ty,ny,c,kx,ky,xx(i),1,yy(nz/2),1,fout,wrk,lwrk,iwrk,kwrk,ier)
!  write(*,'(4e16.8,i3)') xx(i),yy(nz/2),fout,psirz(i,nz/2),ier
!enddo

do i=1,n_tht/2
  angle = 2.d0 * PI * float(i-1)/float(n_tht-1)
  r_bnd(i) = r0 + a0 * cos(angle + tria_u*sin(angle) + quad_u*sin(2.d0*angle))
  z_bnd(i) = z0 + a0 * ellip * sin(angle)
  call bispev(tx,nx,ty,ny,c,kx,ky,r_bnd(i),1,z_bnd(i),1,psi_bnd(i),wrk,lwrk,iwrk,kwrk,ier)
enddo
do i=n_tht/2+1,n_tht
  angle = 2.d0 * PI * float(i-1)/float(n_tht-1)
  r_bnd(i) = r0 + a0 * cos(angle + tria_l*sin(angle) + quad_l*sin(2.d0*angle))
  z_bnd(i) = z0 + a0 * ellip * sin(angle)
  call bispev(tx,nx,ty,ny,c,kx,ky,r_bnd(i),1,z_bnd(i),1,psi_bnd(i),wrk,lwrk,iwrk,kwrk,ier)
enddo


write(*,*) ' plotting results'  
nc = 51
allocate(zc(nc))
call begplt('eqdsk.ps')
call lblbot('eqdsk data',10)

call cplot(22,1,0,xx,yy,nr,nz,1,1,psirz,n_psi,zc,-nc,'fluxcontours',12,'R [m]',5,'Z [m]',5)
call lincol(3)
call lplot6(2,1,rlim,zlim,-limitr,'limiter')
call lincol(1)
call lplot6(2,1,rbnd,zbnd,-nbbs,'boundary')
call lincol(2)
call lplot6(2,1,r_bnd,z_bnd,-n_tht,'JOREK boundary')
call lincol(0)
call lplot6(3,2,psi,p,n_psi,'pressure')
call lplot6(3,3,psi,q,n_psi,'q')

call lplot6(2,2,psi,df2,n_psi,'df2')
call lplot6(3,2,psi,p,n_psi,'pressure')
call lplot6(2,3,psi,f,n_psi,'f')
call lplot6(3,3,psi,q,n_psi,'q')


!---------------------------- write JOREK input files
n_sol = (n_psi-1)/2
n_ext = n_psi + n_sol

write(*,*) ' n_psi, n_sol, n_ext : ',n_psi, n_sol, n_ext

allocate(df2_ext(n_ext),rho_ext(n_ext),T_ext(n_ext),psi_ext(n_ext),p_ext(n_ext))

df2_ext(1:n_psi) = df2(1:n_psi)
rho_ext(1:n_psi) = 1.d0
T_ext(1:n_psi)   = p(1:n_psi)

df2_ext(n_psi-1:n_ext) = df2_ext(n_psi)
rho_ext(n_psi-1:n_ext) = rho_ext(n_psi)
T_ext(n_psi-1:n_ext)   = T_ext(n_psi)

psi_sep = 1.0d0     ! in normalised psi units
sig_sep = 0.005     ! in normalised psi units
rho_bnd = 0.01      ! in jorek units
T_bnd   = 1.d-5     ! in jorek units

psi_ext(1:n_psi) = psi(1:n_psi)
do i=n_psi+1,n_ext
  psi_ext(i) = 1.d0 + 0.5 * float(i-n_psi)/float(n_sol)
enddo

zmu0 = 4.d-7 * PI

do i=1,n_ext
  tanh1 = tanh((psi_ext(i) - psi_sep)/sig_sep)
  df2_ext(i) = df2_ext(i) * (0.5d0 - 0.5d0*tanh1)
  rho_ext(i) = (rho_ext(i) - rho_bnd) * (0.5d0 - 0.5d0*tanh1) + rho_bnd
  T_ext(i)   = T_ext(i)   * (0.5d0 - 0.5d0*tanh1) * zmu0 +T_bnd 
!   T_ext(i)   = T_ext(i) / rho_ext(i) * zmu0 + T_bnd 
  p_ext(i)   = rho_ext(i) * T_ext(i)
enddo

call lplot6(2,2,psi_ext,df2_ext,n_ext,'df2')
call lplot6(3,2,psi_ext,p_ext,n_ext,'pressure')
call lplot6(2,3,psi_ext,rho_ext,n_ext,'density')
call lplot6(3,3,psi_ext,T_ext,n_ext,'T')

call lincol(1)
call lplot6(2,2,psi,df2,-n_psi,'df2')
call lplot6(3,2,psi,p*zmu0,-n_psi,'pressure')

open(21,file='jorek_ffprime')
! We change or not the sign of ff' depending on the sign of Ip because (we assume that) in EQDSK files, 
! psi_axis is always < psi_boundary, whatever the direction of Ip.
if (xip>0) then
  do i=1,n_ext
    write(21,*) psi_ext(i),-df2_ext(i) ! The minus sign is because ff' in JOREK is opposite to the usual ff' for historical reasons.
  enddo  
else
  do i=1,n_ext
    write(21,*) psi_ext(i),df2_ext(i) 
  enddo  
end if
close(21)

open(21,file='jorek_density')
do i=1,n_ext
  write(21,*) psi_ext(i),rho_ext(i)
enddo
close(21)
open(21,file='jorek_temperature')
do i=1,n_ext
  write(21,*) psi_ext(i),T_ext(i)
enddo
close(21)

open(21,file='jorek_namelist')

write(21,*)             '***************************************'
write(21,'(A)')        '*  namelist produced by eqdsk2jorek   *'
write(21,*)             '***************************************'
write(21,*) AA
write(21,'(A,f8.3,A)') '   magnetic field   : ',Bcentr,' [T]'
write(21,'(A,f8.3,A)') '   current          : ',xip/1d6,' [MA]'
write(21,'(A,e14.6,A)')'   central pressure : ',p(1), '[Pa]'
write(21,*)             '***************************************'
write(21,*)
write(21,*) ' &in1'
write(21,*) ' restart = .f.'
write(21,*) ' regrid  = .f.'
write(21,*) ' tstep   = 5.' 
write(21,*) ' nstep   = 0' 
write(21,*) ' nout = 1'
write(21,*)
write(21,*) ' time_evol_scheme= "Gears"'
write(21,*)  
write(21,*) ' tgnum(2) = 0.'
write(21,*) ' tgnum(5) = 0.'
write(21,*) ' tgnum(6) = 0.'
write(21,*) ' tgnum(7) = 0.'
write(21,*)
write(21,*) ' gmres_4     = 1.d0'
write(21,*) ' iter_precon = 21'
write(21,*) ' gmres_tol   = 1.d-8'
write(21,*) ' gmres_m     = 20'
write(21,*)
write(21,*) ' !use_mumps =.true.'
write(21,*) ' !use_pastix = .false.'
write(21,*) ' write_ps = .false.'
write(21,*)
write(21,*) ' !_____________________________________boundary definition'
write(21,*) ' mf = 0'
write(21,*) ' n_boundary = ',n_tht
write(21,*)
! We change or not the sign of psi_bnd depending on the sign of Ip because (we assume that) in EQDSK files, 
! psi_axis is always < psi_boundary, whatever the direction of Ip.
if (xip>0) then
  do j=1,n_tht
    write(21,'(A,i3,A,e16.8,A,i3,A,e16.8,A,i3,A,e16.8,A)'), &
             '  R_boundary(',j,') =',r_bnd(j), &
             ', Z_boundary(',j,') =',z_bnd(j), &
             ', psi_boundary(',j,') =',psi_bnd(j),','
  enddo
else
  do j=1,n_tht
    write(21,'(A,i3,A,e16.8,A,i3,A,e16.8,A,i3,A,e16.8,A)'), &
             '  R_boundary(',j,') =',r_bnd(j), &
             ', Z_boundary(',j,') =',z_bnd(j), &
             ', psi_boundary(',j,') =',-psi_bnd(j),','
  enddo	     
end if
write(21,*)
write(21,*) ' ellip  = ',ellip
write(21,*) ' tria_u = ',tria_u
write(21,*) ' tria_l = ',tria_l
write(21,*) ' quad_u = ',quad_u
write(21,*) ' quad_l = ',quad_l
write(21,*)
write(21,*) ' xampl  = +0.'
write(21,*) ' xpoint = .t.'
write(21,*)
write(21,*) ' freeboundary = .f.'
write(21,*) ' resistive_wall = .f.'
write(21,*)
write(21,*) ' R_geo = ',r0
write(21,*) ' Z_geo = ',z0
if (tokamak_name=='JET') then
  write(21,*) ' F0    = ',-2.96*bcentr ! By convention, the vacuum toroidal field is given at 2.96m in JET eqdsk files. 					
else
  write(21,*) ' F0    = ',-r0*bcentr
end if
write(21,*) ' amin  = 1.d0 ! scale factor for plasma size only'
write(21,*)
write(21,*) ' psi_axis_init  = ', ssimag 
write(21,*)
write(21,*) ' !_____________________________________grid parameters'

write(21,*) ' n_R      = 0'
write(21,*) ' n_Z      = 0'
write(21,*) ' n_radial = 101'
write(21,*) ' n_pol    = 128' 

write(21,*) ' n_flux   = 0'
write(21,*) ' n_tht    = 64'
write(21,*) ' n_open   = 15'
write(21,*) ' n_leg    = 15'
write(21,*) ' n_private = 9'
write(21,*) ' dPSI_open    = 0.04'
write(21,*) ' dPSI_private = 0.02'

write(21,*) ' amix = 0.d0'
write(21,*) ' !_____________________________________physics parameters'

write(21,*) ' eta   = 1.d-7'
write(21,*) ' visco = 1.d-6'
write(21,*) ' visco_par = 1.d-5'
write(21,*)
write(21,*) ' eta_num       = 1.d-12'
write(21,*) ' visco_num     = 1.d-12'
write(21,*) ' visco_par_num = 1.d-12'
write(21,*) ' d_perp_num    = 1.d-12'
write(21,*) ' zk_perp_num   = 1.d-12'
write(21,*)
write(21,*) ' bc_natural_open = .false.'
write(21,*) ' gamma_sheath = 2'
write(21,*)
write(21,*) ' rho_file     = "jorek_density"'
write(21,*) ' T_file       = "jorek_temperature"'
write(21,*) ' ffprime_file = "jorek_ffprime"'
write(21,*)
write(21,*) ' D_par     = 0.d0'
write(21,*) ' D_perp(1) = 1.d-5'
write(21,*) ' D_perp(2) = 0.85d0'
write(21,*) ' D_perp(3) = 0.d0'
write(21,*) ' D_perp(4) = 0.01d0'
write(21,*) ' D_perp(5) = 0.92d0'
write(21,*)
write(21,*) ' ZK_par     = 1.d0'
write(21,*) ' ZK_perp(1) = 1.d-5'
write(21,*) ' ZK_perp(2) = 0.85d0'
write(21,*) ' ZK_perp(3) = 0.d0'
write(21,*) ' ZK_perp(4) = 0.01d0'
write(21,*) ' ZK_perp(5) = 0.92d0'
write(21,*)
write(21,*) ' heatsource     = 1.d-7'
write(21,*) ' particlesource = 5.d-6'
write(21,*)
write(21,*) ' &end'

close(21)

call finplt

!-------------------------------------- eqdsk to vtk (careful VTK expects single precision)
lf   = char(10)
ivtk = 23

write(*,'(A)') ' writing VTK output'

open(unit=ivtk,file='eqdsk.vtk',form='unformatted',access='stream',convert='BIG_ENDIAN')

buffer = '# vtk DataFile Version 3.0'//lf                        ; write(ivtk) trim(buffer)
buffer = 'eqdsk'//lf                                             ; write(ivtk) trim(buffer)
buffer = 'BINARY'//lf                                            ; write(ivtk) trim(buffer)
buffer = 'DATASET RECTILINEAR_GRID'//lf                          ; write(ivtk) trim(buffer)

write(str2(1:24),'(3i8)') nr,nz,1
buffer = 'DIMENSIONS '//str2//lf                             ; write(ivtk) trim(buffer)

write(str1(1:12),'(i12)') nr
buffer = 'X_COORDINATES '//str1//' FLOAT'//lf                    ; write(ivtk) trim(buffer)
write(ivtk) (real(xx(i),4),i=1,nr)

write(str1(1:12),'(i12)') nz
buffer = 'Y_COORDINATES '//str1//' FLOAT'//lf                    ; write(ivtk) trim(buffer)
write(ivtk) (real(yy(j),4),j=1,nz)

write(str1(1:12),'(i12)') 1
buffer = 'Z_COORDINATES '//str1//' FLOAT'//lf                    ; write(ivtk) trim(buffer)
write(ivtk) real(0.d0,4)

! POINT_DATA SECTION
write(str1(1:12),'(i12)') nr*nz
buffer = lf//lf//'POINT_DATA '//str1//lf                             ; write(ivtk) trim(buffer)

buffer = 'SCALARS psi float'//lf                                     ; write(ivtk) trim(buffer)
buffer = 'LOOKUP_TABLE default'//lf                                  ; write(ivtk) trim(buffer)
write(ivtk) ((real(psirz(i,j),4), i=1,nr), j=1,nz)

close(ivtk)

end
 
 

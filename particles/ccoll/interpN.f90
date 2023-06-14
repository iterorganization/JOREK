MODULE interpN
  USE interpolatorss
  IMPLICIT NONE
  INTEGER, PARAMETER :: FIND_DICHO=1,FIND_LPARSE=2,FIND_RPARSE=3
  INTEGER, PARAMETER :: MAXINT=2147483647

CONTAINS

  REAL(kind=8) FUNCTION interp1(x,f,xi,nx)
    INTEGER,INTENT(in)::nx
    REAL*8,INTENT(in)::x(1:nx)
    REAL*8,INTENT(in)::f(1:nx)
    REAL*8,INTENT(in)::xi
    INTEGER :: i
    !     If out of bounds, just release a null value
    IF(xi.GT.x(nx).OR.xi.LT.x(1)) THEN
       interp1 = 0.0D0
    ELSE
       ! find the lowest cardinal point to xi
       i = lower_index(x,xi,nx)
       interp1 = linear_interp(norm_it(xi,x(i:i+1)),f(i:i+1))
    ENDIF
  END FUNCTION interp1

  !--------------------------------------------------

  REAL(kind=8) FUNCTION interp2(x,y,f,xi,yi,nx,ny)
    IMPLICIT NONE
    INTEGER,INTENT(in)::nx,ny
    REAL*8,INTENT(in)::x(1:nx)
    REAL*8,INTENT(in)::y(1:ny)
    REAL*8,INTENT(in)::f(1:nx,1:ny)
    REAL*8,INTENT(in)::xi
    REAL*8,INTENT(in)::yi

    INTEGER :: i,j

    !     If out of bounds, just release a null diffusivity
    IF(xi.GT.x(nx).OR.xi.LT.x(1)&
         .OR.&
         yi.GT.y(ny).OR.yi.LT.y(1)) THEN
       interp2 = 0.0D0
    ELSE
       !     Find the cardinal point of the first cell in x
       i = lower_index(x,xi,nx)
       j = lower_index(y,yi,ny)

       interp2 = bilinear_interp(norm_it(xi,x(i:i+1)),norm_it(yi,y(j:j+1)),f(i:i+1,j:j+1))
    ENDIF
  END FUNCTION interp2

  !--------------------------------------------------

  REAL(kind=8) FUNCTION interp3(x,y,z,f,xi,yi,zi,nx,ny,nz)
    IMPLICIT NONE
    INTEGER,INTENT(in)::nx,ny,nz

    REAL*8,INTENT(in)::x(1:nx)
    REAL*8,INTENT(in)::y(1:ny)
    REAL*8,INTENT(in)::z(1:nz)

    REAL*8,DIMENSION(nx,ny,nz),INTENT(in)::f

    REAL*8,INTENT(in)::xi
    REAL*8,INTENT(in)::yi
    REAL*8,INTENT(in)::zi

    INTEGER::i,j,k

    !     If out of bounds, just release a null diffusivity
    IF(xi.GT.x(nx).OR.xi.LT.x(1)&
         .OR.&
         yi.GT.y(ny).OR.yi.LT.y(1)&
         .OR.&
         zi.GT.z(nz).OR.zi.LT.z(1)) THEN
       interp3 = 0.0D0
    ELSE
       !     Find the cardinal point of the first cell in x
       i = lower_index(x,xi,nx)
       j = lower_index(y,yi,ny)
       k = lower_index(z,zi,nz)
       
       ! BUGGGGGGGGGGGGGGGGGGGGGG, it was like this
!        interp3 = trilinear_interp(&
!             norm_it(xi,x(i:i+1)),&
!             norm_it(yi,y(i:i+1)),&
!             norm_it(zi,z(i:i+1)),&
!             f(i:i+1,j:j+1,k:k+1))

       ! Should be like thisssssssss
       interp3 = trilinear_interp(&
            norm_it(xi,x(i:i+1)),&
            norm_it(yi,y(j:j+1)),&
            norm_it(zi,z(k:k+1)),&
            f(i:i+1,j:j+1,k:k+1))
    END IF
  END FUNCTION interp3

  !--------------------------------------------------

  REAL(kind=8) FUNCTION interp4(x,y,z,v,f,xi,yi,zi,vi,nx,ny,nz,nv)
    IMPLICIT NONE
    INTEGER,INTENT(in)::nx,ny,nz,nv

    REAL*8,INTENT(in)::x(1:nx)
    REAL*8,INTENT(in)::y(1:ny)
    REAL*8,INTENT(in)::z(1:nz)
    REAL*8,INTENT(in)::v(1:nv)

    REAL*8,INTENT(in)::f(1:nx,1:ny,1:nz,1:nv)

    REAL*8,INTENT(in)::xi
    REAL*8,INTENT(in)::yi
    REAL*8,INTENT(in)::zi
    REAL*8,INTENT(in)::vi

    INTEGER::i,j,k,w

    !     If out of bounds, just release a null diffusivity
    IF(xi.GT.x(nx).OR.xi.LT.x(1)&
         .OR.&
         yi.GT.y(ny).OR.yi.LT.y(1)&
         .OR.&
         zi.GT.z(nz).OR.zi.LT.z(1)&
         .OR.&
         vi.GT.v(nv).OR.vi.LT.v(1)) THEN
       interp4 = 0.0D0
    ELSE
       !     Find the cardinal point of the first cell in x
       i = lower_index(x,xi,nx)
       j = lower_index(y,yi,ny)
       k = lower_index(z,zi,nz)
       w = lower_index(v,vi,nv)
  
       interp4=quadrilinear_interp(&
            norm_it(xi,x(i:i+1)),&
            norm_it(yi,y(j:j+1)),&
            norm_it(zi,z(k:k+1)),&
            norm_it(vi,v(w:w+1)),&
            f(i:i+1,j:j+1,k:k+1,w:w+1))
    END IF
  END FUNCTION interp4

  !**************************************************
  ! Here we propose a few algorithms to find 
  ! where xi is located inside x array
  !**************************************************
  INTEGER FUNCTION lower_index(x,xi,nx,ALG)
    IMPLICIT NONE

    INTEGER,INTENT(in) :: nx
    REAL*8, INTENT(in) :: x(1:nx)
    REAL*8, INTENT(in) :: xi
    INTEGER, OPTIONAL, INTENT(IN) :: ALG
    INTEGER::pivot,upper,choice

    choice=FIND_DICHO
    IF(PRESENT(ALG))choice=ALG

    ! Warning : Commented code makes sure that xi belongs 
    !           to the given interval. If you are sure this it is always 
    !           the case, you may keep it commented (for optimisation).
    !
    ! IF(xi.GT.x(nx).OR.xi.LT.x(1)) THEN
    !    print*,'Lower index cannot be found. Interval does not contain xi.'
    !    lower_index=MAXINT+1
    ! ELSE

    ! Dichotomy algorithm : 
    ! Take half-segments until pin-pointed the correct one
    lower_index=1
    upper=nx
    pivot=(nx+1)/2
    DO WHILE((upper-lower_index).GT.1)
       IF (xi.GE.x(pivot)) THEN
          lower_index=pivot
          !if(xi.LE.x(pivot+1)) upper=pivot+1           ! using this peak may be faster (try commenting)
       ELSE
          upper=pivot
          !if(xi.GE.x(pivot-1)) lower_index=pivot-1     ! using this peak may be faster  (try commenting)
       ENDIF
       pivot=(upper+lower_index)/2
    ENDDO
    !ENDIF
  END FUNCTION lower_index

  INTEGER FUNCTION closest_index(x,xi,nx)
    IMPLICIT NONE

    INTEGER,INTENT(in) :: nx
    REAL*8, INTENT(in) :: x(1:nx)
    REAL*8, INTENT(in) :: xi
    INTEGER :: i
    
    IF(nx.EQ.1) THEN
       closest_index = 1
    ELSE
       i = lower_index(x,xi,nx)
       IF(xi-x(i).LE.x(i+1)-xi) THEN
          closest_index = i
       ELSE
          closest_index = i+1
       ENDIF
    ENDIF
  END FUNCTION closest_index

END MODULE interpN

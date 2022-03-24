!> Routines for reading/writing/manipulating 1D data profiles.
!!
!! The module is used to handle numerical input profiles for temperature,
!! density, FFprime, D_perp, and ZK_perp.
module profiles
  
  use tr_module 
  
  implicit none
  
  
  
  private
  
  ! --- Public routines
  public constructProf, destructProf, resizeProf, readProf, writProf, derivProf, interpolProf, &
       constructProfNeo, destructProfNeo, resizeProfNeo, readProfNeo
  
  
  
  contains
  
  
  
  !> Interpolate a given profile to a certain position x0.
  real*8 function interpolProf(x, y, len, x0)
    real*8, allocatable, intent(in) :: x(:)
    real*8, allocatable, intent(in) :: y(:)
    integer,             intent(in) :: len
    real*8,              intent(in) :: x0
    
    integer :: left, mid, right
    real*8  :: aux1, aux2

#if _OPENMP >= 201511
    !$omp declare simd uniform(x,y,len)
#endif
    
    left  = 1
    right = len
    do
      if ( right == left + 1 ) exit
      mid = (left + right) / 2
      if ( x(mid) >= x0 ) then
        right = mid
      else
        left = mid
      end if
    end do
    aux1  = (x0 - x(left)) / (x(right) - x(left))
    aux2  = (1. - aux1)
    interpolProf = y(left) * aux2 + y(right) * aux1
    
  end function interpolProf
  
  
  
  !> Construct a profile.
  recursive subroutine constructProf(x, y, len)

    real*8, allocatable, intent(inout) :: x(:), y(:)
    integer,           intent(inout) :: len
    
    call destructProf(x, y, len)
    
    len = max(len,1) ! at least length 1
    call tr_allocate(x,1,len,"x",CAT_GRID)
    call tr_allocate(y,1,len,"y",CAT_GRID)
    
  end subroutine constructProf
  
  
  
  !> Change the size of a profile.
  recursive subroutine resizeProf(x, y, len, newLen, keep)

    real*8, allocatable, intent(inout) :: x(:), y(:)
    integer,           intent(inout) :: len
    integer,           intent(in)    :: newLen
    logical,           intent(in)    :: keep !< Keep the data in the x and y arrays?
    
    real*8, ALLOCATABLE            :: px(:) ! copy of x (in case keep=.TRUE.)
    real*8, ALLOCATABLE            :: py(:) ! copy of y (in case keep=.TRUE.)
    
    ! --- Recursive call with newLen=1 if newLen < 1.
    if ( newLen < 1 ) then
      call resizeProf(x, y, len, 1, keep)
      return
    end if
    
    ! --- Backup data from profile if keep=.true.
    if ( keep ) then
      call tr_allocate(px,1,len,"px",CAT_GRID)
      call tr_allocate(py,1,len,"py",CAT_GRID)
      if ( allocated(x) ) then
        px(1:len) = x(1:len)
      else
        px = 0.
      end if
      if ( allocated(y) ) then
        py(1:len) = y(1:len)
      else
        py = 0.
      end if
    end if
    
    ! --- Resize x and y.
    if ( allocated(x) ) call tr_deallocate(x,"x",CAT_GRID)
    if ( allocated(y) ) call tr_deallocate(y,"y",CAT_GRID)
    call tr_allocate(x,1,newLen,"x",CAT_GRID)
    call tr_allocate(y,1,newLen,"y",CAT_GRID)
    
    ! --- Restore data to profile if keep=.true.
    if ( keep ) then
      x(1:min(len,newLen)) = px(1:min(len,newLen))
      y(1:min(len,newLen)) = py(1:min(len,newLen))
      call tr_deallocate(px,"px",CAT_GRID)
      call tr_deallocate(py,"py",CAT_GRID)
    end if
    len = newLen
    
  end subroutine resizeProf
  
  
  
  !> Destroy a profile.
  recursive subroutine destructProf(x, y, len)

    real*8, allocatable, intent(inout) :: x(:), y(:)
    integer,           intent(inout) :: len
    
    if ( allocated(x) ) call tr_deallocate(x,"x",CAT_GRID)
    if ( allocated(y) ) call tr_deallocate(y,"y",CAT_GRID)
    len = 0
    
  end subroutine destructprof
  
  
  
  !> Read a profile from a file.
  recursive subroutine readProf(x, y, len, file)
    
    real*8, allocatable, intent(inout) :: x(:), y(:)
    integer,           intent(inout)   :: len
    CHARACTER(LEN=*),  intent(in)      :: file    !< Filename.
    
    integer :: err
    integer :: usedLen
    real*8  :: xx, yy
    
    call destructProf(x, y, len)
    
    ! --- Open the file.
    OPEN(UNIT=42, FILE=file, FORM='FORMATTED', STATUS='OLD', ACTION='READ', IOSTAT=err)
    if ( err /= 0 ) then
      write(*,*) 'ERROR in readProf: Cannot open file '//TRIM(file)//'.'
      return
    end if
    
    ! --- Construct prof with an initial length of 10.
    len = 10
    call constructProf(x, y, len)
    usedLen = 0
    
    ! --- Read profile.
    do
      READ(42, *, IOSTAT=err) xx, yy
      
      if ( err /= 0 ) exit ! end of profile reached
      
      usedLen = usedLen + 1
      
      ! Double the profile length if it becomes to small.
      if ( usedLen > len ) call resizeProf(x, y, len, 2*len, .TRUE.)
      
      x(usedLen) = xx
      y(usedLen) = yy
      
    end do
    
    ! --- Crop the profile to the length that is really used.
    call resizeProf(x, y, len, usedLen, .TRUE.)
    
    ! --- Close the file.
    CLOSE(UNIT=42)

  end subroutine readProf
  
  
  
  !> Write a profile to a file.
  recursive subroutine writProf(x, y, len, file)
    
    real*8, allocatable, intent(in)    :: x(:), y(:)
    integer,           intent(in)    :: len
    CHARACTER(LEN=*),  intent(in)    :: file    !< Filename.
    
    integer :: err
    integer :: i
    
    if ( (.not. allocated(x)) .or. (.not. allocated(y)) ) return
    
    ! --- Open the output file.
    OPEN(UNIT=42, FILE=file, FORM='FORMATTED', STATUS='REPLACE', ACTION='write', IOSTAT=err)
    if ( err /= 0 ) return
    
    ! --- Write data.
    do i = 1, len
      write(42, '(ES22.15,1X,ES22.15)') x(i), y(i)
    end do
    
    ! --- Close the output file.
    CLOSE(UNIT=42)

  end subroutine writProf
  
  
  
  
  !> Determine the derivative of a profile.
  recursive subroutine derivProf(x, y, len, yd)
  
    real*8, allocatable, intent(in)    :: x(:), y(:)
    real*8, allocatable, intent(inout) :: yd(:)
    integer,           intent(inout) :: len
    
    integer :: i       ! Point index of profile.
    real*8  :: d(-2:4) ! Distances between points.
    real*8  :: f(-2:4) ! Function values.
    real*8  :: c(-2:4) ! Coefficients for function values.
    
    ! The derivatives will be determined at the same x-positions as the profile.
    if ( allocated(yd) ) call tr_deallocate(yd,"yd",CAT_GRID)
    call tr_allocate(yd,1,len,"yd",CAT_GRID) 
    
    do i = 1, len
      
      c = 0.
      f = 0.
      f(0) = y(i)
      
      if ( (i==1) .or. (i==len) ) then
        
        if ( i==1 ) then
          f(0) = y(i)
          f(1) = y(2)
          f(2) = y(3)
          
          d(1) = x(2) - x(i)
          d(2) = x(3) - x(i)
        else
          f(0) = y(i)
          f(1) = y(len-1)
          f(2) = y(len-2)
          
          d(1) = x(len-1) - x(i)
          d(2) = x(len-2) - x(i)
        end if
        c(0) = -(1/d(1)) - 1/d(2)
        c(1) = -(d(2)/(d(1)**2 - d(1)*d(2)))
        c(2) = d(1)/((d(1) - d(2))*d(2))
        
      else
        
        f(-1) = y(i-1)
        f(0)  = y(i)
        f(1)  = y(i+1)
        
        d(-1) = x(i-1) - x(i)
        d(1)  = x(i+1) - x(i)
        
        c(-1) = -(d(1)/(d(-1)**2 - d(-1)*d(1)))
        c(0)  = -(1/d(-1)) - 1/d(1)
        c(1)  = d(-1)/((d(-1) - d(1))*d(1))

      end if
      
      yd(i) = sum( c * f )
        
    end do
    
  end subroutine derivProf

  ! Construct a profile.
  recursive subroutine constructProfNeo(x1, x2, x3, len)

    real*8, allocatable, intent(inout) :: x1(:), x2(:), x3(:)
    integer,             intent(inout) :: len
    
    call destructProfNeo(x1, x2, x3, len)
    
    len = max(len,1) ! at least length 1
    call tr_allocate(x1,1,len,"x1",CAT_GRID)
    call tr_allocate(x2,1,len,"x2",CAT_GRID)
    call tr_allocate(x3,1,len,"x3",CAT_GRID)
    
  end subroutine constructProfNeo

 
  ! Change the size of a profile.
  recursive subroutine resizeProfNeo(x1, x2, x3, len, newLen, keep)

    real*8, allocatable, intent(inout) :: x1(:), x2(:), x3(:)
    integer,           intent(inout)   :: len
    integer,           intent(in)      :: newLen
    logical,           intent(in)      :: keep !< Keep the data in the x1, x2, x3 arrays?
    
    real*8, ALLOCATABLE          :: px1(:) ! copy of x1 (in case keep=.TRUE.)
    real*8, ALLOCATABLE          :: px2(:) ! copy of x2 (in case keep=.TRUE.)
    real*8, ALLOCATABLE          :: px3(:) ! copy of x3 (in case keep=.TRUE.)
    
    ! --- Recursive call with newLen=1 if newLen < 1.
    if ( newLen < 1 ) then
      call resizeProfNeo(x1, x2, x3, len, 1, keep)
      return
    end if
    
    ! --- Backup data from profile if keep=.true.
    if ( keep ) then
      call tr_allocate(px1,1,len,"px1",CAT_GRID)
      call tr_allocate(px2,1,len,"px2",CAT_GRID)
      call tr_allocate(px3,1,len,"px3",CAT_GRID)
      if ( allocated(x1) ) then
        px1(1:len) = x1(1:len)
      else
        px1 = 0.
      end if
      if ( allocated(x2) ) then
        px2(1:len) = x2(1:len)
      else
        px2 = 0.
      end if
      if ( allocated(x3) ) then
        px3(1:len) = x3(1:len)
      else
        px3 = 0.
      end if
    end if
    
    ! --- Resize x1, x2 and x3.
    if ( allocated(x1) ) call tr_deallocate(x1,"x1",CAT_GRID)
    if ( allocated(x2) ) call tr_deallocate(x2,"x2",CAT_GRID)
    if ( allocated(x3) ) call tr_deallocate(x3,"x3",CAT_GRID)
    call tr_allocate(x1,1,newLen,"x1",CAT_GRID)
    call tr_allocate(x2,1,newLen,"x2",CAT_GRID)
    call tr_allocate(x3,1,newLen,"x3",CAT_GRID)
    
    ! --- Restore data to profile if keep=.true.
    if ( keep ) then
      x1(1:min(len,newLen)) = px1(1:min(len,newLen))
      x2(1:min(len,newLen)) = px2(1:min(len,newLen))
      x3(1:min(len,newLen)) = px3(1:min(len,newLen))
      call tr_deallocate(px1,"px1",CAT_GRID)
      call tr_deallocate(px2,"px2",CAT_GRID)
      call tr_deallocate(px3,"px3",CAT_GRID)
    end if
    len = newLen
    
  end subroutine resizeProfNeo

  
  ! Destroy a profile.
  recursive subroutine destructProfNeo(x1, x2, x3, len)

    real*8, allocatable, intent(inout) :: x1(:), x2(:), x3(:)
    integer,           intent(inout) :: len
    
    if ( allocated(x1) ) call tr_deallocate(x1,"x1",CAT_GRID)
    if ( allocated(x2) ) call tr_deallocate(x2,"x2",CAT_GRID)
    if ( allocated(x3) ) call tr_deallocate(x3,"x3",CAT_GRID)
    len = 0
    
  end subroutine destructProfNeo
  
  
  ! Read a profile from a file.
  recursive subroutine readProfNeo(x1, x2, x3, len, file)
    
    real*8, allocatable, intent(inout) :: x1(:), x2(:), x3(:)
    integer,           intent(inout)   :: len
    CHARACTER(LEN=*),  intent(in)      :: file    ! Filename.
    
    integer :: err
    integer :: usedLen
    real*8  :: xx1, xx2, xx3

   
    call destructProfNeo(x1, x2, x3, len)
    
    ! --- Open the file.
    OPEN(UNIT=42, FILE=file, FORM='FORMATTED', STATUS='OLD', ACTION='READ', IOSTAT=err)
    if ( err /= 0 ) then
      write(*,*) 'ERROR in readProfNeo: Cannot open file '//TRIM(file)//'.'
      return
    end if

    
    ! --- Construct prof with an initial length of 10.
    len = 10
    call constructProfNeo(x1, x2, x3, len)
    usedLen = 0
    
    ! --- Read profile.
    do
      READ(42, *, IOSTAT=err) xx1, xx2, xx3

!      write (*,*) '********* read neo profiles **********'
!      write (*,*) xx1, xx2, xx3
      
      if ( err /= 0 ) exit ! end of profile reached
      
      usedLen = usedLen + 1
      
      ! Double the profile length if it becomes to small.
      if ( usedLen > len ) call resizeProfNeo(x1, x2, x3, len, 2*len, .TRUE.)
      
      x1(usedLen) = xx1
      x2(usedLen) = xx2
      x3(usedLen) = xx3
      
    end do

    
    ! --- Crop the profile to the length that is really used.
    call resizeProfNeo(x1, x2, x3, len, usedLen, .TRUE.)

    
    ! --- Close the file.
    CLOSE(UNIT=42)

  end subroutine readProfNeo
  
end module profiles

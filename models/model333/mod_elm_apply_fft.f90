module mod_elm_apply_fft
implicit none
contains

!------------------------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------- Compute the FFT for the matrix elements -----------------------------------------------
!------------------------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------------------------
subroutine ELM_apply_fft(RHS, RHS_p, RHS_k, ELM, ELM_p, ELM_n, ELM_k, ELM_kn, tid)
!DEC$ ATTRIBUTES FORCEINLINE :: ELM_apply_fft

  ! --- Modules
  use phys_module
  use mod_parameters    
  use data_structure
  
  implicit none
  
  ! --- Matrix elements and toroidal functions
  integer, intent(in)	      :: tid
  real*8, dimension (:,:)     :: ELM
  real*8, dimension (:)       :: RHS
  real*8, dimension(:,:,:)    :: ELM_p, ELM_n, ELM_k, ELM_kn
  real*8, dimension(:,:)      :: RHS_p, RHS_k 
  real*8		      :: in_fft(1:n_plane)
  complex*16		      :: out_fft(1:n_plane)
  
  ! --- FFT Indexes
  integer    :: i, j, k, l, m, ik, im
  integer    :: index, index_k, index_m
	  
  ! --- RHS_p
  do j=1, n_vertex_max*n_var*n_degrees

    in_fft = RHS_p(1:n_plane,j)

#ifdef USE_FFTW
    call dfftw_execute_dft_r2c(fftw_plan, in_fft, out_fft)
#else
    call my_fft(in_fft, out_fft, n_plane)
#endif

    index = n_tor*(j-1) + 1

    RHS(index) = real(out_fft(1))

    do k=2,(n_tor+1)/2

      index = n_tor*(j-1) + 2*(k-1)

      RHS(index)   =   real(out_fft(k))
      RHS(index+1) = - imag(out_fft(k))

    enddo

  enddo

  ! --- RHS_k
  do j=1, n_vertex_max*n_var*n_degrees

    in_fft = RHS_k(1:n_plane,j)

#ifdef USE_FFTW
    call dfftw_execute_dft_r2c(fftw_plan, in_fft, out_fft)
#else
    call my_fft(in_fft, out_fft, n_plane)
#endif

    index = n_tor*(j-1) + 1
    ik    = 1
  
    RHS(index) = RHS(index) + imag(out_fft(1)) * float(mode(ik))
  
    do k=2,(n_tor+1)/2
  
      ik    = max(2*(k-1),1)
      index = n_tor*(j-1) + 2*(k-1)
  
      RHS(index)   = RHS(index)   + imag(out_fft(k)) * float(mode(ik))
      RHS(index+1) = RHS(index+1) + real(out_fft(k)) * float(mode(ik))

    enddo

  enddo

  ! --- ELM_p
  do i=1,n_vertex_max*n_var*n_degrees
    do j=1, n_vertex_max*n_var*n_degrees

      in_fft =  ELM_p(1:n_plane,i,j)

#ifdef USE_FFTW
      call dfftw_execute_dft_r2c(fftw_plan, in_fft, out_fft)
#else
      call my_fft(in_fft, out_fft, n_plane)
#endif

      do k=1,(n_tor+1)/2

	index_k = n_tor*(i-1) + max(2*(k-1),1)

	do m=1,(n_tor+1)/2

	  index_m = n_tor*(j-1) + max(2*(m-1),1)

	  l = (k-1) + (m-1)

	  if ( (l .ge. 0) .and. (l .le. n_plane/2) ) then

	    ELM(index_k,  index_m  ) = ELM(index_k,  index_m)	+ real(out_fft(l+1))
	    ELM(index_k+1,index_m  ) = ELM(index_k+1,index_m)	- imag(out_fft(l+1))
	    ELM(index_k,  index_m+1) = ELM(index_k,  index_m+1) - imag(out_fft(l+1))
	    ELM(index_k+1,index_m+1) = ELM(index_k+1,index_m+1) - real(out_fft(l+1))

	  elseif ( (l .lt. 0) .and. (abs(l) .le. n_plane/2) ) then

	    ELM(index_k,  index_m  ) = ELM(index_k,  index_m)	+ real(out_fft(abs(l)+1))
	    ELM(index_k+1,index_m  ) = ELM(index_k+1,index_m)	+ imag(out_fft(abs(l)+1))
	    ELM(index_k,  index_m+1) = ELM(index_k,  index_m+1) + imag(out_fft(abs(l)+1))
	    ELM(index_k+1,index_m+1) = ELM(index_k+1,index_m+1) - real(out_fft(abs(l)+1))

	  endif

	  l = (k-1) - (m-1)

	  if ( (l .ge. 0) .and. (l .le. n_plane/2) ) then

	    ELM(index_k,  index_m  ) = ELM(index_k,  index_m)	+ real(out_fft(l+1))
	    ELM(index_k+1,index_m  ) = ELM(index_k+1,index_m)	- imag(out_fft(l+1))
	    ELM(index_k,  index_m+1) = ELM(index_k,  index_m+1) + imag(out_fft(l+1))
	    ELM(index_k+1,index_m+1) = ELM(index_k+1,index_m+1) + real(out_fft(l+1))

	  elseif ( (l .lt. 0) .and. (abs(l) .le. n_plane/2) ) then

	    ELM(index_k,  index_m  ) = ELM(index_k,  index_m)	+ real(out_fft(abs(l)+1))
	    ELM(index_k+1,index_m  ) = ELM(index_k+1,index_m)	+ imag(out_fft(abs(l)+1))
	    ELM(index_k,  index_m+1) = ELM(index_k,  index_m+1) - imag(out_fft(abs(l)+1))
	    ELM(index_k+1,index_m+1) = ELM(index_k+1,index_m+1) + real(out_fft(abs(l)+1))

	  endif

	enddo
      enddo
    enddo
  enddo

  ! --- ELM_n
  do i=1,n_vertex_max*n_var*n_degrees
    do j=1, n_vertex_max*n_var*n_degrees
      if (maxval(abs(ELM_n(1:n_plane,i,j))) .ne. 0.d0) then

	in_fft =  ELM_n(1:n_plane,i,j)

#ifdef USE_FFTW
        call dfftw_execute_dft_r2c(fftw_plan, in_fft, out_fft)
#else
	call my_fft(in_fft, out_fft, n_plane)
#endif

	do k=1,(n_tor+1)/2

	  index_k = n_tor*(i-1) + max(2*(k-1),1)

	  do m=1,(n_tor+1)/2

	    im = max(2*(m-1),1)
	    index_m = n_tor*(j-1) + max(2*(m-1),1)

	    l = (k-1) + (m-1)

	    if ( (l .ge. 0) .and. (l .le. n_plane/2) ) then

	      ELM(index_k,  index_m  ) = ELM(index_k,  index_m)   + imag(out_fft(l+1)) * float(mode(im))
	      ELM(index_k+1,index_m  ) = ELM(index_k+1,index_m)   + real(out_fft(l+1)) * float(mode(im))
	      ELM(index_k,  index_m+1) = ELM(index_k,  index_m+1) + real(out_fft(l+1)) * float(mode(im))
	      ELM(index_k+1,index_m+1) = ELM(index_k+1,index_m+1) - imag(out_fft(l+1)) * float(mode(im))

	    elseif ( (l .lt. 0) .and. (abs(l) .le. n_plane/2) ) then

	      ELM(index_k,  index_m  ) = ELM(index_k,  index_m)   - imag(out_fft(abs(l)+1)) * float(mode(im))
	      ELM(index_k+1,index_m  ) = ELM(index_k+1,index_m)   + real(out_fft(abs(l)+1)) * float(mode(im))
	      ELM(index_k,  index_m+1) = ELM(index_k,  index_m+1) + real(out_fft(abs(l)+1)) * float(mode(im))
	      ELM(index_k+1,index_m+1) = ELM(index_k+1,index_m+1) + imag(out_fft(abs(l)+1)) * float(mode(im))

	    endif

	    l = (k-1) - (m-1)

	    if ( (l .ge. 0) .and. (l .le. n_plane/2) ) then

	      ELM(index_k,  index_m  ) = ELM(index_k,  index_m)   - imag(out_fft(l+1)) * float(mode(im))
	      ELM(index_k+1,index_m  ) = ELM(index_k+1,index_m)   - real(out_fft(l+1)) * float(mode(im))
	      ELM(index_k,  index_m+1) = ELM(index_k,  index_m+1) + real(out_fft(l+1)) * float(mode(im))
	      ELM(index_k+1,index_m+1) = ELM(index_k+1,index_m+1) - imag(out_fft(l+1)) * float(mode(im))

	    elseif ( (l .lt. 0) .and. (abs(l) .le. n_plane/2) ) then

	      ELM(index_k,  index_m  ) = ELM(index_k,  index_m)   + imag(out_fft(abs(l)+1)) * float(mode(im))
	      ELM(index_k+1,index_m  ) = ELM(index_k+1,index_m)   - real(out_fft(abs(l)+1)) * float(mode(im))
	      ELM(index_k,  index_m+1) = ELM(index_k,  index_m+1) + real(out_fft(abs(l)+1)) * float(mode(im))
	      ELM(index_k+1,index_m+1) = ELM(index_k+1,index_m+1) + imag(out_fft(abs(l)+1)) * float(mode(im))

	    endif

	  enddo
	enddo
      endif
    enddo
  enddo

  ! --- ELM_k
  do i=1,n_vertex_max*n_var*n_degrees
    do j=1, n_vertex_max*n_var*n_degrees
      if (maxval(abs(ELM_k(1:n_plane,i,j))) .ne. 0.d0) then

	in_fft =  ELM_k(1:n_plane,i,j)

#ifdef USE_FFTW
        call dfftw_execute_dft_r2c(fftw_plan, in_fft, out_fft)
#else
	call my_fft(in_fft, out_fft, n_plane)
#endif

	do k=1,(n_tor+1)/2

	  ik	  = max(2*(k-1),1)
	  index_k = n_tor*(i-1) + max(2*(k-1),1)

	  do m=1,(n_tor+1)/2

	    index_m = n_tor*(j-1) + max(2*(m-1),1)

	    l = (k-1) + (m-1)

	    if ( (l .ge. 0) .and. (l .le. n_plane/2) ) then

	      ELM(index_k,  index_m  ) = ELM(index_k,  index_m)   + imag(out_fft(l+1)) * float(mode(ik))
	      ELM(index_k+1,index_m  ) = ELM(index_k+1,index_m)   + real(out_fft(l+1)) * float(mode(ik))
	      ELM(index_k,  index_m+1) = ELM(index_k,  index_m+1) + real(out_fft(l+1)) * float(mode(ik))
	      ELM(index_k+1,index_m+1) = ELM(index_k+1,index_m+1) - imag(out_fft(l+1)) * float(mode(ik))

	    elseif ( (l .lt. 0) .and. (abs(l) .le. n_plane/2) ) then

	      ELM(index_k,  index_m  ) = ELM(index_k,  index_m)   - imag(out_fft(abs(l)+1)) * float(mode(ik))
	      ELM(index_k+1,index_m  ) = ELM(index_k+1,index_m)   + real(out_fft(abs(l)+1)) * float(mode(ik))
	      ELM(index_k,  index_m+1) = ELM(index_k,  index_m+1) + real(out_fft(abs(l)+1)) * float(mode(ik))
	      ELM(index_k+1,index_m+1) = ELM(index_k+1,index_m+1) + imag(out_fft(abs(l)+1)) * float(mode(ik))

	    endif

	    l = (k-1) - (m-1)

	    if ( (l .ge. 0) .and. (l .le. n_plane/2) ) then

	      ELM(index_k,  index_m  ) = ELM(index_k,  index_m)   + imag(out_fft(l+1)) * float(mode(ik))
	      ELM(index_k+1,index_m  ) = ELM(index_k+1,index_m)   + real(out_fft(l+1)) * float(mode(ik))
	      ELM(index_k,  index_m+1) = ELM(index_k,  index_m+1) - real(out_fft(l+1)) * float(mode(ik))
	      ELM(index_k+1,index_m+1) = ELM(index_k+1,index_m+1) + imag(out_fft(l+1)) * float(mode(ik))

	    elseif ( (l .lt. 0) .and. (abs(l) .le. n_plane/2) ) then

	      ELM(index_k,  index_m  ) = ELM(index_k,  index_m)   - imag(out_fft(abs(l)+1)) * float(mode(ik))
	      ELM(index_k+1,index_m  ) = ELM(index_k+1,index_m)   + real(out_fft(abs(l)+1)) * float(mode(ik))
	      ELM(index_k,  index_m+1) = ELM(index_k,  index_m+1) - real(out_fft(abs(l)+1)) * float(mode(ik))
	      ELM(index_k+1,index_m+1) = ELM(index_k+1,index_m+1) - imag(out_fft(abs(l)+1)) * float(mode(ik))

	    endif

	  enddo
	enddo
      endif
    enddo
  enddo


  ! --- ELM_kn
  do i=1,n_vertex_max*n_var*n_degrees
    do j=1, n_vertex_max*n_var*n_degrees
      if (maxval(abs(ELM_kn(1:n_plane,i,j))) .ne. 0.d0) then

	in_fft =  ELM_kn(1:n_plane,i,j)

#ifdef USE_FFTW
        call dfftw_execute_dft_r2c(fftw_plan, in_fft, out_fft)
#else
	call my_fft(in_fft, out_fft, n_plane)
#endif

	do k=1,(n_tor+1)/2

	  ik	  = max(2*(k-1),1)
	  index_k = n_tor*(i-1) + max(2*(k-1),1)

	  do m=1,(n_tor+1)/2

	    im      = max(2*(m-1),1)
	    index_m = n_tor*(j-1) + max(2*(m-1),1)

	    l = (k-1) + (m-1)

	    if ( (l .ge. 0) .and. (l .le. n_plane/2) ) then

	       ELM(index_k,  index_m  ) = ELM(index_k,  index_m)   - real(out_fft(l+1)) * float(mode(im)) * float(mode(ik))
	       ELM(index_k+1,index_m  ) = ELM(index_k+1,index_m)   + imag(out_fft(l+1)) * float(mode(im)) * float(mode(ik))
	       ELM(index_k,  index_m+1) = ELM(index_k,  index_m+1) + imag(out_fft(l+1)) * float(mode(im)) * float(mode(ik))
	       ELM(index_k+1,index_m+1) = ELM(index_k+1,index_m+1) + real(out_fft(l+1)) * float(mode(im)) * float(mode(ik))

	    elseif ( (l .lt. 0) .and. (abs(l) .le. n_plane/2) ) then

	       ELM(index_k,  index_m  ) = ELM(index_k,  index_m)   - real(out_fft(abs(l)+1)) * float(mode(im)) * float(mode(ik))
	       ELM(index_k+1,index_m  ) = ELM(index_k+1,index_m)   - imag(out_fft(abs(l)+1)) * float(mode(im)) * float(mode(ik))
	       ELM(index_k,  index_m+1) = ELM(index_k,  index_m+1) - imag(out_fft(abs(l)+1)) * float(mode(im)) * float(mode(ik))
	       ELM(index_k+1,index_m+1) = ELM(index_k+1,index_m+1) + real(out_fft(abs(l)+1)) * float(mode(im)) * float(mode(ik))

	    endif

	    l = (k-1) - (m-1)

	    if ( (l .ge. 0) .and. (l .le. n_plane/2) ) then

	      ELM(index_k,  index_m  ) = ELM(index_k,  index_m)   + real(out_fft(l+1)) * float(mode(im)) * float(mode(ik))
	      ELM(index_k+1,index_m  ) = ELM(index_k+1,index_m)   - imag(out_fft(l+1)) * float(mode(im)) * float(mode(ik))
	      ELM(index_k,  index_m+1) = ELM(index_k,  index_m+1) + imag(out_fft(l+1)) * float(mode(im)) * float(mode(ik))
	      ELM(index_k+1,index_m+1) = ELM(index_k+1,index_m+1) + real(out_fft(l+1)) * float(mode(im)) * float(mode(ik))

	    elseif ( (l .lt. 0) .and. (abs(l) .le. n_plane/2) ) then

	      ELM(index_k,  index_m  ) = ELM(index_k,  index_m)   + real(out_fft(abs(l)+1)) * float(mode(im)) * float(mode(ik))
	      ELM(index_k+1,index_m  ) = ELM(index_k+1,index_m)   + imag(out_fft(abs(l)+1)) * float(mode(im)) * float(mode(ik))
	      ELM(index_k,  index_m+1) = ELM(index_k,  index_m+1) - imag(out_fft(abs(l)+1)) * float(mode(im)) * float(mode(ik))
	      ELM(index_k+1,index_m+1) = ELM(index_k+1,index_m+1) + real(out_fft(abs(l)+1)) * float(mode(im)) * float(mode(ik))

	    endif

	  enddo
	enddo
      endif
    enddo
  enddo

  ELM = 0.5d0 * ELM

  return
end subroutine ELM_apply_fft
















!------------------------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------------------------
!------------------------------------------------ The FFT routine -------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------------------------
subroutine my_fft(in_fft,out_fft,n)
!DEC$ ATTRIBUTES FORCEINLINE :: my_fft
  
  implicit none
    
  real*8     :: in_fft(*)
  complex*16 :: out_fft(*)	
  real*8     :: tmp_fft(2*n+2)
  integer    :: i, n
      
  tmp_fft(1:n) = in_fft(1:n)	  
  call RFT2(tmp_fft,n,1)
      
  do i=1,n
    out_fft(i) = cmplx(tmp_fft(2*i-1),tmp_fft(2*i))
  enddo
      
  return

end subroutine my_fft


end module mod_elm_apply_fft

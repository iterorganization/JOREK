subroutine broadcast_num_profiles(my_id)
!----------------------------------------------------------
! Broadcast the numerical input profiles from MPI thread 0 to the others
!----------------------------------------------------------
use tr_module
use phys_module
use mpi_mod

implicit none


integer, intent(in) :: my_id

integer :: ierr

if ( num_rho ) then
  call MPI_BCAST(num_rho_len,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  if ( my_id /= 0 ) then
     call tr_allocate(num_rho_x,1,num_rho_len,"num_rho_x",CAT_UNKNOWN)
     call tr_allocate(num_rho_y0,1,num_rho_len,"num_rho_y0",CAT_UNKNOWN)
     call tr_allocate(num_rho_y1,1,num_rho_len,"num_rho_y1",CAT_UNKNOWN)
     call tr_allocate(num_rho_y2,1,num_rho_len,"num_rho_y2",CAT_UNKNOWN)
     call tr_allocate(num_rho_y3,1,num_rho_len,"num_rho_y3",CAT_UNKNOWN)
  end if
  call MPI_BCAST(num_rho_x,num_rho_len,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(num_rho_y0,num_rho_len,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(num_rho_y1,num_rho_len,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(num_rho_y2,num_rho_len,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(num_rho_y3,num_rho_len,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
end if

if ( num_rot ) then
  call MPI_BCAST(num_rot_len,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  if ( my_id /= 0 ) then
     call tr_allocate(num_rot_x,1,num_rot_len,"num_rot_x",CAT_UNKNOWN)
     call tr_allocate(num_rot_y0,1,num_rot_len,"num_rot_y0",CAT_UNKNOWN)
     call tr_allocate(num_rot_y1,1,num_rot_len,"num_rot_y1",CAT_UNKNOWN)
     call tr_allocate(num_rot_y2,1,num_rot_len,"num_rot_y2",CAT_UNKNOWN)
     call tr_allocate(num_rot_y3,1,num_rot_len,"num_rot_y3",CAT_UNKNOWN)
  end if
  call MPI_BCAST(num_rot_x,num_rot_len,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(num_rot_y0,num_rot_len,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(num_rot_y1,num_rot_len,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(num_rot_y2,num_rot_len,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(num_rot_y3,num_rot_len,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
end if

if ( num_rhon ) then
  call MPI_BCAST(num_rho_len,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  if ( my_id /= 0 ) then
     call tr_allocate(num_rhon_x,1,num_rhon_len,"num_rhon_x",CAT_UNKNOWN)
     call tr_allocate(num_rhon_y0,1,num_rhon_len,"num_rhon_y0",CAT_UNKNOWN)
     call tr_allocate(num_rhon_y1,1,num_rhon_len,"num_rhon_y1",CAT_UNKNOWN)
     call tr_allocate(num_rhon_y2,1,num_rhon_len,"num_rhon_y2",CAT_UNKNOWN)
     call tr_allocate(num_rhon_y3,1,num_rhon_len,"num_rhon_y3",CAT_UNKNOWN)
  end if
  call MPI_BCAST(num_rhon_x,num_rhon_len,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(num_rhon_y0,num_rhon_len,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(num_rhon_y1,num_rhon_len,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(num_rhon_y2,num_rhon_len,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(num_rhon_y3,num_rhon_len,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
end if


if ( num_T ) then
  call MPI_BCAST(num_T_len,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  if ( my_id /= 0 ) then
     call tr_allocate(num_T_x,1,num_T_len,"num_T_x",CAT_UNKNOWN)
     call tr_allocate(num_T_y0,1,num_T_len,"num_T_y0",CAT_UNKNOWN)
     call tr_allocate(num_T_y1,1,num_T_len,"num_T_y1",CAT_UNKNOWN)
     call tr_allocate(num_T_y2,1,num_T_len,"num_T_y2",CAT_UNKNOWN)
     call tr_allocate(num_T_y3,1,num_T_len,"num_T_y3",CAT_UNKNOWN)
  end if
  call MPI_BCAST(num_T_x,num_T_len,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(num_T_y0,num_T_len,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(num_T_y1,num_T_len,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(num_T_y2,num_T_len,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(num_T_y3,num_T_len,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
end if

if ( num_Ti ) then
  call MPI_BCAST(num_Ti_len,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  if ( my_id /= 0 ) then
     call tr_allocate(num_Ti_x,1,num_Ti_len,"num_Ti_x",CAT_UNKNOWN)
     call tr_allocate(num_Ti_y0,1,num_Ti_len,"num_Ti_y0",CAT_UNKNOWN)
     call tr_allocate(num_Ti_y1,1,num_Ti_len,"num_Ti_y1",CAT_UNKNOWN)
     call tr_allocate(num_Ti_y2,1,num_Ti_len,"num_Ti_y2",CAT_UNKNOWN)
     call tr_allocate(num_Ti_y3,1,num_Ti_len,"num_Ti_y3",CAT_UNKNOWN)
  end if
  call MPI_BCAST(num_Ti_x,num_Ti_len,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(num_Ti_y0,num_Ti_len,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(num_Ti_y1,num_Ti_len,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(num_Ti_y2,num_Ti_len,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(num_Ti_y3,num_Ti_len,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
end if

if ( num_Te ) then
  call MPI_BCAST(num_Te_len,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  if ( my_id /= 0 ) then
     call tr_allocate(num_Te_x,1,num_Te_len,"num_Te_x",CAT_UNKNOWN)
     call tr_allocate(num_Te_y0,1,num_Te_len,"num_Te_y0",CAT_UNKNOWN)
     call tr_allocate(num_Te_y1,1,num_Te_len,"num_Te_y1",CAT_UNKNOWN)
     call tr_allocate(num_Te_y2,1,num_Te_len,"num_Te_y2",CAT_UNKNOWN)
     call tr_allocate(num_Te_y3,1,num_Te_len,"num_Te_y3",CAT_UNKNOWN)
  end if
  call MPI_BCAST(num_Te_x,num_Te_len,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(num_Te_y0,num_Te_len,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(num_Te_y1,num_Te_len,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(num_Te_y2,num_Te_len,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(num_Te_y3,num_Te_len,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
end if

#ifdef fullmhd
  call MPI_BCAST(num_Fprofile_len,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  if ( my_id /= 0 ) then
     if ( allocated(num_Fprofile_x ) ) deallocate( num_Fprofile_x  )
     if ( allocated(num_Fprofile_y0) ) deallocate( num_Fprofile_y0 )
     if ( allocated(num_Fprofile_y1) ) deallocate( num_Fprofile_y1 )
     if ( allocated(num_Fprofile_y2) ) deallocate( num_Fprofile_y2 )
     if ( allocated(num_Fprofile_y3) ) deallocate( num_Fprofile_y3 )
     call tr_allocate(num_Fprofile_x,1,num_Fprofile_len,"num_Fprofile_x",CAT_UNKNOWN)
     call tr_allocate(num_Fprofile_y0,1,num_Fprofile_len,"num_Fprofile_y0",CAT_UNKNOWN)
     call tr_allocate(num_Fprofile_y1,1,num_Fprofile_len,"num_Fprofile_y1",CAT_UNKNOWN)
     call tr_allocate(num_Fprofile_y2,1,num_Fprofile_len,"num_Fprofile_y2",CAT_UNKNOWN)
     call tr_allocate(num_Fprofile_y3,1,num_Fprofile_len,"num_Fprofile_y3",CAT_UNKNOWN)
  end if
  call MPI_BCAST(num_Fprofile_x,num_Fprofile_len,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(num_Fprofile_y0,num_Fprofile_len,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(num_Fprofile_y1,num_Fprofile_len,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(num_Fprofile_y2,num_Fprofile_len,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(num_Fprofile_y3,num_Fprofile_len,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
#endif

if ( num_ffprime ) then
  call MPI_BCAST(num_ffprime_len,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  if ( my_id /= 0 ) then
     if ( allocated(num_ffprime_x ) ) deallocate( num_ffprime_x  )
     if ( allocated(num_ffprime_y0) ) deallocate( num_ffprime_y0 )
     if ( allocated(num_ffprime_y1) ) deallocate( num_ffprime_y1 )
     if ( allocated(num_ffprime_y2) ) deallocate( num_ffprime_y2 )
     call tr_allocate(num_ffprime_x,1,num_ffprime_len,"num_ffprime_x",CAT_UNKNOWN)
     call tr_allocate(num_ffprime_y0,1,num_ffprime_len,"num_ffprime_y0",CAT_UNKNOWN)
     call tr_allocate(num_ffprime_y1,1,num_ffprime_len,"num_ffprime_y1",CAT_UNKNOWN)
     call tr_allocate(num_ffprime_y2,1,num_ffprime_len,"num_ffprime_y2",CAT_UNKNOWN)
  end if
  call MPI_BCAST(num_ffprime_x,num_ffprime_len,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(num_ffprime_y0,num_ffprime_len,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(num_ffprime_y1,num_ffprime_len,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(num_ffprime_y2,num_ffprime_len,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
end if

if ( num_d_perp ) then
  call MPI_BCAST(num_d_perp_len,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  if ( my_id /= 0 ) then
     call tr_allocate(num_d_perp_x,1,num_d_perp_len,"num_d_perp_x")
     call tr_allocate(num_d_perp_y,1,num_d_perp_len,"num_d_perp_y")
  end if
  call MPI_BCAST(num_d_perp_x,num_d_perp_len,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(num_d_perp_y,num_d_perp_len,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
end if

if ( num_d_perp_imp ) then
  call MPI_BCAST(num_d_perp_len_imp,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  if ( my_id /= 0 ) then
    call tr_allocate(num_d_perp_x_imp,1,num_d_perp_len_imp,"num_d_perp_x_imp")
    call tr_allocate(num_d_perp_y_imp,1,num_d_perp_len_imp,"num_d_perp_y_imp")
  end if
  call MPI_BCAST(num_d_perp_x_imp,num_d_perp_len_imp,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(num_d_perp_y_imp,num_d_perp_len_imp,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
end if

if ( num_zk_perp ) then
  call MPI_BCAST(num_zk_perp_len,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  if ( my_id /= 0 ) then
     call tr_allocate(num_zk_perp_x,1,num_zk_perp_len,"num_zk_perp_x")
     call tr_allocate(num_zk_perp_y,1,num_zk_perp_len,"num_zk_perp_y")
  end if
  call MPI_BCAST(num_zk_perp_x,num_zk_perp_len,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(num_zk_perp_y,num_zk_perp_len,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
end if

if ( jorek_model == 400 || jorek_model == 502 ) then
  if ( num_zk_e_perp ) then
    call MPI_BCAST(num_zk_e_perp_len,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    if ( my_id /= 0 ) then
       call tr_allocate(num_zk_e_perp_x,1,num_zk_e_perp_len,"num_zk_e_perp_x")
       call tr_allocate(num_zk_e_perp_y,1,num_zk_e_perp_len,"num_zk_e_perp_y")
    end if
    call MPI_BCAST(num_zk_e_perp_x,num_zk_e_perp_len,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(num_zk_e_perp_y,num_zk_e_perp_len,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  end if
  
  if ( num_zk_i_perp ) then
    call MPI_BCAST(num_zk_i_perp_len,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    if ( my_id /= 0 ) then
       call tr_allocate(num_zk_i_perp_x,1,num_zk_i_perp_len,"num_zk_i_perp_x")
       call tr_allocate(num_zk_i_perp_y,1,num_zk_i_perp_len,"num_zk_i_perp_y")
    end if
    call MPI_BCAST(num_zk_i_perp_x,num_zk_i_perp_len,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(num_zk_i_perp_y,num_zk_i_perp_len,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  end if
end if

if (NEO) then
   if ( num_neo_file ) then
      call MPI_BCAST(num_neo_len,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      if ( my_id /= 0 ) then
         call tr_allocate(num_neo_psi,1, num_neo_len, "num_neo_psi")
         call tr_allocate(num_aki_value,1, num_neo_len, "num_aki_value")
         call tr_allocate(num_amu_value,1, num_neo_len, "num_amu_value")
      endif
      call MPI_BCAST (num_neo_psi,num_neo_len,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST (num_aki_value,num_neo_len,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST (num_amu_value,num_neo_len,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
   end if
endif

return
end subroutine broadcast_num_profiles

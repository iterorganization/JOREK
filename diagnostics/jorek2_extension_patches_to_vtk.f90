!> Program to export extension patches to vtk files.
program jorek2_extension_patches_to_vtk

use phys_module

implicit none

integer :: my_id, k_tor, i_node, ierr, i_ext, i_block, i
integer :: count, count_save, n_all_points, n_points
integer :: i_bnd_beg, i_bnd_end
real*8  :: diff_min_beg, diff_min_end, diff
character*256 :: filename_each
character*50  :: char_tmp


write(*,*) '*************************************'
write(*,*) '* jorek2_extension_patches_to_vtk   *'
write(*,*) '*************************************'

! --- Initialise input parameters and read the input namelist.
my_id     = 0
call initialise_parameters(my_id, "__NO_FILENAME__")








! --- Open a vtk file with the wall contour
open(unit=2,file='wall_contour.vtk', ACTION = 'write')

! --- VTK Header for text files
write(2,'(A)') '# vtk DataFile Version 1.0'
write(2,'(A)') 'Line representation of vtk'
write(2,'(A)') 'ASCII'
write(2,'(A)') ''
write(2,'(A)') 'DATASET POLYDATA'

! --- Count total number of points
write(2,'(A,i8,A)') 'POINTS  ',n_limiter, ' float'

! --- Write all the points into the main VTK file
do i = 1,n_limiter
  write(2,'(A,sp,f21.11,A,sp,f21.11,A,sp,f21.11)') '     ',R_limiter(i),'     ',Z_limiter(i),'    ', 0.d0
enddo

! --- Write all the lines the main VTK file
write(2,'(A)') ''
write(2,'(A,i8,i8)') 'LINES    ',  1,    n_limiter + 2
write(2,'(i8,A)',advance='no') n_limiter+1,' '
count = 0
do i_block = 1,n_limiter
  write(2,'(i8,A)',advance='no') count, ' '
  count = count + 1
enddo
write(2,'(i8)') 0

close(2)









! --- Open a vtk file with all the patches together
open(unit=2,file='extension_patches_all.vtk', ACTION = 'write')

! --- VTK Header for text files
write(2,'(A)') '# vtk DataFile Version 1.0'
write(2,'(A)') 'Line representation of vtk'
write(2,'(A)') 'ASCII'
write(2,'(A)') ''
write(2,'(A)') 'DATASET POLYDATA'

! --- Count total number of points
n_all_points = 0
do i_ext = 1,n_wall_blocks
  n_all_points = n_all_points + n_block_points_left(i_ext) + n_block_points_right(i_ext)
enddo
write(2,'(A,i8,A)') 'POINTS  ',n_all_points, ' float'


! --- Write all the points into the main VTK file
do i_ext = 1,n_wall_blocks
  do i_block = 1,n_block_points_left(i_ext)
    write(2,'(A,sp,f21.11,A,sp,f21.11,A,sp,f21.11)') '     ',R_block_points_left(i_ext,i_block),'     ',Z_block_points_left(i_ext,i_block),'    ', 0.d0
  enddo
  do i_block = n_block_points_right(i_ext),1,-1
    write(2,'(A,sp,f21.11,A,sp,f21.11,A,sp,f21.11)') '     ',R_block_points_right(i_ext,i_block),'     ',Z_block_points_right(i_ext,i_block),'    ', 0.d0
  enddo
enddo

! --- Write all the lines the main VTK file
write(2,'(A)') ''
write(2,'(A,i8,i8)') 'LINES    ',  n_wall_blocks,    n_all_points + 2*n_wall_blocks
count = 0
do i_ext = 1,n_wall_blocks
  write(2,'(i8,A)',advance='no') n_block_points_left(i_ext)+n_block_points_right(i_ext)+1,' '
  count_save = count
  do i_block = 1,n_block_points_left(i_ext)
    write(2,'(i8,A)',advance='no') count, ' '
    count = count + 1
  enddo
  do i_block = n_block_points_right(i_ext),1,-1
    write(2,'(i8,A)',advance='no') count, ' '
    count = count + 1
  enddo
  write(2,'(i8)') count_save
enddo

close(2)









! --- Now we do the same but for each patch independently
do i_ext = 1,n_wall_blocks
  ! --- Create filename
  write(char_tmp,'(i0.2)')i_ext
  write(filename_each,'(A18,A,A4)')'extension_patches_',trim(char_tmp),'.vtk'
  open(unit=2,file=trim(filename_each), ACTION = 'write')

  ! --- VTK Header for text files
  write(2,'(A)') '# vtk DataFile Version 1.0'
  write(2,'(A)') 'Line representation of vtk'
  write(2,'(A)') 'ASCII'
  write(2,'(A)') ''
  write(2,'(A)') 'DATASET POLYDATA'
  
  ! --- Count total number of points
  n_all_points = n_block_points_left(i_ext) + n_block_points_right(i_ext)
  write(2,'(A,i8,A)') 'POINTS  ',n_all_points, ' float'
  
  ! --- Write all the points into the main VTK file
  do i_block = 1,n_block_points_left(i_ext)
    write(2,'(A,sp,f21.11,A,sp,f21.11,A,sp,f21.11)') '     ',R_block_points_left(i_ext,i_block),'     ',Z_block_points_left(i_ext,i_block),'    ', 0.d0
  enddo
  do i_block = n_block_points_right(i_ext),1,-1
    write(2,'(A,sp,f21.11,A,sp,f21.11,A,sp,f21.11)') '     ',R_block_points_right(i_ext,i_block),'     ',Z_block_points_right(i_ext,i_block),'    ', 0.d0
  enddo
  
  ! --- Write all the lines the main VTK file
  write(2,'(A)') ''
  write(2,'(A,i8,i8)') 'LINES    ',  1,    n_all_points + 2
  write(2,'(i8,A)',advance='no') n_all_points+1,' '
  count = 0
  do i_block = 1,n_block_points_left(i_ext)
    write(2,'(i8,A)',advance='no') count, ' '
    count = count + 1
  enddo
  do i_block = n_block_points_right(i_ext),1,-1
    write(2,'(i8,A)',advance='no') count, ' '
    count = count + 1
  enddo
  write(2,'(i8)') 0

  close(2)
enddo






end program jorek2_extension_patches_to_vtk








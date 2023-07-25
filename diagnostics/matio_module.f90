!> Set of subroutines for custom diagnistics
module matio_module
  use hdf5_io_module
  use iso_c_binding, only: c_double
  implicit none

  logical :: fileexist=.false.

  private
  public :: read_matrix_h5, save_mat_h5, timestamp, slurmid

  contains 

!> Reads sparse matrix from HDF5 file
  subroutine read_matrix_h5(fname,n,nnz,rowptr,colptr,val,rhs)
    integer :: i,stat, ierr
    integer(HID_T) fid

    character(len=255),intent(in) :: fname
    integer,intent(out) :: nnz,n
    integer :: indexing

    integer, dimension(:), pointer :: rowptr
    integer, dimension(:), pointer :: colptr
    real(kind=c_double), dimension(:), pointer    :: val
    real(kind=c_double), dimension(:), pointer    :: rhs   

    write(*,*) "Reading ",trim(fname)

    call HDF5_open(trim(fname),fid,ierr)
    if (ierr/=0) return
    call HDF5_integer_reading(fid,n,"n")  
    call HDF5_integer_reading(fid,nnz,"nnz")        
    call HDF5_integer_reading(fid,indexing,"indexing")              
     
    allocate(colptr(nnz),rowptr(nnz))      
    allocate(val(nnz))
    allocate(rhs(n))

    call HDF5_array1D_reading_int(fid,rowptr,"irn")
    call HDF5_array1D_reading_int(fid,colptr,"jcn")      
    call HDF5_array1D_reading(fid,val,"val")            
    call HDF5_array1D_reading(fid,rhs,"rhs")
    call HDF5_close(fid)

! change to C indexing
    if (indexing==1) then
      do i=1,nnz
        rowptr(i)=rowptr(i)-1
        colptr(i)=colptr(i)-1                
      enddo
    endif
    return
  end subroutine read_matrix_h5  

!> Save sparse matrix into HDF5 file
  subroutine save_mat_h5(fname, ng, nl, nnz, irn, jcn, val, l2g, rhs, ind_min, ind_max, block_size)
    use phys_module, only: n_tor, n_var
    integer :: rank, ng, nl, nt, nv, nd, nnz, ierr
    integer(HID_T) fid
    CHARACTER(LEN=10)                                    :: fname
    integer, dimension(:), pointer                       :: irn, jcn
    real(kind=c_double), dimension(:), pointer           :: val
    integer, dimension(:), pointer, optional             :: l2g
    real(kind=c_double), dimension(:), pointer, optional :: rhs
    integer, optional             :: ind_min, ind_max, block_size

    !write(fname,'(A5,I2.2,A3)') "matA_",rank,".h5"
    fname = trim(fname)

    call HDF5_create(filename=fname,file_id=fid,ierr=ierr)
    call HDF5_integer_saving(fid,ng,'ng')
    call HDF5_integer_saving(fid,nl,'nl')
    call HDF5_integer_saving(fid,n_tor,'ntor')
    call HDF5_integer_saving(fid,n_var,'nvar')
    call HDF5_integer_saving(fid,nnz,'nnz')
    call HDF5_integer_saving(fid,1,'indexing')
    call HDF5_array1D_saving_int(fid,irn,nnz,'irn')
    call HDF5_array1D_saving_int(fid,jcn,nnz,'jcn')
    call HDF5_array1D_saving(fid,val,nnz,'val')
    if (present(l2g)) call HDF5_array1D_saving_int(fid,l2g,nl,'loc2glob')
    if (present(rhs)) call HDF5_array1D_saving(fid,rhs,nl,'rhs')
    if (present(ind_min)) call HDF5_integer_saving(fid,ind_min,'ind_min')
    if (present(ind_max)) call HDF5_integer_saving(fid,ind_max,'ind_max')
    if (present(block_size)) call HDF5_integer_saving(fid,block_size,'block_size')
    call HDF5_close(fid)
    return

  end subroutine save_mat_h5

!> Saves 1D array into HDF5 file
  subroutine save_solution_h5(fname,n,x)
      integer(HID_T) fid
      integer :: i,stat, ierr
      integer, intent(in) :: n
      character(len=255),intent(in) :: fname
      real*8, dimension(:), intent(in) :: x

      write(*,*) "Saving ",trim(fname)

      call HDF5_create(trim(fname),fid,ierr)
      call HDF5_integer_saving(fid,n,"n")  
      call HDF5_array1D_saving(fid,x,n,"x")            
      call HDF5_close(fid)        

      return
  end subroutine save_solution_h5

!> Saves SLURM_PROCID variable into file
  subroutine slurmid(rank)
    character(len=12) :: envname="SLURM_PROCID"
    character(len=4) :: val
    integer :: rank

    call get_environment_variable (envname, val)
    open(unit = 101, file = 'procid.out', status='REPLACE', action='WRITE')
    write(101,*) val, rank
    close(101)

  end subroutine slurmid

!> Saves timestamp and message into timeline file
  subroutine timestamp(msg,id)
    character(len=*), intent(in) :: msg
    integer, intent(in), optional :: id

    integer,dimension(8) :: values
    character(len=16) :: fname
    real :: t 
    
    if (present(id)) then
      write (fname, "(A8,(I0.4),A4)") 'timeline', id, '.out'
    else
      write (fname, "(A12)") 'timeline.out'
    endif
    
    inquire(file=fname, exist=fileexist)
    if (.not.fileexist) open(unit = 1001, file = trim(fname), status='REPLACE', action='WRITE')


    call date_and_time(VALUES=values)
    t = values(5)*3600000 + values(6)*60000 + values(7)*1000 + values(8)
    open(unit = 1001, file = trim(fname), status='OLD', position="append", action='WRITE')
    write(1001,*) t, trim(msg) 
    close(1001)
    return

  end subroutine timestamp


end module matio_module

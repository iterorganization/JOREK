!> Contains subroutines to sort/remove duplicates for equilibrium solve
!! and fast column sorting for PC sparse matrix as needed for STRUMPACK
module sorting_module

  use iso_c_binding
  use mod_integer_types
  implicit none
  private
  public remove_duplicates, convert2csr, convert_sorting

#define INTSIZE 8
#define CINT c_int64_t

interface

  subroutine qsort(array, elem_count, elem_size, compar) bind(C,name="qsort")
    import
    type(c_ptr),value       :: array
    integer(C_INT_ALL),value :: elem_count
    integer(C_SIZE_T),value :: elem_size
    type(c_funptr),value    :: compar
  end subroutine qsort

  subroutine convert2csr(indx, n, m, nnz, irn, jcn, val) bind(C)
    use iso_c_binding
    use mod_integer_types
    implicit none
    type(C_PTR) :: irn, jcn, val
    integer(kind=C_INT_ALL), intent(in) :: n, m, indx
    integer(kind=C_INT_ALL), intent(inout) :: nnz
  end subroutine convert2csr

  subroutine sortunique(nnz,ijn) bind(C)
    use iso_c_binding
    use mod_integer_types
    implicit none
    integer(kind=C_INT_ALL) :: nnz
    integer(kind=CINT),dimension(:), pointer, intent(inout) :: ijn
  end subroutine sortunique

end interface

contains

  integer(2) function compar(a, b) bind(C)
    use iso_c_binding
    integer(kind=CINT) a, b

    if ( a .lt. b ) compar = -1
    if ( a .eq. b ) compar = 0
    if ( a .gt. b ) compar = 1
  end function compar

  !> Sort and remove duplicates from 1D array using qsort
  !! replace list with uniquelly sorted entries; return number of unique elements
  subroutine unique_sorted(array,n)

    use, intrinsic :: iso_c_binding
    implicit none

    integer(kind=INTSIZE), dimension(:), pointer :: array
    integer(kind=int_all), intent(inout) :: n
    integer(kind=C_INT_ALL) :: array_len
    integer(C_SIZE_T) array_size
    logical,allocatable :: duplicates(:)
    integer(kind=INTSIZE) :: m, i, j
    integer :: cc, cr
    real t0, t1

    call system_clock(count=cc, count_rate=cr); t0 =  real(cc)/cr

    array_len = n
    array_size = INTSIZE
    call qsort(c_loc(array(1)), array_len, array_size, c_funloc(compar))

    allocate(duplicates(array_len)); duplicates=.false.
    duplicates(1:array_len)=array(1:array_len-1).eq.array(2:array_len)

    m = count(duplicates)
    if (m.gt.0) then
      j = 1
      do i=1, array_len
        if (.not.duplicates(i)) then
          array(j) = array(i)
          j = j + 1
        endif
      enddo
      array_len = j - 1
    endif

    call system_clock(count=cc, count_rate=cr); t1 =  real(cc)/cr
    write(*,*) "Sorting time (s) =",t1-t0

    return

  end subroutine unique_sorted

  !> Find index of element x in the list
  recursive function find_index(list,low,high,x) result(idx)

    use mod_integer_types

    integer(kind=INTSIZE), intent(in) :: x
    integer(kind=int_all), intent(in) :: low, high
    integer(kind=INTSIZE), dimension(:), pointer :: list(:)
    integer(kind=int_all) :: mid
    integer(kind=INTSIZE) :: idx

    if (low.gt.high) then
      idx = 0
      write(*,*) "Error in find_index: element not found", x, low, high
      call exit(0)
    endif

    mid = (low + high)/2

    ! target value is found
    if (x .eq. list(mid)) then
      idx = mid

    ! discard all elements in the right search space
    ! including the mid element
    elseif (x .lt. list(mid)) then
      idx = find_index(list, low,  mid - 1, x)

    ! discard all elements in the left search space
    ! including the mid element
    else
      idx = find_index(list, mid + 1, high, x)
    endif

  end function find_index

  !> Sort and remove duplicates from sparse matrix
  !! by converting to 1D array of ij index
  subroutine remove_duplicates(n,nnz,irn,jcn,val)

    use, intrinsic :: iso_c_binding
    use mod_integer_types

    integer(kind=int_all), intent(in) :: n
    integer(kind=int_all), intent(inout) :: nnz
    integer(kind=C_INT_ALL), dimension(:), pointer  :: irn, jcn
    real(kind=C_DOUBLE), dimension(:), pointer  :: val

    real(kind=C_DOUBLE), allocatable :: val_new(:)
    ! long integer is required for 1d representation of coordinate index
    integer(kind=INTSIZE), dimension(:), pointer :: ij, ij_new, new_ind
    integer(kind=INTSIZE) :: dum, i1, i2, i3
    integer(kind=int_all) :: i, j, nnz0, indmin, indmax

    allocate(ij(nnz), new_ind(nnz))
    do i = 1, nnz
      i1 = int(irn(i)-1,kind=INTSIZE)
      i2 = int(n,kind=INTSIZE)
      i3 = int(jcn(i),kind=INTSIZE)
      ij(i) = i1*i2 + i3
    enddo

    nnz0 = nnz
    call sortunique(nnz,ij)
    !call unique_sorted(ij,nnz)

    if (nnz.ne.nnz0) write(*,*) "Number of nnz changed: nnz_old, nnz_new = ", nnz0, nnz

    ! find index of original element in the new (ordered) list
    i1 = int(n,kind=INTSIZE)
    indmin = 1; indmax = nnz;
    do i = 1, nnz0
      i1 = int(irn(i)-1,kind=INTSIZE)
      i2 = int(n,kind=INTSIZE)
      i3 = int(jcn(i),kind=INTSIZE)
      dum = i1*i2 + i3
      new_ind(i) = find_index(ij,indmin,indmax,dum)
    enddo

    ! sum up possible duplicates
    allocate(val_new(nnz)); val_new = 0.0
    do i =1, nnz0
      val_new(new_ind(i)) = val_new(new_ind(i)) + val(i)
    enddo

    do i = 1, nnz
      irn(i) = int((ij(i)-1)/n) + 1
      jcn(i) = mod(ij(i)-1,n) + 1
      val(i) = val_new(i)
    enddo

    deallocate(ij,new_ind,val_new)

  end subroutine remove_duplicates

  !> Convert to CSR while sorting column-wise
  !! based on matrix being structured in consecutive (non-uniform)
  !! blocks of irn values
  subroutine convert_sorting(nnz,irn,jcn,val,block_size,indx)

    use, intrinsic :: iso_c_binding
    use mod_integer_types

    integer(kind=int_all), intent(in) :: nnz
    integer, intent(in) :: indx, block_size
    integer(kind=int_all), dimension(:), pointer  :: irn, jcn
    real(kind=c_double), dimension(:), pointer  :: val

    integer(kind=int_all), dimension(:), allocatable :: jcn_tmp, indmin, indmax, iblock, iptr
    real(kind=c_double),  dimension(:), allocatable :: val_tmp

    integer(kind=int_all) :: i, nloc, n1, n2, ni, irn0, cnt, idum
    integer :: n_irn_block, ib

    logical :: check

    integer :: cc, cr
    real t0, t1

    call system_clock(count=cc, count_rate=cr); t0 =  real(cc)/cr

    irn0 = minval(irn(1:nnz))
    nloc = maxval(irn(1:nnz)) - irn0 + 1
    !write(*,*) minval(irn(1:nnz)), maxval(irn(1:nnz)), "nloc", nloc
    irn(1:nnz) = irn(1:nnz) - irn0 + 1 ! irn used as index

    write(*,*) "nloc", nloc, "nnz", nnz, "residue", mod(nnz,block_size)
        
    allocate(indmin(nloc), indmax(nloc), iptr(nloc+1))
    iptr = 0

    indmin = nnz
    indmax = 1
    iptr(1) = 1
    do i = 1, nnz, block_size
      iptr(irn(i) + 1) = iptr(irn(i) + 1) + block_size
      indmin(irn(i)) = min(indmin(irn(i)),i)
      indmax(irn(i)) = max(indmax(irn(i)),i+block_size-1)
    enddo

    do i = 2, nloc+1
      iptr(i) = iptr(i) + iptr(i-1)
    enddo
    if ((iptr(nloc+1)-1).ne.nnz) write(*,*) "Warning: iptr(nloc+1)", iptr(nloc+1)-1

    !write(*,*) indmin(nloc), indmax(nloc)

    ! determine number of irn-blocks
    n_irn_block = 1
    do idum = 2, nloc
      if (indmin(idum).gt.indmax(idum-1)) n_irn_block = n_irn_block + 1
    enddo

    allocate(iblock(n_irn_block+1))
    iblock(1) = 1; iblock(n_irn_block+1) = nloc + 1
    ib = 2
    do idum = 2, nloc
      if (indmin(idum).gt.indmax(idum-1)) then
        iblock(ib) = idum ! min irn belonging to block
        ib = ib + 1
      endif
    enddo
    
    ! find maximal block size for temporary buffer allocation
    ni = 0
    do ib = 1, n_irn_block
      cnt = 1
      n1 = indmin(iblock(ib))
      n2 = indmax(iblock(ib+1)-1)
      ni = max(ni, n2 - n1 + 1)
    enddo
    allocate(jcn_tmp(ni),val_tmp(ni))

    do ib = 1, n_irn_block
      cnt = 1
      n1 = indmin(iblock(ib))
      n2 = indmax(iblock(ib+1)-1)
      ni = n2 - n1 + 1
      do idum = iblock(ib),iblock(ib+1)-1
        do i = indmin(idum), indmax(idum), block_size
          if (irn(i).eq.idum) then
            jcn_tmp(cnt:cnt + block_size - 1) = jcn(i:i + block_size - 1)
            val_tmp(cnt:cnt + block_size - 1) = val(i:i + block_size - 1)
            cnt = cnt + block_size
          endif
        enddo
      enddo
      jcn(n1:n2) = jcn_tmp(1:ni)
      val(n1:n2) = val_tmp(1:ni)
    enddo
    
    deallocate(jcn_tmp,val_tmp)

    irn(1:nloc+1) = iptr(1:nloc+1)
    
    deallocate(indmin, indmax, iptr)
    deallocate(iblock)
    
    ! check sorting consistency
    if (.false.) then
      do n1 = 1, nloc
        check = .true.
        ni = irn(n1+1) - irn(n1)
        do i = irn(n1+1)-ni+1, irn(n1+1)-1
          if (jcn(i).le.jcn(i-1)) check = .false.
          idum = n1
        enddo
        if (.not.check) exit
      enddo
      write(*,*) "Consistency:", idum, check
      if (.not.check) then
        ni = irn(n1+1) - irn(n1)
        write(*,*) "idum", n1, "iptr", ni, "jcn", jcn(irn(n1+1)-ni:irn(n1+1)-1)
      endif
    endif

    call system_clock(count=cc, count_rate=cr); t1 =  real(cc)/cr
    write(*,*) "Sorting/csr time (s) =",t1-t0

  end subroutine convert_sorting

end module sorting_module

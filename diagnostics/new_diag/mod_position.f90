!> This module contains datastructures describing lists of poloidal and toroidal positions.
!! 
!! Together with mod_expression, this provides a general diagnostic framework for many applications.
module mod_position
  
  
  
  
  
  use constants
  use mod_parameters
  use equil_info
  use data_structure
  use mod_straight_field_line
  use mod_interp
  
  
  
  
  
  implicit none
  
  
  
  
  
  public
  
  
  
  
  
  ! --- Constants
  character(len=12), parameter, private :: THIS_MOD_NAME = 'mod_position'
  
  
  
  
  
  ! --- Data structures
  
  !> Data structure for a poloidal position
  type t_pol_pos
    logical            :: outside = .false. !< Outside computational domain?
    real*8             :: R, R_s, R_t, R_st, R_ss, R_tt, Z, Z_s, Z_t, Z_st, Z_ss, Z_tt, s, t
    integer            :: ielm
    type(type_element) :: element
    type(type_node)    :: nodes(n_vertex_max)
    ! --- The following quantities will only be available in certain cases of poloidal positions
    real*8             :: theta_star !< Straight field line angle (for flux surfaces)
    real*8             :: r_minor    !< Minor radius from A = r_minor^2 pi (for flux surfaces)
    real*8             :: length     !< Length along a line (for line)
    ! --- Quantities related to boundary elements
    real*8             :: bnd_normal(2) = 0.d0 !< Normal vector to the computational boundary (pointing outside)
  end type t_pol_pos
  
  !> Data structure for a list of poloidal positions
  type t_pol_pos_list
    type(t_pol_pos), allocatable :: pos(:,:)
    integer :: n_pos(2) = 0
    logical :: full_turn !< Do the positions cover a full poloidal turn?
  end type t_pol_pos_list
  
  !> Data structure for a toroidal position
  type t_tor_pos
    real*8  :: phi
  end type t_tor_pos
  
  !> Data structure for a list of toroidal positions
  type t_tor_pos_list
    type(t_tor_pos), allocatable :: pos(:)
    integer :: n_pos = 0
    logical :: full_period !< Do the positions cover a full toroidal period?
  end type t_tor_pos_list
  
  
  
  
  
  contains
  
  
  
  
  
  !> Clean up a poloidal position data structure.
  subroutine cleanup_pol_pos(pos_list)
    
    type(t_pol_pos_list), intent(inout) :: pos_list
    
    pos_list%n_pos = 0
    
    if ( allocated(pos_list%pos) ) deallocate( pos_list%pos )
    
  end subroutine cleanup_pol_pos
  
  
  
  
  
  !> Clean up a toroidal position data structure.
  subroutine cleanup_tor_pos(pos_list)
    
    type(t_tor_pos_list), intent(inout) :: pos_list
    
    pos_list%n_pos       = 0
    pos_list%full_period = .false.
    
    if ( allocated(pos_list%pos) ) deallocate( pos_list%pos )
    
  end subroutine cleanup_tor_pos
  
  
  
  
  
  !> Allocate a poloidal position data structure.
  subroutine alloc_pol_pos(pos_list, n_pos)
    
    type(t_pol_pos_list), intent(inout) :: pos_list
    integer,              intent(in)    :: n_pos(2)
    
    integer :: i, j
    
    call cleanup_pol_pos(pos_list)
    
    pos_list%n_pos(:) = n_pos(:)
    
    allocate( pos_list%pos(n_pos(1),n_pos(2)) )
    
    do j = 1, n_pos(2)
      do i = 1, n_pos(1)
        pos_list%pos(i,j)%theta_star = 0.d0
        pos_list%pos(i,j)%length     = 0.d0
        pos_list%pos(i,j)%r_minor    = 0.d0
      end do
    end do
    
  end subroutine alloc_pol_pos
  
  
  
  
  
  !> Allocate a toroidal position data structure.
  subroutine alloc_tor_pos(pos_list, n_pos)
    
    type(t_tor_pos_list), intent(inout) :: pos_list
    integer,              intent(in)    :: n_pos
    
    call cleanup_tor_pos(pos_list)
    
    pos_list%n_pos = n_pos
    
    allocate( pos_list%pos(n_pos) )
    
  end subroutine alloc_tor_pos
  
  
  
  
  
  !> Simple function to generate poloidal positions (wrapper for routine create_pol_pos).
  function pol_pos(node_list, element_list, eq, R, Z, ielm, s, t, Rmin, Rmax, nR, Zmin, Zmax, nZ,  &
    Rstart, Rend, Zstart, Zend, n, PsiN, nTht, PsiNmin, PsiNmax, nPsiN, nmaxsteps, deltaphi,       &
    nsmallsteps) result(pos_list)
    type(t_pol_pos_list), target :: pos_list
    
    ! --- Routine parameters
    type(type_node_list),         intent(in)    :: node_list
    type(type_element_list),      intent(in)    :: element_list
    type(t_equil_state),          intent(in)    :: eq
    real*8,  optional,            intent(in)    :: R, Z, s, t, Rmin, Rmax, Zmin, Zmax, PsiN,       &
      Rstart, Rend, Zstart, Zend, PsiNmin, PsiNmax, deltaphi
    integer, optional,            intent(in)    :: ielm, nR, nZ, n, nTht, nPsiN, nmaxsteps,        &
      nsmallsteps
    
    ! --- Local variables
    integer :: ierr
    
    call create_pol_pos(pos_list, ierr, node_list, element_list, eq, R, Z, ielm, s, t, Rmin, Rmax, &
      nR, Zmin, Zmax, nZ, Rstart, Rend, Zstart, Zend, n, PsiN, nTht, PsiNmin, PsiNmax, nPsiN,      &
      nmaxsteps, deltaphi, nsmallsteps)
    
  end function pol_pos
  
  
  
  
  
  !> Create a poloidal position data structure (one or several poloidal positions) from:
  !! - (R,Z) or (ielm,s,t)           -> single position
  !! - (Rmin,Rmax,nR,Zmin,Zmax,nZ)   -> 2D array of positions
  !! - (Rstart,Rend,Zstart,Zend,n)   -> Equidistant points along a straight line
  !! - (PsiN,nTht)                   -> flux surface (equidistant points in theta*)
  !! - (PsiNmin,PsiNmax,nPsiN,nTht)  -> flux surfaces (equidistant points in theta*)
  !! To be added: Single node; All nodes; All nodes with subdivision of elements (for vtk) ###
  recursive subroutine create_pol_pos(pos_list, ierr, node_list, element_list, eq, R, Z, ielm, s,  &
    t, Rmin, Rmax, nR, Zmin, Zmax, nZ, Rstart, Rend, Zstart, Zend, n, PsiN, nTht, PsiNmin, PsiNmax,&
    nPsiN, nmaxsteps, deltaphi, nsmallsteps)
    
    character(len=64), parameter :: THIS_ROUTINE_NAME = trim(THIS_MOD_NAME) // ':create_pol_pos'
    
    ! --- Routine parameters
    type(t_pol_pos_list), target, intent(inout) :: pos_list
    integer,                      intent(out)   :: ierr
    type(type_node_list),         intent(in)    :: node_list
    type(type_element_list),      intent(in)    :: element_list
    type(t_equil_state),          intent(in)    :: eq
    real*8,  optional,            intent(in)    :: R, Z, s, t, Rmin, Rmax, Zmin, Zmax, PsiN,       &
      Rstart, Rend, Zstart, Zend, PsiNmin, PsiNmax, deltaphi
    integer, optional,            intent(in)    :: ielm, nR, nZ, n, nTht, nPsiN, nmaxsteps,        &
      nsmallsteps
    
    ! --- Local variables
    type(t_theta_mapping) :: mapping
    type(t_pol_pos), pointer :: pos
    real*8  :: R_out, Z_out, hh, gx, gy, gg, ax, ay, full_length
    integer :: i, j
    real*8, allocatable :: surface(:) !< Poloidal surface inside flux surface (for r_minor)
    
    ierr = 0
    
    ! --- Single position given by R and Z.
    if ( present(R) .and. present(Z) ) then
      
      call alloc_pol_pos(pos_list, (/1,1/))
      pos   => pos_list%pos(1,1)
      pos%R = R
      pos%Z = Z
      call find_RZ(node_list, element_list, R, Z, R_out, Z_out, pos%ielm, pos%s, pos%t, ierr)
      pos%outside = ( ierr /= 0 )
      call fill_pol_pos(pos, node_list, element_list)
      
    ! --- Single position given by ielm, s, and t.
    else if ( present(ielm) .and. present(s) .and. present(t) ) then
      
      call alloc_pol_pos(pos_list, (/1,1/))
      pos      => pos_list%pos(1,1)
      pos%ielm = ielm
      pos%s    = s
      pos%t    = t
      call fill_pol_pos(pos, node_list, element_list)
      
    ! --- Rectangular array of positions in R and Z.
    else if ( present(Rmin) .and. present(Rmax) .and. present(nR) .and. present(Zmin) .and.        &
      present(Zmax) .and. present(nZ) ) then
      
      call alloc_pol_pos(pos_list, (/nR,nZ/))
      
      !$OMP parallel do default(none) firstprivate(R_out,Z_out,ierr) private(pos,i,j)              & 
      !$OMP shared(nR,nZ,node_list,element_list,pos_list,Rmin,Rmax,Zmin,Zmax) schedule(static)
      do i = 1, nR
        do j = 1, nZ
          pos   => pos_list%pos(i,j)
          pos%R = Rmin + (Rmax-Rmin) * real(i-1)/real(nR-1)
          pos%Z = Zmin + (Zmax-Zmin) * real(j-1)/real(nZ-1)
          ierr = 0
          call find_RZ(node_list, element_list, pos%R, pos%Z, R_out, Z_out, pos%ielm, pos%s, pos%t,&
            ierr)
          pos%outside = ( ierr /= 0 )
          call fill_pol_pos(pos, node_list, element_list)
        end do
      end do
      !$OMP end parallel do
      
    ! --- Equidistant points along a straight line in R and Z.
    else if ( present(Rstart) .and. present(Rend) .and. present(Zstart) .and. present(Zend) .and.  &
      present(n) ) then
      
      call alloc_pol_pos(pos_list, (/1,n/))
      full_length = sqrt( (Rend-Rstart)**2 + (Zend-Zstart)**2 )
      do i = 1, n
        pos   => pos_list%pos(1,i)
        pos%R = Rstart + (Rend-Rstart) * real(i-1)/real(n-1)
        pos%Z = Zstart + (Zend-Zstart) * real(i-1)/real(n-1)
        pos%length = full_length * real(i-1)/real(n-1)
        call find_RZ(node_list, element_list, pos%R, pos%Z, R_out, Z_out, pos%ielm, pos%s, pos%t,  &
          ierr)
        pos%outside = ( ierr /= 0 )
        call fill_pol_pos(pos, node_list, element_list)
      end do
      
    ! --- Equidistant points along a straight line in R at fixed Z (e.g., midplane profiles).
    else if ( present(Rstart) .and. present(Rend) .and. present(Z) .and. present(n) ) then
      
      call create_pol_pos(pos_list, ierr, node_list, element_list, eq, Rstart=Rstart, Rend=Rend,   &
        Zstart=Z, Zend=Z, n=n)
      
    ! --- Single (closed) flux surface.
    else if ( present(PsiN) .and. present(nTht) ) then
      
      call create_pol_pos(pos_list, ierr, node_list, element_list, eq, PsiNmin=PsiN, PsiNmax=PsiN, &
        nPsiN=1)
      
    ! --- Several (closed) flux surfaces.
    else if ( present(PsiNmin) .and. present(PsiNmax) .and. present(nPsiN) .and. present(nTht) )   &
      then
      
      if ( (min(PsiNmin,PsiNmax) < 0.d0) .and. (max(PsiNmin,PsiNmax) > 1.d0) ) then
        write(*,*) 'ERROR in '//trim(THIS_ROUTINE_NAME)//': PsiNmin and PsiNmax must be between 0 and 1.'
        ierr = 300
        return
      end if 
      
      call determine_theta_mag(mapping, node_list, element_list, eq, (/PsiNmin,PsiNmax/), nPsiN,   &
        nTht, ierr, nmaxsteps2=nmaxsteps, deltaphi2=deltaphi, nsmallsteps2=nsmallsteps)
        ! (Optional parameters passed on to this routine: nmaxsteps, deltaphi, nsmallsteps)
      if ( ierr /= 0 ) then
        write(*,*) 'ERROR in '//trim(THIS_ROUTINE_NAME)//' calling determine_theta_mag.'
        ierr = 400
        return
      end if 
      
      allocate(surface(nPsiN))
      surface(:) = 0.d0
      call alloc_pol_pos(pos_list, (/nTht,nPsiN/))
      pos_list%full_turn = .true.
      
      do j = 1, nTht
        do i = 1, nPsiN
          pos   => pos_list%pos(j,i)
          pos%R = mapping%rre(i,j-1)
          pos%Z = mapping%zze(i,j-1)
          pos%theta_star = 2.d0 * PI * real(j-1) / real(nTht-1) !######### check if nTht or nTht-1
          call find_RZ(node_list, element_list, pos%R, pos%Z, R_out, Z_out, pos%ielm, pos%s, pos%t,&
            ierr)
          pos%outside = ( ierr /= 0 )
          call fill_pol_pos(pos, node_list, element_list)
          
          ! --- Calculate poloidal surface inside flux surface (for r_minor)
          if ( j /= 1 ) then
            gx = pos_list%pos(j,i)%R - pos_list%pos(j-1,i)%R
            gy = pos_list%pos(j,i)%Z - pos_list%pos(j-1,i)%Z
          else ! if (j == 1) ######### check if this is necessary (it is, if first point /= last point in flux surface!)
            gx = pos_list%pos(1,i)%R - pos_list%pos(nTht,i)%R
            gy = pos_list%pos(1,i)%Z - pos_list%pos(nTht,i)%Z
          end if
          gg = sqrt( gx**2 + gy**2 )
          ax = pos_list%pos(j,i)%R - eq%R_axis
          ay = pos_list%pos(j,i)%Z - eq%Z_axis
          hh = ( gy * ax - gx * ay ) / gg
          surface(i) = surface(i) + abs(hh * gg / 2.d0)
        end do
      end do
      
      ! --- Fill in r_minor
      do j = 1, nTht
        do i = 1, nPsiN
          pos_list%pos(j,i)%r_minor = sqrt( surface(i) / PI )
        end do
      end do
      
      call cleanup_mapping(mapping)
      deallocate(surface)
      
    ! --- Several (closed) flux surfaces.
    else if ( present(nPsiN) .and. present(nTht) ) then
      
      call create_pol_pos(pos_list, ierr, node_list, element_list, eq, nTht=nTht, nPsiN=nPsiN,     &
        psiNmin=0.001d0, psiNmax=0.998d0)
      
    else
      
      ierr = 99
      write(*,*)
      write(*,*) 'ERROR in '//trim(THIS_ROUTINE_NAME)//':'
      write(*,*) 'No valid representation for poloidal position(s) provided.'
      write(*,*)
      
    end if
    
  end subroutine create_pol_pos
  
  
  
  
  
  !> Auxilliary routine used by create_pol_pos: Fill information (R, R_s, ..., Z_tt, element, nodes)
  !! for a single poloidal position. Requires that ielm, s, t are already set to correct values.
  subroutine fill_pol_pos(pos, node_list, element_list)
    
    ! --- Routine parameters
    type(t_pol_pos), pointer, intent(inout) :: pos
    type(type_node_list),     intent(in)    :: node_list
    type(type_element_list),  intent(in)    :: element_list
    
    ! --- Local variables
    integer :: i
    
    if ( (pos%ielm < 0) .or. (pos%ielm > element_list%n_elements) ) pos%outside = .true.
    if ( pos%outside ) return
    
    call interp_RZ(node_list, element_list, pos%ielm, pos%s, pos%t, pos%R, pos%R_s, pos%R_t,       &
      pos%R_st, pos%R_ss, pos%R_tt, pos%Z, pos%Z_s, pos%Z_t, pos%Z_st, pos%Z_ss, pos%Z_tt)
    pos%element = element_list%element(pos%ielm)
    
    do i = 1, n_vertex_max
      pos%nodes(i) = node_list%node(pos%element%vertex(i))
    end do
    
  end subroutine fill_pol_pos
  
  
  


  !> Function to generate positions on the boundary of the JOREK's domain 
  function bnd_pos(node_list, element_list, bnd_node_list, bnd_elm_list, n_elm_pts) result(pos_list)

    type(t_pol_pos_list), target :: pos_list
    
    ! --- Routine parameters
    type(type_node_list),         intent(in)    :: node_list
    type(type_element_list),      intent(in)    :: element_list
    type(type_bnd_node_list),     intent(in)    :: bnd_node_list
    type(type_bnd_element_list),  intent(in)    :: bnd_elm_list
    integer,                      intent(in)    :: n_elm_pts

    ! --- Local variables
    integer                  :: i_bnd, m_bndelem, mv1, m_elm, m_pt
    real*8                   :: s_or_t, s, t
    real*8                   :: R, R_s, R_t, Z, Z_s, Z_t
    real*8                   :: vec_out(2)
    logical                  :: s_const
    type(t_pol_pos), pointer :: pos

    i_bnd = 0  ! index for bnd point
  
    ! --- alllocate position list
    call alloc_pol_pos(pos_list, (/1, n_elm_pts * bnd_elm_list%n_bnd_elements /))

    ! --- For every boundary element do 
    do m_bndelem = 1, bnd_elm_list%n_bnd_elements

      mv1     = bnd_elm_list%bnd_element(m_bndelem)%side
      m_elm   = bnd_elm_list%bnd_element(m_bndelem)%element

      ! --- For every point in the element do
      do m_pt = 1, n_elm_pts

        i_bnd  = i_bnd + 1 
        s_or_t = float(m_pt-1)/float(n_elm_pts)

        ! --- Which s and t values correspond to the current point and is the
        !     boundary element an s=const or t=const side of the 2D element?
        select case (mv1)
        case (1)
          s = s_or_t;  t = 0.d0;    s_const = .false.
        case (2)
          s = 1.d0;    t = s_or_t;  s_const = .true.
        case (3)
          s = s_or_t;  t = 1.d0;    s_const = .false.
        case (4)
          s = 0.d0;    t = s_or_t;  s_const = .true.
        end select

        ! --- Fill in positions
        pos   => pos_list%pos(1,i_bnd)
        pos%ielm = m_elm
        pos%s    = s
        pos%t    = t
        call fill_pol_pos(pos, node_list, element_list)

        ! --- Normal vector to the boundary
        if ( s_const ) then
          pos%bnd_normal = (/ -pos%Z_t, pos%R_t /) / sqrt(pos%R_t**2.d0 + pos%Z_t**2.d0)  
        else
          pos%bnd_normal = (/ -pos%Z_s, pos%R_s /) / sqrt(pos%R_s**2.d0 + pos%Z_s**2.d0)
        end if

        ! --- Correct normal direction to point outwards 
        ! --- Get point inside the element
        if ( s_const ) then
          call interp_RZ(node_list, element_list, m_elm, 0.5d0,     t, R, R_s, R_t, Z, Z_s, Z_t)
        else
          call interp_RZ(node_list, element_list, m_elm,     s, 0.5d0, R, R_s, R_t, Z, Z_s, Z_t)
        endif
        vec_out = (/  pos%R - R, pos%Z  - Z/)    ! vector pointing outside the domain

        pos%bnd_normal = pos%bnd_normal*sign(1.d0, vec_out(1)*pos%bnd_normal(1)+vec_out(2)*pos%bnd_normal(2))

      enddo

    enddo 
   
  end function bnd_pos
  




  
  !> Simple function to generate toroidal positions (wrapper for routine create_tor_pos).
  function tor_pos(phi, iplane, phistart, phiend, nphi) result(pos_list)
    type(t_tor_pos_list), target :: pos_list
    
    ! --- Routine parameters
    real*8,  optional,            intent(in)    :: phi, phistart, phiend
    integer, optional,            intent(in)    :: iplane, nphi
    
    ! --- Local variables
    integer :: ierr
    
    call create_tor_pos(pos_list, ierr, phi, iplane, phistart, phiend, nphi)
    
  end function tor_pos
  
  
  
  
  
  !> Create a toroidal position data structure from the specified information:
  !! - phi or iplane            -> single toroidal position
  !! - (phistart, phiend, nphi) -> several equidistant toroidal positions
  !! - nphi                     -> distribute nphi positions over full toroidal turn (or period)
  subroutine create_tor_pos(pos_list, ierr, phi, iplane, phistart, phiend, nphi)
    
    character(len=64), parameter :: THIS_ROUTINE_NAME = trim(THIS_MOD_NAME) // ':create_tor_pos'
    
    ! --- Routine parameters
    type(t_tor_pos_list), target, intent(inout) :: pos_list
    real*8,  optional,            intent(in)    :: phi, phistart, phiend
    integer, optional,            intent(in)    :: iplane, nphi
    integer,                      intent(out)   :: ierr
    
    ! --- Local variables
    type(t_tor_pos), pointer :: pos
    integer :: i
    
    ierr = 0
    
    if ( present(phi) ) then
      
      call alloc_tor_pos(pos_list, 1)
      pos     => pos_list%pos(1)
      pos%phi = phi
      
    else if ( present(iplane) ) then
      
      call alloc_tor_pos(pos_list, 1)
      pos      => pos_list%pos(1)
      pos%phi = 2.d0 * PI * real(iplane) / real(n_plane * n_period) !###check###
      
    else if ( present(phistart) .and. present(phiend) .and. present(nphi) ) then
      
      call alloc_tor_pos(pos_list, nphi)
      if (nphi==1) then
        pos     => pos_list%pos(1)
        pos%phi = phistart 
      else
        do i = 1, nphi
          pos     => pos_list%pos(i)
          pos%phi = phistart + (phiend-phistart) * real(i-1)/real(nphi-1)
        end do
      endif
     
    else if ( present(nphi) ) then
      
      call alloc_tor_pos(pos_list, nphi)
      do i = 1, nphi
        pos     => pos_list%pos(i)
        pos%phi = (2.d0*PI)/real(n_period) * real(i-1)/real(nphi)
        ! IMPORTANT: Position phi=0 included and position phi=2pi/n_period ommitted!
      end do
      pos_list%full_period = .true.
      
    else
      
      ierr = 99
      write(*,*)
      write(*,*) 'ERROR in '//trim(THIS_ROUTINE_NAME)//':'
      write(*,*) 'No valid representation for toroidal position(s) provided.'
      write(*,*)
      
    end if
    
  end subroutine create_tor_pos
  
  
  
  
  
  !> Output positions to a file for debugging purposes.
  subroutine output_pos(filename, pol_pos_list, tor_pos_list)
    
    character(len=*),               intent(in) :: filename
    type(t_pol_pos_list), optional, intent(in) :: pol_pos_list
    type(t_tor_pos_list), optional, intent(in) :: tor_pos_list
    
    integer :: i, j, k
    
    open(37, file=trim(filename))
    
    if ( present(pol_pos_list) .and. present(tor_pos_list) ) then
      do j = 1, pol_pos_list%n_pos(2)
        do i = 1, pol_pos_list%n_pos(1)
          do k = 1, tor_pos_list%n_pos
            write(37,'(9es25.15)') pol_pos_list%pos(i,j)%R, pol_pos_list%pos(i,j)%Z,               &
              pol_pos_list%pos(i,j)%theta_star, pol_pos_list%pos(i,j)%r_minor,                     &
              pol_pos_list%pos(i,j)%length, tor_pos_list%pos(k)%phi
          end do
        end do
      end do
    else if ( present(pol_pos_list) ) then
      do j = 1, pol_pos_list%n_pos(2)
        do i = 1, pol_pos_list%n_pos(1)
            write(37,'(9es25.15)') pol_pos_list%pos(i,j)%R, pol_pos_list%pos(i,j)%Z,               &
              pol_pos_list%pos(i,j)%theta_star, pol_pos_list%pos(i,j)%r_minor,                     &
              pol_pos_list%pos(i,j)%length
        end do
      end do
    else if ( present(tor_pos_list) ) then
      do k = 1, tor_pos_list%n_pos
        write(37,'(9es25.15)') tor_pos_list%pos(k)%phi
      end do
    else
      write(*,*) 'WARNING: output_pos has no data to output...'
    end if
    
    close(37)
    
  end subroutine output_pos
  
  
  
  
  
end module mod_position

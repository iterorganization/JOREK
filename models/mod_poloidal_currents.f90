!> Contains functions to calculate poloidal currents and the associated FF' 
!!
!!  * Poloidal currents are calculated from teh equilibrium assumption JxB=\grad p 
module mod_poloidal_currents
  
  use constants
  use mod_parameters
  use data_structure
  use gauss
  use basis_at_gaussian
  use tr_module
  use phys_module
  use mod_interp
  use convert_character
  
  implicit none
  
  private
  
  public :: J_pol, normal_bnd_curr, integrated_normal_bnd_curr 
  
  contains
  
  !!-------------------------------------------------------------------
  !> Calculates poloidal currents from equilibrium assumption
  !!
  !!                  JxB = \grad p
  !! 
  !! which gives in the JOREK coordinate system and variables
  !!
  !!      J_R = (-j * B_R - R * dp/dZ) / F0
  !!      J_Z = (-j * B_Z + R * dp/dR) / F0
  !!
  !! Also the local FFprime value is given, which we define as   
  !!
  !!      FFp = F0 * F' = F0 * (J_pol \cdot B_pol) / B_pol^2  
  !!-------------------------------------------------------------------
  subroutine J_pol(node_list, element_list, i_elm, s, t, i_plane, axisym, JR, JZ, FFp)

    implicit none

    type (type_node_list),    intent(in) :: node_list
    type (type_element_list), intent(in) :: element_list
    
    integer, intent(in)       :: i_elm, i_plane ! element index / toroidal plane 
    real*8,  intent(in)       ::   s,  t        ! s-t local coordinates
    real*8,  intent(inout)    ::  JR, JZ, FFp   ! output current density and local FFprime
    logical, intent(in)       :: axisym         ! if true, only calculate axisymmetric component

    ! --- local variables    
    real*8     :: Psi,Ps_s, Ps_t, Ps_st, Ps_ss, Ps_tt
    real*8     :: ZJ ,ZJ_s, ZJ_t, ZJ_st, ZJ_ss, ZJ_tt
    real*8     :: RHO,RHO_s,RHO_t,RHO_st,RHO_ss,RHO_tt
    real*8     :: Ti0,Ti0_s,Ti0_t,Ti0_st,Ti0_ss,Ti0_tt
    real*8     :: Te0,Te0_s,Te0_t,Te0_st,Te0_ss,Te0_tt
    real*8     :: T0,T0_s,T0_t,T0_st,T0_ss,T0_tt
    real*8     :: R,R_s,R_t,R_st,R_ss,R_tt,Z,Z_s,Z_t,Z_st,Z_ss,Z_tt
    real*8     :: xjac, psi_x, psi_y, T0_x, T0_y, RHO_x, RHO_y
    real*8     :: BR, BZ, P0_R, P0_Z, zj_sum, rho_sum, T0_sum, grad_psi
    integer    :: i_tor
   
    call interp_RZ(node_list,element_list,i_elm,s,t,R,R_s,R_t,R_st,R_ss,R_tt,Z,Z_s,Z_t,Z_st,Z_ss,Z_tt)
    xjac  = R_s * Z_t - R_t * Z_s

    psi_x  = 0.d0;   psi_y = 0.d0;    T0_x = 0.d0;     T0_y = 0.d0 
    RHO_x  = 0.d0;   RHO_y = 0.d0;  zj_sum = 0.d0;  rho_sum = 0.d0
    T0_sum = 0.d0 

    do i_tor=1, n_tor

      if ( ( i_tor > 1 ) .and. axisym  ) cycle ! Just include the n=0 mode

      call interp(node_list,element_list,i_elm,1,i_tor,s,t,Psi,Ps_s, Ps_t, Ps_st, Ps_ss, Ps_tt)
      call interp(node_list,element_list,i_elm,3,i_tor,s,t,ZJ ,ZJ_s, ZJ_t, ZJ_st, ZJ_ss, ZJ_tt)
      if (jorek_model .eq. 400 .or. jorek_model .eq. 502) then
        call interp(node_list,element_list,i_elm,5,i_tor,s,t,RHO,RHO_s,RHO_t,RHO_st,RHO_ss,RHO_tt)
        call interp(node_list,element_list,i_elm,6,1,s,t,Ti0,Ti0_s,Ti0_t,Ti0_st,Ti0_ss,Ti0_tt)
        call interp(node_list,element_list,i_elm,n_var,1,s,t,Te0,Te0_s,Te0_t,Te0_st,Te0_ss,Te0_tt) !WARNING, here we always assume the electron temperature is the last variable
        T0    = Ti0   + Te0
        T0_s  = Ti0_s + Te0_s
        T0_t  = Ti0_t + Te0_t
      elseif (jorek_model < 200) then 
        T0 = 0.d0;  T0_s = 0.d0; T0_t = 0.d0; rho = 1.d0; rho_t = 0.d0; rho_s = 0.d0
      else
        call interp(node_list,element_list,i_elm,6,i_tor,s,t,T0,T0_s,T0_t,T0_st,T0_ss,T0_tt) 
        call interp(node_list,element_list,i_elm,5,i_tor,s,t,RHO,RHO_s,RHO_t,RHO_st,RHO_ss,RHO_tt)
      endif

      zj_sum  = zj_sum   +  ZJ * HZ(i_tor,i_plane)
      rho_sum = rho_sum  + RHO * HZ(i_tor,i_plane)
      T0_sum  = T0_sum   +  T0 * HZ(i_tor,i_plane)

      if (abs(xjac) > 1.d-6) then 

        psi_x = psi_x + (   Z_t * Ps_s - Z_s * Ps_t )   / xjac * HZ(i_tor,i_plane)
        psi_y = psi_y + ( - R_t * Ps_s + R_s * Ps_t )   / xjac * HZ(i_tor,i_plane)

        T0_x  = T0_x  + (   Z_t * T0_s - Z_s * T0_t )   / xjac * HZ(i_tor,i_plane)
        T0_y  = T0_y  + ( - R_t * T0_s + R_s * T0_t )   / xjac * HZ(i_tor,i_plane)
 
        RHO_x = RHO_x + (   Z_t * RHO_s - Z_s * RHO_t ) / xjac * HZ(i_tor,i_plane)
        RHO_y = RHO_y + ( - R_t * RHO_s + R_s * RHO_t ) / xjac * HZ(i_tor,i_plane)       

      endif
       
    enddo   ! n_tor loop     

    P0_R = T0_sum * RHO_x + T0_x * RHO_sum
    P0_Z = T0_sum * RHO_y + T0_y * RHO_sum
    BR   =  psi_y/R
    BZ   = -psi_x/R

    JR   = ( -ZJ_sum * BR - R * P0_Z ) / F0
    JZ   = ( -ZJ_sum * BZ + R * P0_R ) / F0

    grad_psi = sqrt(psi_x**2.d0 + psi_y**2.d0)
    if ( grad_psi > 1.d-6 ) then
      FFp    = ZJ_sum + (R**2.d0)*(psi_x*P0_R + psi_y*P0_Z)/(grad_psi**2.0d0)
    else
      FFp    = ZJ_sum !--- difficult to correct the other term when grad_psi=0...
    endif
    
  end subroutine J_pol 





  !!-------------------------------------------------------------------
  !> Calculates poloidal normal current integrated along the wall 
  !!
  !!    I_halo(\phi) = 1/2 * \int |J_n| R dl             
  !! 
  !!          I_halo = \int I_halo (\phi) d\phi
  !!
  !! See documentation in attached file in JIRA issue IMAS-2016
  !!-------------------------------------------------------------------
  subroutine integrated_normal_bnd_curr(node_list, bnd_node_list, bnd_elm_list, I_halo, TPF, print_halo)

    implicit none

    type (type_node_list),        intent(in)    :: node_list
    type (type_bnd_element_list), intent(in)    :: bnd_elm_list   
    type (type_bnd_node_list),    intent(in)    :: bnd_node_list
    real*8,                       intent(inout) :: I_halo, TPF      
    logical,                      intent(in)    :: print_halo      

    ! --- Local variables
    integer               :: m_bndelem, m_pt, m_elm, mv1, mp, in, ms
    integer               :: k_vertex, k_dof, k_node, k_dir, ierr, i_file
    real*8                :: zj(n_gauss, n_plane), psi_s(n_gauss, n_plane)
    real*8                :: rho(n_gauss, n_plane), T0(n_gauss, n_plane)
    real*8                :: rho_s(n_gauss, n_plane), T0_s(n_gauss, n_plane) 
    real*8                :: P0_s(n_gauss, n_plane)
    real*8                :: I_halo_mp(n_plane), phi(n_plane)
    real*8                :: k_size, I_net
    real*8                :: R(n_gauss), Z(n_gauss), R_s(n_gauss), Z_s(n_gauss)    
    type(type_node)       :: node_k
    type(type_bnd_element):: bndelem
    character(len=1024)   :: filename
    character(len=19), parameter :: DIR = './I_halo_phi_files/' 

    ! --- create folder with files for each time with I_halo(phi) profile if
    ! --- n_tor > 1 
    if ((n_tor > 1) .and. print_halo) then
      call system('mkdir -p '//DIR)

      write(filename,'(4a)') DIR, 'I_halo_phi_tnow_', trim(real2str(t_now,'(f12.4)')), '.dat'
      i_file=133
      open(i_file, file=trim(filename), form='formatted', status='replace', access='sequential',  &
          iostat=ierr)    
      write(i_file,'(a)') '#             phi               I_halo(phi) [MA/rad]'
    endif
 
    I_halo    = 0.d0
    I_net     = 0.d0
    I_halo_mp = 0.d0

    !--- go through the boundary elements
    do m_bndelem = 1, bnd_elm_list%n_bnd_elements
    
      bndelem = bnd_elm_list%bnd_element(m_bndelem)
  
      !--- calculate values at gaussian points on the element
      R     = 0.d0;      Z = 0.d0;   R_s = 0.d0;   Z_s = 0.d0      
      zj    = 0.d0;  psi_s = 0.d0;   rho = 0.d0;    T0 = 0.d0
      rho_s = 0.d0;   T0_s = 0.d0;  P0_s = 0.d0
   
      do k_vertex = 1, 2
        do k_dof = 1, 2
          k_node   = bndelem%vertex(k_vertex)
          k_dir    = bndelem%direction(k_vertex,k_dof)
          k_size   = bndelem%size(k_vertex,k_dof)
          node_k   = node_list%node(k_node)
        
          R  (:)   = R  (:)  + node_k%x(1,k_dir,1) * k_size * H1  (k_vertex,k_dof,:)
          Z  (:)   = Z  (:)  + node_k%x(1,k_dir,2) * k_size * H1  (k_vertex,k_dof,:)
          R_s(:)   = R_s(:)  + node_k%x(1,k_dir,1) * k_size * H1_s(k_vertex,k_dof,:)
          Z_s(:)   = Z_s(:)  + node_k%x(1,k_dir,2) * k_size * H1_s(k_vertex,k_dof,:)
        
          do mp=1,n_plane
            do in=1,n_tor
              psi_s(:,mp) = psi_s(:,mp)  + node_k%values(in,k_dir,1) * k_size*H1_s(k_vertex,k_dof,:) * HZ(in,mp)
              zj   (:,mp) = zj   (:,mp)  + node_k%values(in,k_dir,3) * k_size*  H1(k_vertex,k_dof,:) * HZ(in,mp)
              if (jorek_model > 190) then  
                rho  (:,mp) = rho  (:,mp)  + node_k%values(in,k_dir,5) * k_size*  H1(k_vertex,k_dof,:) * HZ(in,mp)
                rho_s(:,mp) = rho_s(:,mp)  + node_k%values(in,k_dir,5) * k_size*H1_s(k_vertex,k_dof,:) * HZ(in,mp)
                T0   (:,mp) = T0   (:,mp)  + node_k%values(in,k_dir,6) * k_size*  H1(k_vertex,k_dof,:) * HZ(in,mp)
                T0_s (:,mp) = T0_s (:,mp)  + node_k%values(in,k_dir,6) * k_size*H1_s(k_vertex,k_dof,:) * HZ(in,mp)
              else
                rho = 1.d0; rho_s = 0.d0; T0 = 0.d0; T0_s = 0.d0 
              endif
            enddo
          enddo
        
        end do
      end do

      do mp=1, n_plane
        P0_s(:,mp)     = T0(:,mp) * rho_s(:,mp) + T0_s(:,mp) * rho(:,mp)
      enddo
  
      !--- integrate wished function over the element
      do ms = 1, n_gauss
        do mp = 1, n_plane
          I_halo = I_halo + wgauss(ms) * abs( ( zj(ms,mp) * psi_s(ms,mp) + R(ms)**2.d0 * P0_s(ms,mp) )/F0 )  
          I_net  = I_net  + wgauss(ms) *      ( zj(ms,mp) * psi_s(ms,mp) + R(ms)**2.d0 * P0_s(ms,mp) )/F0  
   
          I_halo_mp(mp) = I_halo_mp(mp) + wgauss(ms) * abs( ( zj(ms,mp) * psi_s(ms,mp) + R(ms)**2.d0 * P0_s(ms,mp) )/F0 ) 
          phi(mp)       = 2.d0*PI*float(mp-1)/float(n_plane) / float(n_period)
        enddo  
      enddo
  
    enddo !--- bnd elements

    I_halo       = 0.5*I_halo       / mu_zero * 1.d-6 * 2.d0 * PI / n_plane
    I_net        = 0.5*I_net        / mu_zero * 1.d-6 * 2.d0 * PI / n_plane
    I_halo_mp(:) = 0.5*I_halo_mp(:) / mu_zero * 1.d-6 * 2.d0 * PI 
    TPF          = maxval(I_halo_mp)/I_halo 

    if (print_halo) then
      write(*,*) ' I_halo, I_net [MA]:  ', I_halo, I_net
      write(*,*) '               TPF :  ', TPF
    endif

    if ((n_tor > 1) .and. print_halo) then 
      do mp = 1, n_plane 
        write(i_file,'(2es20.10)') phi(mp), I_halo_mp(mp)
      enddo
      close(i_file) 
    endif

  end subroutine integrated_normal_bnd_curr





  !!-------------------------------------------------------------------
  !> Calculates normal current density to the JOREK boundary
  !! and gives it as a function of R,Z in a file 
  !! See documentation in attached file in JIRA issue IMAS-2016
  !!-------------------------------------------------------------------
  subroutine normal_bnd_curr(node_list, element_list, bnd_node_list, bnd_elm_list, i_plane, si_units)

    use mod_basisfunctions

    implicit none

    integer,                      intent(in)    :: i_plane
    logical,                      intent(in)    :: si_units 
    type (type_node_list),        intent(in)    :: node_list
    type (type_element_list),     intent(in)    :: element_list
    type (type_bnd_element_list), intent(in)    :: bnd_elm_list   
    type (type_bnd_node_list),    intent(in)    :: bnd_node_list

    ! --- Local variables
    integer               :: m_bndelem, m_pt, m_elm, mv1, in, ms, i, j
    integer               :: k_vertex, k_dof, k_node, k_dir, i_file, ierr
    real*8                :: zj(n_gauss), psi_s(n_gauss),rho(n_gauss), T0(n_gauss)
    real*8                :: rho_s(n_gauss), T0_s(n_gauss) 
    real*8                :: P0_s(n_gauss), J_normal(n_gauss)
    real*8                :: k_size, fact, R_c, Z_c, vec_inside(2), grad_t(2)
    real*8                :: G(4,4), sign_out
    real*8                :: R(n_gauss), Z(n_gauss), R_s(n_gauss), Z_s(n_gauss)    
    type(type_node)       :: node_k
    type(type_element)    :: elm_k
    type(type_bnd_element):: bndelem
    character(len=1024)   :: filename
    character(len=14), parameter :: DIR = './jnorm_files/' 

    fact = 1.d0
    if (si_units) fact = 1.d0/mu_zero

    call system('mkdir -p '//DIR)
    
    write(filename,'(4a)') DIR, 'Jnorm_RZ_tnow_', trim(real2str(t_now,'(f12.4)')), '.dat'
    i_file=133

    open(i_file, file=trim(filename), form='formatted', status='replace', access='sequential',  &
        iostat=ierr)
    
    write(i_file,'(a)') '#               R               Z              J_norm'
 
    !--- go through the boundary elements
    do m_bndelem = 1, bnd_elm_list%n_bnd_elements
    
      bndelem = bnd_elm_list%bnd_element(m_bndelem)
      elm_k   = element_list%element(bndelem%element)
 
      !--- calculate values at gaussian points on the element
      R     = 0.d0;      Z = 0.d0;   R_s = 0.d0;   Z_s = 0.d0      
      zj    = 0.d0;  psi_s = 0.d0;   rho = 0.d0;    T0 = 0.d0
      rho_s = 0.d0;   T0_s = 0.d0;  P0_s = 0.d0
   
      do k_vertex = 1, 2
        do k_dof = 1, 2
          k_node   = bndelem%vertex(k_vertex)
          k_dir    = bndelem%direction(k_vertex,k_dof)
          k_size   = bndelem%size(k_vertex,k_dof)
          node_k   = node_list%node(k_node)
        
          R  (:)   = R  (:)  + node_k%x(1,k_dir,1) * k_size * H1  (k_vertex,k_dof,:)
          Z  (:)   = Z  (:)  + node_k%x(1,k_dir,2) * k_size * H1  (k_vertex,k_dof,:)
          R_s(:)   = R_s(:)  + node_k%x(1,k_dir,1) * k_size * H1_s(k_vertex,k_dof,:)
          Z_s(:)   = Z_s(:)  + node_k%x(1,k_dir,2) * k_size * H1_s(k_vertex,k_dof,:)
        
          do in=1,n_tor
            psi_s(:) = psi_s(:) + node_k%values(in,k_dir,1) * k_size*H1_s(k_vertex,k_dof,:) * HZ(in,i_plane)
            zj   (:) = zj   (:) + node_k%values(in,k_dir,3) * k_size*  H1(k_vertex,k_dof,:) * HZ(in,i_plane)
            if (jorek_model > 200 ) then
              rho  (:) = rho  (:) + node_k%values(in,k_dir,5) * k_size*  H1(k_vertex,k_dof,:) * HZ(in,i_plane)
              rho_s(:) = rho_s(:) + node_k%values(in,k_dir,5) * k_size*H1_s(k_vertex,k_dof,:) * HZ(in,i_plane)
              T0   (:) = T0   (:) + node_k%values(in,k_dir,6) * k_size*  H1(k_vertex,k_dof,:) * HZ(in,i_plane)
              T0_s (:) = T0_s (:) + node_k%values(in,k_dir,6) * k_size*H1_s(k_vertex,k_dof,:) * HZ(in,i_plane)
            else
              rho = 1.d0; rho_s = 0.d0; T0 = 0.d0; T0_s = 0.d0 
            endif        
          enddo
        
        end do
      end do


      !--- Find out correct sign of the normal (it has to point outwards the domain)
      !---------------------------------------------------------------------------------- 
      ! --- Calculate an inside point on the element to calculate the
      ! direction of bnd normals
      call basisfunctions(xgauss(2),xgauss(2), G)  
      R_c = 0.d0 ;  Z_c = 0.d0 
      do i = 1, n_vertex_max
        do j = 1, n_order+1
          node_k = node_list%node(elm_k%vertex(i)) 
          R_c    = R_c + node_k%x(1,j,1) * elm_k%size(i,j) * G(i,j)
          Z_c    = Z_c + node_k%x(1,j,2) * elm_k%size(i,j) * G(i,j)
        enddo
      enddo  
      vec_inside = (/ R_c - R(2), Z_c - Z(2) /)       ! vector pointing towards the domain
      grad_t     = (/ -Z_s(2) , R_s(2) /)     ! gradient of the coordinate t (normal to the boundary here)
      sign_out   = -1.d0 * sign( 1.d0, ( vec_inside(1)*grad_t(1) + vec_inside(2)*grad_t(2) ) )  
      !--------------------------------------------------------------------------------

      P0_s(:)     = T0(:) * rho_s(:) + T0_s(:) * rho(:)
      J_normal(:) = (zj(:) * psi_s(:) + R(:)**2.d0 * P0_s(:)) &
                    /sqrt(R_s(:)**2.d0 + Z_s(:)**2.d0) / (F0*R(:)) * sign_out

      do ms = 1, n_gauss 
        write(i_file,'(3es20.10)') R(ms), Z(ms), J_normal(ms) * fact
      enddo
  
    enddo !--- bnd elements
    close(i_file)

  end subroutine normal_bnd_curr


  
end module mod_poloidal_currents

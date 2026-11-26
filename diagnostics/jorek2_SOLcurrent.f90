! Calculates SOL currents from the jorek_restart.h5 file of the simulation folder
! the program is run in.
! Solves for Stangeby eq 17.29 by determining the fluxline integrals by line tracing
! (based on jorek2_poincare) and solving the equation using a Newton solver
! 
! requires the input namelist to be piped into the program:
! ./jorek2_SOLcurrent < input
!
! References:
! 
! Stangeby, P.C. (2000). The Plasma Boundary of Magnetic Fusion Devices (1st ed.). 
! CRC Press. https://doi.org/10.1201/9780367801489
! 
! G.M. Staebler and F.L. Hinton 1989 Nucl. Fusion 29 1820 
! https://doi.org/10.1088/0029-5515/29/10/017
!
! F.L. Hinton in chapter 1.5 of Basic Plasma Physics (eds; A. Galeev and R.N. Sudan),
! North-Holland, Amsterdam (1983) ISBN: 9780444866455
!
! *** warning *** currently does not work with full MHD or two temperature model
!
! *** warning *** please check whether for a given fluxline the calculated current density 
! starting at the inner target is the same as that starting at the outer target. It should, but
! for some cases you need to increase n_phi_steps to be able to resolve strong fluctuations
! in the parallel pressure gradient
program jorek2_SOLcurrent

  use constants
  use data_structure
  use phys_module
  use basis_at_gaussian
  use elements_nodes_neighbours
  use mod_neighbours
  use mod_import_restart
  use mod_log_params
  use equil_info, only : get_psi_n, ES
  use mod_interp
  use mod_boundary
  use mod_position
  use mod_model_settings
  use corr_neg
  use mod_plasma_functions, only: coulomb_log_ei
  !$ use omp_lib
  
  implicit none

  type (type_bnd_element_list), pointer :: bnd_elm_list    
  type (type_bnd_node_list),    pointer :: bnd_node_list    

  type(t_pol_pos_list), target :: initial_bnd_pos

  character(len=512) :: s
  real*8, allocatable :: Phi_start_list(:)
  real*8  :: R_b, R_e, Z_b, Z_e, Phi_b, Phi_e  !< start and end points of the traced fieldline
  real*8  :: Te_b, Te_e, ne_b, ne_e            !< quantities at the beginning and end
  real*8  :: Te_h, ratio_T, ne_h, ratio_n      !< hot end of fluxline quantities, and ratios
  real*8  :: c_sh                              !< ion speed of sound at the hot end of the fluxline
  logical :: hot2cold                          !< true if Te_b > Te_e, else false. Tracks whether the fluxline starts at the hot or cold side, to determine the sign of j//
  real*8  :: r0_corr, ne, ne_s, ne_t, ne_phi   !< ne and Te and their derivatives at the tracer
  real*8  :: T_corr, Te, Te_s, Te_t, Te_phi           
  real*8  :: n_norm, T_norm                    !< normalisation constants
  real*8  :: nu_ei                             !< electron ion colission frequency, Stangeby 19.23
  real*8  :: lnA                               !< Coulomb logarithm (ignoring impurities)
  real*8, parameter  :: lambda11=1.975d0       !< Spitzer-Härm coefficient, see Hinton table 1
  real*8  :: resisty                           !< running fluxline average of Stangeby 9.35
  real*8  :: av_sigma                          !< average conductivity as calculated from Staebler eq 12
  real*8  :: L                                 !< connection length (m)
  real*8  :: pe_s, pe_t, pe_R, pe_Z, pe_phi    !< pressure and their derivatives at the tracer
  real*8  :: B, B_s, B_t, B_r, B_z, B_phi      !< magnetic field for \grad_\parallel Pe equation
  real*8  :: par_grad_p                        !< parallel pressure gradient at the tracer
  real*8  :: pressure_drop                     !< Sanity check: integral of parallel pe gradient, i.e. electron pressure drop from hot to cold side (Pa)
  real*8  :: pe_integral                       !< pressure integral \int (1/n) dp_e/ds_\parallel ds_\parallel
  real*8  :: pe_term                           !< pressure contribution to j_hat_par (1/(T_h[eV] e)) \int (1/n) dp_e/ds_\parallel ds_\parallel
  real*8  :: gamma_fact                        !< \gamma factor as given by Stangeby eq 17.24
  real*8  :: alpha                             !< constant factor defined in Stangeby eq 17.17
  real*8  :: j_hat_par                         !< normalised parallel current density along the fluxline, direction from hold to cold target
  real*8  :: j_par                             !< j_\parallel [A/m²], current density along fluxline, direction from hot to cold target
  real*8  :: j_wall                            !< j_{target} [A/m²], current density perpendicular to the target, direction from hot to cold target
  real*8  :: I_loc                             !< [A], current onto the wall corresponding to the area represented by this starting position and - angle
  real*8  :: int_I_wall                        !< \int I over the full wall [A], should be 0 (sanity check)
  logical, parameter :: correct4fluxexpansion = .true. !< whether to correct the j_wall for the fact that j_par is not actually constant due to flux expansion. Only takes major radius into account
  real*8  :: R_average                         !< [m] average major radius of the beginning and end of the traced fieldline (i.e. R_average=(R_b+R_e)/2), used for the flux expansion correction
  integer :: my_id
  integer :: i, j, iside_i, iside_j, ip, np, i_tor, i_harm
  integer :: i_elm, ifail, i_elm_out, i_elm_prev, i_elm_tmp,i_steps, i_phi0
  real*8  :: R_line, Z_line, s_line, t_line, p_line, R_mid, Z_mid, s_mid, t_mid, p_mid
  real*8  :: R, R_s, R_t, Z, Z_s, Z_t 
  real*8,dimension(3) :: P, P_s, P_t, P_phi
  real*8  :: tol, delta_phi, delta_phi_signed, Zjac, psi_s, psi_t, R_in, Z_in, R_out, Z_out, Rmin, Rmax, Zmin, Zmax, delta_s, delta_t, R_keep, Z_keep
  real*8  :: small_delta, small_delta_s, small_delta_t, delta_phi_local
  real*8  :: atmp, cur_pert
  real*8  :: psi_out
  integer :: ierr
  real*8  :: dr,R_old,Z_old   !< integration variables
  logical :: in_domain        !< keeps track of whether tracer position is inside the domain or on its boundary
  real*8  :: B3(3)            !< Magnetic field [T]
  real*8  :: bdotn            !< b (=B/|B|) dot outward pointing normal to the boundary
  real*8  :: inv_st_jac, psi_R, psi_Z
  real*8  :: B_deb(3), grad_p_deb(3)
  character(len=6)   :: t_index_char

  logical, parameter :: debug = .false.            !< useful for debugging the main code 
  logical, parameter :: write_debug_file = .false. !< useful for debugging the Newton solver
  logical, parameter :: solve_j = .true. !< If true, Stangeby eq 17.29 is solved in the routine by a custom Newton solver

  integer, parameter :: n_phi_steps   = 150000 !number of steps in the toroidal direction
  integer, parameter :: n_phi0        = 1!2    !number of toroidal angles to start from for each poloidal bnd point
  integer, parameter :: n_elm_pts     = 1!10   !number of starting points along each bnd element

  character(len=512) :: dir = "./SOLcurrent/" !< subdirectory you want the datafile for this step to write to

  write(*,*) '***************************************'
  write(*,*) '* JOREK2_SOLcurrent                   *'
  write(*,*) '***************************************'

  allocate(bnd_elm_list)
  allocate(bnd_node_list)

  my_id=0

  ! --- Initialize mode and mode_type arrays
  call det_modes()

  call initialise_parameters(my_id,  "__NO_FILENAME__")
  
  call import_restart(node_list, element_list, 'jorek_restart', rst_format, ierr)
  
  call initialise_basis

  call boundary_from_grid(node_list, element_list, bnd_node_list, bnd_elm_list, .false.)

  write(*,*) 'Initialising element neighbours'
  allocate(element_neighbours(4,element_list%n_elements))

  element_neighbours = 0

  ! --- Neighbour search is O(N) with a bounding-box grid (see build_element_neighbours in mod_neighbours)
  call build_element_neighbours(node_list, element_list, element_neighbours)

  ! maximum stepsize delta_phi (stepsize is automatically made smaller if an element boundary is crossed)
  delta_phi = 2.d0 * PI / float(n_period*n_phi_steps)

  allocate(Phi_start_list(n_phi0))
  do i_phi0 = 1, n_phi0
    Phi_start_list(i_phi0) = float(i_phi0-1)/float(n_phi0) * 2.d0 * PI / float(n_period)
  end do

  initial_bnd_pos = bnd_pos(node_list, element_list, bnd_node_list, bnd_elm_list, n_elm_pts)

  np = initial_bnd_pos%n_pos(2)

  write(*,*)
  write(*,*) 'minimal number of steps in toroidal turn:          ',n_phi_steps
  write(*,*) 'number of poloidal starting positions per element: ',n_elm_pts
  write(*,*) 'toroidal angles per starting position:             ',n_phi0
  write(*,*) 'number of boundary elements:                       ',bnd_elm_list%n_bnd_elements
  write(*,*) 'number of fieldlines to trace:                     ',n_phi0*np
  write(*,*)

  mode(1) = 0
  do i=1,(n_tor-1)/2
    mode(2*i)   = i * n_period
    mode(2*i+1) = i * n_period
  enddo
  write(*,*) ' modes   : ',mode
  write(*,*) ' nperiod : ',n_period
  Write(*,*)
  
  if(.not. solve_j) then
    write(*,*) 'j_hat_par will not be solved locally'
    write(*,*)
  end if

  !> Stangeby eq 17.17 \alpha = 1/2 (m_i/(π m_e))^0.5, only depends on ion mass
  alpha = 0.5d0*sqrt(central_mass * ATOMIC_MASS_UNIT/(PI * MASS_ELECTRON)) 

  ! --- Open the output file to which the SOL current data will be written in ascii format
  call system('mkdir -p '//trim(DIR))

  write(t_index_char,'(I6.6)') index_start
  open(21,file=trim(DIR)//'step'//t_index_char//'.dat')
  if(write_debug_file) open(22,file='debug_Newton_solver.dat')
  if(debug) open(23,file='debug2.dat')

  write(21,'(A17, 19A18)') '#            R_b','Z_b','Phi_b','R_e','Z_e','Phi_e','ne_b','ne_e','Te_b','Te_e','L','av_sigma_par','P_h - P_c','pe_term','j_hat_par','j_par','bdotn','j_wall','dl','I_loc'
  
  !normalisation constants to be used
  n_norm = central_density*1.d20                     ! ne[atoms/m³]     = n[jor] * n_norm
  T_norm = 1.d0 / (2.d0 * EL_CHG * MU_ZERO * n_norm) ! Te[eV] = T[eV]/2 = T[jor] * T_norm

  int_I_wall = 0
  
  ! --- Trace the fieldlines
  !$omp parallel do schedule(monotonic:dynamic,1) default(none)                                                                                                  &
  !$omp shared(np, initial_bnd_pos, Phi_start_list, node_list, element_list, n_norm, T_norm, F0, delta_phi, element_neighbours, alpha, central_mass)             &
  !$omp private(i_phi0, s_line, t_line, i_elm, i_elm_prev, Phi_b, p_line, P, P_s, P_t, P_phi, R, R_s, R_t, Z, Z_s, Z_t, ne_b, Te_b, R_b, Z_b, R_old, Z_old,      &
  !$omp         inv_st_jac,psi_R,psi_Z,B3,bdotn,delta_phi_signed,L,resisty,pe_integral,pressure_drop,i_steps,in_domain,s_mid, t_mid, p_mid, delta_s, delta_t,    &
  !$omp         small_delta_s, small_delta_t, small_delta, R_in, Z_in, i_elm_tmp, R_out, Z_out, delta_phi_local, r0_corr, ne, ne_s, ne_t, ne_phi, T_corr, Te,    &
  !$omp         Te_s, Te_t, Te_phi, Zjac, dr, lnA, nu_ei, pe_s, pe_t, pe_phi, pe_R, pe_Z, B_R, B_Z, B_phi, B, par_grad_p, R_e, Z_e, ne_e, Te_e, Phi_e, av_sigma, &
  !$omp         hot2cold, Te_h, ne_h, ratio_T, ratio_n, pe_term, c_sh, gamma_fact, j_hat_par, j_par, j_wall, R_average, I_loc)                                   &
  !$omp reduction(+:int_I_wall)
  do ip=1,np !loop over the starting positions
    
    do i_phi0 = 1, n_phi0 !loop over starting angles
      if(write_debug_file) then
        write(22,*)
        write(22,'(A80)') '-------------------------------------------------------------------------------'
        write(22,'(A3,I5,A7,I3,A1,I3)') 'ip=',ip,' angle ',i_phi0,'/',n_phi0
      end if
      
      s_line = initial_bnd_pos%pos(1,ip)%s
      t_line = initial_bnd_pos%pos(1,ip)%t
      i_elm  = initial_bnd_pos%pos(1,ip)%ielm

      if (.not. ((s_line .eq. 0.d0) .or. (s_line .eq. 1.d0) .or. (t_line .eq. 0.d0) .or. (t_line .eq. 1.d0))) &
        write(*,*) 'WARNING, startpoint is not edge',s_line,t_line,i_elm,i_elm_prev

      Phi_b = Phi_start_list(i_phi0)
      p_line = Phi_b
      
      ! Calculate the starting location and variables
      call interp_PRZ(node_list,element_list,i_elm,[var_psi,var_rho,var_T],3,s_line,t_line,p_line,P,P_s,P_t,P_phi,R,R_s,R_t,Z,Z_s,Z_t)
      ne_b = P(2) * n_norm
      ne_b = max(ne_b, 1.d16)
      Te_b = P(3) * T_norm
      Te_b = max(Te_b, 0.1d0)
      
      R_b   = R
      Z_b   = Z
      R_old = R
      Z_old = Z

      if(debug) write(*,*) 'R,Z',R,Z
      
      ! --- determine sign of delta_phi
      inv_st_jac = 1.d0/(R_s * Z_t - R_t * Z_s)
      psi_R      = (  P_s(1) * Z_t - P_t(1) * Z_s ) * inv_st_jac
      psi_Z      = (- P_s(1) * R_t + P_t(1) * R_s ) * inv_st_jac

      B3         = [+psi_Z, -psi_R, F0] * 1.d0/R
      !> since the phi element of the wall normal vector is 0, the dot product is that of the first two elements 
      bdotn      = dot_product(B3(1:2),initial_bnd_pos%pos(1,ip)%bnd_normal)/norm2(B3)
      delta_phi_signed = sign(abs(delta_phi),-bdotn*F0/R)
      if(debug) write(*,*) 'B:',B3, 'bdotn:',bdotn, 'delta_phi_signed:',delta_phi_signed
      
      !resetting running integrals
      L             = 0.d0
      resisty       = 0.d0
      pe_integral   = 0.d0
      pressure_drop = 0.d0

      i_steps = 0

      in_domain = .true.
      do while (in_domain)
      
        i_steps = i_steps + 1
        if(i_steps .ge. 1000*n_phi_steps) then
          !$omp critical
          write(*,*) 'ERROR: could not find connecting wall', R_b,Z_b,Phi_b,R,Z,p_line
          !$omp end critical
          cycle
        end if
      
        
        call step(i_elm,s_line,t_line,p_line,delta_phi_signed,delta_s,delta_t)

        s_mid = s_line + 0.5d0 * delta_s
        t_mid = t_line + 0.5d0 * delta_t
        p_mid = p_line + 0.5d0 * delta_phi_signed
            
        call step(i_elm,s_mid,t_mid,p_mid,delta_phi_signed,delta_s,delta_t)
        
        small_delta_s = 1.d0
      
        if  (s_line + delta_s .gt. 1.d0) then
      
          small_delta_s = (1.d0 - s_line)/delta_s

        elseif  (s_line + delta_s .lt. 0.d0) then

          small_delta_s = abs(s_line/delta_s)

        endif
      
        small_delta_t = 1.d0

        if  (t_line + delta_t .gt. 1.d0)  then
      
          small_delta_t = (1.d0 - t_line)/delta_t

        elseif  (t_line + delta_t .lt. 0.d0)  then

          small_delta_t = abs(t_line/delta_t)
      
        endif      

        small_delta = min(small_delta_s, small_delta_t)

        if (small_delta .lt. 1.d0)  then ! moving the tracked position by delta_t and delta_s moves it out of the element so correct that the tracked position ends up on the element instead
          
          ! This extra step in jorek2_poincare seems unnecessary and could potentially even lead to problems,
          ! so better to comment it out.
          !s_mid = s_line + 0.5d0 * small_delta * delta_s
          !t_mid = t_line + 0.5d0 * small_delta * delta_t
          !p_mid = p_line + 0.5d0 * small_delta * delta_phi_signed
            
          !call step(i_elm,s_mid,t_mid,p_mid,delta_phi_signed,delta_s,delta_t) ! why is this in here? better accuracy at the cost of getting a possibly getting new position s,t \in [0,1], meaning the position is not updated as it's inside the wrong if statement
      
          if (small_delta_s .lt. small_delta_t) then

            if (s_line + delta_s .gt. 1.d0) then
      
              s_line = 1.d0
              t_line = t_line + small_delta * delta_t
              p_line = p_line + small_delta * delta_phi_signed
        
              call interp_RZ(node_list,element_list,i_elm,s_line,t_line,R_in,Z_in)

              i_elm_prev = i_elm      
              i_elm      = element_neighbours(2,i_elm_prev)
              if ( i_elm == 0 ) then !tracer is at domain boundary
                if(debug) write(*,*) 'end of field line tracing (0)'
                in_domain = .false.
              else
                i_elm_tmp  = element_neighbours(4,i_elm)
          
                if (i_elm_prev .ne. i_elm_tmp) write(*,*) ' WARNING : CHANGE OF ORIENTATION (1)'

                s_line = 0.d0

                call interp_RZ(node_list,element_list,i_elm,s_line,t_line,R_out,Z_out)
          
                if ( (abs(R_in - R_out) .gt. 1.d-8) .or. (abs(Z_in - Z_out) .gt. 1.d-8)) &
                  write(*,'(A,2i6,4f8.4)') ' error in element change (1) ',i_elm_prev,i_elm,R_in,R_out,Z_in,Z_out
              endif
            elseif (s_line + delta_s .lt. 0.d0) then

              s_line = 0.d0
              t_line = t_line + small_delta * delta_t
              p_line = p_line + small_delta * delta_phi_signed

              call interp_RZ(node_list,element_list,i_elm,s_line,t_line,R_in,Z_in)
      
              i_elm_prev = i_elm      
              i_elm      = element_neighbours(4,i_elm_prev)
              if ( i_elm == 0 ) then ! at domain boundary
                if(debug) write(*,*) 'end of field line tracing (1)'
                in_domain = .false.
              else
                i_elm_tmp  = element_neighbours(2,i_elm)
          
                if (i_elm_prev .ne. i_elm_tmp) write(*,*) ' WARNING : CHANGE OF ORIENTATION (2)'

                s_line = 1.d0
                
                call interp_RZ(node_list,element_list,i_elm,s_line,t_line,R_out,Z_out)
          
                if ( (abs(R_in - R_out) .gt. 1.d-8) .or. (abs(Z_in - Z_out) .gt. 1.d-8)) &
                  write(*,'(A,2i6,4f8.4)') ' error in element change (2) ',i_elm_prev,i_elm,R_in,R_out,Z_in,Z_out
              end if

            endif ! after step, is s < 0 or s > 1?

          else ! step crosses t boundary first

            if (t_line + delta_t .gt. 1.d0) then
      
              s_line = s_line + small_delta * delta_s
              t_line = 1.d0
              p_line = p_line + small_delta * delta_phi_signed

              call interp_RZ(node_list,element_list,i_elm,s_line,t_line,R_in,Z_in)

              i_elm_prev = i_elm      
              i_elm      = element_neighbours(3,i_elm_prev)
              if ( i_elm == 0 ) then
                if(debug) write(*,*) 'end of field line tracing (2)'
                in_domain = .false.
              else
                i_elm_tmp  = element_neighbours(1,i_elm)
          
                if (i_elm_prev .ne. i_elm_tmp) write(*,*) ' WARNING : CHANGE OF ORIENTATION (3)'

                t_line = 0.d0

                call interp_RZ(node_list,element_list,i_elm,s_line,t_line,R_out,Z_out)
          
                if ( (abs(R_in - R_out) .gt. 1.d-8) .or. (abs(Z_in - Z_out) .gt. 1.d-8)) &
                  write(*,'(A,2i6,4f8.4)') ' error in element change (3) ',i_elm_prev,i_elm,R_in,R_out,Z_in,Z_out
              end if

            elseif (t_line + delta_t .lt. 0.d0) then

              s_line = s_line + small_delta * delta_s	
              t_line = 0.d0
              p_line = p_line + small_delta * delta_phi_signed

              call interp_RZ(node_list,element_list,i_elm,s_line,t_line,R_in,Z_in)

              i_elm_prev = i_elm      
              i_elm      = element_neighbours(1,i_elm_prev)
              if ( i_elm == 0 ) then
                if(debug) write(*,*) 'end of field line tracing (3)'
                in_domain = .false.
              else
                i_elm_tmp  = element_neighbours(3,i_elm)
          
                if (i_elm_prev .ne. i_elm_tmp) write(*,*) ' WARNING : CHANGE OF ORIENTATION (4)',i_elm_prev,i_elm

                t_line = 1.d0

                call interp_RZ(node_list,element_list,i_elm,s_line,t_line,R_out,Z_out)
          
                if ( (abs(R_in - R_out) .gt. 1.d-8) .or. (abs(Z_in - Z_out) .gt. 1.d-8)) &
                  write(*,'(A,2i6,4f8.4)') ' error in element change (4) ',i_elm_prev,i_elm,R_in,R_out,Z_in,Z_out
              end if

            endif ! after step, is t > 1 or t < 0?
        
          endif ! does step cross s or t boundary first?
          
        else ! step does not cross an element boundary

          s_line = s_line + delta_s
          t_line = t_line + delta_t
          p_line = p_line + delta_phi_signed

          small_delta = 1.d0
      
        endif ! does step cross an element boundary?
        
        delta_phi_local = small_delta * delta_phi_signed
        
        if(i_elm .le. 0) then
          call interp_PRZ(node_list,element_list,i_elm_prev,[var_psi,var_rho,var_T],3,s_line,t_line,p_line,P,P_s,P_t,P_phi,R,R_s,R_t,Z,Z_s,Z_t)
        else
          call interp_PRZ(node_list,element_list,i_elm,[var_psi,var_rho,var_T],3,s_line,t_line,p_line,P,P_s,P_t,P_phi,R,R_s,R_t,Z,Z_s,Z_t)
        end if
        
        r0_corr = corr_neg_dens(P(2))
        ne      = r0_corr  * n_norm
        ne_s    = P_s(2)   * n_norm
        ne_t    = P_t(2)   * n_norm
        ne_phi  = P_phi(2) * n_norm
        T_corr  = corr_neg_temp(P(3))
        Te      = T_corr   * T_norm
        Te_s    = P_s(3)   * T_norm
        Te_t    = P_t(3)   * T_norm
        Te_phi  = P_phi(3) * T_norm
        Zjac    = R_s * Z_t - R_t * Z_s

        dr = sqrt((R-R_old)**2 + (Z-Z_old)**2 + (R*delta_phi_local)**2)
        
        call coulomb_log_ei(P(3),T_corr,P(2),r0_corr, 0.d0, 0.d0, 0.d0, lnA)
        nu_ei = 0.51 * (EL_CHG**4.d0) * lnA * ne / (3.d0*sqrt(MASS_ELECTRON)*(EPS_ZERO**2.d0)*((2.d0*PI*EL_CHG*Te)**(3.d0/2.d0))) ! Stangeby eq 9.23
        
        pe_s   = ne_s   * Te * EL_CHG + ne * Te_s   * EL_CHG !< product rule for pe = ne Te[eV] e
        pe_t   = ne_t   * Te * EL_CHG + ne * Te_t   * EL_CHG
        pe_phi = ne_phi * Te * EL_CHG + ne * Te_phi * EL_CHG
        pe_R   =   ( pe_s   * Z_t - pe_t   * Z_s ) / Zjac
        pe_Z   = - ( pe_s   * R_t - pe_t   * R_s ) / Zjac
        psi_R  =   ( P_s(1) * Z_t - P_t(1) * Z_s ) / Zjac
        psi_Z  = - ( P_s(1) * R_t - P_t(1) * R_s ) / Zjac
        B_R    =   psi_Z / R
        B_Z    = - psi_R / R
        B_phi  =   F0    / R

        !cyclindrical coordinates have an orthonormal basis, 
        !so 3D pythagoras can be used for the length of a vector
        B = sqrt(B_R**2 + B_Z**2 +B_phi**2) 
        
        ! \nabla_parallel P_e = \hat{b} \cdot \nabla P_e = b_R P_R + b_Z P_Z + b_phi P_phi
        par_grad_p = (B_R/B)*pe_R + (B_Z/B)*pe_Z + (B_phi/B)*(1.d0/R)*pe_phi 
        
        !the following statement never wrote so we're good on the formulation of the parallel gradient
        !B_deb      = [B_R, B_Z, B_phi]
        !grad_p_deb = [pe_R,pe_Z,pe_phi]
        !if((abs(dot_product(B_deb/norm2(B_deb),grad_p_deb) - par_grad_p) .gt. 1.d-10) .or. (abs(norm2(B_deb) - B) .gt. 1.d-10 )) &
        !  write(*,*) 'PROBLEM! ', dot_product(B_deb/norm2(B_deb),grad_p_deb), par_grad_p, norm2(B_deb), B 
        !if((ip .eq. 344) .and. (i_steps .lt. 1000) ) write(23,'(A10,I8,10es18.8)') 'debug: ',i_steps,pe_s, pe_t, pe_phi, pe_R, pe_Z, B_R, B_Z, B_phi, B, par_grad_p

        !the direction of the integral is either the direction of B or -B
        !dp_e/ds_\parallel = par_grad_p if \vec{s} = \vec{b}, but dp_e/ds_\parallel = - par_grad_p if \vec{s} = - \vec{b}
        par_grad_p = sign(abs(par_grad_p), B_phi * delta_phi_local * par_grad_p)
        
        L             = L             +  1.d0                 * dr
        resisty       = resisty       + (nu_ei / ne )         * dr !1/(n_e \tau_ei) = \nu_ei / n_e
        pe_integral   = pe_integral   + (1 / ne) * par_grad_p * dr
        pressure_drop = pressure_drop + par_grad_p            * dr !useful sanity check for debugging
        R_old = R
        Z_old = Z

        if ( i_elm == 0 ) then
          if(debug) write(*,*) '(4), do nothing'
        end if
      
      enddo ! traced a fieldline for one bnd point for 1 starting angle phi

      if (.not. ((s_line .eq. 0.d0) .or. (s_line .eq. 1.d0) .or. (t_line .eq. 0.d0) .or. (t_line .eq. 1.d0))) &
        write(*,*) 'WARNING, endpoint   is not edge',s_line,t_line,i_elm,i_elm_prev

      call interp_PRZ(node_list,element_list,i_elm_prev,[var_psi,var_rho,var_T],3,s_line,t_line,p_line,P,P_s,P_t,P_phi,R,R_s,R_t,Z,Z_s,Z_t)
      R_e = R
      Z_e = Z
      
      ne_e = P(2) * n_norm
      ne_e = max(ne_e, 1.d16)
      Te_e = P(3) * T_norm
      Te_e = max(Te_e, 0.1d0)
      Phi_e = p_line
      
      if(resisty .ne. 0.d0) then  
        av_sigma = L * (lambda11 * EL_CHG**2 / MASS_ELECTRON) / resisty
      else
        av_sigma = 0.d0
      endif
      
      if (Te_b > Te_e) then
        hot2cold = .true.
        Te_h     = Te_b
        ne_h     = ne_b
        ratio_T  = Te_b/Te_e
        ratio_n  = ne_b/ne_e
      else
        hot2cold = .false.
        Te_h     = Te_e
        ne_h     = ne_e
        ratio_T  = Te_e/Te_b
        ratio_n  = ne_e/ne_b
        !> the integral should be hot to cold side, so correct the - sign afterwards 
        !> if it was cold to hot
        pe_integral   = - pe_integral 
        pressure_drop = - pressure_drop
      end if

      pe_term    = pe_integral / (Te_h * EL_CHG)
      c_sh       = sqrt(2 * Te_h * EL_CHG / (central_mass * ATOMIC_MASS_UNIT))
      
      if (solve_j) then
        if(L .ne. 0.d0) then
          gamma_fact = av_sigma * Te_h * EL_CHG / (EL_CHG ** 2 * L * ne_h * c_sh)
          j_hat_par = Newton_solver(gamma_fact, ratio_T, ratio_n, pe_term, alpha, write_debug_file)
        else 
          gamma_fact = 0.d0
          j_hat_par  = 0.d0
        end if
        
        if((ne_h .eq. 0.d0) .or. (c_sh .eq. 0.d0)) write(*,*) 'WARNING, ne_h = 0 or c_sh = 0', R_b, Z_b, Phi_b, R_e, Z_e, Phi_e, ne_b, ne_e, Te_b, Te_e, L, av_sigma, pe_term, j_hat_par, ne_h, Te_h, c_sh
        j_par = j_hat_par * (EL_CHG * ne_h * c_sh)
      else
        j_hat_par = 0.d0
        j_par     = 0.d0
      end if

      j_wall = j_par*abs(bdotn)

      ! approximately correct the radial component of the flux expansion dependence of j_par
      if(correct4fluxexpansion) then
        R_average = (R_b + R_e)/2.d0
        j_wall = j_wall*R_average/R_b 
      end if

      if (.not. hot2cold) j_wall = - j_wall 
      !< defining the sign of j_wall as positive when current flows out of the starting wall, and negative 
      !< when current flows into the starting wall; since j_par is defined as positive for hot -> 
      !< cold, if e.g. the tracer started at the hot side and j_par > 0, current is out of the 
      !< starting wall, thus j_wall is positive
      
      ! I_loc = j * A = j_wall * angle*R*dl
      I_loc = j_wall * (TWOPI/n_phi0)*R_b*initial_bnd_pos%pos(1,ip)%dl
      int_I_wall = int_I_wall + I_loc

      !$omp critical
      write(21,'(20es18.8)') R_b, Z_b, Phi_b, R_e, Z_e, Phi_e, ne_b, ne_e, Te_b, Te_e, L, av_sigma, pressure_drop, pe_term, j_hat_par, j_par, bdotn, j_wall, initial_bnd_pos%pos(1,ip)%dl, I_loc
      !$omp end critical
    enddo !traced all starting angles for one bnd point
    
    if (mod(ip,(np/10+1)) == 0) then !< writes some intermediate progress updates
      write(*,'(A17,I10,A4,I10)') 'starting point = ',ip,' of ',np
    endif

  end do !traced all starting points and angles
  !$omp end parallel do

  write(*,'(A,e18.8,A)') 'Sanity check: total current over the whole wall (i.e. discretisation error, usually of order 10 A, should be negligible compared to target current, usually of order 10 kA) = ',int_I_wall, " A"

  close(21)
  if(write_debug_file) close(22)
  if(debug) close(23)
contains

  !> detemines j_hat_par based on Newton's method with error function f
  function Newton_solver(gamma_fact, ratio_T, ratio_n, pe_term, alpha, write_debug_file) result(j_hat_par)
    
    implicit none
    
    real*8, intent(in)            :: gamma_fact, ratio_T, ratio_n, pe_term, alpha !< input parameters from the fluxlines
    real*8                        :: j_hat_par !< normalised parallel current density along the fluxline from hot to cold wall
    real*8                        :: x, x_try, delta_x, delta_x_prev !< guesses of j_hat_par, where through iterations the guess is improved
    real*8                        :: fx !< the error function evaluated at the old guess of j_hat_par
    real*8                        :: tol !< tolerance in the error of j_hat_par
    integer                       :: i, max_iterations, k
    logical,intent(in)            :: write_debug_file

    if(write_debug_file) then
      write(22,*)
      write(22,'(A30, 4es18.8)') 'Newton solver for terms ',gamma_fact, ratio_T, ratio_n, pe_term
    end if

    tol            = 1.d-10 !< note that j_hat_par < 1 due to its normalisation
    max_iterations = 500    !< there are cases where this is not enough to get small machine error in f but it's usually enought to get small error in j
   !x              = 0.1d0  !< simplistic initial guess
    !> initial guess which assumes the term with j on the RHS of Stangeby eq 17.29 is 0
    !> works well if |j_hat| << 1, otherwise the initial guess is set to 0 later on
    x = - gamma_fact * (((1.d0/ratio_T) - 1.d0)*(log(2.d0) - 0.71d0 + log(alpha)) - pe_term)

    !> checking reasonableness of initial guess
    if ((abs(x) .gt. 0.5d0) .or. (1.d0 - ratio_n*sqrt(ratio_T)*x .le. 0.d0) .or. (x .le. -1.d0)) then 
      x = 0.d0
    end if
    
    delta_x_prev = 0.5d0 !limits the first step to 0.5
    do i=1,max_iterations
      fx = f(x,gamma_fact, ratio_T, ratio_n, pe_term, alpha)

      if(write_debug_file) write(22,'(A10, I4, 3es26.16)') 'i,x,dx,fx=',i,x,delta_x_prev,fx
      
      if(fx .eq. -1.d99) then !< f(x) hit a limit in its equation
        write(*,'(A32,4es16.6)') 'Newton solver failed for values ',gamma_fact, ratio_T, ratio_n, pe_term
        j_hat_par = 0.d0
        if(write_debug_file) write(22,*) 'Solver failed due to f(x) limit'
        return
      end if
      
      if (abs(fx) .le. tol) then !< if the error is smaller than the tolerance, we found j_hat_par
        j_hat_par = x
        return !< Newton's method did what it should, leave the function
      end if
      
      ! the step to be made to improve the guess of x is
      delta_x = - fx/dfdj(x, ratio_T, ratio_n) !the full step from Newton's method
      
      if(delta_x*delta_x_prev .le. 0.d0) then !if the step being taken now is in the opposite direction to the previous step
        delta_x = sign(min(abs(delta_x),abs(0.8d0*delta_x_prev)),delta_x) !make sure the step is smaller than last step
      else
        delta_x = sign(min(abs(delta_x),abs(1.1d0*delta_x_prev)),delta_x) !scaled version to not have a much larger step than the previous step
      end if
      
      x_try = x + delta_x !< improving the guess using Newton's method

      !exponentially decreases the stepsize if the step is such that f(x) is not defined
      !due to the logarithm term
      k = 1
      do while ((1.d0 - ratio_n*sqrt(ratio_T)*x_try .le. 0.d0) .or. (x_try .le. -1.d0)) 
        if (k .gt. 100) then
          write (*,*) 'ERROR: no possible step can be made'
          j_hat_par = 0.d0
          return
        end if
        delta_x = delta_x/2.d0 !making the step smaller
        x_try = x + delta_x
        !if(write_debug_file) write(22,'(A10, I4, 2es18.8)') 'k,delta_x,x',k,delta_x,x_try
        k = k + 1
      end do
      
      x = x_try
      delta_x_prev = delta_x
      !if(abs(x) .gt. 1) then
      !  write(*,*) 'serious trouble, x=',x,fx
      !  return
      !end if
    end do
    
    !> When dfdj is very large at the solution for j, sometimes the solution has large error 
    !> (i.e. f(x) is large) but j_hat_par is still known with machine accuracy
    if(abs(delta_x) .le. 1.d-14) then 
      j_hat_par = x
      if(write_debug_file) write(22,'(A50,6es18.8)') 'Newton solver encountered difficult solution', gamma_fact, ratio_T, ratio_n, pe_term, j_hat_par, delta_x
      return
    end if

    !> if we're in this part of the code, the Newton solver wasn't able to find a solution
    write(*,'(A32,4es16.6)') 'Newton solver failed for values ',gamma_fact, ratio_T, ratio_n, pe_term
    if(write_debug_file) write(22,*) 'Solver failed after many iterations'
    j_hat_par = 0.d0
    return
    
  end function Newton_solver
  
  !> calculates Stangeby eq. 17.29 with all terms moved to the LHS so that 
  !> f(\hat{j}_\parallel) := LHS = 0 for the true \hat{j}_\parallel, and thus a 
  !> guess f(\hat{j}_{\parallel,guess}) = error. Used in Newton's method
  function f(j_hat_par, gamma_fact, ratio_T, ratio_n, pe_term, alpha) result(error)
    
    implicit none
    
    real*8, intent(in)  :: gamma_fact, ratio_T, ratio_n, pe_term, alpha !< input parameters from the fluxlines
    real*8, intent(in)  :: j_hat_par   !< normalised current density
    real*8              :: error !< the error made in the equation
    
    !> check for limits of the equation
    if((j_hat_par .le. -1.d0) .or. (ratio_n .eq. 0.d0) .or. (ratio_T .eq. 0.d0) .or. (1.d0 -ratio_n*sqrt(ratio_T)*j_hat_par .le. 0.d0)) then
      write(*,'(A20,7es16.6)') 'limit hit in f(x):', j_hat_par, gamma_fact, ratio_T, ratio_n, pe_term, alpha, 1.d0 -ratio_n*sqrt(ratio_T)*j_hat_par
      error = -1.d99
      return
    end if

    !> 
    error = j_hat_par + gamma_fact * (((1.d0/ratio_T) - 1.d0)*(log(2.d0) - 0.71d0 + log(alpha))        &
            + log((1.d0 + j_hat_par)/((1.d0 - ratio_n*sqrt(ratio_T)*j_hat_par)**(1.d0/ratio_T)))       &
            - pe_term)
    
  end function f
    
  function dfdj(j_hat_par, ratio_T, ratio_n) result(dfdj_eval_j)
    !< calculates the derivative of Stangeby eq. 17.29 with all terms moved to the LHS so that 
    ! f(\hat{j}_\parallel) := LHS = 0. Used in Newton's method
    
    implicit none
    
    real*8, intent(in) :: ratio_T, ratio_n !< ratio's of hot/cold temperature and density
    real*8, intent(in) :: j_hat_par   !< normalised current density
    real*8             :: dfdj_eval_j !< derivative evaluated at j_hat_par
    
    dfdj_eval_j = 1 + 1/(1+j_hat_par) + (ratio_n / sqrt(ratio_T))/(1-ratio_n*sqrt(ratio_T)*j_hat_par)
    
  end function dfdj

end program jorek2_SOLcurrent


!> calculates a single fieldline tracing step in s,t coordinates for a step in the phi direction of delta_p
!> copied from jorek2_poincare
subroutine step(i_elm,s_in,t_in,p_in,delta_p,delta_s,delta_t)
  use mod_parameters
  use elements_nodes_neighbours
  use phys_module
  use mod_interp
  use mod_model_settings

  implicit none

  integer :: i_elm, i_tor, i_harm

  real*8 :: s_in, t_in, p_in, delta_p, delta_s, delta_t
  real*8 :: R,R_s,R_t,Z,Z_s,Z_t
  real*8 :: Pcos,Pcos_s,Pcos_t,Pcos_st,Pcos_ss,Pcos_tt, Psin,Psin_s,Psin_t,Psin_st,Psin_ss,Psin_tt
  real*8 :: P0,P0_s,P0_t,P0_st,P0_ss,P0_tt, psi_s, psi_t, Zjac
  real*8 :: AR0_Z, AR0_p, AR0_s, AR0_t, AZ0_R, AZ0_p, AZ0_s, AZ0_t, A30_R, A30_Z, BR0, BZ0, Bp0, Fprof

  call interp_RZ(node_list,element_list,i_elm,s_in,t_in,R,R_s,R_t,Z,Z_s,Z_t)

  Zjac = (R_s * Z_t - R_t * Z_s)

  call interp(node_list,element_list,i_elm,var_psi,1,s_in,t_in,P0,P0_s,P0_t,P0_st,P0_ss,P0_tt)

  psi_s = P0_s 
  psi_t = P0_t 

#ifdef fullmhd
    call interp(node_list,element_list,i_elm,var_AR,1,s_in,t_in,P0,P0_s,P0_t,P0_st,P0_ss,P0_tt)
    AR0_s = P0_s 
    AR0_t = P0_t 

    call interp(node_list,element_list,i_elm,var_AZ,1,s_in,t_in,P0,P0_s,P0_t,P0_st,P0_ss,P0_tt)
    AZ0_s = P0_s 
    AZ0_t = P0_t 

    AR0_p = 0.d0
    AZ0_p = 0.d0

    call interp(node_list,element_list,i_elm,710,1,s_in,t_in,P0,P0_s,P0_t,P0_st,P0_ss,P0_tt)
    Fprof = P0
#endif

  do i_tor = 1, (n_tor-1)/2

    i_harm = 2*i_tor

    call interp(node_list,element_list,i_elm,var_psi,i_harm,s_in,t_in,Pcos,Pcos_s,Pcos_t,Pcos_st,Pcos_ss,Pcos_tt)

    psi_s = psi_s + Pcos_s * cos(mode(i_harm)*p_in)
    psi_t = psi_t + Pcos_t * cos(mode(i_harm)*p_in)

    call interp(node_list,element_list,i_elm,var_psi,i_harm+1,s_in,t_in,Psin,Psin_s,Psin_t,Psin_st,Psin_ss,Psin_tt)

    psi_s = psi_s + Psin_s * sin(mode(i_harm+1)*p_in)
    psi_t = psi_t + Psin_t * sin(mode(i_harm+1)*p_in)

#ifdef fullmhd
    call interp(node_list,element_list,i_elm,var_AR,i_harm,s_in,t_in,Pcos,Pcos_s,Pcos_t,Pcos_st,Pcos_ss,Pcos_tt)
    AR0_s = AR0_s + Pcos_s * cos(mode(i_harm)*p_in)
    AR0_t = AR0_t + Pcos_t * cos(mode(i_harm)*p_in)
    AR0_p = AR0_p - Pcos   * sin(mode(i_harm)*p_in) * mode(i_harm)
    call interp(node_list,element_list,i_elm,var_AR,i_harm+1,s_in,t_in,Psin,Psin_s,Psin_t,Psin_st,Psin_ss,Psin_tt)
    AR0_s = AR0_s + Psin_s * sin(mode(i_harm+1)*p_in)
    AR0_t = AR0_t + Psin_t * sin(mode(i_harm+1)*p_in)
    AR0_p = AR0_p + Psin   * cos(mode(i_harm+1)*p_in) * mode(i_harm+1)

    call interp(node_list,element_list,i_elm,var_AZ,i_harm,s_in,t_in,Pcos,Pcos_s,Pcos_t,Pcos_st,Pcos_ss,Pcos_tt)
    AZ0_s = AZ0_s + Pcos_s * cos(mode(i_harm)*p_in)
    AZ0_t = AZ0_t + Pcos_t * cos(mode(i_harm)*p_in)
    AZ0_p = AZ0_p - Pcos   * sin(mode(i_harm)*p_in) * mode(i_harm)
    call interp(node_list,element_list,i_elm,var_AZ,i_harm+1,s_in,t_in,Psin,Psin_s,Psin_t,Psin_st,Psin_ss,Psin_tt)
    AZ0_s = AZ0_s + Psin_s * sin(mode(i_harm+1)*p_in)
    AZ0_t = AZ0_t + Psin_t * sin(mode(i_harm+1)*p_in)
    AZ0_p = AZ0_p + Psin   * cos(mode(i_harm+1)*p_in) * mode(i_harm+1)
#endif

  enddo

#ifdef fullmhd
  AR0_Z = ( - R_t * AR0_s  + R_s * AR0_t ) / Zjac
  AZ0_R = (   Z_t * AZ0_s  - Z_s * AZ0_t ) / Zjac
  A30_R = (   Z_t * psi_s  - Z_s * psi_t ) / Zjac
  A30_Z = ( - R_t * psi_s  + R_s * psi_t ) / Zjac

  BR0 = ( A30_Z - AZ0_p )/ R
  BZ0 = ( AR0_p - A30_R )/ R
  Bp0 = ( AZ0_R - AR0_Z )       +   Fprof / R

  ! dR/Rdphi = B_R / B_phi ; dZ/Rdphi = B_Z / B_phi
  ! ds = (Z_t dR - R_t dZ) / Zjac ; dt = ( -Z_s dR + R_s dZ) / Zjac
  delta_s =  ( Z_t*BR0 - R_t*BZ0) / ( Bp0 * Zjac ) * R * delta_p
  delta_t =  (-Z_s*BR0 + R_s*BZ0) / ( Bp0 * Zjac ) * R * delta_p     

#else
  delta_s =   psi_t * R / (Zjac * F0) * delta_p
  delta_t = - psi_s * R / (Zjac * F0) * delta_p
#endif

  return
end subroutine step
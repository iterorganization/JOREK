!> Newton methods for root solutions
module mod_newton_methods

  contains



  subroutine find_variable_minmax(node_list,element_list,i_elm, i_var, var_min, var_max)

    use data_structure
    use mod_interp, only: interp
    use mod_parameters, only: n_vertex_max
    implicit none

    ! --- Routine variables
    type (type_node_list),    intent(in) :: node_list
    type (type_element_list), intent(in) :: element_list
    integer, intent(in)      :: i_elm                          ! element on which we search
    integer, intent(in)      :: i_var                          ! index of the variable involved (-1 for R and -2 for Z)
    real*8,  intent(out)     :: var_min, var_max               ! min max values along edges

    ! --- Internal variables
    real*8    :: psimin, psimax
    integer   :: i_side, n1, n2

    var_min = +1.d15
    var_max = -1.d15

    do i_side=1,4

      ! --- careful with axis
      n1 = element_list%element(i_elm)%vertex(i_side)
      n2 = element_list%element(i_elm)%vertex(mod(i_side,n_vertex_max) + 1)
      if (node_list%node(n1)%axis_node .and. node_list%node(n2)%axis_node) cycle

      call find_variable_minmax_on_edge(node_list,element_list,i_elm, i_side, i_var, psimin, psimax)
      var_min = min(var_min,psimin)
      var_max = max(var_max,psimax)

    enddo

    return
  end subroutine find_variable_minmax



  subroutine find_variable_minmax_on_edge(node_list,element_list,i_elm, i_side, i_var, var_min, var_max)

    use data_structure
    use mod_interp, only: interp, interp_PRZ_combined
    use mod_parameters, only: n_order
    implicit none

    ! --- Routine variables
    type (type_node_list),    intent(in) :: node_list
    type (type_element_list), intent(in) :: element_list
    integer, intent(in)      :: i_elm                          ! element on which we search
    integer, intent(in)      :: i_side                         ! side of the element along which we search
    integer, intent(in)      :: i_var                          ! index of the variable involved (-1 for R and -2 for Z)
    real*8,  intent(out)     :: var_min, var_max               ! min max values along edges

    ! --- Internal variables
    real*8    :: psi, psi_s, psi_t, psi_st, psi_ss, psi_tt
    real*8    :: R,dR_ds,dR_dt,dR_dst,dR_dss,dR_dtt
    real*8    :: Z,dZ_ds,dZ_dt,dZ_dst,dZ_dss,dZ_dtt
    real*8    :: step, s_or_t
    integer   :: i
    integer   :: n_found
    real*8    :: inflex(n_order)

    var_min = +1.d15
    var_max = -1.d15

    if (i_side .eq. 1) s_or_t = 0.d0 ! side t=0
    if (i_side .eq. 2) s_or_t = 1.d0 ! side s=1
    if (i_side .eq. 3) s_or_t = 1.d0 ! side t=1
    if (i_side .eq. 4) s_or_t = 0.d0 ! side s=0

    ! --- First the end points
    do i=1,2

      if (i .eq. 1) step = 0.0d0
      if (i .eq. 2) step = 1.0d0

      if ( mod(i_side,2) .eq. 0 ) then
        call interp_PRZ_combined(node_list,element_list,i_elm,i_var,1,s_or_t,step, psi, psi_s, psi_t, psi_st, psi_ss, psi_tt)
      else
        call interp_PRZ_combined(node_list,element_list,i_elm,i_var,1,step,s_or_t, psi, psi_s, psi_t, psi_st, psi_ss, psi_tt)
      endif

      var_min = min(var_min,psi)
      var_max = max(var_max,psi)

    enddo

    ! --- Then look for inflexions along the edge
    call newton_1D_inflexions(node_list,element_list,i_elm, i_side, i_var, n_found, inflex)
    do i=1,n_found
      step = inflex(i)
      if ( mod(i_side,2) .eq. 0 ) then
        call interp_PRZ_combined(node_list,element_list,i_elm,i_var,1,s_or_t,step, psi, psi_s, psi_t, psi_st, psi_ss, psi_tt)
      else
        call interp_PRZ_combined(node_list,element_list,i_elm,i_var,1,step,s_or_t, psi, psi_s, psi_t, psi_st, psi_ss, psi_tt)
      endif
      var_min = min(var_min,psi)
      var_max = max(var_max,psi)
    enddo

    return
  end subroutine find_variable_minmax_on_edge






  ! --------------------------------------------------------------------------------------------------------
  ! Newton method to find inflexions (df=0) of variable (n=0) along element side
  subroutine newton_1D_inflexions(node_list,element_list,i_elm, i_side, i_var, n_found, inflex)

    use data_structure
    use mod_interp, only: interp, interp_PRZ_combined
    use mod_parameters, only: n_order
    implicit none
  
    ! --- Routine variables
    type (type_node_list),    intent(in) :: node_list
    type (type_element_list), intent(in) :: element_list
    integer, intent(in)      :: i_elm                          ! element on which we search
    integer, intent(in)      :: i_side                         ! side of the element along which we search
    integer, intent(in)      :: i_var                          ! index of the variable involved (-1 for R and -2 for Z)
    integer, intent(out)     :: n_found                        ! number of inflexions found
    real*8,  intent(out)     :: inflex(n_order)                ! inflexions found (coordinate along element side)

    ! --- Internal variables
    real*8    :: psi, psi_s, psi_t, psi_st, psi_ss, psi_tt
    real*8    :: tolf, tolx
    real*8    :: step, s_or_t, jump, df, d2f, errx, errf
    integer   :: ntrial, i, k
  
    ntrial = 20
    tolx = 1.d-8
    tolf = 1.d-16
  
    n_found = 0
  
    if (i_side .eq. 1) s_or_t = 0.d0 ! side t=0
    if (i_side .eq. 2) s_or_t = 1.d0 ! side s=1
    if (i_side .eq. 3) s_or_t = 1.d0 ! side t=1
    if (i_side .eq. 4) s_or_t = 0.d0 ! side s=0

    ! --- We do 3 trials, starting from the edges, and from the middle
    ! --- If Newton get us out of the interval [0,1] each time, we assume there is no inflexion
    do i=1,3
      if (i .eq. 1) step = 0.0d0
      if (i .eq. 2) step = 1.0d0
      if (i .eq. 3) step = 0.5d0
      do k=1,ntrial
   
        ! --- Get 1st and 2nd derivative
        if ( mod(i_side,2) .eq. 0 ) then
          call interp_PRZ_combined(node_list,element_list,i_elm,i_var,1,s_or_t,step, psi, psi_s, psi_t, psi_st, psi_ss, psi_tt)
          df  = psi_t
          d2f = psi_tt
        else
          call interp_PRZ_combined(node_list,element_list,i_elm,i_var,1,step,s_or_t, psi, psi_s, psi_t, psi_st, psi_ss, psi_tt)
          df  = psi_s
          d2f = psi_ss
        endif
   
        errf=abs(df)
   
        if (errf .le. tolf) then
          n_found = n_found + 1
          inflex(n_found) = step
          exit
        endif
     
        if (d2f .ne. 0.d0) then
          jump = - df/d2f
        else
          jump = - df/abs(df) * 0.25
        endif
        errx = abs(jump)
        jump = min(jump,+0.25d0)
        jump = max(jump,-0.25d0)
     
        step = step + jump
     
        ! --- If we get out of the element, it probably means there is no inflexion
        if ( (step .lt. 0.d0) .or. (step .gt. 1.d0) ) exit
   
        if (errx .le. tolx) then
          n_found = n_found + 1
          inflex(n_found) = step
          exit
        endif
   
      enddo ! ntrials
    enddo ! edges_and_middle

    ! --- If these trials didn't find anything, just get out
    if (n_found .eq. 0) return

    ! --- Otherwise, there can be several possibilities, depending on the order of the elements

    if (n_found .eq. 2) then
      ! --- If both trials found the same inflexion point, we assume it's the only one
      if (abs(inflex(1)-inflex(2)) .lt. 2.0*tolx) then
        n_found = 1
        inflex(1) = ( inflex(1) + inflex(2) ) / 2.d0
        return
      ! --- If we found two different inflexion points when starting from the edges
      ! --- and the elements are bi-cubic, then these are our two inflextions, there won't be more
      else
        if (n_order .eq. 3) return
      endif
    endif
 
    if (n_found .eq. 3) then
      ! --- If all trials found the same inflexion point, we assume it's the only one
      if ( (abs(inflex(1)-inflex(2)) .lt. 2.0*tolx) .and. (abs(inflex(2)-inflex(3)) .lt. 2.0*tolx) ) then
        n_found = 1
        inflex(1) = ( inflex(1) + inflex(2) + inflex(3) ) / 3.d0
        return
      ! --- If we found two different inflexion points when starting from the edges
      ! --- and the elements are bi-cubic, then these are our two inflextions, there won't be more
      else
        if (abs(inflex(1)-inflex(2)) .lt. 2.0*tolx) then
          n_found = 2
          inflex(1) = ( inflex(1) + inflex(2) ) / 2.d0
          inflex(2) = inflex(3)
          return
        endif
        if (abs(inflex(2)-inflex(3)) .lt. 2.0*tolx) then
          n_found = 2
          inflex(2) = ( inflex(2) + inflex(3) ) / 2.d0
          return
        endif
        if (abs(inflex(1)-inflex(3)) .lt. 2.0*tolx) then
          n_found = 2
          inflex(1) = ( inflex(1) + inflex(3) ) / 2.d0
          return
        endif
        if (n_order .eq. 3) return
      endif
    endif
 
    ! --- If we found only one inflexion point when starting from the edge
    ! --- and the elements are bi-cubic, then this means the second inflexion point is outside the element
    if ( (n_found .eq. 1) .and. (n_order .eq. 3) ) return

    ! --- If we found two different inflexion points when starting from the edges
    ! --- and the elements are bi-quintic, then there are potentially two more
    ! --- inflextions between these two points

    return
  end subroutine newton_1D_inflexions





  ! --------------------------------------------------------------------------------------------------------
  ! Newton method to find value of variable on the side of an element
  ! Assumes you have already checked the minmax values on element!
  subroutine newton_1D_find_value_on_element_side(node_list,element_list,i_elm, i_side, i_var, var_find, n_found, st_found)

    use data_structure
    use mod_interp, only: interp, interp_PRZ_combined
    use mod_parameters, only: n_order
    implicit none
  
    ! --- Routine variables
    type (type_node_list),    intent(in) :: node_list
    type (type_element_list), intent(in) :: element_list
    integer, intent(in)      :: i_elm                          ! element on which we search
    integer, intent(in)      :: i_side                         ! side of the element along which we search
    integer, intent(in)      :: i_var                          ! index of the variable involved (-1 for R and -2 for Z)
    real*8,  intent(in)      :: var_find                       ! value of variable to locate inside interval
    integer, intent(out)     :: n_found                        ! number of solutions found
    real*8,  intent(out)     :: st_found(n_order)              ! points found

    ! --- Internal variables
    real*8    :: psi, psi_s, psi_t, psi_st, psi_ss, psi_tt
    real*8    :: tolx
    real*8    :: step, s_or_t, st_tmp, jump, ff, df
    real*8    :: st_min, st_max
    integer   :: ntrial, i, j, k, n_inflex, n_tmp
    real*8    :: inflex(n_order), inflex_tmp(n_order), st_single, st_found_tmp(n_order)
    integer   :: i_min
    real*8    :: var_min
    logical   :: found_duplicate
 
    tolx = 1.d-8
    n_found = 0

    ! --- Check if there are any inflexion points
    call newton_1D_inflexions(node_list,element_list,i_elm, i_side, i_var, n_inflex, inflex)

    ! --- If there are no inflexion points, easy
    if (n_inflex .eq. 0) then
      st_min = 0.d0
      st_max = 1.d0
      call newton_1D_find_value_inside_interval(node_list,element_list,i_elm, i_side, i_var, st_min, st_max, var_find, st_single)
      if (st_single .ge. 0.d0) then
        n_found = 1
        st_found(1) = st_single
      endif
      return
    endif

    ! --- If there are inflexion points, we look at each interval
    if (n_inflex .ge. 1) then
      ! --- reorder inflexion points
      if (n_inflex .gt. 1) then
        inflex_tmp(1:n_inflex) = inflex(1:n_inflex)
        do i=1,n_inflex
          var_min = 1.d15
          i_min   = 0
          do j=1,n_inflex
            if (inflex_tmp(j) .lt. var_min) then
              var_min = inflex_tmp(j)
              i_min   = j
            endif
          enddo
          inflex(i) = var_min
          inflex_tmp(i_min) = 999.0
        enddo
      endif
      ! --- Each interval
      do i=1,n_inflex+1
        if (i .eq. 1) then
          st_min = 0.d0
          st_max = inflex(i)
        elseif (i .eq. n_inflex+1) then
          st_min = inflex(n_inflex)
          st_max = 1.d0
        else
          st_min = inflex(i-1)
          st_max = inflex(i)
        endif
        call newton_1D_find_value_inside_interval(node_list,element_list,i_elm, i_side, i_var, st_min, st_max, var_find, st_single)
        if (st_single .ge. 0.d0) then
          n_found = n_found + 1
          st_found(n_found) = st_single
        endif
      enddo
    endif

    ! --- Clean-up solutions in case we found several times the same (shouldn't happen)
    if (n_found .gt. 1) then
      st_found_tmp(1:n_found) = st_found(1:n_found)
      n_tmp = 0
      do i=1,n_found
        found_duplicate = .false.
        do j=i+1,n_found
          if ( abs(st_found_tmp(i)-st_found_tmp(j)) .lt. 2.0*tolx ) found_duplicate = .true.
        enddo
        if (.not. found_duplicate) then
          n_tmp = n_tmp + 1
          st_found(n_tmp) = st_found_tmp(i)
        endif
      enddo
      n_found = n_tmp
    endif

    return
  end subroutine newton_1D_find_value_on_element_side






  ! --------------------------------------------------------------------------------------------------------
  ! Newton method to find value of variable between two point on the side of an element
  subroutine newton_1D_find_value_inside_interval(node_list,element_list,i_elm, i_side, i_var, st_min, st_max, var_find, st_found)

    use data_structure
    use mod_interp, only: interp, interp_PRZ_combined
    use mod_parameters, only: n_order
    implicit none
  
    ! --- Routine variables
    type (type_node_list),    intent(in) :: node_list
    type (type_element_list), intent(in) :: element_list
    integer, intent(in)      :: i_elm                          ! element on which we search
    integer, intent(in)      :: i_side                         ! side of the element along which we search
    integer, intent(in)      :: i_var                          ! index of the variable involved (-1 for R and -2 for Z)
    real*8,  intent(inout)   :: st_min, st_max                 ! min/max values of interval (should be inside [0,1])
    real*8,  intent(in)      :: var_find                       ! value of variable to locate inside interval
    real*8,  intent(out)     :: st_found                       ! point found

    ! --- Internal variables
    real*8    :: psi, psi_s, psi_t, psi_st, psi_ss, psi_tt
    real*8    :: tolf, tolx, var_min, var_max, var_tmp
    real*8    :: step, s_or_t, st_tmp, jump, ff, df, errx, errf
    integer   :: ntrial, i, k
 
    st_found = -999.0 ! ie. outisde
  
    if (st_min .gt. st_max) then
      st_tmp = st_min
      st_min = st_max
      st_max = st_tmp
    endif

    ntrial = 20
    tolx = 1.d-8
    tolf = 1.d-16
  
    if (i_side .eq. 1) s_or_t = 0.d0 ! side t=0
    if (i_side .eq. 2) s_or_t = 1.d0 ! side s=1
    if (i_side .eq. 3) s_or_t = 1.d0 ! side t=1
    if (i_side .eq. 4) s_or_t = 0.d0 ! side s=0

    ! --- Sanity check
    var_min = +1.d15
    var_max = -1.d15
    step = st_min
    if ( mod(i_side,2) .eq. 0 ) then
      call interp_PRZ_combined(node_list,element_list,i_elm,i_var,1,s_or_t,step, var_tmp, psi_s, psi_t, psi_st, psi_ss, psi_tt)
    else
      call interp_PRZ_combined(node_list,element_list,i_elm,i_var,1,step,s_or_t, var_tmp, psi_s, psi_t, psi_st, psi_ss, psi_tt)
    endif
    var_min = min(var_min,var_tmp)
    var_max = max(var_max,var_tmp)
    step = st_max
    if ( mod(i_side,2) .eq. 0 ) then
      call interp_PRZ_combined(node_list,element_list,i_elm,i_var,1,s_or_t,step, var_tmp, psi_s, psi_t, psi_st, psi_ss, psi_tt)
    else
      call interp_PRZ_combined(node_list,element_list,i_elm,i_var,1,step,s_or_t, var_tmp, psi_s, psi_t, psi_st, psi_ss, psi_tt)
    endif
    var_min = min(var_min,var_tmp)
    var_max = max(var_max,var_tmp)
    if ( (var_find .lt. var_min) .or. (var_find .gt. var_max) ) return

    ! --- Start between the two points
    step = 0.5 * (st_min + st_max)
    do k=1,ntrial
   
      ! --- Get value and 1st derivative
      if ( mod(i_side,2) .eq. 0 ) then
        call interp_PRZ_combined(node_list,element_list,i_elm,i_var,1,s_or_t,step, psi, psi_s, psi_t, psi_st, psi_ss, psi_tt)
        ff = psi - var_find
        df = psi_t
      else
        call interp_PRZ_combined(node_list,element_list,i_elm,i_var,1,step,s_or_t, psi, psi_s, psi_t, psi_st, psi_ss, psi_tt)
        ff = psi - var_find
        df = psi_s
      endif
   
      errf=abs(ff)
   
      if (errf .le. tolf) then
        st_found = step
        return
      endif
    
      if (df .ne. 0.d0) then
        jump = - ff/df
      else
        jump = - ff/abs(ff) * 0.25 * (st_max-st_min)
      endif
      errx = abs(jump)
      jump = min(jump,+0.25d0*(st_max-st_min))
      jump = max(jump,-0.25d0*(st_max-st_min))
    
      step = step + jump
    
      ! --- If we get out of the element, it probably means there is no inflexion
      if ( (step .lt. st_min) .or. (step .gt. st_max) ) exit
   
      if (errx .le. tolx) then
        st_found = step
        return
      endif
   
    enddo ! ntrials

  end subroutine newton_1D_find_value_inside_interval




end module mod_newton_methods

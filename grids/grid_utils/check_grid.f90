subroutine check_grid(my_id, node_list, element_list)

use mod_parameters
use data_structure
use mod_basisfunctions
use mod_export_restart

implicit none

! --- Routine parameters
integer,                    intent(in) :: my_id
type(type_node_list),       intent(in) :: node_list
type(type_element_list),    intent(in) :: element_list

! --- Hard-coded parameters: where to check inside an element
integer, parameter  :: NST = 5
real*8,  parameter  :: ST(NST) = (/ 0.01d0, 0.25d0, 0.50d0, 0.75d0, 0.99d0 /)

! --- Local variables
type (type_element) :: element
type (type_node)    :: node
integer             :: ielm, ivert, idof, ms, mt
real*8              :: x_g, x_s, x_t, y_g, y_s, y_t, xjac, xjac_min, xjac_max, s, t
logical             :: problem_found
real*8, dimension(4,n_degrees,NST,NST) :: H, H_s, H_t

problem_found = .false.

! Remark: only do this check on my_id=0 for the moment, can be MPI parallelized later on
if ( my_id /= 0 ) return

! --------------------------------------------------------------------------------------------------
! --- Check whether the Jacobian changes sign inside an element ------------------------------------
! --------------------------------------------------------------------------------------------------

do ms = 1, NST
  s = ST(ms)
  do mt = 1, NST
    t = ST(mt)
    call basisfunctions(s, t, H(:,:,ms,mt), H_s(:,:,ms,mt), H_t(:,:,ms,mt))
  end do
end do

do ielm = 1, element_list%n_elements
  element = element_list%element(ielm)
  xjac_min = +1.d99
  xjac_max = -1.d99
  do ms = 1, NST
    do mt = 1, NST
      x_g  = 0.d0; x_s  = 0.d0; x_t  = 0.d0; y_g  = 0.d0; y_s  = 0.d0; y_t  = 0.d0
      do ivert = 1, n_vertex_max
        node    = node_list%node(element%vertex(ivert))
        do idof = 1, n_degrees
          x_g = x_g  + node%x(1,idof,1) * element%size(ivert,idof) * H  (ivert,idof,ms,mt)
          x_s = x_s  + node%x(1,idof,1) * element%size(ivert,idof) * H_s(ivert,idof,ms,mt)
          x_t = x_t  + node%x(1,idof,1) * element%size(ivert,idof) * H_t(ivert,idof,ms,mt)
          y_g  = y_g + node%x(1,idof,2) * element%size(ivert,idof) * H  (ivert,idof,ms,mt)
          y_s  = y_s + node%x(1,idof,2) * element%size(ivert,idof) * H_s(ivert,idof,ms,mt)
          y_t  = y_t + node%x(1,idof,2) * element%size(ivert,idof) * H_t(ivert,idof,ms,mt)
        end do
      end do
      xjac     = x_s*y_t  - x_t*y_s
      xjac_min = min( xjac_min, xjac )
      xjac_max = max( xjac_max, xjac )
    end do
  end do
  if ( xjac_min * xjac_max < 0.d0 ) then
    write(*,*) 'ERROR in check_grid: Jacobian changes sign inside element!'
    write(*,'(1x,a,i7,a,f9.4,a,f9.4)') 'Element ', ielm, ' near R=',x_g, ' Z=', y_g
    write(*,*) xjac_min, xjac_max
    problem_found = .true.
  end if
end do

if ( problem_found ) then
  call export_restart(node_list, element_list, 'jorek_stopped')
  stop
else
  write(*,*) 'Routine check_grid did not find issues.'
end if

end subroutine check_grid

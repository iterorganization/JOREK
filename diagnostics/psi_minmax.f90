subroutine psi_minmax(node_list,element_list,i_elm,psimin,psimax)

use data_structure
use mod_newton_methods
use mod_parameters, only: n_order

implicit none

type (type_node_list)    :: node_list
type (type_element_list) :: element_list
type (type_surface_list) :: surface_list

real*8  :: psimin, psimax, psma, psmi, psmima, psim, psimr, psip, psipr
real*8  :: aa, bb, cc, det, r, dummy
real*8,external :: root
integer :: i_elm, iv, n, im, n1, n2
real*8  :: s,t,P,P_s,P_t,P_st,P_ss,P_tt
integer :: k

! --- For n_order>3, we need to use Newton methods (not exactly true, should implement quartic root finder) 
if (n_order .ge. 5) then
  call find_variable_minmax(node_list,element_list,i_elm, var_psi, psimin, psimax)
  return
endif

! --- Continue for bi-cubic elements

psimin = 1d10
psimax =-1d10

do iv= 1, n_vertex_max

  im = mod(iv,n_vertex_max) + 1
  n1 = element_list%element(i_elm)%vertex(iv)
  n2 = element_list%element(i_elm)%vertex(im)

  if (node_list%node(n1)%axis_node .and. node_list%node(n2)%axis_node) cycle

  if ((iv .eq. 1) .or. (iv .eq. 3)) THEN

    PSIM  =  node_list%node(n1)%values(1,1,1) * element_list%element(i_elm)%size(iv,1)             ! PSI(1,n1)
    PSIMR =  node_list%node(n1)%values(1,2,1) * element_list%element(i_elm)%size(iv,2) * 3.d0/2.d0 ! PSI(2,n1)
    PSIP  =  node_list%node(n2)%values(1,1,1) * element_list%element(i_elm)%size(im,1)             ! PSI(1,n2)
    PSIPR = -node_list%node(n2)%values(1,2,1) * element_list%element(i_elm)%size(im,2) * 3.d0/2.d0 ! PSI(2,n2)

  elseif ((iv .eq. 2) .or. (iv .eq. 4)) then

    PSIM  =   node_list%node(n1)%values(1,1,1) * element_list%element(i_elm)%size(iv,1)             ! PSI(1,n1)
    PSIMR =   node_list%node(n1)%values(1,3,1) * element_list%element(i_elm)%size(iv,3) * 3.d0/2.d0 ! PSI(3,n1)
    PSIP  =   node_list%node(n2)%values(1,1,1) * element_list%element(i_elm)%size(im,1)             ! PSI(1,n2)
    PSIPR = - node_list%node(n2)%values(1,3,1) * element_list%element(i_elm)%size(im,3) * 3.d0/2.d0 ! PSI(3,n2)

  endif

  PSMA = MAX(PSIM,PSIP)
  PSMI = MIN(PSIM,PSIP)
  AA =  3.d0 * (PSIM + PSIMR - PSIP + PSIPR ) / 4.d0
  BB =  ( - PSIMR + PSIPR ) / 2.d0
  CC =  ( - 3.d0*PSIM - PSIMR + 3.d0*PSIP - PSIPR) / 4.d0
  DET = BB**2 - 4.d0*AA*CC
  IF (DET .GE. 0.d0) THEN
    R = ROOT(AA,BB,CC,DET,1.d0)
    IF (ABS(R) .GT. 1.d0) THEN
      R = ROOT(AA,BB,CC,DET,-1.d0)
    ENDIF
    IF (ABS(R) .LE. 1.d0) THEN
      CALL CUB1D(PSIM,PSIMR,PSIP,PSIPR,R,PSMIMA,DUMMY)
      psma = max(psma,psmima)
      psmi = min(psmi,psmima)
    ENDIF
  ENDIF
  psimin = min(psimin,psmi)
  psimax = max(psimax,psma)

ENDDO

RETURN
END

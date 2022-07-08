# #!/usr/bin/python3
import sys
import numpy
import pylab
import math
import random
from scipy import interpolate
import sympy as sym



# ------------------------------------------------------
# ------------------------------------------------------
# ------- Section: Global variables for all functions
# ------------------------------------------------------
# ------------------------------------------------------



# --- Order of Bezier polynomial basis
n_order  = 7

# --- Do we want to use the same vector names as in the paper?
# --- ie. (u,v,w,i,j,m,n,k)
vectors_like_paper = True

# --- Derived sizes
n_vertex = 4
n_degrees = int((n_order+1)**2/4)
n_vectors = int((n_order+1)/2)

# --- This is important. In principle, we should never need basis functions derivatives
# --- higher than _st, _ss, _tt. There is only one place: when aligning a grid to the psi-contours
# --- However, at n_order>5, instead of using interp(psi) to get the derivatives, we will use the
# --- surface spline of surface_list
n_deriv = 3 #(n_order+1)/2


# --- The ordering of nodes as (u,v,w,i,j,m,n,k) is a bit tricky
count = 0
basis_index  = [[0 for j in range(n_vectors)] for i in range(n_vectors)] # nodes index from 1 to 4
# --- We do it square by square
for k in range(n_vectors):
    # --- For each square, we do i row and j column
    for l in range(k+1):
        # --- First the row
        j=l
        i = k
        basis_index[i][j] = count
        count = count + 1
        # --- Then the column
        i=l
        j = k
        if (i==j): continue # don't record corners twice
        basis_index[i][j] = count
        count = count + 1
# --- Printout nodes?
if (False):
    for j in range(n_vectors-1,-1,-1):
        for i in range(n_vectors):
            sys.stdout.write(" %3d " % (basis_index[i][j]))
            if (i == n_vectors-1): sys.stdout.write('\n\n')

# --- Name of Fortran file
filename = "mod_basisfunctions.f90"


# --- Formalism (python or Fortran)
python_formalism = False
if (python_formalism):
    newline = "\\"
    Lbrack  = "["
    Rbrack  = "]"
    array_index = 0
    filename = "python_basis.txt"
else:
    newline = "&"
    Lbrack  = "("
    Rbrack  = ")"
    array_index = 1







# ------------------------------------------------------
# ------------------------------------------------------
# ------- Section: Bezier Functions
# ------------------------------------------------------
# ------------------------------------------------------




def fact(i): # --- The factorial
    return math.factorial(i)
def b_coef(i,n): # --- The binomial coef for Bezier polynomials
    return float(fact(n)) / float(fact(i)*fact(n-i))
def bezier_coef(s, t, i, j, n_order): # --- The Bezier formula
    return b_coef(i,n_order) * s**i * (1.0-s)**(n_order-i) * b_coef(j,n_order) * t**j * (1.0-t)**(n_order-j)






# ------------------------------------------------------
# ------------------------------------------------------
# ------- Section: Write basis functions and derivatives
# ------------------------------------------------------
# ------------------------------------------------------


# --- Write 2D basis functions, with derivatives between deriv_min and deriv_max
# --- where i_deriv=0 means no derivatives (ie. the value)
def write_basis_function(F_basis, deriv_min, deriv_max_tmp, transpose, st_only):
    deriv_max = deriv_max_tmp + 1
    if (st_only):
        deriv_min = 1
        deriv_max = 2
    basis_file = open(filename,'a')
    # --- Write basis functions for each derivative
    for d_s in range(0,deriv_max):
        for d_t in range(0,deriv_max):
            if (st_only):
                if ( (d_s !=1 ) or (d_t != 1) ): continue
            else:
                if (d_s + d_t < deriv_min): continue # only up to _ss, _tt, _st
                if (d_s + d_t > deriv_max-1): continue # only up to _ss, _tt, _st
            # --- Define the derivative string (eg. _st) and offset for nicely aligned terms
            offset = "              "
            deriv = ""
            if ( (d_s!=0) or (d_t!=0) ):
                deriv = "_"
                offset = offset + " "
                for i in range(d_s):
                    deriv = deriv + "s"
                    offset = offset + " "
                for i in range(d_t):
                    deriv = deriv + "t"
                    offset = offset + " "
            # --- Write basis functions for each node
            for k in range(n_vertex):
                if ( (d_s==0) and (d_t==0) ):
                    basis_file.write("\n  ! --- Main values on node %d" %(k+1))
                else:
                    basis_file.write("\n  ! --- Derivatives %s on node %d" %(deriv,k+1))
                # --- Loop over each degree of the basis function
                for i in range(n_vectors):
                    for j in range(n_vectors):
                        i_deg = basis_index[i][j]
                        # --- Name of basis function
                        basis_file.write("\n")
                        if (transpose):
                            basis_file.write("  H%s%s%2d,%2d%s =  " % (deriv,Lbrack,i_deg+array_index,k+array_index,Rbrack))
                        else:
                            basis_file.write("  H%s%s%2d,%2d%s =  " % (deriv,Lbrack,k+array_index,i_deg+array_index,Rbrack))
                        # --- Write all terms
                        basis_string = str(F_basis[k][i_deg][d_s][d_t])
                        basis_terms = cleanup_string(basis_string)
                        n_max = 2 ; count = 0 # allow up to 2 terms per line (1 for derivatives)
                        for term in (basis_terms):
                            if (term.strip() != ""):
                                if (count == n_max):
                                    basis_file.write(" "+newline+"\n")
                                    basis_file.write(offset)
                                    count = 0
                                basis_file.write(term)
                                count = count + 1
                basis_file.write("\n")
    basis_file.close()








# --- Write 1D basis functions, with derivatives between deriv_min and deriv_max
def write_basis_function_1D(F_basis, deriv_min, deriv_max_tmp):
    deriv_max = deriv_max_tmp + 1
    basis_file = open(filename,'a')
    # --- Write basis functions for each derivative
    for d_s in range(0,deriv_max):
        d_t = 0
        if (d_s + d_t < deriv_min): continue # only up to _ss, _tt, _st
        if (d_s + d_t > deriv_max-1): continue # only up to _ss, _tt, _st
        # --- Define the derivative string (eg. _st) and offset for nicely aligned terms
        offset = "              "
        deriv = ""
        if ( (d_s!=0) or (d_t!=0) ):
            deriv = "_"
            offset = offset + " "
            for i in range(d_s):
                deriv = deriv + "s"
                offset = offset + " "
        # --- Write basis functions for each node
        for k in range(int(n_vertex/2)):
            if ( (d_s==0) and (d_t==0) ):
                basis_file.write("\n  ! --- Main values on node %d" %(k+1))
            else:
                basis_file.write("\n  ! --- Derivatives %s on node %d" %(deriv,k+1))
            # --- Loop over each degree of the basis function
            for i in range(n_vectors):
                j = 0
                i_deg = i
                # --- Name of basis function
                basis_file.write("\n")
                basis_file.write("  H%s%s%2d,%2d%s =  " % (deriv,Lbrack,k+array_index,i_deg+array_index,Rbrack))
                basis_string = str(F_basis[k][i_deg][d_s])
                if (basis_string == ""): continue
                # --- Write all terms
                basis_terms = cleanup_string(basis_string)
                n_max = 2 ; count = 0 # allow up to 2 terms per line (1 for derivatives)
                for term in (basis_terms):
                    if (term.strip() != ""):
                        if (count == n_max):
                            basis_file.write(" "+newline+"\n")
                            basis_file.write(offset)
                            count = 0
                        basis_file.write(term)
                        count = count + 1
            basis_file.write("\n")
    basis_file.close()




# --- Starting from a Sympy string, get a clean list of terms with spaced out *
def cleanup_string(basis_string):
    basis_terms = []
    count = 0
    term_tmp = ""
    c = basis_string[0]
    count = count + 1
    while(c != ""):
        c = basis_string[count]
        cm1 = basis_string[count-1]
        if (count < len(basis_string) - 1): cp1 = basis_string[count+1]
        if (c == "("):
            found_additional_parenthesis = 0
            term_tmp = term_tmp + c
            count = count + 1
            while(True):
                c = basis_string[count]
                if (count < len(basis_string) - 1): cp1 = basis_string[count+1]
                if ( (c == ".") and (cp1 == "0") ):
                    term_tmp = term_tmp + ".d"
                else:
                    term_tmp = term_tmp + c
                count = count + 1
                if (c == "("):
                    found_additional_parenthesis = found_additional_parenthesis + 1
                if (c == ")"):
                    if (found_additional_parenthesis == 0):
                        break
                    else:
                        found_additional_parenthesis = found_additional_parenthesis - 1
        else:
            if ( (c == "*") and (cm1 != "*") and (cp1 != "*") ):
                term_tmp = term_tmp + " * "
            elif ( (c == "+") or (c == "-") ):
                basis_terms.append(term_tmp)
                term_tmp = "  "+c+"  "
            elif ( (c == ".") and (cp1 == "0") ):
                term_tmp = term_tmp + ".d"
            else:
                term_tmp = term_tmp + c
            count = count + 1
        if (count > len(basis_string) - 1):
            basis_terms.append(term_tmp)
            break
    return basis_terms






# ------------------------------------------------------
# ------------------------------------------------------
# ------- Section: Each Fortran routine in the module
# ------------------------------------------------------
# ------------------------------------------------------



# --- The headers of the Fortran module
def write_basis_function_headers():
    basis_file = open(filename,'w')
    
    basis_file.write("!> Module providing the different implementations of the basis functions derived from a\n")
    basis_file.write("!> mixed Bezier/Cubic finite element representation.\n")
    basis_file.write("!> \see ::basisfunctions and ::basisfunctions1\n")
    basis_file.write("!> IMPORTANT: auto-generated code for n_order=%d\n" % (n_order))
    basis_file.write("module mod_basisfunctions\n")
    basis_file.write("\n")
    basis_file.write("use mod_parameters, only: n_order, n_degrees\n")
    basis_file.write("\n")
    basis_file.write("implicit none\n")
    basis_file.write("\n")
    basis_file.write("private\n")
    basis_file.write("public :: basisfunctions1, basisfunctions\n")
    basis_file.write("public :: basisfunctions_T !< Transposed version, for faster interp_PRZ\n")
    basis_file.write("public :: basisfunctions3\n")
    basis_file.write("\n")
    basis_file.write("!> One-dimensional basisfunctions with derivatives of order n\n")
    basis_file.write("interface basisfunctions1\n")
    basis_file.write("  module procedure basisfunctions_1D_0, basisfunctions_1D_1, basisfunctions_1D_2\n")
    basis_file.write("end interface basisfunctions1\n")
    basis_file.write("\n")
    basis_file.write("! Two-dimensional basisfunctions with derivatives of order n\n")
    basis_file.write("! basisfunctions_2D_1 has only s and t derivatives\n")
    basis_file.write("! basisfunctions_2D_1p includes the st cross-derivative too\n")
    basis_file.write("! basisfunctions_2D_2 includes the ss and tt derivatives additionally\n")
    basis_file.write("interface basisfunctions\n")
    basis_file.write("  module procedure basisfunctions_2D_0, basisfunctions_2D_1, basisfunctions_2D_1p, basisfunctions_2D_2\n")
    basis_file.write("end interface basisfunctions\n")
    basis_file.write("\n")
    basis_file.write("!> Two dimensional basisfunction with derivatives of order n\n")
    basis_file.write("!> and transposed matrix for better vectorisation\n")
    basis_file.write("!> basisfunctions_2D_1_T: first order derivatives in s and t\n")
    basis_file.write("!> basisfunctions_2D_2_T: first and second order derivatives in s and t\n")
    basis_file.write("interface basisfunctions_T\n")
    basis_file.write("  module procedure basisfunctions_2D_1_T, basisfunctions_2D_2_T\n")
    basis_file.write("end interface basisfunctions_T\n")
    basis_file.write("\n")
    basis_file.write("contains\n")
    basis_file.write("\n")
    
    basis_file.close()


# --- Close the Fortran module
def write_basis_function_tail():
    basis_file = open(filename,'a')
    basis_file.write("\n\n\nend module\n\n\n")
    basis_file.close()




# --- basisfunctions_1D_0
def write_basis_function_basisfunctions_1D_0(F_basis):
    
    basis_file = open(filename,'a')
    basis_file.write("!> Subroutine which defines the basis functions in one dimension with no derivatives\n")
    basis_file.write("!> IMPORTANT: auto-generated code for n_order=%d\n" % (n_order))
    basis_file.write("pure subroutine basisfunctions_1D_0(s,H)\n")
    basis_file.write("  real*8, intent(in)  :: s                      !< s-coordinate in the element (in [0,1])\n")
    basis_file.write("  real*8, intent(out) :: H(2,(n_order+1)/2)     !< Basis functions\n")
    basis_file.write("\n")
    basis_file.close()
    
    write_basis_function_1D(F_basis, 0, 0)
    
    basis_file = open(filename,'a')
    basis_file.write("end subroutine basisfunctions_1D_0\n\n")
    basis_file.close()



# --- basisfunctions_1D_1
def write_basis_function_basisfunctions_1D_1(F_basis):
    
    basis_file = open(filename,'a')
    basis_file.write("!> Subroutine which defines the basis functions in one dimension with first derivatives\n")
    basis_file.write("!> IMPORTANT: auto-generated code for n_order=%d\n" % (n_order))
    basis_file.write("pure subroutine basisfunctions_1D_1(s,H,H_s)\n")
    basis_file.write("  real*8, intent(in)  :: s                      !< s-coordinate in the element (in [0,1])\n")
    basis_file.write("  real*8, intent(out) :: H  (2,(n_order+1)/2)   !< Basis functions\n")
    basis_file.write("  real*8, intent(out) :: H_s(2,(n_order+1)/2)   !< Basis functions derived with respect to s\n")
    basis_file.write("\n")
    basis_file.write("  call basisfunctions_1D_0(s,H)\n")
    basis_file.write("\n")
    basis_file.close()
    
    write_basis_function_1D(F_basis, 1, 1)
    
    basis_file = open(filename,'a')
    basis_file.write("end subroutine basisfunctions_1D_1\n\n")
    basis_file.close()



# --- basisfunctions_1D_2
def write_basis_function_basisfunctions_1D_2(F_basis):
    
    basis_file = open(filename,'a')
    basis_file.write("!> Subroutine which defines the basis functions in one dimension with first and second derivatives\n")
    basis_file.write("!> IMPORTANT: auto-generated code for n_order=%d\n" % (n_order))
    basis_file.write("pure subroutine basisfunctions_1D_2(s,H,H_s,H_ss)\n")
    basis_file.write("  real*8, intent(in)  :: s                         !< s-coordinate in the element (in [0,1])\n")
    basis_file.write("  real*8, intent(out) :: H   (2,(n_order+1)/2)  !< Basis functions\n")
    basis_file.write("  real*8, intent(out) :: H_s (2,(n_order+1)/2)  !< Basis functions derived with respect to s\n")
    basis_file.write("  real*8, intent(out) :: H_ss(2,(n_order+1)/2)  !< Basis functions derived two times with respect to s\n")
    basis_file.write("\n")
    basis_file.write("  call basisfunctions_1D_1(s,H,H_s)\n")
    basis_file.write("\n")
    basis_file.close()
    
    write_basis_function_1D(F_basis, 2, 2)
    
    basis_file = open(filename,'a')
    basis_file.write("end subroutine basisfunctions_1D_2\n\n")
    basis_file.close()



# --- basisfunctions_2D_0
def write_basis_function_basisfunctions_2D_0(F_basis):
    
    basis_file = open(filename,'a')
    basis_file.write("!> Basisfunctions in 2D, value only.\n")
    basis_file.write("!> IMPORTANT: auto-generated code for n_order=%d\n" % (n_order))
    basis_file.write("pure subroutine basisfunctions_2D_0(s, t, H)\n")
    basis_file.write("  implicit none\n")
    basis_file.write("  real*8, intent(in)  :: s                  !< s-coordinate in the element [0,1]\n")
    basis_file.write("  real*8, intent(in)  :: t                  !< t-coordinate in the element [0,1]\n")
    basis_file.write("  real*8, intent(out) :: H(4,n_degrees)     !< Basis functions\n")
    basis_file.write("\n")
    basis_file.close()
    
    write_basis_function(F_basis, 0, 0, False, False)
    
    basis_file = open(filename,'a')
    basis_file.write("end subroutine basisfunctions_2D_0\n\n")
    basis_file.close()



# --- basisfunctions_2D_1
def write_basis_function_basisfunctions_2D_1(F_basis):
    
    basis_file = open(filename,'a')
    basis_file.write("!> Basisfunctions in 2D, with 1st derivatives.\n")
    basis_file.write("!> IMPORTANT: auto-generated code for n_order=%d\n" % (n_order))
    basis_file.write("pure subroutine basisfunctions_2D_1(s, t, H, H_s, H_t)\n")
    basis_file.write("  implicit none\n")
    basis_file.write("  real*8, intent(in)  :: s                  !< s-coordinate in the element [0,1]\n")
    basis_file.write("  real*8, intent(in)  :: t                  !< t-coordinate in the element [0,1]\n")
    basis_file.write("  real*8, intent(out) :: H  (4,n_degrees)   !< Basis functions\n")
    basis_file.write("  real*8, intent(out) :: H_s(4,n_degrees)   !< Basis functions derived with respect to s\n")
    basis_file.write("  real*8, intent(out) :: H_t(4,n_degrees)   !< Basis functions derived with respect to t\n")
    basis_file.write("\n")
    basis_file.write("  call basisfunctions_2D_0(s, t, H)\n")
    basis_file.write("\n")
    basis_file.close()
    
    write_basis_function(F_basis, 1, 1, False, False)
    
    basis_file = open(filename,'a')
    basis_file.write("end subroutine basisfunctions_2D_1\n\n")
    basis_file.close()



# --- basisfunctions_2D_1_T
def write_basis_function_basisfunctions_2D_1_T(F_basis):
    
    basis_file = open(filename,'a')
    basis_file.write("!> Transposed from the normal usage for easier vector operations (i.e.  basisfunctions_2D_1.T)\n")
    basis_file.write("!> IMPORTANT: auto-generated code for n_order=%d\n" % (n_order))
    basis_file.write("pure subroutine basisfunctions_2D_1_T(s, t, H, H_s, H_t)\n")
    basis_file.write("  implicit none\n")
    basis_file.write("  real*8, intent(in)  :: s                  !< s-coordinate in the element [0,1]\n")
    basis_file.write("  real*8, intent(in)  :: t                  !< t-coordinate in the element [0,1]\n")
    basis_file.write("  real*8, intent(out) :: H  (n_degrees,4)   !< Basis functions\n")
    basis_file.write("  real*8, intent(out) :: H_s(n_degrees,4)   !< Basis functions derived with respect to s\n")
    basis_file.write("  real*8, intent(out) :: H_t(n_degrees,4)   !< Basis functions derived with respect to t\n")
    basis_file.write("\n")
    basis_file.close()
    
    write_basis_function(F_basis, 0, 1, True, False)
    
    basis_file = open(filename,'a')
    basis_file.write("end subroutine basisfunctions_2D_1_T\n\n")
    basis_file.close()




# --- basisfunctions_2D_1p
def write_basis_function_basisfunctions_2D_1p(F_basis):
    
    basis_file = open(filename,'a')
    basis_file.write("!> Basisfunctions in 2D with first and cross-derivative\n")
    basis_file.write("!> IMPORTANT: auto-generated code for n_order=%d\n" % (n_order))
    basis_file.write("pure subroutine basisfunctions_2D_1p(s,t,H,H_s,H_t,H_st)\n")
    basis_file.write("  implicit none\n")
    basis_file.write("  real*8, intent(in)  :: s                   !< s-coordinate in the element [0,1]\n")
    basis_file.write("  real*8, intent(in)  :: t                   !< t-coordinate in the element [0,1]\n")
    basis_file.write("  real*8, intent(out) :: H   (4,n_degrees)   !< Basis functions\n")
    basis_file.write("  real*8, intent(out) :: H_s (4,n_degrees)   !< Basis functions derived with respect to s\n")
    basis_file.write("  real*8, intent(out) :: H_t (4,n_degrees)   !< Basis functions derived with respect to t\n")
    basis_file.write("  real*8, intent(out) :: H_st(4,n_degrees)   !< Basis functions derived with respect to s and t\n")
    basis_file.write("\n")
    basis_file.write("  call basisfunctions_2D_1(s, t, H, H_s, H_t)\n")
    basis_file.write("\n")
    basis_file.close()
    
    write_basis_function(F_basis, 1, 1, False, True)
    
    basis_file = open(filename,'a')
    basis_file.write("end subroutine basisfunctions_2D_1p\n\n")
    basis_file.close()




# --- basisfunctions_2D_2
def write_basis_function_basisfunctions_2D_2(F_basis):
    
    basis_file = open(filename,'a')
    basis_file.write("!> Basisfunctions with second derivatives in 2D\n")
    basis_file.write("!> IMPORTANT: auto-generated code for n_order=%d\n" % (n_order))
    basis_file.write("pure subroutine basisfunctions_2D_2(s, t, H, H_s, H_t, H_st, H_ss, H_tt)\n")
    basis_file.write("  implicit none\n")
    basis_file.write("  real*8, intent(in)  :: s                  !< s-coordinate in the element [0,1]\n")
    basis_file.write("  real*8, intent(in)  :: t                  !< t-coordinate in the element [0,1]\n")
    basis_file.write("  real*8, intent(out) :: H   (4,n_degrees)  !< Basis functions\n")
    basis_file.write("  real*8, intent(out) :: H_s (4,n_degrees)  !< Basis functions derived with respect to s\n")
    basis_file.write("  real*8, intent(out) :: H_t (4,n_degrees)  !< Basis functions derived with respect to t\n")
    basis_file.write("  real*8, intent(out) :: H_st(4,n_degrees)  !< Basis functions derived with respect to s and t\n")
    basis_file.write("  real*8, intent(out) :: H_ss(4,n_degrees)  !< Basis functions derived two times with respect to s\n")
    basis_file.write("  real*8, intent(out) :: H_tt(4,n_degrees)  !< Basis functions derived two times with respect to t\n")
    basis_file.write("\n")
    basis_file.write("  call basisfunctions_2D_1(s, t, H, H_s, H_t)\n")
    basis_file.write("\n")
    basis_file.close()
    
    write_basis_function(F_basis, 1, 1, False, True)
    write_basis_function(F_basis, 2, 2, False, False)
    
    basis_file = open(filename,'a')
    basis_file.write("end subroutine basisfunctions_2D_2\n\n")
    basis_file.close()



# --- basisfunctions_2D_2_T
def write_basis_function_basisfunctions_2D_2_T(F_basis):
    
    basis_file = open(filename,'a')
    basis_file.write("!> Basisfunctions with second derivatives in 2D and transposed matrix for better vectorisation\n")
    basis_file.write("!> IMPORTANT: auto-generated code for n_order=%d\n" % (n_order))


    basis_file.write("pure subroutine basisfunctions_2D_2_T(s, t, H, H_s, H_t, H_st, H_ss, H_tt)\n")
    basis_file.write("  implicit none\n")
    basis_file.write("  real*8, intent(in)  :: s                  !< s-coordinate in the element [0,1]\n")
    basis_file.write("  real*8, intent(in)  :: t                  !< t-coordinate in the element [0,1]\n")
    basis_file.write("  real*8, intent(out) :: H   (n_degrees,4)  !< Basis functions\n")
    basis_file.write("  real*8, intent(out) :: H_s (n_degrees,4)  !< Basis functions derived with respect to s\n")
    basis_file.write("  real*8, intent(out) :: H_t (n_degrees,4)  !< Basis functions derived with respect to t\n")
    basis_file.write("  real*8, intent(out) :: H_st(n_degrees,4)  !< Basis functions derived with respect to s and t\n")
    basis_file.write("  real*8, intent(out) :: H_ss(n_degrees,4)  !< Basis functions derived two times with respect to s\n")
    basis_file.write("  real*8, intent(out) :: H_tt(n_degrees,4)  !< Basis functions derived two times with respect to t\n")
    basis_file.write("\n")
    basis_file.write("  call basisfunctions_2D_1_T(s, t, H, H_s, H_t)\n")
    basis_file.write("\n")
    basis_file.close()
    
    write_basis_function(F_basis, 1, 1, True, True)
    write_basis_function(F_basis, 2, 2, True, False)
    
    basis_file = open(filename,'a')
    basis_file.write("end subroutine basisfunctions_2D_2_T\n\n")
    basis_file.close()



# --- basisfunctions3
def write_basis_function_basisfunctions3(F_basis):
    
    basis_file = open(filename,'a')
    
    basis_file.write("!> IMPORTANT: auto-generated code for n_order=%d\n" % (n_order))
    basis_file.write("pure subroutine basisfunctions3(s, t, H, H_s, H_t)\n")
    basis_file.write("  implicit none\n")
    basis_file.write("  real*8, intent(in)  :: s                  !< s-coordinate in the element\n")
    basis_file.write("  real*8, intent(in)  :: t                  !< t-coordinate in the element\n")
    basis_file.write("  real*8, intent(out) :: H  (4,n_degrees)   !< Basis functions\n")
    basis_file.write("  real*8, intent(out) :: H_s(4,n_degrees)   !< Basis functions derived with respect to s\n")
    basis_file.write("  real*8, intent(out) :: H_t(4,n_degrees)   !< Basis functions derived with respect to t\n")
    
    basis_file.close()
    
    write_basis_function(F_basis, 0, 1, False, False)
    
    basis_file = open(filename,'a')
    basis_file.write("end subroutine basisfunctions3\n\n")
    basis_file.close()
    




# --- gauss.f90
def write_gauss_points():
    
    n_gauss = n_order+1
    
    gauss_numpy = numpy.polynomial.legendre.leggauss(n_gauss)
    xgauss = gauss_numpy[0]
    wgauss = gauss_numpy[1]
    # --- Rescale because numpy is in interval [-1,+1], and JOREK is in [0,1]
    for i in range(n_gauss):
        xgauss[i] = 0.5 + 0.5 * xgauss[i]
        wgauss[i] = 0.5 * wgauss[i]
    
    gauss_file = open('gauss.f90','w')
    
    gauss_file.write("!> IMPORTANT: auto-generated code for n_order=%d\n" % (n_order))
    gauss_file.write("!> Contains positions (xgauss) and weights (wgauss) of\n")
    gauss_file.write("!! Gaussian points for Gaussian integration.\n")
    gauss_file.write("!!\n")
    gauss_file.write("!! The values are valid for normalised coordinates in the range [0,1]\n")
    gauss_file.write("!!\n")
    gauss_file.write("!! Taken from numpy.polynomial.legendre.leggauss function\n")
    gauss_file.write("!! converted to [0,1], with weights normalized to sum 1\n")
    gauss_file.write("module gauss\n")
    gauss_file.write("\n")
    gauss_file.write("  integer, parameter :: n_gauss   = %d    !< Number of Gaussian points\n" % (n_gauss))
    gauss_file.write("\n")
    gauss_file.write("  real*8,  parameter :: Xgauss(n_gauss) = &\n")
    gauss_file.write("                        [ &\n")
    for i in range(n_gauss):
        gauss_file.write("                         %e" % (xgauss[i]))
        if (i < n_gauss - 1):
            gauss_file.write(", &\n")
        else:
            gauss_file.write(" &\n")
    gauss_file.write("                        ]\n")
    gauss_file.write("\n")
    gauss_file.write("  real*8,  parameter :: Wgauss(n_gauss) = &\n")
    gauss_file.write("                        [ &\n")
    for i in range(n_gauss):
        gauss_file.write("                         %e" % (wgauss[i]))
        if (i < n_gauss - 1):
            gauss_file.write(", &\n")
        else:
            gauss_file.write(" &\n")
    gauss_file.write("                        ]\n")
    gauss_file.write("\n")
    gauss_file.write("integer, parameter :: n_gauss_2 = n_gauss * n_gauss  !< Square of n_gauss\n")
    gauss_file.write("\n")
    gauss_file.write("end module gauss\n")
    
    gauss_file.close()
    










# ------------------------------------------------------
# ------------------------------------------------------
# ------- Section: Add new term to control point
# ------------------------------------------------------
# ------------------------------------------------------




# --- Sympy is actually very bad with long complicated terms
# --- Therefore, it is safer to separate each component into array entries
# --- But this is a bit painful, you need to loop over all components...
def add_terms_to_array(P_ij,my_coef,P_kl,vec,P_ref):
    # --- First the P_ref
    for ic in range(len(P_kl)):
        if (str(P_ref) in str(P_kl[ic])):
            none_found = True
            for ix in range(len(P_ij)):
                if (str(P_ref) in str(P_ij[ix])):
                    P_ij[ix] = P_ij[ix] + my_coef*P_kl[ic]
                    none_found = False
            if (none_found):
                P_ij.append(my_coef*P_kl[ic])
    # --- Check each vector
    for i in range(n_vectors):
        for j in range(n_vectors):
            for ic in range(len(P_kl)):
                if (str(vec[i][j]) in str(P_kl[ic])):
                    none_found = True
                    for ix in range(len(P_ij)):
                        if (str(vec[i][j]) in str(P_ij[ix])):
                            P_ij[ix] = P_ij[ix] + my_coef*P_kl[ic]
                            none_found = False
                    if (none_found):
                        P_ij.append(my_coef*P_kl[ic])
    return P_ij








# ------------------------------------------------------
# ------------------------------------------------------
# ------- Section: the main function
# ------------------------------------------------------
# ------------------------------------------------------



# --- Generate basis functions automatically
def main():

    # -----------------------------------------------------------------
    # ------- Sub-section: Define basis vectors and control points
    # -----------------------------------------------------------------
    
    # --- We need a reference node for any control point
    # --- For example, the reference node of P_11 is P_00
    # --- While the reference node of P_22 (in cubic) is P_33
    # --- And the reference node of P_21 is P_30
    i_refnode  = [[0 for j in range(n_order+1)] for i in range(n_order+1)] # nodes index from 1 to 4
    j_refnode  = [[0 for j in range(n_order+1)] for i in range(n_order+1)] # nodes index from 1 to 4
    ij_ref = ['00',str(n_order)+'0',str(n_order)+str(n_order),'0'+str(n_order)]
    # --- We also need to know which control point we mean when we loop on n_vectors at each node
    node_ref   = [[0 for j in range(n_order+1)] for i in range(n_order+1)] # nodes index from 1 to 4
    i_ctpt = [ [ [0 for j in range(n_vectors)] \
               for i in range(n_vectors)] \
             for k in range(n_vertex)]
    j_ctpt = [ [ [0 for j in range(n_vectors)] \
               for i in range(n_vectors)] \
             for k in range(n_vertex)]
    i_vect = [ [0 for j in range(n_order+1)] \
               for i in range(n_order+1)]
    j_vect = [ [0 for j in range(n_order+1)] \
               for i in range(n_order+1)]
    for i in range(n_order+1):
        for j in range(n_order+1):
            if ( (i < n_vectors) and (j < n_vectors) ):
                node_ref[i][j]  = 0
                i_refnode[i][j] = 0
                j_refnode[i][j] = 0
                i_vec = i
                j_vec = j
            if ( (i >= n_vectors) and (j < n_vectors) ):
                node_ref[i][j]  = 1
                i_refnode[i][j] = n_order
                j_refnode[i][j] = 0
                i_vec = n_order - i
                j_vec = j
            if ( (i < n_vectors) and (j >= n_vectors) ):
                node_ref[i][j]  = 3
                i_refnode[i][j] = 0
                j_refnode[i][j] = n_order
                i_vec = i
                j_vec = n_order - j
            if ( (i >= n_vectors) and (j >= n_vectors) ):
                node_ref[i][j]  = 2
                i_refnode[i][j] = n_order
                j_refnode[i][j] = n_order
                i_vec = n_order - i
                j_vec = n_order - j
            i_vect[i][j] = i_vec
            j_vect[i][j] = j_vec
            i_ctpt[ node_ref[i][j] ][ i_vec ][ j_vec ] = i
            j_ctpt[ node_ref[i][j] ][ i_vec ][ j_vec ] = j
    
    
    # --- Bezier control points P(n_order+1,n_order+1)
    P = [[[] for j in range(n_order+1)] for i in range(n_order+1)]
    P_ref = [[sym.Symbol('') for j in range(n_order+1)] for i in range(n_order+1)]
    
    # --- Define vectors as u(k)_ij, with k the node index(4), i and j the x,y indices of vectors
    # --- eg. in cubic elements, at node P_00, we have 
    # ---     u = u(0)_10
    # ---     v = u(0)_01
    # ---     w = u(0)_11
    # --- note: this also defines an unused vector u(0)_00
    u = [ [ [sym.Symbol('u('+str(k)+')_'+str(i)+str(j)) \
            for j in range(n_vectors)] \
          for i in range(n_vectors)] \
        for k in range(n_vertex)]
    
    # --- Special case if we want the same vector names as in the paper
    # --- ie. u,v,w,i,j,m,n,k
    if (vectors_like_paper):
        if (n_order > 5):
            print("warning: vectors_like_paper only valid for n_order<=5")
        else:
            for k in range(n_vertex):
                u[k][1][0] = sym.Symbol('u_'+ij_ref[k])
                u[k][0][1] = sym.Symbol('v_'+ij_ref[k])
                u[k][1][1] = sym.Symbol('w_'+ij_ref[k])
                if (n_order == 5):
                    u[k][2][0] = sym.Symbol('i_'+ij_ref[k])
                    u[k][0][2] = sym.Symbol('j_'+ij_ref[k])
                    u[k][2][1] = sym.Symbol('m_'+ij_ref[k])
                    u[k][1][2] = sym.Symbol('n_'+ij_ref[k])
                    u[k][2][2] = sym.Symbol('k_'+ij_ref[k])
    
    # --- Define reference nodes
    for i in range(n_order+1):
        for j in range(n_order+1):
            i_ref  = i_refnode[i][j]
            j_ref  = j_refnode[i][j]
            P_ref[i][j] = sym.Symbol('P_'+str(i_ref)+str(j_ref))
    # --- Define vertex nodes plus u and v vectors
    # --- General formulation:
    #     P_ij = h_ij * u_ij &
    #            + \sum_k^i \sum_l^j (-1)^{1+i+j+k+l}  (1 - \delta_ki \delta_lj)  \binom_k^i \binom_l^j  P_kl
    # --- Note, the (1-delta) function means that the RHS is a function of all P_kl with k<=i and l<=j, except
    # --- for P_ij itself. Therefore, we must work our way up the chain of control points, each time.
    # --- Reference point
    for kn in range(n_vertex):
        # --- Reference point
        i = 0 ; j = 0
        ii = i_ctpt[kn][i][j]
        jj = j_ctpt[kn][i][j]
        P[ii][jj].append(P_ref[ii][jj])
        for Iter in range(n_vectors):
            if (Iter != 0):
                i = Iter
                j = Iter
                ii = i_ctpt[kn][i][j]
                jj = j_ctpt[kn][i][j]
                # --- The vector itself
                P[ii][jj].append(u[kn][i][j])
                # --- The previous points
                for k in range(i+1):
                    for l in range(j+1):
                        if ( (k==i) and (l==j) ): continue
                        my_coef = (-1)**(1+i+j+k+l) * b_coef(k,i) * b_coef(l,j)
                        kk = i_ctpt[kn][k][l]
                        ll = j_ctpt[kn][k][l]
                        # --- Add each component separately (safer for Sympy with large terms)
                        P[ii][jj] = add_terms_to_array(P[ii][jj],my_coef,P[kk][ll],u[kn],P_ref[ii][jj])
            j = Iter
            for i in range(Iter+1,n_vectors):
                if ( (i==0) and (j==0) ): continue
                ii = i_ctpt[kn][i][j]
                jj = j_ctpt[kn][i][j]
                # --- The vector itself
                P[ii][jj].append(u[kn][i][j])
                # --- The previous points
                for k in range(i+1):
                    for l in range(j+1):
                        if ( (k==i) and (l==j) ): continue
                        my_coef = (-1)**(1+i+j+k+l) * b_coef(k,i) * b_coef(l,j)
                        kk = i_ctpt[kn][k][l]
                        ll = j_ctpt[kn][k][l]
                        # --- Add each component separately (safer for Sympy with large terms)
                        P[ii][jj] = add_terms_to_array(P[ii][jj],my_coef,P[kk][ll],u[kn],P_ref[ii][jj])
            i = Iter
            for j in range(Iter+1,n_vectors):
                if ( (i==0) and (j==0) ): continue
                ii = i_ctpt[kn][i][j]
                jj = j_ctpt[kn][i][j]
                # --- The vector itself
                P[ii][jj].append(u[kn][i][j])
                # --- The previous points
                for k in range(i+1):
                    for l in range(j+1):
                        if ( (k==i) and (l==j) ): continue
                        my_coef = (-1)**(1+i+j+k+l) * b_coef(k,i) * b_coef(l,j)
                        kk = i_ctpt[kn][k][l]
                        ll = j_ctpt[kn][k][l]
                        # --- Add each component separately (safer for Sympy with large terms)
                        P[ii][jj] = add_terms_to_array(P[ii][jj],my_coef,P[kk][ll],u[kn],P_ref[ii][jj])
    
    # --- Clean-up zeros
    for k in range(n_vertex):
        for i in range(n_vectors):
            for j in range(n_vectors):
                ii = i_ctpt[k][i][j]
                jj = j_ctpt[k][i][j]
                for ijk in range(len(P[ii][jj])-1,-1,-1):
                    if (P[ii][jj][ijk] == 0*sym.Symbol('')):
                        P[ii][jj].pop(ijk)
                for ijk in range(len(P[ii][jj])-1,-1,-1):
                    none_found = True
                    if (str(P_ref[ii][jj]) in str(P[ii][jj][ijk])):
                        none_found = False
                    for iv in range(n_vectors):
                        for jv in range(n_vectors):
                            if (str(u[k][iv][jv]) in str(P[ii][jj][ijk])):
                                none_found = False
                                break
                    if (none_found):
                        P[ii][jj].pop(ijk)
    
    # --- Print the nodal formulation?
    if (True):
        for k in range(n_vertex):
            if (k > 0): continue
            for i in range(n_vectors):
                for j in range(n_vectors):
                    ii = i_ctpt[k][i][j]
                    jj = j_ctpt[k][i][j]
                    print("%s%s" % ('P_'+str(ii)+str(jj)+' = ',P[ii][jj]) )
    
    
    # -----------------------------------------------------------------
    # ------- Sub-section: Define basis functions for 2D element
    # -----------------------------------------------------------------
    
    
    # --- Declare array for each vector component
    FF = {}
    for k in range(n_vertex):
        for i in range(n_vectors):
            for j in range(n_vectors):
                ii = i_ctpt[k][i][j]
                jj = j_ctpt[k][i][j]
                if not (str(P_ref[ii][jj]) in FF.keys()):
                    FF[str(P_ref[ii][jj])] = []
                if ( (i==0) and (j==0) ): continue # avoid irrelevant vectors u(0)_00
                if not (str(u[k][i][j]) in FF.keys()):
                    FF[str(u[k][i][j])] = []
    
    # --- Build array for each vector component
    s = sym.Symbol('s')
    t = sym.Symbol('t')
    for i in range(n_order+1):
        for j in range(n_order+1):
            F = bezier_coef(s, t, i, j, n_order)
            k = node_ref[i][j]
            for ijk in range(len(P[i][j])):
                if (str(P_ref[i][j]) in str(P[i][j][ijk])):
                    FF[str(P_ref[i][j])].append(F * P[i][j][ijk])
                # --- Check all vectors on that node
                for i_vec in range(n_vectors):
                    for j_vec in range(n_vectors):
                        if ( (i_vec==0) and (j_vec==0) ): continue # avoid irrelevant vectors u(0)_00
                        if (str(u[k][i_vec][j_vec]) in str(P[i][j][ijk])):
                            FF[str(u[k][i_vec][j_vec])].append(F * P[i][j][ijk])
    
    # --- Create the basis function for JOREK
    # --- with the derivatives, we need (n_order-2)**2 of them
    F_basis = [[ [[ sym.Symbol("") for m in range(n_deriv)] for n in range(n_deriv)] \
                for i in range(n_degrees)] for k in range(n_vertex)] # the nodes and degrees
    # --- The node contributions
    for d_s in range(n_deriv):
        for d_t in range(n_deriv):
            if (d_s + d_t > n_deriv-1): continue # only up to _ss, _tt, _st
            for i_node in range(n_vertex):
                if (i_node == 0):
                    i = 0 ; j = 0
                if (i_node == 1):
                    i = 0 ; j = n_order
                if (i_node == 2):
                    i = n_order ; j = n_order
                if (i_node == 3):
                    i = n_order ; j = 0
                k = node_ref[i][j]
                for mno in range(len(FF[str(P_ref[i][j])])):
                    deriv = FF[str(P_ref[i][j])][mno] / P_ref[i][j]
                    for dd in range(d_s):
                        deriv = sym.diff(deriv,s)
                    for dd in range(d_t):
                        deriv = sym.diff(deriv,t)
                    F_basis[k][0][d_s][d_t] = F_basis[k][0][d_s][d_t] + deriv
    # --- The vector contributions
    for d_s in range(n_deriv):
        for d_t in range(n_deriv):
            if (d_s + d_t > n_deriv-1): continue # only up to _ss, _tt, _st
            for k in range(n_vertex):
                for i in range(n_vectors):
                    for j in range(n_vectors):
                        if ( (i==0) and (j==0) ): continue # avoid irrelevant vectors u(0)_00
                        ii = i_ctpt[k][i][j]
                        jj = j_ctpt[k][i][j]
                        i_deg = basis_index[i][j]
                        if (str(u[k][i][j]) in FF.keys()):
                            for mno in range(len(FF[str(u[k][i][j])])):
                                deriv = FF[str(u[k][i][j])][mno] / u[k][i][j]
                                for dd in range(d_s):
                                    deriv = sym.diff(deriv,s)
                                for dd in range(d_t):
                                    deriv = sym.diff(deriv,t)
                                F_basis[k][i_deg][d_s][d_t] = F_basis[k][i_deg][d_s][d_t] + deriv
    
    
    
    # -----------------------------------------------------------------
    # ------- Sub-section: Define basis functions for 1D curve
    # -----------------------------------------------------------------
    
    
    # --- Declare array for each vector component in 1D
    FF_1D = {}
    for k in range(n_vertex):
        for i in range(n_vectors):
            j = 0
            ii = i_ctpt[k][i][j]
            jj = j_ctpt[k][i][j]
            if not (str(P_ref[ii][jj]) in FF_1D.keys()):
                FF_1D[str(P_ref[ii][jj])] = []
            if ( (i==0) and (j==0) ): continue # avoid irrelevant vectors u(0)_00
            if not (str(u[k][i][j]) in FF_1D.keys()):
                FF_1D[str(u[k][i][j])] = []
    
    # --- Build array for each vector component in 1D
    s = sym.Symbol('s')
    t = 0
    for i in range(n_order+1):
        j = 0
        F = bezier_coef(s, t, i, j, n_order)
        k = node_ref[i][j]
        for ijk in range(len(P[i][j])):
            if (str(P_ref[i][j]) in str(P[i][j][ijk])):
                FF_1D[str(P_ref[i][j])].append(F * P[i][j][ijk])
            # --- Check all vectors on that node
            for i_vec in range(n_vectors):
                for j_vec in range(n_vectors):
                    if ( (i_vec==0) and (j_vec==0) ): continue # avoid irrelevant vectors u(0)_00
                    if (str(u[k][i_vec][j_vec]) in str(P[i][j][ijk])):
                        FF_1D[str(u[k][i_vec][j_vec])].append(F * P[i][j][ijk])
    
    # --- Create the basis function for JOREK
    # --- with the derivatives, we need (n_order-2)**2 of them
    F_basis_1D = [[ [ sym.Symbol("") for m in range(n_deriv)] \
                for i in range(n_degrees)] for k in range(int(n_vertex/2))] # the nodes and degrees
    # --- The node contributions
    for d_s in range(n_deriv):
        d_t = 0
        if (d_s + d_t > n_deriv-1): continue # only up to _ss, _tt, _st
        for i_node in range(int(n_vertex/2)):
            if (i_node == 0):
                i = 0 ; j = 0
            if (i_node == 1):
                i = n_order ; j = 0
            k = node_ref[i][j]
            for mno in range(len(FF_1D[str(P_ref[i][j])])):
                deriv = FF_1D[str(P_ref[i][j])][mno] / P_ref[i][j]
                for dd in range(d_s):
                    deriv = sym.diff(deriv,s)
                for dd in range(d_t):
                    deriv = sym.diff(deriv,t)
                F_basis_1D[k][0][d_s] = F_basis_1D[k][0][d_s] + deriv
    # --- The vector contributions
    for d_s in range(n_deriv):
        d_t = 0
        if (d_s + d_t > n_deriv-1): continue # only up to _ss, _tt, _st
        for k in range(int(n_vertex/2)):
            for i in range(n_vectors):
                j = 0
                if ( (i==0) and (j==0) ): continue # avoid irrelevant vectors u(0)_00
                ii = i_ctpt[k][i][j]
                jj = j_ctpt[k][i][j]
                i_deg = i
                if (str(u[k][i][j]) in FF_1D.keys()):
                    for mno in range(len(FF_1D[str(u[k][i][j])])):
                        deriv = FF_1D[str(u[k][i][j])][mno] / u[k][i][j]
                        for dd in range(d_s):
                            deriv = sym.diff(deriv,s)
                        for dd in range(d_t):
                            deriv = sym.diff(deriv,t)
                        F_basis_1D[k][i_deg][d_s] = F_basis_1D[k][i_deg][d_s] + deriv
    
    
    
    # -----------------------------------------------------------------
    # ------- Sub-section: Print basis_function Fortran module
    # -----------------------------------------------------------------
    
    
    if (python_formalism):
        basis_file = open(filename,'w')
        basis_file.write("\n")
        write_basis_function(F_basis, 0, 0, False, False)
        basis_file.close()
    else:
        write_basis_function_headers()
        write_basis_function_basisfunctions_1D_0(F_basis_1D)
        write_basis_function_basisfunctions_1D_1(F_basis_1D)
        write_basis_function_basisfunctions_1D_2(F_basis_1D)
        write_basis_function_basisfunctions_2D_0(F_basis)
        write_basis_function_basisfunctions_2D_1(F_basis)
        write_basis_function_basisfunctions_2D_1_T(F_basis)
        write_basis_function_basisfunctions_2D_1p(F_basis)
        write_basis_function_basisfunctions_2D_2(F_basis)
        write_basis_function_basisfunctions_2D_2_T(F_basis)
        write_basis_function_basisfunctions3(F_basis)
        write_basis_function_tail()
        #if (n_order > 7):
        #    write_gauss_points()
        write_gauss_points()
    








# --- Execute as a script
main()


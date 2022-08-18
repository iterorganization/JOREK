"""
jorek_read_h5.py

Read JOREK hdf5 restart files and define functions for interpolation
on the exported fields.

Before reading, set:
    - Variable numbers
Before interpolating, set:
    - Without_n0_mode (default false)

Created by Daan van Vugt on 2017-01-18

prec controls the precision of the calculation. Values: np.float32 or np.float64
Calculation time is not really affected, render time maybe?
"""
import h5py
import numpy as np
prec=np.float32
vtk_prec=vtk.VTK_FLOAT



"""
Class encapsulating some read data from a restart file

input arguments:
    Variables: a list of all variables you want
    i_plane: which plane to interpolate at. If -1 make a 3D field.
"""
class fields(object):
    def read(self, filename, variables='', file_prev=None, interp_fraction=None):
        if (isinstance(variables, str)):
            self.vars = [int(x)-1 for x in variables.split(',')]
        else:
            raise "Error: expected a string of variables"

        with h5py.File(filename, 'r') as hf:
            self.n_var        = hf.get('n_var')[0]
            self.n_period     = hf.get('n_period')[0]
            self.n_tor        = hf.get('n_tor')[0]
            self.n_vertex_max = hf.get('n_vertex_max')[0]
            n_elements   = hf.get('n_elements')[0]

        # Assume the grid not to change. Important!
        if (not (hasattr(self, 'n_elements') and self.n_elements == n_elements)):
            self.vertex   = read_mmap_or_h5py(filename, 'vertex')
            self.x        = read_mmap_or_h5py(filename, 'x', type_out=prec)
            self.size     = read_mmap_or_h5py(filename, 'size', type_out=prec)
        self.n_elements = n_elements
        self.values   = read_mmap_or_h5py(filename, 'values', type_out=prec)

        if (file_prev is not None):
            with h5py.File(file_prev, 'r') as hf2:
                if (self.n_elements != hf2.get('n_elements')[0]):
                    raise "Error: Files with different numbers of elements read! Refinement is not supported"
            # Interpolate before making grid etc
            self.values = f*self.values + (1.0-f)*read_mmap_or_h5py(file_prev, 'values', type_out=prec)


    """
    Create VTK objects from points and connectivity matrix
    input:
        n_sub: number of subdivisions per element
        n_plane: number of planes in toroidal direction
        phi: range (2 elements) or single value of phi
        without_n0_mode: do not include n0 mode if true

    returns:
        vtkUnstructuredGrid
    """
    def to_vtk(self, n_sub=2, phi=[0,2*np.pi], n_plane=8, without_n0_mode=False, 
               force_remake_grid=False, quadratic=True, output=None):
        import vtk
        from vtk.util import numpy_support as npvtk

        # If we do not have a range in phi make only one plane
        periodic=False
        if (isinstance(phi,int)):
            n_plane = 1
            phis = np.asarray([phi])
        elif (n_plane == 1):
            phis = np.asarray([phi[0]])
        else:
            if (quadratic):
                n_plane = n_plane*2-1
            periodic = (np.mod(phi[0]-phi[1],2*np.pi) < 1e-9)
            phis = np.linspace(phi[0],phi[1],num=n_plane,endpoint=not periodic)

        if (quadratic):
            n_sub = n_sub*2+1

        if (output == None):
            output = vtk.vtkUnstructuredGrid()

        settings = [n_sub, phi, n_plane, without_n0_mode, quadratic, self.n_elements]
        if (not hasattr(self, 'old_settings') or self.old_settings != settings):
            force_remake_grid = True
        self.old_settings = settings

        if (force_remake_grid or not hasattr(self, 'points') or not hasattr(self, 'cells')):
            (xyz, ien) = create_grid(self.x, self.vertex, self.size, self.n_elements,
                                     n_sub, phis, n_plane, periodic, quadratic)

            pcoords = npvtk.numpy_to_vtk(xyz, deep=True, array_type=vtk_prec)
            self.points = vtk.vtkPoints()
            self.points.SetData(pcoords)

            self.cells = vtk.vtkCellArray()
            self.cells.SetCells(ien.shape[0], npvtk.numpy_to_vtk(ien, deep=True, array_type=vtk.VTK_ID_TYPE))

        output.SetPoints(self.points)

        HZ = toroidal_basis(self.n_tor, self.n_period, phis, without_n0_mode)
        val = interp_scalars_3D(self.values[self.vars,:,:,:],
                               self.vertex, self.size, n_sub, HZ).reshape((len(self.vars),-1))
        # Could split here if we run into memory problems and delete each part after use

        for i in range(len(self.vars)):
            tmp = npvtk.numpy_to_vtk(val[i,:], deep=True, array_type=vtk_prec)
            tmp.SetName(str(self.vars[i]))
            output.GetPointData().AddArray(tmp)


        if (quadratic):
            if (n_plane > 1):
                etype = vtk.VTK_QUADRATIC_HEXAHEDRON
            else: # or stay 2D
                etype = vtk.VTK_QUADRATIC_QUAD
        else:
            if (n_plane > 1):
                etype = vtk.VTK_HEXAHEDRON
            else: # or stay 2D
                etype = vtk.VTK_QUAD

        output.SetCells(etype, self.cells)
        return output


"""
Read part of a hdf5 file directly

memmap is roughly 20% faster on my system
"""
def read_mmap_or_h5py(path, h5path, type_out=None):
    with h5py.File(path, 'r') as f:
        ds = f[h5path]
        offset = ds.id.get_offset()
        if (ds.chunks is None and ds.compression is None and offset > 0):
            dtype = ds.dtype
            shape = ds.shape
            arr = np.memmap(path, mode='r', shape=shape, offset=offset, dtype=dtype)
            if (type_out is not None and type_out is not dtype):
                return arr.astype(type_out)
        else:
            if (type_out is None):
                arr = np.array(f.get(h5path))
            else:
                arr = np.array(f.get(h5path), dtype=type_out)
    return arr


"""
Calculate RZ positions of all points
return x[element,is,it,var] (where var is 0->R or 1->Z)
"""
def grid_2D(x, vertex, size, n_sub):
    # Calculate RZ for all of the elements (dimension 0) for each of the s
    # positions (dimension 1) for each of the t positions (dimension 2)
    # Multiply x[var, order, node[vertex, element]] on the last two dimensions
    # with size[order, vertex, element]*bf[order, vertex, s, t]
    # See http://stackoverflow.com/questions/26089893/understanding-numpys-einsum
    #return np.einsum('lijk,ijk,ijmn->kmnl', x[:,:,vertex], size, bf(n_sub))
    # Code below is ~5x faster or so! try again when einsum supports optimize=True

    # First create a temporary array holding: x[order, vertex, element, var]
    tmp = np.zeros((x.shape[0], x.shape[1],vertex.shape[0],vertex.shape[1]), dtype=prec)
    # Fill it with the right x
    for i in range(vertex.shape[0]): # small loop over vertices (hardcode 4 here?)
        tmp[:,:,i,:] = x[:,:,vertex[i,:]-1]
    # multiply by size[order, vertex, element]
    tmp[0,:,:,:] *= size
    tmp[1,:,:,:] *= size
    # Create output array
    out = np.zeros((vertex.shape[1],n_sub+2,n_sub+2,2), dtype=prec)
    out[:,:,:,0] = np.tensordot(tmp[0,:,:,:], bf(n_sub), axes=((0,1),(0,1)))
    out[:,:,:,1] = np.tensordot(tmp[1,:,:,:], bf(n_sub), axes=((0,1),(0,1)))
    return out


"""
Calculate toroidal basis functions
"""
def toroidal_basis(n_tor, n_period, phis, without_n0_mode):
    # Setup toroidal coefficients for each plane and toroidal harmonic
    HZ = np.zeros((n_tor,len(phis)))
    for i in range(n_tor):
        mode = np.floor((i+1)/2)*n_period
        if (i == 0):
            if (not without_n0_mode):
                HZ[i,:] = 1
        elif (i % 2 == 0):
            HZ[i,:] = np.sin(mode*phis)
        elif (i % 2 == 1):
            HZ[i,:] = np.cos(mode*phis)
    return HZ



"""
Interpolate scalars on 2D poloidal plane

returns:
    values: interpolated values, values[var, harmonic, element, is, it]
"""
def interp_scalars(values, vertex, size, n_sub):
    # Multiply values[var,order,harm,vertex,element] with
    # size[order, vertex, element] and bf[order, vertex, s, t]
    return np.einsum('lihjk,ijk,ijmn->lhkmn',
                    values[:,:,:,vertex-1],
                    size, bf(n_sub))

"""
Interpolate scalars on 2D planes * n_planes
"""
def interp_scalars_3D(values, vertex, size, n_sub, HZ):
    vals = interp_scalars(values, vertex, size, n_sub)
    return np.einsum('lhkmn,hp->lpkmn', vals, HZ)


"""
Create a grid of (nsub+2)**2 points per element, at phis positions
return points and connectivity matrix
"""
def create_grid(x, vertex, size, n_elements, n_sub, phis, n_plane, periodic, quadratic):
    RZ     = grid_2D(x, vertex, size, n_sub)

    n_p = n_sub + 2 # number of points per dimension in each element

    # Create connectivity data
    # Calculate 2D connectivity first
    # For each element, calculate the number of the lowest point
    # Create (n_p-1)**2 quadrangles
    n_points = n_elements*(n_p**2) # number of points in one plane
    n_cells  = n_elements*((n_p-1)**2) # Number of cells in one plane
    if (n_plane > 1): # Create a volume
        if (periodic):
            n_cells_tor = n_plane
        else:
            n_cells_tor = n_plane - 1
    else:
        n_cells_tor = 1

    xyz = np.zeros((n_points*n_plane,3))
    for i in range(n_plane):
        xyz[i*n_points:(i+1)*n_points,0] = np.ravel(RZ[:,:,:,0]*np.cos(phis[i]))
        xyz[i*n_points:(i+1)*n_points,1] = np.ravel(RZ[:,:,:,1])
        xyz[i*n_points:(i+1)*n_points,2] = np.ravel(RZ[:,:,:,0]*np.sin(phis[i]))

    # See http://www.vtk.org/doc/nightly/html/classvtkQuadraticHexahedron.html#details
    if (quadratic):
        if (n_plane > 1):
            # The base block for a quadratic element
            face_block = np.zeros(21, dtype=np.int32)
            face_block[1:21] = [0,2*n_p,2*n_p+2,2, # corners front face
                                2*n_points,2*n_points+2*n_p,2*n_points+2*n_p+2,2*n_points+2, # corners back face
                                n_p,2*n_p+1,n_p+2,1,# midedges front face
                                2*n_points+n_p,2*n_points+2*n_p+1,2*n_points+n_p+2,2*n_points+1,# midedges back face
                                n_points,2*n_p+n_points,2*n_p+2+n_points,2+n_points]# midedges middle
            i_start_t = np.arange(0,n_p*(n_p-1),2*n_p,dtype=np.int32)
            i_start_s = np.arange(0,n_p-1,2,dtype=np.int32)

            element_block = i_start_t[:,np.newaxis,np.newaxis] + \
                            i_start_s[np.newaxis,:,np.newaxis] + \
                            face_block[np.newaxis,np.newaxis,:]

            i_start_elm   = np.arange(0,n_points, n_p**2, dtype=np.int32)
            i_start_plane = np.arange(0,n_points*(n_plane-1),2*n_points, dtype=np.int32)
            if (periodic):
                i_start_plane[-1] = 0
            ien = np.reshape(i_start_plane[:,np.newaxis,np.newaxis,np.newaxis,np.newaxis]+
                             i_start_elm[np.newaxis,:,np.newaxis,np.newaxis,np.newaxis]+
                             element_block[np.newaxis,np.newaxis,:,:,:], (-1,21))
            ien[:,0] = 20
        else:
            # The base block for a quadratic element
            face_block = np.zeros(9, dtype=np.int32)
            face_block[1:9] = [0,2*n_p,2*n_p+2,2, # corners
                               n_p,2*n_p+1,n_p+2,1] # midedges
            i_start_t = np.arange(0,n_p*(n_p-1),2*n_p,dtype=np.int32)
            i_start_s = np.arange(0,n_p-1,2,dtype=np.int32)

            element_block = i_start_t[:,np.newaxis,np.newaxis] + \
                            i_start_s[np.newaxis,:,np.newaxis] + \
                            face_block[np.newaxis,np.newaxis,:]

            i_start_elm = np.arange(0,n_points, n_p**2, dtype=np.int32)
            ien = np.reshape(i_start_elm[:,np.newaxis,np.newaxis,np.newaxis]+
                             element_block[np.newaxis,:,:,:], (-1,9))
            ien[:,0] = 8
    else:
        # The base block in a 2D plane
        block = np.zeros((n_p-1,n_p-1,4), dtype=np.int32)
        for j in range(n_p-1):
            for k in range(n_p-1):
                block[j,k,:] = [n_p*j    +k  ,n_p*(j+1)+k,
                                n_p*(j+1)+k+1,n_p*j    +k+1]

        i_start = np.arange(0,n_points, n_p**2, dtype=np.int32)
        ien = np.reshape(i_start[:,np.newaxis,np.newaxis,np.newaxis]+
                         block[np.newaxis,:,:,:], (-1,4))

        # Define only _within_ an element for now
        if (n_plane > 1):
            ien_2D = ien
            ien = np.zeros((n_cells*n_cells_tor,9), dtype=np.int32)
            ien[:,0] = 8
            ien[:,1:9] = np.tile(ien_2D, (n_cells_tor,2))
            for i in range(len(phis)-1):
                ien[i*n_cells:(i+1)*n_cells,1:9] += np.concatenate(([n_points*i]*4,[n_points*(i+1)]*4))
            if (periodic):
                i = len(phis)
                ien[i*n_cells:(i+1)*n_cells,1:9] += np.concatenate(([n_points*i]*4,[0]*4))

            n_cells = n_cells * n_cells_tor
        else:
            ien = np.insert(ien, 0, 4, axis=1)

    return (xyz, ien)


"""
Calculate values of the basis functions at positions s and t
Optionally put many values of s and t at once as numpy arrays.
Dimension 0: order
Dimension 1: vertex
optional dimension 2, 3: position s, t
"""
def basis_functions(s,t):
    return np.asarray([
        [ (-1 + s)**2*(1 + 2*s)*(-1 + t)**2*(1 + 2*t),
         -(s**2*(-3 + 2*s)*(-1 + t)**2*(1 + 2*t)),
          s**2*(-3 + 2*s)*t**2*(-3 + 2*t),
         -((-1 + s)**2*(1 + 2*s)*t**2*(-3 + 2*t))],
        [ 3*(-1 + s)**2*s*(-1 + t)**2*(1 + 2*t),
         -3*(-1 + s)*s**2*(-1 + t)**2*(1 + 2*t),
         3*(-1 + s)*s**2*t**2*(-3 + 2*t),
         -3*(-1 + s)**2*s*t**2*(-3 + 2*t)],
        [ 3*(-1 + s)**2*(1 + 2*s)*(-1 + t)**2*t,
         -3*s**2*(-3 + 2*s)*(-1 + t)**2*t,
          3*s**2*(-3 + 2*s)*(-1 + t)*t**2,
         -3*(-1 + s)**2*(1 + 2*s)*(-1 + t)*t**2],
        [ 9*(-1 + s)**2*s*(-1 + t)**2*t,
         -9*(-1 + s)*s**2*(-1 + t)**2*t,
          9*(-1 + s)*s**2*(-1 + t)*t**2,
         -9*(-1 + s)**2*s*(-1 + t)*t**2]])


"""
Calculate basis functions at (n_sub+2)**2 points
n_sub is the number of subdivisions per element, + 2 for the endpoints
"""
def bf(n_sub):
    # Get the basis functions at each of the points
    lin = np.linspace(0.0, 1.0, n_sub+2)
    s  = np.tensordot(lin, [1]*(n_sub+2), axes=0)
    t  = s.transpose()
    return basis_functions(s, t)

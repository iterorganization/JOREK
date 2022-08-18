"""
jorek_read_h5.py

Read JOREK hdf5 restart files and define functions for interpolation
on the exported fields.

Before reading, set:
    - Variable numbers
Before interpolating, set:
    - Without_n0_mode (default false)

Created by Daan van Vugt <the JOREK team>

See https://www.jorek.eu/wiki/doku.php?id=jorek_read_h5

prec controls the precision of the calculation. Values: np.float32 or np.float64
Calculation time is not really affected, render time maybe?
"""
import h5py
import numpy as np
import vtk
prec=np.float32
vtk_prec=vtk.VTK_FLOAT



"""
Class encapsulating some read data from a restart file

input arguments:
    Variables: a list of all variables you want
    i_plane: which plane to interpolate at. If -1 make a 3D field.
"""
class fields(object):
    var_names = ["psi", "u", "j", "w", "rho", "T", "v_par", 'T_e', 'rho_n'] # rest is ambiguous
    def read(self, filename, variables=[], file_prev=None, interp_fraction=None):
        self.vars = variables
        with h5py.File(filename, 'r') as hf:
            self.n_var        = hf.get('n_var')[0]
            self.n_period     = hf.get('n_period')[0]
            self.n_tor        = hf.get('n_tor')[0]
            self.n_vertex_max = hf.get('n_vertex_max')[0]
            n_elements   = hf.get('n_elements')[0]
            self.tstep = hf['tstep'][0]
            self.t_now = hf['t_now'][0]

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
            self.values = interp_fraction*self.values + (1.0-interp_fraction)*read_mmap_or_h5py(file_prev, 'values', type_out=prec)


    """
    Create VTK objects from points and connectivity matrix
    input:
        n_sub: number of subdivisions per element
        n_plane: number of planes in toroidal direction
        phi: range (2 elements) or single value of phi
        without_n0_mode: do not include n0 mode if true
        quadratic: create quadratic elements instead of quads unless bezier enabled
        bezier: create anisotropic Bezier elements 
           (3rd order in poloidal plane and rational quadratic in toroidal plane). 

    returns:
        vtkUnstructuredGrid
    """
    def to_vtk(self, n_sub=2, phi=[0,360], n_plane=8, without_n0_mode=False, force_remake_grid=False,
               quadratic=True, output=None, bezier=False):
        import vtk
        from vtk.util import numpy_support as npvtk
        if (bezier):
            n_plane = 1 + (n_plane - 1) * 2
            n_sub = 4
            quadratic = False

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
            periodic = (np.mod(phi[0]-phi[1],360) < 1e-9)
            phis = np.linspace(phi[0],phi[1],num=n_plane,endpoint=not periodic)
        # Convert to radians
        phis = phis*np.pi/180


        if (quadratic):
            n_sub = n_sub*2-1

        if (output == None):
            output = vtk.vtkUnstructuredGrid()

        settings = [n_sub, phi, n_plane, without_n0_mode, quadratic, self.n_elements]
        if (not hasattr(self, 'old_settings') or self.old_settings != settings):
            force_remake_grid = True
        self.old_settings = settings

        etype = None

        if (force_remake_grid or not hasattr(self, 'points') or not hasattr(self, 'cells')):
            (xyz, ien) = create_grid(self.x, self.vertex, self.size, self.n_elements,
                                     n_sub, phis, n_plane, periodic, quadratic, bezier)

            # (Control) points from Daan Van Vugt thesis must have following index notation in vtk file
            # For bezier grid number n_planes must be 1+2n for n is natural number
            if (bezier):
                if n_plane == 1:
                    index = np.array([0, 1, 2, 3, 4, 5, 9, 10, 7, 6, 8, 11, 12, 13, 15, 14])
                    step = np.array([16 for i in range(16)])
                    ien2 = np.array([index])
                    for i in range(1, np.shape(xyz)[0] // 16):
                        ien2 = np.concatenate((ien2, np.array([index + i * step])), axis=0)
                    ien = np.insert(ien2, 0, 16, axis=1)
                    etype = vtk.VTK_BEZIER_QUADRILATERAL

                else:
                    alpha = (phi[1] - phi[0]) / (n_plane - 1)
                    s = np.shape(xyz)[0] // n_plane
                    w = np.cos(np.deg2rad(alpha))
                    w1 = np.ones((np.shape(xyz)[0]))
                    ien = None
                    index = np.array([0, 1, 2, 3, 0 + 2 * s, 1 + 2 * s, 2 + 2 * s, 3 + 2 * s, 4, 5, 9, 10, 7, 6, 8, 11,
                                      4 + 2 * s, 5 + 2 * s, 9 + 2 * s, 10 + 2 * s, 7 + 2 * s, 6 + 2 * s, 8 + 2 * s,
                                      11 + 2 * s, s, 1 + s, 3 + s, 2 + s, 8 + s, 11 + s, 9 + s, 10 + s, 4 + s, 5 + s,
                                      7 + s, 6 + s, 12, 13, 15, 14, 12 + 2 * s, 13 + 2 * s, 15 + 2 * s, 14 + 2 * s,
                                      12 + s, 13 + s, 15 + s, 14 + s])
                    for i in range((n_plane - 1) // 2):
                        index2 = index + i * np.array([s * 2 for k in range(48)])
                        w1[s + i * 2 * s: 2 * s + i * 2 * s] = np.array([w for i in range(s)])
                        xyz[s + i * 2 * s: 2 * s + i * 2 * s, 0] = 1 / w * xyz[s + i * 2 * s: 2 * s + i * 2 * s, 0]
                        xyz[s + i * 2 * s: 2 * s + i * 2 * s, 2] = 1 / w * xyz[s + i * 2 * s: 2 * s + i * 2 * s, 2]
                        step = np.array([16 for i in range(48)])
                        ien2 = np.zeros((np.shape(xyz)[0] // n_plane // 16, 48))
                        for j in range(1, np.shape(xyz)[0] // n_plane // 16):
                            ien2[j, :] = np.array([index2 + j * step])

                        if (np.any(ien)):
                            ien = np.concatenate((ien, ien2), axis=0)
                        else:
                            ien = ien2

                    ien = np.insert(ien, 0, 48, axis=1)
                    etype = vtk.VTK_BEZIER_HEXAHEDRON

                    weights = npvtk.numpy_to_vtk(w1, deep=True, array_type=vtk_prec)
                    weights.SetName("RationalWeights")
                    output.GetPointData().SetRationalWeights(weights)
                    n_c = int(np.shape(xyz)[0] / n_plane / 16 * (n_plane - 1) / 2)
                    degrees = npvtk.numpy_to_vtk(np.array([[3, 3, 2] for k in range(n_c)]), deep=True,
                                                 array_type=vtk.VTK_ID_TYPE)
                    degrees.SetName("HigherOrderDegrees")

                    output.GetCellData().SetHigherOrderDegrees(degrees)

            pcoords = npvtk.numpy_to_vtk(xyz, deep=True, array_type=vtk_prec)
            self.points = vtk.vtkPoints()
            self.points.SetData(pcoords)

            self.cells = vtk.vtkCellArray()
            self.cells.SetCells(ien.shape[0], npvtk.numpy_to_vtk(ien, deep=True, array_type=vtk.VTK_ID_TYPE))

            self.xyz = xyz # expose for IMAS
            self.ien = ien

        output.SetPoints(self.points)

        HZ = toroidal_basis(self.n_tor, self.n_period, phis, without_n0_mode)

        val = interp_scalars_3D(self.values[self.vars,:,:,:],
                                    self.vertex, self.size, n_sub, HZ).reshape((len(self.vars),-1))

        if (bezier):
            a = np.shape(val)
            val = val.reshape((a[0], a[1] // 16, 16))
            val[:, :, [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]] = val[:, :, [0, 12, 15, 3, 4, 8, 11, 7, 1,
                                                                                           13, 14, 2, 5, 9, 10, 6]]
            val = val.reshape((a[0], a[1]))

        # Could split here if we run into memory problems and delete each part after use

        for i in range(len(self.vars)):
            tmp = npvtk.numpy_to_vtk(val[i,:], deep=True, array_type=vtk_prec)
            tmp.SetName(self.var_names[self.vars[i]])
            output.GetPointData().AddArray(tmp)

        self.val = val # Expose for IDS

        if not(etype):
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
def grid_2D(x, vertex, size, n_sub, bezier=False):
    # Calculate RZ for all of the elements (dimension 0) for each of the s
    # positions (dimension 1) for each of the t positions (dimension 2)
    # Multiply x[var, order, node[vertex, element]] on the last two dimensions
    # with size[order, vertex, element]*bf[order, vertex, s, t]
    # See http://stackoverflow.com/questions/26089893/understanding-numpys-einsum
    #return np.einsum('lijk,ijk,ijmn->kmnl', x[:,:,vertex], size, bf(n_sub))
    # Code below is ~5x faster or so! try again when einsum supports optimize=True
    # First create a temporary array holding: x[order, vertex, element, var]
    tmp = np.zeros((x.shape[0], x.shape[1],vertex.shape[0],vertex.shape[1]))
    # Fill it with the right x
    if x.ndim == 4: # New format (after Feb 2021)
        for i in range(vertex.shape[0]): # small loop over vertices (hardcode 4 here?)
            tmp[:,:,i,:] = x[:,:,0,vertex[i,:]-1]
    else:
        for i in range(vertex.shape[0]): 
            tmp[:,:,i,:] = x[:,:,vertex[i,:]-1]

    # multiply by size[order, vertex, element]
    tmp[0,:,:,:] *= size
    tmp[1,:,:,:] *= size
    # Create output array
    if (bezier):
        xy = np.zeros((np.shape(tmp)[3] * 16, 2))
        for i in range(np.shape(tmp)[3]):
            x1 = tmp[0, :, :, i]
            y1 = tmp[1, :, :, i]
            x1[3, :] = x1[3, :] + x1[1, :] + x1[2, :]
            y1[3, :] = y1[3, :] + y1[1, :] + y1[2, :]
            x1[1:, :] += np.tile(x1[0, :], (3, 1))
            y1[1:, :] += np.tile(y1[0, :], (3, 1))
            xy[16 * i:16 * (i + 1), 0] = np.ravel(x1)
            xy[16 * i:16 * (i + 1), 1] = np.ravel(y1)
            out = xy
    else:
        out = np.zeros((vertex.shape[1],n_sub,n_sub,2), dtype=prec)
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
Create a grid of nsub**2 points per element, at phis positions
return points and connectivity matrix
"""
def create_grid(x, vertex, size, n_elements, n_sub, phis, n_plane, periodic, quadratic, bezier):
    RZ     = grid_2D(x, vertex, size, n_sub, bezier)

    # Create connectivity data
    # Calculate 2D connectivity first
    # For each element, calculate the number of the lowest point
    # Create (n_sub-1)**2 quadrangles
    n_points = n_elements*(n_sub**2) # number of points in one plane
    n_cells  = n_elements*((n_sub-1)**2) # Number of cells in one plane
    if (n_plane > 1): # Create a volume
        if (periodic):
            n_cells_tor = n_plane
        else:
            n_cells_tor = n_plane - 1
    else:
        n_cells_tor = 1

    if (bezier):
        n_xy = np.shape(RZ)[0]
        xyz = np.zeros((n_xy * n_plane, 3))
        for i in range(n_plane):
            xyz[i * n_xy:(i + 1) * n_xy, 0] = np.ravel(RZ[:, 0] * np.cos(phis[i]))
            xyz[i * n_xy:(i + 1) * n_xy, 1] = np.ravel(RZ[:, 1])
            xyz[i * n_xy:(i + 1) * n_xy, 2] = np.ravel(RZ[:, 0] * np.sin(phis[i]))
    else:
        xyz = np.zeros((n_points*n_plane,3))
        for i in range(n_plane):
            xyz[i*n_points:(i+1)*n_points,0] = np.ravel(RZ[:,:,:,0]*np.cos(phis[i]))
            xyz[i*n_points:(i+1)*n_points,1] = np.ravel(RZ[:,:,:,1])
            xyz[i*n_points:(i+1)*n_points,2] = np.ravel(RZ[:,:,:,0]*np.sin(phis[i]))

    ien = None

    if not(bezier):
        # See http://www.vtk.org/doc/nightly/html/classvtkQuadraticHexahedron.html#details
        if (quadratic):
            if (n_plane > 1):
                # The base block for a quadratic element
                face_block = np.zeros(21, dtype=np.int32)
                face_block[1:21] = [0,2*n_sub,2*n_sub+2,2, # corners front face
                                    2*n_points,2*n_points+2*n_sub,2*n_points+2*n_sub+2,2*n_points+2, # corners back face
                                    n_sub,2*n_sub+1,n_sub+2,1,# midedges front face
                                    2*n_points+n_sub,2*n_points+2*n_sub+1,2*n_points+n_sub+2,2*n_points+1,# midedges back face
                                    n_points,2*n_sub+n_points,2*n_sub+2+n_points,2+n_points]# midedges middle
                i_start_t = np.arange(0,n_sub*(n_sub-1),2*n_sub,dtype=np.int32)
                i_start_s = np.arange(0,n_sub-1,2,dtype=np.int32)

                element_block = i_start_t[:,np.newaxis,np.newaxis] + \
                                i_start_s[np.newaxis,:,np.newaxis] + \
                                face_block[np.newaxis,np.newaxis,:]

                i_start_elm   = np.arange(0,n_points, n_sub**2, dtype=np.int32)
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
                face_block[1:9] = [0,2*n_sub,2*n_sub+2,2, # corners
                                   n_sub,2*n_sub+1,n_sub+2,1] # midedges
                i_start_t = np.arange(0,n_sub*(n_sub-1),2*n_sub,dtype=np.int32)
                i_start_s = np.arange(0,n_sub-1,2,dtype=np.int32)

                element_block = i_start_t[:,np.newaxis,np.newaxis] + \
                                i_start_s[np.newaxis,:,np.newaxis] + \
                                face_block[np.newaxis,np.newaxis,:]

                i_start_elm = np.arange(0,n_points, n_sub**2, dtype=np.int32)
                ien = np.reshape(i_start_elm[:,np.newaxis,np.newaxis,np.newaxis]+
                                 element_block[np.newaxis,:,:,:], (-1,9))
                ien[:,0] = 8
        else:
            # The base block in a 2D plane
            block = np.zeros((n_sub-1,n_sub-1,4), dtype=np.int32)
            for j in range(n_sub-1):
                for k in range(n_sub-1):
                    block[j,k,:] = [n_sub*j    +k  ,n_sub*(j+1)+k,
                                    n_sub*(j+1)+k+1,n_sub*j    +k+1]

            i_start = np.arange(0,n_points, n_sub**2, dtype=np.int32)
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
         -(-1 + s)**2*(1 + 2*s)*t**2*(-3 + 2*t)],
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
Calculate values of the basis functions derived to s at positions s and t
Optionally put many values of s and t at once as numpy arrays.
Dimension 0: order
Dimension 1: vertex
optional dimension 2, 3: position s, t
"""
def basis_functions_s(s,t):
    return np.asarray([
        [ 6*(-1 + s)*s*(-1 + t)**2*(1 + 2*t),
         -6*(-1 + s)*s*(-1 + t)**2*(1 + 2*t),
          6*(-1 + s)*s*t**2*(-3 + 2*t),
         -6*(-1 + s)**2*t**2*(-3 + 2*t)],
        [ 3*(-1 + s)*(-1+3*s)*(-1 + t)**2*(1 + 2*t),
         -3*s*(-2 + 3*s)*(-1 + t)**2*(1 + 2*t),
          3*s*(-2 + 3*s)*t**2*(-3 + 2*t),
         -3*(-1 + 3*s)*(-1 + s)*t**2*(-3 + 2*t)],
        [ 18*(-1 + s)*s*(-1 + t)**2*t,
         -18*(-1 + s)*s*(-1 + t)**2*t,
          18*(-1 + s)*s*(-1 + t)*t**2,
         -18*(-1 + s)*s*(-1 + t)*t**2],
        [ 9*(-1 + s)*(-1+3*s)*(-1 + t)**2*t,
         -9*s*(-2 + 3*s)*(-1 + t)**2*t,
          9*s*(-2 + 3*s)*(-1 + t)*t**2,
         -9*(-1 + 3*s)*(-1 + s)*(-1 + t)*t**2]])

"""
Calculate values of the basis functions derived to t at positions s and t
Optionally put many values of s and t at once as numpy arrays.
Dimension 0: order
Dimension 1: vertex
optional dimension 2, 3: position s, t
"""
def basis_functions_t(s,t):
    return np.asarray([
        [ 6*(-1 + s)**2*(1 + 2*s)*(-1 + t)*t,
         -6*s**2*(-3 + 2*s)*(-1 + t)*t,
          6*s**2*(-3 + 2*s)*(-1 + t)*t,
         -6*(-1 + s)**2*(1 + 2*s)*(-1 + t)*t],
        [ 18*(-1 + s)**2*s*(-1 + t)*t,
         -18*(-1 + s)*s**2*(-1 + t)*t,
          18*(-1 + s)*s**2*(-1 + t)*t,
         -18*(-1 + s)**2*s*(-1 + t)*t],
        [ 3*(-1 + s)**2*(1 + 2*s)*(-1 + t)*(-1 + 3*t),
         -3*s**2*(-3 + 2*s)*(1 - 3*t)*(-1 + t),
          3*s**2*(-3 + 2*s)*t*(-2 + 3*t),
         -3*(-1 + s)**2*(1 + 2*s)*t*(-2 + 3*t)],
        [ 9*(-1 + s)**2*s*(-1 + t)*(-1 + 3*t),
         -9*(-1 + s)*s**2*(-1 + t)*(-1 + 3*t),
          9*(-1 + s)*s**2*t*(-2 + 3*t),
         -9*(-1 + s)**2*s*t*(-2 + 3*t)]])


"""
Calculate basis functions at n_sub**2 points
"""
def bf(n_sub):
    # Get the basis functions at each of the points
    lin = np.linspace(0.0, 1.0, n_sub)
    s  = np.tensordot(lin, [1]*n_sub, axes=0)
    t  = s.transpose()
    return basis_functions(s, t)
def bf_s(n_sub):
    # Get the basis functions at each of the points
    lin = np.linspace(0.0, 1.0, n_sub)
    s  = np.tensordot(lin, [1]*n_sub, axes=0)
    t  = s.transpose()
    return basis_functions_s(s, t)
def bf_t(n_sub):
    # Get the basis functions at each of the points
    lin = np.linspace(0.0, 1.0, n_sub)
    s  = np.tensordot(lin, [1]*n_sub, axes=0)
    t  = s.transpose()
    return basis_functions_t(s, t)

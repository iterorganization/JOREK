#!local/bin/python
#   Name : jorekHDF5toIDS.py
#
#   Description :
#       A script which reads JOREK HDF5 output file(s) and writes its
#       contents (grid geometry as Bezier elements, quad connectivity,
#       data fields with Fourer harmonics) to IMAS MHD IDS.
#  Main features:
#      1. There are two grid_ggd spaces:
#         - Space 1 is two-dimensional (R,Z) space of unstructured grid
#           with geometry_2d associated to nodes and cells.
#         - Space 2 is one-dimensional space with $\phi$ and
#           geometry_type.index = 1 (Fourier) with geometry
#           objects (1, 2, 3, ..., number of harmonics)
#     2. Lossless (Te, n, w, ...), stored under ggd, with the above two
#        spaces form a "structured" (implicitly defined) grid of node
#        values, where explicitly RZ values of first harmonics are saved
#        first, then RZ values of the second harmonics follow up to the
#        last RZ harmonics. Similarly, coefficients on the nodes are saved.
#        This definition follows column major (FORTRAN) notation,
#        meaning that with varying first index (R) the values are close
#        together in the memory and that the last index
#        (phi in Fourier space) defines RZ block of values in memory.
#
#       JOREK HDF5 file contents:
#           https://www.jorek.eu/wiki/doku.php?id=hdf5-tools&s[]=hdf5
#
#   Requirements:
#       - h5py
#       - IMAS
#       - mkdir -p $HOME/public/imasdb/jorek/3/0
#
#   Author :
#       Leon Kos and Miha Radež
#   E-mail :
#       leon.kos@lecad.fs.uni-lj.si
#

import numpy as np
import sys
import getpass
import argparse

import h5py
import vtk

import imas
from imas import imasdef

prec=np.float32
vtk_prec=vtk.VTK_FLOAT

parser = argparse.ArgumentParser(description="Convert JOREK HDF5 file(s) to IMAS (IDSs)",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-s", "--shot", type=int, default=1,
                    help="Shot number")
parser.add_argument("-r", "--run", type=int, default=10,
                    help="Run number")
parser.add_argument("-u", "--user", type=str, default=getpass.getuser(),
                    help="Location of ~$USER/public/imasdb")
parser.add_argument("-d", "--database", type=str, default="jorek",
                    help="Database name under public/imasdb/")
parser.add_argument("-o", "--occurrence", type=int, default=0,
                    help="Occurrence number")
parser.add_argument("-f", "--backend", type=int, default=imasdef.MDSPLUS_BACKEND,
                    help="Database format: 12=MDSPLUS, 13=HDF5")
parser.add_argument("hdf5files", metavar='jorek?????.h5', nargs='*',
                    help="JOREK HDF5 file(s)", default=["/tmp/jorek_restart.h5"])

args = parser.parse_args()

np.set_printoptions(
    formatter={'float': '{: 0.3f}'.format})

# Open Database entry
data_entry = imas.DBEntry(args.backend, args.database, args.shot, args.run,
                          user_name=args.user)
status, idx = data_entry.create()
if status != 0:
    print("Creation of data entry FAILED!")
    sys.exit(status)

comment = "Written results of multiple JOREK output HDF5 files/timeslices."
# Create Grid GGD
mhd_ids = imas.mhd()
mhd_ids.ids_properties.homogeneous_time = 1
mhd_ids.time.resize(1)
mhd_ids.time = np.array([0.0])
mhd_ids.grid_ggd.resize(1)
mhd_grid_ggd = mhd_ids.grid_ggd
# w_ids.data_entry.mhd.time.resize(len(filePathList))
allTimeValues = np.array([0]*len(args.hdf5files))

# Resize GGD to the number of timeslices
mhd_ids.ggd.resize(len(args.hdf5files))
mhd_ids.time.resize(len(args.hdf5files))

# Loop through the list of HDF5 files
for i_slice in range(len(args.hdf5files)):

    if len(args.hdf5files) > 1 :
        print(f"[{i_slice}]", end="", flush=True)

    # Read file
    with h5py.File(args.hdf5files[i_slice], 'r') as hf:
        n_var        = hf.get('n_var')[0]
        n_period     = hf.get('n_period')[0]
        n_tor        = hf.get('n_tor')[0]
        n_vertex_max = hf.get('n_vertex_max')[0]
        n_elements   = hf.get('n_elements')[0]
        vertex       = np.array(hf.get('vertex'))
        x            = np.array(hf.get('x'))
        size         = np.array(hf.get('size'))
        values       = np.array(hf.get('values'))
        tstep        = hf['tstep'][0]
        t_now        = hf['t_now'][0]
        t_norm       = hf.get('t_norm', default=[1])[0]

    # Grid geometry is taken only from the first time slice. All other
    # time slices share the same geometry. No need to write the same grid
    # multiple times as all of them have the same.
    if i_slice == 0:

        #ien for quad
        ien = np.swapaxes(vertex, 1, 0) - 1
        vtk_quad_conn_array = np.insert(ien, 0, 4, axis=1)
        quad_conn_array = ien

        #xyz for quad
        dimensions = np.ndim(x)-1
        if dimensions == 3: # New format 
            r0 = x[0, 0, 0]
            z0 = x[0, 1, 0]
        else: # format before February 2021
            r0 = x[0, 0]
            z0 = x[1, 0]
        rz = np.zeros((np.shape(r0)[0], 2))
        rz[:, 0] = r0
        rz[:, 1] = z0
        

        #val for quad
        val0 = values[:, 0, :, :]
        tor = [1, 1, 0, 1, 0]
        #val = np.einsum('ijk,j->ik', val0, tor)
        
        # # Changing arrays notation from Python (starting with 0) to
        # # Fortran notation (starting with 1)
        obj_2D_list_f90 = np.array(quad_conn_array)
        obj_2D_list_f90 = obj_2D_list_f90 + 1

        mhd_grid_ggd = mhd_ids.grid_ggd

        # Write grid geometry
        mhd_grid_ggd[0].space.resize(2)
        mhd_grid_ggd[0].space[0].objects_per_dimension.resize(4)

        mhd_grid_ggd[0].space[0].coordinates_type.resize(2)
        # Set coordinates type to [R, Z]
        mhd_grid_ggd[0].space[0].coordinates_type = np.array([4, 3])
        # mhd_grid_ggd[0].space[0].coordinates_type = np.array([1, 2, 3])

        mhd_grid_ggd[0].identifier.description = "Mesh JOREK output HDF5 file grid with quantities"
        mhd_grid_ggd[0].identifier.name = "JOREK output HDF5 file grid with quantities"
        mhd_grid_ggd[0].identifier.index = 1

        # Writing nodes
        points = rz[:, :]
        num_points = len(points)
        ids_space = mhd_grid_ggd[0].space[0]
        ids_space.objects_per_dimension[0].object.resize(num_points)

        # For nodes
        ids_dim_0D = ids_space.objects_per_dimension[0]
        # For edges. No edges, so we create a dummy object...
        ids_dim_1D = ids_space.objects_per_dimension[1]
        # For cells
        ids_dim_2D = ids_space.objects_per_dimension[2]

        for i in range(num_points):
            ids_dim_0D.object[i].geometry.resize(2)
            ids_dim_0D.object[i].geometry = points[i]

        # Now to write dummy 1D object
        ids_dim_1D.object.resize(1)
        ids_dim_1D.object[0].nodes.resize(1)
        ids_dim_1D.object[0].nodes[0] = 0

        # Writing cells
        num_obj_2D = len(obj_2D_list_f90)
        ids_dim_2D.object.resize(num_obj_2D)
        cell_size = len(obj_2D_list_f90[0])
        for i in range(num_obj_2D):
            ids_dim_2D.object[i].nodes.resize(cell_size)
            for j in range(cell_size):
                ids_dim_2D.object[i].nodes[j] = obj_2D_list_f90[i][j]


        # Writing grid_subsets
        num_grid_subset = 2 # points, 0D, 2D
        mhd_grid_ggd[0].grid_subset.resize(num_grid_subset)
        gs_index = 0

        # First for points
        gs = mhd_grid_ggd[0].grid_subset[gs_index]
        gs.identifier.name = "Nodes"
        gs.identifier.index = gs_index + 1 # Fortran notation
        gs.identifier.description = """All points/nodes/vertices/0D
                                  objects in the domain."""
        gs.dimension = 1
        gs.element.resize(num_points)
        for j in range(num_points):
            gs.element[j].object.resize(1)
            # Write in Fortran notation (!)
            gs.element[j].object[0].space = 1
            gs.element[j].object[0].index = j + 1 # Fortran notation
            gs.element[j].object[0].dimension = 1

        gs_index += 1

        # Now for cells
        gs = mhd_grid_ggd[0].grid_subset[gs_index]
        gs.identifier.name = "Cells2D"
        gs.identifier.index = gs_index + 1 # Fortran notation
        gs.identifier.description = "All 2D cells/2D objects in the domain."
        gs.dimension = 3
        gs.element.resize(num_obj_2D)
        for j in range(num_obj_2D):
            gs.element[j].object.resize(1)
            # Write in Fortran notation (!)
            gs.element[j].object[0].space = 1
            gs.element[j].object[0].index = j + 1 # Fortran notation
            gs.element[j].object[0].dimension = 3

        mhd_grid_ggd = mhd_ids.grid_ggd
        # w_ids.grid_ggd.array[0].space.resize(2)  # Assure second space is allocated
        toroidal_space = mhd_grid_ggd[0].space[1]  # TODO multiple time slices
        toroidal_space.coordinates_type.resize(1)  # One dimensional (toroidal) space
        toroidal_space.coordinates_type[0] = 5  # Toroidal angle (phi) in radians
        toroidal_space.geometry_type.index = n_period  # Fourier with periodicity
        toroidal_space.identifier.description = "Toroidal space"

        toroidal_space.objects_per_dimension.resize(1)  # We have only one dimension of
        toroidal_space.objects_per_dimension[0].object.resize(n_tor)  # toroidal harmonics

        for tor in range(n_tor):
            toroidal_space.objects_per_dimension[0].object[tor].geometry = np.array([tor+1])

        gr2d = mhd_grid_ggd[0].space[0]
        gr2d.geometry_type.index = 0  # Standard geometry (non Fourier)
        x_shape = np.shape(x)[dimensions]
        # coordinate and derivates (s, t, mixed)
        if dimensions == 3:
            for j in range(x_shape):
                gr2d.objects_per_dimension[0].object[j].geometry_2d = x[:, :, 0, j]
        else:
            for j in range(x_shape):
                gr2d.objects_per_dimension[0].object[j].geometry_2d = x[:, :, j]


        # size 1, d_{uk}, d_{vk}, d{uv}d{vk} as in Daan Van Vugt thesis
        size_shape = np.shape(size)[2]
        for i in range(size_shape):
            gr2d.objects_per_dimension[2].object[i].geometry_2d = size[:, :, i]



        
    quantity_names_list = ["psi", "u", "j", "w", "rho", "T", "v_par"]
    
    # Set empty IDS path for quantity tree node
    IDSQuantityPath = None

    # Set time
    mhd_ids.ggd[i_slice].time = t_now*t_norm
    mhd_ids.time[i_slice] = t_now*t_norm

    for i in range(len(quantity_names_list)):
        label = quantity_names_list[i]
        if label == 'psi': # Flux / poloidal magnetic flux
            mhd_ids.ggd[i_slice].psi.resize(1)
            IDSQuantityPath = mhd_ids.ggd[i_slice].psi[0]
            arr = values[0, :, :, :]
            # Swap axes
        elif label == 'u': # Potential / electric potential
            mhd_ids.ggd[i_slice].phi_potential.resize(1)
            IDSQuantityPath = mhd_ids.ggd[i_slice].phi_potential[0]
            arr = values[1, :, :, :]
        elif label == 'j': # Current * R/ toroidal current density * R
            mhd_ids.ggd[i_slice].j_tor_r.resize(1)
            IDSQuantityPath = mhd_ids.ggd[i_slice].j_tor_r[0]
            arr = values[2, :, :, :]
        elif label == 'w': # Vorticity / R
            mhd_ids.ggd[i_slice].vorticity_over_r.resize(1)
            IDSQuantityPath = mhd_ids.ggd[i_slice].vorticity_over_r[0]
            arr = values[3, :, :, :]
        elif label == 'rho': # Mass density
            mhd_ids.ggd[i_slice].mass_density.resize(1)
            IDSQuantityPath = mhd_ids.ggd[i_slice].mass_density[0]
            arr = values[4, :, :, :]
        elif label == 'T': # Current / toroidal current density
            mhd_ids.ggd[i_slice].electrons.temperature.resize(1)
            IDSQuantityPath = mhd_ids.ggd[i_slice].electrons.temperature[0]
            arr = values[5, :, :, :]
        elif label == 'v_par': # Current / toroidal current density
            mhd_ids.ggd[i_slice].velocity_parallel_over_b_field.resize(1)
            IDSQuantityPath = mhd_ids.ggd[i_slice].velocity_parallel_over_b_field[0]
            arr = values[6, :, :, :]
        else:
            print(f"Unknown label: {label}")
            continue

        # Swap and correct axes
        arr = np.swapaxes(arr, 0, 1)
        shape = np.shape(arr)
        reshape_arr = np.reshape(arr, (shape[0], shape[1] * shape[2]))


        IDSQuantityPath.grid_index = 1
        IDSQuantityPath.grid_subset_index = 1
        #IDSQuantityPath.values = val[i, :]
        IDSQuantityPath.coefficients = reshape_arr
print("[put...]", end="", flush=True)
data_entry.put(mhd_ids)
data_entry.close()
print("[OK]")
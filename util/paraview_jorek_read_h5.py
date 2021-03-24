Name = 'JorekReadH5'
Label = 'JorekH5Reader'
Help = 'Read a jorek restart file in H5 format'
Extension = 'h5'
FileDescription = 'JOREK hdf5 files'

import h5py
import inspect
import numpy as np
from os.path import dirname
from collections import OrderedDict
PythonPaths = [dirname(dirname(inspect.getfile(h5py)))]

NumberOfInputs = 0
OutputDataType = 'vtkUnstructuredGrid'
# Top-level tuples are propertygroups. One level below is a list of dictionaries for the elements
# Names are converted to labels by '_' > ' ' and uppercasing the first letter

# Create a new class to distinguish the list to select
# Keys = string names
# Values = 1 or 0 (enabled by default or not)
class ArraySelectionDomain(OrderedDict):
    pass
class PropertyGroup(OrderedDict):
    pass

Properties = OrderedDict([
    ('variables', ArraySelectionDomain([
        ('Psi', 1), ('u', 1), ('j', 1), ('w', 1), ('rho', 1), ('T', 1),
        ('v_par', 0), ('T_e', 0), #model400
        ('rho_n', 0), # model500
    ])),
    ('interpolation_options', PropertyGroup([
        ('number_of_subdivisions', dict(value=3, min=2, max=6)),
        ('quadratic', False),
        ('exclude_n0_mode', False),
        ('smooth_time_interpolation', True),
    ])),
    ('toroidal_direction', PropertyGroup([
        ('number_of_planes', dict(value=1, min=1, max=360)),
        ('phi_range', dict(value=[0.0, 90.0], widget='double_range', min=0.0, max=360)),
    ])),
])


# from paraview import vtk # is done automatically
def RequestData(self):
    import re, os
    #include 'jorek_read_h5.py' # Don't delete: include content of this file here

    def GetUpdateTimestep(algorithm):
        """Returns the requested time value, or None if not present"""
        executive = algorithm.GetExecutive()
        outInfo = executive.GetOutputInformation(0)
        return outInfo.Get(executive.UPDATE_TIME_STEP()) \
                if outInfo.Has(executive.UPDATE_TIME_STEP()) else None
    # Get the current timestep
    req_time = GetUpdateTimestep(self)

    # Read the timestep info from the files
    t_norm = h5py.File(FileNames[0]).get("t_norm", default=[1])[0]
    # Read xtime from last file and correlate against file numbers
    f = h5py.File(FileNames[-1])
    if "xtime" in f:
        xtime_all = np.insert(f.get("xtime")*t_norm, 0, 0.0)
    else:
        xtime_all = np.asarray([0.0])
    if len(FileNames) > 1:
        xtime = np.asarray([xtime_all[int(re.findall(r'\d+', os.path.basename(fname))[0])] for fname in FileNames])
    else:
        xtime = np.asarray([xtime_all[-1]])

    # 4 possibilities here:
    # After last step: return last file
    # Before first step: return first file
    # Very close match: return that file
    # Between 2 files: interpolate
    interp = False
    try:
        index = np.isclose(xtime, req_time).tolist().index(True)
        # We have a match, read and return data for this filename
    except ValueError:
        # Check for other 3 cases
        if (np.count_nonzero(req_time < xtime) == 0):
            # After last file
            index = len(FileNames)-1
        elif (np.count_nonzero(req_time > xtime) == 0):
            # Before first file
            index = 0
        else:
            interp = True
            index = (xtime > req_time).tolist().index(True)
            interp_fraction = (req_time - xtime[index-1])/(xtime[index] - xtime[index-1])
            if (interp_fraction > 1 or interp_fraction < 0):
                print("ERROR in interp_fraction", interp_fraction, xtime, index)
            # how much of second to take == 1-how much of first to take

    # Make a list of variables to read (0-based) and a list of variables to interpret
    # See https://www.jorek.eu/wiki/doku.php?id=models
    assert self.model < 700, 'Model700 and above restart files are not supported, model=%s'%(self.model,)
    to_read = []
    # Model-agnostic variables first
    if (Psi): to_read.append(0)
    if (u): to_read.append(1)
    if (j): to_read.append(2)
    if (w): to_read.append(3)
    if (rho): to_read.append(4)
    if (T): to_read.append(5)
    if self.model > 199:
        if (v_par): to_read.append(6)
    if self.model >= 400 and self.model < 500:
        if (T_e): to_read.append(7)
    if self.model >= 500 and self.model < 600:
        if (rho_n): to_read.append(7)
    if self.model == 712:
        if (rho_n): to_read.append(7)




    # Read the h5 file
    if (not hasattr(self, 'f')):
        self.f = fields()
    if (smooth_time_interpolation and interp):
        self.f.read(FileNames[index], variables=to_read, file_prev=FileNames[index-1],
               interp_fraction=interp_fraction)
    else:
        self.f.read(FileNames[index], variables=to_read)

    output = self.f.to_vtk(n_sub=number_of_subdivisions, phi=phi_range,
                      n_plane=number_of_planes, without_n0_mode=exclude_n0_mode,
                      output=self.GetUnstructuredGridOutput(), quadratic=quadratic)
    return output

"""
See paraview guide 13.2.2
"""
def RequestInformation(self):
    import re, os
    import numpy as np
    import h5py
    def setOutputTimesteps(algorithm):
        executive = algorithm.GetExecutive()
        outInfo = executive.GetOutputInformation(0)

        t_norm = h5py.File(FileNames[0]).get("t_norm", default=[1])[0]

        # Read xtime from last file and correlate against file numbers
        # Also read the model number and set this
        with h5py.File(FileNames[-1], 'r') as f:
            self.model = f['jorek_model'][0]
            if "xtime" in f:
                xtime_all = np.insert(f.get("xtime")*t_norm, 0, 0.0)
            else:
                xtime_all = np.asarray([0.0])
        if len(FileNames) > 1:
            xtime = [xtime_all[int(re.findall(r'\d+', os.path.basename(fname))[0])] for fname in FileNames]
        else:
            xtime = [xtime_all[-1]]

        outInfo.Remove(executive.TIME_STEPS())
        for i in range(len(FileNames)):
            outInfo.Append(executive.TIME_STEPS(), xtime[i])

        # Remove time range info
        outInfo.Remove(executive.TIME_RANGE())
        outInfo.Append(executive.TIME_RANGE(), xtime[0])
        outInfo.Append(executive.TIME_RANGE(), xtime[-1])

    setOutputTimesteps(self)

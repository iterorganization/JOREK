paraview_plugin_name = "Read JOREK files in H5 format"
paraview_plugin_version = "1.1"

## TODO 
#    @smhint.xml("""
#        <PropertyGroup label="Interpolation options"
#           <Property name="User Quadratic elements" />
#           <Property name="Exclude n0 mode" />
#        </PropertyGroup>
#        """)

from paraview.util.vtkAlgorithm import smproxy, smproperty, smdomain, smhint
from vtkmodules.util.vtkAlgorithm import VTKPythonAlgorithmBase

import jorek.jorek_read_h5 as jorek

import h5py
import inspect
import numpy as np
from os.path import dirname

def createModifiedCallback(anobject):
    import weakref
    weakref_obj = weakref.ref(anobject)
    anobject = None
    def _markmodified(*args, **kwars):
        o = weakref_obj()
        if o is not None:
            o.Modified()
    return _markmodified

@smproxy.reader(name="JOREK HDF5 restart file reader", label="JOREK HDF5 Reader",
                extensions="h5", file_description="HDF5 files")
class PythonJorekHDF5Reader(VTKPythonAlgorithmBase):
    """A reader that reads a HDF5 file. If the HDF5 has a "time" column, then
    the data is treated as a temporal dataset"""
    def __init__(self):
        VTKPythonAlgorithmBase.__init__(self, 
                                        nInputPorts=0, 
                                        nOutputPorts=1, 
                                        outputType='vtkUnstructuredGrid')
        self._filenames = []
        self._number_of_values = None
        self._fields = None
        self._timesteps = None
        
        self.__n_plane = 2
        self.__phi_start = 0.0
        self.__phi_end = 90.0
        self.__bezier = True
        self.__quadratic = True
        self.__exclude_n0_mode = False
        self.__smooth_time_interpolation = True
        self.__number_of_subdivisions = 3



        from vtkmodules.vtkCommonCore import vtkDataArraySelection
        self._arrayselection = vtkDataArraySelection()
        self._arrayselection.AddObserver("ModifiedEvent", createModifiedCallback(self))

    def _get_raw_data(self, requested_time=None):
        if self._ndata is not None:
            if requested_time is not None:
                return self._ndata[self._ndata["time"]==requested_time]
            return self._ndata

        if self._filename is None:
            # Note, exceptions are totally fine!
            raise RuntimeError("No filename specified")

        import numpy
        self._ndata = numpy.genfromtxt(self._filename, dtype=None, names=True, delimiter=',', autostrip=True)
        self._timesteps = None
        if "time" in self._ndata.dtype.names:
            self._timesteps = numpy.sort(numpy.unique(self._ndata["time"]))

        for aname in self._ndata.dtype.names:
            # note, this doesn't change MTime on the array selection, which is
            # good!
            self._arrayselection.AddArray(aname)
        return self._get_raw_data(requested_time)

    def _get_timesteps(self):
        import re, os
        import numpy as np
        import h5py
        t_norm = float(h5py.File(self._filenames[0]).get("t_norm", default=[1])[0])

        # Read xtime from last file and correlate against file numbers
        # Also read the model number and set this
        with h5py.File(self._filenames[-1], 'r') as f:
            self.model = f['jorek_model'][0]
            # print(f'jorek_model = {self.model}')
            if "xtime" in f:
                xtime = f.get("xtime").astype(float)[:]
                xtime_all = np.insert(f.get("xtime").astype(float)[:]*t_norm, 0, 0.0)
            else:
                xtime_all = np.asarray([0.0])
            self._number_of_values = f.get('values').len()

        if len(self._filenames) > 1:
            xtime = [xtime_all[int(re.findall(r'\d+', os.path.basename(fname))[0])] for fname in self._filenames]
        else:
            xtime = [xtime_all[-1]]

        return xtime

    def _get_array_selection(self):
        return self._arrayselection

    @smproperty.xml("""
        <StringVectorProperty name="FileNames"
            animateable="0"
            clean_command="RemoveAllFileNames"
            command="AddFileName"          
            number_of_elements="0"
            panel_visibility="never"
            repeat_command="1">
            <FileListDomain name="files" />
            <Documentation>The list of files to be read by the reader.</Documentation>
        </StringVectorProperty>
      """)
    def AddFileName(self, name):
        self._filenames.append(name)

    def RemoveAllFileNames(self):
        self._filenames = []

    @smproperty.xml("""
        <IntVectorProperty name="Use Bezier elements with nonlinear subdivisions"
            animateable="1"
            command="SetBezier"          
            number_of_elements="1"
            default_values="1">
            <BooleanDomain name="bool" />
            <Documentation>Bezier elements do not require subdivision since 
            the subdivision in controlled by "Nonlinear Subdivision Level" 
            by ParaView under Miscellaneous Properties after Apply. 
            If unchecked Quadratic or Hexa elements are used with
            number of subdivisions specified.</Documentation>
        </IntVectorProperty>
      """)
    def SetBezier(self, val):
        if bool(val)  != self.__bezier:
            self.__bezier = bool(val)
            self.Modified()

    @smproperty.xml("""
        <IntVectorProperty name="Use Quadratic elements"
            animateable="1"
            command="SetQuadratic"          
            number_of_elements="1"
            default_values="1">
            <BooleanDomain name="bool" />
        </IntVectorProperty>
      """)
    def SetQuadratic(self, val):
        if bool(val) != self.__quadratic:
            self.__quadratic = bool(val)
            self.Modified()

    @smproperty.xml("""
        <IntVectorProperty name="Exclude n0 mode"
            animateable="1"
            command="SetN0Mode"          
            number_of_elements="1"
            default_values="0">
            <BooleanDomain name="bool" />
        </IntVectorProperty>
      """)
    def SetN0Mode(self, val):
        if bool(val) != self.__exclude_n0_mode:
            self.__exclude_n0_mode = bool(val)
            self.Modified()


    @smproperty.xml("""
        <IntVectorProperty name="Smooth time interpolation"
            animateable="1"
            command="SetSmoothTimeInterpolation"          
            number_of_elements="1"
            default_values="1">
            <BooleanDomain name="bool" />
            <Documentation>When set the slice is interpolated between two nearest slices</Documentation>
        </IntVectorProperty>
      """)
    def SetSmoothTimeInterpolation(self, val):
        if bool(val) != self.__smooth_time_interpolation:
            self.__smooth_time_interpolation = bool(val)
            self.Modified()


    @smproperty.intvector(name="Number of subdivisions", default_values=3)
    @smdomain.intrange(min=0, max=6)
    def SetNumberOfSubdivisions(self, val):
        if val != self.__number_of_subdivisions:
            self.__number_of_subdivisions = val
            self.Modified()


    @smproperty.doublevector(name="TimestepValues", information_only="1", si_class="vtkSITimeStepsProperty")
    def GetTimestepValues(self):
        return self._get_timesteps()

    # Array selection API is typical with readers in VTK
    # This is intended to allow ability for users to choose which arrays to
    # load. To expose that in ParaView, simply use the
    # smproperty.dataarrayselection().
    # This method **must** return a `vtkDataArraySelection` instance.
    @smproperty.dataarrayselection(name="Arrays")
    def GetDataArraySelection(self):
        return self._get_array_selection()


    @smproperty.intvector(name="N plane", default_values=2)
    def SetNPlane(self, val):
        if val != self.__n_plane:
            self.__n_plane = val
            self.Modified()

    @smproperty.doublevector(name="Phi Range", default_values=[0.0, 90.0])
    @smdomain.doublerange(min=0, max=360.0)
    def SetPhiRange(self, val, val2):
        if val != self.__phi_start:
            self.__phi_start = val
            self.Modified()
        if val2 != self.__phi_end:
            self.__phi_end = val2
            self.Modified()


    def RequestInformation(self, request, inInfoVec, outInfoVec):
        #executive = self.GetExecutive()
        #outInfo = executive.GetOutputInformation(0)

        executive = self.GetExecutive()
        outInfo = outInfoVec.GetInformationObject(0)

        self._timesteps = self._get_timesteps()

        outInfo.Remove(executive.TIME_STEPS())
        for i in range(len(self._filenames)):
            outInfo.Append(executive.TIME_STEPS(), self._timesteps[i])

        # Remove time range info
        outInfo.Remove(executive.TIME_RANGE())
        outInfo.Append(executive.TIME_RANGE(), self._timesteps[0])
        outInfo.Append(executive.TIME_RANGE(), self._timesteps[-1])

        if self._fields is None:
            self._fields = jorek.fields()

        for i in range(self._number_of_values):
            name = self._fields.var_names[i]
            self._arrayselection.AddArray(name)

        return 1

    def RequestData(self, request, inInfoVec, outInfoVec):
        from vtkmodules.vtkCommonDataModel import vtkUnstructuredGrid

        def GetUpdateTimestep(algorithm):
            """Returns the requested time value, or None if not present"""
            executive = algorithm.GetExecutive()
            outInfo = executive.GetOutputInformation(0)
            return outInfo.Get(executive.UPDATE_TIME_STEP()) \
                    if outInfo.Has(executive.UPDATE_TIME_STEP()) else None
        # Get the current timestep
        req_time = GetUpdateTimestep(self)

        if self._get_timesteps is None:
            self._timesteps = self._get_timesteps()

        # 4 possibilities here:
        # After last step: return last file
        # Before first step: return first file
        # Very close match: return that file
        # Between 2 files: interpolate
        interp = False
        interp_fraction = 0
        try:
            index = np.isclose(self._timesteps, req_time).tolist().index(True)
            # We have a match, read and return data for this filename
        except ValueError:
            # Check for other 3 cases
            if req_time > self._timesteps[-1]:
                # After last file
                index = len(self._filenames)-1
            elif req_time < self._timesteps[0]:
                # Before first file
                index = 0
            else:
                interp = True
                index = (self._timesteps > req_time).tolist().index(True)
                interp_fraction = (req_time - self._timesteps[index-1])/(self._timesteps[index] - self._timesteps[index-1])
                if (interp_fraction > 1 or interp_fraction < 0):
                    print("ERROR in interp_fraction", interp_fraction, self._timesteps, index)
     
        to_read = []
        if self._fields is None:
            self._fields = jorek.fields()

        for i in range(self._number_of_values):
            name = self._fields.var_names[i]
            if self._arrayselection.ArrayIsEnabled(name):
                to_read.append(i)  

        output = vtkUnstructuredGrid.GetData(outInfoVec)

        if (self.__smooth_time_interpolation and interp):
            self._fields.read(self._filenames[index], variables=to_read, file_prev=self._filenames[index-1],
               interp_fraction=interp_fraction)
        else:
            self._fields.read(self._filenames[index], variables=to_read)

        self._fields.to_vtk(n_sub=self.__number_of_subdivisions, n_plane=self.__n_plane, 
            phi=[self.__phi_start, self.__phi_end], without_n0_mode=self.__exclude_n0_mode,
            output=output, quadratic=self.__quadratic, bezier=self.__bezier,
             force_remake_grid=True)
        return 1
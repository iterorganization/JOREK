# script-version: 2.0
# Catalyst state generated using paraview version 5.11.0

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# ----------------------------------------------------------------
# setup views used in the visualization
# ----------------------------------------------------------------

# get the material library
materialLibrary1 = GetMaterialLibrary()

# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.ViewSize = [928, 789]
renderView1.InteractionMode = '2D'
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.CenterOfRotation = [10.0, 0.0, 0.0]
renderView1.StereoType = 'Crystal Eyes'
renderView1.CameraPosition = [10.0, 0.0, 5.463685892985357]
renderView1.CameraFocalPoint = [10.0, 0.0, 0.0]
renderView1.CameraFocalDisk = 1.0
renderView1.CameraParallelScale = 1.4141059655625834
renderView1.BackEnd = 'OSPRay raycaster'
renderView1.OSPRayMaterialLibrary = materialLibrary1

# init the 'GridAxes3DActor' selected for 'AxesGrid'
renderView1.AxesGrid.Visibility = 1

SetActiveView(None)

# ----------------------------------------------------------------
# setup view layouts
# ----------------------------------------------------------------

# create new layout object 'Layout #1'
layout1 = CreateLayout(name='Layout #1')
layout1.AssignView(0, renderView1)
layout1.SetSize(928, 789)

# ----------------------------------------------------------------
# restore active view
SetActiveView(renderView1)
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# setup the data processing pipelines
# ----------------------------------------------------------------

# create a new 'XML Partitioned Dataset Reader'
grid = XMLPartitionedDatasetReader(registrationName='grid', FileName=['/home/mr618/rds/rds-intel-vis-0hKvbqKpczQ/mr618/jorek/model199/test_intear_catalyst/datasets/grid_000000.vtpd', '/home/mr618/rds/rds-intel-vis-0hKvbqKpczQ/mr618/jorek/model199/test_intear_catalyst/datasets/grid_000001.vtpd', '/home/mr618/rds/rds-intel-vis-0hKvbqKpczQ/mr618/jorek/model199/test_intear_catalyst/datasets/grid_000002.vtpd', '/home/mr618/rds/rds-intel-vis-0hKvbqKpczQ/mr618/jorek/model199/test_intear_catalyst/datasets/grid_000003.vtpd', '/home/mr618/rds/rds-intel-vis-0hKvbqKpczQ/mr618/jorek/model199/test_intear_catalyst/datasets/grid_000004.vtpd', '/home/mr618/rds/rds-intel-vis-0hKvbqKpczQ/mr618/jorek/model199/test_intear_catalyst/datasets/grid_000005.vtpd', '/home/mr618/rds/rds-intel-vis-0hKvbqKpczQ/mr618/jorek/model199/test_intear_catalyst/datasets/grid_000006.vtpd', '/home/mr618/rds/rds-intel-vis-0hKvbqKpczQ/mr618/jorek/model199/test_intear_catalyst/datasets/grid_000007.vtpd', '/home/mr618/rds/rds-intel-vis-0hKvbqKpczQ/mr618/jorek/model199/test_intear_catalyst/datasets/grid_000008.vtpd', '/home/mr618/rds/rds-intel-vis-0hKvbqKpczQ/mr618/jorek/model199/test_intear_catalyst/datasets/grid_000009.vtpd', '/home/mr618/rds/rds-intel-vis-0hKvbqKpczQ/mr618/jorek/model199/test_intear_catalyst/datasets/grid_000010.vtpd', '/home/mr618/rds/rds-intel-vis-0hKvbqKpczQ/mr618/jorek/model199/test_intear_catalyst/datasets/grid_000011.vtpd', '/home/mr618/rds/rds-intel-vis-0hKvbqKpczQ/mr618/jorek/model199/test_intear_catalyst/datasets/grid_000012.vtpd', '/home/mr618/rds/rds-intel-vis-0hKvbqKpczQ/mr618/jorek/model199/test_intear_catalyst/datasets/grid_000013.vtpd', '/home/mr618/rds/rds-intel-vis-0hKvbqKpczQ/mr618/jorek/model199/test_intear_catalyst/datasets/grid_000014.vtpd', '/home/mr618/rds/rds-intel-vis-0hKvbqKpczQ/mr618/jorek/model199/test_intear_catalyst/datasets/grid_000015.vtpd', '/home/mr618/rds/rds-intel-vis-0hKvbqKpczQ/mr618/jorek/model199/test_intear_catalyst/datasets/grid_000016.vtpd', '/home/mr618/rds/rds-intel-vis-0hKvbqKpczQ/mr618/jorek/model199/test_intear_catalyst/datasets/grid_000017.vtpd', '/home/mr618/rds/rds-intel-vis-0hKvbqKpczQ/mr618/jorek/model199/test_intear_catalyst/datasets/grid_000018.vtpd', '/home/mr618/rds/rds-intel-vis-0hKvbqKpczQ/mr618/jorek/model199/test_intear_catalyst/datasets/grid_000019.vtpd', '/home/mr618/rds/rds-intel-vis-0hKvbqKpczQ/mr618/jorek/model199/test_intear_catalyst/datasets/grid_000020.vtpd', '/home/mr618/rds/rds-intel-vis-0hKvbqKpczQ/mr618/jorek/model199/test_intear_catalyst/datasets/grid_000021.vtpd', '/home/mr618/rds/rds-intel-vis-0hKvbqKpczQ/mr618/jorek/model199/test_intear_catalyst/datasets/grid_000022.vtpd', '/home/mr618/rds/rds-intel-vis-0hKvbqKpczQ/mr618/jorek/model199/test_intear_catalyst/datasets/grid_000023.vtpd', '/home/mr618/rds/rds-intel-vis-0hKvbqKpczQ/mr618/jorek/model199/test_intear_catalyst/datasets/grid_000024.vtpd', '/home/mr618/rds/rds-intel-vis-0hKvbqKpczQ/mr618/jorek/model199/test_intear_catalyst/datasets/grid_000025.vtpd', '/home/mr618/rds/rds-intel-vis-0hKvbqKpczQ/mr618/jorek/model199/test_intear_catalyst/datasets/grid_000026.vtpd', '/home/mr618/rds/rds-intel-vis-0hKvbqKpczQ/mr618/jorek/model199/test_intear_catalyst/datasets/grid_000027.vtpd', '/home/mr618/rds/rds-intel-vis-0hKvbqKpczQ/mr618/jorek/model199/test_intear_catalyst/datasets/grid_000028.vtpd', '/home/mr618/rds/rds-intel-vis-0hKvbqKpczQ/mr618/jorek/model199/test_intear_catalyst/datasets/grid_000029.vtpd', '/home/mr618/rds/rds-intel-vis-0hKvbqKpczQ/mr618/jorek/model199/test_intear_catalyst/datasets/grid_000030.vtpd', '/home/mr618/rds/rds-intel-vis-0hKvbqKpczQ/mr618/jorek/model199/test_intear_catalyst/datasets/grid_000031.vtpd', '/home/mr618/rds/rds-intel-vis-0hKvbqKpczQ/mr618/jorek/model199/test_intear_catalyst/datasets/grid_000032.vtpd', '/home/mr618/rds/rds-intel-vis-0hKvbqKpczQ/mr618/jorek/model199/test_intear_catalyst/datasets/grid_000033.vtpd', '/home/mr618/rds/rds-intel-vis-0hKvbqKpczQ/mr618/jorek/model199/test_intear_catalyst/datasets/grid_000034.vtpd', '/home/mr618/rds/rds-intel-vis-0hKvbqKpczQ/mr618/jorek/model199/test_intear_catalyst/datasets/grid_000035.vtpd', '/home/mr618/rds/rds-intel-vis-0hKvbqKpczQ/mr618/jorek/model199/test_intear_catalyst/datasets/grid_000036.vtpd', '/home/mr618/rds/rds-intel-vis-0hKvbqKpczQ/mr618/jorek/model199/test_intear_catalyst/datasets/grid_000037.vtpd', '/home/mr618/rds/rds-intel-vis-0hKvbqKpczQ/mr618/jorek/model199/test_intear_catalyst/datasets/grid_000038.vtpd', '/home/mr618/rds/rds-intel-vis-0hKvbqKpczQ/mr618/jorek/model199/test_intear_catalyst/datasets/grid_000039.vtpd', '/home/mr618/rds/rds-intel-vis-0hKvbqKpczQ/mr618/jorek/model199/test_intear_catalyst/datasets/grid_000040.vtpd', '/home/mr618/rds/rds-intel-vis-0hKvbqKpczQ/mr618/jorek/model199/test_intear_catalyst/datasets/grid_000041.vtpd', '/home/mr618/rds/rds-intel-vis-0hKvbqKpczQ/mr618/jorek/model199/test_intear_catalyst/datasets/grid_000042.vtpd', '/home/mr618/rds/rds-intel-vis-0hKvbqKpczQ/mr618/jorek/model199/test_intear_catalyst/datasets/grid_000043.vtpd', '/home/mr618/rds/rds-intel-vis-0hKvbqKpczQ/mr618/jorek/model199/test_intear_catalyst/datasets/grid_000044.vtpd', '/home/mr618/rds/rds-intel-vis-0hKvbqKpczQ/mr618/jorek/model199/test_intear_catalyst/datasets/grid_000045.vtpd', '/home/mr618/rds/rds-intel-vis-0hKvbqKpczQ/mr618/jorek/model199/test_intear_catalyst/datasets/grid_000046.vtpd', '/home/mr618/rds/rds-intel-vis-0hKvbqKpczQ/mr618/jorek/model199/test_intear_catalyst/datasets/grid_000047.vtpd', '/home/mr618/rds/rds-intel-vis-0hKvbqKpczQ/mr618/jorek/model199/test_intear_catalyst/datasets/grid_000048.vtpd', '/home/mr618/rds/rds-intel-vis-0hKvbqKpczQ/mr618/jorek/model199/test_intear_catalyst/datasets/grid_000049.vtpd', '/home/mr618/rds/rds-intel-vis-0hKvbqKpczQ/mr618/jorek/model199/test_intear_catalyst/datasets/grid_000050.vtpd'])
grid.TimeArray = 'timestep'

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from grid
gridDisplay = Show(grid, renderView1, 'UnstructuredGridRepresentation')

# get 2D transfer function for 'rho'
rhoTF2D = GetTransferFunction2D('rho')

# get color transfer function/color map for 'rho'
rhoLUT = GetColorTransferFunction('rho')
rhoLUT.TransferFunction2D = rhoTF2D
rhoLUT.RGBPoints = [0.9999920725822449, 4.05432e-07, 0.0, 5.90122e-06, 0.9999932431688905, 0.0, 0.120401, 0.302675, 0.9999944137555361, 0.0, 0.216583, 0.524574, 0.9999955843515098, 0.0552475, 0.345025, 0.6595, 0.9999967549381554, 0.128047, 0.492588, 0.720288, 0.999997925524801, 0.188955, 0.641309, 0.792092, 0.9999990961114467, 0.327673, 0.784935, 0.873434, 1.0000002666980923, 0.60824, 0.892164, 0.935547, 1.0000014372900548, 0.881371, 0.912178, 0.818099, 1.0000026078807116, 0.951407, 0.835621, 0.449279, 1.0000037784673572, 0.904481, 0.690489, 0.0, 1.0000049490540028, 0.85407, 0.510864, 0.0, 1.0000061196406484, 0.777093, 0.33018, 0.00088199, 1.0000072902366222, 0.672862, 0.139087, 0.00269398, 1.0000084608232678, 0.508815, 0.0, 0.0, 1.0000096314099134, 0.299417, 0.000366289, 0.000547829, 1.0000107288360596, 0.0157519, 0.00332021, 4.55569e-08]
rhoLUT.ColorSpace = 'Lab'
rhoLUT.NanColor = [0.25, 0.0, 0.0]
rhoLUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'rho'
rhoPWF = GetOpacityTransferFunction('rho')
rhoPWF.Points = [0.9999920725822449, 0.0, 0.5, 0.0, 1.0000107288360596, 1.0, 0.5, 0.0]
rhoPWF.ScalarRangeInitialized = 1

# trace defaults for the display properties.
gridDisplay.Representation = 'Surface'
gridDisplay.ColorArrayName = ['POINTS', 'rho']
gridDisplay.LookupTable = rhoLUT
gridDisplay.SelectTCoordArray = 'None'
gridDisplay.SelectNormalArray = 'None'
gridDisplay.SelectTangentArray = 'None'
gridDisplay.OSPRayScaleArray = 'Psi'
gridDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
gridDisplay.SelectOrientationVectors = 'None'
gridDisplay.ScaleFactor = 0.2
gridDisplay.SelectScaleArray = 'None'
gridDisplay.GlyphType = 'Arrow'
gridDisplay.GlyphTableIndexArray = 'None'
gridDisplay.GaussianRadius = 0.01
gridDisplay.SetScaleArray = ['POINTS', 'Psi']
gridDisplay.ScaleTransferFunction = 'PiecewiseFunction'
gridDisplay.OpacityArray = ['POINTS', 'Psi']
gridDisplay.OpacityTransferFunction = 'PiecewiseFunction'
gridDisplay.DataAxesGrid = 'GridAxesRepresentation'
gridDisplay.PolarAxes = 'PolarAxesRepresentation'
gridDisplay.ScalarOpacityFunction = rhoPWF
gridDisplay.ScalarOpacityUnitDistance = 0.1302402729886959
gridDisplay.OpacityArrayName = ['POINTS', 'Psi']
gridDisplay.SelectInputVectors = [None, '']
gridDisplay.WriteLog = ''

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
gridDisplay.ScaleTransferFunction.Points = [-0.20266033709049225, 0.0, 0.5, 0.0, -2.015789535292356e-09, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
gridDisplay.OpacityTransferFunction.Points = [-0.20266033709049225, 0.0, 0.5, 0.0, -2.015789535292356e-09, 1.0, 0.5, 0.0]

# setup the color legend parameters for each legend in this view

# get color legend/bar for rhoLUT in view renderView1
rhoLUTColorBar = GetScalarBar(rhoLUT, renderView1)
rhoLUTColorBar.Title = 'rho'
rhoLUTColorBar.ComponentTitle = ''

# set color bar visibility
rhoLUTColorBar.Visibility = 1

# show color legend
gridDisplay.SetScalarBarVisibility(renderView1, True)

# ----------------------------------------------------------------
# setup color maps and opacity mapes used in the visualization
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# setup extractors
# ----------------------------------------------------------------

# create extractor
pNG1 = CreateExtractor('PNG', renderView1, registrationName='PNG1')
# trace defaults for the extractor.
pNG1.Trigger = 'TimeStep'

# init the 'PNG' selected for 'Writer'
pNG1.Writer.FileName = 'rho_{timestep:06d}{camera}.png'
pNG1.Writer.ImageResolution = [1920, 1080]
pNG1.Writer.Format = 'PNG'
pNG1.Writer.ResetDisplay = 1

# ----------------------------------------------------------------
# restore active source
SetActiveSource(pNG1)
# ----------------------------------------------------------------

# ------------------------------------------------------------------------------
# Catalyst options
from paraview import catalyst
options = catalyst.Options()
options.GlobalTrigger = 'TimeStep'
options.CatalystLiveTrigger = 'TimeStep'

# ------------------------------------------------------------------------------
if __name__ == '__main__':
    from paraview.simple import SaveExtractsUsingCatalystOptions
    # Code for non in-situ environments; if executing in post-processing
    # i.e. non-Catalyst mode, let's generate extracts using Catalyst options
    SaveExtractsUsingCatalystOptions(options)

---
title: "util/jorek_read_h5.py"
nav_order: 36
render_with_liquid: false
parent: "General"
---

## util/jorek_read_h5.py
`jorek_read_h5.py` is a python script that makes **.vtk** object from jorek-output data. This object is then called by `test_jorek_read_h5.py`, which writes this data into **.vtk** or .**vtu** file. To properly run python script we must understand basic parameters and writer properties in *test* file.

### util/test_jorek_read.h5.py

First, we must choose variables and jorek file with desired data.  
```
f.read('/tmp/jorek_restart.h5', variables=[2, 5])
```
To select variables we specify their number in list of desired variables, where *0 = psi, 1 = u, 2 = j, 3 = w, 4 = rho, 5 = T, 6 = v_par*.

Secondly, we decide how geometry of our grid will look like. More elaborate description of grid's cells is presented in two following sections, here just input parameters are explained.  
```
grid = f.to_vtk(n_sub=3, phi=[0,180], n_plane=3, quadratic=False, bezier=False)
```
**n_sub** or number of subdivisions specifies how many points are to be contained in each edge of a previous on jorek cell in RZ plane (min = 2). Number of jorek cells is fixed by jorek file, here we just "increase resolution".  
**phi** is an interval of grid's polar angle in degrees (min = 0, max = 360). If start and end of interval are the same, grid is 2D.  
**n_plane** tells us how many RZ planes are used in our grid in toroidal coordinate, between which 3D cells can be made. If number of planes equals "1" our grid is 2D.  
**quadratic** and **bezier** are two bool parameters which describe properties of cells (explained in geometry).

Lastly, we choose format of our output file with `writer=vtk.vtkXMLUnstructuredGridWriter()` for .vtu  or `writer=vtk.vtkUnstructuredGridWriter()` for .vtk, if desired change encodint to "ascii" with `writer.SetDataModeToAscii()` and choose file name with `writer.SetFileName('jorek_2d.vtu')`.

### Basic geometry
Basic geometry is used when `bezier=False`. Cells used here for grid use only 4 points in 2D and 8 points in 3D which are connected with straight lines, so their vtk types are [quad](https/*vtk.org/doc/nightly/html/classvtkquad.html.md) and [hexahedron](https/*vtk.org/doc/nightly/html/classvtkhexahedron.html.md).  
If `quadratic=True` cells are a bit more complex. For 2D [quadratic quad](https/*vtk.org/doc/nightly/html/classvtkquadraticquad.html.md), which has 4 nodes in corners and 4 mid-edge nodes. In 3D the same principle of mid-edge nodes is used for [quadratic hexahedron](https/*vtk.org/doc/nightly/html/classvtkquadratichexahedron.html.md) which gives us 20 nodes. Mid-edge points are calculated in the same way as for higher subdivisions, with Bezier interpolation explained in [Daan Van Vugt thesis](https://research.tue.nl/en/publications/nonlinear-coupled-mhd-kinetic-particle-simulations-of-heavy-impur).

### Bezier geometry
Bezier geometry is the most complex, but it enables us to increase nonlinear subdivision level for same number of points within visualisation platform (such as **paraview**), because it uses Bezier polynomials to interpolate more points from initial control points.  
Calculation of control points in RZ plane used here is explained in *thesis*. For 2D these control points are used to create [bezier quadrilateral](https/*vtk.org/doc/nightly/html/classvtkbezierquadrilateral.html.md)s. As a result each cell is made of 4x4 nodes. Noding used for 2D and for 3D cells is explained in [this pdf](https/*staging.coreform.com/papers/implementation-of-rational-bezier-cells-into-vtk-report.pdf.md).

When we move to 3D, [bezier hexahedron](https/*vtk.org/doc/nightly/html/classvtkbezierhexahedron.html.md)s are made of 3 planes of 4x4 points, which gives us 48 nodes per cell. Because our cell is not 4x4x4 as in *pdf* we just skip nodes that are not there. Anisotropic structure of cell's shape is specified by `SetHigherOrderDegrees()` with shape `[3,3,2]` in *jorek_read_h5.py//. (Even though cells are made of 4x4x3 points, one number less must be given for each dimension.) Below pictures of example cell, it's top and side view and it's points are displayed and provided as ![ VTK file for ParaView](vtk_bezier_hexahedron.zip?linkonly).

![](points_and_surface.png?800)
![](side_view.png?430)
![](top_view.png?430)

Control points in RZ plane are still calculated with previous method in toroidal direction other method is used. Because circle arcs can only be displayed with rational quadratic Bezier curves, these control points must have rational weights. As previously explained cells are made of 3 planes, of which middle plane must have rational weights. If angle between first and third plane equals *2α* (second plane is at α) weight is calculated as `w=cos(α)` and new radius of points in second plane (still at the same polar angle) as `R=r/cos(α)`. We can see the position of control points with weights in top-view picture. See [this guide](https://ctan.math.illinois.edu/macros/latex/contrib/lapdf/rcircle.pdf) for rational quadratic Bezier curves.

For Bezier cells *.vtu* file format must be used because *.vtk* writer still has some problems with writing weights and higher order derees. VTK version 9 needs to be used for Bezier elements that are included also in SMITER 1.6.3 with patched *ParaView* version *5.8.0* or *ParaView* 5.9.1.

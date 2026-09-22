---
title: "Jorek2_fieldlines_vtk_newdiag"
nav_order: 26
render_with_liquid: false
parent: "General"
---

## Jorek2_fieldlines_vtk_newdiag

This diagnostics evaluates arbitrary physical expressions (using the new_diag framework) along field lines and writes this information to both vtk and txt files. Run it as
   jorek2_fieldlines_vtk_newdiag < jorek.in

It uses the file jorek_restart.h5 to load the simulation state, and a text file called stpts containing required information about the field lines along which physical expressions are to be to evaluated, and the physical expressions to be evaluated (see below). 
Output is generated in the files 
- field_lines.vtk is a binary file containing the expressions along field lines (using a coordinate system consistent with jorek2vtk) and can be read with visualization tools such as Paraview
- field_lines.txt is an ASCII file containing the same data (using the coordinate system in https://www.jorek.eu/wiki/doku.php?id=coordinates)

### Program settings
If a file stpts exists, it is read in to set the required parameters. The stpts file must have this format:
  # n_scalars
  3
  # expr  name
  Te    Te_eV
  ne    ne_m-3
  nimp  nimp_m-3
  # n_lines
  100
  # nr    R_start   Z_start   phi_start   n_turns
  1    3.003     1.606     6.0900      1
  2    2.988     1.558     6.0900      1
  3    2.990     1.560     6.0900      1
  4    2.946     1.576     6.0900      1
  5    2.989     1.571     6.0900      1
  ...
  100  2.997     1.620     6.0900      1
Here, n_scalars is the number of scalar expressions to be evaluated, expr and name are the expressions (as defined in the new_diag framework) and its name, n_lines are the number of field lines to be evaluated (in both backward and forward direction so the final number of computed field lines will be twice this number). Moreover, nr is the field line number, R_start, Z_start and phi_start its initial coordinate, and n_turns the number of toroidal turns along which the field line should be followed. As in the jorek2_poincare diagnostics, the user can either provide the initial coordinates of all the field lines, or specify groups of field lines where only the first and last lines are explicitly given, and their initial coordinates are assumed to be uniformily distributed between the first and the last.

### Output files
field_lines.txt contains the evaluated expressions along field lines using the coordinate system as in https://www.jorek.eu/wiki/doku.php?id=coordinates
The format of this file is like this:
            line               x               y               z             Phi    Te_eV           ne_m-3          nimp_m-3
               1  0.29471371E+01  0.57653373E+00  0.16059999E+01  0.60900002E+01  0.32168159E+01  0.25254070E+21  0.11506390E+21
               1  0.29454236E+01  0.59543729E+00  0.16052808E+01  0.60837169E+01  0.32214899E+01  0.25050444E+21  0.11412951E+21
               1  0.29435921E+01  0.61434269E+00  0.16045547E+01  0.60774336E+01  0.32317951E+01  0.24841616E+21  0.11315785E+21
               1  0.29416425E+01  0.63324934E+00  0.16038219E+01  0.60711503E+01  0.32482059E+01  0.24628058E+21  0.11214984E+21
               1  0.29395747E+01  0.65215647E+00  0.16030822E+01  0.60648675E+01  0.32712436E+01  0.24410277E+21  0.11110648E+21

field_lines.vtk is a binary file that can be read within Paraview. It contains the same data as field_line.txt, but the coordinate system is consistent with jorek2vtk. As example of the kind of plot that can be produced with this diagnostics is shown here.
![](jorek2_fieldlines_vtk.png)

Jorek plugins for ParaView
==========================

This project contains a ParaView plugin that reads JOREK HDF5 restart files (timed)
and shows them in the ParaView pipeline. 


Developer instructions
----------------------
Load the pre-requisite modules, create and activate a virtual environment and install
the project files. Virtual environment is required to install `h5py` package that is not provided by ParaView. Plugins are installed under the same virtual environment.

After installation, export the installation path of plugins to `PV_PLUGIN_PATH`. 
For ParaView to find the libraries `PYTHONPATH` is needed to find additional site-packages.

The following are instructions with ParaView built with EasyBuild that provides a Python module as a dependency used by ParaView and therefore can be extended with *venv*. 

~~~ bash
module load ParaView/5.10.0-intel-2020b-mpi
python -m venv --system-site-packages --clear --prompt paraview-jorek-plugin local 
source local/bin/activate
source install.sh local
# Either launch paraview and test the plugins where it can be reloaded if modified
paraview
# Or open up your IDE/code editor and begin development.
~~~

Upon modifying the source, run `source install.sh local` again.

Separation of `jorek` and `jorek/plugins` site-packages is required by ParaView
to prevent loading of non-plugin Python files pointed by `PV_PLUGIN_PATH`.

User instructions
-----------------

Subsequent use of installed local package can be with using virtual environment

~~~ bash
module load ParaView/5.10.0-intel-2020b-mpi
source local/bin/activate
source install.sh local
paraview
~~~

or by exporting install directory and plugin directly.

~~~ bash
module load ParaView/5.10.0-intel-2020b-mpi
export PYTHONPATH=${PWD}/local/lib64/python3.8/site-packages:${PYTHONPATH}
export PV_PLUGIN_PATH=${PWD}/local/lib64/python3.8/site-packages/jorek/plugins
~~~

Open `jorek*.h5` file(s), select Arrajs and set parameters before Apply. 
If using *Bezier* elements then 'Nonlinear Subdivision Level' that defaults to 1 can 
be increased up to 3 and Paraview will do the subdivision. The easiest to find this
parameter is to type 'non' in *Properties* Search.

Standalone ParaView use
-----------------------

If not having ParaView built with separate Python then one can just install
pre-built ParaView downloaded from https://www.paraview.org/download/ 
Versions from 5.10.1 include `h5py` module and just export to current and 
plugins directory is required.

~~~ bash
tar xzf ~/ParaView-5.10.1-MPI-Linux-Python3.9-x86_64.tar.gz
export PYTHONPATH=${PWD}
export PV_PLUGIN_PATH=${PWD}/jorek/plugins
./ParaView-5.10.1-MPI-Linux-Python3.9-x86_64/bin/paraview
~~~



TODO
----

Projections and particles in HDF5 are not supported and stored under 'old/'.
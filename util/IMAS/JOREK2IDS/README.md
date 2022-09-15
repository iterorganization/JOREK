# Writing timed IMAS MHD IDS from HDF5 files and generating MHD and Radiation VTK files

IMAS needs to be set up with

    module load IMAS # at least 3.37 required
    imasdb jorek
    imasdb

## Preparing python environment


	python3 -m venv local
	local/bin/python -m pip install --upgrade pip
	local/bin/python -m pip install h5py VTK pytz 
	local/bin/python jorekHDF5toIDS.py --help

Sample run can be saved with

    local/bin/python3 jorekHDF5toIDS.py --shot=303 --run=1 --user=${USER} --database=jorek --occurrence=0 jorek00000.h5


and created for visualisation with ParaView using

	local/bin/python3 IDS_to_VTK.py --shot=303 --run=1 --user=${USER} --database=jorek --occurrence=0
	paraview jorek..vtu


MHD IDS is written under `~/public/imasdb/jorek/3/0/` as 

    ids_3030001.characteristics  ids_3030001.datafile  ids_3030001.tree

## Generating MHD and Radiation VTK

    local/bin/python IDS_to_VTK.py --help
    module load ParaView
    paraview jorek..vtu

If you start typing sub in ParaView Advanced Properties pane, the 
`NonLinear Subdivision Level` can be increased from 1 to 3 
to get recomputted at finer FEM.
     
 

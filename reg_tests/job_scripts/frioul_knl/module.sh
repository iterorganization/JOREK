module load intel/17.1  intelmpi/2017.1.132
. /opt/software/common/intel/compilers_and_libraries_2017.1.132/linux/bin/ifortvars.sh intel64
export PASTIX_LIB=/home/glatu/jorek_tools/pastix_3184_knl_17_mt/install
export MUMPS_HOME=/home/glatu/jorek_tools/mumps
export HDF5_HOME=/home/glatu/jorek_tools/hdf5
export PATH=${HDF5_HOME}/bin:$PATH
export LD_LIBRARY_PATH=${HDF5_HOME}/lib:$LD_LIBRARY_PATH

# Author: Kos Leon
# Description: Install paraview_jorek_plugins on any system with pip.

import setuptools
from jorek._version import __version__

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="jorek_plugins",
    version =  __version__,
    author="Leon Kos",
    author_email="leon.kos@lecad.fs.uni-lj.si",
    description="Read JOREK HDF5 with ParaView plugins",
    long_description=long_description,
    url="https://www.jorek.eu/wiki/doku.php?id=paraview-plugin",
    classifiers=[
        "Programming Language :: Python :: 3",
        # "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    #packages=setuptools.find_packages(),
    packages=['jorek'],
    package_dir={'jorek':'jorek'},
    package_data={'jorek': ['plugins/*.py']},
    python_requires=">=3.8",
    install_requires=[
        "h5py >= 3.5.0",
        "numpy >= 1.19.4",
        "vtk >= 9.1.0"
    ],
)
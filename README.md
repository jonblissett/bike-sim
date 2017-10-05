This software enables the simulation of an electric motorcycle accounting for numerous electrical and mechanical 
effects. It is particularly suited to lap time simulation and system optimisation of such motorcycles.

## Installation

Use Python2.7, requires the following packages:

matplotlib, numpy, scipy, cython, os, itertools, time

Commonnly used python functions are compiled in C using Cython for speed. These .pyx files must first be compiled 
with Cython uing the following commands:

python2.7 Cythonise_bike_fcn.py build_ext --inplace
python2.7 Cythonise_interp_fcn.py build_ext --inplace 
python2.7 Cythonise_losses_fcn.py build_ext --inplace 
python2.7 Cythonise_pmsm_fcn.py build_ext --inplace 

## Usage

Saved data from the simulation will be storted in the output directory './data_export'

Filenames beginning 'main' are a good place to start.

If running using parallel processing, in the working directory execute command 'ipcluster start -n ?' where ? 
specfies the number of cores. This is not nessecary for 'one off' simulations.

This software enables the simulation of an electric motorcycle accounting for numerous electrical and mechanical 
effects. It is particularly suited to lap time simulation and system optimisation of such motorcycles.

## Prerequisites

Use Python2.7, requires the following packages:

matplotlib, numpy, scipy, cython, ipyparallel

Commonnly used python functions are compiled in C using Cython for speed. These .pyx files must first be compiled 
with Cython by executing the following commands, or the shell script "compile_cython.sh":

```sh
python2.7 Cythonise_bike_fcn.py build_ext --inplace
python2.7 Cythonise_interp_fcn.py build_ext --inplace 
python2.7 Cythonise_losses_fcn.py build_ext --inplace 
python2.7 Cythonise_pmsm_fcn.py build_ext --inplace 
```

## Usage

Saved data from the simulation will be storted in the output directory './data_export'

Filenames beginning 'main' are a good place to start.

If running using parallel processing, in the working directory execute command `ipcluster start -n ?` where ? 
specfies the number of cores. This is not required for single simulations.

## Windows installation

If you must use Windows, the following steps are required:

0. Install [Python2.7](https://www.python.org/downloads/windows/).

1. To edit and debug the code you will want an IDE. I use 
[PyCharm](https://www.jetbrains.com/pycharm/). 

2. You should now install the prerequesit packages, likely within the IDE. You may find it easier to install Numpy 
and SciPy from [here](http://www.lfd.uci.edu/~gohlke/pythonlibs/#numpy) and 
[here](http://www.lfd.uci.edu/~gohlke/pythonlibs/#scipy).

3. Install a suitable C compiler, such as [Microsoft Visual C++ Compiler for Python 2.7
](https://www.microsoft.com/en-gb/download/details.aspx?id=44266).

4. Compile the C code using Cython. First ensure your IDE is not running. Typically:

```
C:\Python27\python.exe Cythonise_bike_fcn.py build_ext --inplace --compiler=msvc
C:\Python27\python.exe Cythonise_interp_fcn.py build_ext --inplace --compiler=msvc
C:\Python27\python.exe Cythonise_losses_fcn.py build_ext --inplace --compiler=msvc
C:\Python27\python.exe Cythonise_pmsm_fcn.py build_ext --inplace --compiler=msvc
```

If you see the error 'Unable to find vcvarsall.bat', open the command prompt for the C++ Compiler and execute the following 
commands first:

```
SET DISTUTILS_USE_SDK=1
SET MSSdk=1
```

echo "Compiling four files with Cython"
python2.7 Cythonise_bike_fcn.py build_ext --inplace
python2.7 Cythonise_interp_fcn.py build_ext --inplace 
python2.7 Cythonise_losses_fcn.py build_ext --inplace 
python2.7 Cythonise_pmsm_fcn.py build_ext --inplace 
echo "Ready to run"

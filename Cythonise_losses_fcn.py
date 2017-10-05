from distutils.core import setup
from Cython.Build import cythonize
import numpy

setup(
  name = 'C_losses',
  ext_modules = cythonize("C_losses.pyx"),
  include_dirs=[numpy.get_include()]

)
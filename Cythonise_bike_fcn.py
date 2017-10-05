from distutils.core import setup
from Cython.Build import cythonize
import numpy

setup(
  name = 'motorbike_C',
  ext_modules = cythonize("C_motorbike.pyx"),
  include_dirs=[numpy.get_include()]

)
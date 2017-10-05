from distutils.core import setup
from Cython.Build import cythonize
import numpy

setup(
  name = 'C_pmsm',
  ext_modules = cythonize("C_pmsm.pyx"),
  include_dirs=[numpy.get_include()]

)
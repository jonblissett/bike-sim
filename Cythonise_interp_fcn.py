from distutils.core import setup
from Cython.Build import cythonize
import numpy

setup(
  name='C_interp',
  ext_modules=cythonize("C_interp.pyx"),
  include_dirs=[numpy.get_include()]
)
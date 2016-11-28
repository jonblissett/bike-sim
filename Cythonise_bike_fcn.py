from distutils.core import setup
from Cython.Build import cythonize

setup(
  name = 'motorbike_C',
  ext_modules = cythonize("C_motorbike.pyx"),
)
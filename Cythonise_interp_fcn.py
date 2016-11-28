from distutils.core import setup
from Cython.Build import cythonize

setup(
  name='C_interp',
  ext_modules=cythonize("C_interp.pyx"),
)
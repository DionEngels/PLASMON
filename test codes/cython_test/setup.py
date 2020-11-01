"""from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
from Cython.Build import cythonize
import numpy

#  ext_modules = [Extension("cython_test", ["cython_test.pyx"], extra_compile_args=["-ffast-math"])]
#  setup(name="cython_test", cmdclass={"build_ext": build_ext}, ext_modules=ext_modules)

setup(name="cython_test", ext_modules=cythonize('cython_test.pyx'), include_dirs=[numpy.get_include()])
"""
"""

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy
ext_modules = [Extension("cython_test", sources=["cython_test.pyx"], include_dirs=[numpy.get_include()], extra_compile_args=["-w"])]
setup(name='cython_test', cmdclass={'build_ext': build_ext}, ext_modules=ext_modules)

"""
from setuptools import setup
from Cython.Build import cythonize

setup(
    ext_modules = cythonize("helloworld.pyx")
)
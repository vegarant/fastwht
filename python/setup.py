from distutils.core import setup, Extension
import numpy
import os

os.environ['CC'] = 'g++';
setup(name='hadamardKernel_test', version='1.0', ext_modules =[Extension('_hadamardKernel',
 ['../hadamard.cpp', 'hadamardKernel.i'], include_dirs = [numpy.get_include(),'.'], extra_compile_args=['-std=c++11', '-O3'] )])


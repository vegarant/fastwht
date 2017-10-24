#! /bin/sh

rm -f hadamardKernel_wrap* _hadamardKernel.so hadamardKernel.py*
rm -fR build
rm -f *.pyc
swig -c++ -python hadamardKernel.i

echo "Running python setup.py build_ext --inplace"

python setup.py build_ext --inplace


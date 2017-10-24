#! /bin/sh

rm -f hadamardKernel_wrap* _hadamardKernel.so hadamardKernel.py*
rm -fR build

swig -c++ -python hadamardKernel.i

echo "Running python2 setup.py build_ext --inplace"

python2 setup.py build_ext --inplace


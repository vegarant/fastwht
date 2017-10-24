#! /bin/sh

rm -f hadamardKernel_wrap* _hadamardKernel.so hadamardKernel.py*
rm -fR build
rm -f *.pyc
swig -c++ -python hadamardKernel.i

echo "Running python setup.py build_ext --inplace"

PYTHON_VERSION_STR=`python -V 2> /dev/stdout` 
SUBSTRING=$(echo $PYTHON_VERSION_STR | cut -d' ' -f 2)
VERSION=$(echo $SUBSTRING | cut -d'.' -f 1)

if [ "$VERSION" -eq "2" ]
then
    python3 setup.py build_ext --inplace
else
    python setup.py build_ext --inplace
fi




# Walsh-Hadamard Transform

## Introduction

Hadamard is a python 2 module for the *Walsh-Hadamard transform*. It preforms
this transform in O(N log N) operations, using one of the tree orderings
*sequency*, *hadamard* and *dyadic*. It does also provide functionality to 
produce the values to the functions *WAL* and *PAL*.   


## Usage

```
>>> from hadamard import *
>>> from numpy import *
>>> 
>>> R = 2**3
>>> U = zeros([R,R])
>>> 
>>> for i in range(R):
...     x = zeros(R) 
...     x[i] = 1
...     U[:,i] = fastwht(x);
... 
>>> print U*R
[[ 1.  1.  1.  1.  1.  1.  1.  1.]
 [ 1.  1.  1.  1. -1. -1. -1. -1.]
 [ 1.  1. -1. -1. -1. -1.  1.  1.]
 [ 1.  1. -1. -1.  1.  1. -1. -1.]
 [ 1. -1. -1.  1.  1. -1. -1.  1.]
 [ 1. -1. -1.  1. -1.  1.  1. -1.]
 [ 1. -1.  1. -1. -1.  1. -1.  1.]
 [ 1. -1.  1. -1.  1. -1.  1. -1.]]
>>> 
```

## Installation

The hadamard module imports its core modules from code written in C++. The 
wrapper code used to interface python with C++ is auto generated using Swig. 
and placed in the module *hadamardKernel*. The Hadamard module import all its 
functionality from this kernel module. 

The easiest way to install this module is to download the precompiled binaries [Download](http://folk.uio.no/vegarant/fastwht_matlab.zip). If this fails, one 
can try to compile it using the setup.sh script. 

In order to make the python code callable from other directories add the
directory to your python path.
```export PYTHONPATH=$PYTHONPATH:/path/to/directory/fastwht/python```




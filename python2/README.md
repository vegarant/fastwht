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
>>> N = 2**3
>>> U = zeros([N,N])
>>>
>>> for i in range(N):
...     x = zeros(N)
...     x[i] = 1
...     U[:,i] = fastwht(x);
... 
>>> print N*U
[[ 1.  1.  1.  1.  1.  1.  1.  1.]
 [ 1.  1.  1.  1. -1. -1. -1. -1.]
 [ 1.  1. -1. -1. -1. -1.  1.  1.]
 [ 1.  1. -1. -1.  1.  1. -1. -1.]
 [ 1. -1. -1.  1.  1. -1. -1.  1.]
 [ 1. -1. -1.  1. -1.  1.  1. -1.]
 [ 1. -1.  1. -1. -1.  1. -1.  1.]
 [ 1. -1.  1. -1.  1. -1.  1. -1.]]
>>>
>>> N = 4;
>>> 
>>> U_wal = zeros([N,N]);
>>> U_pal = zeros([N,N]);
>>> for n in range(N):
...     for t in range(N):
...         U_wal[n,t] = WAL(N, n, t);
...         U_pal[n,t] = PAL(N, n, t);
... 
>>> # Sequency ordered Hadamard matrix of size 4
... print U_wal;
[[ 1.  1.  1.  1.]
 [ 1.  1. -1. -1.]
 [ 1. -1. -1.  1.]
 [ 1. -1.  1. -1.]]
>>> # Paley ordered Hadamard matrix of size 4
... print U_pal;
[[ 1.  1.  1.  1.]
 [ 1.  1. -1. -1.]
 [ 1. -1.  1. -1.]
 [ 1. -1. -1.  1.]]
>>> 
```

## Installation

The Hadamard module imports its core modules from code written in C++. The
wrapper code used to interface python with C++ is auto generated using Swig. This
code creates the module *hadamardKernel*. The Hadamard module import all its
functionality from this kernel module.

The easiest way to install this module is to download the precompiled binaries [Download](http://folk.uio.no/vegarant/fastwht_python.zip). If this fails, one
can try to compile it using the setup.sh script.

In order to make the python code callable from other directories add the
directory to your python path.
```export PYTHONPATH=$PYTHONPATH:/path/to/directory/fastwht/python```




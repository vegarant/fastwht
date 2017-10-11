# Fastwht
## A fast Walsh-Hadamard transform for Matlab and Python 2

`fastwht` is an C++ implementation (using [fxt](http://www.jjj.de/fxt/)) of the
fast Walsh-Hadamard transform with bindings to Matlab and Python 2. The
algorithm perform the transform in-place in O( N log(N) ) operations. Current
version are order of magnitude faster than Matlab's own implementation `fwht`.
As Python's NumPy package does not implement this transform, no such comparison
have been made for the Python implementation. 

### Matlab

The Matlab interface supports the two functions `fastwht` and `wal`. The 
`fastwht` function is a faster version of Matlab's `fwht` function, and for its
usage we refer to Matlab's own 
[documentation](http://se.mathworks.com/help/signal/ref/fwht.html).
The [`wal` function](src/master/matlab/wal.md) is the function generating the matrix entries in a 2^n x 2^n 
sequency ordered Hadamard matrix. The current version of this function is 
vectorized, so that it can handle vector input. 

A test of performance for the `fastwht` function can be seen in the image
below. The test where performed on an Lenovo ThinkPad T440s with Intel Core
i7-4600U CPU and 8 GB of RAM. The computer where using the operating system
Arch Linux. The time measurements where performed using Matlab's `timeit()`
function.
![alt text](matlab/compare_performance.png =200x150)

### Python
The interface between C++ and Python 2 is auto generated using
[Swig](http://www.swig.org) version 3.0.8. On top of this python interface,
there have been built an extra layer of wrapper code, to support complex arrays
and matrix objects. None of swig's auto generated code is submitted in this
repository.

In addition to the `fastwht` function, one fins the `WAL` and `PAL` funcitons, but 
these have not yet been extended to support vectorized input. 

## Install
You may either pull the required code directly from this
repository and compile it yourself or you can download these precompiled binaries:
[Matlab](http://folk.uio.no/vegarant/fastwht_matlab.zip) and
[Python](http://folk.uio.no/vegarant/fastwht_python.zip).
Remember to update your Python and Matlab path after you have complied the
code. The current release have only been tested on one machine running Arch
Linux.

## License
The project is published under GNU General Public License version 3.

## Contact
For any questions concerning the code please contact Vegard Antun at
vegarant@math.uio.no

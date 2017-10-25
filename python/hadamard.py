# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 3 of the License.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# Copyright 2016 Vegard Antun
#

"""
Hadamard is a python module providing an implementation of the Walsh-Hadamard
transform.
"""
from hadamardKernel import *;
import sys
import numpy as np


def fastwht(x, N=0, order='sequency'):
    '''
    fastwht - Preform the fast Walsh-Hadamard transform on the array or matrix `x`.

    Input:

    x - NumPy array or matrix.
    N - The size of the transform. Must be a power of 2.
    order - One of the tree orderings, *sequency*, *hadamard* or *dyadic*

    Background:

    The Hadamard matrices are orthogonal matrices of size NxN  with N = 2^R,
    consisting of the elements {1/N, -1/N}. They can be ordered in tree different
    ways, known as *sequency*, *hadamard* and *dyadic* ordering. Due to the
    structure of these matrices their matrix product with the vector `x` can be
    computed in-place with O( N log N ) operations. This is how fastwht is
    implemented.

    Functionality:

    If the length of `x` is not a power of 2, it will be zero expanded with the
    appropriate number of zeros, such that the transform is defined. The returned
    array will necessarily be longer than the input array in such cases.

    If N is non-zero, the input array will be zero-extended or shrunken to have
    length N. The returned array will then have size N, while the input array
    will be returned untouched.

    The function is only defined for one-dimensional NumPy arrays and matrices
    of type numpy.float64 and numpy.complex64.
    '''

    if isinstance(x, (np.ndarray, np.matrix)):
        if (x.dtype == np.complex128):
            y_real = fastwht(x.real, N, order);
            y_imag = fastwht(x.imag, N, order);
            y = y_real + 1j*y_imag;
            if (isinstance(x, np.matrix)):
                y = np.matrix(y);
                if x.shape[0] == 1:
                    y.shape = (1, max(y.shape));
                else:
                    y.shape = (max(y.shape),1);

            return y;

        if (x.dtype != np.float64 and x.dtype != np.int64):
            raise ValueError('fastwht: Only numpy.int64, numpy.float64 and numpy.complex128 supported');

        if (len(x.shape) > 2):
            raise TypeError('Only one or two dimensional arrays/matrices');

    else:
        raise TypeError('fastwht: Input must be of type numpy.ndarray og numpy.matrix');

    if (order != 'sequency' and order != 'dyadic' and order != 'hadamard'):
        raise ValueError('fastwht: Unknown ordering');
    if (len(x.shape) > 2):
        raise TypeError("fastwht: x must be one or two dimensional");
    if (not isinstance(N, int)):
        raise TypeError("fastwht: N must be integer");
    if (N < 0):
        raise ValueError('fastwht: N must be non-negative integer');
    

    originalShape = x.shape;
    if(x.shape[0] != 1):
        originalLength = x.shape[0];
    else:
        originalLength = max(x.shape);

    # N has been provided by user 
    if (N != 0):
        if (not _is_power_of_2(N)):
            raise ValueError('fastwht: N must be a power of 2');

        if (N == originalLength):
            y = np.asarray(x.copy());
        elif ( N > originalLength ):
            if (len(x.shape) == 1):
                y = np.zeros(N);
                y[:len(x)] = np.asarray(x.copy());
            else:
                if (x.shape[0] == 1):
                    y = np.zeros([1,N]);
                    y[:, :x.shape[1]] = np.asarray(x.copy());
                else: 
                    y = np.zeros([N, x.shape[1]]);
                    y[:x.shape[0],:] = np.asarray(x.copy());
        elif ( N < originalLength ):
            if (len(x.shape) == 1):
                y = np.asarray(np.copy(x[:N]));
            else:
                if (x.shape[0] == 1):
                    y = np.asarray(np.copy(x[:,:N]));
                else:
                    y = np.asarray(np.copy(x[:N,:]));
    else:
        if (not _is_power_of_2(originalLength) ):
            N = _make_power_of_2(originalLength);
            if (len(x.shape) == 1):
                y = np.zeros(N);
                y[:len(x)] = np.asarray(x.copy());
            else:
                y = np.zeros([N, x.shape[1]]);
                y[:x.shape[0], :] = np.asarray(x.copy());
        else:
            N = originalLength;
            y = np.asarray(x.copy());

    # y is the correct shape and we can apply the transform
    if (len(y.shape) == 1):
        if (order == 'hadamard'):
            fwhtKernelOrdinary(y);
        elif (order == 'sequency'):
            fwhtKernelSequency(y);
        else:
            fwhtKernelPaley(y);
    else:
        if (y.shape[0] != 1):
            if (order == 'hadamard'):
                fwhtKernel2dOrdinary(y);
            elif (order == 'sequency'):
                fwhtKernel2dSequency(y);
               
            else: 
                fwhtKernel2dPaley(y);
        else: # shape is size 2 and shape[0] == 1
            y.shape = (max(y.shape), );
            if (order == 'hadamard'):
                fwhtKernelOrdinary(y);
            elif (order == 'sequency'):
                fwhtKernelSequency(y);
            else: 
                fwhtKernelPaley(y);
            y.shape = (1, y.shape[0]);
    
    if y.dtype == np.int64:
        y = y.astype('float64');
    y /= N;
    if (N == 0 or N == originalLength):
        y.shape = originalShape;
    if isinstance(x, np.matrix):
        y = np.matrix(y);

    return y;


def fastwht2(X, A=[0,0], order='sequency'):
    if (len(A) != 2):
        raise ValueError('A must contain two elements');
    Y = fastwht(X, A[0], order);
    Z = fastwht(Y.transpose(), A[1], order);
    return Z.transpose();


def WAL(N, n, t):
    """
    The sequency ordered Walsh-function.

    INPUT:

    N - Dimension of the sequency ordered Hadamard matrix. Must be N = 2**r for
        some positive integer r
    n - Walsh-function number.
    t - Input to the Walsh function.

    All inputs must be non-negative integers. This function is only defined for
    0 <= n,t < N. As the n'th order Walsh-function's input is only defined for
    input in the interval [0,1). The function WAL(N,n,t) can be interpreted as
    w_n(t/N).

    """

    if (not _is_power_of_2(N)) :
        raise ValueError('N must equal 2**r for some positive integer r');
    if (t < 0 or t > N):
        raise ValueError("Illegal t-value: Must be in the interval 0 <= t < N");
    if (n < 0 or n > N):
        raise ValueError("Illegal n-value: Must be in the interval 0 <= n < N");

    return WAL_kernel(N, n, t);


def PAL(N, n, t):
    """
    The Paley ordered Walsh-function.

    INPUT:

    N - Dimension of the Paley ordered Hadamard matrix. Must be N = 2**r for
        some positive integer r
    n - Walsh-function number.
    t - Input to the Walsh function.

    All inputs must be non-negative integers. This function is only defined for
    0 <= n,t < N. As the n'th order Walsh-function's input is only defined for
    input in the interval [0,1). The function PAL(N,n,t) can be interpreted as
    w_n(t/N).

    """

    if (not _is_power_of_2(N)) :
        raise ValueError('N must equal 2**r for some positive integer r');
    if (t < 0 or t > N):
        raise ValueError("Illegal t-value: Must be in the interval 0 <= t < N");
    if (n < 0 or n > N):
        raise ValueError("Illegal n-value: Must be in the interval 0 <= n < N");

    return PAL_kernel(N, n, t);


def _is_power_of_2(N):
    k = int(np.log2(N));
    return (2**k) == N;


def _make_power_of_2(N):
    k = int(np.log2(N));

    N_new = 2**k;
    if (N_new == N):
        return N;
    else:
        return 2*N_new;


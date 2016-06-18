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

from numpy import *;
from hadamard import *;

Had_seq8 = array([[1,  1,  1,  1,  1,  1,  1,  1], \
[1,  1,  1,  1, -1, -1, -1, -1], \
[1,  1, -1, -1, -1, -1,  1,  1], \
[1,  1, -1, -1,  1,  1, -1, -1], \
[1, -1, -1,  1,  1, -1, -1,  1], \
[1, -1, -1,  1, -1,  1,  1, -1], \
[1, -1,  1, -1, -1,  1, -1,  1], \
[1, -1,  1, -1,  1, -1,  1, -1]])

Had_seq16 = array([[1,  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],\
[1,  1, 1, 1, 1, 1, 1, 1,-1,-1,-1,-1,-1,-1,-1,-1],\
[1,  1, 1, 1,-1,-1,-1,-1,-1,-1,-1,-1, 1, 1, 1, 1],\
[1,  1, 1, 1,-1,-1,-1,-1, 1, 1, 1, 1,-1,-1,-1,-1],\
[1,  1,-1,-1,-1,-1, 1, 1, 1, 1,-1,-1,-1,-1, 1, 1],\
[1,  1,-1,-1,-1,-1, 1, 1,-1,-1, 1, 1, 1, 1,-1,-1],\
[1,  1,-1,-1, 1, 1,-1,-1,-1,-1, 1, 1,-1,-1, 1, 1],\
[1,  1,-1,-1, 1, 1,-1,-1, 1, 1,-1,-1, 1, 1,-1,-1],\
[1, -1,-1, 1, 1,-1,-1, 1, 1,-1,-1, 1, 1,-1,-1, 1],\
[1, -1,-1, 1, 1,-1,-1, 1,-1, 1, 1,-1,-1, 1, 1,-1],\
[1, -1,-1, 1,-1, 1, 1,-1,-1, 1, 1,-1, 1,-1,-1, 1],\
[1, -1,-1, 1,-1, 1, 1,-1, 1,-1,-1, 1,-1, 1, 1,-1],\
[1, -1, 1,-1,-1, 1,-1, 1, 1,-1, 1,-1,-1, 1,-1, 1],\
[1, -1, 1,-1,-1, 1,-1, 1,-1, 1,-1, 1, 1,-1, 1,-1],\
[1, -1, 1,-1, 1,-1, 1,-1,-1, 1,-1, 1,-1, 1,-1, 1],\
[1, -1, 1,-1, 1,-1, 1,-1, 1,-1, 1,-1, 1,-1, 1,-1]])

Had_pal8 = array([[1,  1,  1,  1,  1,  1,  1,  1],\
[1,  1,  1,  1, -1, -1, -1, -1],\
[1,  1, -1, -1,  1,  1, -1, -1],\
[1,  1, -1, -1, -1, -1,  1,  1],\
[1, -1,  1, -1,  1, -1,  1, -1],\
[1, -1,  1, -1, -1,  1, -1,  1],\
[1, -1, -1,  1,  1, -1, -1,  1],\
[1, -1, -1,  1, -1,  1,  1, -1]]);

Had_pal16 = array([[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],\
[1, 1, 1, 1, 1, 1, 1, 1,-1,-1,-1,-1,-1,-1,-1,-1],\
[1, 1, 1, 1,-1,-1,-1,-1, 1, 1, 1, 1,-1,-1,-1,-1],\
[1, 1, 1, 1,-1,-1,-1,-1,-1,-1,-1,-1, 1, 1, 1, 1],\
[1, 1,-1,-1, 1, 1,-1,-1, 1, 1,-1,-1, 1, 1,-1,-1],\
[1, 1,-1,-1, 1, 1,-1,-1,-1,-1, 1, 1,-1,-1, 1, 1],\
[1, 1,-1,-1,-1,-1, 1, 1, 1, 1,-1,-1,-1,-1, 1, 1],\
[1, 1,-1,-1,-1,-1, 1, 1,-1,-1, 1, 1, 1, 1,-1,-1],\
[1,-1, 1,-1, 1,-1, 1,-1, 1,-1, 1,-1, 1,-1, 1,-1],\
[1,-1, 1,-1, 1,-1, 1,-1,-1, 1,-1, 1,-1, 1,-1, 1],\
[1,-1, 1,-1,-1, 1,-1, 1, 1,-1, 1,-1,-1, 1,-1, 1],\
[1,-1, 1,-1,-1, 1,-1, 1,-1, 1,-1, 1, 1,-1, 1,-1],\
[1,-1,-1, 1, 1,-1,-1, 1, 1,-1,-1, 1, 1,-1,-1, 1],\
[1,-1,-1, 1, 1,-1,-1, 1,-1, 1, 1,-1,-1, 1, 1,-1],\
[1,-1,-1, 1,-1, 1, 1,-1, 1,-1,-1, 1,-1, 1, 1,-1],\
[1,-1,-1, 1,-1, 1, 1,-1,-1, 1, 1,-1, 1,-1,-1, 1]] )

Had_ord8 = array([[ 1,  1,  1,  1,  1,  1,  1,  1], \
[1, -1,  1, -1,  1, -1,  1, -1],\
[1,  1, -1, -1,  1,  1, -1, -1],\
[1, -1, -1,  1,  1, -1, -1,  1],\
[1,  1,  1,  1, -1, -1, -1, -1],\
[1, -1,  1, -1, -1,  1, -1,  1],\
[1,  1, -1, -1, -1, -1,  1,  1],\
[1, -1, -1,  1, -1,  1,  1, -1]]);

Had_ord16 = array([[1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1],\
[1, -1,  1, -1,  1, -1,  1, -1,  1, -1,  1, -1,  1, -1,  1, -1],\
[1,  1, -1, -1,  1,  1, -1, -1,  1,  1, -1, -1,  1,  1, -1, -1],\
[1, -1, -1,  1,  1, -1, -1,  1,  1, -1, -1,  1,  1, -1, -1,  1],\
[1,  1,  1,  1, -1, -1, -1, -1,  1,  1,  1,  1, -1, -1, -1, -1],\
[1, -1,  1, -1, -1,  1, -1,  1,  1, -1,  1, -1, -1,  1, -1,  1],\
[1,  1, -1, -1, -1, -1,  1,  1,  1,  1, -1, -1, -1, -1,  1,  1],\
[1, -1, -1,  1, -1,  1,  1, -1,  1, -1, -1,  1, -1,  1,  1, -1],\
[1,  1,  1,  1,  1,  1,  1,  1, -1, -1, -1, -1, -1, -1, -1, -1],\
[1, -1,  1, -1,  1, -1,  1, -1, -1,  1, -1,  1, -1,  1, -1,  1],\
[1,  1, -1, -1,  1,  1, -1, -1, -1, -1,  1,  1, -1, -1,  1,  1],\
[1, -1, -1,  1,  1, -1, -1,  1, -1,  1,  1, -1, -1,  1,  1, -1],\
[1,  1,  1,  1, -1, -1, -1, -1, -1, -1, -1, -1,  1,  1,  1,  1],\
[1, -1,  1, -1, -1,  1, -1,  1, -1,  1, -1,  1,  1, -1,  1, -1],\
[1,  1, -1, -1, -1, -1,  1,  1, -1, -1,  1,  1,  1,  1, -1, -1],\
[1, -1, -1,  1, -1,  1,  1, -1, -1,  1,  1, -1,  1, -1, -1,  1]]);


def test_correctness():
    """
    The correct Hadamard matrices have been hard coded into the source code. The
    fastwht reproduces all of these matrices, and computes the Frobenius norm of
    the difference between the constructed matrix and the hard coded matrix.
    If this norm is non-zero the test fails.
    """
    eps = 1e-8;

    for N in [8,16]:
        U_seq = zeros([N,N]);
        U_pal = zeros([N,N]);
        U_ord = zeros([N,N]);
        for i in range(N):
            x = zeros(N);
            x[i] = 1;

            y_seq = fastwht(x, order='sequency');
            y_pal = fastwht(x, order='dyadic');
            y_ord = fastwht(x, order='hadamard');

            U_seq[:,i] = y_seq;
            U_pal[:,i] = y_pal;
            U_ord[:,i] = y_ord;

        if (N == 8):
            zeroNorm_seq = linalg.norm(Had_seq8 - U_seq*N, 'fro');
            zeroNorm_pal = linalg.norm(Had_pal8 - U_pal*N, 'fro');
            zeroNorm_ord = linalg.norm(Had_ord8 - U_ord*N, 'fro');
            assert zeroNorm_seq < eps, "Wrong matrix N = 8" ;
            assert zeroNorm_pal < eps, "Wrong matrix N = 8" ;
            assert zeroNorm_ord < eps, "Wrong matrix N = 8" ;

        if (N == 16):
            zeroNorm_seq = linalg.norm(Had_seq16 - U_seq*N, 'fro');
            zeroNorm_pal = linalg.norm(Had_pal16 - U_pal*N, 'fro');
            zeroNorm_ord = linalg.norm(Had_ord16 - U_ord*N, 'fro');
            assert zeroNorm_seq < eps, "Wrong matrix N = 16" ;
            assert zeroNorm_pal < eps, "Wrong matrix N = 16" ;
            assert zeroNorm_ord < eps, "Wrong matrix N = 16" ;


def test_zero_expansion():
    """
    Various tests to verify that objects keep their shape and structure.
    """
    N = 16;
    eps = 1e-8;

    x = 5*random.randn(N,1)+5;
    y = fastwht(x);
    zeroNorm = linalg.norm(x-y);

    assert zeroNorm > eps, "The computations changed x";
    assert len(y.shape) == 2, "The shape of y changed";
    assert y.shape[0] == x.shape[0] and y.shape[1] == x.shape[1], "The shape changed";

    x = 5*random.randn(1,N)+5;
    y = fastwht(x);
    zeroNorm = linalg.norm(x-y);

    assert zeroNorm > eps, "The computations changed x";
    assert len(y.shape) == 2, "The shape of y changed";
    assert y.shape[0] == x.shape[0] and y.shape[1] == x.shape[1], "The shape changed";

    N = 6;
    x = zeros(N);
    x[1] = 1;
    y = fastwht(x);

    assert len(y.shape) == 1, "The shape of y changed";
    assert y.shape[0] == 8 and x.shape[0] == 6, "The shape changed";

    N = 8;
    x = zeros(N);
    x[1] = 1;
    y = fastwht(x,2);

    assert len(y.shape) == 1, "The shape of y changed";
    assert y.shape[0] == 2 and x.shape[0] == N, "The shape changed";

    success = False;
    try:
        fastwht(65);
    except TypeError:
        success = True;

    assert success, "Did accept non-array or non-matrix argument";

    success = False;
    try:
        x = zeros([N,N,N])
        fastwht(x);
    except TypeError:
        success = True;

    assert success, "Did accept more than 2-dimensional argument";

    success = False;
    try:
        x = zeros([N,N])
        fastwht(x);
    except TypeError:
        success = True;

    assert success, "Did accept more than 2-dimensional argument";

    N = 32;
    x = matrix(zeros(N));
    y = fastwht(x);
    assert isinstance(y, matrix), "No longer a matrix object";


def test_complex_numbers():
    N = 16;
    eps = 1e-10;

    x = zeros(N, dtype=complex128);
    for i in range(N):
        x[i] = i;
    x[3] = -4;
    x[7] = 20;

    x_copy = x.copy();
    y = fastwht(x);

    assert y.dtype == complex128, "Converted complex array into real"
    assert linalg.norm(y.imag) < eps, "Nonzero imaginary norm"
    assert linalg.norm(x - x_copy) < eps, "x has been changed"

    x = x - 1j*x;

    y = fastwht(x);

    assert linalg.norm(y.real + y.imag) < eps, "The two arrays should cancel"

    x = matrix(x);
    y = fastwht(x);

    assert isinstance(y, matrix), "Is not a matrix"
    assert y.shape[0] == x.shape[0] and y.shape[1] == x.shape[1], \
                                                "The shape changed";
    x.shape = (max(x.shape), min(x.shape));
    y = fastwht(x);
    assert y.shape[0] == x.shape[0] and y.shape[1] == x.shape[1], \
                                                "The shape changed";


def test_correctness_WAL_PAL():
    eps = 1e-8;

    N = 16;
    U_PAL = zeros([N,N]);
    U_WAL = zeros([N,N]);
    for n in range(N):
        for t in range(N):
            U_PAL[n,t] = PAL(N, n, t);
            U_WAL[n,t] = WAL(N, n, t);

    zeroNorm_seq = linalg.norm(Had_seq16 - U_WAL, 'fro');
    zeroNorm_pal = linalg.norm(Had_pal16 - U_PAL, 'fro');
    assert zeroNorm_seq < eps, "Wrong matrix N = 16" ;
    assert zeroNorm_pal < eps, "Wrong matrix N = 16" ;


if __name__ == "__main__":
    test_zero_expansion();
    test_correctness();
    test_complex_numbers();
    test_correctness_WAL_PAL();


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
            zeroNorm_seq = linalg.norm(Had_seq8 - U_seq*N);
            zeroNorm_pal = linalg.norm(Had_pal8 - U_pal*N);
            zeroNorm_ord = linalg.norm(Had_ord8 - U_ord*N);
            assert zeroNorm_seq < eps, "Wrong matrix N = 8" ;
            assert zeroNorm_pal < eps, "Wrong matrix N = 8" ;
            assert zeroNorm_ord < eps, "Wrong matrix N = 8" ;
        
        if (N == 16):
            zeroNorm_seq = linalg.norm(Had_seq16 - U_seq*N);
            zeroNorm_pal = linalg.norm(Had_pal16 - U_pal*N);
            zeroNorm_ord = linalg.norm(Had_ord16 - U_ord*N);
            assert zeroNorm_seq < eps, "Wrong matrix N = 16" ;
            assert zeroNorm_pal < eps, "Wrong matrix N = 16" ;
            assert zeroNorm_ord < eps, "Wrong matrix N = 16" ;

def test_zero_expansion():
    
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
    except ValueError:
        success = True;
    
    assert success, "Did accept non-array or non-matrix argument";
    
    success = False;
    try:
        x = zeros([N,N,N])
        fastwht(x);
    except ValueError:
        success = True;
    
    assert success, "Did accept more than 2-dimensional argument";
    
    success = False;
    try:
        x = zeros([N,N])
        fastwht(x);
    except ValueError:
        success = True;
    
    assert success, "Did accept more than 2-dimensional argument";
    
    
    
    
if __name__ == "__main__":
    test_zero_expansion();
    test_correctness(); 








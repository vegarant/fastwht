// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, version 3 of the License.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright 2017 Vegard Antun
//

#include "mex.h"
#include "../hadamard.h"
#include <iostream>
#include <cstring>
#include <complex>

/* Input Arguments */

#define	N_IN	prhs[0]
#define	FREQ_IN	prhs[1]
#define	K_IN	prhs[2]

/* Output Arguments */

#define	X_OUT	plhs[0]


inline bool isPowOf2(ulong x);



/*

Quick and dirty conversion from C++ to Matlab code. It is not thoroughly
tested yet. 

N = N_IN
n = FREQ_IN
k = K_IN

wal(N,n,k) computes the value w_n(k/N) for N = 2^nu for some nu >= 0, 
where 0 <= n,k < N

*/
void mexFunction( const int nlhs, mxArray *plhs[],
        		  const int nrhs, const mxArray *prhs[] ) {
    
    double N_double, n_double, k_double;
    /* Check for proper number of arguments */
    if (nlhs > 1) {
	    mexErrMsgIdAndTxt( "MATLAB:wal:maxlhs",
                "Too many output arguments.");
    }

    if (nrhs != 3 ) {
	    mexErrMsgIdAndTxt( "MATLAB:wal:invalidNumInputs",
                "Three input arguments are required");
    }
    
    /* Test the input type */
    if ( mxIsDouble(N_IN) and (!mxIsComplex(N_IN)) ) {
        N_double = mxGetScalar(N_IN);
    } else {
	    mexErrMsgIdAndTxt( "MATLAB:wal:unsupportedType",
                "N must be of data type 'double'");
    }
    
    if ( mxIsDouble(FREQ_IN) and (!mxIsComplex(FREQ_IN)) ) {
        n_double = mxGetScalar(FREQ_IN);
    } else {
	    mexErrMsgIdAndTxt( "MATLAB:wal:unsupportedType",
                "n must be of data type 'double'");
    }
    
    if ( mxIsDouble(FREQ_IN) and (!mxIsComplex(FREQ_IN)) ) {
        k_double = mxGetScalar(K_IN);
    } else {
	    mexErrMsgIdAndTxt( "MATLAB:wal:unsupportedType",
                "k must be of data type 'double'");
    }
    
    unsigned long N, n, k;
    double eps = 1e-15*N;
    
    // To avoid that N is truncated to N-1, due to round off errors.
    // std::round is only supported in C++11
    N = (unsigned long) N_double+eps;
    n = (unsigned long) n_double+eps;
    k = (unsigned long) k_double+eps;

    if (!isPowOf2(N)) {
	    mexErrMsgIdAndTxt( "MATLAB:wal:invalidInputArgument",
                           "N must be a power of 2, i.e. N = 2^x where x is positive integer");
    }
    
    if (k_double < 0 || n_double < 0 || N_double < 0) {
	    mexErrMsgIdAndTxt( "MATLAB:wal:invalidInputArgument",
                           "All arguments must be non-negative");
    }
    
    if (k >= N || n >= N) {
	    mexErrMsgIdAndTxt( "MATLAB:wal:invalidInputArgument",
                           "0 <= k,n < N");
    }
    
    int val = WAL(N,n,k);
    
    X_OUT       = mxCreateDoubleMatrix((mwSize) 1, (mwSize) 1, mxREAL);
    double * xd = mxGetPr(X_OUT);
    xd[0] = (double) val;
    
}

/*
   Is `x` on the form 2^nu;
*/
inline bool isPowOf2(ulong x)
{
    return  !(x & (x-1));
}

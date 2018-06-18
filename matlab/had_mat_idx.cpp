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
//



#include "mex.h"
#include "../hadamard.h"
#include "wal_comp_core.cpp"
#include <iostream>
#include <cstring>
#include <complex>

/* Input Arguments */

#define	N_IN	prhs[0]
#define	FREQ_IN	prhs[1]
#define	K_IN	prhs[2]

/* Output Arguments */

#define	X_OUT	plhs[0]

#if !defined(MAX)
#define	MAX(A, B)	((A) > (B) ? (A) : (B))
#endif


inline bool isPowOf2(ulong x);
inline void zeroOut(double * x, unsigned long N);


/*

N = N_IN
n = FREQ_IN
k = K_IN

Extract elements from a sequency ordered hadamard matrix.

>> W = N*fwht(eye(N));
>> W(n,k) % This is what the function computes

n and k can be arrays.

*/
void mexFunction( const int nlhs, mxArray *plhs[],
        		  const int nrhs, const mxArray *prhs[] ) {

    double N_double;
    double * n_double;
    double * k_double;

    unsigned long n_rows, n_cols, k_rows, k_cols = 0;
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
    if ( mxIsDouble(N_IN) && (!mxIsComplex(N_IN)) ) {
        N_double = mxGetScalar(N_IN);
    } else {
	    mexErrMsgIdAndTxt( "MATLAB:wal:unsupportedType",
                "N must be of data type 'double'");
    }


    if ( mxIsDouble(FREQ_IN) && (!mxIsComplex(FREQ_IN)) ) {
        n_rows = mxGetM(FREQ_IN);
        n_cols = mxGetN(FREQ_IN);
    } else {
	    mexErrMsgIdAndTxt( "MATLAB:wal:unsupportedType",
                           "n must be of data type 'double'");
    }

    if ( mxIsDouble(K_IN) && (!mxIsComplex(K_IN)) ) {
        k_rows = mxGetM(K_IN);
        k_cols = mxGetN(K_IN);
    } else {
	    mexErrMsgIdAndTxt( "MATLAB:wal:unsupportedType",
                "k must be of data type 'double'");
    }
    
    unsigned long N = (unsigned long) N_double;

    if (!isPowOf2(N)) {
	    mexErrMsgIdAndTxt( "MATLAB:wal:invalidInputArgument",
                           "N must be a power of 2, i.e. N = 2^x where x is positive integer");
    }
    
    if (k_rows != 1 && k_cols != 1) {
            mexErrMsgIdAndTxt( "MATLAB:fastwht:inputNotVector",
                               "Input must be a vector.");
    }
    
    if (n_rows != 1 && n_cols != 1) {
            mexErrMsgIdAndTxt( "MATLAB:fastwht:inputNotVector",
                               "Input must be a vector.");
    }
    
    const unsigned long n_size = MAX(n_cols, n_rows);
    const unsigned long k_size = MAX(k_cols, k_rows);
    
    //std::cout << "n_size: " << n_size << ", k_size: " << k_size << std::endl;
    
    /* Import matrix data  */
    n_double = mxGetPr(FREQ_IN);
    k_double = mxGetPr(K_IN);
    
    unsigned long *n_long = new unsigned long[n_size]; 
    unsigned long *k_long = new unsigned long[k_size]; 

    for (long i = 0; i < n_size; i++) {

        if (n_double[i] < 0 || n_double[i] > N) {
	        mexErrMsgIdAndTxt( "MATLAB:wal:invalidInputArgument",
                               "0 < k,n <= N");
        }

        if ( fabs(n_double[i]-floor(n_double[i])) > 1e-14*n_double[i] ) {
	        mexErrMsgIdAndTxt( "MATLAB:wal:invalidInputArgument",
                               "Walsh frequencies n, must be integers");
        }

        n_long[i] = (unsigned long) n_double[i]-1; // Convert from matlab to cpp indices

    }

    for (long i = 0; i < k_size; i++) {

        if (k_double[i] <= 0 || k_double[i] > N) {
	        mexErrMsgIdAndTxt( "MATLAB:wal:invalidInputArgument",
                               "0 < k,n <= N");
        }

        if ( fabs(k_double[i]-floor(k_double[i])) > 1e-14*k_double[i] ) {
            mexErrMsgIdAndTxt( "MATLAB:wal:invalidInputArgument",
                               "matrix indices must be integers");
        }

        k_long[i] = (unsigned long) k_double[i]-1; // Convert from matlab to cpp indices  

    }


    X_OUT       = mxCreateDoubleMatrix((mwSize) n_size, (mwSize) k_size, mxREAL);

    double * xd = mxGetPr(X_OUT);

    /* Get Walsh coefficients */
    wal_comp_core(n_long, n_size, k_long, k_size, xd, N, SEQUENCY);

    delete [] n_long;
    delete [] k_long;

}


































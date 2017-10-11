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


//                                  WARNING
//      THIS FUNCTION IS UNDER CONSTRUCTION AND IS NOT THOROUGHLY TESTED

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

#if !defined(MAX)
#define	MAX(A, B)	((A) > (B) ? (A) : (B))
#endif


inline bool isPowOf2(ulong x);
inline void zeroOut(double * x, unsigned long N);


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
    if ( mxIsDouble(N_IN) and (!mxIsComplex(N_IN)) ) {
        N_double = mxGetScalar(N_IN);
    } else {
	    mexErrMsgIdAndTxt( "MATLAB:wal:unsupportedType",
                "N must be of data type 'double'");
    }

    unsigned long N = (unsigned long) N_double;

    if (!isPowOf2(N)) {
	    mexErrMsgIdAndTxt( "MATLAB:wal:invalidInputArgument",
                           "N must be a power of 2, i.e. N = 2^x where x is positive integer");
    }
    
    if ( mxIsDouble(FREQ_IN) and (!mxIsComplex(FREQ_IN)) ) {
        n_rows = mxGetM(FREQ_IN);
        n_cols = mxGetN(FREQ_IN);
    } else {
	    mexErrMsgIdAndTxt( "MATLAB:wal:unsupportedType",
                "n must be of data type 'double'");
    }
    
    if ( mxIsDouble(K_IN) and (!mxIsComplex(K_IN)) ) {
        k_rows = mxGetM(K_IN);
        k_cols = mxGetN(K_IN);
    } else {
	    mexErrMsgIdAndTxt( "MATLAB:wal:unsupportedType",
                "k must be of data type 'double'");
    }
    
    if (k_rows != 1 and k_cols != 1) {
            mexErrMsgIdAndTxt( "MATLAB:fastwht:inputNotVector",
                               "Input must be a vector.");
    }
    
    if (n_rows != 1 and n_cols != 1) {
            mexErrMsgIdAndTxt( "MATLAB:fastwht:inputNotVector",
                               "Input must be a vector.");
    }
    
    const unsigned long n_size = MAX(n_cols, n_rows);
    const unsigned long k_size = MAX(k_cols, k_rows);
    
    //std::cout << "n_size: " << n_size << ", k_size: " << k_size << std::endl;
    
    /* Import matrix data  */
    n_double = mxGetPr(FREQ_IN);
    k_double = mxGetPr(K_IN);
    
    unsigned long n_long[n_size]; 
    
    for (long i = 0; i < n_size; i++) {

        if (n_double[i] < 0 or n_double[i] >= N) {
	        mexErrMsgIdAndTxt( "MATLAB:wal:invalidInputArgument",
                               "0 <= k,n < N");
        }

        if ( fabs(n_double[i]-floor(n_double[i])) > 1e-14*n_double[i] ) {
	        mexErrMsgIdAndTxt( "MATLAB:wal:invalidInputArgument",
                               "Walsh frequencies n, must be integers");
        }
        
        n_long[i] = (unsigned long) n_double[i];    
        
    }
    
    double k_max = 0;    
    double k_min = 0;    
    for (long i = 0; i < k_size; i++) {
        k_max =  (k_double[i] > k_max) ? k_double[i] : k_max; 
        k_min =  (k_double[i] < k_min) ? k_double[i] : k_min; 
    }
    
    bool is_in_0_1 = (k_max < 1.0);
    
    if (k_max >= N_double or k_min < 0) {
	    mexErrMsgIdAndTxt( "MATLAB:wal:invalidInputArgument",
                           "0 <= k,n < N");
    }

    unsigned long k_long[k_size];

    if( is_in_0_1) {
        for (unsigned long i = 0; i < k_size; i++) {
            k_long[i] = (unsigned long) floor(N_double*k_double[i]);
        }
    } else {
        
        for (unsigned long i = 0; i < k_size; i++) {
            if ( fabs(k_double[i]-floor(k_double[i])) > 1e-14*k_double[i] ) {
	            mexErrMsgIdAndTxt( "MATLAB:wal:invalidInputArgument",
                                   "Walsh indices k, must be integers or in [0,1]");
            }
        
            k_long[i] = (unsigned long)k_double[i];
        }
    
    }
    
    X_OUT       = mxCreateDoubleMatrix((mwSize) n_size, (mwSize) k_size, mxREAL);

    double * xd = mxGetPr(X_OUT);
    
    if (k_size == 1 and n_size == 1) {
        
        int val = WAL( N, n_long[0], k_long[0] );
        xd[0] = (double) val;
        
        //std::cout << "n_size: " << n_size << ", k_size: " << k_size << std::endl;
        //std::cout << "n_long[0]: " << n_long[0] << std::endl;
        //std::cout << "k_long[0]: " << k_long[0] << std::endl;

        //std::cout << "n_double[0]: " << n_double[0] << std::endl;
        //std::cout << "k_double[0]: " << k_double[0] << std::endl;

    } else if (k_size == 1) {

        double * x  = new double[N];
        zeroOut(x,N);

        x[k_long[0]] = 1;
        hadamardTransform<double>(x, N, SEQUENCY);

        for (unsigned long i = 0; i < n_size; i++) {
            xd[i] = x[n_long[i]];
        }

        delete [] x;    
        
    } else if (n_size == 1) {
        
        double * x  = new double[N];
        zeroOut(x,N);
        
        x[n_long[0]] = 1;
        hadamardTransform<double>(x, N, SEQUENCY);

        for (unsigned long i = 0; i < k_size; i++) {
            xd[i] = x[k_long[i]];
        }

        delete [] x;    
        
    } else { // n_size > 1 and k_size > 1

        double * x  = new double[N];
        zeroOut(x,N);
        
        if (n_size > k_size) {
            
            for (long j = 0; j < k_size; j++ ) {
                x[k_long[j]] = 1;
                hadamardTransform<double>(x, N, SEQUENCY);
                for (unsigned long i = 0; i < n_size; i++) {
                    xd[j*n_size + i] = x[n_long[i]];
                }
                zeroOut(x,N);
            }
            
        } else { // k_size > n_size
            
            for (long j = 0; j < n_size; j++ ) {
                x[n_long[j]] = 1;
                hadamardTransform<double>(x, N, SEQUENCY);
                for (unsigned long i = 0; i < k_size; i++) {
                    xd[i*n_size + j] = x[k_long[i]];
                }
                zeroOut(x,N);
            }
        
        }

        delete [] x;    

    }

}


/*
   Is `x` on the form 2^nu;
*/
inline bool isPowOf2(ulong x)
{
    return  !(x & (x-1));
}

inline void zeroOut(double * x, unsigned long N) {
    for (unsigned long i = 0; i < N; i++) {
        x[i] = 0;
    }
}
































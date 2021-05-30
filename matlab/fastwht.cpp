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
// Copyright 2016 Vegard Antun
//

#include "mex.h"
#include "../hadamard.h"
#include <iostream>
#include <cstring>
#include <complex>
// Uncomment the line below as well as the lines in the function hadamardTransformColumnWise to enable openMP (remember to compile correctly).
//#include <omp.h>

/* Input Arguments */

#define	X_IN	prhs[0]
#define	LENGTH_IN	prhs[1]
#define	ORDER_IN	prhs[2]

/* Output Arguments */

#define	X_OUT	plhs[0]

#if !defined(MAX)
#define	MAX(A, B)	((A) > (B) ? (A) : (B))
#endif

#if !defined(MIN)
#define	MIN(A, B)	((A) < (B) ? (A) : (B))
#endif

enum MatlabDataType {INT32, DOUBLE, COMPLEX};
inline unsigned long getLegalLength(const unsigned long N);
inline bool isPowOf2(ulong x);

void copyData(double * x, const double *data, const unsigned long M_data,
                                              const unsigned long N_data, 
                                              const unsigned long M, 
                                              const unsigned long N,
                                              const bool useOneDimTransform);

void hadamardTransformColumnWise(double * x, const unsigned long M, 
                                             const unsigned long N, 
                                             const HadamardOrder order);


void mexFunction( const int nlhs, mxArray *plhs[],
        		  const int nrhs, const mxArray *prhs[] ) {

    /* Set default order */
    HadamardOrder order = SEQUENCY;

    /* Define variables */
    MatlabDataType dataType;
    const char *inputBuf;
    bool useOneDimTransform = true;

    // Only one of these set of pointers will be used.
    double *xd, *xd_out;
    double *xd_im, *xd_out_im;
    int *xi, *xi_out;

    double input_N = 1;

    /* Check for proper number of arguments */
    if (nlhs > 1) {
	    mexErrMsgIdAndTxt( "MATLAB:fastwht:maxlhs",
                "Too many output arguments.");
    }

    if (0 == nrhs || nrhs > 3 ) {
	    mexErrMsgIdAndTxt( "MATLAB:fastwht:invalidNumInputs",
                "one, two, or three input arguments are required");
    }

    const unsigned long M_data = mxGetM(X_IN);
    const unsigned long N_data = mxGetN(X_IN);
    unsigned long M = M_data;
    unsigned long N = N_data;
    
    /* Test the input type */
    if ( mxIsDouble(X_IN) && (!mxIsComplex(X_IN)) ) {
        dataType = DOUBLE;
    } else if ( mxIsComplex(X_IN) ) {
        dataType = COMPLEX;
    } else {
	    mexErrMsgIdAndTxt( "MATLAB:fastwht:unsupportedType",
                "The vector must be of type double or complex");
    }

    if (nrhs >= 2) {

        const unsigned long numElementsLength = mxGetNumberOfElements(LENGTH_IN);

        switch (numElementsLength) 
        {

            case 1:
            {
                if ( !mxIsDouble(LENGTH_IN) || mxIsComplex(LENGTH_IN) ) {
                    mexErrMsgIdAndTxt("MyToolbox:fastwht:notDouble",
                                      "Input length must be type double.");
                }

                if (N_data == 1) {
                    M = mxGetScalar(LENGTH_IN);
                    N = 1;
                } else if (M_data == 1) {
                    N = mxGetScalar(LENGTH_IN);
                    M = 1;
                } else {
                    M = mxGetScalar(LENGTH_IN);
                    N = N_data;
                }
                if ((! isPowOf2(M)) || (! isPowOf2(N)) ) {
                    mexErrMsgIdAndTxt("MyToolbox:fastwht:notPowerOf2",
                  "Input must be a non-negative integer which is a power of 2");
                }
                break;
            }
            case 0: 
            {
                M = M_data;
                N = N_data;
                break;
            }
            case 2:
            {
                if ( !mxIsDouble(LENGTH_IN) || mxIsComplex(LENGTH_IN) ) {
                    mexErrMsgIdAndTxt("MyToolbox:fastwht:notDouble",
                                      "Input length must be type double.");
                }
                useOneDimTransform = false;
                const double *pLengthIn = mxGetPr(LENGTH_IN);
                M = (long) pLengthIn[0];
                N = (long) pLengthIn[1];
                
                if ((! isPowOf2(M)) || (! isPowOf2(N)) ) {
                    mexErrMsgIdAndTxt("MyToolbox:fastwht:notPowerOf2",
                  "Input must be a non-negative integer which is a power of 2");
                }
                
                break;
            }
            default :
	        {
                mexErrMsgIdAndTxt( "MATLAB:fastwht:invalidNumInputs",
                  "one, two, or three input arguments are required");
            }
        }
    }
    
    M = getLegalLength(M);
    N = getLegalLength(N);
    
    // Deduce the order of the transform 
    if (nrhs == 3) {
        /* Input must be a string */
        if ( mxIsChar( ORDER_IN ) != 1) {
            mexErrMsgIdAndTxt( "MATLAB:fastwht:inputNotString",
                             "Input must be a string.");
        }
        /* Input must be a row vector */
        if (mxGetM( ORDER_IN ) != 1) {
            mexErrMsgIdAndTxt( "MATLAB:fastwht:inputNotVector",
                             "Input must be a row vector.");
        }

        /* copy the string data from prhs[2] into a C string input_ buf.    */
        inputBuf = mxArrayToString(ORDER_IN);

        if( inputBuf == NULL ) {
            mexErrMsgIdAndTxt( "MATLAB:fastwht:conversionFailed",
                  "Could not convert input to string.");
        }

        if (std::strcmp(inputBuf, "sequency") == 0) {
            order = SEQUENCY;
        } else if (std::strcmp(inputBuf, "dyadic") == 0) {
            order = PALEY;
        } else if (std::strcmp(inputBuf, "hadamard") == 0) {
            order = ORDINARY;
        } else {
            mexErrMsgIdAndTxt( "MATLAB:fastwht:unknownOrder",
                  "Did not recognize the order option.");
        }
    }

    if (dataType == DOUBLE) {

        X_OUT  = mxCreateDoubleMatrix((mwSize) M, (mwSize) N, mxREAL);

    } else if (dataType == COMPLEX) {

        X_OUT  = mxCreateDoubleMatrix((mwSize) M, (mwSize) N, mxCOMPLEX);

        xd_im     = mxGetPi(X_IN);
        xd_out_im = mxGetPi(X_OUT);

        copyData(xd_out_im, xd_im, M_data, N_data, M, N, useOneDimTransform);    

        if ( useOneDimTransform ) {
            hadamardTransformColumnWise(xd_out_im, M, N, order); 
        } else {
            hadamardTransform2dColumn<double>(xd_out_im, M, N, order);
        }

    } else {
	    mexErrMsgIdAndTxt( "MATLAB:fastwht:unsupportedType",
                "The vector must be of type double or complex");
    }

    xd        = mxGetPr(X_IN);
    xd_out    = mxGetPr(X_OUT);

    copyData(xd_out, xd, M_data, N_data, M, N, useOneDimTransform);

    if ( useOneDimTransform ) {
        hadamardTransformColumnWise(xd_out, M, N, order); 
    } else {
        hadamardTransform2dColumn<double>(xd_out, M, N, order);
    }

    return;
}


inline unsigned long getLegalLength(const unsigned long N) {
    unsigned long M = 1UL;
    while( M < N ) {
        M = M << 1;
    }
    return M;
}


void copyData(double * x, const double *data, const unsigned long M_data,
                                              const unsigned long N_data, 
                                              const unsigned long M, 
                                              const unsigned long N,
                                              const bool useOneDimTransform)
{
    double factor = 1;
    if (useOneDimTransform) {
        const unsigned long G = (M > N) ? M : N;
        factor = 1.0/G;
    } else {
        factor = 1.0/(M*N);
    }
    
    const double constFactor = factor;
    const unsigned long minN = (N < N_data) ? N : N_data;
    const unsigned long minM = (M < M_data) ? M : M_data;

    for (int j = 0; j < minN; j++) {
        for (int i = 0; i < minM; i++) {
            x[M*j + i] = data[M_data*j + i]*constFactor;
        }
    }
    
    
}

void hadamardTransformColumnWise(double* x, const unsigned long M, const unsigned long N, 
                                          const HadamardOrder order) 
{
    if (M == 1 || N == 1) {
        const unsigned long G = (M > N) ? M : N;
        hadamardTransform<double>(x, G, order);
    } else {
// Uncomment these lines to use OpenMP
//        #pragma omp parallel
//        {
//        #pragma omp for
        for (int i = 0; i < N; i++) {
            
            hadamardTransform<double>(&x[i*M], M, order);
        }
//        }
    }
}

/*
   Is `x` on the form 2^nu;
*/
inline bool isPowOf2(ulong x)
{
    return  !(x & (x-1));
}




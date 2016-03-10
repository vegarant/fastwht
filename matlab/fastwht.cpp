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

enum Type {INT32, DOUBLE, COMPLEX};
enum Order {ORDINARY, SEQUENCY, PALEY};


static void fastwht(int  N, double *arr) {
     hadamardOrdinary<double>(arr, N);    
} 

void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray *prhs[] )
{ 
    
    Type dataType;
    Order hadamardOrder = SEQUENCY;
    const char *inputBuf;    
    
    // Only one of these set of pointers will be used.
    double *xd, *xd_out;  
    double *xd_im, *xd_out_im; 
    int *xi, *xi_out; 
    
    double input_N = 1;
    //int sizeOfDataType = 0;
    size_t M,N; 
    
    /* Check for proper number of arguments */
    if (nlhs > 1) {
	    mexErrMsgIdAndTxt( "MATLAB:fastwht:maxlhs",
                "Too many output arguments."); 
    }
 
    if (nrhs > 3) { 
	    mexErrMsgIdAndTxt( "MATLAB:fastwht:invalidNumInputs",
                "No more that tree input arguments are allowed."); 
    }
    
    /* Check the dimensions of X.  X can be N X 1 or 1 X N. */ 
    M = mxGetM(X_IN); 
    N = mxGetN(X_IN);
    
    if ( (MIN(M,N) != 1) ) {  //mxIsComplex(Y_IN) || 
	    mexErrMsgIdAndTxt( "MATLAB:fastwht:toManyDimensions",
                "The vector can only be one dimensional"); 
    } 
    

    /* Test the input type */
    if ( mxIsDouble(X_IN) and not mxIsComplex(X_IN) ) {
        dataType = DOUBLE;
    //    sizeOfDataType = sizeof(double);
    } else if ( mxIsComplex(X_IN) ) {
        dataType = COMPLEX;
    //    sizeOfDataType = sizeof(std::complex<double>);
    } else if ( mxIsInt32(X_IN) ) {
        dataType = INT32;
    //    sizeOfDataType = sizeof(int);
    } else {
	    mexErrMsgIdAndTxt( "MATLAB:fastwht:unsupportedType",
                "The vector must be of type int32, double or complex"); 
    }
     
    size_t originalVectorLength = MAX(M,N);
     
    if (nrhs >= 2) {
         
        if ( !mxIsDouble(LENGTH_IN) || mxIsComplex(LENGTH_IN) ) {
            mexErrMsgIdAndTxt("MyToolbox:fastwht:notDouble",
                              "Input length must be type double.");
        }
        
        if ( mxGetNumberOfElements(LENGTH_IN) > 1 ) {
            mexErrMsgIdAndTxt("MyToolbox:fastwht:notScalar",
                              "The secound input must be scalar");
            
        }
        
        if (mxGetNumberOfElements(LENGTH_IN) == 0) {
            
            
            
        }  else {
            
            /* Get the value of the scalar input  */
            input_N = mxGetScalar(LENGTH_IN);
             
            size_t k = 1;
             
            while ( (input_N - k) > 1e-5 ) {
                k *= 2; 
            }
             
            if (N > M) {
                N = k; 
            } else {  
                M = k; 
            }
        }
         
    }
    
    
    
    if (nrhs == 3) {
        /* input must be a string */
        if ( mxIsChar( ORDER_IN ) != 1) {
            mexErrMsgIdAndTxt( "MATLAB:fastwht:inputNotString",
                             "Input must be a string.");
        }
        /* input must be a row vector */
        if (mxGetM( ORDER_IN )!=1) {
            mexErrMsgIdAndTxt( "MATLAB:fastwht:inputNotVector",
                             "Input must be a row vector.");
        }
        ///* get the length of the input string */
        //buflen = ( mxGetM( ORDER_IN ) * mxGetN( ORDER_IN) ) + 1;
         
        ///* allocate memory for output string */
        //output_buf=mxCalloc(buflen, sizeof(char));
         
        /* copy the string data from prhs[0] into a C string input_ buf.    */
        inputBuf = mxArrayToString(ORDER_IN);
        
        if(inputBuf == NULL) {
            mexErrMsgIdAndTxt( "MATLAB:fastwht:conversionFailed",
                  "Could not convert input to string.");
        } 
        
        if (std::strcmp(inputBuf, "sequency") == 0) {
            hadamardOrder = SEQUENCY;    
        } else if (std::strcmp(inputBuf, "dyadic") == 0) {
            hadamardOrder = PALEY;    
        } else if (std::strcmp(inputBuf, "hadamard") == 0) {
            hadamardOrder = ORDINARY;    
        } else {
            mexErrMsgIdAndTxt( "MATLAB:fastwht:unknownOrder",
                  "Did not recognize the order option.");
        }
    }
     
    /* create the output matrix */
    size_t newVectorLength = MAX(M,N);
    const double newVectorLengthDouble = (double) newVectorLength;
    originalVectorLength = (originalVectorLength > newVectorLength) ? newVectorLength : originalVectorLength;

    if (dataType == DOUBLE) {
        X_OUT  = mxCreateDoubleMatrix((mwSize) M, (mwSize)N ,mxREAL);
        xd     = mxGetPr(X_IN); 
        xd_out = mxGetPr(X_OUT); 
         
        for (int i = 0; i < originalVectorLength; i++) {
            xd_out[i] = xd[i]/newVectorLengthDouble;
        }
        for (int i = originalVectorLength; i < newVectorLength; i++) {
            xd_out[i] = 0;
        }
        
        if (hadamardOrder == ORDINARY) {
           hadamardOrdinary<double>(xd_out, newVectorLength);
        } else if (hadamardOrder == PALEY) {
           hadamardPaley<double>(xd_out, newVectorLength);
        } else if (hadamardOrder == SEQUENCY) {
           hadamardSequency<double>(xd_out, newVectorLength);
        }  

    } else if (dataType == COMPLEX) {
        X_OUT     = mxCreateDoubleMatrix((mwSize) M, (mwSize)N ,mxCOMPLEX);
         
        xd        = mxGetPr(X_IN); 
        xd_out    = mxGetPr(X_OUT); 
         
        xd_im     = mxGetPi(X_IN); 
        xd_out_im = mxGetPi(X_OUT); 
         
        std::complex<double> *c_vector =   new std::complex<double> [newVectorLength]; 
        for (int i = 0; i < originalVectorLength; i++) {
            c_vector[i] = std::complex<double>( xd[i], xd_im[i] )/newVectorLengthDouble;
        }
        for (int i = originalVectorLength; i < newVectorLength; i++) {
            c_vector[i] = 0;
        }
         
        if (hadamardOrder == ORDINARY) {
           hadamardOrdinary<std::complex<double> >(c_vector, newVectorLength);
        } else if (hadamardOrder == PALEY) {
           hadamardPaley<std::complex<double> >(c_vector, newVectorLength);
        } else if (hadamardOrder == SEQUENCY) {
           hadamardSequency<std::complex<double> >(c_vector, newVectorLength);
        } 
         
        for (int i = 0; i < newVectorLength; i++) {
             
            xd_out[i] = c_vector[i].real();
            xd_out_im[i] = c_vector[i].imag();
             
        }
         
         
        delete [] c_vector;
    } else if (dataType == INT32) {
         
        mwSize dim[2] = {M, N};
         
        X_OUT = mxCreateNumericArray(2, dim, mxINT32_CLASS, mxREAL); 
         
        xi     = (int*) mxGetData(X_IN); 
        xi_out = (int*) mxGetData(X_OUT); 
        
        // Zero pad vector  
        std::memcpy(xi_out, xi, sizeof(int)*newVectorLength);
        if ( originalVectorLength - newVectorLength ) {
            std::memset(&xi_out[originalVectorLength], 0, 
                        sizeof(int)*(newVectorLength-originalVectorLength));
        }
        
        
        if (hadamardOrder == ORDINARY) {
            
           hadamardOrdinary<int>(xi_out, newVectorLength);
            
        } else if (hadamardOrder == PALEY) {
             
            hadamardPaley<int>(xi_out, newVectorLength);
             
        } else if (hadamardOrder == SEQUENCY) {
             
            hadamardSequency<int>(xi_out, newVectorLength);
             
        } 
         
    } 
     
     
    return;
     
}







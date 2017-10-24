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

/*

This file contain the interface to the computational core of the Hadamard transform.

*/

#include "hadamard.h"
#include "python/hadamardKernel.h"


/*
    The Hadamard transform

    INPUT:
    x - The vector one would like to transform
    N - Length of the vector
    order - The order of the Hadamard transform

*/
template <typename T>
void hadamardTransform(T* x, const unsigned long N, const HadamardOrder order) 
{
    
    const unsigned long ldn = ld(N);
    
    switch (order) 
    {
        case SEQUENCY :
            walsh_wal<T>(x,ldn);
            break;
        
        case PALEY :
            revbin_permute<T>(x, N);
            walsh_wak(x, ldn);
            break;
        
        case ORDINARY :
            walsh_wak(x, ldn);
            break; 
    }
}

template void hadamardTransform<>(short * x, const unsigned long N, const HadamardOrder order); 
template void hadamardTransform<>(int * x, const unsigned long N, const HadamardOrder order); 
template void hadamardTransform<>(long * x, const unsigned long N, const HadamardOrder order); 
template void hadamardTransform<>(float * x, const unsigned long N, const HadamardOrder order); 
template void hadamardTransform<>(double * x, const unsigned long N, const HadamardOrder order); 
template void hadamardTransform<>(long double * x, const unsigned long N, const HadamardOrder order); 
template void hadamardTransform<>(std::complex<float> * x, const unsigned long N, const HadamardOrder order); 
template void hadamardTransform<>(std::complex<double> * x, const unsigned long N, const HadamardOrder order); 
template void hadamardTransform<>(std::complex<long double> * x, const unsigned long N, const HadamardOrder order); 


/*
    A two-dimensional Hadamard transform using a tensor product expansion

    This function will apply a Hadamard transform to every column and row of
    input matrix. 

    INPUT:
    x - Two-dimensional matrix, with column-major order
    M - Number of rows
    N - Number of columns
    order - The order of the Hadamard transform

*/
template<typename T> 
void hadamardTransform2dColumn(T* x, const unsigned long M , 
                                     const unsigned long N, 
                                     const HadamardOrder order) 
{
    
    T * y = new T[N]; // Intermediate array
    
    // For each column perform the Hadamard transform
    for (unsigned int k = 0; k < N; k++) {
        hadamardTransform<T>(&x[k*M], M, order);
    }

    // For each row apply the Hadamard transform
    for (unsigned int k = 0; k < M; k++) {
        // Copy out the row column
        for (unsigned int i = 0; i < N; i++) {
            y[i] = x[i*M+k];
        } 
         
        hadamardTransform<T>(y, N, order);
         
        // Copy the calculated elements back in the array
        for (unsigned int i = 0; i < N; i++) {
            x[i*M+k] = y[i];
        } 
    }
     
    delete [] y;
     
}

template void hadamardTransform2dColumn<>(double* x, const unsigned long M , 
                                               const unsigned long N, 
                                               const HadamardOrder order); 

/*

Various kernels used for the python binding.

*/
void fwhtKernelSequency(int n, double *arr) 
{
     hadamardTransform<double>(arr, n, SEQUENCY);
}

void fwhtKernelOrdinary(int n, double *arr) 
{
     hadamardTransform<double>(arr, n, ORDINARY);
}

void fwhtKernelPaley(int n, double *arr) 
{
     hadamardTransform<double>(arr, n, PALEY);
}

int PAL_kernel(unsigned int N, unsigned int n, unsigned int x) 
{
    return PAL(N,n,x);
}

int WAL_kernel(unsigned int N, unsigned int n, unsigned int x) 
{
    return WAL(N,n,x);
}

/*

Computes the element in the sequency ordered Walsh-Hadamard matrix.

N - Dimension of the matrix N × N
n - row number i.e., Walsh-Hadamard function ψ_n
t - column i.e., the input ψ_n(t/N)

*/
int WAL(unsigned long N, unsigned long n, unsigned long t) 
{
    const unsigned long ldn = ld(N);
    unsigned long s = 0;
    unsigned long n_pr = 1;
    unsigned long t_r  = 1;
    n_pr = n_pr << (ldn-1);
    unsigned long n_r = 0;
    unsigned long t_s = 0;

    for (unsigned long r = 0; r < ldn; r++) 
    {
        n_r = (((n_pr >> r) & n) >> (ldn - r-1));
        t_s = (((t_r << r) & t) >> r) - (((t_r << (r+1)) & t) >> (r+1));
        s += n_r*t_s;
    }

    if ( s%2 ) { // s is odd
        return -1;
    } else {
        return 1;
    }
}

/*

Computes the element in the Paley ordered Walsh-Hadamard matrix.

N - Dimension of the matrix N × N
n - row number i.e., Walsh-Hadamard function ψ_n
t - column i.e., the input ψ_n(t/N)

*/
int PAL(unsigned long N, unsigned long n, unsigned long x) 
{
    const unsigned int ldn = ld(N); 
    unsigned long s = 0;
    unsigned long ONE = 1;
    unsigned long n_j = 0;
    unsigned long x_jp1 = 0;

    for (unsigned long j = 0; j < ldn; j++) 
    {
        n_j = (n & (ONE << j)) >> j ;
        x_jp1 = (x & (ONE << (ldn - j-1))) >> (ldn - j-1);

        s += n_j*x_jp1;
    }

    if ( s%2 ) { // s is odd
        return -1;
    } else {
        return 1;
    }
}










//******************************************************************************
//******************************************************************************
//***                           Outdated functions                           ***
//******************************************************************************
//******************************************************************************
// These functions are primarily included for the use with the verification code

/*

Perform the fast Walsh-Hadamard transform in sequency order.

*/
template <typename T>
void hadamardSequency(T * x, const unsigned long N) {
    hadamardTransform<T>(x, N, SEQUENCY);
}

template void hadamardSequency<>(short * x, const unsigned long N);
template void hadamardSequency<>(int * x, const unsigned long N);
template void hadamardSequency<>(long * x, const unsigned long N);
template void hadamardSequency<>(float * x, const unsigned long N);
template void hadamardSequency<>(double * x, const unsigned long N);
template void hadamardSequency<>(long double * x, const unsigned long N);
template void hadamardSequency<>(std::complex<double>* x, const unsigned long N);

/*

Perform the fast Walsh-Hadamard transform in Paley order.

*/
template <typename T>
void hadamardPaley(T * x, const unsigned long N) {
    hadamardTransform<T>(x, N, PALEY);
}

template void hadamardPaley<>(short * x, const unsigned long N);
template void hadamardPaley<>(int * x, const unsigned long N);
template void hadamardPaley<>(long * x, const unsigned long N);
template void hadamardPaley<>(float * x, const unsigned long N);
template void hadamardPaley<>(double * x, const unsigned long N);
template void hadamardPaley<>(long double * x, const unsigned long N);
template void hadamardPaley<>(std::complex<double>* x, const unsigned long N);


/*

Performs the Walsh-Hadamard transform using the ordinary order.
The calculations happens in-place.

INPUT:
x - The vector one would like to transform
N - Length of the vector. N must be a power of 2.

*/
template <typename T>
void hadamardOrdinary(T *x, const unsigned long N) 
{
    hadamardTransform<T>(x, N, ORDINARY);
}

// Specialization
template void hadamardOrdinary<>(int* x, const unsigned long N);
template void hadamardOrdinary<>(short* x, const unsigned long N);
template void hadamardOrdinary<>(long* x, const unsigned long N);
template void hadamardOrdinary<>(float* x, const unsigned long N);
template void hadamardOrdinary<>(double* x, const unsigned long N);
template void hadamardOrdinary<>(long double* x, const unsigned long N);
template void hadamardOrdinary<>(std::complex<double>* x, const unsigned long N);


void fwhtKernelSequency(int n, long *arr) 
{
     hadamardTransform<long>(arr, n, SEQUENCY);
}

void fwhtKernelOrdinary(int n, long *arr) 
{
     hadamardTransform<long>(arr, n, ORDINARY);
}

void fwhtKernelPaley(int n, long *arr) 
{
     hadamardTransform<long>(arr, n, PALEY);
}

// DIM 1 == n, DIM2 == m
template <typename T>
void fwhtKernel2dColumn(int n, int m, T *arr, const HadamardOrder order){
    T * x_tmp = new T[n];
    for (int j = 0; j < m; j++) {
        for (int i = 0; i < n; i++) {
            x_tmp[i] = arr[i*m+j];
        }
        
        hadamardTransform<T>(x_tmp, n, order);            
        
        for (int i = 0; i < n; i++) {
            arr[i*m+j] = x_tmp[i];
        }
    }
    delete [] x_tmp;
}
template void fwhtKernel2dColumn<>(int n, int m, double *arr, const HadamardOrder order);
template void fwhtKernel2dColumn<>(int n, int m, long   *arr, const HadamardOrder order);


void fwhtKernel2dSequency(int n, int m, double *arr) {
    fwhtKernel2dColumn<double>(n, m, arr, SEQUENCY); 
}

void fwhtKernel2dOrdinary(int n, int m, double *arr) {
    fwhtKernel2dColumn<double>(n, m, arr, ORDINARY); 
}
void fwhtKernel2dPaley(int n, int m, double *arr) {
    fwhtKernel2dColumn<double>(n, m, arr, PALEY); 
}

void fwhtKernel2dSequency(int n, int m, long *arr) {
    fwhtKernel2dColumn<long>(n, m, arr, SEQUENCY); 
}

void fwhtKernel2dOrdinary(int n, int m, long *arr) {
    fwhtKernel2dColumn<long>(n, m, arr, ORDINARY); 
}

void fwhtKernel2dPaley(int n, int m, long *arr) {
    fwhtKernel2dColumn<long>(n, m, arr, PALEY); 
}






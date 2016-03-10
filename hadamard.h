/*

This header file is used to interface with any C++ implementation

For the file interfacing with python see python/hadamard.h

*/

#include <iostream> 
#include <iomanip>
#include <cstring>
#include <cmath>
#include <complex>
#include <vector>

unsigned int sequency(unsigned int a, unsigned int N);

int findMostSignificantBit(unsigned int a);

/*

Computes the element in the sequency ordered Walsh-Hadamard matrix. 

N - Dimension of the matrix N × N
n - row number i.e., Walsh-Hadamard function ψ_n
t - column i.e., the input ψ_n(t/N)

*/
int WAL(unsigned int N, unsigned int n, unsigned int t);
int PAL(unsigned int N, unsigned int n, unsigned int t);

template <typename T>
void hadamardOrdinary(T *x, unsigned int N);

template <typename T> 
void hadamardSequency(T * x, unsigned int N);

template <typename T> 
void hadamardPaley(T * x, unsigned int N);

unsigned int reverseBitSequence(const unsigned int N, unsigned int x);

unsigned int binaryToGrayCode( unsigned int x );
unsigned int grayCodeToBinary(unsigned int x);













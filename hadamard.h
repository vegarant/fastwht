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

This header file is used to interface with any C++ implementation

For the file interfacing with python see python/hadamard.h

*/

#include <iostream> 
#include <iomanip>
#include <cstring>
#include <cmath>
#include <complex>
#include <vector>


unsigned int powDyadic(const unsigned int k); 
int findMostSignificantBit(unsigned int a);

/*

Computes the element in the sequency ordered Walsh-Hadamard matrix. 

N - Dimension of the matrix N × N
n - row number i.e., Walsh-Hadamard function ψ_n
t - column i.e., the input ψ_n(t/N)

*/
int WAL( unsigned int N, unsigned int n, unsigned int t);
int PAL(unsigned int N, unsigned int n, unsigned int t);

template <typename T>
void hadamardOrdinary(T *x, const unsigned int N);

template <typename T> 
void hadamardSequency(T * x, const unsigned int N);

template <typename T> 
void hadamardPaley(T * x, const unsigned int N);

unsigned int reverseBit(const unsigned int N, const unsigned int x);
unsigned int reverseBitSequence(const unsigned int N, unsigned int x);

unsigned int binaryToGrayCode( unsigned int x );
unsigned int grayCodeToBinary(unsigned int x);













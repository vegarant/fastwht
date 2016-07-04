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

#ifndef FASTWHT_HADAMARD_H_
#define FASTWHT_HADAMARD_H_

#include <iostream>
#include <iomanip>
#include <cstring>
#include <cmath>
#include <complex>
#include <vector>
#include <cstdint>
#include <thread>


uint32_t powDyadic(const uint32_t k);
unsigned int findMostSignificantBit(uint32_t a);

/*

Computes the element in the sequency ordered Walsh-Hadamard matrix.

N - Dimension of the matrix N × N
n - row number i.e., Walsh-Hadamard function ψ_n
t - column i.e., the input ψ_n(t/N)

*/
int WAL(unsigned int N, unsigned int n, unsigned int t);
int PAL(unsigned int N, unsigned int n, unsigned int t);

template <typename T>
void hadamardOrdinary(T *x, const uint32_t N);

template <typename T>
void hadamardSequency(T * x, const uint32_t N);

template <typename T>
void hadamardPaley(T * x, const uint32_t N);

template <typename T>
void hadamardRecursive(T *x, const unsigned int N);

template <typename T>
void hadamardDepthFirst(T *x, const unsigned int N);

template <typename Type>
void hadamardArndt(Type *f, const uint32_t ldn);

template <typename T>
void hadamardArndtOneStep(T *x, const uint32_t N, const uint32_t ldm );

template <typename T>
void hadamardParallel(T* x, const uint32_t N);

uint32_t reverseBitSequence(const uint32_t N, uint32_t x);

uint32_t idxFromOrdinaryToSequency(uint32_t a, uint32_t N);
uint32_t binaryToGrayCode(uint32_t x );
uint32_t grayCodeToBinary(uint32_t x);





#endif
